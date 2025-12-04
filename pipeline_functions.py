import matplotlib.pyplot as plt
import numpy as np
import wfdb
import pandas as pd
from scipy.signal import firwin, filtfilt, resample_poly
import os
from adaptivePLICanceller import AdaptivePLICanceller
import math
from MECGCanceller import MECGCanceller
from FetalECGExtractor import FetalECGExtractor
from sklearn.decomposition import PCA
import os
import numpy as np

def load_and_interpolate(record_path):
    """
    Carica segnale WFDB e interpola eventuali NaN.
    """
    try:
        signals, fields = wfdb.rdsamp(record_path)
    except Exception as e:
        raise FileNotFoundError(f"Impossibile leggere {record_path}: {e}")

    df = pd.DataFrame(signals)
    if df.isnull().values.any():
        df = df.interpolate(method='linear', limit_direction='both')
        signals = df.values
    return signals, fields

def design_baseline_filter(fs, fc_high=3.0, taps=1001):
    """
    Progetta il filtro FIR passa-alto per il baseline wander
    """
    if taps % 2 == 0: taps += 1
    b_fir = firwin(taps, fc_high, fs=fs, pass_zero=False, window='hamming')
    a_fir = [1.0]
    return b_fir, a_fir

def resample_ecg(signal, fs_old, fs_new):
    """
    Ricampiona un segnale ECG multicanale dalla frequenza fs_old a fs_new.
    Usa un filtro polifase (resample_poly) per evitare aliasing e distorsioni ai bordi.
    """
    # Se le frequenze sono identiche, restituisci il segnale originale (o una copia)
    if fs_old == fs_new:
        print(f"  -> FS già a {fs_new}Hz. Nessun ricampionamento applicato.")
        return signal.copy()
    
    # Calcola il massimo comune divisore per trovare i fattori interi minimi
    # Esempio: fs_old=400, fs_new=2000 -> gcd=400 -> up=5, down=1
    gcd = math.gcd(fs_old, fs_new)
    up = fs_new // gcd
    down = fs_old // gcd
    
    print(f"  -> Resampling: {fs_old}Hz -> {fs_new}Hz (up={up}, down={down})")
    
    # axis=0 indica che ricampioniamo lungo la dimensione temporale
    resampled_signal = resample_poly(signal, up, down, axis=0)
    
    return resampled_signal

def preprocess_signal(raw, fs_native, target_fs):
    """
    Esegue filtraggio baseline, PLI cancellation, upsampling e cancellazione MECG.
    """
    _, n_channels = raw.shape
    
    # 1. Filtro Baseline
    b_base, a_base = design_baseline_filter(fs_native, fc_high=3.0)
    
    # 2. Loop sui canali per filtraggio base e PLI
    filtered_sigs = np.zeros_like(raw)
    for ch in range(n_channels):
        sig = raw[:, ch]
        sig_no_base = filtfilt(b_base, a_base, sig)
        pli_canceller = AdaptivePLICanceller(fs_native, f_line=50.0)
        filtered_sigs[:, ch] = pli_canceller.apply(sig_no_base)
        
    # 3. Upsampling
    up_factor = int(target_fs / fs_native)
    if up_factor > 1:
        sig_upsampled = resample_poly(filtered_sigs, up_factor, 1, axis=0)
    else:
        sig_upsampled = filtered_sigs
    
    # S4
    S4 = sig_upsampled.copy() 

    mecg_canceller = MECGCanceller(target_fs)
    S5 = mecg_canceller.apply(sig_upsampled)
    
    return S5, S4, up_factor

def load_reference_peaks(record_path, up_factor=1):
    """
    Carica le annotazioni .fqrs (Ground Truth) e le scala in base all'upsampling.
    """
    try:
        ann = wfdb.rdann(record_path, extension='fqrs')
        original_peaks = ann.sample
        upsampled_peaks = original_peaks * up_factor
        
        return upsampled_peaks
    except Exception as e:
        print(e)
        return None

def load_ground_truth(full_path, total_samples, up_factor):
    """
    Carica e allinea i picchi di riferimento
    """
    all_gt_peaks = load_reference_peaks(full_path, up_factor=up_factor)
    segment_gt_aligned = np.array([])
    gt_status = "NO FILE"

    if all_gt_peaks is not None:
        limit_up = total_samples * up_factor
        segment_gt_aligned = all_gt_peaks[all_gt_peaks < limit_up]
        gt_status = f"{len(segment_gt_aligned)} peaks"
        
    return segment_gt_aligned, gt_status

def extract_and_analyze(fetal_matrix, extractor, target_fs):
    """
    Applica PCA per ridurre le dimensioni e rileva i picchi fetali.
    """
    pca = PCA(n_components=1)
    best_signal_fet = pca.fit_transform(fetal_matrix).flatten()

    f_peaks = extractor.detect_fetal_peaks(best_signal_fet, dist_sec=0.25)
    
    return best_signal_fet, f_peaks
    
def calculate_paper_metrics(f_peaks, fs, total_samples):
    """ 
    Calcola le metriche di reliability e BPM medio.
    """
    if len(f_peaks) < 5:
        return 0.0, 0.0, False  # Traccia insufficiente

    rr_sec = np.diff(f_peaks) / fs
    fhr_trace = 60.0 / rr_sec
    
    fhr_times = f_peaks[1:] / fs
    duration_sec = total_samples / fs
    
    n_outliers = 0
    total_points = len(fhr_trace)
    
    for t_start in range(0, int(duration_sec), 10):
        t_end = t_start + 10
        
        mask = (fhr_times >= t_start) & (fhr_times < t_end)
        block_fhr = fhr_trace[mask]
        
        if len(block_fhr) == 0:
            continue
        
        block_median = np.median(block_fhr)
        
        outliers_in_block = np.sum(np.abs(block_fhr - block_median) > 10)
        n_outliers += outliers_in_block

    reliability = 1.0 - (n_outliers / total_points)
    
    mean_bpm = np.mean(fhr_trace)
    
    is_successful = (60 <= mean_bpm <= 220) and (reliability > 0.60)
    
    return reliability, mean_bpm, is_successful

def calculate_snr_sir(s4, s5, s6):
    """
    Implementazione esatta Eq (4) e (5)
    """
    # Evita divisioni per zero
    eps = 1e-10
    
    # Calcolo Potenza del segnale (Mean Squared)
    def power(x):
        return np.mean(x**2, axis=0)
    
    P_F = power(s6)
    mecg_est = s4 - s5
    P_M = power(mecg_est)
    noise_est = s5 - s6
    P_N = power(noise_est)
    
    snr_linear = P_F / (P_N + eps)
    sir_linear = P_F / (P_M + eps)
    
    snr_db = 10 * np.log10(snr_linear)
    sir_db = 10 * np.log10(sir_linear)
    
    return np.mean(snr_db), np.mean(sir_db)

def plot_results(rec_name, signal, detected_peaks, gt_peaks, extractor, target_fs):
    
    plt.figure(figsize=(15, 10))
    ax1 = plt.subplot(2, 1, 1)
    start_s, duration_s = 20, 10
    s_idx = start_s * target_fs
    e_idx = (start_s + duration_s) * target_fs
    if e_idx > len(signal): e_idx = len(signal)

    t_zoom = np.arange(s_idx, e_idx) / target_fs
    sig_zoom = signal[s_idx:e_idx]
    
    plt.plot(t_zoom, sig_zoom, color='green', alpha=0.5, label='Segnale Residuo PCA')
    
    p_alg = detected_peaks[(detected_peaks >= s_idx) & (detected_peaks < e_idx)]
    if len(p_alg) > 0:
        plt.scatter(p_alg/target_fs, signal[p_alg], color='red', marker='x', s=100, lw=2, label='Rilevati', zorder=5)
    
    if len(gt_peaks) > 0:
        gt_zoom = gt_peaks[(gt_peaks >= s_idx) & (gt_peaks < e_idx)]
        if len(gt_zoom) > 0:
            plt.scatter(gt_zoom/target_fs, signal[gt_zoom], facecolors='none', edgecolors='blue', marker='o', s=150, lw=2, label='GT')

    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)

    ax2 = plt.subplot(2, 1, 2)
    avg_beat, _ = extractor.compute_average_beat(signal, detected_peaks)
    
    t_ms = np.linspace(0, 250, len(avg_beat))
    
    plt.plot(t_ms, avg_beat, color='blue', lw=2)
    plt.title("Morfologia Fetale Media")
    plt.xlabel("Tempo (ms)")
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
def debug_plot_overlay(raw, fs_raw, processed, fs_proc, title="Sovrapposizione Raw vs Processed", start_s=0, duration_s=10, channel_idx=0):

    plt.figure(figsize=(12, 6))
    s_idx_raw = int(start_s * fs_raw)
    e_idx_raw = int((start_s + duration_s) * fs_raw)
    e_idx_raw = min(e_idx_raw, len(raw))
    
    # Processed
    s_idx_proc = int(start_s * fs_proc)
    e_idx_proc = int((start_s + duration_s) * fs_proc)
    e_idx_proc = min(e_idx_proc, len(processed))
    
    # Necessari perché fs è diversa
    t_raw = np.arange(s_idx_raw, e_idx_raw) / fs_raw
    t_proc = np.arange(s_idx_proc, e_idx_proc) / fs_proc

    sig_raw = raw[s_idx_raw:e_idx_raw, channel_idx] / 10
    sig_proc = processed[s_idx_proc:e_idx_proc, channel_idx] / 10

    plt.plot(t_raw, sig_raw, color='grey', alpha=0.6, label=f'Raw', linewidth=1.5)
    plt.plot(t_proc, sig_proc, color='#007acc', label=f'MECG Removed', linewidth=2)
    
    plt.title(title)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude (uV)")
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

def get_average_beat_per_channel(signal, peaks, window_sec=0.25, fs=2000):

    n_samples, n_channels = signal.shape
    win_samples = int(window_sec * fs)
    half_win = win_samples // 2
    
    avg_beats = np.zeros((win_samples, n_channels))
    
    for ch in range(n_channels):
        beats = []
        for p in peaks:
            start = p - half_win
            end = p + half_win
            if start >= 0 and end < n_samples:
                beats.append(signal[start:end, ch])
        
        if len(beats) > 0:
            avg_beats[:, ch] = np.mean(beats, axis=0)
            
    return avg_beats

def construct_s6_signal(s5_signal, peaks, fs=2000):

    n_samples, n_channels = s5_signal.shape
    s6 = np.zeros_like(s5_signal)
    
    window_sec = 0.25
    win_samples = int(window_sec * fs)
    half_win = win_samples // 2
    
    avg_templates = get_average_beat_per_channel(s5_signal, peaks, window_sec, fs)
    
    for ch in range(n_channels):
        template = avg_templates[:, ch]
        for p in peaks:
            start = p - half_win
            end = p + half_win
            if start >= 0 and end < n_samples:
                s6[start:end, ch] += template
                
    return s6

