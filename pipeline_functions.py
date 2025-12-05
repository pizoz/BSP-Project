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
    
    if fs_old == fs_new:
        return signal.copy()
    
    gcd = math.gcd(fs_old, fs_new)
    up = fs_new // gcd
    down = fs_old // gcd
    
    resampled_signal = resample_poly(signal, up, down, axis=0)
    
    return resampled_signal

def preprocess_signal(raw, fs_native, target_fs):
    """
    Esegue filtraggio baseline, PLI cancellation, upsampling e cancellazione MECG.
    """
    
    _, n_channels = raw.shape
    
    b_base, a_base = design_baseline_filter(fs_native, fc_high=3.0)
    
    filtered_sigs = np.zeros_like(raw)
    for ch in range(n_channels):
        sig = raw[:, ch]
        sig_no_base = filtfilt(b_base, a_base, sig)
        pli_canceller = AdaptivePLICanceller(fs_native, f_line=50.0)
        filtered_sigs[:, ch] = pli_canceller.apply(sig_no_base)
        
    up_factor = int(target_fs / fs_native)
    if up_factor > 1:
        sig_upsampled = resample_poly(filtered_sigs, up_factor, 1, axis=0)
    else:
        sig_upsampled = filtered_sigs
    
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
    
    pca = PCA(n_components=1)
    best_signal_fet = pca.fit_transform(fetal_matrix).flatten()

    if np.abs(np.min(best_signal_fet)) > np.abs(np.max(best_signal_fet)):
        best_signal_fet = -best_signal_fet

    f_peaks = extractor.detect_fetal_peaks(best_signal_fet, dist_sec=0.25)
    average_beat, _  = extractor.compute_average_beat(best_signal_fet, f_peaks, window_sec=0.125)
    
    return best_signal_fet, average_beat, f_peaks

def extract_and_analyze_multichannel(fetal_matrix, extractor, target_fs):

    pca = PCA(n_components=1)
    detection_signal = pca.fit_transform(fetal_matrix).flatten()

    if np.abs(np.min(detection_signal)) > np.abs(np.max(detection_signal)):
        detection_signal = -detection_signal

    f_peaks = extractor.detect_fetal_peaks(detection_signal, dist_sec=0.25)

    average_beat_matrix, _ = extractor.compute_average_beat(fetal_matrix, f_peaks, window_sec=0.125)
    
    return fetal_matrix, average_beat_matrix, f_peaks

def evaluate_performance(ref_peaks, det_peaks, fs, tol_ms=50):
    """
    Confronta i picchi rilevati con quelli di riferimento.
    Calcola TP (True Positive), FP (False Positive), FN (False Negative).
    """
    tolerance_samples = int((tol_ms / 1000.0) * fs)
    
    tp = 0
    fp = 0
    fn = 0
    
    det_peaks_list = list(det_peaks)
    
    for ref in ref_peaks:
        matches = [d for d in det_peaks_list if abs(d - ref) <= tolerance_samples]
        
        if matches:
            tp += 1
            det_peaks_list.remove(matches[0]) 
        else:
            fn += 1
            
    fp = len(det_peaks_list)
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
    accuracy = (tp) / (tp + fp + fn) if (tp + fp + fn) > 0 else 0
    
    
    return {
        "Accuracy": accuracy * 100,
        "Se": sensitivity * 100,
        "PPV": precision * 100,
        "F1": f1_score * 100
    }

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

def calculate_snr_sir(s4, s5, s6, noise):
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
    P_N = power(noise)
    
    snr_linear = P_F / (P_N + eps)
    sir_linear = P_F / (P_M + eps)
    
    snr_db = 10 * np.log10(snr_linear)
    sir_db = 10 * np.log10(sir_linear)
    
    return np.mean(snr_db), np.mean(sir_db)

def subtract_average_beat(s5_signal, peaks, average_beat_matrix):
    """
    Sottrae il template medio dal segnale S5 per ottenere il rumore.
    """
    
    win_samples, _ = average_beat_matrix.shape
    half_win = win_samples // 2
    
    n_samples_total = s5_signal.shape[0]

    noise_residual = s5_signal.copy()
    
    for p in peaks:
        start = p - half_win
        end = start + win_samples 
        
        if start < 0:
            offset = -start
            if end > 0:
                 noise_residual[0:end, :] -= average_beat_matrix[offset:, :]
                 
        elif end > n_samples_total:
            excess = end - n_samples_total
            valid_len = win_samples - excess
            if valid_len > 0:
                noise_residual[start:n_samples_total, :] -= average_beat_matrix[:valid_len, :]
                
        else:
            noise_residual[start:end, :] -= average_beat_matrix

    return noise_residual

"""
    Plot functions
"""
def plot_results(rec_name, signal, detected_peaks, gt_peaks, extractor, target_fs):
    
    plt.figure(figsize=(15, 12))
    
    ax1 = plt.subplot(2, 1, 1)
    
    start_s, duration_s = 20, 10
    s_idx = int(start_s * target_fs)
    e_idx = int((start_s + duration_s) * target_fs)
    
    if e_idx > len(signal): 
        e_idx = len(signal)
        s_idx = max(0, e_idx - int(duration_s * target_fs))

    t_zoom = np.arange(s_idx, e_idx) / target_fs
    sig_zoom = signal[s_idx:e_idx]

    if signal.ndim > 1:
        n_channels = signal.shape[1]
        colors = plt.cm.viridis(np.linspace(0, 1, n_channels))
        for ch in range(n_channels):
            plt.plot(t_zoom, sig_zoom[:, ch], color=colors[ch], alpha=0.6, label=f'Ch {ch+1}')
        
        ref_full_signal = signal[:, 0]
        plt.title(f"Segnale Multicanale (Zoom {duration_s}s) - {rec_name}")
    else:
        plt.plot(t_zoom, sig_zoom, color='green', alpha=0.8, label='Segnale Residuo')
        ref_full_signal = signal
        plt.title(f"Segnale Residuo (Zoom {duration_s}s) - {rec_name}")

    p_alg = detected_peaks[(detected_peaks >= s_idx) & (detected_peaks < e_idx)]
    if len(p_alg) > 0:
        plt.scatter(p_alg/target_fs, ref_full_signal[p_alg], 
                   color='red', marker='x', s=100, lw=2, label='Rilevati (Alg)', zorder=5)
    
    if len(gt_peaks) > 0:
        gt_zoom = gt_peaks[(gt_peaks >= s_idx) & (gt_peaks < e_idx)]
        if len(gt_zoom) > 0:
            plt.scatter(gt_zoom/target_fs, ref_full_signal[gt_zoom], 
                       facecolors='none', edgecolors='blue', marker='o', s=150, lw=2, label='GT')

    plt.legend(loc='upper right', fontsize='small', ncol=2)
    plt.grid(True, alpha=0.3)
    plt.ylabel("Ampiezza")

    ax2 = plt.subplot(2, 1, 2)
    
    avg_beat, num_beats = extractor.compute_average_beat(signal, detected_peaks)
    
    t_ms = np.linspace(0, 250, len(avg_beat))
    
    if avg_beat.ndim > 1:
        n_channels = avg_beat.shape[1]
        colors = plt.cm.viridis(np.linspace(0, 1, n_channels))
        
        for ch in range(n_channels):
            plt.plot(t_ms, avg_beat[:, ch], color=colors[ch], lw=2, alpha=0.8, label=f'Ch {ch+1}')
        
        plt.legend(loc='upper right')
        plt.title(f"Morfologia Fetale Media (Multicanale) - Basata su {num_beats} battiti")
    else:
        # Plot Monocanale
        plt.plot(t_ms, avg_beat, color='blue', lw=2)
        plt.title(f"Morfologia Fetale Media - Basata su {num_beats} battiti")

    plt.xlabel("Tempo (ms)")
    plt.ylabel("Ampiezza Media")
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
def debug_plot_overlay(raw, fs_raw, processed, fs_proc, title="Sovrapposizione Raw vs Processed", start_s=0, duration_s=10, channel_idx=0):

    plt.figure(figsize=(12, 6))
    s_idx_raw = int(start_s * fs_raw)
    e_idx_raw = int((start_s + duration_s) * fs_raw)
    e_idx_raw = min(e_idx_raw, len(raw))
    
    s_idx_proc = int(start_s * fs_proc)
    e_idx_proc = int((start_s + duration_s) * fs_proc)
    e_idx_proc = min(e_idx_proc, len(processed))
    
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

def plot_pipeline_stages(s1_raw, s3_cleaned, s5_cancelled, s6_avg_beat, fs, start_sec=10, duration_sec=4):

    color_signal = "#A50000"
    color_beat = '#FDB813'
    line_width = 1.2
    
    idx_start = int(start_sec * fs)
    idx_end = int((start_sec + duration_sec) * fs)
    
    if idx_end > len(s1_raw):
        idx_end = len(s1_raw)
        idx_start = max(0, idx_end - int(duration_sec * fs))
    
    _, axes = plt.subplots(nrows=4, ncols=4, figsize=(18, 10), constrained_layout=True)
    
    cols_titles = ["S1 raw abdominal", "S3 cleaned", "S5 MECG cancelled", "S6 average fetal complex"]
    
    for ch in range(4):
        
        d1 = s1_raw[idx_start:idx_end, ch]
        d3 = s3_cleaned[idx_start:idx_end, ch]
        d5 = s5_cancelled[idx_start:idx_end, ch]
        d6 = s6_avg_beat[:, ch]
        
        all_data_ch = np.concatenate([d1, d3, d5]) 
        y_min = np.min(all_data_ch)
        y_max = np.max(all_data_ch)
        
        margin = (y_max - y_min) * 0.1
        y_lims = (y_min - margin, y_max + margin)

        ax = axes[ch, 0]
        ax.plot(d1, color=color_signal, lw=line_width)
        ax.set_ylabel(f"AECG {ch+1}", fontsize=12, fontweight='bold')
        ax.set_ylim(y_lims)
        
        ax = axes[ch, 1]
        ax.plot(d3, color=color_signal, lw=line_width)
        ax.set_ylim(y_lims)
        
        ax = axes[ch, 2]
        ax.plot(d5, color="#C77415", lw=line_width)
        ax.set_ylim(y_lims)
        
        ax = axes[ch, 3]
        ax.plot(d6, color=color_beat, lw=2)
        
        for col in range(4):
            ax_curr = axes[ch, col]
            ax_curr.spines['top'].set_visible(False)
            ax_curr.spines['right'].set_visible(False)
            ax_curr.set_xticks([])
            ax_curr.set_yticks([]) # Rimuove i numeri per pulizia
            
            if ch == 0:
                ax_curr.set_title(cols_titles[col], fontsize=14, pad=15)

    plt.show()

def plot_reliability_correlations(df, x_fit, y_fit, fit_label):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))

    ax1.scatter(df['SNR'], df['Reliability_Pct'], marker='x', color='black', label='Dati')
    
    if len(x_fit) > 0:
        ax1.plot(x_fit, y_fit, color='red', linewidth=2, label=fit_label)
        ax1.legend(loc='lower right', fontsize='small')

    ax1.set_xlabel('SNR [dB]', fontsize=11)
    ax1.set_ylabel('SA reliability [%]', fontsize=11)
    ax1.tick_params(direction='in')
    ax1.grid(True, alpha=0.3)

    ax2.scatter(df['SIR'], df['Reliability_Pct'], marker='x', color='black')
    ax2.set_xlabel('SIR [dB]', fontsize=11)
    ax2.set_ylabel('SA reliability [%]', fontsize=11)
    ax2.tick_params(direction='in')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()