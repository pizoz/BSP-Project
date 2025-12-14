import numpy as np
from scipy.signal import find_peaks, hilbert
from scipy.stats import skew
from sklearn.decomposition import PCA


class FetalECGExtractor:
    def __init__(self, fs):
        self.fs = fs

    def detect_fetal_peaks(self, residual_signal, dist_sec=0.25, prominence_factor=0.3):
        """
        Individua i complessi QRS fetali nel segnale residuo.
        """
        
        energy_signal = residual_signal ** 3
        
        min_dist = int(dist_sec * self.fs)
        
        signal_level = np.percentile(energy_signal, 98)
        prominence = signal_level * prominence_factor

        peaks, _ = find_peaks(energy_signal, distance=min_dist, prominence=prominence)
        
        return peaks
    
    def detect_fetal_peaks_2(self, residual_matrix):
        """
        Robust Fetal QRS detection using PCA and Hilbert envelope.
        """
        if residual_matrix.ndim == 1:
            residual_matrix = residual_matrix.reshape(-1, 1)
        
        centered = residual_matrix - np.mean(residual_matrix, axis=0)
        
        pca = PCA(n_components=1)
        fetal_comp = pca.fit_transform(centered).flatten()

        if skew(fetal_comp) < 0:
            fetal_comp = -fetal_comp

        env = np.abs(hilbert(fetal_comp))
        
        win_sec = 0.015
        win_samples = int(win_sec * self.fs)
        env_smooth = np.convolve(env, np.ones(win_samples)/win_samples, mode='same')
        
        thresh = np.mean(env_smooth) + 1.0 * np.std(env_smooth)
        

        min_dist = int(0.25 * self.fs) 
            
        peaks, _ = find_peaks(env_smooth, height=thresh, distance=min_dist)
            
        return peaks

    def compute_average_beat(self, signal, peaks, window_sec=0.125):
        """
        Calcola la morfologia media del battito cardiaco sincronizzando 
        il segnale sugli indici dei picchi forniti
        """
        
        w_samples = int(window_sec * self.fs)
        half_w = w_samples // 2
        
        if signal.ndim == 1:
            n_channels = 1
        else:
            n_channels = signal.shape[1]
        
        beats = []
        
        for p in peaks:
            start = p - half_w
            end = p + half_w
            
            if start < 0 or end > len(signal):
                continue
            segment = signal[start:end]
            segment = segment - np.mean(segment, axis=0)
            
            beats.append(segment)
            
        if len(beats) == 0:
            if n_channels == 1:
                return np.zeros(w_samples), 0
            else:
                return np.zeros((w_samples, n_channels)), 0
            
        beats_array = np.array(beats) 
        averaged_beat = np.mean(beats_array, axis=0)
        
        return averaged_beat, len(beats)