import numpy as np
from scipy.signal import find_peaks

class FetalECGExtractor:
    """
    Classe per il rilevamento e l'estrazione dell'ECG fetale da segnali residui
    (dopo la cancellazione dell'ECG materno).
    """
    def __init__(self, fs):
        self.fs = fs

    def detect_fetal_peaks(self, residual_signal, dist_sec=0.25, prominence_factor=0.3):
        
        # 1. Enfatizzazione dell'energia locale
        energy_signal = residual_signal ** 3
        
        # 2. Conversione parametri temporali in campioni
        min_dist = int(dist_sec * self.fs)
        
        signal_level = np.percentile(energy_signal, 98)
        prominence = signal_level * prominence_factor

        # 4. Rilevamento picchi con Scipy
        peaks, _ = find_peaks(energy_signal, distance=min_dist, prominence=prominence)
        
        return peaks

    def compute_average_beat(self, residual_signal, peaks, window_sec=0.125):
        
        w_samples = int(window_sec * self.fs)
        half_w = w_samples // 2
        
        beats = []
        
        for p in peaks:
            start = p - half_w
            end = p + half_w
            
            # Scarta battiti ai bordi del segnale
            if start < 0 or end > len(residual_signal):
                continue
                
            segment = residual_signal[start:end]
            
            # Detrending locale: sottrae la media per allineare lo zero
            segment = segment - np.mean(segment)
            
            beats.append(segment)
            
        if len(beats) == 0:
            return np.zeros(w_samples), 0
            
        beats_matrix = np.array(beats)
        averaged_beat = np.mean(beats_matrix, axis=0)
        
        return averaged_beat, len(beats)