import numpy as np
from scipy.signal import find_peaks
from sklearn.decomposition import PCA

class FetalECGExtractor:
    def __init__(self, fs):
        self.fs = fs

    def detect_fetal_peaks(self, residual_signal, dist_sec=0.25, prominence_factor=0.3):
        """
        Individua i complessi QRS fetali nel segnale residuo.
        """
        
        energy_signal = residual_signal ** 3
        
        # 2. Conversione parametri temporali in campioni
        min_dist = int(dist_sec * self.fs)
        
        # 3. Calcolo soglia dinamica
        signal_level = np.percentile(energy_signal, 98)
        prominence = signal_level * prominence_factor

        # 4. Rilevamento picchi con Scipy
        peaks, _ = find_peaks(energy_signal, distance=min_dist, prominence=prominence)
        
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
            
            # Detrending locale: 
            # Se è multicanale (axis=0), sottrae la media di ogni colonna indipendentemente
            segment = segment - np.mean(segment, axis=0)
            
            beats.append(segment)
            
        if len(beats) == 0:
            # Restituisce array di zeri con la forma corretta
            if n_channels == 1:
                return np.zeros(w_samples), 0
            else:
                return np.zeros((w_samples, n_channels)), 0
            
        beats_array = np.array(beats) 
        # beats_array shape: (Num_Battiti, Campioni_Finestra, Canali)
        
        # Media lungo l'asse dei battiti (axis=0)
        averaged_beat = np.mean(beats_array, axis=0)
        
        return averaged_beat, len(beats)