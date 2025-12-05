import numpy as np
from scipy.signal import butter, filtfilt, hilbert, find_peaks, correlate
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

class MECGCanceller:

    def __init__(self, fs):
        self.fs = fs
        
        self.win_pre = 0.2
        self.win_post = 0.4
        self.n_avg = 10
        self.qrs_win = 0.04
        self.lam = 1e-3

    def _bandpass_detect(self, X):
        b, a = butter(3, [5/(self.fs/2), 40/(self.fs/2)], btype='band')
        return filtfilt(b, a, X, axis=0)

    def _detect_maternal_qrs(self, signal_matrix):
        """
        Rileva i complessi QRS materni utilizzando un approccio a due stadi per massima precisione temporale.

        Fase 1:
        - Fusione dei canali tramite PCA per isolare la componente materna dominante.
        - Calcolo dell'inviluppo di Hilbert per enfatizzare i picchi di energia.
        - Sogliatura dinamica per trovare i candidati R-peak iniziali.

        Fase 2:
        - Creazione di un template mediano dai battiti rilevati.
        - Allineamento fine di ogni battito tramite cross-correlazione col template (metodo Abboud & Beker).
        """
        
        bp = self._bandpass_detect(signal_matrix)
        centered = bp - np.mean(bp, axis=0)
        pca_signal = PCA(n_components=1).fit_transform(centered).flatten()

        if np.max(pca_signal) < np.abs(np.min(pca_signal)):
            pca_signal = -pca_signal

        env = np.abs(hilbert(pca_signal))
        win = int(0.02 * self.fs)
        env_smooth = np.convolve(env, np.ones(win)/win, mode='same')
        
        thresh = np.mean(env_smooth) + 1.5*np.std(env_smooth)
        min_dist = int(0.35 * self.fs)
        
        rough_peaks, _ = find_peaks(env_smooth, height=thresh, distance=min_dist)

        window_sec = 0.08
        w_half = int((window_sec * self.fs) / 2)
        
        segments = []
        valid_rough = []
        
        for p in rough_peaks:
            if p - w_half < 0 or p + w_half >= len(pca_signal): continue
            segments.append(pca_signal[p - w_half : p + w_half])
            valid_rough.append(p)
            
        if not segments: return rough_peaks, pca_signal

        maternal_template = np.median(np.array(segments), axis=0)
        
        final_peaks = []
        search_window = 10
        
        for p in valid_rough:
            
            start = p - w_half - search_window
            end = p + w_half + search_window
            
            if start < 0 or end >= len(pca_signal): continue
            
            segment = pca_signal[start:end]
            cc = correlate(segment, maternal_template, mode='valid')
            shift_idx = np.argmax(cc)
            true_peak = start + shift_idx + w_half
            
            final_peaks.append(true_peak)

        return np.array(final_peaks), pca_signal

    def _get_segments(self, win_len):
        """
        Suddivide la finestra del battito in tre regioni morfologiche: onda P, QRS e onda T
        """
        
        center = int(self.win_pre * self.fs)
        qrs_rad = int(self.qrs_win * self.fs)

        qrs_start = max(0, center - qrs_rad)
        qrs_end = min(win_len, center + qrs_rad)

        p_start = 0
        p_end = qrs_start

        t_start = qrs_end
        t_end = win_len

        return (p_start, p_end), (qrs_start, qrs_end), (t_start, t_end)

    def _extract_beat(self, sig, peak, win_len, n_samples):
        """
        Estrae un segmento di segnale centrato (o quasi) sul picco.
        """
        
        start = peak - int(self.win_pre * self.fs)
        end = start + win_len
        
        if start < 0 or end > n_samples:
            return None
        
        return sig[start:end]

    def apply(self, signal_matrix):
        
        """
        Applica l'algoritmo di cancellazione dell'ECG materno (MECG) tramite sottrazione di template adattivo.

        Il metodo esegue le seguenti operazioni:
        1. Padding del segnale per gestire i bordi.
        2. Rilevamento dei complessi QRS materni.
        3. Per ogni canale, costruisce un template medio dinamico (running average) del battito materno.
        4. Adatta il template al battito corrente (scaling dei segmenti P, QRS, T) utilizzando i minimi quadrati regolarizzati (Ridge Regression).
        5. Sottrae la stima materna per isolare il residuo fetale.
        """
        
        _, n_channels = signal_matrix.shape
        
        pad_pre = int(self.win_pre * self.fs) + 100
        pad_post = int(self.win_post * self.fs) + 100
        
        sig_padded = np.pad(signal_matrix, ((pad_pre, pad_post), (0, 0)), mode='edge')
        
        n_samples_pad = sig_padded.shape[0]
        fetal_ecg_padded = np.zeros_like(sig_padded)

        m_peaks, _ = self._detect_maternal_qrs(sig_padded)

        win_len = int((self.win_pre + self.win_post) * self.fs)
        
        for ch in range(n_channels):
            sig = sig_padded[:, ch]
            residue = sig.copy()
            
            init_beats = []
            for peak in m_peaks:
                beat = self._extract_beat(sig, peak, win_len, n_samples_pad)
                if beat is not None:
                    init_beats.append(beat)
                if len(init_beats) >= self.n_avg: break
            
            templates = []
            if len(init_beats) > 0:
                initial_avg = np.mean(init_beats, axis=0)
                templates = [initial_avg for _ in range(self.n_avg)]
            
            for peak in m_peaks:
                curr = self._extract_beat(sig, peak, win_len, n_samples_pad)
                if curr is None: continue
                
                if len(templates) > 0:
                    avg_template = np.mean(templates, axis=0)
                else:
                    avg_template = curr 

                templates.append(curr)
                if len(templates) > self.n_avg: templates.pop(0)

                M = np.zeros((win_len, 3))
                (p_s, p_e), (q_s, q_e), (t_s, t_e) = self._get_segments(win_len)
                M[p_s:p_e, 0] = avg_template[p_s:p_e]
                M[q_s:q_e, 1] = avg_template[q_s:q_e]
                M[t_s:t_e, 2] = avg_template[t_s:t_e]

                try:
                    A = M.T @ M + self.lam * np.eye(3)
                    b = M.T @ curr
                    coeffs = np.linalg.solve(A, b)
                    fitted = M @ coeffs
                except:
                    fitted = 0

                start = peak - int(self.win_pre * self.fs)
                end = start + win_len
                residue[start:end] = curr - fitted

            fetal_ecg_padded[:, ch] = residue

        fetal_ecg_final = fetal_ecg_padded[pad_pre:-pad_post, :]
        
        return fetal_ecg_final