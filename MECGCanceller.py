import numpy as np
from scipy.signal import butter, filtfilt, hilbert, find_peaks, correlate
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

class MECGCanceller:

    def __init__(self, fs):
        self.fs = fs
        
        # Finestre
        self.win_pre = 0.2        # secondi prima del picco
        self.win_post = 0.4       # secondi dopo il picco
        self.n_avg = 10           # template da mediare
        self.qrs_win = 0.04       # finestra centrale stretta per regressione

        # Regolarizzazione LMS (Ridge)
        self.lam = 1e-3

    def _bandpass_detect(self, X):
        b, a = butter(3, [5/(self.fs/2), 40/(self.fs/2)], btype='band')
        return filtfilt(b, a, X, axis=0)

    def _detect_maternal_qrs(self, signal_matrix):
        
        # --- FASE 1: Rilevamento Robusto (IL TUO CODICE) ---
        # Questa parte è ottima per non perdere nessun battito, anche se rumoroso.
        
        bp = self._bandpass_detect(signal_matrix)
        centered = bp - np.mean(bp, axis=0)
        pca_signal = PCA(n_components=1).fit_transform(centered).flatten()

        # Importante: Hilbert distrugge la polarità. 
        # Salviamo il segnale PCA con polarità corretta per dopo.
        if np.max(pca_signal) < np.abs(np.min(pca_signal)):
            pca_signal = -pca_signal  # Rendo il picco R positivo per comodità

        # Calcolo Envelope per detection (come facevi tu)
        env = np.abs(hilbert(pca_signal))
        win = int(0.02 * self.fs)
        env_smooth = np.convolve(env, np.ones(win)/win, mode='same')
        
        thresh = np.mean(env_smooth) + 1.5*np.std(env_smooth)
        min_dist = int(0.35 * self.fs)
        
        # Questi sono i picchi "grezzi" (affetti da jitter dovuto allo smoothing)
        rough_peaks, _ = find_peaks(env_smooth, height=thresh, distance=min_dist)

        # --- FASE 2: Raffinamento alla Abboud & Beker (ALIGNMENT) ---
        # Usiamo i picchi grezzi solo come indizi per trovare il VERO centro nel segnale raw.
        
        # 1. Costruiamo il Template Reale (Mediana dei battiti estratti dal segnale PCA pulito)
        window_sec = 0.08 # 80ms bastano per il nucleo del QRS
        w_half = int((window_sec * self.fs) / 2)
        
        segments = []
        valid_rough = []
        
        for p in rough_peaks:
            if p - w_half < 0 or p + w_half >= len(pca_signal): continue
            segments.append(pca_signal[p - w_half : p + w_half])
            valid_rough.append(p)
            
        if not segments: return rough_peaks, pca_signal

        # Template MEDIO del paziente
        maternal_template = np.median(np.array(segments), axis=0)
        
        # 2. Cross-Correlazione per correggere il Jitter
        final_peaks = []
        search_window = 10 # Cerchiamo +/- 10 samples intorno al picco grezzo
        
        for p in valid_rough:
            start = p - w_half - search_window
            end = p + w_half + search_window
            
            if start < 0 or end >= len(pca_signal): continue
            
            # Segmento locale raw
            segment = pca_signal[start:end]
            
            # Correlazione con il template
            cc = correlate(segment, maternal_template, mode='valid')
            
            # Il picco della correlazione è il punto di allineamento perfetto
            shift_idx = np.argmax(cc)
            
            # Calcolo la posizione assoluta corretta
            # start + shift_idx ci porta all'inizio del match.
            # Aggiungiamo w_half per trovare il centro del QRS.
            true_peak = start + shift_idx + w_half
            
            final_peaks.append(true_peak)

        return np.array(final_peaks), pca_signal

    def _get_segments(self, win_len):

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
        """Estrae un segmento di segnale centrato (o quasi) sul picco."""
        start = peak - int(self.win_pre * self.fs)
        end = start + win_len
        
        # Controllo bordi
        if start < 0 or end > n_samples:
            return None
        
        return sig[start:end]

    def apply(self, signal_matrix):
        
        # Dimensioni originali
        n_samples_orig, n_channels = signal_matrix.shape
        
        # 1) Calcoliamo quanto padding serve
        # Aggiungiamo un po' di margine oltre alle finestre definite
        pad_pre = int(self.win_pre * self.fs) + 100
        pad_post = int(self.win_post * self.fs) + 100
        
        # 2) Eseguiamo il padding del segnale (modalità 'edge' ripete l'ultimo valore)
        # axis=((pad_pre, pad_post), (0,0)) significa: pad temporale, nessun pad sui canali
        sig_padded = np.pad(signal_matrix, ((pad_pre, pad_post), (0, 0)), mode='edge')
        
        n_samples_pad = sig_padded.shape[0]
        fetal_ecg_padded = np.zeros_like(sig_padded)

        # 3) Rileviamo i picchi SUL SEGNALE PADDATO
        # (Così i picchi originali saranno spostati in avanti di 'pad_pre')
        m_peaks, pca_sig = self._detect_maternal_qrs(sig_padded)

        # 4) Cancellazione canale-per-canale (logica identica a prima)
        win_len = int((self.win_pre + self.win_post) * self.fs)
        
        for ch in range(n_channels):
            sig = sig_padded[:, ch]
            residue = sig.copy()
            
            # --- SEEDING (Inizializzazione buffer) ---
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
                if curr is None: continue # Non dovrebbe succedere grazie al padding
                
                if len(templates) > 0:
                    avg_template = np.mean(templates, axis=0)
                else:
                    avg_template = curr 

                templates.append(curr)
                if len(templates) > self.n_avg: templates.pop(0)

                # LMS
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

                # Sottrazione (attenzione agli indici padded)
                start = peak - int(self.win_pre * self.fs)
                end = start + win_len
                residue[start:end] = curr - fitted

            fetal_ecg_padded[:, ch] = residue

        fetal_ecg_final = fetal_ecg_padded[pad_pre:-pad_post, :]
        
        return fetal_ecg_final