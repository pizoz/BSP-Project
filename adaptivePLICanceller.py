import numpy as np
import pandas as pd
from scipy.signal import butter

def lfilter_step(b, a, x, zi):
    
    """Esegue un singolo passo di filtro IIR mantenendo lo stato zi."""

    y = b[0]*x + zi[0]
    for i in range(len(zi)-1):
        zi[i] = b[i+1]*x - a[i+1]*y + zi[i+1]
    zi[-1] = b[-1]*x - a[-1]*y
    return y, zi

def freqz_scalar(b, a, f, fs):
    
    """Calcola risposta in frequenza a una singola frequenza f."""
    
    w = 2 * np.pi * f / fs
    zm1 = np.exp(-1j * w)
    num = np.polyval(b, zm1)
    den = np.polyval(a, zm1)
    return w, num/den

class AdaptivePLICanceller:
    def __init__(self, fs, f_line=50.0):
        self.fs = fs
        self.f_line = f_line
        self.w_n = 2 * np.pi * f_line / fs
        
        tau = 0.13 
        self.K_a = 1.0 / (fs * tau)
        
        zeta = 1.0
        ratio_wn_wp = 0.04
        omega_n = self.w_n * ratio_wn_wp
        
        self.K_dw = omega_n ** 2
        self.K_phi = 2 * zeta * omega_n
        
        cutoff_hz = 80.0
        nyq = 0.5 * fs
        b, a = butter(2, cutoff_hz / nyq, btype='high')
        
        w, h = freqz_scalar(b, a, f_line, fs)
        gain_at_50 = np.abs(h)
        self.b_err = b / gain_at_50
        self.a_err = a 
        
    def _apply_comb_filter_and_detect_blocking(self, d_k):
        """
        Sezione II-E: Adaptation Blocking
        Usa un filtro a pettine (Comb) per stimare l'energia del segnale (QRS)
        senza l'interferenza di rete, per decidere quando bloccare l'adattamento.
        """

        beta = int(round(self.fs / self.f_line))
        
        pad = np.zeros(beta)
        d_shifted = np.concatenate((pad, d_k[:-beta]))
        d_H = d_k - d_shifted
        
        win_len = int(self.fs)
        d_H_series = pd.Series(d_H)
        sigma = d_H_series.rolling(window=win_len, center=True).std().fillna(0).values
        
        chi = np.sqrt(2) * sigma
        
        raw_mask = np.abs(d_H) > chi
        expansion = int(0.05 * self.fs)
        blocking_mask = np.convolve(raw_mask.astype(int), np.ones(2*expansion+1), mode='same') > 0
        
        return blocking_mask

    def apply(self, signal):
        """
        Esegue il loop di cancellazione adattiva (Sample-by-sample).
        Implementa lo schema di Figura 4 e le Equazioni 30-32.
        """
        n_samples = len(signal)
        e_output = np.zeros(n_samples)
        
        theta_a = 0.0
        theta_phi = 0.0
        theta_dw = 0.0
        
        zi_e = np.zeros(max(len(self.a_err), len(self.b_err)) - 1)
        zi_y_sin = np.zeros_like(zi_e)
        zi_y_cos = np.zeros_like(zi_e)
        
        blocking_mask = self._apply_comb_filter_and_detect_blocking(signal)
        
        for k in range(n_samples):

            arg = self.w_n * k + theta_phi
            ref_sin = np.sin(arg)
            ref_cos = np.cos(arg)
            
            x_hat = theta_a * ref_sin
            d_val = signal[k]
            e_val = d_val - x_hat
            e_output[k] = e_val
            
            e_w, zi_e = lfilter_step(self.b_err, self.a_err, e_val, zi_e)
            
            y_sin_w, zi_y_sin = lfilter_step(self.b_err, self.a_err, ref_sin, zi_y_sin)
            y_cos_w, zi_y_cos = lfilter_step(self.b_err, self.a_err, ref_cos, zi_y_cos)
            
            if not blocking_mask[k]:

                alpha = 1.0 / theta_a if theta_a > 1e-6 else 1.0
                
                eta_a = e_w * y_sin_w
                eta_phi = e_w * (alpha * y_cos_w)
                
                theta_a_next = theta_a + self.K_a * eta_a
                
                if theta_a_next < 0: theta_a_next = 0
                
                theta_dw_next = theta_dw + self.K_dw * eta_phi
                
                max_dw = 2 * np.pi * 4.0 / self.fs
                theta_dw_next = np.clip(theta_dw_next, -max_dw, max_dw)
                
                theta_phi_next = theta_phi + self.K_phi * eta_phi + theta_dw
                
                theta_a = theta_a_next
                theta_dw = theta_dw_next
                theta_phi = theta_phi_next
                
        return e_output