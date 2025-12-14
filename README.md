# Robust Fetal ECG Detection for Abdominal Recordings

## Overview

This project replicates the sequential analysis method proposed by Martens et al. for detecting fetal electrocardiogram (fECG) signals from non-invasive abdominal recordings. The implementation addresses the challenge of extracting weak fetal cardiac signals overwhelmed by maternal ECG (mECG), baseline wander, power-line interference, and electromyographic (EMG) noise.

## Clinical Motivation

Monitoring fetal cardiac activity provides critical information for assessing fetal well-being during pregnancy and labor. While Doppler ultrasound is commonly used, abdominal fECG offers several advantages:

- **Non-invasive** monitoring throughout pregnancy
- **Rich morphological information** (P wave, QRS complex, T wave, ST segment)
- **Long-term ambulatory monitoring** potential
- **Multiple parameters** for fetal distress assessment

## Methodology

The project implements a **sequential analysis pipeline** that systematically removes interference signals using a priori knowledge about their characteristics:

### Pipeline Stages

```
Raw Signal (S1) 
    ↓
1. Baseline Wander Removal (S2)
    ↓
2. Power-Line Interference Cancellation (S3)
    ↓
3. Upsampling to 2000 Hz (S4)
    ↓
4. Maternal QRS Detection → MECG Cancellation (S5)
    ↓
5. Fetal QRS Detection → FECG Extraction (S6)
```

### 1. Baseline Wander Remover

**Problem:** Respiratory movement and patient motion cause low-frequency drift (< 3 Hz) that overlaps with fECG spectrum.

**Solution:** 
- Linear-phase FIR high-pass filter (1000 taps, 3 Hz cutoff, Hamming window)
- Preserves fECG morphology while attenuating baseline drift
- Zero-phase filtering via `filtfilt` prevents signal distortion

### 2. Adaptive Power-Line Interference Canceller

**Problem:** 50 Hz power-line noise and harmonics corrupt recordings.

**Solution:**
- Adaptive cancellation with real-time estimation of amplitude, frequency, and phase
- **Blocking detection** protects QRS complexes during adaptation:
  - Comb filter removes periodic interference
  - Threshold-based detection identifies signal regions
  - Mask prevents adaptation during high-energy cardiac events
- LMS-based parameter updates maintain >30 dB suppression

### 3. Sampling Rate Increaser

**Rationale:** Accurate MECG template matching requires high temporal resolution.

**Implementation:**
- Upsampling to 2000 Hz using polyphase filtering
- Prevents aliasing while enabling precise R-peak alignment
- Critical for effective MECG subtraction (reduces residual maternal artifacts)

### 4. MECG Canceller

1. **Maternal QRS Detection (Two-Stage Approach):**
   - PCA extracts primary cardiac component
   - Hilbert envelope with dynamic thresholding locates R-peaks
   - Cross-correlation refinement ensures sample-accurate alignment

2. **Morphology-Aware Cancellation:**
   - Template averaging over N=10 beats (balances fECG rejection and MECG tracking)
   - **Segmented scaling** independently fits P wave, QRS complex, and T wave
   - Ridge regression (`λ = 1e-3`) handles time-varying morphology
   
3. **Mathematical Formulation:**
   ```
   Template Matrix M = [P_wave | QRS | T_wave]
   Fitted beat = M · α
   where α = (M^T M + λI)^(-1) M^T · current_beat
   Residual = current_beat - fitted_beat
   ```

### 5. Fetal ECG Extractor

**Two Detection Strategies:**

#### Energy-Based Detection (Primary)
- Cubic transformation enhances R-peaks: `energy = signal³`
- Percentile-based adaptive thresholding (98th percentile)
- Minimum distance constraint: 0.25s (supports 120-240 bpm range)

#### PCA-Hilbert Detection (Alternative)
- Multi-channel PCA identifies dominant fetal component
- Hilbert envelope extracts instantaneous amplitude
- Gaussian smoothing (15 ms window) + dynamic thresholding

**Beat Averaging:**
- Synchronization on detected R-peaks (250 ms window)
- Ensemble averaging over all beats improves SNR by factor √N
- Baseline-corrected segments ensure morphological fidelity

## Project Structure

```
├── pipeline.ipynb              # Main analysis notebook
├── pipeline_functions.py       # Core processing functions
├── adaptivePLICanceller.py    # Power-line interference removal
├── MECGCanceller.py           # Maternal ECG cancellation
├── FetalECGExtractor.py       # Fetal QRS detection & averaging
├── Paper3.pdf                 # Reference paper (Martens et al.)
└── BSP_Project_Pizzocri.pdf   # Project presentation
```

## Implementation Details

### Dependencies

```python
numpy          # Numerical computing
scipy          # Signal processing (butter, filtfilt, hilbert, resample_poly)
sklearn        # PCA, FastICA for decomposition
pandas         # Data manipulation and rolling statistics
matplotlib     # Visualization
wfdb           # PhysioNet database access
```

### Key Classes

**`AdaptivePLICanceller`**
- Real-time adaptive filtering with blocking detection
- Comb filter for interference-free energy estimation
- Sample-by-sample LMS updates with convergence safeguards

**`MECGCanceller`**
- Two-stage maternal R-peak detection (rough → refined)
- Adaptive template with P-QRS-T segmentation
- Ridge regression for robust fitting under morphological variation

**`FetalECGExtractor`**
- Dual detection modes (energy-based + PCA-Hilbert)
- Synchronized beat averaging with outlier rejection
- Multi-channel morphology extraction

## Performance Metrics

### Detection Success Rate
- **Sensitivity (Se):** TP / (TP + FN) × 100%
- **Precision (PPV):** TP / (TP + FP) × 100%
- **F1 Score:** 2 · (PPV · Se) / (PPV + Se)
- **Tolerance:** ±50 ms window for peak matching

### Signal Quality (Paper Metrics)
```
SNR = 10 log₁₀(P_FECG / P_noise)       [Eq. 4]
SIR = 10 log₁₀(P_FECG / P_MECG)        [Eq. 5]
```

### FHR Reliability
```
Reliability = 1 - (outliers / total_points)
Outlier: |FHR - median₁₀ₛ| > 10 bpm
```

## Usage Example

```python
from pipeline_functions import load_and_interpolate, preprocess_signal, extract_and_analyze
from FetalECGExtractor import FetalECGExtractor

# Load data
raw_signal, fields = load_and_interpolate('path/to/record')
fs_native = fields['fs']

# Process through pipeline
S5, S4, up_factor = preprocess_signal(raw_signal, fs_native, target_fs=2000)

# Extract fetal ECG
extractor = FetalECGExtractor(fs=2000)
fetal_signal, avg_beat, f_peaks = extract_and_analyze(S5, extractor, 2000)

print(f"Detected {len(f_peaks)} fetal heartbeats")
```

## References

**Primary Paper:**
Martens, S. M. M., Rabotti, C., Mischi, M., & Sluijter, R. J. (2007). A robust fetal ECG detection method for abdominal recordings. *Physiological Measurement*, 28(4), 373-388.

**Dataset:**
PhysioNet Non-Invasive Fetal ECG Database (https://physionet.org/content/nifecgdb/)

## License

This project is an academic implementation for educational purposes. Please refer to the original paper for methodological details and cite appropriately.
