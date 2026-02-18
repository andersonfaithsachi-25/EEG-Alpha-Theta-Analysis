## Overview
This project implements a MATLAB-based EEG spectral analysis pipeline for computing:
- Posterior alpha power (8–12 Hz)
- Frontal theta power (4–8 Hz)
- Theta/Alpha ratio biomarker

## Dataset
MIPDB preprocessed EEG dataset.

## Methods
- FFT-based power spectral density (no Signal Processing Toolbox required)
- Automatic frontal/posterior channel detection
- Regional band power averaging

## Outputs
- Alpha power per channel
- Theta power per channel
- Posterior alpha (regional)
- Frontal theta (regional)
- Theta/Alpha ratio

## Example Output
Theta/Alpha ratio ≈ 3.63 (example subject)
