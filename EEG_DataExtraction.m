%% EEG Band Power + Frontal Theta / Posterior Alpha Ratio (MIPDB preprocessed .mat)
% Robust version:
%  - Uses FFT-based PSD (no pwelch / no Signal Processing Toolbox)
%  - Tries to find frontal/posterior channels by standard 10–20 labels
%  - If labels don't match (e.g., E1..E111), auto-selects frontal/posterior
%    using chanlocs Y-coordinates (front/back of head)
%  - PSD plot falls back to channel 1 if requested label doesn't exist
%
% Outputs:
%  - alpha_power: 1 x nCh absolute alpha power (8–12 Hz) per channel
%  - theta_power: 1 x nCh absolute theta power (4–8 Hz) per channel
%  - posterior_alpha: average alpha power over posterior channels
%  - frontal_theta: average theta power over frontal channels
%  - ratio: frontal_theta / posterior_alpha

clear; clc;

%% 1) Load data
filename = 'bp_A00054432006.mat';   % <-- change if needed
S = load(filename);

% Many MIPDB preprocessed files store data in an EEGLAB-like struct called EEG
EEG = S.EEG;

% Extract signal + sampling rate
signal = double(EEG.data);    % [nCh x nSamples]
fs = EEG.srate;               % sampling rate (Hz)
labels = {EEG.chanlocs.labels};

% Quick sanity checks
[nCh, nSamples] = size(signal);
dur_sec = nSamples / fs;

fprintf('Loaded %s\n', filename);
fprintf('Channels: %d | Samples: %d | fs: %g Hz | Duration: %.1f sec (%.2f min)\n', ...
    nCh, nSamples, fs, dur_sec, dur_sec/60);

%% 2) Compute PSD using FFT (NO pwelch / NO toolboxes)
% Output: pxx = [nFreq x nCh], f = [nFreq x 1]
x = signal.';                % [nSamples x nCh]
n = size(x,1);

% Remove DC offset per channel (recommended)
x = x - mean(x,1);

% FFT length (next power of 2 for speed)
nfft = 2^nextpow2(n);

% FFT along time dimension
X = fft(x, nfft, 1);

% Two-sided power (unnormalized)
P2 = (abs(X)/n).^2;

% One-sided power spectrum
pxx = P2(1:nfft/2+1, :);
pxx(2:end-1,:) = 2*pxx(2:end-1,:);

% Frequency vector
f = fs*(0:(nfft/2))'/nfft;

%% 3) Compute absolute band power per channel (mean power in band)
alpha_idx = (f >= 8) & (f <= 12);
theta_idx = (f >= 4) & (f <= 8);

alpha_power = mean(pxx(alpha_idx, :), 1);   % 1 x nCh
theta_power = mean(pxx(theta_idx, :), 1);   % 1 x nCh

%% 4) Pick posterior + frontal channel groups
% First attempt: standard 10–20 labels
posterior_names = {'O1','O2','Oz','Pz','P3','P4','POz'};
frontal_names   = {'Fz','F3','F4','Fp1','Fp2','FCz'};

posterior_idx = find(ismember(labels, posterior_names));
frontal_idx   = find(ismember(labels, frontal_names));

% If label matching fails, auto-pick using chanlocs Y-coordinates
% (front/back axis). Posterior = lowest Y, Frontal = highest Y.
if isempty(posterior_idx) || isempty(frontal_idx)
    warning('10-20 label match failed. Attempting AUTO channel selection using chanlocs coordinates...');

    chanlocs = EEG.chanlocs;

    % Extract Y coordinate robustly (handles numeric or string)
    Y = nan(1, numel(chanlocs));
    for k = 1:numel(chanlocs)
        if isfield(chanlocs(k),'Y')
            yk = chanlocs(k).Y;
            if ischar(yk) || isstring(yk)
                yk = str2double(yk);
            end
            Y(k) = yk;
        end
    end

    valid = ~isnan(Y);

    if sum(valid) < 10
        warning('Not enough valid Y-coordinates in chanlocs to auto-select regions. Results may be NaN.');
        posterior_idx = [];
        frontal_idx = [];
    else
        % Pick top/bottom 10% (at least 5 channels)
        kN = max(5, round(0.10 * sum(valid)));

        validIdx = find(valid);
        [~, sortIdx] = sort(Y(valid), 'ascend');   % ascending: low = posterior, high = frontal

        posterior_idx = validIdx(sortIdx(1:kN));
        frontal_idx   = validIdx(sortIdx(end-kN+1:end));
    end
end

%% 5) Compute regional averages + ratio (guard against empty selections)
if isempty(posterior_idx)
    posterior_alpha = NaN;
else
    posterior_alpha = mean(alpha_power(posterior_idx));
end

fprintf('posterior_alpha :   (%d)\n', posterior_alpha)

if isempty(frontal_idx)
    frontal_theta = NaN;
else
    frontal_theta = mean(theta_power(frontal_idx));
end
 
fprintf('frontal theta:   (%d)\n', frontal_theta)

ratio = frontal_theta / posterior_alpha;
fprintf('ratio:   (%d)\n', ratio)
fprintf('\nPosterior channels found (%d)\n', numel(posterior_idx));
if ~isempty(posterior_idx), disp(labels(posterior_idx)); 
end

fprintf('Frontal channels found   (%d)\n', numel(frontal_idx))