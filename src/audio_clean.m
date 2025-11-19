% ============================
% CS4302 Practical 2 - Audio
% File: audio_clean.m
% ============================
% What this script does:
% 1) Loads each noisy audio file
% 2) Estimates noise frequency/frequencies from the magnitude spectrum
% 3) Designs notch filters and removes the noise
% 4) Shows frequency/time plots before & after
% 5) Saves cleaned audio as: noise_removed<1/2/3>.wav
%

clear; close all; clc;

% ----- Parameters -----
peakProminence = 0.02;   % relative prominence for spectral peaks
minFreqHz = 50;     % ignore very low frequency/DC drift
maxPeaks = 3;      % maximum peaks to label as noise
Q = 35;    % notch filter quality factor

% Input file names
inFiles = { ...
    '../audio/audio_in_noise1.wav', ...
    '../audio/audio_in_noise2.wav', ...
    '../audio/audio_in_noise3.wav'  ...
};

% Output folder
outDir = '../out';
if ~exist(outDir, 'dir'), mkdir(outDir); end

for k = 1:numel(inFiles)
    [x, fs] = audioread(inFiles{k});
    x = mean(x, 2); % convert to mono if stereo

    % ----- Frequency analysis (magnitude spectrum) -----
    Nfft = 2^nextpow2(numel(x));
    X = fft(x, Nfft);
    f = (0:Nfft-1)*(fs/Nfft);
    mag = abs(X)/max(abs(X)); % normalize for easier peak selection

    % consider only up to Nyquist and ignore the DC bin
    halfIdx = 1:floor(Nfft/2);
    fHalf = f(halfIdx);
    magHalf = mag(halfIdx);

    % ignore very low frequencies (below minFreqHz)
    valid = fHalf >= minFreqHz;
    fCand = fHalf(valid);
    mCand = magHalf(valid);

    % find spectral peaks that stand out
    [pks, locs] = findpeaks(mCand, 'MinPeakProminence', peakProminence);
    % sort by peak height and take the top few
    [~, order] = sort(pks, 'descend');
    order = order(1:min(maxPeaks, numel(order)));
    noiseFreqs = fCand(locs(order));
    noiseFreqs = round(noiseFreqs, 1); % neat formatting

    % ----- Plot spectrum with detected noise frequencies -----
    fig1 = figure('Name', sprintf('Spectrum_%d', k), 'Color', 'w');
    plot(fHalf, magHalf, 'LineWidth', 1); grid on;
    xlabel('Frequency (Hz)'); ylabel('Normalized magnitude');
    title(sprintf('Magnitude Spectrum (%s)', inFiles{k}), 'Interpreter', 'none');
    hold on;
    stem(noiseFreqs, interp1(fHalf, magHalf, noiseFreqs, 'linear', 'extrap'), 'r', 'LineWidth', 1.2);
    legend('Spectrum', 'Detected noise peak(s)');
    saveas(fig1, fullfile(outDir, sprintf('spectrum_%d.png', k)));

    % ----- Build cascaded notch filters for each detected noise frequency -----
    xFilt = x;
    for f0 = noiseFreqs(:).'
        bwHz = f0 / Q;               % approximate notch width
        f1 = (f0 - bwHz/2) / (fs/2); % normalized lower band edge
        f2 = (f0 + bwHz/2) / (fs/2); % normalized upper band edge

        % clamp edges to valid range
        f1 = max(f1, 0.0001);
        f2 = min(f2, 0.9999);

        % 2nd order bandstop
        [b, a] = butter(2, [f1, f2], 'stop');

        % zero phase filtering
        xFilt = filtfilt(b, a, xFilt);
    end

    % ----- Time domain comparison plot -----
    nShow = min(3*fs, numel(x)); % show up to 3 seconds
    t = (0:nShow-1)/fs;
    fig2 = figure('Name', sprintf('Time_%d', k), 'Color', 'w');
    subplot(2,1,1);
    plot(t, x(1:nShow), 'LineWidth', 1); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Original (first %.2f s) — %s', nShow/fs, inFiles{k}), 'Interpreter', 'none');

    subplot(2,1,2);
    plot(t, xFilt(1:nShow), 'LineWidth', 1); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Filtered (same window)');
    saveas(fig2, fullfile(outDir, sprintf('time_compare_%d.png', k)));

    % ----- Frequency domain comparison after filtering -----
    Xf = fft(xFilt, Nfft);
    magHalfFilt = abs(Xf(halfIdx)); magHalfFilt = magHalfFilt/max(magHalfFilt);

    fig3 = figure('Name', sprintf('SpectrumAfter_%d', k), 'Color', 'w');
    plot(fHalf, magHalf, 'LineWidth', 1); hold on;
    plot(fHalf, magHalfFilt, 'LineWidth', 1);
    grid on; xlabel('Frequency (Hz)'); ylabel('Normalized magnitude');
    title('Before vs After Filtering (magnitude spectra)');
    legend('Before', 'After');
    saveas(fig3, fullfile(outDir, sprintf('spectrum_compare_%d.png', k)));

    % ----- Save the cleaned file -----
    outName = sprintf('noise_removed%d.wav', k);
    audiowrite(outDir + "/" + outName, xFilt, fs);

    % ----- Console summary -----
    fprintf('File %d: %s\n', k, inFiles{k});
    fprintf('  Detected noise frequencies (Hz): %s\n', mat2str(noiseFreqs));
    fprintf('  Saved cleaned audio to: %s\n\n', outName);
end
