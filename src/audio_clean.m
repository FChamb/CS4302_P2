infiles = {'../audio/audio_in_noise1.wav','../audio/audio_in_noise2.wav','../audio/audio_in_noise3.wav'};
outdir  = '../out'; if ~exist(outdir,'dir'), mkdir(outdir); end

function [f, Xdb, pkHz] = analyse_noise_tones(x, fs)
    Nfft = 2^nextpow2(max(length(x), 4*fs));
    w = hann(min(length(x), 4*fs));
    xw = x(1:length(w)).*w;
    X = fft(xw, Nfft);
    f = (0:Nfft/2)' * fs/Nfft;
    Xmag = abs(X(1:Nfft/2+1));
    Xdb  = 20*log10(Xmag/max(Xmag)+eps);

    [~,locs] = findpeaks(Xdb, 'MinPeakProminence', 6, 'MinPeakDistance', round(0.5*length(f)/ (fs/2)*100)); 
    pkHz = f(locs);
end

for k = 1:numel(infiles)
    [x, fs] = audioread(infiles{k});
    x = mean(x,2);
    x = x / max(abs(x)+eps);

    [f, Xdb, pkHz] = analyse_noise_tones(x, fs);

    figure('Name',sprintf('Spectrum_%d',k),'Visible','off');
    plot(f, Xdb); grid on; xlim([0 fs/2]);
    xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)');
    title(sprintf('Magnitude Spectrum for %s', infiles{k}));
    hold on; yL = ylim; 
    for p = pkHz(:)', line([p p], yL, 'LineStyle','--'); end
    legend('Spectrum','Detected peaks');
    saveas(gcf, fullfile(outdir, sprintf('audio_spectrum_%d.png',k)));

    fprintf('File %d: detected noise tones at ~ [ %s ] Hz\n', k, sprintf('%.1f ', pkHz));
end

for k = 1:numel(infiles)
    [x, fs] = audioread(infiles{k}); x = mean(x,2); x = x/max(abs(x)+eps);
    [~, ~, pkHz] = analyse_noise_tones(x, fs);

    y = x; sos_all = []; g_all = 1;
    for f0 = pkHz(:)'
        w0 = f0/(fs/2); Q = 45;
        [b,a] = iirnotch(w0, w0/Q);
        [sos,g] = tf2sos(b,a);
        sos_all = [sos_all; sos]; g_all = g_all*g;
        y = filtfilt(b,a,y);
    end

    outf = fullfile(outdir, sprintf('noise_removed%d.wav', k));
    audiowrite(outf, y, fs, 'BitsPerSample',16);

    t = (0:numel(x)-1)/fs; Tplot = t<=0.08;
    figure('Name',sprintf('Time_%d',k),'Visible','off');
    plot(t(Tplot), x(Tplot)); hold on; plot(t(Tplot), y(Tplot));
    grid on; xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Time Domain Before/After (file %d)',k));
    legend('Original','Filtered');
    saveas(gcf, fullfile(outdir, sprintf('audio_time_%d.png',k)));

    Nfft = 2^nextpow2(4*fs);
    [f, Xdb] = local_fft_db(x, fs, Nfft);
    [~, Ydb] = local_fft_db(y, fs, Nfft);
    figure('Name',sprintf('Freq_%d',k),'Visible','off');
    plot(f,Xdb); hold on; plot(f,Ydb); grid on; xlim([0 fs/2]);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title(sprintf('Spectrum Before/After (file %d)',k));
    legend('Original','Filtered');
    saveas(gcf, fullfile(outdir, sprintf('audio_freq_%d.png',k)));
end

function [f, Xdb] = local_fft_db(x, fs, Nfft)
    x = x - mean(x);
    X = fft(hann(min(length(x),4*fs)).*x(1:min(length(x),4*fs)), Nfft);
    f = (0:Nfft/2)' * fs/Nfft;
    X = abs(X(1:Nfft/2+1));
    Xdb = 20*log10(X/max(X)+eps);
end