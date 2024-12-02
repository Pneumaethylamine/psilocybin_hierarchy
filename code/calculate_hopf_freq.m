function [f_diff, fce] = calculate_hopf_freq(scanner_table, TR, Tmax, make_plots)
% Calculate Hopf frequencies from timeseries data
%
% Inputs:
% scanner_table: table containing timeseries data in ts_SchaeferTian232
% TR: Repetition Time in seconds
% Tmax: Number of timepoints
%
% Outputs:
% f_diff: Peak frequencies for each ROI
% fce: Mean functional connectivity matrix

NPARCELLS = 200;  % After removing ROIs 201-232
NSUB = height(scanner_table);

% Filter settings
fnq = 1/(2*TR);
flp = 0.008;
fhi = 0.08;
Wn = [flp/fnq fhi/fnq];
k = 2;
[bfilt,afilt] = butter(k,Wn);

% Initialize arrays
fce1 = zeros(NSUB, NPARCELLS, NPARCELLS);
PowSpect = zeros(floor(Tmax/2), NPARCELLS, NSUB);

% Process each subject
for sub = 1:NSUB
    fprintf('Processing subject %d of %d\n', sub, NSUB)

    % Get subject timeseries and remove ROIs 201-232
    ts = scanner_table.ts_SchaeferTian232{sub};
    ts = ts';
    ts = ts(1:200, :);

    TT = Tmax;
    tss = zeros(NPARCELLS, Tmax);
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs = length(freq);

    % Process each ROI
    x = zeros(NPARCELLS, Tmax);  % Preallocate x
    for seed = 1:NPARCELLS
        x(seed,:) = detrend(ts(seed,:)-mean(ts(seed,:)));
        tss(seed,:) = filtfilt(bfilt,afilt,x(seed,:));
        pw = abs(fft(tss(seed,:)));
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    fce1(sub,:,:) = corrcoef(tss','rows','pairwise');
end

% Calculate mean FC
fce = squeeze(mean(fce1));

% Calculate power spectra and frequencies
Power_Areas = squeeze(mean(PowSpect,3));
for seed = 1:NPARCELLS
    Power_Areas(:,seed) = gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

% Get peak frequencies
[maxpowdata,index] = max(Power_Areas);
f_diff = freq(index);
f_diff(f_diff==0) = mean(f_diff(f_diff~=0));

if make_plots
    % Create figure for power spectra
    figure('Position', [100 100 1200 800]);

    % Plot 1: All power spectra
    subplot(2,2,1)
    plot(freq, Power_Areas, 'LineWidth', 1, 'Color', [0.7 0.7 0.7 0.3])
    hold on
    plot(freq, mean(Power_Areas,2), 'k-', 'LineWidth', 2)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title('Power Spectra for All ROIs')
    grid on

    % Plot 2: Mean power spectrum with peak
    subplot(2,2,2)
    mean_power = mean(Power_Areas,2);
    [max_power, max_idx] = max(mean_power);
    plot(freq, mean_power, 'b-', 'LineWidth', 2)
    hold on
    plot(freq(max_idx), max_power, 'r*', 'MarkerSize', 15)
    text(freq(max_idx), max_power, sprintf(' Peak: %.3f Hz', freq(max_idx)), 'VerticalAlignment', 'bottom')
    xlabel('Frequency (Hz)')
    ylabel('Mean Power')
    title('Mean Power Spectrum with Peak')
    grid on

    % Plot 3: Histogram of peak frequencies
    subplot(2,2,3)
    histogram(f_diff, 30)
    xlabel('Peak Frequency (Hz)')
    ylabel('Count')
    title('Distribution of Peak Frequencies Across ROIs')
    grid on

    % Plot 4: Brain regions colored by frequency
    subplot(2,2,4)
    scatter(1:NPARCELLS, f_diff, 50, f_diff, 'filled')
    colorbar
    xlabel('ROI Index')
    ylabel('Peak Frequency (Hz)')
    title('Peak Frequencies by ROI')
    grid on

    sgtitle(sprintf('Hopf Frequency Analysis (TR=%.1fs, Timepoints=%d)', TR, Tmax))

    % Save plot
    print(sprintf('hopf_freq_plots_TR%.1f_TP%d', TR, Tmax), '-dpng', '-r300')
end

end