% ---- paths: set these environment variables (see README) ----
project_dir = getenv('UPF_PROJECT_DIR');   % project root (data/, outputs/, ROIs/)
np2_dir     = getenv('NP2_DATA_DIR');      % NeuroPharm2 derived data (restricted access)
bct_dir     = getenv('BCT_DIR');           % Brain Connectivity Toolbox
roi_dir     = getenv('ROI_DIR');           % Schaefer atlas files

function plot_autocorrelations(mean_autocorr, se_autocorr, time, TR)
    % Plot autocorrelations for a single scanner
    %
    % Inputs:
    %   mean_autocorr: Mean autocorrelation for scanner
    %   se_autocorr: Standard error of autocorrelation
    %   time: Time vector
    %   scanner_name: Name of the scanner (e.g., 'Scanner P')
    %   TR: TR value in seconds
    
    figure('Position', [100 100 400 400]);
    errorbar(time, mean_autocorr, se_autocorr, '-', 'LineWidth', 2);
    xlabel('Time lag (seconds)');
    ylabel('Autocorrelation');
    title(sprintf('(TR = %.1fs)', TR));
    grid on;
    ylim([-0.2 1]);
    xlim([0 max(time)]);
end