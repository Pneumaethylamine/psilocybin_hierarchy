% ---- paths: set these environment variables (see README) ----
project_dir = getenv('UPF_PROJECT_DIR');   % project root (data/, outputs/, ROIs/)
np2_dir     = getenv('NP2_DATA_DIR');      % NeuroPharm2 derived data (restricted access)
bct_dir     = getenv('BCT_DIR');           % Brain Connectivity Toolbox
roi_dir     = getenv('ROI_DIR');           % Schaefer atlas files

function [mean_autocorr, se_autocorr, time] = calculate_autocorrelations(data_table, ts_column_name, maxlag, TR)
    % Calculate autocorrelations for time series data
    %
    % Inputs:
    %   data_table: Table containing the time series data
    %   ts_column_name: Name of the column containing time series cell arrays
    %   maxlag: Maximum lag for autocorrelation calculation
    %   TR: Repetition time in seconds
    %
    % Outputs:
    %   mean_autocorr: Mean autocorrelation across all subjects
    %   se_autocorr: Standard error of autocorrelation
    %   time: Time vector corresponding to lags
    
    % Input validation
    if ~istable(data_table)
        error('First argument must be a table');
    end
    if ~ischar(ts_column_name) && ~isstring(ts_column_name)
        error('Second argument must be a string specifying the column name');
    end
    if ~any(strcmp(data_table.Properties.VariableNames, ts_column_name))
        error('Specified column name not found in table');
    end
    
    % Initialize matrix for autocorrelations
    autocorr_matrix = zeros(height(data_table), maxlag + 1);
    
    % Calculate autocorrelation for each subject
    for i = 1:height(data_table)
        ts = data_table.(ts_column_name){i};
        mean_ts = mean(ts, 2);  % Average across ROIs
        autocorr_matrix(i, :) = autocorr_manual(mean_ts, maxlag);
    end
    
    % Calculate statistics
    mean_autocorr = mean(autocorr_matrix, 1);
    se_autocorr = std(autocorr_matrix, 0, 1) / sqrt(size(autocorr_matrix, 1));
    
    % Create time vector
    time = (0:maxlag) * TR;
end