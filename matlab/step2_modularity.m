% ---- paths: set these environment variables (see README) ----
project_dir = getenv('UPF_PROJECT_DIR');   % project root (data/, outputs/, ROIs/)
np2_dir     = getenv('NP2_DATA_DIR');      % NeuroPharm2 derived data (restricted access)
bct_dir     = getenv('BCT_DIR');           % Brain Connectivity Toolbox
roi_dir     = getenv('ROI_DIR');           % Schaefer atlas files

%% Standalone Modularity Analysis Script
% This script calculates modularity metrics for FC and forward MI matrices
% and saves the full null distributions for comparison with GCOR
%
% Key features:
% - Handles both signed (FC) and unsigned (forward MI) matrices
% - Saves full null distribution for each subject
% - Calculates normalised and absolute modularity
% - Exports results to R-friendly format

clear all
close all

%% Configuration
cd(project_dir);
addpath(genpath(project_dir));
addpath(bct_dir)

% Modularity parameters
n_null = 100;        % Number of null model iterations
resolution = 1.0;    % Resolution parameter
ncores = 32;         % Number of parallel cores

% Output directory
out_dir = fullfile(project_dir, 'outputs');
timestamp = datestr(now, 'yyyy-mm-dd');
out_dir = fullfile(out_dir, timestamp);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
    fprintf('Created output directory: %s\n', out_dir);
end

%% Load data
fprintf('Loading master_tbl_combined...\n');
load(fullfile(project_dir, 'outputs/master_tbl_combined.mat'));

n_subjects = height(master_tbl_combined);
fprintf('Loaded data for %d subjects\n', n_subjects);

%% Initialize storage for modularity results

% We'll calculate modularity for:
% 1. FC_matrix (static FC - signed, undirected)
% 2. forward_MI (directed MI - unsigned, directed)

% Initialize result tables
results_fc = table();
results_fmi = table();

% Add subject identifiers
results_fc.CIMBI_ID = master_tbl_combined.CIMBI_ID;
results_fc.scanner = master_tbl_combined.scanner;
results_fc.psi_conc = master_tbl_combined.psi_conc;
results_fc.time_since_admin = master_tbl_combined.time_since_admin;

results_fmi.CIMBI_ID = master_tbl_combined.CIMBI_ID;
results_fmi.scanner = master_tbl_combined.scanner;
results_fmi.psi_conc = master_tbl_combined.psi_conc;
results_fmi.time_since_admin = master_tbl_combined.time_since_admin;

% Initialize columns for storing distributions (as cells containing vectors)
results_fc.actual_mod_distribution = cell(n_subjects, 1);
results_fc.null_mod_distribution = cell(n_subjects, 1);
results_fc.modularity_abs = zeros(n_subjects, 1);
results_fc.modularity_norm = zeros(n_subjects, 1);

results_fmi.actual_mod_distribution = cell(n_subjects, 1);
results_fmi.null_mod_distribution = cell(n_subjects, 1);
results_fmi.modularity_abs = zeros(n_subjects, 1);
results_fmi.modularity_norm = zeros(n_subjects, 1);

%% Start parallel pool
if isempty(gcp('nocreate'))
    parpool(ncores);
end

%% Calculate modularity for each subject

start_time = tic;
fprintf('\n=== CALCULATING MODULARITY ===\n');
fprintf('Starting calculations for %d subjects...\n\n', n_subjects);

for subj = 1:n_subjects
    if subj > 1
        elapsed_time = toc(start_time);
        estimated_total_time = elapsed_time / (subj-1) * n_subjects;
        time_remaining = estimated_total_time - elapsed_time;
        fprintf('\rProcessing subject %d/%d - Est. time remaining: %.1f min', ...
            subj, n_subjects, time_remaining/60);
    else
        fprintf('Processing subject %d/%d - Calculating...', subj, n_subjects);
    end
    
    %% 1. Static FC (signed, undirected)
    FC_mat = master_tbl_combined.FC_matrix{subj};
    [actual_mods_fc, null_mods_fc] = calculate_modularity_with_nulls(...
        FC_mat, n_null, resolution, 'signed', false);
    
    results_fc.actual_mod_distribution{subj} = actual_mods_fc;
    results_fc.null_mod_distribution{subj} = null_mods_fc;
    results_fc.modularity_abs(subj) = max(actual_mods_fc);
    results_fc.modularity_norm(subj) = max(actual_mods_fc) / mean(null_mods_fc);
    
    %% 2. Forward MI (unsigned, directed)
    FMI_mat = master_tbl_combined.forward_MI{subj};
    [actual_mods_fmi, null_mods_fmi] = calculate_modularity_with_nulls(...
        FMI_mat, n_null, resolution, 'unsigned', true);
    
    results_fmi.actual_mod_distribution{subj} = actual_mods_fmi;
    results_fmi.null_mod_distribution{subj} = null_mods_fmi;
    results_fmi.modularity_abs(subj) = max(actual_mods_fmi);
    results_fmi.modularity_norm(subj) = max(actual_mods_fmi) / mean(null_mods_fmi);
end

fprintf('\n\nCompleted all subjects! Total time: %.1f min\n', toc(start_time)/60);

%% Save MATLAB results
fprintf('\nSaving MATLAB results...\n');
save([out_dir, '/modularity_results_fc.mat'], 'results_fc');
save([out_dir, '/modularity_results_fmi.mat'], 'results_fmi');
fprintf('Saved .mat files\n');

%% Export to R-friendly format

fprintf('\nExporting to R-friendly format...\n');

% Create summary table for R (without full distributions)
summary_fc = table(...
    results_fc.CIMBI_ID, ...
    results_fc.scanner, ...
    results_fc.psi_conc, ...
    results_fc.time_since_admin, ...
    results_fc.modularity_abs, ...
    results_fc.modularity_norm, ...
    'VariableNames', {'CIMBI_ID', 'scanner', 'psi_conc', 'time_since_admin', ...
                      'modularity_abs', 'modularity_norm'});

summary_fmi = table(...
    results_fmi.CIMBI_ID, ...
    results_fmi.scanner, ...
    results_fmi.psi_conc, ...
    results_fmi.time_since_admin, ...
    results_fmi.modularity_abs, ...
    results_fmi.modularity_norm, ...
    'VariableNames', {'CIMBI_ID', 'scanner', 'psi_conc', 'time_since_admin', ...
                      'modularity_abs', 'modularity_norm'});

% Save summary tables
writetable(summary_fc, [out_dir, '/modularity_summary_fc.csv']);
writetable(summary_fmi, [out_dir, '/modularity_summary_fmi.csv']);

% Export null distributions in long format for R
fprintf('Exporting null distributions in long format...\n');

% FC null distributions
export_null_distributions(results_fc, [out_dir, '/null_distributions_fc.csv'], 'FC');
% Forward MI null distributions
export_null_distributions(results_fmi, [out_dir, '/null_distributions_fmi.csv'], 'FMI');

fprintf('\n=== EXPORT COMPLETE ===\n');
fprintf('Files saved to: %s\n', out_dir);
fprintf('\nSummary files:\n');
fprintf('  - modularity_summary_fc.csv\n');
fprintf('  - modularity_summary_fmi.csv\n');
fprintf('\nNull distribution files:\n');
fprintf('  - null_distributions_fc.csv\n');
fprintf('  - null_distributions_fmi.csv\n');

%% Helper Functions

function [actual_mods, null_mods] = calculate_modularity_with_nulls(...
    conn_matrix, n_null, resolution, network_type, is_directed)
    % Calculate modularity with full null distribution
    %
    % Inputs:
    %   conn_matrix: connectivity matrix
    %   n_null: number of null iterations
    %   resolution: resolution parameter
    %   network_type: 'signed' or 'unsigned'
    %   is_directed: true/false
    %
    % Outputs:
    %   actual_mods: vector of actual modularity values
    %   null_mods: vector of null modularity values
    
    % Calculate actual network modularity
    actual_mods = zeros(1, n_null);
    parfor i = 1:n_null
        if strcmp(network_type, 'signed')
            [~, Qa] = community_louvain(conn_matrix, resolution, [], 'negative_asym');
        else
            [~, Qa] = community_louvain(conn_matrix, resolution);
        end
        actual_mods(i) = Qa;
    end
    
    % Calculate null model modularity
    null_mods = zeros(1, n_null);
    parfor i = 1:n_null
        % Generate null model
        if is_directed
            if strcmp(network_type, 'signed')
                null_matrix = null_model_dir_sign(conn_matrix);
            else
                % For unsigned directed, use signed null but matrix is positive
                null_matrix = null_model_dir_sign(conn_matrix);
            end
        else
            if strcmp(network_type, 'signed')
                null_matrix = null_model_und_sign(conn_matrix);
            else
                % For unsigned undirected, use signed null but matrix is positive
                null_matrix = null_model_und_sign(conn_matrix);
            end
        end
        
        % Calculate modularity for this null
        null_mod_reps = zeros(1, n_null);
        for j = 1:n_null
            if strcmp(network_type, 'signed')
                [~, Qn] = community_louvain(null_matrix, resolution, [], 'negative_asym');
            else
                [~, Qn] = community_louvain(null_matrix, resolution);
            end
            null_mod_reps(j) = Qn;
        end
        % Take the maximum across repetitions
        null_mods(i) = max(null_mod_reps);
    end
end

function export_null_distributions(results_table, filename, network_name)
    % Export null distributions in long format for R
    % Each row is one modularity value from one subject
    
    n_subjects = height(results_table);
    
    % Pre-allocate cell arrays for efficiency
    all_cimbi_ids = {};
    all_scanners = {};
    all_psi_conc = [];
    all_time_since_admin = [];
    all_distribution_type = {};
    all_modularity_values = [];
    
    for subj = 1:n_subjects
        if mod(subj, 20) == 0
            fprintf('  Processing subject %d/%d for %s...\n', subj, n_subjects, network_name);
        end
        
        % Get subject info
        cimbi_id = char(results_table.CIMBI_ID(subj));
        scanner = results_table.scanner{subj};
        psi_conc = results_table.psi_conc(subj);
        time_since = results_table.time_since_admin(subj);
        
        % Get actual distribution
        actual_dist = results_table.actual_mod_distribution{subj};
        n_actual = length(actual_dist);
        
        all_cimbi_ids = [all_cimbi_ids; repmat({cimbi_id}, n_actual, 1)];
        all_scanners = [all_scanners; repmat({scanner}, n_actual, 1)];
        all_psi_conc = [all_psi_conc; repmat(psi_conc, n_actual, 1)];
        all_time_since_admin = [all_time_since_admin; repmat(time_since, n_actual, 1)];
        all_distribution_type = [all_distribution_type; repmat({'actual'}, n_actual, 1)];
        all_modularity_values = [all_modularity_values; actual_dist(:)];
        
        % Get null distribution
        null_dist = results_table.null_mod_distribution{subj};
        n_null = length(null_dist);
        
        all_cimbi_ids = [all_cimbi_ids; repmat({cimbi_id}, n_null, 1)];
        all_scanners = [all_scanners; repmat({scanner}, n_null, 1)];
        all_psi_conc = [all_psi_conc; repmat(psi_conc, n_null, 1)];
        all_time_since_admin = [all_time_since_admin; repmat(time_since, n_null, 1)];
        all_distribution_type = [all_distribution_type; repmat({'null'}, n_null, 1)];
        all_modularity_values = [all_modularity_values; null_dist(:)];
    end
    
    % Create table
    export_table = table(...
        all_cimbi_ids, ...
        all_scanners, ...
        all_psi_conc, ...
        all_time_since_admin, ...
        all_distribution_type, ...
        all_modularity_values, ...
        'VariableNames', {'CIMBI_ID', 'scanner', 'psi_conc', 'time_since_admin', ...
                          'distribution_type', 'modularity'});
    
    % Save
    writetable(export_table, filename);
    fprintf('  Exported %d rows to %s\n', height(export_table), filename);
end