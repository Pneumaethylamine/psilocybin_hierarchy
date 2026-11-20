% ---- paths: set these environment variables (see README) ----
project_dir = getenv('UPF_PROJECT_DIR');   % project root (data/, outputs/, ROIs/)
np2_dir     = getenv('NP2_DATA_DIR');      % NeuroPharm2 derived data (restricted access)
bct_dir     = getenv('BCT_DIR');           % Brain Connectivity Toolbox
roi_dir     = getenv('ROI_DIR');           % Schaefer atlas files

%% Setup and data loading
clear all
cd(project_dir);
addpath(genpath(project_dir));
addpath(bct_dir)

load(fullfile(project_dir, 'data/data_09-Mar-2023_SchaeferTian232_ts.mat'));
age_sex = readtable(fullfile(np2_dir, 'critical_files/CIMBI_age_sex_clean.xlsx'));

% Filter data
tbl = SDI_PPL_rsfmri_tbl(strcmpi(SDI_PPL_rsfmri_tbl.include, 'yes'), :);
%Remove PA PED scans
tbl = tbl(~contains(lower(tbl.rsfmri_location), 'pa'), :);

scanner_p_idx = startsWith(tbl.rsfmri_raw_ID, 'p');  % 2.0s TR scans
scanner_n_idx = startsWith(tbl.rsfmri_raw_ID, 'n');  % 0.8s TR scans
tbl_scanner_p = tbl(scanner_p_idx, :);
tbl_scanner_n = tbl(scanner_n_idx, :);

master_tbl_p = table(tbl_scanner_p.CIMBI_ID, tbl_scanner_p.psi_conc, ...
    tbl_scanner_p.time_since_admin, tbl_scanner_p.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

master_tbl_n = table(tbl_scanner_n.CIMBI_ID, tbl_scanner_n.psi_conc, ...
    tbl_scanner_n.time_since_admin, tbl_scanner_n.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

% Add age and sex
age_sex = age_sex(:, {'Cimbi', 'Age', 'Sex'});
missing_dude = table(55992, 33.0, "Male", 'VariableNames', {'Cimbi', 'Age', 'Sex'});
age_sex = [age_sex; missing_dude];
age_sex = sortrows(age_sex, 'Cimbi');
age_sex.Properties.VariableNames{'Cimbi'} = 'CIMBI_ID';

% Convert variables to appropriate types
age_sex.CIMBI_ID = categorical(age_sex.CIMBI_ID);
master_tbl_p.CIMBI_ID = categorical(master_tbl_p.CIMBI_ID); 
master_tbl_n.CIMBI_ID = categorical(master_tbl_n.CIMBI_ID);

% Join tables
master_tbl_p = join(master_tbl_p, age_sex, 'Keys', 'CIMBI_ID');
master_tbl_n = join(master_tbl_n, age_sex, 'Keys', 'CIMBI_ID');

master_tbl_p.Age = double(master_tbl_p.Age);
master_tbl_n.Age = double(master_tbl_n.Age);
master_tbl_p.Sex = categorical(master_tbl_p.Sex);
master_tbl_n.Sex = categorical(master_tbl_n.Sex);

out_dir = fullfile(project_dir, 'outputs');

% Create date-specific output directory
timestamp = datestr(now, 'yyyy-mm-dd');
out_dir = fullfile(out_dir, timestamp);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
    fprintf('Created output directory: %s\n', out_dir);
end

%% Estimate autocorrelations

[mean_autocorr_p, se_autocorr_p, time_p] = calculate_autocorrelations(...
    tbl_scanner_p, 'ts_SchaeferTian232', 150, 2.0);
plot_autocorrelations(mean_autocorr_p, se_autocorr_p, time_p, 2.0);

[mean_autocorr_n, se_autocorr_n, time_n] = calculate_autocorrelations(...
    tbl_scanner_n, 'ts_SchaeferTian232', 300, 0.8);
plot_autocorrelations(mean_autocorr_n, se_autocorr_n, time_n, 0.8);

%% Calculate directionalities

% Scanner P - tau 3 (6 seconds)
global_scores = zeros(height(tbl_scanner_p), 1); %INSIDEOUT
master_tbl_p.FC_matrix = cell(height(tbl_scanner_p), 1);
master_tbl_p.FCtf = cell(height(tbl_scanner_p), 1);
master_tbl_p.MI_matrix = cell(height(tbl_scanner_p), 1);
master_tbl_p.forward_MI = cell(height(tbl_scanner_p), 1);
master_tbl_p.reverse_MI = cell(height(tbl_scanner_p), 1);
master_tbl_p.insideout_regional_dirscore = cell(height(tbl_scanner_p), 1); %Regional INSIDEOUT
master_tbl_p.directionality_scores = cell(height(tbl_scanner_p), 1); %NEW: Directionality scores

for i = 1:height(tbl_scanner_p)
    [directionality_scores, FCtf, MI_matrix, FC_matrix, global_score, forward_MI, reverse_MI, insideout_regional_dirscore] = compute_insideout(tbl_scanner_p.ts_SchaeferTian232{i}, 3, 100);
    
    % Store results
    global_scores(i) = global_score; %INSIDEOUT
    master_tbl_p.directionality_scores{i} = directionality_scores; %NEW
    master_tbl_p.forward_MI{i} = forward_MI;
    master_tbl_p.reverse_MI{i} = reverse_MI;
    master_tbl_p.insideout_regional_dirscore{i} = insideout_regional_dirscore; %Regional INSIDEOUT
    master_tbl_p.FC_matrix{i} = FC_matrix;
    master_tbl_p.MI_matrix{i} = MI_matrix;
    master_tbl_p.FCtf{i} = FCtf;
    if mod(i, 5) == 0
        fprintf('Processed %d/%d participants\n', i, height(tbl_scanner_p));
    end
end
master_tbl_p.insideout = global_scores;

% Scanner N - tau 8 (6.4 seconds)  
global_scores = zeros(height(tbl_scanner_n), 1); %INSIDEOUT
master_tbl_n.FC_matrix = cell(height(tbl_scanner_n), 1);
master_tbl_n.FCtf = cell(height(tbl_scanner_n), 1);
master_tbl_n.MI_matrix = cell(height(tbl_scanner_n), 1);
master_tbl_n.forward_MI = cell(height(tbl_scanner_n), 1);
master_tbl_n.reverse_MI = cell(height(tbl_scanner_n), 1);
master_tbl_n.insideout_regional_dirscore = cell(height(tbl_scanner_n), 1); %Regional INSIDEOUT
master_tbl_n.directionality_scores = cell(height(tbl_scanner_n), 1); %NEW: Directionality scores

for i = 1:height(tbl_scanner_n)
    [directionality_scores, FCtf, MI_matrix, FC_matrix, global_score, forward_MI, reverse_MI, insideout_regional_dirscore] = compute_insideout(tbl_scanner_n.ts_SchaeferTian232{i}, 8, 100);
    
    % Store results
    global_scores(i) = global_score; %INSIDEOUT
    master_tbl_n.directionality_scores{i} = directionality_scores; %NEW
    master_tbl_n.forward_MI{i} = forward_MI;
    master_tbl_n.reverse_MI{i} = reverse_MI;
    master_tbl_n.insideout_regional_dirscore{i} = insideout_regional_dirscore; %Regional INSIDEOUT
    master_tbl_n.FC_matrix{i} = FC_matrix;
    master_tbl_n.MI_matrix{i} = MI_matrix;
    master_tbl_n.FCtf{i} = FCtf;
    if mod(i, 5) == 0
        fprintf('Processed %d/%d participants\n', i, height(tbl_scanner_n));
    end
end
master_tbl_n.insideout = global_scores;

timestamp = datestr(now, 'yyyy-mm-dd');
save([out_dir, '/master_tbl_p_', timestamp, '.mat'], 'master_tbl_p');
save([out_dir, '/master_tbl_n_', timestamp, '.mat'], 'master_tbl_n');

% Combine tables and add scanner identifier
master_tbl_p.scanner = repmat({'P'}, height(master_tbl_p), 1);
master_tbl_n.scanner = repmat({'N'}, height(master_tbl_n), 1);
master_tbl_combined = [master_tbl_p; master_tbl_n];
save([out_dir, '/master_tbl_combined_', timestamp, '.mat'], 'master_tbl_combined');
save([out_dir, '/master_tbl_combined.mat'], 'master_tbl_combined');

%% Quality Control - Check for NaN and Inf values in connectivity matrices
fprintf('Running quality control checks on connectivity matrices...\n\n');

% Define the matrix columns to check
matrix_columns = {'FC_matrix', 'FCtf', 'MI_matrix', 'forward_MI'};
total_subjects = height(master_tbl_combined);

% Initialize tracking variables
subjects_with_issues = [];
issue_summary = struct();

for col_idx = 1:length(matrix_columns)
    current_col = matrix_columns{col_idx};
    fprintf('Checking %s...\n', current_col);
    
    subjects_with_issues_this_col = [];
    
    for i = 1:total_subjects
        matrix = master_tbl_combined.(current_col){i};
        
        num_nans = sum(isnan(matrix(:)));
        num_infs = sum(isinf(matrix(:)));
        total_elements = numel(matrix);
        
        if num_nans > 0 || num_infs > 0
            subjects_with_issues_this_col = [subjects_with_issues_this_col, i];
            subjects_with_issues = [subjects_with_issues, i];
            
            percent_nans = (num_nans / total_elements) * 100;
            percent_infs = (num_infs / total_elements) * 100;
            
            fprintf('  Subject %d: %d NaNs (%.2f%%), %d Infs (%.2f%%)\n', ...
                i, num_nans, percent_nans, num_infs, percent_infs);
            
            if num_infs > 0
                fprintf('    Matrix range: Min=%.4f, Max=%s\n', ...
                    min(matrix(~isinf(matrix))), num2str(max(matrix(:))));
            end
        end
    end
    
    % Store summary for this column
    issue_summary.(current_col).subjects_with_issues = subjects_with_issues_this_col;
    issue_summary.(current_col).num_subjects_affected = length(subjects_with_issues_this_col);
    
    if isempty(subjects_with_issues_this_col)
        fprintf('  ✓ No NaN or Inf values found in %s\n', current_col);
    else
        fprintf('  ⚠ Found NaN/Inf values in %d subjects for %s\n', length(subjects_with_issues_this_col), current_col);
    end
    fprintf('\n');
end

% Overall summary
subjects_with_issues = unique(subjects_with_issues);
fprintf('=== SUMMARY ===\n');
fprintf('Total subjects checked: %d\n', total_subjects);
fprintf('Subjects with NaN/Inf values in at least one matrix: %d\n', length(subjects_with_issues));

if ~isempty(subjects_with_issues)
    fprintf('Subject IDs with issues: %s\n', mat2str(subjects_with_issues));
    fprintf('\nRecommendation: Handle NaN/Inf values before modularity analysis.\n');
else
    fprintf('✓ All matrices are clean - no NaN or Inf values detected!\n');
end
fprintf('\n');

%% Calculate modularity for all matrices
% res = 1.0;
% nreps = 100;
% ncores = 32;
% 
% % Initialize arrays for all 4 modularity metrics
% master_tbl_combined.modularity_fc = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_fc_abs = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_forward_corr = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_forward_corr_abs = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_MI = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_MI_abs = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_forward_MI = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.modularity_forward_MI_abs = zeros(height(master_tbl_combined), 1);
% 
% % Initialize community count columns (mean values across repetitions)
% master_tbl_combined.communities_fc_actual = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_fc_null = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_forward_corr_actual = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_forward_corr_null = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_MI_actual = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_MI_null = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_forward_MI_actual = zeros(height(master_tbl_combined), 1);
% master_tbl_combined.communities_forward_MI_null = zeros(height(master_tbl_combined), 1);
% 
% % Initialize community data storage for histograms
% all_actual_comm_fc = [];
% all_null_comm_fc = [];
% all_actual_comm_fcorr = [];
% all_null_comm_fcorr = [];
% all_actual_comm_mi = [];
% all_null_comm_mi = [];
% all_actual_comm_fmi = [];
% all_null_comm_fmi = [];
% 
% % Progress tracking
% total_subjects = height(master_tbl_combined);
% start_time = tic;
% fprintf('Starting modularity calculations for %d scans ...\n', total_subjects);
% 
% for i = 1:height(master_tbl_combined)
% 
%     if i > 1
%         elapsed_time = toc(start_time);
%         estimated_total_time = elapsed_time / (i-1) * total_subjects;
%         time_remaining = estimated_total_time - elapsed_time;
% 
%         fprintf('\rProcessing scan %d/%d - Est. time remaining: %.1f min', ...
%             i-1, total_subjects, time_remaining/60);
%     else
%         fprintf('\rProcessing scan %d/%d - Calculating...', i, total_subjects);
%     end
% 
%     % Static FC matrix (undirected)
%     [norm_mod, abs_mod, actual_comm_fc, null_comm_fc] = calculate_modularity_metrics3(master_tbl_combined.FC_matrix{i}, nreps, res, ncores, false);
%     master_tbl_combined.modularity_fc(i) = norm_mod;
%     master_tbl_combined.modularity_fc_abs(i) = abs_mod;
%     master_tbl_combined.communities_fc_actual(i) = mean(actual_comm_fc);
%     master_tbl_combined.communities_fc_null(i) = mean(null_comm_fc);
%     all_actual_comm_fc = [all_actual_comm_fc, actual_comm_fc];
%     all_null_comm_fc = [all_null_comm_fc, null_comm_fc];
% 
%     % Forward correlation matrix (directed)
%     [norm_mod, abs_mod, actual_comm_fcorr, null_comm_fcorr] = calculate_modularity_metrics3(master_tbl_combined.FCtf{i}, nreps, res, ncores, true);
%     master_tbl_combined.modularity_forward_corr(i) = norm_mod;
%     master_tbl_combined.modularity_forward_corr_abs(i) = abs_mod;
%     master_tbl_combined.communities_forward_corr_actual(i) = mean(actual_comm_fcorr);
%     master_tbl_combined.communities_forward_corr_null(i) = mean(null_comm_fcorr);
%     all_actual_comm_fcorr = [all_actual_comm_fcorr, actual_comm_fcorr];
%     all_null_comm_fcorr = [all_null_comm_fcorr, null_comm_fcorr];
% 
%     % Static MI matrix (undirected)
%     [norm_mod, abs_mod, actual_comm_mi, null_comm_mi] = calculate_modularity_metrics3(master_tbl_combined.MI_matrix{i}, nreps, res, ncores, false);
%     master_tbl_combined.modularity_MI(i) = norm_mod;
%     master_tbl_combined.modularity_MI_abs(i) = abs_mod;
%     master_tbl_combined.communities_MI_actual(i) = mean(actual_comm_mi);
%     master_tbl_combined.communities_MI_null(i) = mean(null_comm_mi);
%     all_actual_comm_mi = [all_actual_comm_mi, actual_comm_mi];
%     all_null_comm_mi = [all_null_comm_mi, null_comm_mi];
% 
%     % Forward MI matrix (directed)
%     [norm_mod, abs_mod, actual_comm_fmi, null_comm_fmi] = calculate_modularity_metrics3(master_tbl_combined.forward_MI{i}, nreps, res, ncores, true);
%     master_tbl_combined.modularity_forward_MI(i) = norm_mod;
%     master_tbl_combined.modularity_forward_MI_abs(i) = abs_mod;
%     master_tbl_combined.communities_forward_MI_actual(i) = mean(actual_comm_fmi);
%     master_tbl_combined.communities_forward_MI_null(i) = mean(null_comm_fmi);
%     all_actual_comm_fmi = [all_actual_comm_fmi, actual_comm_fmi];
%     all_null_comm_fmi = [all_null_comm_fmi, null_comm_fmi];
% end
% 
% % Final progress update
% fprintf('\rProgress: [====================] 100%% (%d/%d) - Complete! Total time: %.1f min\n', ...
%     total_subjects, total_subjects, toc(start_time)/60);
% 
% save([out_dir, '/master_tbl_combined.mat'], 'master_tbl_combined');

%% Remove some old cols from master_tbl_combined

% Delete specified columns from master_tbl_combined
columns_to_delete = {'insideout_regional_FCtf', 'insideout_regional_diff', ...
                     'insideout_regional_dirscore_old', 'comparison', ...
                     'insideout_old', 'asymmetry'};

% Check which columns actually exist before trying to delete them
existing_vars = master_tbl_combined.Properties.VariableNames;
columns_to_delete_existing = intersect(columns_to_delete, existing_vars);

if ~isempty(columns_to_delete_existing)
    fprintf('Deleting columns: %s\n', strjoin(columns_to_delete_existing, ', '));
    master_tbl_combined = removevars(master_tbl_combined, columns_to_delete_existing);
    fprintf('Deleted %d columns. New table has %d columns.\n', ...
            length(columns_to_delete_existing), width(master_tbl_combined));
else
    fprintf('No specified columns found to delete.\n');
end

% Show remaining column names
fprintf('\nRemaining columns:\n');
disp(master_tbl_combined.Properties.VariableNames');

%% Calculate Global and Regional GCOR and MI-based GCOR for combined table

fprintf('Calculating Global and Regional GCOR and MI-based GCOR...\n');

num_regions = height(master_tbl_combined.FC_matrix{1});

master_tbl_combined.global_GCOR = zeros(height(master_tbl_combined), 1);
master_tbl_combined.global_MI_GCOR = zeros(height(master_tbl_combined), 1);
master_tbl_combined.regional_GCOR = cell(height(master_tbl_combined), 1);
master_tbl_combined.regional_MI_GCOR = cell(height(master_tbl_combined), 1);

for i = 1:height(master_tbl_combined)
    if mod(i, 20) == 0
        fprintf('Processing subject %d/%d\n', i, height(master_tbl_combined));
    end
    
    FC_mat = master_tbl_combined.FC_matrix{i};
    
    % Calculate standard GCOR (mean of upper triangle)
    upper_tri = triu(FC_mat, 1);
    % Set diagonal to 0 
    upper_tri(logical(eye(size(upper_tri)))) = 0;
    % mean of non-zero elements
    master_tbl_combined.global_GCOR(i) = mean(upper_tri(upper_tri ~= 0));
    
    MI_mat = master_tbl_combined.MI_matrix{i};
    % Set diagonal to 0 to ignore it later
    MI_mat(logical(eye(size(MI_mat)))) = 0;
    
    % Calculate MI-based GCOR (mean of upper triangle)
    upper_tri_MI = triu(MI_mat, 1);
    master_tbl_combined.global_MI_GCOR(i) = mean(upper_tri_MI(upper_tri_MI ~= 0));
    
    % Calculate regional GCOR (mean connectivity per region, excluding diagonal)
    regional_GCOR = zeros(num_regions, 1);
    regional_MI_GCOR = zeros(num_regions, 1);
    for reg = 1:num_regions
        regional_GCOR(reg) = mean(FC_mat(reg, [1:reg-1, reg+1:end]));
        regional_MI_GCOR(reg) = mean(MI_mat(reg, [1:reg-1, reg+1:end]));
    end
    
    % Store as 232x1 vectors
    master_tbl_combined.regional_GCOR{i} = regional_GCOR;
    master_tbl_combined.regional_MI_GCOR{i} = regional_MI_GCOR;
end

fprintf('Global and Regional GCOR/MI-GCOR calculation complete.\n');

timestamp = datestr(now, 'yyyy-mm-dd');
save([out_dir, '/master_tbl_combined_', timestamp, '.mat'], 'master_tbl_combined');
save([out_dir, '/master_tbl_combined.mat'], 'master_tbl_combined');

%% Load pre-computed modularity scores
fprintf('\nLoading pre-computed modularity scores from July 13, 2025...\n');

% Load the scanner-specific tables with modularity
load(fullfile(project_dir, 'outputs/master_tbl_p_2025-07-13_post_modularity.mat'));
load(fullfile(project_dir, 'outputs/master_tbl_n_2025-07-13_post_modularity.mat'));

% Store the loaded tables temporarily
master_tbl_p_mod = master_tbl_p;
master_tbl_n_mod = master_tbl_n;

% Reload the current combined table (to ensure we have it fresh)
load([out_dir, '/master_tbl_combined.mat']);

% Extract only the modularity columns we need
modularity_cols = {'modularity_fc', 'modularity_forward_MI'};

% Find indices for P and N scanners in current combined table
p_idx = strcmp(master_tbl_combined.scanner, 'P');
n_idx = strcmp(master_tbl_combined.scanner, 'N');

% Copy modularity values
for col = modularity_cols
    col_name = col{1};
    master_tbl_combined.(col_name)(p_idx) = master_tbl_p_mod.(col_name);
    master_tbl_combined.(col_name)(n_idx) = master_tbl_n_mod.(col_name);
end

fprintf('Modularity scores successfully loaded and merged.\n');
fprintf('  - Scanner P: %d scans\n', sum(p_idx));
fprintf('  - Scanner N: %d scans\n', sum(n_idx));

% Save the updated combined table with modularity
save([out_dir, '/master_tbl_combined_', timestamp, '.mat'], 'master_tbl_combined');
save([out_dir, '/master_tbl_combined.mat'], 'master_tbl_combined');
fprintf('Updated master_tbl_combined saved with modularity scores.\n');


%% Save out and export to R for stats and plots
% Create export table from the combined table
export_tbl = table(...
    master_tbl_combined.CIMBI_ID, ...
    master_tbl_combined.scanner, ...
    master_tbl_combined.psi_conc, ...
    master_tbl_combined.time_since_admin, ...
    master_tbl_combined.SDI_score, ...
    master_tbl_combined.Age, ...
    master_tbl_combined.Sex, ...
    master_tbl_combined.insideout, ...
    master_tbl_combined.modularity_fc, ...
    master_tbl_combined.modularity_forward_MI, ...
    master_tbl_combined.global_GCOR, ...
    'VariableNames', {...
        'CIMBI_ID', ...
        'scanner', ...
        'psi_conc', ...
        'time_since_admin', ...
        'SDI_score', ...
        'Age', ...
        'Sex', ...
        'insideout', ...
        'modularity_fc', ...
        'modularity_forward_MI', ...
        'global_GCOR'});
        
% Save to Excel
writetable(export_tbl, [out_dir, '/combined_scanner_metrics.xlsx']);

%% Export subject-level regional scores for R analysis

fprintf('\n=== EXPORTING SUBJECT-LEVEL REGIONAL SCORES FOR R ===\n');

n_subjects = height(master_tbl_combined);
n_regions = 232;

subject_ids = master_tbl_combined.CIMBI_ID;
psi_conc = master_tbl_combined.psi_conc;
time_since_admin = master_tbl_combined.time_since_admin;

% Directionality scores
all_dirscore = zeros(n_subjects, n_regions);
for subj = 1:n_subjects
    all_dirscore(subj, :) = master_tbl_combined.insideout_regional_dirscore{subj};
end

dirscore_table = array2table(all_dirscore);
dirscore_table.CIMBI_ID = subject_ids;
dirscore_table.psi_conc = psi_conc;
dirscore_table.time_since_admin = time_since_admin;
dirscore_table = [dirscore_table(:,end-2:end), dirscore_table(:,1:end-3)];

for reg = 1:n_regions
    dirscore_table.Properties.VariableNames{3+reg} = sprintf('Region_%d', reg);
end

writetable(dirscore_table, [out_dir, '/subject_directionality_scores.csv']);

% Regional GCOR scores
all_gcor = zeros(n_subjects, n_regions);
for subj = 1:n_subjects
    all_gcor(subj, :) = master_tbl_combined.regional_GCOR{subj};
end

gcor_table = array2table(all_gcor);
gcor_table.CIMBI_ID = subject_ids;
gcor_table.psi_conc = psi_conc;
gcor_table.time_since_admin = time_since_admin;
gcor_table = [gcor_table(:,end-2:end), gcor_table(:,1:end-3)];

for reg = 1:n_regions
    gcor_table.Properties.VariableNames{3+reg} = sprintf('Region_%d', reg);
end

writetable(gcor_table, [out_dir, '/subject_regional_GCOR_scores.csv']);

% Forward MI scores
all_forward_scores = zeros(n_subjects, n_regions);
for subj = 1:n_subjects
    forward_matrix = master_tbl_combined.forward_MI{subj};
    for reg = 1:n_regions
        all_forward_scores(subj, reg) = mean(forward_matrix(reg,:));
    end
end

forward_table = array2table(all_forward_scores);
forward_table.CIMBI_ID = subject_ids;
forward_table.psi_conc = psi_conc;
forward_table.time_since_admin = time_since_admin;
forward_table = [forward_table(:,end-2:end), forward_table(:,1:end-3)];

for reg = 1:n_regions
    forward_table.Properties.VariableNames{3+reg} = sprintf('Region_%d', reg);
end

writetable(forward_table, [out_dir, '/subject_forward_MI_scores.csv']);

% Reverse MI scores
all_reverse_scores = zeros(n_subjects, n_regions);
for subj = 1:n_subjects
    forward_matrix = master_tbl_combined.forward_MI{subj};
    for reg = 1:n_regions
        all_reverse_scores(subj, reg) = mean(forward_matrix(:,reg));
    end
end

reverse_table = array2table(all_reverse_scores);
reverse_table.CIMBI_ID = subject_ids;
reverse_table.psi_conc = psi_conc;
reverse_table.time_since_admin = time_since_admin;
reverse_table = [reverse_table(:,end-2:end), reverse_table(:,1:end-3)];

for reg = 1:n_regions
    reverse_table.Properties.VariableNames{3+reg} = sprintf('Region_%d', reg);
end

writetable(reverse_table, [out_dir, '/subject_reverse_MI_scores.csv']);

fprintf('Exported subject-level regional scores:\n');
fprintf('- subject_directionality_scores.csv (%d subjects x %d regions)\n', n_subjects, n_regions);
fprintf('- subject_regional_GCOR_scores.csv (%d subjects x %d regions)\n', n_subjects, n_regions);
fprintf('- subject_forward_MI_scores.csv (%d subjects x %d regions)\n', n_subjects, n_regions);
fprintf('- subject_reverse_MI_scores.csv (%d subjects x %d regions)\n', n_subjects, n_regions);

