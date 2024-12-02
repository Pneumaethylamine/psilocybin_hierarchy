cd('/indirect/staff/drummondmcculloch/UPF_project');
addpath(genpath('/indirect/staff/drummondmcculloch/UPF_project'));
addpath("/mrdata/proj-doccls/code")

load('/indirect/staff/drummondmcculloch/UPF_project/data/data_09-Mar-2023_SchaeferTian232_ts.mat');
age_sex = readtable('/mrdata/np2/p3/entropy/critical_files/CIMBI_age_sex_clean.xlsx');

tbl = SDI_PPL_rsfmri_tbl(strcmpi(SDI_PPL_rsfmri_tbl.include, 'yes'), :);
%remove PA scans
tbl = tbl(~contains(lower(tbl.rsfmri_location), 'pa'), :);

%% First split the data by scanner
scanner_p_idx = startsWith(tbl.rsfmri_raw_ID, 'p');  % 2.0s TR scans
scanner_n_idx = startsWith(tbl.rsfmri_raw_ID, 'n');  % 0.8s TR scans

% Create separate tables
tbl_scanner_p = tbl(scanner_p_idx, :);
tbl_scanner_n = tbl(scanner_n_idx, :);

%% Create master tables with required columns
master_tbl_p = table(tbl_scanner_p.CIMBI_ID, tbl_scanner_p.psi_conc, ...
    tbl_scanner_p.time_since_admin, tbl_scanner_p.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

master_tbl_n = table(tbl_scanner_n.CIMBI_ID, tbl_scanner_n.psi_conc, ...
    tbl_scanner_n.time_since_admin, tbl_scanner_n.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

%% Calculate INSIDEOUT for different taus
% For Scanner P
taus_scanner_p = [1, 2, 3, 4];  % will represent 2, 4, 6, 8 seconds
for tau = taus_scanner_p
    global_scores = zeros(height(tbl_scanner_p), 1);
    for i = 1:height(tbl_scanner_p)
        [~, global_score, ~, ~, asymmetry] = computeDirectionality(tbl_scanner_p.ts_SchaeferTian232{i}, tau);
        global_scores(i) = global_score;
        asymmetries(i) = asymmetry;
    end
    master_tbl_p.(['insideout_tau_' num2str(tau)]) = global_scores;
    master_tbl_p.(['asymmetry_tau_' num2str(tau)]) = asymmetries';
end

clear global_scores asymmetries

% For Scanner N
taus_scanner_n = [2, 4, 6, 8, 10];  % will represent 1.6, 3.2, 4.8, 6.4, 8 seconds
for tau = taus_scanner_n
    global_scores = zeros(height(tbl_scanner_n), 1);
    for i = 1:height(tbl_scanner_n)
        [~, global_score, ~, ~, asymmetry] = computeDirectionality(tbl_scanner_n.ts_SchaeferTian232{i}, tau);
        global_scores(i) = global_score;
        asymmetries(i) = asymmetry;
    end
    master_tbl_n.(['insideout_tau_' num2str(tau)]) = global_scores;
    master_tbl_n.(['asymmetry_tau_' num2str(tau)]) = asymmetries';
end


% Set up plot parameters
maxlag_p = 10; % For scanner P (adjust based on your needs)
maxlag_n = 25; % For scanner N (adjust based on your needs)
TR_p = 2.0; % TR for scanner P in seconds
TR_n = 0.8; % TR for scanner N in seconds

% Initialize matrices to store autocorrelations
autocorr_p = zeros(height(tbl_scanner_p), maxlag_p + 1);
autocorr_n = zeros(height(tbl_scanner_n), maxlag_n + 1);

% Calculate autocorrelation for each scan in scanner P
for i = 1:height(tbl_scanner_p)
    ts = tbl_scanner_p.ts_SchaeferTian232{i};
    % Average across ROIs
    mean_ts = mean(ts, 2);
    autocorr_p(i, :) = autocorr_manual(mean_ts, maxlag_p);
end

% Calculate autocorrelation for each scan in scanner N
for i = 1:height(tbl_scanner_n)
    ts = tbl_scanner_n.ts_SchaeferTian232{i};
    % Average across ROIs
    mean_ts = mean(ts, 2);
    autocorr_n(i, :) = autocorr_manual(mean_ts, maxlag_n);
end

% Calculate mean and standard error
mean_autocorr_p = mean(autocorr_p, 1);
se_autocorr_p = std(autocorr_p, 0, 1) / sqrt(size(autocorr_p, 1));
mean_autocorr_n = mean(autocorr_n, 1);
se_autocorr_n = std(autocorr_n, 0, 1) / sqrt(size(autocorr_n, 1));

% Create time vectors
time_p = (0:maxlag_p) * TR_p;
time_n = (0:maxlag_n) * TR_n;

% Create the plot
figure('Position', [100 100 800 400]);

% Scanner P
subplot(1,2,1)
errorbar(time_p, mean_autocorr_p, se_autocorr_p, 'b-', 'LineWidth', 2);
xlabel('Time lag (seconds)');
ylabel('Autocorrelation');
title('Scanner P (TR = 2.0s)');
grid on;
ylim([-0.2 1]);
xlim([0 max(time_p)]);

% Scanner N
subplot(1,2,2)
errorbar(time_n, mean_autocorr_n, se_autocorr_n, 'r-', 'LineWidth', 2);
xlabel('Time lag (seconds)');
ylabel('Autocorrelation');
title('Scanner N (TR = 0.8s)');
grid on;
ylim([-0.2 1]);
xlim([0 max(time_n)]);

title('Average Autocorrelation Functions by Scanner', 'Position', [0.5 1.05]);

%% Add age and sex 

% Select only needed columns
age_sex = age_sex(:, {'Cimbi', 'Age', 'Sex'});

% Add missing participant
missing_dude = table(55992, 33.0, "Male", 'VariableNames', {'Cimbi', 'Age', 'Sex'});
age_sex = [age_sex; missing_dude];

% Sort by CIMBI ID
age_sex = sortrows(age_sex, 'Cimbi');

% Rename Cimbi to CIMBI_ID to match your existing tables
age_sex.Properties.VariableNames{'Cimbi'} = 'CIMBI_ID';

% Convert CIMBI_ID to categorical
age_sex.CIMBI_ID = categorical(age_sex.CIMBI_ID);
master_tbl_p.CIMBI_ID = categorical(master_tbl_p.CIMBI_ID);
master_tbl_n.CIMBI_ID = categorical(master_tbl_n.CIMBI_ID);

% Use leftjoin instead of outerjoin
master_tbl_p = join(master_tbl_p, age_sex, 'Keys', 'CIMBI_ID');
master_tbl_n = join(master_tbl_n, age_sex, 'Keys', 'CIMBI_ID');

% Convert variables to appropriate types
master_tbl_p.Age = double(master_tbl_p.Age);
master_tbl_n.Age = double(master_tbl_n.Age);
master_tbl_p.Sex = categorical(master_tbl_p.Sex);
master_tbl_n.Sex = categorical(master_tbl_n.Sex);

%% Stats - Setup results storage
% Define measures to analyze
insideout_measures_p = {'insideout_tau_1', 'insideout_tau_2', 'insideout_tau_3', 'insideout_tau_4'};
insideout_measures_n = {'insideout_tau_2', 'insideout_tau_4', 'insideout_tau_6', 'insideout_tau_8', 'insideout_tau_10'};

asymmetry_measures_p = {'asymmetry_tau_1', 'asymmetry_tau_2', 'asymmetry_tau_3', 'asymmetry_tau_4'};
asymmetry_measures_n = {'asymmetry_tau_2', 'asymmetry_tau_4', 'asymmetry_tau_6', 'asymmetry_tau_8', 'asymmetry_tau_10'};


% Create results tables structure
results_p = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});
results_n = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});

%% Run models for Scanner P
fprintf('\nProcessing Scanner P models:\n');

% INSIDEOUT measures
for i = 1:length(insideout_measures_p)
    measure = insideout_measures_p{i};
    fprintf('Running model for %s\n', measure);
    
    % Fit model
    lme_model = fitlme(master_tbl_p, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
    
    % Extract results for psi_conc effect
    coef_table = lme_model.Coefficients;
    psi_idx = strcmp(coef_table.Name, 'psi_conc');
    
    % Store results
    new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
        coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
    results_p = [results_p; new_row];
end

for i = 1:length(asymmetry_measures_p)
    measure = asymmetry_measures_p{i};
    fprintf('Running model for %s\n', measure);
    
    % Fit model
    lme_model = fitlme(master_tbl_p, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
    
    % Extract results for psi_conc effect
    coef_table = lme_model.Coefficients;
    psi_idx = strcmp(coef_table.Name, 'psi_conc');
    
    % Store results
    new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
        coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
    results_p = [results_p; new_row];
end


%% Run models for Scanner N
fprintf('\nProcessing Scanner N models:\n');

% INSIDEOUT measures
for i = 1:length(insideout_measures_n)
    measure = insideout_measures_n{i};
    fprintf('Running model for %s\n', measure);
    
    % Fit model
    lme_model = fitlme(master_tbl_n, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
    
    % Extract results for psi_conc effect
    coef_table = lme_model.Coefficients;
    psi_idx = strcmp(coef_table.Name, 'psi_conc');
    
    % Store results
    new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
        coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
    results_n = [results_n; new_row];
end

for i = 1:length(asymmetry_measures_n)
    measure = asymmetry_measures_n{i};
    fprintf('Running model for %s\n', measure);
    
    % Fit model
    lme_model = fitlme(master_tbl_n, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
    
    % Extract results for psi_conc effect
    coef_table = lme_model.Coefficients;
    psi_idx = strcmp(coef_table.Name, 'psi_conc');
    
    % Store results
    new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
        coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
    results_n = [results_n; new_row];
end

% Save results
save('lme_results.mat', 'results_p', 'results_n');

% Display results
disp('Scanner P Results:');
disp(results_p);
disp('Scanner N Results:');
disp(results_n);

%% Partial residual plots

figure('Position', [100 100 1500 800]);

% Get overall min/max for axes
all_resid_p = [];
all_resid_n = [];
for i = 1:4
   mdl_p = fitlme(master_tbl_p, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', insideout_measures_p{i}));
   resid_p = residuals(mdl_p, 'Conditional', true);
   partial_resid_p = resid_p + mdl_p.Coefficients.Estimate(2)*master_tbl_p.psi_conc;
   all_resid_p = [all_resid_p; partial_resid_p];
end

for i = 1:5
   mdl_n = fitlme(master_tbl_n, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', insideout_measures_n{i}));
   resid_n = residuals(mdl_n, 'Conditional', true);
   partial_resid_n = resid_n + mdl_n.Coefficients.Estimate(2)*master_tbl_n.psi_conc;
   all_resid_n = [all_resid_n; partial_resid_n];
end

y_min = min([all_resid_p; all_resid_n]);
y_max = max([all_resid_p; all_resid_n]);
x_min = min([master_tbl_p.psi_conc; master_tbl_n.psi_conc]);
x_max = max([master_tbl_p.psi_conc; master_tbl_n.psi_conc]);

% Scanner P plots
for i = 1:4
   subplot(2,5,i)
   measure = insideout_measures_p{i};
   mdl = fitlme(master_tbl_p, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure));
   
   resid = residuals(mdl, 'Conditional', true);
   partial_resid = resid + mdl.Coefficients.Estimate(2)*master_tbl_p.psi_conc;
   
   scatter(master_tbl_p.psi_conc, partial_resid, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
   hold on
   
   valid_idx = ~isnan(partial_resid) & ~isnan(master_tbl_p.psi_conc);
   coef = polyfit(master_tbl_p.psi_conc(valid_idx), partial_resid(valid_idx), 1);
   x = linspace(min(master_tbl_p.psi_conc), max(master_tbl_p.psi_conc), 100);
   y = polyval(coef, x);
   plot(x, y, 'r-', 'LineWidth', 2)
   
   xlabel('PSI Concentration')
   ylabel('Partial Residuals')
   title(sprintf('Scanner P - Tau %ds', i*2))
   grid on
   xlim([x_min x_max])
   ylim([y_min y_max])
end

% Scanner N plots
for i = 1:5
   subplot(2,5,i+5)
   measure = insideout_measures_n{i};
   mdl = fitlme(master_tbl_n, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure));
   
   resid = residuals(mdl, 'Conditional', true);
   partial_resid = resid + mdl.Coefficients.Estimate(2)*master_tbl_n.psi_conc;
   
   scatter(master_tbl_n.psi_conc, partial_resid, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
   hold on
   
   valid_idx = ~isnan(partial_resid) & ~isnan(master_tbl_n.psi_conc);
   coef = polyfit(master_tbl_n.psi_conc(valid_idx), partial_resid(valid_idx), 1);
   x = linspace(min(master_tbl_n.psi_conc), max(master_tbl_n.psi_conc), 100);
   y = polyval(coef, x);
   plot(x, y, 'r-', 'LineWidth', 2)
   
   xlabel('PSI Concentration')
   ylabel('Partial Residuals')
   title(sprintf('Scanner N - Tau %.1fs', i*1.6))
   grid on
   xlim([x_min x_max])
   ylim([y_min y_max])
end

subplot(2,5,5)
axis off

savefig('partiout_dir0al_residual_plots.fig')
print('partial_residual_plots','-dpng','-r300')

%% 

figure('Position', [100 100 1500 800]);

% Create colormap with distinct colors for participants - using HSV instead of parula
colors_p = hsv(length(unique(master_tbl_p.CIMBI_ID)));
colors_n = hsv(length(unique(master_tbl_n.CIMBI_ID)));

% First calculate all deltas to find overall min/max
all_deltas = [];

% Scanner P
for i = 1:length(insideout_measures_p)
    measure = insideout_measures_p{i};
    subjects = unique(master_tbl_p.CIMBI_ID);
    for subj = subjects'
        subj_data = master_tbl_p(master_tbl_p.CIMBI_ID == subj, :);
        baseline_idx = subj_data.time_since_admin < 3;
        baseline_value = subj_data.(measure)(baseline_idx);
        all_deltas = [all_deltas; subj_data.(measure) - baseline_value];
    end
end

% Scanner N
for i = 1:length(insideout_measures_n)
    measure = insideout_measures_n{i};
    subjects = unique(master_tbl_n.CIMBI_ID);
    for subj = subjects'
        subj_data = master_tbl_n(master_tbl_n.CIMBI_ID == subj, :);
        baseline_idx = subj_data.time_since_admin < 3;
        baseline_value = subj_data.(measure)(baseline_idx);
        all_deltas = [all_deltas; subj_data.(measure) - baseline_value];
    end
end

y_min = min(all_deltas);
y_max = max(all_deltas);

% Scanner P plots
for i = 1:length(insideout_measures_p)
    measure = insideout_measures_p{i};
    subjects = unique(master_tbl_p.CIMBI_ID);
    
    subplot(2,5,i)
    hold on
    
    for subj_idx = 1:length(subjects)
        subj = subjects(subj_idx);
        subj_data = master_tbl_p(master_tbl_p.CIMBI_ID == subj, :);
        baseline_idx = subj_data.time_since_admin < 3;
        baseline_value = subj_data.(measure)(baseline_idx);
        
        deltas = subj_data.(measure) - baseline_value;
        times = subj_data.time_since_admin;
        
        % Sort by time to ensure line connects points in temporal order
        [sorted_times, sort_idx] = sort(times);
        sorted_deltas = deltas(sort_idx);
        
        % Plot line and points
        plot(sorted_times, sorted_deltas, '-', 'Color', [colors_p(subj_idx,:), 0.3])
        scatter(sorted_times, sorted_deltas, 50, colors_p(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
    end
    
    xlabel('Time since administration (min)')
    ylabel('Δ from baseline')
    title(sprintf('Scanner P - Tau %ds', i*2))
    ylim([y_min y_max])
    grid on
    hold off
end

% Scanner N plots
for i = 1:length(insideout_measures_n)
    measure = insideout_measures_n{i};
    subjects = unique(master_tbl_n.CIMBI_ID);
    
    subplot(2,5,i+5)
    hold on
    
    for subj_idx = 1:length(subjects)
        subj = subjects(subj_idx);
        subj_data = master_tbl_n(master_tbl_n.CIMBI_ID == subj, :);
        baseline_idx = subj_data.time_since_admin < 3;
        baseline_value = subj_data.(measure)(baseline_idx);
        
        deltas = subj_data.(measure) - baseline_value;
        times = subj_data.time_since_admin;
        
        % Sort by time to ensure line connects points in temporal order
        [sorted_times, sort_idx] = sort(times);
        sorted_deltas = deltas(sort_idx);
        
        % Plot line and points
        plot(sorted_times, sorted_deltas, '-', 'Color', [colors_p(subj_idx,:), 0.3])
        scatter(sorted_times, sorted_deltas, 50, colors_p(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
    end
    
    xlabel('Time since administration (min)')
    ylabel('Δ from baseline')
    title(sprintf('Scanner N - Tau %.1fs', i*1.6))
    ylim([y_min y_max])
    grid on
    hold off
end

subplot(2,5,5)
axis off
subplot(2,5,10)
axis off

savefig('delta_from_baseline_plots.fig')
print('delta_from_baseline_plots','-dpng','-r300')