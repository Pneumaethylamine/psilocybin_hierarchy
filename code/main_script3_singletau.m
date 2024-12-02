%% Setup and data loading
cd('/indirect/staff/drummondmcculloch/UPF_project');
addpath(genpath('/indirect/staff/drummondmcculloch/UPF_project'));
addpath("/mrdata/proj-doccls/code")

load('/indirect/staff/drummondmcculloch/UPF_project/data/data_09-Mar-2023_SchaeferTian232_ts.mat');
age_sex = readtable('/mrdata/np2/p3/entropy/critical_files/CIMBI_age_sex_clean.xlsx');

% Filter data
tbl = SDI_PPL_rsfmri_tbl(strcmpi(SDI_PPL_rsfmri_tbl.include, 'yes'), :);
tbl = tbl(~contains(lower(tbl.rsfmri_location), 'pa'), :);

%% Split data by scanner
scanner_p_idx = startsWith(tbl.rsfmri_raw_ID, 'p');  % 2.0s TR scans
scanner_n_idx = startsWith(tbl.rsfmri_raw_ID, 'n');  % 0.8s TR scans
tbl_scanner_p = tbl(scanner_p_idx, :);
tbl_scanner_n = tbl(scanner_n_idx, :);

%% Create master tables
master_tbl_p = table(tbl_scanner_p.CIMBI_ID, tbl_scanner_p.psi_conc, ...
    tbl_scanner_p.time_since_admin, tbl_scanner_p.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

master_tbl_n = table(tbl_scanner_n.CIMBI_ID, tbl_scanner_n.psi_conc, ...
    tbl_scanner_n.time_since_admin, tbl_scanner_n.SDI_score, ...
    'VariableNames', {'CIMBI_ID', 'psi_conc', 'time_since_admin', 'SDI_score'});

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

%% Calculate INSIDEOUT for specific taus
% Scanner P - tau 3 (6 seconds)
global_scores = zeros(height(tbl_scanner_p), 1);
for i = 1:height(tbl_scanner_p)
    [~, global_score, ~, ~, asymmetry] = computeDirectionality(tbl_scanner_p.ts_SchaeferTian232{i}, 3, 100);
    global_scores(i) = global_score;
    asymmetries(i) = asymmetry;
end
master_tbl_p.insideout = global_scores;
master_tbl_p.asymmetry = asymmetries';

clear global_scores asymmetries

% Scanner N - tau 8 (6.4 seconds)
global_scores = zeros(height(tbl_scanner_n), 1);
for i = 1:height(tbl_scanner_n)
    [~, global_score, ~, ~, asymmetry] = computeDirectionality(tbl_scanner_n.ts_SchaeferTian232{i}, 8, 100);
    global_scores(i) = global_score;
    asymmetries(i) = asymmetry;
end
master_tbl_n.insideout = global_scores;
master_tbl_n.asymmetry = asymmetries';

%% Calculate modularity for different resolutions
% For Scanner P
fprintf('\nProcessing Scanner P data:\n');
res = 1.0;
mod_scores = zeros(height(tbl_scanner_p), 1);
for i = 1:height(tbl_scanner_p)
    if mod(i,2)==0
        fprintf('Processing subject %d of %d\n', i, height(tbl_scanner_p));
    end
    conn_matrix = corr(tbl_scanner_p.ts_SchaeferTian232{i});
    [norm_mod, ~] = calculate_modularity_metrics2(conn_matrix, 2, res, 16);
    mod_scores(i) = norm_mod;
end
master_tbl_p.modularity_res_1_0 = mod_scores;

% For Scanner N
fprintf('\nProcessing Scanner N data:\n');
mod_scores = zeros(height(tbl_scanner_n), 1);
for i = 1:height(tbl_scanner_n)
    if mod(i,2)==0 
        fprintf('Processing subject %d of %d\n', i, height(tbl_scanner_n));
    end
    conn_matrix = corr(tbl_scanner_n.ts_SchaeferTian232{i});
    [norm_mod, ~] = calculate_modularity_metrics2(conn_matrix, 2, res, 16);
    mod_scores(i) = norm_mod;
end
master_tbl_n.modularity_res_1_0 = mod_scores;

%% Add age and sex
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

%% Stats
% Create results tables
results_p = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});
results_n = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});

% Scanner P models
fprintf('\nProcessing Scanner P models (tau = 3):\n');

% INSIDEOUT
fprintf('Running model for insideout\n');
lme_model = fitlme(master_tbl_p, 'insideout ~ psi_conc + Age + Sex + (1|CIMBI_ID)', 'FitMethod', 'ML');
coef_table = lme_model.Coefficients;
psi_idx = strcmp(coef_table.Name, 'psi_conc');
new_row = {'insideout', coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
    coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
results_p = [results_p; new_row];

% Asymmetry
fprintf('Running model for asymmetry\n');
lme_model = fitlme(master_tbl_p, 'asymmetry ~ psi_conc + Age + Sex + (1|CIMBI_ID)', 'FitMethod', 'ML');
coef_table = lme_model.Coefficients;
psi_idx = strcmp(coef_table.Name, 'psi_conc');
new_row = {'asymmetry', coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
    coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
results_p = [results_p; new_row];

% Scanner N models
fprintf('\nProcessing Scanner N models (tau = 8):\n');

% INSIDEOUT
fprintf('Running model for insideout\n');
lme_model = fitlme(master_tbl_n, 'insideout ~ psi_conc + Age + Sex + (1|CIMBI_ID)', 'FitMethod', 'ML');
coef_table = lme_model.Coefficients;
psi_idx = strcmp(coef_table.Name, 'psi_conc');
new_row = {'insideout', coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
    coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
results_n = [results_n; new_row];

% Asymmetry
fprintf('Running model for asymmetry\n');
lme_model = fitlme(master_tbl_n, 'asymmetry ~ psi_conc + Age + Sex + (1|CIMBI_ID)', 'FitMethod', 'ML');
coef_table = lme_model.Coefficients;
psi_idx = strcmp(coef_table.Name, 'psi_conc');
new_row = {'asymmetry', coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
    coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
results_n = [results_n; new_row];

% Save results
save('lme_results.mat', 'results_p', 'results_n');

% Display results
disp('Scanner P Results (tau = 3):');
disp(results_p);
disp('Scanner N Results (tau = 8):');
disp(results_n);

%% Partial residual plots
figure('Position', [100 100 1600 400]);

% Scanner P plots
% INSIDEOUT
subplot(1,4,1)
mdl_p = fitlme(master_tbl_p, 'insideout ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
resid_p = residuals(mdl_p, 'Conditional', true);
partial_resid_p = resid_p + mdl_p.Coefficients.Estimate(2)*master_tbl_p.psi_conc;

scatter(master_tbl_p.psi_conc, partial_resid_p, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
valid_idx = ~isnan(partial_resid_p) & ~isnan(master_tbl_p.psi_conc);
coef = polyfit(master_tbl_p.psi_conc(valid_idx), partial_resid_p(valid_idx), 1);
x = linspace(min(master_tbl_p.psi_conc), max(master_tbl_p.psi_conc), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

xlabel('PSI Concentration')
ylabel('Partial Residuals')
title('Scanner P - Inside-out (Tau 6s)')
grid on

% ASYMMETRY
subplot(1,4,2)
mdl_p = fitlme(master_tbl_p, 'asymmetry ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
resid_p = residuals(mdl_p, 'Conditional', true);
partial_resid_p = resid_p + mdl_p.Coefficients.Estimate(2)*master_tbl_p.psi_conc;

scatter(master_tbl_p.psi_conc, partial_resid_p, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
valid_idx = ~isnan(partial_resid_p) & ~isnan(master_tbl_p.psi_conc);
coef = polyfit(master_tbl_p.psi_conc(valid_idx), partial_resid_p(valid_idx), 1);
x = linspace(min(master_tbl_p.psi_conc), max(master_tbl_p.psi_conc), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

xlabel('PSI Concentration')
ylabel('Partial Residuals')
title('Scanner P - Asymmetry (Tau 6s)')
grid on

% Scanner N plots
% INSIDEOUT
subplot(1,4,3)
mdl_n = fitlme(master_tbl_n, 'insideout ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
resid_n = residuals(mdl_n, 'Conditional', true);
partial_resid_n = resid_n + mdl_n.Coefficients.Estimate(2)*master_tbl_n.psi_conc;

scatter(master_tbl_n.psi_conc, partial_resid_n, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
valid_idx = ~isnan(partial_resid_n) & ~isnan(master_tbl_n.psi_conc);
coef = polyfit(master_tbl_n.psi_conc(valid_idx), partial_resid_n(valid_idx), 1);
x = linspace(min(master_tbl_n.psi_conc), max(master_tbl_n.psi_conc), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

xlabel('PSI Concentration')
ylabel('Partial Residuals')
title('Scanner N - Inside-out (Tau 6.4s)')
grid on

% ASYMMETRY
subplot(1,4,4)
mdl_n = fitlme(master_tbl_n, 'asymmetry ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
resid_n = residuals(mdl_n, 'Conditional', true);
partial_resid_n = resid_n + mdl_n.Coefficients.Estimate(2)*master_tbl_n.psi_conc;

scatter(master_tbl_n.psi_conc, partial_resid_n, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
valid_idx = ~isnan(partial_resid_n) & ~isnan(master_tbl_n.psi_conc);
coef = polyfit(master_tbl_n.psi_conc(valid_idx), partial_resid_n(valid_idx), 1);
x = linspace(min(master_tbl_n.psi_conc), max(master_tbl_n.psi_conc), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

xlabel('PSI Concentration')
ylabel('Partial Residuals')
title('Scanner N - Asymmetry (Tau 6.4s)')
grid on

savefig('partial_residual_plots.fig')
print('partial_residual_plots','-dpng','-r300')

%% Delta from baseline plots
figure('Position', [100 100 1600 400]);

% Create colormap
colors_p = hsv(length(unique(master_tbl_p.CIMBI_ID)));
colors_n = hsv(length(unique(master_tbl_n.CIMBI_ID)));

% Calculate all deltas for y-axis limits
all_deltas_insideout = [];
all_deltas_asymmetry = [];

% Scanner P
subjects = unique(master_tbl_p.CIMBI_ID);
for subj = subjects'
    subj_data = master_tbl_p(master_tbl_p.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    
    baseline_insideout = subj_data.insideout(baseline_idx);
    baseline_asymmetry = subj_data.asymmetry(baseline_idx);
    
    all_deltas_insideout = [all_deltas_insideout; subj_data.insideout - baseline_insideout];
    all_deltas_asymmetry = [all_deltas_asymmetry; subj_data.asymmetry - baseline_asymmetry];
end

% Scanner N
subjects = unique(master_tbl_n.CIMBI_ID);
for subj = subjects'
    subj_data = master_tbl_n(master_tbl_n.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    
    baseline_insideout = subj_data.insideout(baseline_idx);
    baseline_asymmetry = subj_data.asymmetry(baseline_idx);
    
    all_deltas_insideout = [all_deltas_insideout; subj_data.insideout - baseline_insideout];
    all_deltas_asymmetry = [all_deltas_asymmetry; subj_data.asymmetry - baseline_asymmetry];
end

% Scanner P plots
% INSIDEOUT
subplot(1,4,1)
hold on
subjects = unique(master_tbl_p.CIMBI_ID);
for subj_idx = 1:length(subjects)
    subj = subjects(subj_idx);
    subj_data = master_tbl_p(master_tbl_p.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    baseline_value = subj_data.insideout(baseline_idx);
    
    deltas = subj_data.insideout - baseline_value;
    times = subj_data.time_since_admin;
    
    [sorted_times, sort_idx] = sort(times);
    sorted_deltas = deltas(sort_idx);
    
    plot(sorted_times, sorted_deltas, '-', 'Color', [colors_p(subj_idx,:), 0.3])
    scatter(sorted_times, sorted_deltas, 50, colors_p(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
end

xlabel('Time since administration (min)')
ylabel('Δ from baseline')
title('Scanner P - Inside-out (Tau 6s)')
ylim([min(all_deltas_insideout) max(all_deltas_insideout)])
grid on
hold off

% ASYMMETRY
subplot(1,4,2)
hold on
for subj_idx = 1:length(subjects)
    subj = subjects(subj_idx);
    subj_data = master_tbl_p(master_tbl_p.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    baseline_value = subj_data.asymmetry(baseline_idx);
    
    deltas = subj_data.asymmetry - baseline_value;
    times = subj_data.time_since_admin;
    
    [sorted_times, sort_idx] = sort(times);
    sorted_deltas = deltas(sort_idx);
    
    plot(sorted_times, sorted_deltas, '-', 'Color', [colors_p(subj_idx,:), 0.3])
    scatter(sorted_times, sorted_deltas, 50, colors_p(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
end

xlabel('Time since administration (min)')
ylabel('Δ from baseline')
title('Scanner P - Asymmetry (Tau 6s)')
ylim([min(all_deltas_asymmetry) max(all_deltas_asymmetry)])
grid on
hold off

% Scanner N plots
% INSIDEOUT
subplot(1,4,3)
hold on
subjects = unique(master_tbl_n.CIMBI_ID);
for subj_idx = 1:length(subjects)
    subj = subjects(subj_idx);
    subj_data = master_tbl_n(master_tbl_n.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    baseline_value = subj_data.insideout(baseline_idx);
    
    deltas = subj_data.insideout - baseline_value;
    times = subj_data.time_since_admin;
    
    [sorted_times, sort_idx] = sort(times);
    sorted_deltas = deltas(sort_idx);
    
    plot(sorted_times, sorted_deltas, '-', 'Color', [colors_n(subj_idx,:), 0.3])
    scatter(sorted_times, sorted_deltas, 50, colors_n(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
end

xlabel('Time since administration (min)')
ylabel('Δ from baseline')
title('Scanner N - Inside-out (Tau 6.4s)')
ylim([min(all_deltas_insideout) max(all_deltas_insideout)])
grid on
hold off

% ASYMMETRY
subplot(1,4,4)
hold on
for subj_idx = 1:length(subjects)
    subj = subjects(subj_idx);
    subj_data = master_tbl_n(master_tbl_n.CIMBI_ID == subj, :);
    baseline_idx = subj_data.time_since_admin < 3;
    baseline_value = subj_data.asymmetry(baseline_idx);
    
    deltas = subj_data.asymmetry - baseline_value;
    times = subj_data.time_since_admin;
    
    [sorted_times, sort_idx] = sort(times);
    sorted_deltas = deltas(sort_idx);
    
    plot(sorted_times, sorted_deltas, '-', 'Color', [colors_n(subj_idx,:), 0.3])
    scatter(sorted_times, sorted_deltas, 50, colors_n(subj_idx,:), 'filled', 'MarkerFaceAlpha', 0.8)
end

xlabel('Time since administration (min)')
ylabel('Δ from baseline')
title('Scanner N - Asymmetry (Tau 6.4s)')
ylim([min(all_deltas_asymmetry) max(all_deltas_asymmetry)])
grid on
hold off

savefig('delta_from_baseline_plots.fig')
print('delta_from_baseline_plots','-dpng','-r300')

%% Plot correlations between inside-out and asymmetry
figure('Position', [100 100 800 400]);

% Scanner P correlation
subplot(1,2,1)
scatter(master_tbl_p.insideout, master_tbl_p.asymmetry, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on

% Add regression line
valid_idx = ~isnan(master_tbl_p.insideout) & ~isnan(master_tbl_p.asymmetry);
coef = polyfit(master_tbl_p.insideout(valid_idx), master_tbl_p.asymmetry(valid_idx), 1);
x = linspace(min(master_tbl_p.insideout), max(master_tbl_p.insideout), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

% Calculate correlation and p-value
[r_p, p_p] = corr(master_tbl_p.insideout, master_tbl_p.asymmetry, 'rows', 'complete');

% Add correlation info to plot
text(0.05, 0.95, sprintf('r = %.3f\np = %.3e', r_p, p_p), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

xlabel('Inside-out')
ylabel('Asymmetry')
title('Scanner P (Tau 6s)')
grid on

% Scanner N correlation
subplot(1,2,2)
scatter(master_tbl_n.insideout, master_tbl_n.asymmetry, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.5)
hold on

% Add regression line
valid_idx = ~isnan(master_tbl_n.insideout) & ~isnan(master_tbl_n.asymmetry);
coef = polyfit(master_tbl_n.insideout(valid_idx), master_tbl_n.asymmetry(valid_idx), 1);
x = linspace(min(master_tbl_n.insideout), max(master_tbl_n.insideout), 100);
y = polyval(coef, x);
plot(x, y, 'r-', 'LineWidth', 2)

% Calculate correlation and p-value
[r_n, p_n] = corr(master_tbl_n.insideout, master_tbl_n.asymmetry, 'rows', 'complete');

% Add correlation info to plot
text(0.05, 0.95, sprintf('r = %.3f\np = %.3e', r_n, p_n), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

xlabel('Inside-out')
ylabel('Asymmetry')
title('Scanner N (Tau 6.4s)')
grid on

savefig('metric_correlations.fig')
print('metric_correlations','-dpng','-r300')

% Display correlations
fprintf('\nCorrelations between Inside-out and Asymmetry:\n');
fprintf('Scanner P: r = %.3f, p = %.3e\n', r_p, p_p);
fprintf('Scanner N: r = %.3f, p = %.3e\n', r_n, p_n);