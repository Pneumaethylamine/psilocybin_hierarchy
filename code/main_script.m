cd('/indirect/staff/drummondmcculloch/UPF_project');
addpath(genpath('/indirect/staff/drummondmcculloch/UPF_project'));
addpath("/mrdata/proj-doccls/code")

load('/indirect/staff/drummondmcculloch/UPF_project/data/data_09-Mar-2023_SchaeferTian232_ts.mat');

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

%% Calculate modularity for different resolutions
% resolutions = [0.5, 1.0, 1.5, 2];
% 
% % For Scanner P
% fprintf('\nProcessing Scanner P data:\n');
% for res = resolutions
%     fprintf('\nStarting resolution %.1f:\n', res);
%     mod_scores = zeros(height(tbl_scanner_p), 1);
%     for i = 1:height(tbl_scanner_p)
%         if mod(i,2)==0
%             fprintf('Processing subject %d of %d\n', i, height(tbl_scanner_p));
%         end
%         conn_matrix = corr(tbl_scanner_p.ts_SchaeferTian232{i});
%         [norm_mod, ~] = calculate_modularity_metrics2(conn_matrix, 2, res, 16);
%         mod_scores(i) = norm_mod;
%     end
%     master_tbl_p.(['modularity_res_' num2str(res,'%.1f')]) = mod_scores;
%     fprintf('Completed resolution %.1f\n', res);
% end
% 
% % For Scanner N
% fprintf('\nProcessing Scanner N data:\n');
% for res = resolutions
%     fprintf('\nStarting resolution %.1f:\n', res);
%     mod_scores = zeros(height(tbl_scanner_n), 1);
%     for i = 1:height(tbl_scanner_n)
%         if mod(i,2)==0 
%             fprintf('Processing subject %d of %d\n', i, height(tbl_scanner_n));
%         end
%         conn_matrix = corr(tbl_scanner_n.ts_SchaeferTian232{i});
%         [norm_mod, ~] = calculate_modularity_metrics2(conn_matrix, 2, res, 16);
%         mod_scores(i) = norm_mod;
%     end
%     master_tbl_n.(['modularity_res_' num2str(res,'%.1f')]) = mod_scores;
%     fprintf('Completed resolution %.1f\n', res);
% end

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

%% Add age and sex 
age_sex = readtable('/mrdata/np2/p3/entropy/critical_files/CIMBI_age_sex_clean.xlsx');
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
% modularity_measures = {'modularity_res_0.5', 'modularity_res_1.0'};
insideout_measures_p = {'insideout_tau_1', 'insideout_tau_2', 'insideout_tau_3', 'insideout_tau_4'};
insideout_measures_n = {'insideout_tau_2', 'insideout_tau_4', 'insideout_tau_6', 'insideout_tau_8', 'insideout_tau_10'};

% Create results tables structure
results_p = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});
results_n = table('Size', [0, 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Measure', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'});

%% Run models for Scanner P
fprintf('\nProcessing Scanner P models:\n');

% % Modularity measures
% for i = 1:length(modularity_measures)
%     measure = modularity_measures{i};
%     fprintf('Running model for %s\n', measure);
% 
%     % Fit model
%     lme_model = fitlme(master_tbl_p, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
% 
%     % Extract results for psi_conc effect
%     coef_table = lme_model.Coefficients;
%     psi_idx = strcmp(coef_table.Name, 'psi_conc');
% 
%     % Store results
%     new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
%         coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
%     results_p = [results_p; new_row];
% end

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

%% Run models for Scanner N
fprintf('\nProcessing Scanner N models:\n');

% % Modularity measures
% for i = 1:length(modularity_measures)
%     measure = modularity_measures{i};
%     fprintf('Running model for %s\n', measure);
% 
%     % Fit model
%     lme_model = fitlme(master_tbl_n, sprintf('%s ~ psi_conc + Age + Sex + (1|CIMBI_ID)', measure), 'FitMethod', 'ML');
% 
%     % Extract results for psi_conc effect
%     coef_table = lme_model.Coefficients;
%     psi_idx = strcmp(coef_table.Name, 'psi_conc');
% 
%     % Store results
%     new_row = {measure, coef_table.Estimate(psi_idx), coef_table.SE(psi_idx), ...
%         coef_table.tStat(psi_idx), coef_table.DF(psi_idx), coef_table.pValue(psi_idx)};
%     results_n = [results_n; new_row];
% end

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

% Save results
save('lme_results.mat', 'results_p', 'results_n');

% Display results
disp('Scanner P Results:');
disp(results_p);
disp('Scanner N Results:');
disp(results_n);

%% Visualization

% Create figure for insideout results comparison
figure('Position', [100 100 1200 600]);

% Scanner P subplot
subplot(1,2,1)
errorbar(1:4, results_p.Estimate, results_p.SE, 'ko-', 'LineWidth', 1.5)
title('Scanner P (TR=2.0s)')
xticks(1:4)
xticklabels({'2s', '4s', '6s', '8s'})
xlabel('Tau (seconds)')
ylabel('Estimate (psi_conc effect)')
grid on

% Scanner N subplot
subplot(1,2,2)
errorbar(1:5, results_n.Estimate, results_n.SE, 'ko-', 'LineWidth', 1.5)
title('Scanner N (TR=0.8s)')
xticks(1:5)
xticklabels({'1.6s', '3.2s', '4.8s', '6.4s', '8.0s'})
xlabel('Tau (seconds)')
ylabel('Estimate (psi_conc effect)')
grid on

% Add significance markers
for i = 1:4
   subplot(1,2,1)
   if results_p.pValue(i) < 0.05
       text(i, results_p.Estimate(i) + 2*results_p.SE(i), '*', 'FontSize', 14, 'HorizontalAlignment', 'center')
   end
end

for i = 1:5
   subplot(1,2,2) 
   if results_n.pValue(i) < 0.05
       text(i, results_n.Estimate(i) + 2*results_n.SE(i), '*', 'FontSize', 14, 'HorizontalAlignment', 'center')
   end
end

savefig('insideout_results.fig')
print('insideout_results','-dpng','-r300')

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

savefig('partial_residual_plots.fig')
print('partial_residual_plots','-dpng','-r300')

%% Regional effects

% For scanner P (tau=2) and scanner N (tau=6)
regional_scores_p = zeros(232, 232, height(tbl_scanner_p));
regional_scores_n = zeros(232, 232, height(tbl_scanner_n));

% Compute regional scores
for i = 1:height(tbl_scanner_p)
   [regional_scores_p(:,:,i), ~, ~, ~] = computeDirectionality(tbl_scanner_p.ts_SchaeferTian232{i}, 2);
end

for i = 1:height(tbl_scanner_n)
   [regional_scores_n(:,:,i), ~, ~, ~] = computeDirectionality(tbl_scanner_n.ts_SchaeferTian232{i}, 6);
end

% Mass univariate analysis for each scanner
pvals_p = zeros(232, 232);
estimates_p = zeros(232, 232);
pvals_n = zeros(232, 232); 
estimates_n = zeros(232, 232);

% Scanner P analysis
for i = 1:232
   for j = 1:232
       if i ~= j
           conn_scores = squeeze(regional_scores_p(i,j,:));
           tmp_tbl = master_tbl_p;
           tmp_tbl.conn_score = conn_scores;
           mdl = fitlme(tmp_tbl, 'conn_score ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
           pvals_p(i,j) = mdl.Coefficients.pValue(2);
           estimates_p(i,j) = mdl.Coefficients.Estimate(2);
       end
   end
end

% Scanner N analysis
for i = 1:232
   for j = 1:232
       if i ~= j
           conn_scores = squeeze(regional_scores_n(i,j,:));
           tmp_tbl = master_tbl_n;
           tmp_tbl.conn_score = conn_scores;
           mdl = fitlme(tmp_tbl, 'conn_score ~ psi_conc + Age + Sex + (1|CIMBI_ID)');
           pvals_n(i,j) = mdl.Coefficients.pValue(2);
           estimates_n(i,j) = mdl.Coefficients.Estimate(2);
       end
   end
end

% FDR correction and visualization
[~, ~, adj_p_p] = fdr_bh(pvals_p(:));
[~, ~, adj_p_n] = fdr_bh(pvals_n(:));

sig_conn_p = reshape(adj_p_p < 0.05, 232, 232);
sig_conn_n = reshape(adj_p_n < 0.05, 232, 232);

% Plot results
figure('Position', [100 100 1200 500]);

subplot(1,2,1)
imagesc(estimates_p .* sig_conn_p);
colorbar;
title('Scanner P (tau=2) - Significant PSI Effects');
axis square

subplot(1,2,2)
imagesc(estimates_n .* sig_conn_n);
colorbar;
title('Scanner N (tau=6) - Significant PSI Effects');
axis square

savefig('regional_effects.fig')
print('regional_effects','-dpng','-r300')