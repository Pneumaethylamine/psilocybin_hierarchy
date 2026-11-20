% ---- paths: set these environment variables (see README) ----
project_dir = getenv('UPF_PROJECT_DIR');   % project root (data/, outputs/, ROIs/)
np2_dir     = getenv('NP2_DATA_DIR');      % NeuroPharm2 derived data (restricted access)
bct_dir     = getenv('BCT_DIR');           % Brain Connectivity Toolbox
roi_dir     = getenv('ROI_DIR');           % Schaefer atlas files

function [directionality_scores, FCtf, MI_matrix, FC_matrix, global_score, forward_MI, reverse_MI, insideout_regional_dirscore] = compute_insideout(ts, Tau, edge_percent)
    % Compute directionality scores from time series data
    %
    % Inputs:
    %   ts: Time series data (time x regions double)
    %   Tau: Time delay
    %   edge_percent: Percentage of edges to keep (default=100)
    %
    % Outputs:
    %   MI_matrix: Mutual information matrix
    %   FC_matrix: Functional connectivity (correlation) matrix
    %   global_score: Global directionality score
    %   forward_MI: Forward mutual information matrix
    %   reverse_MI: Reverse mutual information matrix
    %   insideout_regional_dirscore: Regional directionality scores
    
    if nargin < 3
        edge_percent = 100;
    end

    % Get dimensions
    Tm = size(ts, 1);
    num_regions = size(ts, 2);

    FC_matrix = corr(ts);
    MI_matrix = -0.5*log(1-FC_matrix.*FC_matrix);
    MI_matrix(logical(eye(size(MI_matrix)))) = 0;
    
    FCtf = corr(ts(1:Tm-Tau,:), ts(1+Tau:Tm,:));
    FCtr = FCtf';
    forward_MI = -0.5*log(1-FCtf.*FCtf);
    reverse_MI = -0.5*log(1-FCtr.*FCtr);
 
    
    % Compute directionality scores
    directionality_scores = (forward_MI - reverse_MI).^2;
    
    % Compute global score
    all_scores = directionality_scores(:);
    if edge_percent < 100
        threshold = prctile(all_scores, 100 - edge_percent);
        global_score = mean(all_scores(all_scores >= threshold));
    else
        global_score = mean(all_scores);
    end
    
    % Compute regional directionality scores
    insideout_regional_dirscore = zeros(num_regions, 1);
    for i = 1:num_regions
        insideout_regional_dirscore(i) = mean(directionality_scores(i,:));
    end
    
end