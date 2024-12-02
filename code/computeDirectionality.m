function [directionality_scores, global_score, forward_MI, reverse_MI, assymmetry] = computeDirectionality(ts, Tau, edge_percent)
    % Compute directionality scores from time series data
    %
    % Inputs:
    %   ts: Time series data (time x regions double)
    %   Tau: Time delay (default=2)
    %   edge_percent: Percentage of edges to keep (default=50)
    %
    % Outputs:
    %   directionality_scores: Matrix of directional differences (regions x regions)
    %   global_score: Single number summarizing overall directionality
    %   forward_MI: Forward mutual information matrix
    %   reverse_MI: Reverse mutual information matrix
    %   assymmetry: Mean absolute asymmetry of temporal correlations
    
    if nargin < 2
        Tau = 2;
    end
    if nargin < 3
        edge_percent = 50;
    end
    
    % Get dimensions
    Tm = size(ts, 1);
    
    % Compute forward and reverse correlations
    FCtf = corr(ts(1:Tm-Tau,:), ts(1+Tau:Tm,:));
    FCtr = corr(ts(Tm:-1:Tau+1,:), ts(Tm-Tau:-1:1,:));
    
    % Convert to mutual information
    forward_MI = -0.5*log(1-FCtf.*FCtf);
    reverse_MI = -0.5*log(1-FCtr.*FCtr);
    diff_mat = FCtf - FCtf';
    assymmetry = mean(mean(abs(diff_mat)));
    
    % Compute directionality scores (keeping matrix structure)
    directionality_scores = (forward_MI - reverse_MI).^2;
    
    % Compute global score using specified percentage of edges
    all_scores = directionality_scores(:);
    threshold = prctile(all_scores, 100 - edge_percent);
    global_score = mean(all_scores(all_scores >= threshold));
end