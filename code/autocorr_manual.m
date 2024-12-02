function r = autocorr_manual(x, maxlag)
    n = length(x);
    x = x - mean(x); % Remove mean
    % Normalize by variance
    x = x / std(x);
    r = zeros(maxlag + 1, 1);
    for k = 0:maxlag
        r(k + 1) = sum(x(1:n-k) .* x(k+1:n)) / n;
    end
end