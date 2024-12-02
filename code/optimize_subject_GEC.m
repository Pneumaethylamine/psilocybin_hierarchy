function Ceff = optimize_subject_GEC(subject_ts, f_diff, TR, sc_matrix, Tau, use_time_lag)
    % Function to optimize Global Effective Connectivity (GEC) for a single subject
    % Inputs:
    % subject_ts: single subject timeseries [time x regions]
    % f_diff: region frequencies from hopf frequency calculation
    % TR: repetition time (2.0s for P, 0.8s for N)
    % sc_matrix: structural connectivity matrix
    % scanner_type: 'P' or 'N'
    % use_time_lag: boolean to include time-lagged covariance in optimization

       % Set epsFCtau based on use_time_lag
    if use_time_lag
        epsFCtau = 0.0005;
    else
        epsFCtau = 0;
    end
   
   sigma = 0.01;
   N = 200;  % ROIs
   epsFC = 0.0002;
   maxC = 0.1;
   Isubdiag = find(tril(ones(N),-1));
   
   % Bandpass filter settings
   fnq = 1/(2*TR);                 % Nyquist frequency
   flp = 0.008;                    % lowpass frequency (Hz)
   fhi = 0.08;                     % highpass
   Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
   k = 2;                          % 2nd order butterworth filter
   [bfilt,afilt] = butter(k,Wn);   % construct the filter
   
   % Preprocess timeseries
   subject_ts = subject_ts';
   subject_ts = subject_ts(1:N,:);
   ts2 = zeros(size(subject_ts));
   for seed = 1:N
       x = detrend(subject_ts(seed,:) - mean(subject_ts(seed,:)));
       ts2(seed,:) = filtfilt(bfilt,afilt,x);
   end
   
   % Calculate empirical FC(0)
   FCemp = corrcoef(ts2');
   COVemp = cov(ts2');
   
   % Calculate empirical COV(tau)
   tst = ts2';
   for i = 1:N
       for j = 1:N
           sigratio(i,j) = 1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
           [clag, lags] = xcov(tst(:,i), tst(:,j), Tau);
           indx = find(lags == Tau);
           COVtauemp(i,j) = clag(indx)/size(tst,1);
       end
   end
   COVtauemp = COVtauemp.*sigratio;
   
   % Optimize GEC
   olderror = 100000;
   Cnew = sc_matrix;
   
   for iter = 1:5000
       fprintf('Iteration %d/5000\n', iter);
       
       % Get simulated FC and covariances
       [FCsim, COVsim, COVsimtotal, A] = hopf_int(Cnew, f_diff, sigma);
       COVtausim = expm((Tau*TR)*A)*COVsimtotal;
       COVtausim = COVtausim(1:N,1:N);
       
       % Calculate simulation ratios
       for i = 1:N
           for j = 1:N
               sigratiosim(i,j) = 1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
           end
       end
       COVtausim = COVtausim.*sigratiosim;
       
       % Calculate errors
       errorFC(iter) = mean(mean((FCemp-FCsim).^2));
       errorCOVtau(iter) = mean(mean((COVtauemp-COVtausim).^2));
       
       % Check convergence
       if mod(iter,100) < 0.1
           errornow = mean(mean((FCemp-FCsim).^2)) + mean(mean((COVtauemp-COVtausim).^2));
           if (olderror-errornow)/errornow < 0.001
               break;
           end
           if olderror < errornow
               break;
           end
           olderror = errornow;
       end
       
       % Update connectivity
       for i = 1:N
           for j = 1:N
               if (sc_matrix(i,j)>0 || i==j+100 || j==i+100)
                   Cnew(i,j) = Cnew(i,j) + epsFC*(FCemp(i,j)-FCsim(i,j)) + ...
                       epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                   if Cnew(i,j) < 0
                       Cnew(i,j) = 0;
                   end
               end
           end
       end
       Cnew = Cnew/max(max(Cnew))*maxC;
   end
   
   % Final calculation
   Ceff = Cnew;
   [FCsim, COVsim, COVsimtotal, A] = hopf_int(Ceff, f_diff, sigma);
   
   % Calculate final fits
   fittFC = corr2(FCemp(Isubdiag), FCsim(Isubdiag));
   COVtausim = expm((Tau*TR)*A)*COVsimtotal;
   COVtausim = COVtausim(1:N,1:N);
   for i = 1:N
       for j = 1:N
           sigratiosim(i,j) = 1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
       end
   end
   COVtausim = COVtausim.*sigratiosim;
   fittCVtau = corr2(COVtauemp(Isubdiag), COVtausim(Isubdiag));
   
   fprintf('Final FC fit: %.3f\n', fittFC);
   fprintf('Final COVtau fit: %.3f\n', fittCVtau);
end