function subjectGEC(s,cnd)

load('XX.mat'); % the data acordingly the subject (s) and condition (brain state)
load('XX.mat'); % region freq from the Compute_Hopf_Freq.m

NSUB=length(BOLD_sbj); % change with the name that you give to the ts 

Tau=1; %CHANGE ACCORDING TO YOUR DATA AND THE SCAN
sigma=0.01;
N=1000; % CHANGE
NPARCELLS=N;
epsFC=0.0002;%0.0002
epsFCtau=0.0005; %0.00004
maxC=0.1;

Isubdiag = find(tril(ones(N),-1));


ts=BOLD_sbj{s,cnd};  % fMRI data organized in row subjet and column the
% condition....change with the name that you give to the ts 


% Parameters of the data
TR=2;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter



load('sc_schaefer_MK.mat') % CHANGE FOR THE PROPER IN SCH200
C = sc_schaefer;  %% anatomie..DTI tractography
C = C/max(max(C))*0.2;


% empirical
clear signal_filt;
for seed=1:N
    ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
end
% FC(0)
ts2=ts(:,:);
Tm=size(ts2,2);
FCemp=corrcoef(ts2');
FC(:,:)=FCemp;
COVemp=cov(ts2');
% COV(tau)
tst=ts2';
for i=1:N
    for j=1:N
        sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
        [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
        indx=find(lags==Tau);
        COVtauemp(i,j)=clag(indx)/size(tst,1);
    end
end
COVtauemp=COVtauemp.*sigratio;


% compute GEC
olderror=100000;
Cnew=C;
for iter=1:5000
	iter
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff(:,cnd),sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));
    
    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end
    
    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || i==j+100 || j==i+100)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
Ceff=Cnew;

[FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff(:,cnd),sigma);
fittFC=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
COVtausim=expm((Tau*TR)*A)*COVsimtotal;
COVtausim=COVtausim(1:N,1:N);
for i=1:N
    for j=1:N
        sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
    end
end
COVtausim=COVtausim.*sigratiosim;
fittCVtau=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));


save(sprintf('LSD_cndCOVtau_%03d_%03d.mat',cnd,s),'Ceff','fittCVtau','fittFC');

