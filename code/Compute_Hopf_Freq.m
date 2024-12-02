clear all;

NSUB=15;  % change depending on your dataset
NPARCELLS=90; 
Tmax=198;

%%%
% Parameters of the data
TR=2;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

%%%%
fce1=zeros(NSUB,NPARCELLS,NPARCELLS);

%PowSpect=zeros(600,NPARCELLS,12);
load('DataSleepW_N3.mat') % load you timeseries file
NSUB=15;
for sub=1:NSUB
    sub
   
    ts=TS_W{1,sub}(:,1:Tmax); % change with your data the name of this variable

    %[Ns, Tmax]=size(ts);
    TT=Tmax;
    tss=zeros(NPARCELLS,Tmax);
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    
    for seed=1:NPARCELLS
        x(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        tss(seed,:)=filtfilt(bfilt,afilt,x(seed,:));
        pw = abs(fft(tss(seed,:)));
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    fce1(sub,:,:)=corrcoef(tss','rows','pairwise');
end
fce=squeeze(mean(fce1));

Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NPARCELLS
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

save hopf_freq.mat f_diff fce;


