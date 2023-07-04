function [f_diff,fce]=CoWBA_Empirical_freq(Cfg,SC,TSk);


NPARCELLS=Cfg.nNodes;


TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                    % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter



xs=TSk;
NSUB=size(find(~cellfun(@isempty,xs)),1);

for sub=1:NSUB
    
    ts=xs{sub,1};
    clear signal_filt 
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        
        if sum(isnan(ts(seed,:)))<1
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        else
            ts(seed,:)=ts(seed,:);
            signal_filt(seed,:) = ts(seed,:);
        end
    end
    
    fce=corrcoef(signal_filt');
end

ts=xs{sub,1};
[Ns, Tmax]=size(ts);
TT=Tmax;
Ts = TT*TR;
freq = (0:TT/2-1)/Ts;
nfreqs=length(freq);

for seed=1:NPARCELLS
    
    if sum(isnan(ts(seed,:)))<1
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        tss(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    else
        tss(seed,:)=ts(seed,:);
    end
    
    
    pw = abs(fft(tss(seed,:)));
    PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
end


Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NPARCELLS
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));



end