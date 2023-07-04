function [output]=empirical_metaesta(Cfg,SC,ts_out,RSN);

NPARCELLS=Cfg.nNodes;
NSUB=size(find(~cellfun(@isempty,ts_out)),1);


TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                    % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));


for nsub=1:NSUB
    ts=ts_out{nsub};
    if Cfg.Tmax>0
        Tmax=Cfg.Tmax;                                   % Timepoints
    else
        Tmax= size(ts{1},2);
    end
    ts = ts(:,1:Tmax);
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        
        if sum(isnan(ts(seed,:)))<1
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        else
            ts(seed,:)=ts(seed,:);
            signal_filt(seed,:) = ts(seed,:);
        end
        
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        PhasesE(seed,:) = angle(Xanalytic);
    end
    
    
     gKOM=nansum(complex(cos(PhasesE),sin(PhasesE)))/NPARCELLS;
    enstrophy1=abs(gKOM);
   % Meta and Synchro
    Metaemp_sub(nsub)=nanstd(enstrophy1(:)); %metastability
    Syncrhoemp_sub(nsub)=nanmean(enstrophy1(:)); % synchronization
    
    % FC and GBC
    FCemp2(nsub,:,:)=corrcoef(signal_filt(:,:)');
    
end
FCemp = squeeze(mean(FCemp2,1));

output.Metaemp_sub=Metaemp_sub;
output.Syncrhoemp_sub=Syncrhoemp_sub;
output.FCemp=FCemp;


