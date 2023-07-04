function [Ceffgroup CeffPB]=CoWBA_eff_conn_linearmodel(Cfg,TSk,SC,f_diff,RSN);
%%  Read the empirical data



xs=TSk;

NPARCELLS=Cfg.nNodes;
N=NPARCELLS;

NSUB=size(find(~cellfun(@isempty,xs)),1);



Isubdiag = find(tril(ones(NPARCELLS),-1));

C=zeros(NPARCELLS,NPARCELLS);


%%%
% Parameters of the data
TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                    % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

% Parameters HOPF
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;



Tau=1;
sigma=0.01;

epsFC=0.0002;
epsFCtau=0.00004;
maxC=0.1;




for nsub=1:NSUB;
    ts=xs{nsub};
    Tmax=size(ts,2);
    tsdata=ts(:,1:Tmax);
    FCdata(nsub,:,:)=corrcoef(squeeze(tsdata)');
end


C = SC;  %% anatomie..DTI tractography
C = C/max(max(C))*maxC;

% global e subject-level

for nsub=1:NSUB
    ts1=xs{nsub};
    Tmax=size(ts1,2);
    ts=ts1(:,1:Tmax);
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(:,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FC(nsub,:,:)=FCemp; 
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:NPARCELLS
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtau(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FC));
COVtauemp=squeeze(mean(COVtau));
Cnew=C;
olderror=100000;
for iter=1:5000
   % iter
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma);
    TR = 2.2;
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:NPARCELLS
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

    for i=1:NPARCELLS  %% learning
        for j=1:NPARCELLS
            if (C(i,j)>0 )
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
    Cnew_iter(iter,:,:)=Cnew;
end

% [a b]=min(errorFC);
% Ceff=squeeze(Cnew_iter(b,:,:));
% clear Cnew_iter errorFC

Ceffgroup=Cnew;

%% Individual 
for nsub=1:NSUB
    fprintf('\n \n GEC subject: %f \n',nsub)    ; 
    ts1=xs{nsub};
    Tmax=size(ts1,2);
    ts=ts1(:,1:Tmax);
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(:,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FC(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=Ceffgroup;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:NPARCELLS
            for j=1:NPARCELLS
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

        for i=1:NPARCELLS  %% learning
            for j=1:NPARCELLS
                if (C(i,j)>0)
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
    CeffPB(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff,sigma);
    fittFC_PB(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_PB(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end


