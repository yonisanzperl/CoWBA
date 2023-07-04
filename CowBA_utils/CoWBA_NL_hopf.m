function [model_output]=CoWBA_NL_hopf(Cfg,TS,SC,f_diff,RSN);


xs=TS;

NPARCELLS=Cfg.nNodes;

NSUBSIM=Cfg.NSIM;

Tmax=Cfg.TmaxSim;

if NSUBSIM==0
    NSUBSIM=NSUB;
end


G_range=Cfg.Glower:Cfg.Gstep:Cfg.Gupper;




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

%%






C=SC/max(max(SC))*0.2;

for nG=1:size(G_range,2)
    G = G_range(nG);
    fprintf('\n \n Coupling Strength G: %f \n',G)    ;
    for sub=1:NSUBSIM

        wC = G*C;
        sumC = repmat(sum(wC,2),1,2);
        
        %% Hopf Simulation
        a=-0.02*ones(NPARCELLS,2);
        xs=zeros(Tmax,NPARCELLS);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
        z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 2000 time steps
        for t=0:dt:2000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts_out{sub}=xs';
        
    end

    [model_output(nG)]= CoWBA_empirical_analysis(Cfg,ts_out);
    
end
    
