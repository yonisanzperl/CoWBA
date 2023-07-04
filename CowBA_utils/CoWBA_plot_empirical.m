function CoWBA_plot_empirical(output,Cfg,metadata)



make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.015], [0.12 0.01], [0.12 0.01]);
if ~make_it_tight,  clear subplot;  end


% Metastability

groups =metadata.group;
stattest = metadata.stattest;

NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;

C = cell(1,Cfg.nBrainStates);
f1=figure('Name','Metastability');

cont=1;
for ii=1:Cfg.nBrainStates
    
    C{1,ii}= output(ii).meta;
end


p = swarm(C, groups, tlt=sprintf('Metastability'), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_turbuStats_lambda %.2f.txt',ii));
%    title(sprintf('lambda %.2f',LAMBDA(i)));
%    movefile(sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)


% syncrho

C = cell(1,Cfg.nBrainStates);
f1=figure('Name','Synchronization');

cont=1;
for ii=1:Cfg.nBrainStates
    
    C{1,ii}= output(ii).sync;
end


p = swarm(C, groups, tlt=sprintf('Synchronization'), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_turbuStats_lambda %.2f.txt',ii));
%    title(sprintf('lambda %.2f',LAMBDA(i)));
%    movefile(sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)




%f1.Position = [100 100 740 600];



% Edge Meta
f1=figure('Name','Edge Meta');

cont=1;
for ii=1:Cfg.nBrainStates
    
    C{1,ii}= output(ii).EdgeMeta;
end


p = swarm(C, groups, tlt=sprintf('Edge Meta'), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_turbuStats_lambda %.2f.txt',ii));
%    title(sprintf('lambda %.2f',LAMBDA(i)));
%    movefile(sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)

% Complexity


f1=figure('Name','LZ complexity');

cont=1;
for ii=1:Cfg.nBrainStates
    
    C{1,ii}= output(ii).C_mean;
end


p = swarm(C, groups, tlt=sprintf('LZ complexity'), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_turbuStats_lambda %.2f.txt',ii));
%    title(sprintf('lambda %.2f',LAMBDA(i)));
%    movefile(sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)





%f1.Position = [100 100 740 600];


end


