% patient 117 hippocampus analysis 
p1   = ('/Users/elinemertens/Data/ephys/Hippocampus/2017_hipp_ephys/2017_03_29_IV');
fn1  = '2017_03_29_0076.abf';
fp1  = fullfile(p1,fn1);
ss  = load('/Users/elinemertens/Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf1 = Abffile(fp1,ss);

%% Plot all sweeps and analyse 
abf1.channels(1, 2).analogins(1, 1).sweeps.plot


%% Find AP sweep and plot
 plot(abf1.channels(1, 2).analogins(1, 1).sweeps(1, 13))


 %% plot only the AP
abf1.channels(1, 2).analogins(1, 1).sweeps(1, 14).aps(1).plot('superimpose','peak')

legend show
figure(1)
title('AP traces Par')
xlabel('time (ms)')
ylabel('mV')
xlim([-2 6])
 
 
%% s1 = s1.analysesweep;



%%
