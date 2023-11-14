
%% Plot all sweeps and analyse 
abf1.channels(1, 2).analogins(1, 1).sweeps.plot
 lot(abf1.channels(1, 2).analogins(1, 1).sweeps(1, 13))


 %% plot only the AP
abf1.channels(1, 2).analogins(1, 1).sweeps(1, 13).aps(1).plot('superimpose','peak')

legend show
figure(1)
title('AP traces Par')
xlabel('time (ms)')
ylabel('mV')
xlim([-2 6])
 
 
%% resonance

abf1.channels(1, 2).analogins.plot