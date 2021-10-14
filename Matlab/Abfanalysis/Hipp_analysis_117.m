% patient 117 hippocampus analysis 
p1   = ('//Users/elinemertens/Data/ephys/Hippocampus/H17.29.117.21/resonance_all');
fn1  = '2017_03_30_0179.abf';
fp1  = fullfile(p1,fn1);
ss  = load('/Users/elinemertens/Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf1 = Abffile(fp1,ss);

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

%% For djai different channels 
signal = abf1.channels.analogins(1, 1).sweeps.avtrace
plot(signal)
signal = signal.low_pass_filter(1000);
data = abfload_pro(fp1);
stim = squeeze(data(:,2,1));
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;

%% thijs
signal = abf1.channels(1, 2).analogins.sweeps.avtrace
plot(signal)
signal = signal.low_pass_filter(1000);
data = abfload_pro(fp1);
stim = squeeze(data(:,4,1));
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;


%%
 Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;

    plot(signal)
    grid off
    
                 ylabel('Voltage (mV)')
              xlabel('Time (sec)')
    
    
    %%
    signal_fft=fft(signal);
    stim_fft=fft(stim);

    %plot

%     %raw
    plot(f, abs((signal_fft)./abs(stim_fft)));
    xlim([1.2 25])
    hold on
    
    
        %smoothed
    smoothed=smooth(abs((signal_fft)./abs(stim_fft)),23,  'lowess');
    plot(f, smoothed);
    xlim([1.2 25])
    
    %% plot only smoothed 
    
    smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([1.2 25])
    grid off   
            
 ylabel('MOhm')
xlabel('Frequency (Hz)')
set(gca, 'TickDir', 'out')
    % Set the remaining axes properties

   
   %% resonance freq
   
   [max_imp, loc] = max(smoothed(f>1.2 & f<25));
    f2 = f(f>1.2 & f<25);
    res_freq = f2(loc);
    
   %% 3db cutoff
   normalized = smoothed(f>0.5 & f<25)./max_imp;
   % when does normalized frequency response go below square root of 0.5?
   % This is the definition of 3 dB cutoff!
   cutoff_3db=f2(find(normalized>=sqrt(0.5), 1, 'last'));
   figure
   plot(f2, normalized)
   hold on
   line(xlim,[sqrt(0.5), sqrt(0.5)],'linewidth',2,'lineStyle',':','color','k')
   %line([cutoff_3db, cutoff_3db],ylim,'linewidth',2,'lineStyle',':','color','k')
 ylim([0.1 1.2])
 set(gca, 'TickDir', 'out')

