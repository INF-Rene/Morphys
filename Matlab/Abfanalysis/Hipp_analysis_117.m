% patient 117 hippocampus analysis  for resonance
p1   = ('/Volumes/WD ELEMENTS/Human ephys/H20.29.182_F/cell3');
fn1  = '2020_09_23_0019.abf';
fp1  = fullfile(p1,fn1);
ss  = load('/Users/elinemertens/Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
zabf1 = Abffile(fp1,ss);



%% For djai different channels 
signal = zabf1.channels.analogins(1, 1).sweeps.avtrace
plot(signal)
signal = signal.low_pass_filter(1000);
data = abfload_pro(fp1);
stim = squeeze(data(:,2,1));
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;

%% thijs
signal = zabf1.channels(1, 2).analogins.sweeps.avtrace
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
              
  %% for tijs you do channels (1,2) and for djai not 
               
time = zabf1.channels(1,2).analogins(1,1).sweeps.Time ;  
 aastart_rmp = mean(signal(time>100 & time<140));

[first_wave, loc] = max(signal(time>140 & time<700));
 adelta_first_wave = (aastart_rmp) - (first_wave) ; 

[max_ampl, loc] = max(signal(time>300 & time<10000));
t2 = time(time>300 & time<10000);
  time_max_amp = t2(loc);
  cdelta_max_wave = (aastart_rmp) - (max_ampl) ;

cend_rmp = mean(signal(time>10000 & time<10200)) ;
    
    
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
   normalized = smoothed(f>1.2 & f<25)./max_imp;
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

