
signal_stim = abfload_pro('/Users/elinemertens/Desktop/ting/Cell5_Fres3.abf');
signal = signal_stim(:,1);
stim = signal_stim(:,2);
signal = signal.low_pass_filter(1000);
stim = squeeze(data(:,2,1));
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
   % xlim([1.2 25])
    hold on
    
    
        %smoothed
    smoothed=smooth(abs((signal_fft)./abs(stim_fft)),23,  'lowess');
  %  plot(f, smoothed);
    xlim([1.2 15])
    
    %% plot only smoothed 
    
    smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    %xlim([1.2 16])
    grid off   
            
 ylabel('MOhm')
xlabel('Frequency (Hz)')
set(gca, 'TickDir', 'out')
    % Set the remaining axes properties

   [max_imp, loc] = max(smoothed(f>1.2 & f<16));
    f2 = f(f>1.2 & f<15);
    res_freq = f2(loc);
    
   %% 3db cutoff
   normalized = smoothed(f>1.2 & f<15)./max_imp;
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








%%

for i=1:numel(fn)
    data=abfload_pro(fn(i));
 
   
    signal=mean(squeeze(data(:,1,:))');
    stim=mean(squeeze(data(:,2,:))');

     Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;

    signal_fft=fft(signal);
    stim_fft=fft(stim);

    %plot

%     %raw
    plot(f, abs(signal_fft)./abs(stim_fft));
    xlim([0.5 35])
    hold on
        %smoothed
    smoothed=smooth(abs(signal_fft)./abs(stim_fft),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 35])
    
    %resonance freq
    [~, loc] = max(smoothed(f>0.5 & f<35));
    f2 = f(f>0.5 & f<35);
    res_freq = f2(loc);

end