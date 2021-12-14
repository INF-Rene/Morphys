%% Resonance analysis  single files (raw nwb files)

fn= '/Users/elinemertens/Data/ephys/Hippocampus/H21.29.198.21/H21.29.198.21.11.01.nwb'

nwb = NWBfile(fn,[{'X8_C'}])


%% Zorg dat hiervoor de traces individueel al bekeken zijn of ze vergelijkbaar zijn 
%(indien niet, geen average trace nemen)
% Bij errors 'expected input to be finite' moet je 
% even de losse sweeps bekijken en dan de juiste sweeps selecteren
% signal = nwb.getstimset.getnwbchannel.getsweep(1:3).avtrace
% or signal = nwb.getstimset.getnwbchannel.getsweep([1:4 5:6]).avtrace;
 signal = nwb.getstimset.getnwbchannel.getsweep(3).avtrace;
 
 plot(signal)
 
 %looks good? combine them 

%%
signal = nwb.getstimset.getnwbchannel.getsweep.avtrace;

plot(signal)
%%
signal = signal.low_pass_filter(1000);
%check if you take correct stimwave(#) (after wash in this is diff)
stim = nwb.getstimset.getnwbchannel.getstimwave(2).Data ;
% if the first sweep is bad or incomplete, take getstimwave(2)
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;


%%
%   signal=normalize(mean(squeeze(data(:,1,:))'), 'range');
%     stim=normalize(mean(squeeze(data(:,2,:))'), 'range');

    %%
    Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;
   
    plot(signal)
    grid off
set(gca, 'TickDir', 'out')
    ylabel('Voltage (mV)','FontSize',12)
              xlabel('Time (sec)', 'FontSize',12)
   
 %% Values resonance 
 
time = nwb.getstimset.getnwbchannel.getstimwave(2).Time ;  
 aastart_rmp = mean(signal(time>100 & time<500));

[first_wave, loc] = max(signal(time>500 & time<1500));
 adelta_first_wave = (aastart_rmp) - (first_wave) ; 

[max_ampl, loc] = max(signal(time>1500 & time<20000));
t2 = time(time>1500 & time<20000);
  time_max_amp = t2(loc);
  cdelta_max_wave = (aastart_rmp) - (max_ampl) ;

end_rmp = mean(signal(time>21000 & time<23000)) ;
          
    
    %%
    signal_fft=fft(signal);
    stim_fft=fft(stim);

    %plot

%     %raw
    plot(f, abs((signal_fft)./abs(stim_fft)*1000));
    xlim([0.5 25])
    hold on
        %smoothed
    smoothed=smooth(abs((signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 25])
    
    %% plot only smoothed 

     smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 25])
    grid off  ;
set(gca, 'TickDir', 'out')
    box off         
                 ylabel('Impedance (MΩ)','FontSize',12)
              xlabel('Frequency (Hz)', 'FontSize',12)
    
   %% resonance freq
    [max_imp, loc] = max(smoothed(f>0.5 & f<25));
    f2 = f(f>0.5 & f<25);
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
 box off
 set(gca, 'TickDir', 'out')
 ylabel('Impedance (normalized)','FontSize',12)
 xlabel('Frequency (Hz)', 'FontSize',12)
   
%%
clear all 
