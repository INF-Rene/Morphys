%% Resonance analysis single files (raw nwb files): first select only correct protocol
% We use the script NWBfile to select this 
fn= '/Volumes/Expansion/Ephys 173-224/H22.29.217_F/H22.29.217.11.01.04-compressed.nwb';
nwb = NWBfile(fn,[{'CHIRP'}]) ; 
%% Look at the sweeps, are they comparable? 
 signal = nwb.getstimset.getnwbchannel.getsweep([2]);
 plot(signal)
%% Only average the sweeps that look comparable
signal = nwb.getstimset.getnwbchannel.getsweep([1:3]).avtrace;
plot(signal)
stim = nwb.getstimset.getnwbchannel.getstimwave(3).Data ;
%signal = signal.low_pass_filter(1000);
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;
%signal = lowpass(signal,1000,10000);
    Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;        
    plot(signal)
    grid off      
    %% Fast fourier transformation: plot only smoothed signal 
     signal_fft=fft(signal);
    stim_fft=fft(stim);
figure(4)
 smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([1 25])
    %ylim([0 100])
    grid off  ;
set(gca, 'TickDir', 'out')
    box off         
                 ylabel('Impedance (MΩ)','FontSize',12)
              xlabel('Frequency (Hz)', 'FontSize',12)
             
    [max_imp, loc] = max(smoothed(f>1.1 & f<25));
    f2 = f(f>1.1 & f<25);
    res_freq = f2(loc);
    
   normalized = smoothed(f>1.1 & f<25)./max_imp;
   
   % when does normalized frequency response go below square root of 0.5?
   % This is the definition of 3 dB cutoff!
   cutoff_3db=f2(find(normalized>=sqrt(0.5), 1, 'last'));
   figure (2)
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

%%  
plot(stim) 
hold on 
ylim([-300,300])
