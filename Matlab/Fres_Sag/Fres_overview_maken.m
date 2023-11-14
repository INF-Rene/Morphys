%% Resonance analysis  single files (raw nwb files)

fn= '/Users/elinemertens/Data/Projects/NEW_Hippocampus_2022_07/Data/mouse/M23.29.45992.11.13.nwb';

nwb = NWBfile(fn,[{'CHIRP'}]);
 
%%
%signal = nwb.getstimset.getnwbchannel.getsweep([1:2]).avtrace;
signal = nwb.stimsets(1, 1).nwbchannels.sweeps([3:4 6 8 10]).avtrace ;   

plot(signal)
%check if you take correct stimwave(#) (after wash in this is diff)
stim = nwb.getstimset.getnwbchannel.getstimwave(2).Data ;
signal = signal.low_pass_filter(1000);
% if the first sweep is bad or incomplete, take getstimwave(2)
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;

    Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;

    %%
    signal_fft=fft(signal);
    stim_fft=fft(stim);
   
figure(4)
     smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 15])
    ylim([80 350])
    grid off  ;
set(gca, 'TickDir', 'out')
    box off         
                 ylabel('Impedance (MΩ)','FontSize',12)
              xlabel('Frequency (Hz)', 'FontSize',12)

              
    [max_imp, loc] = max(smoothed(f>0.5 & f<25));
    f2 = f(f>0.5 & f<25);
    res_freq = f2(loc);
    
    % Values resonance 
 
time = nwb.getstimset.getnwbchannel.getstimwave(2).Time ;  
 aastart_rmp = mean(signal(time>100 & time<500));

[max_ampl, loc] = max(signal(time>1500 & time<20000));
t2 = time(time>1500 & time<20000);
  time_max_amp = t2(loc);
  cdelta_max_wave = (aastart_rmp) - (max_ampl) ;
  
  normalized = smoothed(f>1.1 & f<25)./max_imp;
   
   % when does normalized frequency response go below square root of 0.5?
   % This is the definition of 3 dB cutoff!
   cutoff_3db=f2(find(normalized>=sqrt(0.5), 1, 'last'));

%%
%   signal=normalize(mean(squeeze(data(:,1,:))'), 'range');
%     stim=normalize(mean(squeeze(data(:,2,:))'), 'range');
