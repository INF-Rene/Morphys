%% Resonance analysis  single files (raw nwb files)
fn= '/Volumes/Ephys2/Ephys hipp/healthy/cluster1/H21.29.198.21.11.04.nwb'
nwb = NWBfile(fn,[{'CHIRP'}]);
%%
 signal = nwb.getstimset.getnwbchannel.getsweep(1:5);
 plot(signal)
%%
signal = nwb.getstimset.getnwbchannel.getsweep(1:5).avtrace;
plot(signal)
signal = signal.low_pass_filter(1000);
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;
%signal = signal.Data ; 
stim= nwb.getstimset.getnwbchannel.getstimwave(2).Data ;
%%
  Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;
   
   signal_fft=fft(signal);
    stim_fft=fft(stim);

        %smoothed
    smoothed=smooth(abs((signal_fft)./abs(stim_fft)),23,  'lowess');
    plot(f, smoothed);
    xlim([1 25])
    
    %%
    clear