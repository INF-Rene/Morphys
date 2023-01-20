%% Resonance analysis  single files (raw nwb files)
fn= '/Volumes/Ephys2/Ephys hipp/healthy/cluster2/H21.29.198.21.01.01.nwb'
nwb = NWBfile(fn,[{'CHIRP'}]);
%%
 signal = nwb.getstimset.getnwbchannel.getsweep(2:4);
 plot(signal)
%%
signal = nwb.getstimset.getnwbchannel.getsweep(2:4).avtrace;
plot(signal)
signal = signal.low_pass_filter(1000);
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;
time = nwb.getstimset.getnwbchannel.getsweep(2) ;
time = time.low_pass_filter(1000) ; 
time=time.resample(0:0.1:time.TimeInfo.End).Time;
%%
clear