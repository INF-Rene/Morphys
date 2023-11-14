 signal = obj.stimsets(1).nwbchannels.sweeps(1) ; 
 signal = signal.Data ;
 %Fs = 1000
%signal = signal.low_pass_filter(1000);
 signal = lowpass(signal,1000,10000);
 plot(signal) ; 
 hold on
 ylim([-75 -68]) ; 