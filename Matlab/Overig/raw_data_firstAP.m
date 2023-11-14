Data = (obj.stimsets(1).nwbchannels.sweeps(16).epochs(3).Data)   ; 
Time = (obj.stimsets(1).nwbchannels.sweeps(16).epochs(3).Time) ;
Data = Data.low_pass_filter(50000);
Time = Time.low_pass_filter(50000);
plot(Time, Data);

%%
testAP = table(Data, Time) ;
 writetable(testAP, 'H22.29.209.21.61.02.nwb'.xlsx');







%%
Data3= (obj.stimsets(1, 1).nwbchannels.sweeps(1, 16).aps.ts.Data) ; 
%%
'H21.29.198.21.11.02.nwb'
ts.Data:1
ts.Time