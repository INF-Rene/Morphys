
close all, clear all

%% Load NWB file
fn= 'C:\Users\tamar\Documents\Master\Internship I\Patch data\Human\H219\H22.29.219.11.72.03-compressed.nwb'
%% 
%%nwb = NWBfile(fn)
%% Get triple protocol
nwb = NWBfile(fn, [{'X6SQ_C2SSTRIPLE_DA_0'}])

%% Plot all sweeps in different windows
f1 = figure;
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,1));
f2 = figure;
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,2));
f3 = figure;
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,3));
f4 = figure;
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,4));
f5 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,5));
f6 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,6));
f7 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,7));
f8 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,8));
f9 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,9));
f10 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,10));
f11 = figure
plot(nwb.stimsets(1,1).nwbchannels.sweeps(1,11));








