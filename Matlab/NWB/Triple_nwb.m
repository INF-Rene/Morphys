
close all, clear all

%% Load NWB file
fn= '/Users/elinemertens/Data/ephys/Hippocampus/H21.29.198.21/H21.29.198.21.01.01.nwb'
%% 
%%nwb = NWBfile(fn)
%% Get triple protocol
nwb = NWBfile(fn, [{'X6SQ_C2SSTRIPLE_DA_0'}])
obj=nwb.analyseNWB


%% need to figure this out 
figure(1)
        for j = 1:length(obj.stimsets.nwbchannels.sweeps)
          if obj.stimsets.nwbchannels.sweeps(j).aps > 2 
                obj.getstimset.getnwbchannel.getsweep.getepoch.getap(3).plot('superimpose','peak');
                legend(filelist)
                xlim([-5 10])
              %  title('First AP')
            %   ylabel('mV')
            %    xlabel('ms')
      end
        end


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








