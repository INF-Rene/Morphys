

%plot(nwb.stimsets(1, 1).nwbchannels.stimwaves)
obj.stimsets(1, 3).nwbchannels.sweeps.plot   
grid off
set(gca, 'TickDir', 'out')
set(0, 'DefaultFigureRenderer', 'painters');