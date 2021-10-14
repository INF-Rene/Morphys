   
%plot((obj.stimsets(1, 1).nwbchannels.sweeps(1,1:27)), 'k') 


%als je alleen -150pA wil dan neem je stimwave 1 tot met 30 
plot((obj.stimsets(1, 1).nwbchannels.stimwaves(1,1:27)), 'k')


%stel je wilt enkel de laatste sweeps (met nieuwe scaling dan doe je dit:
%plot((obj.stimsets(1, 1).nwbchannels.stimwaves(1,31:end)), 'k')
%plot((obj.stimsets(1, 1).nwbchannels.sweeps(1,31:end)), 'k')

%als je alles wilt doe je dit 
%plot((obj.stimsets(1, 1).nwbchannels.stimwaves), 'k')  

                xlim([0 1800])
                title('Sag')
                grid off
                 set(gca, 'TickDir', 'out')
                 ylabel('Membrane voltage (mV)')
              xlabel('Time (ms)')