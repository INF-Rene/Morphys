p1   = ('/Users/elinemertens/Data/ephys/Hippocampus/H17.29.117.21/H17.29.117.21_all_thijs');
fn1  = '2017_03_29_0041.abf';
fp1  = fullfile(p1,fn1);
ss  = load('/Users/elinemertens/Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf = Abffile(fp1,ss);
abf.plot
%%
abf = abf.analyseabf ; 

%% if you want to specify the correct channels you do this 
%   abf.channels = abf.channels(1, 2) ; 
%   abf.channels.analogins.sweeps(1:8).plot
% De pA klopt niet, deze wordt nu uit het output kanaal gehaald, 
% indien van belang opnieuw runnen en juiste kanaal pakken

%%
filter_freq = 1000;

%%  when you don't run the part up here, and you keep all channels, you run this 
 abf.channels(1, 2).analogins.sweeps(1:8).plot
 
 %% Djai only run this and the loop 
 abf.channels.analogins(1, 1).sweeps(1:8).plot ; 
 swps = abf.channels.analogins(1, 1).sweeps(1:8) ; 
 

%%
% nwb swps = stimset.getnwbchannel.getsweep();
swps = abf.channels(1, 2).analogins(1, 1).sweeps(1:8) ;

%%
swp1 = abf.channels(1, 2).analogins(1, 1).sweeps(1) ; 
%s = swp1.analysesweep ; 
close all ; 
swps.plotanalysis

%%
 table_initialize = false ;
 
for i=1:numel(swps)
   stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
        pA = [swps.getepoch(stepepoch).amplitude];
         vstep = [swps.getepoch(stepepoch).vstep];
        rmp = [swps.getepoch(stepepoch-1).steadystate];
        v_delta = vstep - rmp;
        ss = [swps.getepoch(stepepoch).steadystate];
        ss_delta = ss - rmp;
        delta = ss - vstep;
        sag_ratio = delta ./ -v_delta;
end
        

%%
      t = table(pA', rmp', vstep', ss', ss_delta', sag_ratio', 'VariableNames', {'pA', 'RMP', 'sagvolt', 'ss', 'ss_delta', 'sag_ratio'});
      %  t.Protocol = repmat(stimset.name, height(t), 1);
      t.Filename = repmat(fn1, height(t), 1);
      %t.Protocol = repmat({stimset.name},height(t),1) %add a ; if you dont want it printed in command
        
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
        
        
%%

stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
pA = [swps.getepoch(stepepoch).amplitude];