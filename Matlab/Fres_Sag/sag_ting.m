p1   = ('/Users/elinemertens/Data/ephys/Hippocampus/Ting');
fn1  = 'Cell6_ccsteps.abf';
fp1  = fullfile(p1,fn1);
%ss  = load('/Users/elinemertens/Downloads/Morphys-master_2022/Setupsettings_Allen.mat');
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

%%
% nwb swps = stimset.getnwbchannel.getsweep();
swps =  abf.channels(1, 1).analogins.sweeps(1:4) ;
swps.plot ; 
swps = swps.analysesweep ; 
swps.plot

swps.plotanalysis ; 

%%
 table_initialize = false ;
 
for i=1:numel(swps)
   stepepoch = find([swps(1).getepoch(5).amplitude]~=0, 1);
        pA = [swps.getepoch(5).amplitude];
         vstep = [swps.getepoch(5).vstep];
        rmp = [swps.getepoch(4).steadystate];
        v_delta = vstep - rmp;
        ss = [swps.getepoch(5).steadystate];
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

%% 
writetable(t, 'summary_sag_Ting.xlsx'); 
