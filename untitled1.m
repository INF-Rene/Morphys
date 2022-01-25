%checking for the reason of error int he file
clear all
clc
fn = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/H19.29.162.11.11.03.nwb'

dir_fixed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/fixed'
%stimsets for APs
nwb = NWBfile(fn,[{'X2LP_Search'} {'teps'}])


%bulk conversion 
a = NWBfile(fullfile(dir_abfs,fnabf), [{'hresh'} {'teps'} {'LSFINEST'} {'LSCOARSE'}])

%all stimsets

b = NWBfile(fn, {})

%what I need
c = NWBfile(fn, [{'teps'} {'hres'} {'Search'} {'Rheo'}])

X1PS_SubThresh
X2LP_Search_DA_0
X3LP_Rheo_DA_0
X4PS_SupraThresh_DA_0
X4PT_C2NSD1SHORT_DA_0 
X4PU_C2NSD2SHORT_DA_0
X5SP_Search_DA_0
X6SP_Rheo_DA_0
X6SQ_C2SSTRIPLE_DA_0
X7_Ramp_DA_0
X8_CHIRP_DA_
X9_C1QCAPCHK_

a = c.analyseNWB

fnmat = [fn '.mat'];
a.saveme(fnmat)

%%
clear all
clc

fn = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/H19.29.162.11.11.03.nwb'


%% settings
filter_freq = 1000; % low pass filter frequency for sag traces (Hz)

%% loop over files

table_initialize = false ;
nwb = NWBfile(fn,[{'X1PS_SubT'} {'teps'}])
swps = nwb.stimsets.nwbchannels.sweeps();
stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
pA = [swps.getepoch(stepepoch).amplitude];
for k = 1:numel(swps)
            loc = find(isnan(swps(k).Data), 1); %find from where swp is NaN
            if ~isempty(loc) ,swps(k) = swps(k).getsamples(1:loc-1);end
        end
    
        vstep = [swps.getepoch(stepepoch).vstep];
        rmp = [swps.getepoch(stepepoch-1).steadystate];
        v_delta = vstep - rmp;
        ss = [swps.getepoch(stepepoch).steadystate];
        ss_delta = ss - rmp;
        delta = ss - vstep;
        sag_ratio = delta ./ -v_delta;
        
t = table(pA', vstep', rmp', ss', ss_delta', delta', sag_ratio', 'VariableNames', {'pA', 'sagvolt', 'RMP', 'ss', 'ss_delta', 'delta', 'sag_ratio'});
      %  t.Protocol = repmat(stimset.name, height(t), 1);
      t.Filename = repmat({filelist{i}}, height(t), 1);
      t.Protocol = repmat({stimset.name},height(t),1) %add a ; if you dont want it printed in command
        
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
        %Anna's edition 2021-Dec-1
        writetable(t_all,'sag_ratios.xlsx');%if you add an extra argument "Sheet",
        %1/2/etc - then you can add multiplie MAtLab tables into one xlsx
    end
    catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end        
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages=[errorMessages; filelist{i} {ME}];%append the error message to a list
        end
    end
end
save('/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/fixed/errorMessages_011221.mat','errorMessages');
