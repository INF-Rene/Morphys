%%Load in analyzed mat files 
basedir = '/Users/elinemertens/Data/ephys/Hippocampus/nwb2 analyzed/198' ;
savedir = '/Users/elinemertens/Data/ephys/Summary/Human' ;
savename = 'Summary' ;

%% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};

%% settings
filter_freq = 1000; % low pass filter frequency for sag traces (Hz)

%% loop over files
%t_all=table();
table_initialize = false ;

for i=1:numel(filelist)
    %load file
    fn = fullfile(basedir, filelist{i});
    nwb = load(fn);
    nwb = nwb.obj;
    
    % loop over stimsets
    for j = 1:nwb.nrofstimsets
        stimset = nwb.getstimset(j);
%         stimset_filter = nwb.getstimset.name
%         ({'Sub'} {'CCSteps_DA_0'});
        swps = stimset.getnwbchannel.getsweep();
        
        % extract from each sweep:
        % step size (pA)
        stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
        pA = [swps.getepoch(stepepoch).amplitude];
        
        %select only negative injection swps
        if all(pA>=0), continue;end
        swps = swps(pA<0);
        pA = pA(pA<0);
        
    
      %  swptime = swps.Time ;
      %  if all(swptime<300000), continue;end
        
        for k = 1:numel(swps)
            loc = find(isnan(swps(k).Data), 1); %find from where swp is NaN
            if ~isempty(loc) ,swps(k) = swps(k).getsamples(1:loc-1);end
        end
        
        % peak and ss voltage
        
%         swps_f = swps.low_pass_filter(1000);
        
        vstep = [swps.getepoch(stepepoch).vstep];
        rmp = [swps.getepoch(stepepoch-1).steadystate];
        v_delta = vstep - rmp;
        ss = [swps.getepoch(stepepoch).steadystate];
        ss_delta = ss - rmp;
        delta = ss - vstep;
        sag_ratio = delta ./ -v_delta;
        
        % peak time
        
        % time constant of sag activation
        
        % rebound amplitude
%         reb_start_time = swps(1).getepoch(stepepoch).Time(end);
%         reb_ss_time = swps(1).Time(end);
%         rebound_window=300; % in ms
%         reb_traces = swps.getsampleusingtime(reb_start_time, reb_start_time+rebound_window);
%         reb_volts = NaN(1, numel(reb_traces));
%         reb_ss=reb_volts;
%         for k = 1:numel(reb_traces)
%             reb_volts(k)=max(reb_traces(k).medianfilter(0.5,'truncate').getdata);
%             reb_ss(k) = swps(k).getsampleusingtime(reb_ss_time-200, reb_ss_time).median;
%         end
%         reb_amps = reb_volts-reb_ss;
        
        % results to table
     % Protocol_column = repmat(stimset.name, height(t), 1);
     % Filename_column = repmat(filelist{i}, height(t), 1); 
      
        t = table(pA', vstep', rmp', ss', ss_delta', sag_ratio', 'VariableNames', {'pA', 'sagvolt', 'RMP', 'ss', 'ss_delta', 'sag_ratio'});
      %  t.Protocol = repmat(stimset.name, height(t), 1);
      t.Filename = repmat(filelist{i}, height(t), 1);
      t.Protocol = repmat({stimset.name},height(t),1) %add a ; if you dont want it printed in command
        
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
        
        
    end
    
end