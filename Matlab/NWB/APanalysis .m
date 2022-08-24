%%Load in analyzed mat files 
basedir = '/Users/elinemertens/Data/ephys/Hippocampus/nwb2_analyzed/Analysis AP' ;
savedir = '/Users/elinemertens/Data/ephys/Summary/Human' ;
savename = 'Summary_209' ;

%% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};


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
        % select the one with cc step 
        stimset = nwb.getstimset(1);
       % stimset = nwb.getstimset.name({'CC'}, {'CCSteps_DA_0'});
        swps = stimset.getnwbchannel.getsweep();
        
        % extract from each sweep:
        % step size (pA)
        stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
        pA = [swps.getepoch(stepepoch).amplitude];
        
        %select only positive injection swps
        if all(pA<=0), continue;end
        swps = swps(pA>0);
        pA = pA(pA>0);
        
    
      %  swptime = swps.Time ;
      %  if all(swptime<300000), continue;end
        
        for k = 1:numel(swps)
            loc = find(isnan(swps(k).Data), 1); %find from where swp is NaN
            if ~isempty(loc) ,swps(k) = swps(k).getsamples(1:loc-1);end
        end
        
        % peak and ss voltage
        
%         swps_f = swps.low_pass_filter(1000);
    nroAPS = {};
    nrofAPS = [nroAPS, swps.nrofaps] ; 
    sweepname = swps.Name ; 
%     currinj = {}
%     currinj = [currinj, swps.epochs.amplitude] ; 
%      
%      
%        
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
      t = table(nrofAPS', pA', 'VariableNames', {'nrofAPS', 'pA'}); 
       t.Protocol_column = repmat(stimset.name, height(t), 1);
       t.Filename_column = repmat(filelist{i}, height(t), 1); 
%       % let op dat alle namen hetzelfde zijn (dus geen compressed.nwb), 
      % dit snapt het script niet 
      
    
      
      
%       index = 1 ; 
% for i = 1:length(swps)
% %        Summary(index).File               = fileinfo.name ;
% %         Summary(index).Date               = nwb.filetimestart ;
%          Summary(index).stimset             = stimset.name ;
%          Summary(index).sweepname    = sweepname ;
%         Summary(index).nrofAPs = nrofAPS ;
%         Summary(index).pA = pA ; 
%        index =  index + 1 ;
        
%           
%         t = arrayfun(nrofAPS', 'VariableNames', {'nrofAPS'});
%        t = arrayfun(nrofAPS', sweepname', currinj', 'VariableNames', {'nrofAPS', 'sweepname', 'currinj'});
%        %  t.Protocol = repmat(stimset.name, height(t), 1);
%       t.Filename = repmat(filelist{i}, height(t), 1);
%       t.Protocol = repmat({stimset.name},height(t),1) %add a ; if you dont want it printed in command
%         
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
        
        
    end
    
    end

%end