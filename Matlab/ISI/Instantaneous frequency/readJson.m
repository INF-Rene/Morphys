% This function extract information on AP kinetics features from json files
% and calculates instantaneous frequency for each spike for one cell

function [data_all, info]=readJson(file2analyse, cell_id, path2saveSpikes)

%define the bins for instantaneous frequency
edges=[0:10:100];  
% read json file 
str=fileread(file2analyse);
data_full = jsondecode(str);
data = data_full.feature_extraction.sweep_features;
sweeps_to_extract=data_full.sweep_extraction.sweep_features;

% get the information on the protocols and quality of the sweeps
for ll=1:size(data_full.qc.sweep_states,1)
    sweep_qc(ll,1)=data_full.qc.sweep_states(ll).sweep_number;
    sweep_qc(ll,2)=data_full.qc.sweep_states(ll).passed;
end

for l=1:size(sweeps_to_extract,1)
    if iscell(sweeps_to_extract)
    info(l).sweep_numbers=sweeps_to_extract{l}.sweep_number;
    info(l).stim_code=sweeps_to_extract{l}.stimulus_code;
    info(l).stim_name=sweeps_to_extract{l}.stimulus_name;
    elseif isstruct(sweeps_to_extract)
    info(l).sweep_numbers=sweeps_to_extract(l).sweep_number;
    info(l).stim_code=sweeps_to_extract(l).stimulus_code;
    info(l).stim_name=sweeps_to_extract(l).stimulus_name;   
    end
    idx=find(sweep_qc(:,1)==info(l).sweep_numbers);
    if ~isempty(idx)
        info(l).qc=sweep_qc(idx,2);     
    else
        info(l).qc=NaN;
    end
end

% extract data for sweeps only if protocol is "Long Square" and sweep
% quality check was "pass"
kk=1;
for l=1:size(sweeps_to_extract,1)
    if strcmp(info(l).stim_name,'Long Square') && info(l).qc==1
        sweeps(kk,1)= info(l).sweep_numbers;
        kk=kk+1;
    end
end

% get the names of the sweeps
sweepnames=fields(data);
varnames={'sweep_nr','spike_nr','peak_t','interval', 'inst_freq'...
    'upstroke','downstroke','peak_v','width', 'threshold','upstroke_downstroke_ratio', 'freq_bin'};
     
% get the data for all spikes in one table (data_all)
data_all=table();

for n= 1:size(sweepnames,1)
spikes=struct();
data_tb=table();
sweep_nr=[];
      spike_nr=[];
      inst_freq=[];
      peak_t=[];
      interval=[];
      inst_freq=[];
      upstroke=[];
      downstroke=[];
      peak_v=[];
      width=[];
      threshold=[];
      upstroke_downstroke_ratio=[];
sweep_nr=data.(sweepnames{n}).sweep_number;
% get the data for spikes
spikes=data.(sweepnames{n}).spikes;
% only spikes for sweeps with spikes and sweeps that passed qualty check
    if isfield(spikes,'peak_t') && any(sweeps==data.(sweepnames{n}).sweep_number)
 % get all the features
 feature={'peak_t','upstroke', 'downstroke','peak_v','width',...
     'threshold_v','upstroke_downstroke_ratio'};

        for i=1:size(spikes,1)
             features.sweep_nr(i,1)=data.(sweepnames{n}).sweep_number;
             features.spike_nr(i,1)=i;
            for kk=1:length(feature)
                % check if the feature exists, otherwise it gets 'NaN'
                if ~isempty(spikes(i). (feature{kk}))
                features.(feature{kk})(i,1)=spikes(i). (feature{kk});
                else 
                features.(feature{kk})(i,1)=NaN;
                end
            end
        end
      
      % calculate interval to the previous spike and instantaneous
      % frequency
      features.interval(1,1)=0;
       features.inst_freq(1,1)=0;
          for nn=2:size( features.peak_t,1)
               features.interval(nn,1)= features.peak_t(nn)- features.peak_t(nn-1);
               features.inst_freq(nn,1)=1/ features.interval(nn);
          end
      % assign each spike a bin for instantaneous frequency bin 1 (0-10Hz), 2(10-20Hz) etc   
       features.freq_bin = discretize( features.inst_freq,edges);
      % make a table with all features
      data_tb=table(features.sweep_nr, features.spike_nr, features.peak_t,...
           features.interval, features.inst_freq,...
       features.upstroke, features.downstroke, features.peak_v, features.width,...
        features.threshold_v,  features.upstroke_downstroke_ratio,  features.freq_bin,...
      'VariableNames', varnames);
    end
    % combine tabels for different sweeps. 
    data_all=[data_all;data_tb]; 
end
name=sprintf('%s',cell_id,'_spikes.mat');
file2save=fullfile(path2saveSpikes,name);
save(file2save, 'data_all', 'info');
end




