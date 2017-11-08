function abf = SubsetTable2struct(abfs,channels,ins,outs,sweeps,epochs,aps)
% returns nested struct containing subset of whichever abfs are entered as
% argument.

% example: >>> abf = SubsetTable2struct(abfs,channels,ins,outs,sweeps,epochs,aps) 
%           : returns nested structure of all abfs in 'abfs'table.

% example: >>> abf = SubsetTable2struct(abfs(1,:),channels,ins,outs,sweeps,epochs,aps) 
%           : returns nested structure of first abf in 'abfs'table.

%   Written by Djai Heyer

%% Get subsets
channels = channels(ismember(channels.guid_abf,abfs.guid),:) ;
ins = ins(ismember(ins.guid_channel,channels.guid),:) ;
outs = outs(ismember(outs.guid_channel,channels.guid),:) ;
sweeps = sweeps(ismember(sweeps.guid_in,ins.guid),:) ;
epochs = epochs(ismember(epochs.guid_swp,sweeps.guid),:) ;
aps = aps(ismember(aps.parent_guid,epochs.guid),:) ;

%% Tables to structs
abfs = table2struct(abfs) ;
channels = table2struct(channels) ;
ins = table2struct(ins) ;
outs = table2struct(outs) ;
sweeps = table2struct(sweeps) ;
epochs = table2struct(epochs) ;
aps = table2struct(aps) ;

%% Create nested structure
for i = 1:length(epochs)    
        epochs(i).ap = aps(ismember({aps.parent_guid},epochs(i).guid)) ;    
end

for i = 1:length(sweeps)    
        sweeps(i).epoch = epochs(ismember({epochs.guid_swp},sweeps(i).guid)) ;    
end

for i = 1:length(ins)    
        ins(i).sweep = sweeps(ismember({sweeps.guid_in},ins(i).guid)) ;    
end

for i = 1:length(channels)    
        channels(i).in = ins(ismember({ins.guid_channel},channels(i).guid)) ;
        channels(i).out = outs(ismember({outs.guid_channel},channels(i).guid)) ;
end

for i = 1:length(abfs)    
        abfs(i).channel = channels(ismember({channels.guid_abf},abfs(i).guid)) ;    
end

abf = abfs ;
end




