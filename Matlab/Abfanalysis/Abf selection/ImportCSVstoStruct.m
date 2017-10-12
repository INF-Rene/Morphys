function abf = ImportCSVstoStruct(mainfolder, savename)
%ImportCSVstoStruct import csv metadatatables created by
%   xAbfMetadataParallel script. convert tables to structures and restore
%   nested hierarchical order.
%
%   mainfolder is the path to the folder containing the necessary
%   subfolders: Abffiles, Channels, Analogins, Sweeps, Epochs and
%   Actionpotentials. which all contain the corresponding CSV file with similar name.
%   for instance: mainfolder contains a subfolder Abffiles which contains
%   Abffiles.csv

%   savename is the filename of the .mat file containing the resulting
%   struct to be saved

%   Written by Djai Heyer

%% import CSV files
abfs = table2struct(readtable(fullfile(mainfolder, 'Abffiles','Abffiles.txt'))) ;
channels = table2struct(readtable(fullfile(mainfolder, 'Channels','Channels.txt'))) ;
ins = table2struct(readtable(fullfile(mainfolder, 'Analogins','Analogins.txt'))) ;
sweeps = table2struct(readtable(fullfile(mainfolder, 'Sweeps','Sweeps.txt'))) ;
epochs = table2struct(readtable(fullfile(mainfolder, 'Epochs','Epochs.txt'))) ;
aps = table2struct(readtable(fullfile(mainfolder, 'Actionpotentials','Actionpotentials.txt'))) ;

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
end

for i = 1:length(abfs)    
        abfs(i).channel = channels(ismember({channels.guid_abf},abfs(i).guid)) ;    
end

abf = abfs ;
%% Save struct
save(fullfile(mainfolder, savename), 'abf') ;

end




