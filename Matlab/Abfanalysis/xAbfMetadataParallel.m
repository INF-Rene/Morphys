%% bulk data extraction from Abffile objects
% When nr. of files is too large to be practical in one Abfbatch object, this function comes in handy. It chucks the list
% of files into smaller batches of files, keeping the total number of MBs per batch roughy constant. Then creates ABFBATCH 
% objects for them using parallel for looping (parfor) to obtain the metadata tables. Parfor allows you to use more of your 
% CPU
%   ---------------------------------------------------------------------------------------------------------------------
%   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
%   Created:      2017-08-18       
%   Modifications (Date/Name/Description):
%   ---------------------------------------------------------------------------------------------------------------------

%% presets
close all,clear all

% Global username must be defined!
global USERNAME
id = USERNAME;
if isempty(id), 
    s = Sharedpaths; % this is just to pop up the dialog box for username selection, so that there exists a 
    clear s
end

% set paths
dir_mats_analysed = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Natalia\Selection\onrap30\analyzed';% where the mat files of Abffile objects are stored
dir_batch_base    = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Natalia\Selection\onrap30';

% get list of filenames and make all paths
fileinfo  = dir(fullfile(dir_mats_analysed,'*.mat'));
filelist  = {fileinfo.name};
pathnames = arrayfun(@(x) fullfile(dir_mats_analysed,filelist{x}),1:numel(filelist),'UniformOutput',false)';

% make batches, each with an approximate megabyte limit
batchsizelimit = 500e6 ; % so that's the number of MBs per batch.
[batchnrs,~,indexes] = unique(floor(cumsum([fileinfo.bytes])/batchsizelimit));
batchlists = arrayfun(@(x) pathnames(indexes==x),batchnrs+1,'UniformOutput',false);

% storage variable to keep track
indexesbatched = zeros(numel(batchlists),1);

%% load abfs into batches first.  
parfor i=1:numel(batchlists)   
    
    % make the Abfbatch object and store
    fprintf('\nBatching loop nr %d, %d files in total:\n\n',i,numel(batchlists{i}))
    bb = Abfbatch(batchlists{i});
    try
        fprintf('Saving Abffile table...')
        p = fullfile(dir_batch_base,'Abffiles'); if ~isdir(p),mkdir(p); end
        bb.savetableas('abfs',fullfile(dir_batch_base,'Abffiles'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Channel table...')
        p = fullfile(dir_batch_base,'Channels'); if ~isdir(p),mkdir(p); end
        bb.savetableas('channels',fullfile(dir_batch_base,'Channels'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Analogout table...')
        p = fullfile(dir_batch_base,'Analogouts'); if ~isdir(p),mkdir(p); end
        bb.savetableas('analogouts',fullfile(dir_batch_base,'Analogouts'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Analogin table...')
        p = fullfile(dir_batch_base,'Analogins'); if ~isdir(p),mkdir(p); end
        bb.savetableas('analogins',fullfile(dir_batch_base,'Analogins'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Sweep table...')
        p = fullfile(dir_batch_base,'Sweeps'); if ~isdir(p),mkdir(p); end
        bb.savetableas('sweeps',fullfile(dir_batch_base,'Sweeps'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Epoch table...')
        p = fullfile(dir_batch_base,'Epochs'); if ~isdir(p),mkdir(p); end
        bb.savetableas('epochs',fullfile(dir_batch_base,'Epochs'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    try
        fprintf('Saving Actionpotential table...')
        p = fullfile(dir_batch_base,'Actionpotentials'); if ~isdir(p),mkdir(p); end
        bb.savetableas('aps',fullfile(dir_batch_base,'Actionpotentials'),'csv')
        fprintf('done.\n')
    catch
        fprintf('FAILED.\n')
    end
    %bb.saveme(fullfile(dir_batch_base,'Batches'))
    indexesbatched(i)=i;
end
