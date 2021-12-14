clear all
clc
%% Bulk conversion for NWB2 files
% set path to find .NWB files
dir_abfs          = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/test for error handling';
% set path to save analyzed files
dir_mats_analysed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human/test for error handling/mat'
% set path to save figures
dir_figs_analysed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/PatchSeq analysed/HUMAN/NOT WORKING/new nwb2/figures'
% get a list of all .nwb files in the specified path
fileinfo  = dir(fullfile(dir_abfs,'*.nwb'));
filelist  = {fileinfo.name};


%% this is optional: exclude files from the list that you have already analyzed
fileinfo2  = dir(fullfile(dir_mats_analysed,'*.mat'));
filelist2  = {fileinfo2.name};
filelist2= cellfun(@(x) [x(1:end-4) '.nwb'],filelist2,'UniformOutput',false) ;
filelist = filelist(~ismember(filelist,filelist2)) ;

%% Convert and analyze
% loop through filelist
errorMessages=[];%keep a list of any errors, Anna 01-Dec
verbose=1%keep a list of any errors, Anna 01-Dec
for i=1:size(filelist,2)
    try%try to process data for a given session
        
    % add the right extensions to the filename for saving and loading
    fn    = filelist{i}(1:end-4);
    fnabf = [fn '.nwb'];
    fnmat = [fn '.mat'];
    
    % load file and convert to NWBfile object. 
    % use the final argument to specify which stimsets to include
    %a = NWBfile(fullfile(dir_abfs,fnabf), [{'LP'} {'hresh'} {'LSFINEST'} {'LSCOARSE'}]) ;
    a = NWBfile(fullfile(dir_abfs,fnabf), [{'hresh'} {'teps'} {'LSFINEST'} {'LSCOARSE'}]) ;
    %a = NWBfile(fullfile(dir_abfs,fnabf), [{'Rheo'} {'CCSteps'}]) ;
    % analyse NWBfile object
    a = a.analyseNWB ;
    % save the object to .mat file to be used for further analysis
    a.saveme(dir_mats_analysed,fnmat);
    %a.plotanalysis(dir_figs_analysed)
    
    %Anna_edition 1-Dec - handling errors
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

%General end error message - Anna 1-Dec
if ~isempty(errorMessages)
    for i=1:size(errorMessages,1)
        fprintf('\nError at session %d\n',errorMessages{i,1});%print problematic session number to screen
        errorMessages{i,2}.message%display error message
        errorMessages{i,2}.stack%print function and line number at which error occurred
    end
end
    
    
%% Get stimset lists 
% set path to find .NWB files
dir_abfs          = 'C:\Users\DBHeyer\Documents\PhD\Human Database\nwb2\2019\H19.29.158\H19.29.158.11.11';
% set path to save analyzed files
dir_mats_analysed = 'C:\Users\DBHeyer\Documents\PhD\Human Database\nwb2\analyzed';

% get a list of all .nwb files in the specified path
fileinfo  = dir(fullfile(dir_abfs,'*.nwb'));
filelist  = {fileinfo.name};

% Create struct
stimsets = struct();
% loop through filelist
for i=1:size(filelist,2) 
    % get filename
    fn    = filelist{i};
    % Get nwb file structure
    info=h5info(fullfile(dir_abfs,fn));  
    % Get sweep names
    grloc = strcmp({info.Groups.Name}, '/acquisition');
    swps={info.Groups(grloc).Groups.Name};
    % Get protocol/stimset names for all sweeps
    protocols=cell(numel(swps), 1);
    for j=1:numel(swps)
        protocols(j) = h5readatt(fullfile(dir_abfs,fn), [swps{j}], 'stimulus_description');
    end
    % remove duplicates
    protocols = unique(protocols);
    % Add list to struct
    stimsets(i).file = fn;
    stimsets(i).protocols = protocols;
end        
    
% Save struct
save(fullfile(dir_mats_analysed, 'AllStimsets'), 'stimsets') ;    
    

  