%% 
%%how to load specific stimsets of the file - stimsets containing APs (and
%%supposed to contain APs)
clear all
clc

%fn = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/Test smth that workes all/H21.29.201.11.11.01.nwb'

%nwb = NWBfile(fn,[{'X2LP_Search'} {'teps'}])
basedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human' ;
savedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/ephys human' ;
fileinfo  = dir(fullfile(basedir,'*.nwb'));
filelist  = {fileinfo.name};
table_initialize = false ; % Anna 2-Dec-2021
errorMessages3=[];%keep a list of any errors, Anna 01-Dec
verbose=1%keep a list of any errors, Anna 01-Dec
for i=1:numel(filelist)
    try%try to process data for a given session
    fn = fullfile(basedir, filelist{i});
    nwb = NWBfile(fn, {'X2LP_Search','teps'});
    %getting the sweeps I need
    swps = nwb.getstimset.getnwbchannel.getsweep;
    %analysis of the sweeps
    swps = swps.analysesweep;
    % getting th eAPs from the sweeps which have them
    APs = swps.getap;
    %%getting the frequency of every AP
    inst_freqs = [APs.freq];
    % dividing the APS into the bins to get the ones we need
    binedges = [0,0.00001, 10, 20, 30, 40, 500];
    bins = discretize(inst_freqs, binedges);
    binlabels = {'FirstAP', 's0to10 Hz', 's10to20 Hz','s20to30 Hz' 's30to40 Hz', 's40toINF Hz'};
    
    upstrokes = [APs.maxdvdt];
    
    %scatter(categorical(binlabels(bins)), upstrokes)
    %choose specific APs and get data per bin
    mean_binned_upstrokes=NaN(1, numel(binlabels));
    median_binned_upstrokes=NaN(1, numel(binlabels));
    for j=1:numel(binlabels)
        mean_binned_upstrokes(j) = nanmean(upstrokes(bins==j) )  ;
        median_binned_upstrokes(j) = nanmedian(upstrokes(bins==j) )  ;
    end
    %boxplot(upstrokes, bins)
    %hold on
    %plot(swps)
    %scatter(categorical(binlabels(bins)), upstrokes)
    Mean_upstroke = array2table(mean_binned_upstrokes, 'VariableNames', {'M_FirstAP', 'M_s0to10 Hz', 'M_s10to20 Hz','M_s20to30 Hz' 'M_s30to40 Hz', 'M_s40toINF Hz'});
    Median_upstroke = array2table(median_binned_upstrokes, 'VariableNames', {'Md_FirstAP', 'Md_s0to10 Hz', 'Md_s10to20 Hz','Md_s20to30 Hz' 'Md_s30to40 Hz', 'Md_s40toINF Hz'});
    Summary_per_cell = [Mean_upstroke Median_upstroke];
    Summary_per_cell.Filename = repmat({filelist{i}}, height(Summary_per_cell), 1);
    %append results to overview
        if table_initialize 
            t_all = [t_all; Summary_per_cell];
        else 
            t_all = Summary_per_cell ; 
            table_initialize = true
        end
        %Anna_edition 1-Dec - handling errors
    catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end        
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages3=[errorMessages3; filelist{i} {ME}];%append the error message to a list
        end
    end
        
        
end  

% save the files into the excel - first row - mean, second row - median
%Mean_upstroke = array2table(mean_binned_upstrokes, 'VariableNames', {'M_FirstAP', 'M_s0to10 Hz', 'M_s10to20 Hz','M_s20to30 Hz' 'M_s30to40 Hz', 'M_s40toINF Hz'});
%Median_upstroke = array2table(median_binned_upstrokes, 'VariableNames', {'Md_FirstAP', 'Md_s0to10 Hz', 'Md_s10to20 Hz','Md_s20to30 Hz' 'Md_s30to40 Hz', 'Md_s40toINF Hz'});
%Summary_per_cell = [Mean_upstroke Median_upstroke];
writetable(t_all,'APs.xlsx')

%General end error message - Anna 1-Dec
if ~isempty(errorMessages3)
    for i=1:size(errorMessages3,1)
        fprintf('\nError at session %d\n',errorMessages3{i,1});%print problematic session number to screen
        errorMessages3{i,2}.message%display error message
        errorMessages3{i,2}.stack%print function and line number at which error occurred
    end
end

    
%% 
%getting the sweeps I need
%swps = nwb.getstimset.getnwbchannel.getsweep;  
%analysis of the sweeps
%swps = swps.analysesweep;
% getting th eAPs from the sweeps which have them 
%APs = swps.getap;
%%getting the frequency of every AP
%inst_freqs = [APs.freq];
% dividing the APS into the bins to get the ones we need 
%binedges = [0,0.00001, 10, 20, 30, 40, 500];

%bins = discretize(inst_freqs, binedges);
%binlabels = {'FirstAP', 's0to10 Hz', 's10to20 Hz','s20to30 Hz' 's30to40 Hz', 's40toINF Hz'};

%upstrokes = [APs.maxdvdt];

%scatter(categorical(binlabels(bins)), upstrokes) 

%choose specific APs and get data per bin
%mean_binned_upstrokes=NaN(1, numel(binlabels));
%median_binned_upstrokes=NaN(1, numel(binlabels));
%for i=1:numel(binlabels)
 %   mean_binned_upstrokes(i) = nanmean(upstrokes(bins==i) )  ; 
  %  median_binned_upstrokes(i) = nanmedian(upstrokes(bins==i) )  ;
%end

%boxplot(upstrokes, bins)
%hold on
%plot(swps)
%scatter(categorical(binlabels(bins)), upstrokes)

%% save the files into the excel - first row - mean, second row - median
%Mean_upstroke = array2table(mean_binned_upstrokes, 'VariableNames', {'M_FirstAP', 'M_s0to10 Hz', 'M_s10to20 Hz','M_s20to30 Hz' 'M_s30to40 Hz', 'M_s40toINF Hz'});
%Median_upstroke = array2table(median_binned_upstrokes, 'VariableNames', {'Md_FirstAP', 'Md_s0to10 Hz', 'Md_s10to20 Hz','Md_s20to30 Hz' 'Md_s30to40 Hz', 'Md_s40toINF Hz'});
%Summary_per_cell = [Mean_upstroke Median_upstroke];
%writetable(Summary_per_cell,'APs.xlsx')
