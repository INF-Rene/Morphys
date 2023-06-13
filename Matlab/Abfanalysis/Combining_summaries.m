%% Combine summaries Eline
%load in all the summaries and name them
%summary_patientnummer_regionletter
%Summary_177_F for example 

basedir = '/Users/elinemertens/Data/ephys/abf/Human/2023/matlab/good' ;
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};

for i = 1:length(filelist)
    load(fullfile(basedir,filelist{i})) ;
    for j = 1:length(Summary)
        Summary(j).subject_ID = ['H0' filelist{i}(9:11)] ;
        Summary(j).patcher = filelist{i}(13) ;
    end
    temp = struct2table(Summary) ;
    
    if i < 2
        Summ = temp ;
    else
        Summ = [Summ; temp] ;
    end
end

save(fullfile(basedir, 'MasterSummary'),'Summary') ;

%% only works on windows
fileName = 'Summary_excel.xlsx' ;
xlwrite(Summ,fileName);

%%
 writetable(Summ, 'Summary_all.xlsx');


%%
save(fullfile(basedir, fileName)) 

% check dit: Summ.region = categorical(Summ.region);

%%
 
gscatter(Summ.vmbaseM,Summ.vmbaseSD,Summ.region)

%%
xlwrite('Mastersummary_excel.xls', 'Summ')

%%
save(fullfile(basedir, 'MasterSummary.xlsx'), 'Summ') ;