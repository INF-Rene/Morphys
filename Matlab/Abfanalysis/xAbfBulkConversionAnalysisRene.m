% Abffile bulk conversion script
    %   
    %   Goal:   - convert and analyse abf files
    %           - skip certain protocols that are not analyzable at the moment
    %           - only analyse channels with something meaning ful associated to it: so check labbook list to see if channel
    %             was used.
    %
    %   NOTE please ensure that a global USERNAME is defined. First section in presets is a quick attempt to ensure this. 
    %   Reason is that the dialog box that would pop up in case no globale USERNAME is defined is incompatible with the parfor 
    %   loop below.
    %   
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    %        
%% presets
% set a USERNAME
listofuserids = {'AKS','EJM', 'DBH','DRU','GTS','IKS','JDZ','JOR','JSR','MBV','NAG','RBP','RWS','SHT','THK','TKN'};
global USERNAME
if isempty(USERNAME)
    idx = listdlg('PromptString','Select a userid:','SelectionMode','single','ListString',listofuserids);
    if isempty(idx)
        disp('No user name selected, script cannot run further unless the parfor loop below is changed into a for loop.')
        return
    else
        USERNAME = listofuserids{idx};
    end
end

% set paths
dir_abfs          = 'D:\Courses\LCI 2018 week 1\SpikeProfiles';
dir_mats_converts = 'D:\Courses\LCI 2018 week 1\SpikeProfiles\Converted';
dir_mats_analysed = 'D:\Courses\LCI 2018 week 1\SpikeProfiles\Analyzed';

% load table and keep only relevant columns
todo = {'2018_10_15_0003.abf';'2018_10_15_0007.abf';'2018_10_16_0036.abf';'2018_10_16_0039.abf';'2018_10_17_0001.abf';'2018_10_17_0005.abf';'2018_10_17_0007.abf';'2018_10_17_0012.abf';'2018_10_17_0016.abf'}
%% setup settings
ss = load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
ss = ss.obj;

%% To run if files are partly run already
% Reanalyze everything? reanalyze=1 
%Convert and analyze only nonconverted files? reanalyze=0

% if reanalyze==0
%     F = dir(dir_mats_analysed);
%     F={F.name};
%     F=cellfun(@(x) strsplit(x, '.'), F, 'UniformOutput', false);
%     F=cellfun(@(x) x{1}, F, 'UniformOutput', false)';
%     ttname=tt;
%     ttname=cellfun(@(x) strsplit(x, '.'), ttname{:,1}, 'UniformOutput', false);
%     ttname=cellfun(@(x) x{1}, ttname, 'UniformOutput', false);
%     todo=find(~ismember(ttname,F));
% else
%     todo=1:numel(tt);
% end
% %% load abffile objects and analyse
% conversion_fails={};
% analysis_fails={};
% analysis_error_msg={};
% conversion_error_msg={};
for i=1:numel(todo)
    % load the Abf xfile object, analyse and save
%     myrow = fn2chnnl(i,:);
%     fn    = myrow{1};
    fnabf = todo(i);
%     fnabf=fnabf{1};
    fnmat = strsplit(fnabf, '.');
    fnmat = [fnmat{1}, '.mat'];
    try
        fprintf('\n#%05d: converting %s\n',todo(i),fnabf)
        a = Abffile(fullfile(dir_abfs,fnabf),ss);
        
        % keep only amplifier channels listed in the info file
%         maxchannelset  = 1:4;
%         channels2ditch = maxchannelset(~ismember(maxchannelset, myrow{2}));
%         for ii=1:numel(channels2ditch)
%             a = a.removechannel('number',channels2ditch(ii));
%         end        
        
        % save it
        a.saveme(dir_mats_converts,fnmat)        
        
        % attempt analysis
        if ~ismember(a.proname,protocols2skip4analysis)
            try
                fprintf('\n#%05d: analysing %s\n',todo(i),fnabf)
                % analyse and save
                a = a.analyseabf;
                a.saveme(dir_mats_analysed,fnmat);
            catch err
                fprintf('#%05d, %s $$$ ANALYSIS FAIL $$$: %s\n',todo(i),fnabf,err.message);
                analysis_fails{end+1}=fnabf;
                analysis_error_msg{end+1}=err.message;
            end
        end
    catch err
        fprintf('#%05d, %s *** CONVERSION FAIL ***: %s\n',todo(i),fnabf,err.message);
        conversion_fails{end+1}=fnabf;
        conversion_error_msg{end+1}=err.message;
    end
end