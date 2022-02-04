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
listofuserids = {'AKS','DBH','DRU','GTS','IKS','JDZ','JOR','JSR','MBV','NAG','RBP','RWS','SHT','THK','TKN'};
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
dir_abfs          = 'D:\Morphys\Data\Electrophysiology\Abffiles\MBV\abf';
dir_mats_converts = 'D:\Morphys\Data\Electrophysiology\Abffiles\MBV\abf\Converted';
dir_mats_analysed = 'D:\Morphys\Data\Electrophysiology\Abffiles\MBV\abf\Analyzed';

% leave protocols
protocols2skip4analysis = { 'eCode_1_BridgeBalance'
                            'BridgeBalance'
                            'BB'
                            'BB2'
                            'Resonance_100pA'
                            'eCode_2_Resonance'
                            'eCode_2_Spontaneous'
                            'eCode_1_Delta'
                         };

%% check info file for meaningful channels
dir_info = 'D:\Morphys\Data\Labbooks\MBV\Preps';
fn_info  = 'Labbook_2017_10_02.xlsx';

% load table and keep only relevant columns
tt = readtable(fullfile(dir_info,fn_info),'Sheet','Channel');
tt(:,~ismember(tt.Properties.VariableNames,{'ABF_name','ampChannel'}))=[];

% make list of channels to keep for every recfile
[C,IA,IC] = unique(tt.ABF_name);
fn2chnnl  = cat(2,C,arrayfun(@(x) tt.ampChannel(IC==x),1:numel(C),'UniformOutput',false)');

%% setup settings
ss = load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_MBV.mat');
ss = ss.obj;

%% load abffile objects and analyse
conversion_results={};
analysis_results={};

reanalyze=0;
if reanalyze==0
    F = dir(dir_mats_analysed);
    F={F.name};
    F=cellfun(@(x) strsplit(x, '.'), F, 'UniformOutput', false);
    F=cellfun(@(x) x{1}, F, 'UniformOutput', false)';
    ttname=fn2chnnl;
    ttname=cellfun(@(x) strsplit(x, '.'), ttname(:,1), 'UniformOutput', false);
    ttname=cellfun(@(x) x{1}, ttname, 'UniformOutput', false);
    todo=find(~ismember(ttname,F));
else
    todo=1:numel(fn2chnnl);
end

%%
parfor i=1
    
    % load the Abf xfile object, analyse and save
    myrow = fn2chnnl(todo(i),1:2);
    fn    = myrow{1};
    fnabf = [fn '.abf'];
    fnmat = [fn '.mat'];
    try
        fprintf('\n#%05d: converting %s\n',i,fnabf)
        a = Abffile(fullfile(dir_abfs,fnabf),ss);
        
        % keep only amplifier channels listed in the info file
        maxchannelset  = 1:4;
        channels2ditch = maxchannelset(~ismember(maxchannelset, myrow{2}));
        for ii=1:numel(channels2ditch)
            a = a.removechannel('number',channels2ditch(ii));
        end        
        
        % save it
        a.saveme(dir_mats_converts,fnmat)

        
        % attempt analysis
        if ~ismember(a.proname,protocols2skip4analysis)
            try
                fprintf('\n#%05d: analysing %s\n',i,fn)
                % analyse and save
                a = a.analyseabf;
                a.saveme(dir_mats_analysed,fnmat);

            catch err
                fprintf('#%05d, %s $$$ ANALYSIS FAIL $$$: %s\n',i,fn,err.message);

            end
        end
    catch err
        fprintf('#%05d, %s *** CONVERSION FAIL ***: %s\n',i,fn,err.message);
    end
end