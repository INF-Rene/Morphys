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
dir_abfs          = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Natalia\Selection\abfs';
dir_mats_converts = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Natalia\Selection\onrap30\converted';
dir_mats_analysed = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Natalia\Selection\onrap30\analyzed';

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
% files = dir(dir_abfs) ;
% files = struct2table(files) ;
% files = files(files.isdir==0,:) ;
% files = table2struct(files) ;

files = table2struct(readtable(fullfile(dir_abfs,'NAGselectOverview.csv'))) ;
 
%% setup settings
ss = load('C:\Users\DBHeyer\Documents\PhD\Human Database\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
ss = ss.obj;

% for abfs from AKS
ss2 = Setupsettings ;
ss2 = ss2.addchannel('number',1,'dacnum',0,'primary',2) ;

%% load abffile objects and analyse
parfor i=1:size(files,1)
    
    % load the Abf xfile object, analyse and save
    myrow = files(i);
    fn    = myrow.name(1:end-4);
    fnabf = [fn '.abf'];
    fnmat = [fn '.mat'];
  
    try
        fprintf('\n#%05d: converting %s\n',i,fnabf)
        if myrow.aks == 0
            a = Abffile(fullfile(dir_abfs,fnabf),ss);
        elseif myrow.aks == 1
            a = Abffile(fullfile(dir_abfs,fnabf),ss2);
        end
        % keep only amplifier channels listed in the info file
        maxchannelset  = 1:4;
        channels2ditch = maxchannelset(~ismember(maxchannelset, myrow.channel));
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