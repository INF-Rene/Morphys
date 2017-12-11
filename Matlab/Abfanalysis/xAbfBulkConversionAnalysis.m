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
dir_abfs          = 'C:\Users\DBHeyer\Documents\PhD\Human Database\YairData\Edata';
dir_mats_converts = 'C:\Users\DBHeyer\Documents\PhD\Human Database\YairData\Eresults\onrap30\Converted';
dir_mats_analysed = 'C:\Users\DBHeyer\Documents\PhD\Human Database\YairData\Eresults\onrap30\Analyzed';

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
files = dir(dir_abfs) ;
files = struct2table(files) ;
files = files(files.isdir==0,:) ;
files = table2struct(files) ;

%% setup settings
ss = Setupsettings ;
ss = ss.addchannel('number',1,'dacnum',1,'primary',1) ;
ss.setupsettingsname = 'custom' ;
%% load abffile objects and analyse
parfor i=1:size(files,1)
    
    % load the Abf xfile object, analyse and save
    myrow = files(i);
    fn    = myrow.name;

    try
        fprintf('\n#%05d: converting %s\n',i,fn)
        a = Abffile(fullfile(dir_abfs,fn),ss);
        
        % keep only amplifier channels listed in the info file
%         maxchannelset  = 1:4;
%         channels2ditch = maxchannelset(~ismember(maxchannelset, myrow{2}));
%         for ii=1:numel(channels2ditch)
%             a = a.removechannel('number',channels2ditch(ii));
%         end        
%         
        % save it
        a.saveme(dir_mats_converts,fn)        
        
        % attempt analysis
        if ~ismember(a.proname,protocols2skip4analysis)
            try
                fprintf('\n#%05d: analysing %s\n',i,fn)
                % analyse and save
                a = a.analyseabf;
                a.saveme(dir_mats_analysed,fn);
            catch err
                fprintf('#%05d, %s $$$ ANALYSIS FAIL $$$: %s\n',i,fn,err.message);
            end
        end
    catch err
        fprintf('#%05d, %s *** CONVERSION FAIL ***: %s\n',i,fn,err.message);
    end
end