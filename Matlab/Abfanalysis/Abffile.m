classdef Abffile < Sharedpaths & Setupsettings
    
    %   ABFFILE: A Matlab class for the axon binary file (ABF)
    %   The axon binary file is the standard data format storing electrophysiological recordings obtained using pClamp, the 
    %   Axon Instruments data acquisition software. 
    %   
    %   Notes:
    %   The class constructor for this object uses the abfload_pro.m function to read data from an abf file. The abfload_pro 
    %   is based on abf2load.m, available on matlab file exchange. It was adapted by Tim Kroon and me to add fields called
    %   'EpochSec' and 'stringSection' to the 'h' header output, if available. 
    %   -   h.EpochSec provides information on the pClamp analog waveform that was injected via each OUT (DAC) channel.
    %   -   h.stringSection provides among other things the names of pClamp protocol files and external stimulus files used   
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################

    properties (Access = public)
        guid_batch          % Globally unique identifier of parent Abfbatch object.
        fileconversiondate  % date of conversion of ABFfile
        filesystemdate      % system date of abf file storage
        filename            % original file name of ABFfile
        filedirectory       % original storage directory of file
        filetype            % axon binary file
        fileversion         % abf file version
        filesize            % size in bytes of original ABFfile
        fileduration        % total duration of recording
        filetimestart       % start date and time of recording (ms precision)
        filetimeend         % end date and time of recording (ms precision)
        savename            % name to use when saving this file
        
        dataptscheck        % check number of data points
        datasize            % size of data matrix
        datadimorderstr     % string indicating what data is stored in what dimensions in data matrix
        datadimcount        % number of dimensions
        
        proname             % name of the pClamp protocol file that was run to obtain this abf file.
        prodirectory        % directory where pClamp protocol file was stored 
        operationmode       % integer indicating mode of acquisition (e.g. 3 = Gap free, 5 = Episodic stimulation mode)
        operationmodestr    % string indicating mode of acquisition
        sampleint           % interval between samples in microseconds (1e-6)
        samplefreq          % sampling frequency in Hz
        stopwatchtime       % time of pClamp stop watch at onset of recording.
        
        nrofsweeps  = 0     % number of sweeps actually recorded (can be less than specified in protocol file)
        updatephase = 1     % phase of abf conversion (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)
          
    end
    
%####################################################### METHODS ############################################################
    
    methods
        
        % ----------------------------------------- CLASS CONSTRUCTOR--------------------------------------------------------
        function obj = Abffile(path2file, setup) 
            % This function instantiates an object of class "Abffile".
            % If a valid path is provided as input to this function, all info relevant for analysis of the axon binary file
            % (ABF) will be extracted and stored as properties of this object. Meta data extends to protocol information such 
            % as epoch definitions (as defined in analog waveform tab in pClamp), or to the atf used as stimulus files. 
            % INPUT:
            % "path2file" = a valid path to an Abffile location, with or without ".abf" extension, cannot be empty or have a
            %                 different extension.
            % "setup"     = a Setupsettings object containing info required to match DACs to ADC channels. This object serves
            %               as a template to construct the setupsettings superclass of Abffile. 
                       
%             % check setup settings object. No empty dacs allowed, no duplicates
%             if numel(setup.getchannel.get('dacnum')) < setup.nrofchannels, 
%                 error('Setup settings error: Not all channels have a DAC number specified.'); 
%             %elseif unique(cell2mat(setup.getchannel.get('dacnum'))) < setup.nrofchannels, %bugfix by RWS 21-09-2017        
%             elseif length(unique(cell2mat(setup.getchannel.get('dacnum')))) < setup.nrofchannels,
%                 error('Setup settings error: Mulitple channels seem to share the same DAC number.'); 
%             end
            
            % if setup checks passed, proceed to superclass constructor methods. Note that Abffile will also inherit 
            % properties from the Sharedmethods class via its Setupsettings superclass.
            obj = obj@Sharedpaths;
            obj = obj@Setupsettings(setup);
            
            % return empty if no input arguments
            if nargin == 0, return; end  
            
            % check path input
            if nargin == 0 || isempty(path2file)
                error('No pathname provided.')
            elseif ~ischar(path2file)
                error('Only string input accepted.')
            else
                [fileDir,~,ext] = fileparts(path2file);
                if isempty(ext), 
                    path2file = [path2file '.abf']; % add ".abf" if required
                elseif ~strcmp(ext,'.abf'),         % if extension doesn't match, error
                    error('Only ''.abf'' files supported, provided ''%s''',ext)
                end
                if ~exist(path2file,'file'),        % if path doesn't exist, error
                    error('Cannot find specified file: \n%s',path2file), 
                end    
            end
            
            % load ABF file
            obj.fileconversiondate = char(datetime(datestr(now()),'Format',obj.datetimefmt));     

            % load ABFfile (more loading options from abload_pro could be added here...)
            [dataMtx,si,h] = abfload_pro(path2file);
            % assignin('base','h',h)

            % Add file info
            fileinfo         = dir(path2file);
            obj.filename     = fileinfo.name;
            obj.filedirectory= fileDir;
            obj.filetype     = 'Axon Binary File';
            obj.fileversion  = h.fFileSignature;   
            obj.filesize     = fileinfo.bytes;
            
            % create a savename
            obj.savename = sprintf('Abf_%s_%s.mat',obj.filename(1:end-4),obj.guid);
            
            % Get dates and durations
            obj.filesystemdate= fileinfo.date;
            obj.fileduration  = duration(0,0,0,diff(h.recTime*1e3),'Format',obj.durationfmt);
            obj.filetimestart = datetime(datevec(num2str(h.uFileStartDate),'yyyymmdd') + datevec(duration(0,0,0,h.uFileStartTimeMS)),'Format',obj.datetimefmt);
            obj.filetimeend   = obj.filetimestart+obj.fileduration;

            % Check nr. of datapoints (moved here from abfload_pro). Still not quite sure why this check is nesessary and why/how it sometimes
            % fails... ABF files with failed check (so where obj.dataptscheck=0) seem perfectly fine and can be analysed without apparent issues.
            if rem(h.dataPts,h.nADCNumChannels)>0 || rem(h.dataPtsPerChan,h.lActualEpisodes)>0, obj.dataptscheck = 0;
            else obj.dataptscheck = 1;
            end  

            % General protocol info
            path_and_file    = obj.extractstringsectioninfo(h.stringSection, 1);
            obj.proname      = path_and_file{2};
            
            % HACK ALERT!!! This is just for 2017-08-16 session, do not analyse bridge balance abfs!!!
            if ismember(obj.proname,{'eCode_1_BridgeBalance','BridgeBalance','BB','BB2'})
                error('Bridge balance error!')
            end
            
            obj.prodirectory     = path_and_file{1};
            obj.operationmode    = h.nOperationMode;
            obj.operationmodestr = h.nOperationModeStr;
            if ~any(obj.operationmode == obj.implrecmode), 
                error('ABFfile operation mode unknown. Only Episodic and Gapfree recordings please.')
            end
            
            obj.sampleint     = si; % sample interval is in microseconds (= abfload default)
            obj.samplefreq    = 1/(obj.sampleint*1e-6); 
            obj.stopwatchtime = duration(0,0,h.uStopwatchTime);

            % Data dimensions
            % Data dimensions of the data matrix outputted by abfload_pro.m can vary. If recording only one channel or sweep, 
            % the number of dimensions of the data matrix is 2, whereas recording more channels and sweeps results in a data
            % matrix of 3 dimensions. 
            
            % Data dimensions
            obj.datasize     = size (dataMtx);
            obj.datadimcount = ndims(dataMtx);
            
            % Find out number of input channels and sweeps:
            nrofanalogins_init = numel(h.ADCSec);           % initial number of analogins
            switch obj.operationmode
                case 3, obj.nrofsweeps = 1;                 % initial number of sweeps = 1 in case of Gapfree
                case 5, obj.nrofsweeps = h.lActualEpisodes; % initial number of sweeps
            end            

            % Analog inputs
            % Basically go through the Analog inputs, check what amplifier Channel it belongs to based on setup settings,
            % then add an Analogin to the relevant Channel.
            fprintf('Adding analog inputs...\n')
            for i=1:nrofanalogins_init
                infostruct = struct;
                % get the Analogin number (so that is the number of the Analog input of th digidata)
                infostruct.number = h.ADCSec(i).nADCNum;
                for ii=1:obj.nrofchannels
                    if ismember(infostruct.number, obj.getchannel(ii).getanaloginputnrs)
                        
                        % make a struct with necessary info to make an Analogin object
                        infostruct.signal                = obj.getchannel(ii).getsignalfromanaloginnr(infostruct.number);
                        infostruct.adcusername           = h.recChNames{i};
                        infostruct.units                 = h.recChUnits{i};
                        infostruct.samplefreq            = 1/(h.si*1e-6); 
                        infostruct.telegraphenabled      = h.ADCSec(i).nTelegraphEnable;
                        infostruct.telegraphinstrument   = h.ADCSec(i).nTelegraphInstrument;
                        infostruct.gain                  = h.ADCSec(i).fTelegraphAdditGain;
                        infostruct.lowpassfilter         = h.ADCSec(i).fTelegraphFilter;
                        infostruct.instrumentscalefactor = h.ADCSec(i).fInstrumentScaleFactor;
                        infostruct.h = h; % only reuired for making sweeps now, for adddates. Remove soon;
                        
                        % Select channel data and make 2D array (samples x sweeps)
                        if     obj.nrofsweeps     == 1, analogindata = dataMtx(:,i); obj.datadimorderstr = 'samples x analogins';
                        elseif nrofanalogins_init == 1, analogindata = dataMtx(:,:); obj.datadimorderstr = 'samples x sweeps';
                        else   analogindata = squeeze(dataMtx(:,i,:)); % squeeze to 2D
                               obj.datadimorderstr = 'samples x analogins x sweeps';
                        end
                        
                        % add Analogin object to channel in list
                        obj.channels(ii) = obj.getchannel(ii).addin(infostruct,analogindata);
                    end
                end
            end 
            % since data has now been passed on, we can delete the data matrix. Use obj.getdata to retrieve data from
            % downstream channels/analogins/sweeps.
            clear dataMtx          
                        
            % Digital-to-Analog (DAC) Channels
            %   Abffiles can have multiple of these channels, which can either be defined in "h.EpochSec" or, in case a 
            %   stimulus file (".atf") was used for generating the anaolog waveform, the details must be extracted from 
            %   "h.stringSection". 
            %   NOTE: Overlap can occur, the stringsection lists all the atf files associated with an analog waveform tab in 
            %   pClamp, regardless of whether the analog out was enabled or not, and regardless of whether the stimfile was 
            %   actually used. If EpochSec AND stringsection both define an analog OUT (so share same dacnum), then EpochSec 
            %   takes precedence, since epoch section only includes enabled channels where the waveform epoch table was used
            %   as output. (Update RWS 21-11-2018): Tested this and turns out not to be true: if a stimulusfile is used on a
            %   file that previously had an EpochSec, than the EpochSec information stays in the file. Therefore changed
            %   precedence to stimfile.
            %   NOTE: To prevent DAC channels that do have a stimulus file but where DAC was disabled still initialise, 
            %   only DAC channels defined in the provided Setupsettings object will be taken into account.
            fprintf('Adding analog outputs...\n')
%             if contains(obj.proname, 'endurance master9')
%                 h.stringSection = [h.stringSection ' C:\Rene\Protocols\Endurance\endurance master9.atf '];
%             end
            tmpoutinfo = obj.extractstringsectioninfo(h.stringSection,2); 
            if isfield(h,'EpochSec') 
                dacnumberset       = unique([h.EpochSec.nDACNum]);
                
                for i=dacnumberset
                    if ~ismember(i, tmpoutinfo.number) %precedence of stimfile over epochsec
                        infostruct  = struct;                                                       % start empty
                        epochstruct = h.EpochSec([h.EpochSec.nDACNum] == i);                        % make table of epoch information for this OUT channel
                        infostruct.number       = epochstruct(1).nDACNum;                           % collect OUT nr
                        infostruct.name         = '';
                        infostruct.units        = '';
                        infostruct.path2pro     = obj.prodirectory;
                        infostruct.profilename  = obj.proname;
                        infostruct.path2stim    = '';
                        infostruct.stimfilename = '';
                        
                        epochstruct = rmfield(epochstruct,'nDACNum');                               % remove DAC number from table (now redundant)
                        epochstruct = obj.standardiseepochstruct(epochstruct,obj.sampleint*1e-3);   % standardise epochStruct
                        t = struct2table(epochstruct);
                        
                        % now an ugly fix for cases where epochstruct has only 1 epoch. In these cases, struct2table converts
                        % cells into char, which leads to errors later. Next measures are to keep idxstr,timespan and typestr the
                        % same type throughout.
                        if isscalar(epochstruct)
                            if ischar(t.idxstr),   t.idxstr   = {t.idxstr}; end
                            if ischar(t.timespan), t.timespan = {t.timespan}; end
                            if ischar(t.typestr),  t.typestr  = {t.typestr}; end
                        end
                        
                        infostruct.analogwaveformtable = Analogwaveformtab(t);                      % convert to analog waveform tab object
                        infostruct.analogwaveformtable.samplefreq = obj.samplefreq;                 % assign sampleFrequency of this file to analogWaveformTab.
                        infostruct.analogwaveformtable.nrofsweeps = obj.nrofsweeps;                 % assign number of sweeps in this file to analogWaveformTab.
                        
                        % add to Analogout object to relevant channel.
                        for ii=1:obj.nrofchannels
                            if infostruct.number == obj.getchannel(ii).get('dacnum')
                                obj.channels(ii) = obj.getchannel(ii).addout(infostruct);
                            end
                        end
                    end
                end
            end

            % Analogouts defined in string section to employ external stimulus files (atf files, listed in h.stringSection)
            
            if ~isempty(tmpoutinfo)                
                outchannels = tmpoutinfo{:,'number'}';
                
                % Ignore DACs in stringSection that are already defined by protocol section. These MUST be inactive atf
                % files and therefore not interesting.
                if numel(tmpoutinfo.stimfilename{1}) < 5 | ~strcmp('eCode',tmpoutinfo.stimfilename{1}(1:5));
                    outchannels = outchannels(~ismember(outchannels,obj.getanalogoutputnrs_current));
                end
                
                % Ignore DACs in stringSection that are not assigned in setup
                outchannels = outchannels(ismember(outchannels,obj.getanalogoutputnrs));

                for i=outchannels
                    infostruct = table2struct(tmpoutinfo(tmpoutinfo{:,'number'}==i,:));
                    infostruct.path2pro    = obj.prodirectory;
                    infostruct.profilename = obj.proname;
                    
                    % check if the atf file has an existing Analogwaveformtab object associated with it
                    if exist(fullfile(obj.path2analogwaveforms,[infostruct.stimfilename '.mat']),'file')==2
                        infostruct.analogwaveformtable = load(fullfile(obj.path2analogwaveforms,infostruct.stimfilename)); 
                        infostruct.analogwaveformtable = infostruct.analogwaveformtable.obj;
                        if ~isempty(infostruct.analogwaveformtable)
                            % assign sampleFrequency and number of sweeps in this file to analogWaveformTab.
                            infostruct.analogwaveformtable.samplefreq = obj.samplefreq; 
                            infostruct.analogwaveformtable.nrofsweeps = obj.nrofsweeps;    
                        end
                    else % if not, skip.
                        infostruct.analogwaveformtable = [];
                        warning('%s, OUT #%d: Stimulus file "%s" not found among Analog Waveform collection.',obj.filename,i,infostruct.stimfilename)
                    end
                    
                    % add to Analogout object to relevant channel.
                    for ii=1:obj.nrofchannels
                        if infostruct.number == obj.getchannel(ii).get('dacnum')
                            if ~isempty(obj.getchannel(ii).getout)
                                obj.channels(ii) = obj.getchannel(ii).removeout ;
                            end
                            obj.channels(ii) = obj.getchannel(ii).addout(infostruct);
                            if ~isempty(obj.getchannel(ii).getin('signal','secondary'))
                                sf = obj.getchannel(ii).getout.getscalefactor(obj.getchannel(ii).getin('signal','secondary'));
                                obj.channels(ii).analogouts = obj.getchannel(ii).getout.set('scalefactor',sf);
                                if ~isempty(obj.getchannel(ii).getout.analogwaveformtable)
                                obj.channels(ii).analogouts.analogwaveformtable = obj.getchannel(ii).getout.analogwaveformtable.set('scalefactor',sf);
                                end
                            end
                        end
                    end
                end
            end
            
            % Ditch Channels with no Analoginputs and update guid and filename fields
            %for i=1:obj.nrofchannels  %Bugfix RWS 21-09-2017
            for i=obj.nrofchannels:-1:1  
                if isempty(obj.getchannel(i).nrofanalogins) || obj.getchannel(i).nrofanalogins == 0
                    obj = obj.removechannel(i); 
                else
                    obj.channels(i) = obj.getchannel(i).set('guid_abf',obj.guid,'filename',obj.filename);
                end
            end

            % Make Epochs using Analogwaveformtab of Analogout
            fprintf('Adding epochs...\n')
            obj = obj.makeepochs;
            
            % Get holding current or voltage
            for i=1:obj.nrofchannels
                if ~isempty(obj.getchannel(i).getin('signal','secondary')) && ~isempty(obj.getchannel(i).getout) && ~isempty(obj.getchannel(i).getout.analogwaveformtable) 
                    [holdI, holdV] = obj.getchannel(i).getout.getholdingIorV(obj.getchannel(i));
                    obj.channels(i).analogouts = obj.getchannel(i).getout.set('holdingI',holdI);
                    obj.channels(i).analogouts = obj.getchannel(i).getout.set('holdingV',holdV);
                end      
            end     
        end
        
        % ------------------------------------------ HELPER METHODS ---------------------------------------------------------
        function identify(obj)
            % This dissplays a small summary of the contents of the ABF file object. 
            if ~isscalar(obj), error('Object must be scalar.'); end
            disp(obj.idmessage)
        end
        
        function id = idmessage(obj)
            % This function returns a small text summary of the contents of the ABF file object. 
            if ~isscalar(obj), error('Object must be scalar.'); end
            bfilts = [obj.getchannel.getin.lowpassfilter];
            if range(bfilts)~=0 % so if different bessel filters were used for different channels...
                bfiltstr = sprintf('%d to %d kHz',min(bfilts)/1000,max(bfilts)/1000);
            else
                bfiltstr = sprintf('%d kHz',bfilts(1)/1000);
            end
            id = sprintf([   'Rec date:     %s\n',...
                             'Protocol:     %s\n',...
                             'AmpChannels:  %d\n',...
                             'AnalogInputs: %d\n'...
                             'Sweeps:       %d\n',...
                             'Sample freq:  %d kHz\n',...
                             'Bessel filt:  %s\n'],...
                             obj.filesystemdate, obj.proname, obj.nrofchannels, sum([obj.getchannel.nrofanalogins]),...
                             obj.nrofsweeps,floor(obj.samplefreq/1000),bfiltstr);
        end
        
        function id = idmessagetab(obj)
            % This function returns a small text summary of the contents of the ABF file object. Tabs used for spacing
            if ~isscalar(obj), error('Object must be scalar.'); end
            bfilts = [obj.getchannel.getin.lowpassfilter];
            if range(bfilts)~=0 % so if different bessel filters were used for different channels...
                bfiltstr = sprintf('%d to %d kHz',min(bfilts)/1000,max(bfilts)/1000);
            else
                bfiltstr = sprintf('%d kHz',bfilts(1)/1000);
            end
            id = sprintf([   'Rec date:\t\t%s\n',...
                             'Protocol:\t\t%s\n',...
                             'AmpChannels:\t%d\n',...
                             'AnalogInputs:\t%d\n'...
                             'Sweeps:\t\t\t%d\n',...
                             'Sample freq:\t%d kHz\n',...
                             'Bessel filt:\t%s\n'],...
                             obj.filesystemdate, obj.proname, obj.nrofchannels, sum([obj.getchannel.nrofanalogins]),...
                             obj.nrofsweeps,floor(obj.samplefreq/1000),bfiltstr);
        end
               
        function data = getdata(obj)
            % retrieve all data (now stored downstream in sweeps). Data returned as N-D matrix. 
            % ABFFILE arrays will result in cellarray of N-D matrixes.
            data  = cell(size(obj));
            for i = 1:numel(obj); data{i} = [obj(i).getchannel.getin.getdata]; end
            if numel(data)==1, data=data{:}; end
        end
        
        function obj = selectsweep(obj,varargin)
            % make an ABFFILE object including only specified SWEEP(s).
            % Note, using the 'end' keyword as in 1:4:end does not work for this function. 
            % --------------------
            % See also SELECTITEM
            if nargin == 1, return; end
            for i = 1:numel(obj)
                for ii = 1:obj(i).nrofchannels
                    for iii=1:obj(i).getchannel(ii).nrofanalogins
                        % replace analogin with one filtered based on varargin
                        obj(i).channels(ii).analogins(iii) = obj(i).getchannel(ii).getin(iii).selectsweep(varargin{:});
                    end
                end
                obj(i).nrofsweeps=numel(obj(i).getchannel(ii).getin(iii).getsweep);
            end
        end
        
        function obj = selectin(obj,varargin)
            % make an ABFFILE object including only specified SWEEP(s).
            % Note, using the 'end' keyword as in 1:4:end does not work for this function. 
            % --------------------
            % See also SELECTITEM
            if nargin == 1, return; end
            for i = 1:numel(obj)
                for ii = 1:obj(i).nrofchannels
                        % replace channel with one filtered based on varargin
                        obj(i).channels(ii) = obj(i).getchannel(ii).selectin(varargin{:});
                end
            end
        end
        
        function obj = selectepoch(obj,varargin)
            % make an ABFFILE object including only specified epoch(s).
            % Note, using the 'end' keyword as in 1:4:end does not work for this function. 
            % --------------------
            % See also SELECTITEM
            if nargin == 1, return; end
            for i = 1:numel(obj)
                for ii = 1:obj(i).nrofchannels
                    for iii=1:obj(i).getchannel(ii).nrofanalogins
                        for iiii=1:obj(i).nrofsweeps
                            % replace sweep with one filtered based on varargin
                            obj(i).channels(ii).analogins(iii).sweeps(iiii) = obj(i).getchannel(ii).getin(iii).getsweep(iiii).selectepoch(varargin{:});
                        end
                    end
                end
            end
        end
        
        function aplist = apsweeps(obj)
            % Returns a list with the amount of AP's detected in each sweep of the
            % primary signal of each channel
            % 
            % --------------------
            % See also SELECTITEM
            if ~isscalar(obj), error('Object must be scalar.'); end
            
            % method dependent on whether abf is already analyzed
            if obj.updatephase < 3
                obj=obj.selectin('signal', 'primary');
                aplist=zeros(obj.nrofchannels, obj.nrofsweeps);
                for i = 1:obj.nrofchannels
                    for ii=1:obj.nrofsweeps
                        %find nr of APs in every sweep
                        swp=obj.getchannel(i).getin(1).getsweep(ii).findaps;
                        swp=swp.updateapstats;
                        aplist(i, ii) = swp.nrofaps;
                    end
                end
            else
                obj=obj.selectin('signal', 'primary');
                aplist=zeros(obj.nrofchannels, obj.nrofsweeps);
                for i = 1:obj.nrofchannels
                    for ii=1:obj.nrofsweeps
                        aplist(i, ii) = obj.getchannel(i).getin(1).getsweep(ii).nrofaps;
                    end
                end
            end
        end
        function firstapsweep = firstapsweep(obj)
            % Return the first sweep in which an AP is found in any
            % channel
            % 
            % --------------------
            % See also SELECTITEM
            if ~isscalar(obj), error('Object must be scalar.'); end
            firstapsweep=[];
            % method dependent on whether abf is already analyzed
            if obj.updatephase < 3
                obj=obj.selectin('signal', 'primary');
                for i = 1:obj.nrofchannels
                    for ii=1:obj.getchannel(i).nrofanalogins
                        for iii=1:obj.nrofsweeps
                            %find nr of APs in every sweep
                            swp=obj.getchannel(i).getin(ii).getsweep(iii).findaps;
                            swp=swp.updateapstats;
                            if swp.nrofaps > 0, firstapsweep=iii; return, end
                        end
                    end
                end
            else
                obj=obj.selectin('signal', 'primary');
                for i = 1:obj.nrofchannels
                    for ii=1:obj.getchannel(i).nrofanalogins
                        for iii=1:obj.nrofsweeps
                            if obj.getchannel(i).getin(ii).getsweep(iii).nrofaps > 0, firstapsweep=iii; return, end
                        end
                    end
                end
            end
        end
 
        % ----------------------------------------- UPDATING METHODS --------------------------------------------------------        
        function obj = analyseabf(obj)
            % analyse primary inputs of Channels in list of Abffile object
            % checks for action potentials and calculates passive properties for steps.
            for i = 1:numel(obj), 
                switch obj(i).updatephase
                    case 1
                        error('Please add settings to the Channels of Abffile object "%s" first.',obj(i).filename);
                    case 2
                        fprintf('analysing abf %s\n',obj(i).filename)
                        obj(i).channels    = obj(i).getchannel.analysechannel; 
                        obj(i).updatephase = 3;
                    case 3
                        error('Abffile object "%s" has already been analysed.',obj(i).filename);
                end
            end
        end
  
        function obj = makeepochs(obj)
            % Add Epochs to Sweeps. For every Channel, uses the Analogwaveformtab object of every Analogout
            % to add Epochs to all Sweeps of Analogin objects belonging to same Channel (which thus share the same Epochs).
            if ~isscalar(obj), error('Abffile object must be scalar.'); end
            switch obj.updatephase
                case 1
                    for i = 1:obj.nrofchannels                                      
                        obj.channels(i) = obj.getchannel(i).makeepochs;
                    end
                    % update class phase to indicate that setup configuration info has been added succesfully.
                    obj.updatephase = 2;
                case 2 
                    disp('Epoch information already added! No action taken.')
                case 3
                    disp('Epoch information already added and file has been analysed. No action taken.')
            end
            obj = obj.updatechannelstats;
        end
                
        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function [hfig,ax_handles] = plot(obj,varargin)
            % Plots the Abffile object. Returns figure handle.
            % --------------------
            if isscalar(obj)
                nrofsubplots = numel(obj.getchannel.getin);
                ax_handles   = zeros(nrofsubplots,1);
                subplotcount = 1;
                for i=1:obj.nrofchannels
                    
                    for ii=1:obj.getchannel(i).nrofanalogins
                        ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                        subplotcount  = subplotcount+1;
                        
                        analogin2plot = obj.getchannel(i).getin(ii);
                        analogin2plot.plot(varargin{:})
                        if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
                            configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.getchannel(i).number);
                        else
                            configstring = '';
                        end
                        
                        title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
                        ylabel(analogin2plot.units)
                    end
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
                hfig = gcf;
            else
                arrayfun(@(x) obj(x).plot(varargin{:}),1:numel(obj),'UniformOutput',false)
            end
        end       
               
        function hfig = plotanalysis(obj)
            % Plots the Abffile object with Epoch boundaries and results of basic analysis. Returns figure handle.
            % --------------------
            if isscalar(obj)
                nrofsubplots = numel(obj.getchannel.getin);
                ax_handles   = zeros(nrofsubplots,1);
                subplotcount = 1;
                for i=1:obj.nrofchannels
                    for ii=1:obj.getchannel(i).nrofanalogins
                        ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                        subplotcount  = subplotcount+1;
                        
                        analogin2plot = obj.getchannel(i).getin(ii);
                        analogin2plot.plotanalysis;
                        if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
                            configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.getchannel(i).number);
                        else
                            configstring = '';
                        end
                        
                        title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
                        ylabel(analogin2plot.units)
                    end
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
                hfig = gcf;
            else
                arrayfun(@(x) obj(x).plotanalysis,1:numel(obj),'UniformOutput',false)
            end
        end             
    end
    
    methods (Static)
        
        function [infoout] = extractstringsectioninfo(stringsection, request)
            % ABF stringSection info extractor.
            %
            % A lot of information is locked in the "stringSection" of the abf file, including the protocol name, 
            % inputs and outputs, units, stimulus files used... etc. See example stringSection: 
            % Clampex C:\Data\Thijs\Protocols\eCode_PRO_Scaled\eCode_2_NoisePP.pro IN 1 pA IN 2 mV IN 3 mV IN 5 pA Cmd 0 mV Cmd 1 mV C:\Data\Thijs\Protocols\eCode_ATF_Scaled_M1\eCode_2_NoisePP.atf Cmd 2 pA C:\Data\Thijs\Protocols\eCode_ATF_Scaled_M2\eCode_2_NoisePP.atf Cmd 3 mV AO #4 mV AO #5 mV AO #6 mV AO #7 mV 
            %
            % This function allows extraction of info from this big string.
            % stringSection = the stringSection as stored in "h.stringSection" in the "h" struct returned by "abfload_pro.m".
            % request == 1: return protocol file (.pro) directory and filename as 1x2 cell
            % request == 2: return stimulus file (.atf) directory, filename and associated DAC name and number (name/number of the Analog Output channel) as table
            % request == 3: return full input-output table as table
            %
            % Assumption: only .pro and .atf files are listed in string section

            %########################################################################################################################
            % replace all \ with / to be compatable with macbook files 
            stringsection = strrep(stringsection, '\', '/') ;


            % Find all paths first:
            stringsection = stringsection(stringsection~='#'); % remove annoying hashtags
            stringsection(uint8(stringsection)==0)=' ';        % fix weird undetectable spaces
            pathsIdx = regexp(stringsection,':|.pro |.atf ');  % find path strings (assuming here that only .pro and .atf are possible!)
            for i=1:numel(pathsIdx)/2
                pathsidxfull(i,1:2) = [pathsIdx(i*2-1)-1, pathsIdx(i*2-1+1)+3];
                pathsidxstr {i,1}   = stringsection(pathsIdx(i*2-1)-1: pathsIdx(i*2-1+1)+3);
            end    

           
            switch request
                case 1 % Return protocol (assuming there can be only one for an abf file):
                    [filepath,filename] = fileparts(pathsidxstr{[cellfun(@(x) strcmp('.pro',x(end-3:end)),pathsidxstr)],:});
                    infoout = {filepath,filename};    

                case {2,3} % 
                    default_section_nrs = 3; % with no stimfiles defined, sections in stringSection are composed of 3 parts: 'IN 1 pA' for example
                    iotable = cell2table(cell(0,6),'VariableNames',{'IOtype','number','units','path2stim','stimfilename','name'});
                    for i=1:size(pathsidxfull,1)
                        if i<size(pathsidxfull,1)
                            io = strsplit(stringsection(pathsidxfull(i,2)+1:pathsidxfull(i+1,1)-1),' ');
                        else
                            io = strsplit(stringsection(pathsidxfull(i,2)+1:end),' ');
                        end
                        io = io(cellfun(@(x) ~isempty(x),io));                                % remove empties
                        if ismember('RawOutput', io)
                            n=find(strcmp(io, 'RawOutput'));
                            io(n+1:end+1)=io(n:end);
                            io{n+1}='2';
                        end
                        io = reshape(io,default_section_nrs,size(io,2)/default_section_nrs)'; % reshape into discrete sections
                        iotable = [iotable; cell2table(cat(2,io,repmat({''},size(io,1),2),...
                            arrayfun(@(x) strjoin(io(x,1:2)),1:size(io,1),'UniformOutput',false)'),... % add IO names (e.g. Cmd 1)
                            'VariableNames',{'IOtype','number','units','path2stim','stimfilename','name'})];
                        iotable.name=strrep(iotable.name, 'RawOutput 4', 'RawOutput');

                        % split StimFile into path and file and add to IOtable
                        if i+1<=numel(pathsidxstr)
                            [p,fn] = fileparts(char(pathsidxstr(i+1,1)));
                            iotable{end,{'path2stim','stimfilename'}} = {p,fn}; %pathsIdx_str(i+1,1);
                        end       
                    end

                    % fix numbers...
                    number_tmp = table(arrayfun(@(x) str2num(iotable{x,'number'}{:}),1:size(iotable,1))','VariableNames',{'number'});
                    iotable(:,'number')=[];        % delete old
                    iotable=[iotable,number_tmp]; % append new (now as number instead of string)

                    if request == 2 % return only identified stimfiles, with additional DACinfo as table
                        infoout = iotable(cellfun(@(x) ~isempty(x),iotable{:,'stimfilename'}),{'name','number','units','path2stim','stimfilename'});
                    else % return full table
                        infoout=iotable; 
                    end
            end
        end   

        function epochstruct = standardiseepochstruct(epochstruct,idx2ms)
            % This function standardises the struct with epoch information, to ensure that the epoch struct from h.EpochSec 
            % (case where epoch info comes from pClamp waveform tab) and the epoch struct from ATF files are similar and can 
            % be used interchangebly.
            % --------------------
            
            % for ease of use first convert to table
            epochtable = struct2table(epochstruct);
            
            % change variable names
            variablenameconversion = {
                'nEpochNum'         , 'number'
                'nEpochType'        , 'type'
                'fEpochInitLevel'   , 'firstlevel'
                'fEpochLevelInc'    , 'deltalevel'
                'lEpochInitDuration', 'timespan'
                'lEpochDurationInc' , 'deltaduration'
                'lEpochPulsePeriod' , 'pulseperiod'
                'lEpochPulseWidth'  , 'pulsewidth'
                };
            vars2convert = ismember(epochtable.Properties.VariableNames,variablenameconversion(:,1));
            epochtable.Properties.VariableNames(vars2convert) = variablenameconversion(vars2convert,2);
            
            % add variables
            epochtypes             = {'step';'ramp';'train';'trngl';'cos';'';'biphsc';'chirp';'noisepp';'noisespiking';'truenoise'}; 
            epochtable.idx         = (1:height(epochtable))';               % add indexes
            epochtable.idxstr      = cellstr(char(epochtable{:,'idx'}+64)); % add epoch letters (as in pClamp waveform tab)
            epochtable.typestr     = epochtypes(epochtable{:,'type'});      % add strings to specify epoch type
            epochtable.maxfrequency= zeros(height(epochtable),1);           % frequencies (for chirps, not available in pClamp waveform tab)
            
            % change variable type; timespan can be string format when epoch table is based on atf info file, so
            % standardise to string everywhere
            newvals = epochtable.timespan*idx2ms;   % convert duration to milliseconds (was samples)
            newvals = strsplit(num2str(newvals'));  % make into cell array of strings
            epochtable.timespan = [];               % delete old variable column in table
            epochtable.timespan = newvals';         % replace with new.
            epochstruct = table2struct(epochtable); % convert table back to struct
            epochstruct = orderfields(epochstruct); % sort fieldnames
        end 
    end             
end

