classdef NWBfile < Sharedpaths & Sharedmethods
    
    %   NWBfile: 
    %   This uses matlab's built in hdf5 functions to load electrophysiology data from a 
    %   NWB file (Neurodata Without Borders, see https://neurodatawithoutborders.github.io/)
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       René Wilbers (renewilbers@gmail.com)
    %   Created:      29-03-2019
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################

    properties (Access = public)
        fileconversiondate  % date of conversion of ABFfile
        filesystemdate      % system date of abf file storage
        filecreatedate
        filename            % original file name of ABFfile
        filedirectory       % original storage directory of file
        filetype            % axon binary file
        fileversion         % abf file version
        filesize            % size in bytes of original ABFfile
        fileduration        % total duration of recording
        filetimestart       % start date and time of recording (ms precision)
        filetimeend         % end date and time of recording (ms precision)
        savename            % name to use when saving this file
        
        labbooktext
        labbooktext_keys
        labbooknum
        labbooknum_keys
        labbook_timestamps
        
        activeHS            % active Headstages
        stimsetfilters
        stimsets
        
        nrofstimsets = 0
        nrofsweeps  = 0     % number of sweeps actually recorded (can be less than specified in protocol file)
        updatephase = 1     % phase of abf conversion (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)
          
    end
    
%####################################################### METHODS ############################################################
    
    methods
        
        % ----------------------------------------- CLASS CONSTRUCTOR--------------------------------------------------------
        function obj = NWBfile(fn, stimsetfilters, sweepselect) 
            % stimsetfilters is either a string or a cell array of strings to specify the search term for the stimsets that 
            % should be loaded in the NWBfile object
            
            % This function instantiates an object of class "Abffile".
            obj = obj@Sharedpaths;
            % return empty if no input arguments
            if nargin == 0, return; end  
            

            % check path input
            if nargin == 0 || isempty(fn)
                error('No pathname provided.')
            elseif ~ischar(fn)
                error('Only string input accepted.')
            else
                [fileDir,~,ext] = fileparts(fn);
                if isempty(ext), 
                    fn = [fn '.nwb']; % add ".abf" if required
                elseif ~strcmp(ext,'.nwb'),         % if extension doesn't match, error
                    error('Only ''.nwb'' files supported, provided ''%s''',ext)
                end
                if ~exist(fn,'file'),        % if path doesn't exist, error
                    error('Cannot find specified file: \n%s',fn), 
                end    
            end
            
            
            % load NWB info
            obj.fileconversiondate = char(datetime(datestr(now()),'Format',obj.datetimefmt));
            
            info=h5info(fn);
            if any(strcmp({info.Datasets.Name}, 'nwb_version'))
                nwbv=h5read(fn, '/nwb_version');
            else
                nwbv=h5readatt(fn, '/', 'nwb_version');
            end
         %   nwbv=nwbv{1};
                        
%             if ~strcmp(nwbv, 'NWB-1.0.5')
%                 warnmsg=sprintf('This script was designed for NWB-1.0.5. The loaded file has version %s \n', nwbv{1});
%                 warning(warnmsg)
%             end
            
            if strcmp(nwbv, 'NWB-1.0.5')
                grloc = strcmp({info.Groups(1).Groups.Name}, '/acquisition/timeseries');
                swps={info.Groups(1).Groups(grloc).Groups.Name};
            elseif strcmp(nwbv, '2.2.4')
                grloc = strcmp({info.Groups.Name}, '/acquisition');
                swps={info.Groups(grloc).Groups.Name};
            end
            % dialog option if no stimset or sweep filters are specified
            if ~exist('stimsetfilters', 'var'), stimsetfilters={};end
            if ~exist('sweepselect', 'var'), sweepselect=[];end
            if isempty(stimsetfilters) && isempty(sweepselect)
                protocols=cell(numel(swps), 1);
                for i=1:numel(swps)
                    if strcmp(nwbv, 'NWB-1.0.5')
                        protocols(i) = h5read(fn, [swps{i} '/stimulus_description']);
                    else
                        protocols(i) = h5readatt(fn, swps{i}, 'stimulus_description');
                    end
                end
                list = unique(protocols);
                [indx,tf] = listdlg('ListString',list);
                if tf==1
                    stimsetfilters=list(indx);
                end
            end
            
            grloc=strcmp({info.Groups.Name}, '/stimulus');
            if ~isempty(info.Groups(grloc).Groups(1).Groups)
                stimswps={info.Groups(grloc).Groups(1).Groups.Name};
                swpstimnames={};
                sweeptable=table();
                for i=1:numel(swps)
                    tmp = strsplit(swps{i}, '_');
                    stimdataloc = stimswps{contains(stimswps, tmp{2}) & endsWith(stimswps, ['DA' tmp{3}(end)] ) };
                    if any(contains(stimswps, tmp{2}) & contains(stimswps, 'TTL') )
                        TTLdataloc = stimswps{contains(stimswps, tmp{2}) & contains(stimswps, 'TTL') };
                    else
                        TTLdataloc = '';
                    end
                    % warning: will probably give trouble when multiple TTL waves are used
                    if strcmp(nwbv, 'NWB-1.0.5')
                        protocol = h5read(fn, [swps{i} '/stimulus_description']);
                    else
                        protocol = h5readatt(fn, swps{i}, 'stimulus_description');
                    end
                    sweeptable(i,{'dataloc', 'sweepnr', 'ADname','ADnr', 'protocol', 'stimdataloc', 'TTLdataloc'})=...
                        {swps{i}, str2num(tmp{2})+1, tmp{3},str2double(tmp{3}(end))+1, protocol,...
                        stimdataloc, TTLdataloc};
%                     info.Groups(1).Groups(2).Groups(1).Datasets(3).Dataspace.Size
                end
            else
                error('No stimulus data in nwb. Creating stimulus trace based on description')
            end
           
            if exist('sweepselect', 'var') && ~isempty(sweepselect)
                sweeptable(~ismember(sweeptable.sweepnr,sweepselect),:)=[];
            end

%             for i=1:numel(swps)
%                 swpstimnames(i)=h5read(fn, [swps{i} '/aibs_stimulus_name']);
%             end

            % load Lab notebook
            grloc = strcmp({info.Groups.Name}, '/general');
            grloc2 = strcmp({info.Groups(grloc).Groups.Name}, '/general/labnotebook');
            lbloc = info.Groups(grloc).Groups(grloc2).Groups.Name;
            
            obj.labbooknum=h5read(fn, [lbloc '/numericalValues']);
            obj.labbooknum_keys=h5read(fn, [lbloc '/numericalKeys']);
            
            % find out active headstages
            obj.activeHS = find(any(squeeze(obj.labbooknum(:, strcmp('Headstage Active', obj.labbooknum_keys),:))'));
             
            % LBN=squeeze(LBN(1,:,:))';
            % LBN=cell2table(num2cell(LBN));
%             LBN.Properties.VariableNames=genvarname(LBN_keys(:,1));
%             LBN.Properties.VariableUnits=LBN_keys(:,2);
            obj.labbooktext=h5read(fn, [lbloc '/textualValues']);
%             LBT=squeeze(LBT(1,:,:))';
            obj.labbooktext_keys=h5read(fn, [lbloc '/textualKeys']);
%             LBT=cell2table(LBT);
%             LBT.Properties.VariableNames=genvarname(LBT_keys(:,1));
            obj.labbook_timestamps = squeeze(datetime(obj.labbooknum(1, strcmp('TimeStamp', obj.labbooknum_keys),:), ...
                'ConvertFrom', 'epochtime', 'Epoch', '1904-01-01'));
            obj.labbooknum(:, strcmp('SweepNum', obj.labbooknum_keys),:) = ...
                obj.labbooknum(:, strcmp('SweepNum', obj.labbooknum_keys),:)+1; % since MIES/Igor Pro counts from sweep 0 and we count from sweep 1
               
            % make overview of swps and channels and stimsets
            
                

            
            

            % Add file info
            fileinfo         = dir(fn);
            obj.filename     = fileinfo.name;
            obj.filedirectory= fileDir;
            obj.filetype     = 'Neurodata Without Borders file';
            obj.fileversion  = nwbv;
            obj.filesize     = fileinfo.bytes;
            
            % create a savename
            obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
            
            % Get dates and durations
            obj.filesystemdate= fileinfo.date;
            tmp=h5read(fn, '/file_create_date');
%             obj.filecreatedate= datetime(tmp{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''' ) + duration(1, 0, 0);
%            obj.filecreatedate= datetime(tmp{1}(1:end-6), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS' ) + duration(1, 0, 0);
            
            firstLBentry  = find(any(obj.labbooknum(obj.activeHS, strcmp('SweepNum', obj.labbooknum_keys),:)==1),1,'last'); %last entry sweeps are often re-recorded when new cells are patched
            %firstLBentry = find(LBN.SweepNum==1,1, 'last'); 
%             obj.filetimestart = LBN.TimeStamp(firstLBentry);
%             obj.filetimeend = LBN.TimeStamp(find(~isnat(LBN.TimeStamp),1, 'last'));
%             obj.fileduration  = duration(obj.filetimeend-obj.filetimestart,'Format',obj.durationfmt);
            obj.filetimestart = obj.labbook_timestamps(firstLBentry);
            obj.filetimeend = obj.labbook_timestamps(find(~isnat(obj.labbook_timestamps),1, 'last'));
            obj.fileduration  = duration(obj.filetimeend-obj.filetimestart,'Format',obj.durationfmt);

            % load data from NWBfile according to stimset filters
            if ~iscell(stimsetfilters), stimsetfilters={stimsetfilters}; end
            obj.stimsetfilters=stimsetfilters;
            
            stimselect=contains(sweeptable.protocol,stimsetfilters);
            selectedswps=unique(sweeptable.sweepnr(stimselect));
            [stimsetnms, ~, ic] = unique(sweeptable.protocol(stimselect));
            
            % gather data for each unique stimset and save the sweeps
            fprintf('Adding stimsets and loading sweep data...')
            for i=1:numel(stimsetnms)
                stimsetswpt=sweeptable(stimselect,:);
                stimsetswpt = stimsetswpt(ic==i & ismember(stimsetswpt.sweepnr, selectedswps),:) ;

                obj=obj.addstimset(stimsetnms{i}, stimsetswpt);
%                 obj.stimsets(i).sweepnrs = stimsetswpnrs;
%                 obj.stimsets(i).datalocs = swps(stimsetswpnrs);
%                 obj.stimsets(i).stimdatalocs = stimswps(stimsetswpnrs);
%                 obj.stimsets(i).wavename = stimsetnms{i};
                obj.stimsets(i) = obj.stimsets(i).addnwbchannels; % add ADs and sweep and epochs etc
            end

            
            obj.nrofsweeps   = numel(swps);
            obj.nrofstimsets = numel(stimsetnms);
            obj.updatephase=2;
            
            % Get holding current or voltage
%             for i=1:obj.nrofchannels
%                 if ~isempty(obj.getchannel(i).getin('signal','secondary')) && ~isempty(obj.getchannel(i).getout) && ~isempty(obj.getchannel(i).getout.analogwaveformtable) 
%                     [holdI, holdV] = obj.getchannel(i).getout.getholdingIorV(obj.getchannel(i));
%                     obj.channels(i).analogouts = obj.getchannel(i).getout.set('holdingI',holdI);
%                     obj.channels(i).analogouts = obj.getchannel(i).getout.set('holdingV',holdV);
%                 end      
%             end

        end
        
        % ------------------------------------------ HELPER METHODS ---------------------------------------------------------
        function identify(obj)
            % This dissplays a small summary of the contents of the NWB file object. 
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
        
        %make select functions here later for selecting stimsets and sweeps
%         function obj = selectsweep(obj,varargin)
%             % make an ABFFILE object including only specified SWEEP(s).
%             % Note, using the 'end' keyword as in 1:4:end does not work for this function. 
%             % --------------------
%             % See also SELECTITEM
%             if nargin == 1, return; end
%             for i = 1:numel(obj)
%                 for ii = 1:obj(i).nrofchannels
%                     for iii=1:obj(i).getchannel(ii).nrofanalogins
%                         % replace analogin with one filtered based on varargin
%                         obj(i).channels(ii).analogins(iii) = obj(i).getchannel(ii).getin(iii).selectsweep(varargin{:});
%                     end
%                 end
%                 obj(i).nrofsweeps=numel(obj(i).getchannel(ii).getin(iii).getsweep);
%             end
%         end
        function ap = getstimset(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'stimsets',0,varargin{:});
        end
        
 
        % ----------------------------------------- UPDATING METHODS --------------------------------------------------------        
        function obj = analyseNWB(obj)
%             analyse sweeps of stimsets in the NWBfile object
%             checks for action potentials and calculates passive properties for steps.
            for i = 1:numel(obj), 
                switch obj(i).updatephase
                    case 1
                        error('The stimsets, sweeps and epochs have not been loaded on the NWBfile "%s" .',obj(i).filename);
                    case 2
                        fprintf('analysing NWB %s\n',obj(i).filename)
                        obj(i).stimsets    = obj(i).getstimset.analysestimset; 
                        obj(i).updatephase = 3;
                    case 3
                        error('NWBfile object "%s" has already been analysed.',obj(i).filename);
                end
            end
        end
        function obj =addstimset(obj,name, swpt) %stimsetnms, stimsetswpnrs
            if ~isscalar(obj), error('NWBfile object must be scalar.'); end
            
            stimsetswpnrs = unique(swpt.sweepnr);
            
            stimsetidx = squeeze(ismember(obj.labbooknum(1, strcmp('SweepNum', obj.labbooknum_keys),:),stimsetswpnrs));
            
            firstLBtime = obj.labbook_timestamps(find(obj.labbooknum(1, strcmp('SweepNum', obj.labbooknum_keys),:)==min(stimsetswpnrs),1,'first'));
            firstLBtime = firstLBtime - seconds(10);
            
            timeidx = squeeze(obj.labbook_timestamps>=firstLBtime);
            stimsetidx = stimsetidx & timeidx;
            stimsetLB = obj.labbooknum(:,:,stimsetidx);
            
            activeHS = find(any(squeeze(stimsetLB(:, strcmp('Headstage Active', obj.labbooknum_keys),:))'));
            associated_stimsets = unique(swpt.protocol);
            
            
            starttime = nanmin(obj.labbook_timestamps(stimsetidx));
            endtime = nanmax(obj.labbook_timestamps(stimsetidx));
            stim2add = Stimset('name', name,'sweeptable', swpt,'activeHS', activeHS, 'filename', obj.filename,...
                'filedirectory', obj.filedirectory,'associated_stimsets', associated_stimsets,'sweepnrs',stimsetswpnrs,...
                'datetimestart', starttime,'datetimeend', endtime, 'labbooknum', stimsetLB, 'labbooknum_keys',...
                obj.labbooknum_keys, 'labbook_timestamps', obj.labbook_timestamps(stimsetidx),...
                'guid_NWB', obj.guid);
            if isempty(obj.stimsets) 
                obj.stimsets        = stim2add;
            else
                obj.stimsets(end+1) = stim2add;
            end
%             obj = obj.updatestimstats;
        end
        
        function ap = getstimsets(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'stimsets',0,varargin{:});
        end
      
        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function [hfig,ax_handles] = plot(obj,varargin)
            %popup dialog asking which stimset to plot
            obj.getstimsets.plot;
        end       
               
        function hfig = plotanalysis(obj)
            %popup dialog asking which stimset to plot
            obj.getstimsets.plotanalysis;
        end             
    end
               
end