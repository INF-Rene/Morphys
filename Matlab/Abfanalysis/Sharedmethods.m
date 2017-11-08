 classdef Sharedmethods
    %SHAREDMETHODS A superclass with universal methods
    %   From ABFBatch, through ABFFiles, AnalogInputs, and Sweeps down to Epochs, a number of methods are shared by all.
    %   Using this Sharedmethods class, all levels can share these methods using class inheritance. Also used as a central 
    %   point for storing fixed preset values, paths to common folders etc...
    %   
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    properties (Access = public)
        created 	% date of object creation
        guid    	% assign globally unique identifier
    end
    
    properties (Hidden = true) % formatting
        datetimefmt = 'yyyy-MM-dd HH:mm:ss.SSSSSS'; % formatting of datetime objects, with microsecond precision
        durationfmt = 'hh:mm:ss.SSSSSSS';           % formatting of duration objects, with microsecond precision
    end
    
    properties (Hidden = true) % abf file recording modes currently supported (gap free and episodic)
        implrecmode = [3,5];   % numbers indicating ABF recording modes that currently can be handled
    end
    
    properties (Hidden = true)      % default analysis settings
        apMinPeakHeightVm   = -5;    % (mV) minimum AP height to be detected
        apMinPeakHeightDvdt = 30;   % (mV/ms) for detecting peaks in dvdt trace 
        apThreshDvdt        = 15;   % (mV/ms) dVdt threshold for determining AP threshold. Value of 10 mV/ms usually does well, but perhaps ~20 best for recordings with more high-frequency noise...    
        apThreshRapidity    = 30;   % (mV/ms) onset rapidity will be measured as slope of phase plane when it crosses this threshold.    
        apMinPeakDistance   = 2;    % (ms) refractory period; minimum time between peaks for peak detection.
        apMinPeakProminence = 10;   % (mV) The amount of mV the peak has to stand out from the minimum between the neighboring peaks and/or the beginning/end of sweep
        plateaucorrection   = 1e-6; % small value to fix 'plateau' peaks to make them detectable by findpeaks function.
        ahpSnakeWin         = 3;    % (ms) length of running window to find AHP after AP peak.
        adpSnakeWin         = 15;   % (ms) length of running window to find ADP after AP peak.
        postAPtime          = 100;  % maximum time after AP to be taken into waveform
        preAPtime           = 1;    % time before AP threshold ot be taken into waveform
        dvdtmedfilt         = 0.1;  % size of window for median filter of dvdt traces (ms). Used to remove the capacative peaks from your signal.
        onsetrapfitwin      = 5;    % number of samples (in effect microseconds) to include in the onset rapidity fit, before and after the crossing of the onset rapidity dvdt threshold (apThreshRapidity)
        upsample            = 1e-3; % degree of upsampling for sections of trace (used in halfwidth and onset rapidity calculations)
        minrange4fitting    = 10;   % minimum duration of section of trace that will be fit using gettau function of Epoch objects.
    end
    
%####################################################### METHODS ############################################################
    methods (Access = public)
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Sharedmethods
            % make a SHAREDMETHODS object. 
            
            % set some defaults
            obj.created = datestr(now());                           % date of object creation
            obj.guid    = upper(char(java.util.UUID.randomUUID())); % assign globally unique identifier
                       
        end
        
        % ------------------------------------------- SHARED METHODS --------------------------------------------------------
    
        function prop = get(obj,propertyname)
            % Get object properties. 
            % V = get(obj,'PropertyName') returns the value of the specified property for the object 'obj' (scalar or non-scalar).
            if nargin < 2, error('Too little input')
            elseif ~ischar(propertyname) %|| ~ismember(propertyname,properties(obj))
                error('Please provide a valid property name (as string).')
            end
            if isscalar(obj)
                prop = obj.(propertyname);
            else
                prop = arrayfun(@(x) obj(x).get(propertyname),1:numel(obj),'UniformOutput',false);
                prop = reshape(prop,size(obj));
            end
        end
        
        function obj = set(obj,varargin)
            % Function assigns 'value' to any property matching string 'PropertyName' from object.
            % Supports Property/value pair input arguments.
            if nargin < 3 || mod(nargin-1,2)~=0, error('Property/value pairs must come in even number.'); end
            prop2val = cat(1,varargin(1:2:end),varargin(2:2:end))';
            for i=1:size(prop2val,1)
                propertyname = prop2val{i,1};
                value        = prop2val{i,2};
                for ii = 1:numel(obj), obj(ii).(propertyname) = value; end
            end
        end
        
        function item = getitem(obj,listname,inverselect,varargin)
            % get an item from a list. 
            % Returns item as a 1xN object with N equal to numel(idx). If no input provided, returns all items in list. This 
            % function is intended for use by subclasses which have lists of objects, such as Sweep objects, which have a 
            % list of epochs. The listname would then be 'epochs'.
            % inverselect is a switch to select the inverse of whatever specified. if 0, normal selection, if 1, inverse. 
            % Varargin may contain indexes (real positive integers or logicals)
            % Usage: 
            % The following examples will all return the 2nd and 3rd sweep of the analogin:
            %   s = analogin.getitem('sweeps',0,[2,3])) 
            %   s = analogin.getitem('sweeps',0,[false true true]) 
            %   s = analogin.getitem('sweeps',1,logical([1 0 0]))
            % Items can be also be selected using property name/value pairs, where the property must match a property name of
            % the object in the list. Items will be returned that have a corresponding value, but selection can only be 
            % done on one criterion at a time. This example will return all Sweeps in the Analogin that have 3
            % Actionpotentials:
            %   s = analogin.getitem('sweeps',0,'nrofaps',3))
            
            if nargin < 3, error('please provide the name of the list to select from and choose selection method (0 = normal, 1 = inverse)'); end
            if isscalar(obj), 
                % if list is empty, return empty
                if isempty(obj.(listname)), item = []; return; end
                % if no selection criteria provided (so no varargin), get all 
                if nargin == 3;
                    logicalidxs = true(1,numel(obj.(listname)));
                % if 1 element in varargin, must be logical or index array.
                elseif numel(varargin) == 1
                    indxes = 1:numel(obj.(listname));
                    i2keep = indxes(varargin{1});
                    logicalidxs = logical(ismember(indxes,i2keep));
                % if 2 elements in varargin, must be name/value pair.
                elseif numel(varargin) == 2
                    propname  = varargin{1}; if ~ischar(propname), error('Property names must be strings'); end
                    propvalue = varargin{2};
                    if isnumeric(propvalue) && ~isempty(propvalue)
                        proplist = obj.getitem(listname,0).get(propname);
                        if iscell(proplist), proplist = cell2mat(obj.getitem(listname,0).get(propname)); end
                        logicalidxs = proplist == propvalue;
                    elseif ischar(propvalue)
                        logicalidxs = strcmp(obj.getitem(listname,0).get(propname),propvalue);
                    else
                        error('only properties with numeric or character (string) data types can be selected this way.') 
                    end
                % anything else is an error
                else
                    % selecting multiple properties should be possible to implement. for now however:
                    error('Only one property can be selected at a time when getting an item from list.')
                end
                % now either keep or invert selection.
                if inverselect == 1
                    item = obj.(listname)(~logicalidxs);
                else
                    item = obj.(listname)(logicalidxs);
                end
            else
                % if object is not scalar, do 'em one by one. Note dimensions are lost, and items are returned in a vector.
                if nargin < 3,
                    item = arrayfun(@(x) obj(x).getitem(listname,inverselect),1:numel(obj),'UniformOutput',false);
                else
                    item = arrayfun(@(x) obj(x).getitem(listname,inverselect,varargin{:}),1:numel(obj),'UniformOutput',false);
                end
                item = [item{:}];
            end
        end
                
        function obj = removeitem(obj,listname,varargin)
            % remove item(s) from list. 
            %
            % see also GETITEM, SHAREDMETHODS
            if isscalar(obj)
                items2keep = obj.getitem(listname,1,varargin{:});
                if ~isempty(items2keep), obj.(listname) = items2keep;
                else obj.(listname) = {};
                end
            else
                error('removeitem only handles scalar objects')
            end
        end
        
        function obj = selectitem(obj,listname,varargin)
            % select item(s) from list. 
            %
            % see also GETITEM, SHAREDMETHODS
            if isscalar(obj)
                items2keep = obj.getitem(listname,0,varargin{:});
                if ~isempty(items2keep), obj.(listname) = items2keep;
                else obj.(listname) = [];
                end
            else
                error('selectitem only handles scalar objects')
            end
        end
                
        function obj = updateitem(obj,listname,idx,varargin)
            % update properties of item at location IDX in list
            for i = 1:numel(obj)
                obj(i).(listname)(idx) = obj(i).getitem(listname,0,idx).set(varargin{:});   
            end
        end
        
        function obj = sortitems(obj,listname,property,descending)
            % sort items in list by PROPERTY
            % Returns object with list (LISTNAME) sorted according to value of item PROPERTY. 
            % Default order is ascending, unless DESCENDING is set to 1.
            %
            % Example:
            %   epoch.sortitems('aps')                --> sort aps in ascending order by time of peak
            %   epoch.sortitems('aps','thresh_time')  --> sort aps in ascending order by time of threshold
            %   epoch.sortitems('aps','amp',1)        --> sort aps in descending order by their amplitude 
            %
            % See also GETITEM, REMOVEITEM, UPDATEITEM, SELECTITEM, SORT
            
            if ~isscalar(obj), error('object must be scalar'); end
            if nargin<4, descending = 0; end
            if nargin == 2,
                switch listname
                    case 'abfs',  property = 'fileName';
                    case {'channels','analogouts','analogins','sweeps','epochs'},
                                  property = 'number';
                    case 'aps',   property = 'peak_time';
                    otherwise,    error('listname not recognised')
                end
            end
            nrofitems = numel(obj.(listname));
            if nrofitems > 0, 
                vals = obj.getitem(listname,0).get(property);
                if isempty(vals), warning('no values to sort on, returning list unsorted'); return; end % return object as is, no vaules to sort on.
                if numel(vals)==1 && ~iscell(vals), vals = {vals}; end  % in case only one item is returned, the vals can be of any type... convert to cell for next step.
                obj = obj.removeitem(listname, cellfun(@isempty,vals)); % remove any item with no value for feature
                if nrofitems > numel(obj.(listname)), warning('%d items(s) removed during sorting.',nrofitems-numel(obj.(listname))); end
                [~,i] = sort([obj.getitem(listname,0).(property)]); 
                obj.(listname) = obj.(listname)(i); 
                if descending == 1, obj.(listname) = fliplr(obj.(listname)); end
            end
        end
        
        function s = metadata(obj,metadataproperties)
            % make a struct of object properties relevant to be included in a metadata overview. METADATAPROPERTIES is a
            % cell array of strings that should contain property names of the object. If none provided, the default 
            % metadataproperties of this object are used.
            s = obj.structme;
            if nargin<2,
                switch class(obj)
                    case 'Abfbatch', metadataproperties = { 'guid';'savename';'nrofabfs'};
                    case 'Abffile',  metadataproperties = { 'guid'
                                                            'guid_batch'
                                                            'userid'
                                                            'fileconversiondate' 
                                                            'filesystemdate'
                                                            'filename'
                                                            'filedirectory'
                                                            'filetype'
                                                            'fileversion'
                                                            'filesize'
                                                            'fileduration'
                                                            'filetimestart'
                                                            'filetimeend'
                                                            'datadimorderstr'
                                                            'proname'
                                                            'prodirectory'
                                                            'operationmodestr'
                                                            'samplefreq'
                                                            'stopwatchtime'
                                                            'nrofsweeps'
                                                            'nrofchannels'
                                                            'setupsettingsname'
                                                          };
                    case 'Channel',  metadataproperties = { 'guid'
                                                            'guid_abf'
                                                            'filename'
                                                            'number'
                                                            'nrofanalogins'
                                                            'nrofanalogouts'
                                                          };
                    case 'Analogout',metadataproperties = { 'guid'
                                                            'guid_channel'
                                                            'number'        
                                                            'dacusername'
                                                            'units'
                                                            'epochinfosource'
                                                            'path2pro'
                                                            'protocolfile'
                                                            'path2stim'
                                                            'stimulusfile'
                                                            'scalefactor'
                                                            'holdingI'
                                                            };
                    case 'Analogin', metadataproperties = { 'guid'
                                                            'guid_channel'
                                                            'number' 
                                                            'adcusername'
                                                            'units'
                                                            'telegraphenabled'
                                                            'telegraphinstrument'
                                                            'gain'
                                                            'lowpassfilter'
                                                            'instrumentscalefactor'
                                                            'protocollag'
                                                            'samplefreq'
                                                            'nrofsweeps'
                                                            'signal'
                                                            }; 
                    case 'Sweep',    metadataproperties = { 'guid'
                                                            'guid_in'
                                                            'number'
                                                            'Name'
                                                            'sampleFreq'
                                                            'timespan'
                                                            'datetimestart'
                                                            'datetimeend'
                                                            'nrofepochs'
                                                            'nrofaps'
                                                            'sweeplagpnts'
                                                            'sweeplagtime'
                                                            };
                    case 'Epoch',    metadataproperties = { 'guid'
                                                            'guid_swp'
                                                            'number'
                                                            'idxstr'
                                                            'Name'
                                                            'typestr'
                                                            'datetimestart'
                                                            'datetimeend'
                                                            'samplefreq'          
                                                            'timespan'            
                                                            'nrofaps'
                                                            'stepdiff'
                                                            'vstep'
                                                            'tau'
                                                            'steadystate'
                                                            'steadystate_diff'
                                                            'sag' 
                                                            'rinput'
                                                            };
                    case 'Actionpotential', metadataproperties = {'guid'
                                                            'parent_guid'
                                                            'start_time' 
                                                            'end_time'
                                                            'peak'
                                                            'peak_time'
                                                            'ahp'
                                                            'ahp_time'
                                                            'relahp'
                                                            'adp'
                                                            'adp_time'
                                                            'reladp'
                                                            'thresh'
                                                            'thresh_time'
                                                            'amp'
                                                            'halfwidth'
                                                            'halfwidth_strt_time'
                                                            'halfwidth_end_time'
                                                            'maxdvdt'
                                                            'maxdvdt_time'
                                                            'mindvdt'     
                                                            'mindvdt_time'
                                                            'onsetrapidity'  
                                                            'onsetrapfit'    
                                                            'onsetrapvm'
                                                            'isi'
                                                            'freq'
                                                            'number'
                                                            };
                    otherwise, return
                end
            end
            f = fieldnames(s);
            s = rmfield(s,f(~ismember(f,metadataproperties))); % remove fields that aren't needed
        end
        
        function s = structme(obj,flag)
            % Transform objects in list into an object of type struct. If FLAG==1, make structs of all
            % non-Matlab objects present in 'downstream' lists.
            if nargin < 2, flag=0; end
            if isscalar(obj), 
                s = cell2struct(cellfun(@(x) obj.(x),properties(obj),'UniformOutput',false),properties(obj));
                if flag==1
                    switch class(obj)
                        case 'Abfbatch', lists = {'abfs'};
                        case 'Abffile',  lists = {'channels'};
                        case 'Channel',  lists = {'analogins','analogouts'};
                        case 'Analogin', lists = {'sweeps'};
                        case 'Analogout',lists = {'analogwaveformtable'};
                        case 'Sweep',    lists = {'epochs','aps'};
                        case 'Epoch',    lists = {'aps'};
                        otherwise, return
                    end      
                    for i=1:numel(lists)
                        tmp = s.(lists{i});
                        s   = rmfield(s,lists{i});
                        for ii = 1:numel(tmp)
                            s.(lists{i})(ii) = tmp(ii).structme(flag);
                        end
                    end
                end
            else 
                s = arrayfun(@(x) obj(x).structme(flag),1:numel(obj),'UniformOutput',false);
                s = [s{:}];
            end
        end
        
        function xml = xmlme(obj,destination,fn)
            % This function writes all meta data of this object into an XML file. 
            % check inputs
            error('this function is not ready yet. certain objects cannot be converted to xml yet.')
            if ~isscalar(obj), error('Object must be scalar.'); end
            if nargin<3, fn = [obj.guid '.xml']; end
            if nargin<2, destination=cd(); end
            xml = xml_write(fullfile(destination,fn),obj.structme(1));
        end
        
        function saveme(obj,destination,fn)
            % save the object as .mat file
            % check inputs
            if ~isscalar(obj), error('Object must be scalar.'); end
            if nargin<3, fn = obj.guid; end
            if nargin<2, destination=cd(); end
            fprintf('Saving object as %s.mat ... ',fn)
            save(fullfile(destination,fn),'obj')
            fprintf('done.\n')
        end
    end
end

