classdef Spike2Channel < Sharedmethods
    %Channel Object describing amplifier channel
    %   Each amplifier channel has analog outputs and analog inputs. Here, the mapping of these inputs and outputs onto the
    %   BNCs of the axon digidata is specified. Also contains the lists of corresponding Analogin and Analogout objects.
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ########################################################## 
    properties
        guid_abf            % Globally unique identifier of parent Abffile object.
        filename            % filename of parent abffile
        name                % Channel name (e.g. "Channel 2")
        number              % number of this Channel (make it to correspond with the experimenter's numbering of amplifier channels)
        
        primary             % Number of the primary analogin in the spike2 generated mat file
        secondary           % Number of the secondary analogin in the spike2 generated mat file
        
        trigger             % Name of a trigger channel. Not necessarily unique associated with this amplifier channel, but good to have the option, for instance, if one records the Master8 output triggering an EPSP in one of the Analogins.
        updatephase = 1     % phase of class; (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)   
        
        analogins           % list of Analogins  (= Analogin  class instances)
        analoginnrs         % numbers of analog input  channels (corresponding to numbers of digidata)
        nrofanalogins  = 0  % number of Analog Input channels recorded 
    end
    
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Spike2Channel(varargin)
            % call to superclass constructor. Takes a struct or name/value pair arguments as inputs.
            obj = obj@Sharedmethods;
                        
            % check/process inputs
            if     numel(varargin) == 0, return % returns empty Channel object   
            elseif numel(varargin) == 1
                s = varargin{:};
                if isstruct(s)
                    flds = fieldnames(s);
                    for i=1:numel(flds)
                        obj = obj.set(flds{i},s.(flds{i}));
                    end
                else
                    error('Please provide a struct or name/value pairs as inputs')
                end
            elseif mod(numel(varargin),2)~=0
                error('Uneven number of name/value pairs.')
            elseif ~all(cellfun(@ischar,varargin(1:2:end)))
                error('All property names must be strings')
            else
                for i=1:2:numel(varargin)
                    obj = obj.set(varargin{i},varargin{i+1});
                end
            end
            if ~isempty(obj.number),
                obj.name = sprintf('Channel %d',obj.number);
            end

        end
        % -------------------------------------------- OTHER METHODS --------------------------------------------------------
        function n = getanaloginputnrs(obj)
            % get list of all POSSIBLE/DEFAULT analog input numbers used in this CHANNEL (sorted).
            if isscalar(obj)
                n = [obj.get('primary'),obj.get('secondary'),obj.get('trigger')];
            elseif isvector(obj)
                n = arrayfun(@(x) obj(x).getanaloginputnrs,1:numel(obj),'UniformOutput',false);
                n = [n{:}];
            else
                error('Object must be scalar or vector')
            end
            n = sort(n);
        end        
        
        function obj = addin(obj,infostruct,channeldata)
            % Instantiates an Analogin object and appends this to the analog IN list of the Channel object.
            if ~isscalar(obj),error('Object must be scalar.'); end
            
            % find signal for this in
            in2add = Spike2Analogin(infostruct,channeldata,obj.guid);
            if isempty(obj.analogins), 
                 obj.analogins        = in2add;
            else obj.analogins(end+1) = in2add;
            end
            obj = obj.updateinstats;
        end
        
        function obj = addout(obj,outinfo)
            % Instantiates an Analogout object using information in OUTinfo and appends this to the analog OUT list of the Channel object.
            if ~isscalar(obj),     error('Object must be scalar.'); end
            if ~isstruct(outinfo), error('Input must be a struct.'); end
            
            out2add = Analogout(outinfo,obj.guid);
            if isempty(obj.analogouts), obj.analogouts  = out2add;
            else                  obj.analogouts(end+1) = out2add;
            end
            obj = obj.updateoutstats;
        end     
        
        function obj = updatein(obj,idx,varargin)
            % update properties of ANALOGIN in list
            %
            % See also UPDATEIO, ANALOGIN
            obj = updateio(obj,'analogins',idx,varargin{:});
        end
        
        function obj = updateout(obj,idx,varargin)
            % update properties of ANALOGOUT in list
            %
            % See also UPDATEIO, ANALOGOUT
            obj = updateio(obj,'analogouts',idx,varargin{:});
        end
        
        function obj = updateio(obj,listname,idx,varargin)
            % update properties of ANALOGIN/OUT in list. 
            % If no index and varargin provided, just calls the updateinstats/updateoutstats methods.
            %
            % See also ANALOGIN, ANALOGOUT
            if ~isscalar(obj), error('object must be scalar'); end
            if nargin > 2
                obj = obj.updateitem(listname,idx,varargin{:});
            end
            switch listname
                case 'analogins', obj = obj.updateinstats; % update stats for analog ins
                case 'analogouts',obj = obj.updateoutstats;
                otherwise, error('unknown list')
            end
        end

        function obj = sortins(obj,varargin)
            % sort ANALOGIN in list
            % See also SORTIO, ANALOGIN.
            obj = sortios(obj,'analogins',varargin{:});
        end
        
        function obj = sortouts(obj,varargin)
            % sort ANALOGOUT in list
            % See also SORTIO, ANALOGOUT.
            obj = sortios(obj,'analogouts',varargin{:});
        end        
        
        function obj = sortios(obj,listname,varargin)
            % sort ANALOGIN/OUT in list
            % Returns object with list sorted according to sweep FEATURE. 
            % Default feature is number, default order is ascending, unless DESCENDING is set to 1.
            % Usage
            %   obj = sortsweeps(OBJ,FEATURE,DESCENDING)
            %
            % Example:
            %   obj.sortios()               --> sort ANALOGIN/analogouts in ascending order by ADC/DACnumber
            %   obj.sortios('number')       --> sort ANALOGIN/analogouts in ascending order by ADC/DACnumber
            %   obj.sortios('guid',1)       --> sort ANALOGIN/analogouts in descending order by their guid 
            %
            % See also ANALOGIN, ANALOGOUT, SORTITEM, SORT
            
            if ~isscalar(obj), error('object must be scalar'); end
            if obj.nrofaps > 0, 
                 obj = obj.sortitems(listname,varargin{:});
            end
        end
        
        function out = getout(obj,varargin)
            % get ANALOGOUT(s) from list. 
            % Returns ANALOGOUT as a 1xN ANALOGOUT object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also ANALOGOUT, GETITEM, SHAREDMETHODS
            out = getitem(obj,'analogouts',0,varargin{:});
        end        
        
        function in = getin(obj,varargin)
            % get ANALOGIN(s) from list. 
            % Returns ANALOGIN as a 1xN ANALOGIN object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also ANALOGIN, GETITEM, SHAREDMETHODS
            in = getitem(obj,'analogins',0,varargin{:});
        end
        
        function obj = removeio(obj,listname,varargin)
            % remove ANALOGIN/OUT(s) from list. 
            %
            % see also ANALOGIN, ANALOGOUT, GETITEM, REMOVEITEM, SELECTEPOCH, SHAREDMETHODS
            for i = 1:numel(obj), obj(i) = obj(i).removeitem(listname,varargin{:}).updateio(listname); end
        end
        
        function obj = removeout(obj,varargin)
            % remove ANALOGOUT(s) from list. 
            %
            % see also ANALOGOUT, GETITEM, REMOVEITEM, SELECTEPOCH, SHAREDMETHODS
            obj = removeio(obj,'analogouts',varargin{:});
        end        
        
        function obj = removein(obj,varargin)
            % remove ANALOGIN(s) from list. 
            %
            % see also ANALOGIN, GETITEM, REMOVEITEM, SELECTEPOCH, SHAREDMETHODS
            obj = removeio(obj,'analogins',varargin{:});
        end
        
        function obj = selectio(obj,listname,varargin)
            % select ANALOGIN/OUT(s) from list. 
            %
            % see also ANALOGIN, ANALOGOUT, GETITEM, SELECTITEM, SHAREDMETHODS
            for i = 1:numel(obj), obj(i) = obj(i).selectitem(listname,varargin{:}).updateio(listname); end
        end
        
        function obj = selectout(obj,varargin)
            % select ANALOGOUT(s) from list. 
            %
            % see also SELECTIO, ANALOGOUT, GETITEM, SELECTITEM, SHAREDMETHODS
            obj = selectio(obj,'analogouts',varargin{:});
        end        
        
        function obj = selectin(obj,varargin)
            % select ANALOGIN(s) from list. 
            %
            % see also SELECTIO, ANALOGIN, GETITEM, SELECTITEM, SHAREDMETHODS
            obj = selectio(obj,'analogins',varargin{:});
        end
        
        function obj = selectsweep(obj,varargin)
            % make a CHANNEL object including only specified SWEEP(s).
            % Note, using the 'end' keyword as in 1:4:end does not work for this function. 
            % --------------------
            % See also SELECTITEM
            if nargin == 1, return; end
            for i = 1:numel(obj)
                for ii=1:obj(i).nrofanalogins
                    % replace analogin with one filtered based on varargin
                    obj(i).analogins(ii) = obj(i).getin(ii).selectsweep(varargin{:});
                end
            end
        end
        
        function obj = makeepochs(obj)
            %
            if ~isscalar(obj), error('Abffile object must be scalar.'); end
            fprintf('- Updating channel %d...\n',obj.number)           
            if ~isempty(obj.getout) && ~isempty(obj.getout.analogwaveformtable)                                                       
                % update every Analogin in Channel
                for i = 1:obj.nrofanalogins
                    % replace old IN object in list with updated one. If epochinfosource is 'pClampProFile', a lag will be
                    % added to protocol. 
                    addlag = strcmp(obj.getout.epochinfosource,'pClampProFile');
                    obj.analogins(i) = obj.getin(i).makeepochs(obj.getout.analogwaveformtable, addlag);                        
                end
                obj.updatephase = 2;
            else
                warning('No Analogout object present in channel, no protocol/epoch information to update.')
            end
        end
        
        function obj = analysechannel(obj)
            % analyse Analogins with signal type "primary" of Channel.
            % checks for action potentials and calculates passive properties for steps.
            for i = 1:numel(obj), 
                switch obj(i).updatephase
                    case 1
                        warning('Channel in update phase 1 cannot be analysed. Most likely protocol information is not available.');
                    case 2
                        fprintf('\tanalysing channel %d\n',obj(i).number)
                        obj(i).analogins   = obj(i).getin.analysein; 
                        obj(i).updatephase = 3;
                    case 3
                        error('Channel object "%s" has already been analysed.',obj(i).name);
                    otherwise
                        error('Unknown option for updatephase')
                end
            end
        end
        
        function s = getsignalfromanaloginnr(obj,analoginnr)
            % get the signal type (primary, secondary, ... ) of Analogin object using Analogin number (number of the ADC
            % channel of digidata). 
            % NOTE: it may be better to have the Analogin array made below as a property of the Setupsettings object...
            analoginarray = {   'primary'   ,obj.primary;
                                'secondary' ,obj.secondary;           
                                'scope'     ,obj.scope;        
                                'mode'      ,obj.mode;
                                'trigger'   ,obj.trigger;
                            };
            analoginarray = analoginarray(~cellfun(@isempty,analoginarray(:,2)),:);            
            s = analoginarray{ismember(cell2mat(analoginarray(:,2)),analoginnr),1};
        end
        
        function obj = updateinstats(obj)
            % function updates the count of analog inputs associated with this abffile object and updates the list
            % of analog input numbers.
            % --------------------
            for i = 1:numel(obj)
                obj(i).nrofanalogins = numel(obj(i).analogins); % update IN count
                obj(i).analoginnrs = cell2mat(arrayfun(@(x) obj(i).getin(x).number,1:obj(i).nrofanalogins,'UniformOutput',false)); % update list of IN numbers
            end
        end
        
        function obj = updateoutstats(obj)
            % function updates the count of analog outputs associated with this abffile object and updates the list
            % of analog output numbers.
            % --------------------
            for i = 1:numel(obj)
                obj(i).nrofanalogouts = numel(obj(i).analogouts); % update OUT count
                obj(i).analogoutnrs = cell2mat(arrayfun(@(x) obj(i).getout(x).number,1:obj(i).nrofanalogouts,'UniformOutput',false)); % update list of OUT numbers
            end
        end

        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function [hfig,ax_handles] = plot(obj,varargin)
            % Plots the Channel object. Returns figure handle.
            % --------------------
            if isscalar(obj)
                ax_handles = zeros(obj.nrofanalogins,1);
                for i=1:obj.nrofanalogins
                    ax_handles(i) = subplot(obj.nrofanalogins,1,i); 
                    analogin2plot = obj.getin(i);
                    analogin2plot.plot(varargin{:})
                    if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
                        configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.number);
                    else
                        configstring = '';
                    end
                    title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
                    ylabel(analogin2plot.units)
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
                hfig = gcf;
                set(hfig, 'color', [1 1 1])
            else
                arrayfun(@(x) obj(x).plot(varargin{:}),1:numel(obj),'UniformOutput',false)
            end
        end       

        function hfig = plotanalysis(obj)
            % Plots the Channel object with Epoch boundaries and results of basic analysis. Returns figure handle.
            % --------------------
            if isscalar(obj)
                ax_handles = zeros(obj.nrofanalogins,1);
                for i=1:obj.nrofanalogins
                    ax_handles(i) = subplot(obj.nrofanalogins,1,i); 
                    analogin2plot = obj.getin(i);
                    analogin2plot.plotanalysis
                    if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
                        configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.number);
                    else
                        configstring = '';
                    end
                    title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
                    ylabel(analogin2plot.units)
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
                hfig = gcf;
                set(hfig, 'color', [1 1 1])
            else
                arrayfun(@(x) obj(x).plotanalysis,1:numel(obj),'UniformOutput',false)
            end
        end
    end
end

