classdef Stimwave < Sharedmethods & Trace
    %Stimwave object
    %   A Stimwave object is a group of digital to analog out waveforms that are used in NWB stimsets
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       René Wilbers (renewilbers@gmail.com)
    %   Created:      24-01-2020
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
    %###################################################### PROPERTIES ##########################################################
    properties
        guid_nwbchannel          % Globally unique identifier of parent NWB object.
        filename            % filename of parent NWB
        filedirectory       % directory of parent NWB
        name                % Stimset name
        
        sweepnrs   % to which sweep nrs this wave was applied
        dataloc    % the location of the wave data in the NWB file
        units
        
    end
    
    %####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Stimwave(varargin)
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
            %             if ~isempty(obj.number),
            %                 obj.name = sprintf('Channel %d',obj.number);
            %             end
            
        end
        %% --------------------------------------- GET/SET/SELECT METHODS ----------------------------------------------------
        function q = info(obj)
            disp(obj)
        end
        
        function prop = get(obj,propertyName)
            % Function returns any property matching string 'PropertyName' from object (scalar or non-scalar). 
            % Note this is calling the 'get' method of timeseries superclass. Thus, we ignore the superclass 'get' method 
            % from Sharedmethods, which would otherwise give conflicting definitions.
            prop = get@timeseries(obj,propertyName); 
        end
        
        function obj = set(obj,varargin)
            % Set any property matching string 'PropertyName' from object (scalar or non-scalar). 
            % Note this is calling the 'set' method of timeseries superclass. Thus, we ignore the superclass 'set' method 
            % from Sharedmethods, which would otherwise give conflicting definitions.
            obj = set@timeseries(obj,varargin{:}); 
        end

        %% ------------------------------------------ HELPER METHODS --------------------------------------------------------
        function obj = addts(obj,ts)
            % add data and units to sweep timeseries. These functions could be improved, very similar to timeseries method
            % 'addsample'.
            obj = obj.set('Data',ts.Data,'Time',ts.Time);  % Time is set here to default values to keep length of Data and Time matching
            for i=1:numel(obj)
                obj(i).DataInfo.Units = ts.DataInfo.Units;
                obj(i).TimeInfo.Units = ts.TimeInfo.Units;
            end
        end
        
        function obj = addtime(obj,samplefreq)
            % add time and units to sweep timeseries. Note, only 'within' sweep time is set here, starting from 0 with the 
            % default unit of milliseconds. Absolute start times of sweeps are set in 'adddates' method using h.
            if nargin < 2, error('Please provide sampling frequency.'); end
            for i=1:numel(obj); 
                if isempty(obj.Data), error('Please add data to sweep timeseries first before specifying time.'); end
                obj(i) = obj(i).set('samplefreq',samplefreq).setuniformtime('StartTime',0,'Interval',1e3/samplefreq);
                obj(i).TimeInfo.Units = 'milliseconds';
            end
        end
        % -------------------------------------------- OTHER METHODS --------------------------------------------------------
        
        function ap = getstimwave(obj,varargin)
            % get SWEEP(s) from list.
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps.
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use
            % variable arguments in for selecting items.
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'stimwaves',0,varargin{:});
        end
        
        
        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function plot(obj,varargin)
            % NOTE when asking a sweep to plot itself, it will always plot its entire timeseries, regardless of wether
            % certain epochs have been removed from list...
            hold on
            for i=1:numel(obj)
                plot@timeseries(obj(i),varargin{:})
                grid on
            end
            yl = ylim ;
            yrange = yl(2)-yl(1) ;
            yrange2add = yrange*0.1 ;
            ylim([yl(1)-yrange2add yl(2)+yrange2add])
            hold off
            % some formatting...
            title('Time Series Plot')
            datainfo = [obj.DataInfo]; unitvals = unique({datainfo.Units});
            timeinfo = [obj.TimeInfo]; timevals = unique({timeinfo.Units});
            if numel(unitvals)==1, ylabel(unitvals{:}); end
            if numel(timevals)==1, xlabel(timevals{:}); end

        end   
    end
end

