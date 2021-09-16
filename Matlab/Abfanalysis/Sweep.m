classdef Sweep < Sharedmethods & Trace
    %SWEEP A vector of continuously recorded samples 
    %   Sweep objects are a subclass of timeseries with access to all timeseries methods. Further, 
    %   Sweeps inherit a few methods from Sharedmethods. The default time unit is milliseconds. Note time series collections
    %   (tscollection) cannot be made from Sweeps, these behave very unpredictably.
    %   NOTE: when using timeseries methods such as GETSAMPLEUSINGTIME, DELSAMPLE, and ADDSAMPLE, the Sweep properites of
    %   timespan, nrofepochs, aps etc do not update accordingly. 
    %   
    %   TODO:
    %   - as with epochs and APs, keep epochs sorted 
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%% ###################################################### PROPERTIES ##########################################################
    properties (Access = public)
        number        = 0;  % number of sweep
        guid_in             % guid of parent Analogin
        guid_AD
        nrofepochs    = 0;  % number of elements in the epoch list
        sweeplagpnts  = 0;  % the annoying pClamp lag in points (=1/64th of trace length; occurs only when using analog waveform tab, not when using stimulus files)
        sweeplagtime  = duration(0,0,0,0,'Format','hh:mm:ss.SSSSSSS'); % lag in milliseconds, as duration object;
        epochs              % list of Epoch objects. Only populated if the setup configuration info has been added to the parent ABFFile/AnalogIN objects (link between analog OUT and analog IN is required) 
        updatedconfig = 0;  % indication of whether configuration info has been added. 0 = no, 1 = yes;
        labbooknum
        nrinset
        stimwavemode        % 1 if the stimwave is always the same, 2 if the stimwave changes with every sweep
        stimwavename
        stimdataloc
    end
    
%% ##################################################### METHODS ############################################################
    methods
        %% ---------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Sweep(varargin)
            % This is the Sweep constructor function. 

            % call superclass constructors
            obj = obj@Sharedmethods;
            obj = obj@Trace;
            
            % default actions and settings
            obj.datetimestart = datetime; % just a place holder
            obj.samplefreq    = 10000;    % default sample frequency
            
            % check/process inputs
            if nargin == 0, return % returns empty sweep
            elseif mod(numel(varargin),2)~=0
                error('Uneven number of name/value pairs.')
            elseif ~all(cellfun(@ischar,varargin(1:2:end)))
                error('All property names must be strings')
            else
                for i=1:2:numel(varargin)
                    obj = obj.set(varargin{i},varargin{i+1});
                end
            end
            
            if obj.number ~= 0; obj.Name = sprintf('Sweep %d',obj.number); end
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
        function obj = adddata(obj,data,units)
            % add data and units to sweep timeseries. These functions could be improved, very similar to timeseries method
            % 'addsample'.
            obj = obj.set('Data',data,'Time',1:numel(data));  % Time is set here to default values to keep length of Data and Time matching
            for i=1:numel(obj); obj(i).DataInfo.Units = units; end
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
        
        function obj= addNWBepochs(obj, epochtable,fn)
            
            scalefactor=obj.labbooknum.StimScaleFactor;
            scalefactor=unique(scalefactor(~isnan(scalefactor)));
            if ~isempty(scalefactor) && ~isscalar(scalefactor)
                warning('Multiple scalefactors found for this sweep')
                scalefactor=scalefactor(end);
            end
            if ~isempty(scalefactor)
                epochtable.firstlevel=epochtable.firstlevel.*scalefactor;
            end
            
            obj.nrinset=obj.nrinset(end);
            if any(any(epochtable.deltas(:,:)))>0
                epochtable.duration=epochtable.duration+epochtable.deltas(:,1).*(obj.nrinset-1);
                epochtable.firstlevel=epochtable.firstlevel+epochtable.deltas(:,2).*(obj.nrinset-1);
                epochtable.pulseperiod=epochtable.pulseperiod+1000./epochtable.deltas(:,3).*(obj.nrinset-1);
                epochtable.pulsewidth=epochtable.pulsewidth+epochtable.deltas(:,4).*(obj.nrinset-1);
            end
            
            
            for i=1:height(epochtable)
                epochtab=table2struct(epochtable(i,:));
                epochtab.strttime=sum(epochtable.duration(1:i-1));
                epochtab.endtime=epochtab.strttime+epochtable.duration(i);
                
                if obj.Time(end)> epochtab.strttime
                
                    %epoch time info
                    epochtab.datetimestart = obj.datetimestart + duration(0,0,0,epochtab.strttime,'Format',obj.durationfmt);  
                    epochtab.timespan      = duration(0,0,0,epochtab.duration,'Format',obj.durationfmt);
                    epochtab.datetimeend   = epochtab.datetimestart + epochtab.timespan;   

                    %get the data
                    epochts=obj.getsampleusingtime(epochtab.strttime,epochtab.endtime);
                    epochtab.data=epochts.Data;
                    epochtab.units=obj.DataInfo.Units;
                    
                    %get stim data to read input current for "step" epoch;
                    %hotfix for EM data
                    if strcmp(epochtab.idxstr, 'B')
                        stimdata=h5read(fn, [obj.stimdataloc '/data'], double(epochtab.strttime*1e-3*obj.samplefreq),double(epochtab.duration*1e-3*obj.samplefreq));
                        epochtab.firstlevel = nanmedian(stimdata);
                    end

                    %add the epoch to the sweep
                    obj=obj.addepoch(epochtab);
                else
                    break
                end
            end
            
        end
        
        function obj = adddates(obj,h)
            % function to update Sweep properties with info stored in header info ('h')
            if obj.number == 0, error('Please set the number property of Sweep. The original index of sweep (sweep number) is required to calculate datetime start.'); end
            fileDatetimeStart = datetime(datevec(num2str(h.uFileStartDate),'yyyymmdd') + datevec(duration(0,0,0,h.uFileStartTimeMS)),'Format',obj.datetimefmt);
            switch h.nOperationMode
                case 5, % Episodic stimulation mode
                    obj.timespan      = duration(0,0,0,h.sweepLengthInPts*h.si*1e-3,'Format',obj.durationfmt); % duration of sweeps in episodic recording mode
                    obj.datetimestart = fileDatetimeStart + milliseconds(h.sweepStartInPts(obj.number)*h.si*1e-3);
                    obj.datetimeend   = obj.datetimestart+obj.timespan; % end time
                case 3, % Gap-free mode
                    obj.timespan      = duration(0,0,0,h.lActualAcqLength*h.si*1e-3,'Format',obj.durationfmt); % in gap free mode, duration of 'sweep' is equal to entire length of recording
                    obj.datetimestart = fileDatetimeStart; % start time
                    obj.datetimeend   = obj.datetimestart+obj.timespan; % end time
            end
            %obj.TimeInfo.StartDate = fileDatetimeStart;
            %obj.TimeInfo.Format    = obj.datetimefmt;
        end
                
        function obj = addepoch(obj,epochtab)
            % Instantiates an Epoch object using info in the epochStruct and appends this to the epoch list.
            % --------------------
            if ~isscalar(obj),    error('Sweep object must be scalar.'); end
            if numel(epochtab)>1, error('Epoch table may only have one row.'), end
            epochtab.guid_swp   = obj.guid;         % add guid of parent sweep
            epochtab.samplefreq = obj.samplefreq;   % add sampling frequency (actually redundant, since epoch.TimeInfo.Increment has similar information...)
            
            switch lower(epochtab.typestr)
                case {'step','chirp','noisepp','noisespiking','truenoise'}
                            if isempty(obj.epochs), obj.epochs        = Epoch(epochtab);
                            else                    obj.epochs(end+1) = Epoch(epochtab);
                            end
                case lower({'Square pulse','Ramp','Noise','Sin','Saw tooth','Pulse train','PSC','Load custom wave','Combine', 'MIES testpulse'}) %NWB epoch types
                            if isempty(obj.epochs), obj.epochs        = Epoch(epochtab);
                            else                    obj.epochs(end+1) = Epoch(epochtab);
                            end
                case {'ramp', 'train'}                 %RWS
                            if isempty(obj.epochs), obj.epochs        = Epoch(epochtab,obj.getstartendepochwave(obj.nrofepochs,'end'));
                            else                    obj.epochs(end+1) = Epoch(epochtab,obj.getstartendepochwave(obj.nrofepochs,'end'));
                            end
                otherwise, error('Epoch type "%s" unknown or not implemented yet',epochtab.typestr)
                % so these will include the "pulse train", "biphasic train", "triangle train", and "cosine train" epoch types...
            end
            obj = updateepochstats(obj);
        end
        
        %% -------------------------------------- GET/SET/SELECT METHODS ----------------------------------------------------
        function obj = updateepoch(obj,idx,varargin)
            % update properties of an EPOCH in list
            %
            % See also EPOCH 
            for i=1:numel(obj), obj(i) = obj(i).updateitem('epochs',idx,varargin{:}).updateepochstats; end
        end

        function obj = sortepochs(obj,varargin)
            % sort EPOCHs in list
            % 
            % See also EPOCH, SORTITEM, SORT
            for i=1:numel(obj), if obj.nrofaps > 0, obj(i) = obj(i).sortitems('epochs',varargin{:}); end; end
        end
        
        function ap = getepoch(obj,varargin)
            % get EPOCH(s) from list. 
            % Returns EPOCH as a 1xN EPOCH object with N equal to numel(obj). If no input provided, returns all epochs. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also EPOCH, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'epochs',0,varargin{:});
        end
        
        function obj = removeepoch(obj,varargin)
            % remove epoch(s) from list. 
            %
            % see also EPOCH, GETITEM, REMOVEITEM, SELECTEPOCH, SHAREDMETHODS
            for i=1:numel(obj), obj(i) = obj(i).removeitem('epochs',varargin{:}).updateepochstats; end
        end
        
        function obj = selectepoch(obj,varargin)
            % select epoch(s) from list. 
            %
            % see also EPOCH, GETITEM, SELECTITEM, REMOVEEPOCH, SHAREDMETHODS
            for i=1:numel(obj), obj(i) = obj(i).selectitem('epochs',varargin{:}).updateepochstats; end
        end
        
        %% -------------------------------------- ANALYSE SWEEP METHODS -----------------------------------------------------
        function obj = analysesweep(obj,epochsaswell)
            % analyse sweep
            % Wrapper function that analyses ACTIONPOTENTIALs, then assigns them to EPOCHs based on peak time. Then adds
            % the step differences between step epochs.
            % See also ACTIONPOTENTIAL, EPOCHBORDERS, ANALYSEAPS, APS2EPOCHS, EPOCHSTEPS
            if nargin<2, epochsaswell=1; end
            fprintf('\t\t\tanalysing sweeps\n')
            obj = obj.analyseaps;
            if obj(1).nrofepochs>0
                obj = obj.aps2epochs;
                obj = obj.epochsteps;
                if epochsaswell, obj = obj.analyseepochs(0); end
            end
        end
        
        function obj = analyseepochs(obj,apsaswell)
            if nargin<2, apsaswell=1; end
            for i=1:numel(obj)
                prev_steadystate = [];
                prev_amp = [];
                for ii=1:obj(i).nrofepochs
                    if ii>1, prev_steadystate = obj(i).getepoch(ii-1).steadystate; prev_amp = obj(i).getepoch(ii-1).amplitude; end % collect previous baseline
                    obj(i).epochs(ii) = obj(i).getepoch(ii).analyseepoch(apsaswell,prev_steadystate,prev_amp);
                end
            end
        end
        
        function b = epochborders(obj)
            % get EPOCH border times
            % See also EPOCH, TIMESERIES
            if~isscalar(obj), error('object must be scalar'); end
            b = [obj.getepoch.getstarttime, obj.getendtime];
        end
        
        function obj = aps2epochs(obj)
            % assign ACTIONPOTENTIALs in list to EPOCHs in list
            % See also ACTIONPOTENTIAL, EPOCH, TIMESERIES
            for i=1:numel(obj)
                if obj(i).nrofaps>0
                    epochnumbers = discretize([obj(i).getap.peak_time],obj(i).epochborders);
                    if any(isnan(epochnumbers)), error('action potentials found that do not fall within an epoch boundary... Time for debuggin'); end
                    for ii=1:obj(i).nrofaps
                        obj(i).epochs(epochnumbers(ii)) = obj(i).getepoch(epochnumbers(ii)).addap(obj(i).getap(ii), 'parent_guid', obj(i).epochs(epochnumbers(ii)).guid); % added guid argument RWS 20171003
                    end
                end
            end
        end
        
        function obj = epochsteps(obj)
            % determine the size and direction of the step between adjacent step epochs and add this value to epochs in list
            % Used for automated analysis of step epoch based protocols. 
            %       stepdiff_inj:   difference in injected signal between steps, uses analog waveform tab info.
            % See also UPDATEEPOCH, EPOCH.
            for i=1:numel(obj)
                for ii=2:obj(i).nrofepochs
                    if all(ismember(obj(i).getepoch([ii-1,ii]).get('typestr'),'step')) % if both adjacent epochs are step epochs
                        obj(i) = obj(i).updateepoch(ii,'stepdiff' , obj(i).getepoch(ii).amplitude - obj(i).getepoch(ii-1).amplitude);
                    end
                end
            end
        end
        

        %% -------------------------------------- MISCELLANEOUS METHODS -----------------------------------------------------
        function value = getstartendepochwave(obj,idx,inputstr)
            % function to calculate the start or end level of an Epoch (at position 'idx' in Epoch list). That is, what is 
            % injected during the first or last sample of this epoch? In practice, this function is used to determine the 
            % 'initial level' for an EpochRamp object, by taking the last level of previous epoch. 
            if ~isscalar(obj), error('Sweep object must be scalar.'); end
            if idx == 0, value = 0; return, end
            wf = obj.getepoch(idx).waveform.Data;
            switch inputstr
                case 'start', value = wf(1);
                case 'end'  , value = wf(end);
                otherwise
                    error('Only "end" and "start" are valid inputs')   
            end
        end
        
        function wf = waveform(obj)
            % function to make analog output waveform as calculated from the epoch collection in epoch list. Basically
            % concatenates all individual epoch waveforms.
            % See also TIMESERIES, APPEND
            if ~isscalar(obj), error('Sweep object must be scalar.'); end
            if obj.nrofepochs == 0 || isempty(obj.nrofepochs)
                disp('no epoch, no waveform.')
                wf=[];
            else
                for i=1:obj.nrofepochs
                    if i==1 
                        wf = obj.getepoch(i).waveform;
                    else
                        nextwf = obj.getepoch(i).waveform;
                        nextwf = nextwf.set('Time',nextwf.Time+obj.getepoch(i).TimeInfo.Start);
                        wf = wf.append(nextwf);
                    end
                end
                wf.TimeInfo.Units = 'milliseconds';
            end
        end        
        
        %% ----------------------------------------- UPDATING METHODS -------------------------------------------------------
        function obj = updateepochstats(obj)
            for i=1:numel(obj), obj(i).nrofepochs = numel(obj(i).epochs); end
        end
        
        function obj = makeepochs(obj,swptab,swplag)
            % funtion that adds an Epoch object for every epoch in sweepEpochList to the list of epochs of this sweep. swplag
            % is optional, specifies number of samples lag before protocol start (1/64th of total sweep length in points.
            % Does not occur when using stimulus files). Then, duration of leading and lagging epochs can be determined; the
            % leading epoch is the time spent waiting for the analog waveform to begin, the lagging epoch is the time
            % recorded in a sweep after the analog waveform out has stopped. These two epochs are then added to the swptab
            % and added to the list of epochs for this sweep as any other, but only if their duration is greater than zero.
            % The epoch string (idxstr) given to the leading Epoch is "<", the lagging Epoch is ">". 
            
            if ~isscalar(obj), error('Sweep object must be scalar.'); end
            if nargin == 2, 
                warning('No sweep lag provided, assuming no lag.'), 
            else
                obj.sweeplagpnts   = swplag;
                obj.sweeplagtime   = duration(0,0,obj.sweeplagpnts/obj.samplefreq,0,'Format',obj.durationfmt);
            end
            
            % make epoch template
            template = {'idx',          0;
                        'number',       -1;
                        'idxstr',       '-';
                        'type',         1;
                        'typestr',      'step';
                        'timespan',     0;
                        'firstlevel',   0; % assumption that lead and lagging epochs inject 0, not checked but pretty sure. Unless intersweep holding is set differently perhaps...
                        'pulseperiod',  0;
                        'pulsewidth',   0;
                        'maxfrequency', 0;
                        'deltascaled',  0;
                        };
            template = struct2table(cell2struct(template(:,2),template(:,1)));

            % add leading epoch
            if obj.sweeplagpnts > 0;
                t = template;
                t.timespan = milliseconds(obj.sweeplagtime);
                t.idxstr   = '<'; 
                swptab = [t;swptab]; % append
            end
            
            % add lagging epoch
            unallocatedtime = round(milliseconds(obj.timespan) - sum(swptab.timespan),3); % so thats to microsecond resolution...
            unallocatedtime = duration(0,0,0,unallocatedtime);
            if milliseconds(unallocatedtime)>0
                t = template;
                t.timespan = milliseconds(unallocatedtime);
                t.idxstr   = '>'; 
                t.idx      = max(swptab.idx)+1;
                t.number   = max(swptab.number)+1;
                swptab = [swptab;t]; % append
            end

            epochdatetimestart = obj.datetimestart;
            for i = 1:height(swptab)               
                if i>1,
                    epochdatetimestart = epochdatetimestart + duration(0,0,0,swptab(i-1,:).timespan);    % increase delay
                end
                epochtab = table2struct(swptab(i,:));
                if  epochtab.timespan > 0 % this condition is necessary to avoid creating 'empty' Epochs (with duration = 0)
                    epochtab.datetimestart = epochdatetimestart;                                         % add start time
                    epochtab.timespan      = duration(0,0,0,epochtab.timespan,'Format',obj.durationfmt); % convert to duration
                    epochtab.datetimeend   = epochtab.datetimestart + epochtab.timespan;                 % add end time

                    strttime = milliseconds(epochdatetimestart-obj.datetimestart);
                    endtime  = strttime + milliseconds(epochtab.timespan);
                    
                    swpsection          = obj.getsampleusingtime(strttime,endtime);
                    epochtab.data       = swpsection.Data';
                    epochtab.units      = obj.DataInfo.Units;
                    epochtab.samplefreq = obj.samplefreq;
                    epochtab.strttime   = strttime;

                    obj = obj.addepoch(epochtab); % add epoch object to list
                end
            end
            obj.updatedconfig = 1;
        end

        %% ----------------------------------------- sTING METHODS -------------------------------------------------------
        function plot(obj,varargin)
            % NOTE when asking a sweep to plot itself, it will always plot its entire timeseries, regardless of wether
            % certain epochs have been removed from list...
            hold on
            for i=1:numel(obj)
                plot@timeseries(obj(i),varargin{:})
                grid on
            end
            hold off
            % some formatting...
            title('Time Series Plot')
            datainfo = [obj.DataInfo]; unitvals = unique({datainfo.Units});
            timeinfo = [obj.TimeInfo]; timevals = unique({timeinfo.Units});
            if numel(unitvals)==1, ylabel(unitvals{:}); end
            if numel(timevals)==1, xlabel(timevals{:}); end

        end   
        
        function plotanalysis(obj)
            % plot SWEEP with analysis details (e.g. baselines, ap peak events, tau fits and sags etc...). Note to plot 
            % these, the SWEEP must be analysed first with the MAKEEPOCHS method. Only adds epoch borders for first Sweep
            % in array.
            % NOTE: this method does not work on sections of Sweep selected using the getsampleusingtime method of 
            % timeseries, it will simply plot the whole sweep because plotting is based on Epochs in list, which are 
            % unaffected by getsampleusingtime. 
            % See also MAKEPOCHS
            for i = 1:numel(obj)
                obj(i).getepoch.plotanalysis
                hold on
                    % add baseline
                    bl = obj(i).getepoch(1).getsteadystate;
                    line(xlim,[bl bl],'color','k','linestyle','-')
                hold off
            end
            obj(1).addepochborders2plot;
        end
        
        function plotwithborders(obj,varargin)
            % plot sweep with epoch borders.
            for i=1:numel(obj)
                obj(i).plot(varargin{:})
                obj(i).addepochborders2plot;
            end
        end 
        
        function addepochborders2plot(obj)
            % plot sweep with epoch borders. Assumes all sweeps have same borders.
            hold on
                b = obj(1).epochborders;
                line([b;b],repmat(ylim',1,numel(b)),'color',[0.5 0.5 0.5],'linestyle','-.','linewidth',1)
            hold off
        end 

    end
    
end

