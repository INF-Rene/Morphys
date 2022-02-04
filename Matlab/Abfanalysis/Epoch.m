classdef Epoch < Sharedmethods & Trace
    %EPOCH An Epoch object
    %   Epochs are sections of sweeps and describe the analog output signal within that section. Most common is the epoch 
    %   "step", which means a constant current or voltage of a given amplitude was applied to the cell during this epoch. 
    %   Epochs inherit properties from Sharedmethods, and importantly, from the Matlab timeseries data type. Note time series 
    %   collections (tscollection) cannot be made from Epochs, these behave very unpredictably.
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    properties (Access = public) % general properties
        guid_swp        % guid of parent sweep
        idx             % index of epoch 
        number          % number of epoch (differs from index in starting at 0, defined by pClamp)
        idxstr          % string indicating epoch. As in pClamp, where epochs are labelled 'A','B',...,'J'
        metadataproperties = {  'guid'        % properties that hold relevant metadata
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
                                'gof'
                                'steadystate'
                                'steadystate_diff'
                                'sag' 
                                'rinput'
                                'amplitude'
                                'jitter'
                                'meansignal'
                                };
    end
    
    properties (Access = public) % epochs specifics
        type                % epoch type as a number. 
        typestr             % epoch type as string (e.g. "step" or "ramp")
        amplitude           % amplitude of injected current. For steps: constant, for ramps: final value at end of ramp, for chirps: absolute deviation from zero line (so peak-to-peak amp = 2*amplitude). 
        maxfrequency        % maximum frequency of chirp (chirps always starts at 0!).
        pulsewidth          % Width of train pulses          %RWS
        pulseperiod         % Period (start to start interval) of train pulses          %RWS
        meansignal          % average of raw data %DBH
        jitter              % std of raw data, can be used as measure of jitter if epoch=step with amplitude=0 %DBH
        initlevel           % starting point of ramp. (Determined by previous epoch) ##NOTE!! not sure where and if to update with scaling factor...    
        extfilename  = '';  % name of text file with noise waveform
        extfilescale = [];  % scaling of text file
    end
    
    properties (Access = public) % step analysis props       
        steadystatewin = 20;   % for EpochSteps. Percentage of time from end of epoch to be used to calculate voltage deflection caused by step. 
        runningavwin   = 5;    % size of running average window in ms for filtering trace to get max Vm deflection.
        tau_win        = 80;   % time in ms after onset pulse to use in fitting tau... Now fixed, but could be flexible, for instance until peak hyppol/depol...    
    end
    
    properties (Access = public) % step analysis result properties   
        stepdiff        = nan;   % relevant for step epochs. The difference in output signal units between this epoch and last (only if both epochs are step epochs!); Also, thi is based on protocol, does not incorporate any scaling
        tau_time        = nan;   % time range used for fitting tau
        fitobject       = nan;   % cfit object resulting from fit
        tau             = nan;   % membrane time constant extracted from single exponential fit.
        gof             = nan;   % goodness-of-fit for tau
        steadystate     = nan;   % voltage at steadystate
        steadystate_diff= nan;   % difference in steadystate with last epoch
        sswin           = nan; 
        sag             = nan; 
        rinput          = nan;   % input resisitance, in MOhm.
        vstep           = nan;   % voltage at max deflection
    end
    
%####################################################### METHODS ############################################################
    methods
        %% ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Epoch(epochtab,initlevel)
            
            % call to superclass constructor
            obj = obj@Sharedmethods;
            obj = obj@Trace;
            
            % return empty if no input arguments
            if nargin == 0, return; end
            
            % fix some timeseries details
            obj = obj.set('Data',epochtab.data,'Time',1:numel(epochtab.data));  % Time is set here to default values to keep length of Data and Time matching
            obj.Data = squeeze(obj.Data); % don't know why, but without squeezing end up with 1x1xN matrix...
            if ~isempty(obj.Data)
                obj = setuniformtime(obj,'StartTime',epochtab.strttime,'Interval',1e3/epochtab.samplefreq);
                obj.TimeInfo.Units = 'milliseconds';
                %obj.TimeInfo.StartDate = epochStruct.strttime;
                obj.DataInfo.Units = epochtab.units;
            end           
            
            % general info
            obj.guid_swp     = epochtab.guid_swp;
            obj.timespan     = epochtab.timespan;
            obj.idx          = epochtab.idx;
            obj.number       = epochtab.number;
            obj.idxstr       = epochtab.idxstr;
            obj.Name         = sprintf('Epoch %s',obj.idxstr);
            obj.datetimestart= epochtab.datetimestart;
            obj.datetimeend  = epochtab.datetimeend;
            obj.samplefreq   = epochtab.samplefreq;
            
            % epoch specification
            obj.type         = epochtab.type;  
            obj.typestr      = epochtab.typestr;
            obj.amplitude    = epochtab.firstlevel;
            obj.maxfrequency = epochtab.maxfrequency;
            obj.pulsewidth   = epochtab.pulsewidth; %RWS
            obj.pulseperiod  = epochtab.pulseperiod;   %RWS
            obj.jitter       = nanstd(obj.Data); %DBH
            obj.meansignal   = nanmean(obj.Data); %DBH
            if nargin == 2, obj.initlevel = initlevel; end
            
            if ~isempty(obj.typestr)
                switch lower(obj.typestr)
                    case 'noisepp',      obj.extfilename  = 'NoiseWave_interpol.mat';
                                         obj.extfilescale = 1.5; % default pre-scaling for noise pp
                    case 'noisespiking', obj.extfilename  = 'NoiseWave_interpol.mat';
                                         obj.extfilescale = 3; % default pre-scaling for noise spiking 
                    case 'truenoise',    obj.extfilename  = 'TrueNoise_interpol.mat';
                                         obj.extfilescale = 1; % default pre-scaling for noise spiking 
                end
            end
        end
        
        %% --------------------------------------- GET/SET/SELECT METHODS ----------------------------------------------------
        function prop = get(obj,propertyname)
            % Function returns any property matching string 'PropertyName' from object (scalar or non-scalar). Note this is
            % calling the 'get' method of timeseries superclass. Thus, we ignore the superclass 'get' method from 
            % Sharedmethods, which would otherwise give conflicting definitions.
            prop = cell(size(obj));
            for i = 1:numel(obj); prop{i} = get@timeseries(obj(i),propertyname); end
            if numel(prop) == 1, prop = prop{:}; end
        end
        
        function obj = set(obj,varargin)
            % Set any property matching string 'PropertyName' from object (scalar or non-scalar). Note this is
            % calling the 'set' method of timeseries superclass. Thus, we ignore the superclass 'set' method from 
            % Sharedmethods, which would otherwise give conflicting definitions.
            for i = 1:numel(obj), obj(i) = set@timeseries(obj(i),varargin{:}); end
        end

        function wf = waveform(obj)
            % plot epoch waveform 
            for i = 1:numel(obj)
                switch lower(obj(i).typestr)
                    case 'step',                                wf = waveform_step(obj(i));
                    case 'ramp',                                wf = waveform_ramp(obj(i));
                    case 'train',                               wf = waveform_train(obj(i));    %RWS      
                    case 'chirp',                               wf = waveform_chirp(obj(i));
                    case {'noisepp','noisespiking','truenoise'},wf = waveform_noise(obj(i));
                    otherwise, error('Epoch type "%s" unknown or not implemented yet',obj(i).typestr)
                    % so these will include the "pulse train", "biphasic train", "triangle train", and "cosine train" epoch types...
                end
            end
        end
        
        function wf = waveform_step(obj)
            % get waveform of step epoch as a timeseries object
            wf = linspace(obj.amplitude,obj.amplitude,milliseconds(obj.timespan)*obj.samplefreq*1e-3)';
            wf = timeseries(wf);
            wf = setuniformtime(wf,'StartTime',0,'Interval',1e3/obj.samplefreq);
        end

        function wf = waveform_train(obj)   %Added by RWS 22-09-2017
            % get waveform of train epoch as a timeseries object
            pulse = linspace(obj.amplitude,obj.amplitude,obj.pulsewidth)';
            gap = linspace(obj.initlevel,obj.initlevel,obj.pulseperiod-obj.pulsewidth)';
            pulsenr = round((milliseconds(obj.timespan)*obj.samplefreq*1e-3-obj.pulsewidth)/obj.pulseperiod);
            wf=[repmat([pulse; gap],pulsenr,1);pulse];
            wf = timeseries(wf);
            wf = setuniformtime(wf,'StartTime',0,'Interval',1e3/obj.samplefreq);
        end
        
        function wf = waveform_ramp(obj)
            % get waveform of ramp epoch object as a timeseries
            wf = linspace(obj.initlevel,obj.amplitude,milliseconds(obj.timespan)*obj.samplefreq*1e-3)';
            wf = timeseries(wf);
            wf = setuniformtime(wf,'StartTime',0,'Interval',1e3/obj.samplefreq);
        end
        
        function wf = waveform_chirp(obj, samplefreq)
            % get waveform of chirp epoch object as a timeseries object.
            wf = obj.chirp(samplefreq);
            wf = timeseries(wf);
            wf = setuniformtime(wf,'StartTime',0,'Interval',1e3/obj.samplefreq);
        end
        
        function wf = waveform_noise(obj)
            % get waveform of noise epoch object as a timeseries
            if obj.samplefreq~=8000, warning('warning, noise waveforms are specifically made for an 8 kHz sampling frequency'); end
            load(fullfile(obj.path2noise, obj.extfilename));
            wf = dataI*obj.extfilescale; 
            wf = timeseries(wf);
            wf = setuniformtime(wf,'StartTime',0,'Interval',1e3/obj.samplefreq);
        end
        
        function [data, time] = chirp(obj)
            % make a sine sweeping through frequencies (Chirp). Chirp always starts at 0 Hz.
            time  = 0:1/obj.samplefreq:seconds(obj.timespan);
            const = obj.maxfrequency/seconds(obj.timespan);    % sets the slope of linear progression through frequencies
            data  = sin(const*pi*time.^2);                     % make data
            data  = (data*obj.amplitude)';                     % scale data
        end
        
        %% -------------------------------------- ANALYSE EPOCH METHODS -----------------------------------------------------
        function obj = analyseepoch(obj,apsaswell,prev_steadystate,prev_amp)
            % analyse epoch
            % Wrapper function that analyses EPOCHs. If apsaswell==1, analyses ACTIONPOTENTIALs as well.
            
            % See also ACTIONPOTENTIAL, GETTAU, GETSAG, GETSTEADYSTATE.
            
            if ~isscalar(obj),error('Object must be scalar'); end
            if nargin < 2, apsaswell=1; prev_steadystate = []; prev_amp = []; end
            if nargin == 2, prev_steadystate = []; prev_amp = []; end
            
            % analyse action potentials
            if apsaswell == 1,
                obj = obj.analyseaps; 
            end

            % 211116 removed this part due to errors with analyzing abf 
%             if ~isempty(prev_amp)
%                 obj.stepdiff = obj.amplitude - prev_amp ;
%             end
            
            % analyse passive properties
            if obj.nrofaps == 0 
                [obj.steadystate, obj.sswin] = obj.getsteadystate;
                if ~isempty(prev_steadystate)
                    obj.steadystate_diff = obj.getsteadystatediff(prev_steadystate);
                end
                if obj.stepdiff ~= 0 && obj.number > 0 % obtain rinput and fit tau to all deflections
                    obj.rinput = obj.getrinput;
                    [obj.tau, obj.fitobject, obj.tau_time, obj.gof] = gettau(obj);
                end
                if obj.stepdiff < 0 % get sag only from negative-going steps
                    [obj.sag, obj.vstep] = obj.getsag;
                end
            end
        end
        
        function ssdiff = getsteadystatediff(obj,prev_steadystate)
            if ~isscalar(obj),error('Object must be scalar'); end
            ssdiff = obj.getsteadystate - prev_steadystate;
        end
        
        function ri = getrinput(obj)
            if ~isscalar(obj),error('Object must be scalar'); end
            ri = [];
            if ~isempty(obj.steadystate_diff)
                switch obj.DataInfo.units
                    case 'pA', ri = abs(obj.stepdiff / obj.steadystate_diff * 1e3);
                    case 'mV', ri = abs(obj.steadystate_diff / obj.stepdiff * 1e3);
                    otherwise, error('unable to calculate ri with units: %s',obj.DataInfo.units);
                end
            else
                %disp('steadystate_diff is empty, ri cannot be calculated')
            end
        end
        
        function [ss, sswin] = getsteadystate(obj)
            % get a default steadystate value and window, using steadystatewin. Intended for step epochs 
            % NOTE: for more flexible selection of epoch sections, use the timeseries method 'getsampleusingtime'.
            for i = 1:numel(obj)
                if ~strcmpi(obj(i).typestr,'step') &  ~strcmpi(obj(i).typestr,'square pulse')  & ~strcmpi(obj(i).typestr,'MIES testpulse')
                    warning('the getsteadystate method was intended for use with step epochs and was not tested for other epoch types. Be careful when interpreting results')
                end
                winstart = obj(i).TimeInfo.Start + range([obj(i).TimeInfo.End,obj(i).TimeInfo.Start])*(100-obj(i).steadystatewin)/100;
                winend   = obj(i).TimeInfo.End;
                sswin    = [winstart,winend];
                ss = mean(obj(i).getsampleusingtime(winstart,winend));
            end
        end
        
        function [sag,vstep] = getsag(obj)
            % calculate sag 
            % Sag is defined as mV difference between initial max deflection and vm at end of current pulse. Assumes pulse is
            % hyperpolarising. 
            sag   = zeros(size(obj));
            vstep = zeros(size(obj));
            for i = 1:numel(obj)
                vstep(i) = min(obj(i).medianfilter(0.5,'truncate').getdata);
                sag(i)   = vstep(i) - obj(i).getsteadystate; % use medianfilter to filter out capacitive peak and use truncate option to avoid the 'zero padding'
                sag(i)   = abs(sag);
            end
        end
        
        function [tau, fitresult, tau_time, gof] = gettau(obj, fitstart, fitend)
            %CREATEFIT(XXX,YYY)
            %  Create a fit.
            %
            %  Data for 'untitled fit 1' fit:
            %      X Input : xxx
            %      Y Output: yyy
            %  Output:
            %      fitresult : a fit object representing the fit.
            %      gof : structure with goodness-of fit info.
            %
            %  See also FIT, CFIT, SFIT.

            %  Auto-generated by MATLAB on 21-Jun-2017 14:46:55
            %  Further modified by Thijs
            
            if isscalar(obj)
                if ~strcmpi(obj.typestr,'step') &  ~strcmpi(obj.typestr,'square pulse')  & ~strcmpi(obj.typestr,'MIES testpulse')
                    warning('the gettau method was intended for use with step epochs and was not tested for other epoch types. Be careful when interpreting results')
                end

                % check inputs
                if nargin == 1
                    fitstart = obj.TimeInfo.Start+0.1;
                    if obj.stepdiff<0
                        [~,fitborder] = min(obj.medianfilter(0.5,'truncate').getdata);
                        fitborder=(obj.Time(fitborder)-obj.TimeInfo.Start)*0.9;
                        if fitborder<obj.tau_win && obj.stepdiff<0
                            fitend=obj.TimeInfo.Start+fitborder;  %end of fit can never be too close to the max hyperpol, otherwise Ih current may influence the result
                        else
                            fitend = obj.TimeInfo.Start+obj.tau_win;
                        end
                    else
                        fitend = obj.TimeInfo.Start+obj.tau_win;
                    end
                    
                elseif nargin == 2, error('please provide a start and end time for fitting'),           
                end
                
                % set defaults
                f_status = 'Failed'; fitresult = nan; tau=nan; tau_time = nan; gof=nan;

                % select time range
                if fitend   > obj.TimeInfo.End,   fitend   = obj.TimeInfo.End; end
                if fitstart < obj.TimeInfo.Start, fitstart = obj.TimeInfo.Start+0.1; end
                
                % check if time range is long enough to make sense for fitting
                if diff([fitstart,fitend]) < obj.minrange4fitting; return; end
                
                % select data using time
                fitts = obj.getsampleusingtime(fitstart,fitend);
                [xdata, ydata] = prepareCurveData( double(fitts.Time-fitts.TimeInfo.Start), double(fitts.Data) );
           
                % attempt fit Set up fittype and options.
                try
                ft = fittype( 'a*exp(b*-x)+c', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display     = 'Off';
                opts.Lower       = [-Inf -100 -200];
                opts.Upper       = [Inf  100 200];
                opts.MaxFunEvals = 1000;
                opts.MaxIter     = 1000;
                opts.StartPoint  = [ydata(1) 0 fitts.getsteadystate];

                % Fit model to data.
                [fitresult, gof] = fit( xdata, ydata, ft, opts ); f_status = 'Passed';
                catch err;
                end
                if strcmp(f_status,'Passed'),
                    tau      = 1/fitresult.b;
                    gof=gof.rsquare;
                    if tau < 5 || tau > 100 || gof<0.8        
                        tau=NaN;
                    end
                    tau_time = xdata;
                end 
            else
                if numel(unique(arrayfun(@(x) obj(x).TimeInfo.Start,1:numel(obj))))>1 % if not all start times are the same
                    error('When using non-scalar 1xN epoch array input, all epochs must have the same start time')
                end
                if nargin == 1,
                    tau = arrayfun(@(x) obj(x).gettau,1:numel(obj),'UniformOutput',false);
                else
                    tau = arrayfun(@(x) obj(x).gettau(fitstart,fitend),1:numel(obj),'UniformOutput',false);
                end
                tau = [tau{:}]; 
            end
        end
        
        %% -------------------------------------- ANALYSE EPOCH METHODS -----------------------------------------------------
        function plotwithfit(obj,varargin)
            % plot tau fits of epoch
            % varargin are passed on to the plot of the regular trace. Fits are plotted red with linewidth 2.
            % See also FITME, EPOCH, ADDTAUFIT2PLOT.
            for i = 1:numel(obj)
                hold on
                    obj(i).plot(varargin{:})
                    obj(i).addtau2plot
                hold off
            end
            grid on
        end
        
        function addtau2plot(obj)
            % add fit data to an existing plot
            for i = 1:numel(obj)
                hold on
                    if ~isempty(obj(i).fitobject) && ~isnan(obj(i).tau)
                        plot(obj(i).tau_time+obj(i).TimeInfo.Start,feval(obj(i).fitobject,obj(i).tau_time),'color','r','linewidth',2)
                    end
                hold off
            end
        end
        
        function addsag2plot(obj)
            % add sag to an existing plot
            for i = 1:numel(obj)
                hold on
                    if ~isempty(obj(i).sag) && ~isempty(obj(i).vstep)
                        line(obj(i).sswin,[obj(i).vstep, obj(i).vstep],            'color','r','linewidth',1,'linestyle','-')
                        line(obj(i).sswin,[obj(i).sag,   obj(i).sag]+obj(i).vstep, 'color','r','linewidth',1,'linestyle','--')
                    end
                hold off
            end
        end
        
        function addss2plot(obj)
            % add steadystate to an existing plot
            for i = 1:numel(obj)
                hold on
                     line(obj(i).sswin,[obj(i).steadystate, obj(i).steadystate],'color','b','linewidth',2,'linestyle','-')
                hold off
            end
        end
        
        function plotanalysis(obj)
            % plot EPOCH with analysis details (e.g. baselines, ap peak events, tau fits and sags etc...). Note to plot 
            % these, the EPOCH must be analysed first with the ANALYSEEPOCH method.
            % See also ANALYSEEPOCH
            for i = 1:numel(obj)
                hold on
                    obj(i).plotapevents
                    obj(i).addtau2plot
                    %obj.addsag2plot
                    obj(i).addss2plot
                hold off
            end
        end
        
    end
end

