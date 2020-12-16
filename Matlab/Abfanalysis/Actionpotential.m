classdef Actionpotential < Sharedmethods
    %ACTIONPOTENTIAL Action potential object
    %   Detailed explanation goes here
    
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
    %###################################################### PROPERTIES ##########################################################
    properties (Access = public)
        parent_guid         % guid of parent sweep OR epoch. depends on where ap was created.
        ts                  % timeseries
        start_time          = nan;% time of AP start (ms)
        end_time            = nan;% time of AP end (ms)
        peak                = nan;% absolute peak amplitude (mV)
        peak_time           = nan;% time of peak (ms)
        ahp                 = nan;% fast component of after-hyperpolarising potential (mV)
        ahp_time            = nan;% time of fast ahp (ms)
        ahp_slow            = nan;% slow component of after-hyperpolarising potential (mV)
        ahp_slow_time       = nan;% time of slow ahp(ms)
        relahp              = nan;% relative fast ahp (ahp-thresh; mV)
        relahp_slow         = nan;% relative slow ahp (ahp-thresh; mV)
        adp                 = nan;% after-depolarising potential (mV)
        adp_time            = nan;% adp time (ms)
        reladp              = nan;% relative adp (adp-thresh; mV)
        thresh              = nan;% threshold potential (mV)
        thresh_time         = nan;% threshold time (ms)
        amp                 = nan;% amplitude (peak-thresh; mV)
        halfwidth           = nan;% width of AP at half maximal amplitude (amp*0.5+thresh; ms)
        halfwidth_strt_time = nan;% time when halfheight is first reached
        halfwidth_end_time  = nan;% time when halfheight is last visited
        maxdvdt             = nan;% maximum dVdt during rising phase of AP (mV/ms)
        maxdvdt_time        = nan;% time of maximum dvdt
        mindvdt             = nan;% minimum dVdt during repolarisation phase of AP (mV/ms)
        mindvdt_time        = nan;% time of minimum dvdt
        upstroke            = nan;% rising slope from 30 to 70% between amplitude and threshold
        downstroke          = nan;% falling slope from 70 to 30% between ampltiude and threshold voltage
        onsetrapidity       = nan;% onset rapidity, slope of phase-plane when crossing dvdt threshold (defined in apThreshRapidity), usually set at 20 or 30 mV/ms (1/ms)
        onsetrapfit   = [nan,nan];% coefficients of linear fit. Use these in polyval to plot.
        onsetrapvm          = nan;% vm in centre of onset rapidity fit window
        isi                 = nan;% inter-spike-interval (ms)
        freq                = nan;% instantaneous firing frequency (Hz)
        number              = nan;% number of AP in set
    end
    
    %####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Actionpotential(varargin)
            obj = obj@Sharedmethods;
            if nargin == 0, return % returns empty AP object
            elseif mod(numel(varargin),2)~=0
                error('Uneven number of name/value pairs.')
            elseif ~all(cellfun(@ischar,varargin(1:2:end)))
                error('All property names must be strings')
            else
                for i=1:2:numel(varargin)
                    obj = obj.set(varargin{i},varargin{i+1});
                end
            end
        end
        
        function e = eventme(obj, feature)
            % transform Actionpotential object into tsdata.event type object, ready to be added as an event to any timeseries
            % object (e.g. Sweep or Epoch).
            %
            % See also TSDATA.EVENT, ACTIONPOTENTIAL, SWEEP, EPOCH
            
            if isscalar(obj)
                switch feature
                    case 'peak',     t = 'peak_time';
                    case 'ahp',      t = 'ahp_time';
                    case 'ahp_slow', t = 'ahp_slow_time';
                    case 'adp',      t = 'adp_time';
                    case 'thresh',   t = 'thresh_time';
                    otherwise, error('Unknown AP feature requested')
                end
                if ~isempty(obj.(t)) && ~isnan(obj.(t))
                    e = tsdata.event(feature,obj.(t));
                else
                    e = tsdata.event;
                end
                e.EventData = obj.(feature);
                e.Units = 'milliseconds';
            else
                e = arrayfun(@(x) obj(x).eventme(feature),1:numel(obj),'UniformOutput',false);
                e = [e{:}];
            end
        end
        
        function data = getdata(obj)
            % get ap trace data from its timeseries
            %
            % See also TIMESERIES
            data = obj.ts.Data;
        end
        
        function data = gettime(obj)
            % get ap trace time from its timeseries
            %
            % See also TIMESERIES
            data = obj.ts.Time;
        end
        
        function dvdt = getdvdt(obj,n)
            % make nth derivative of AP waveform.
            % Based on matlab function 'diff'. obj.getDvdt(n) returns nth difference, so n=2 would give second order
            % differential for example. When n not provided, returns 1st derivative. Note length(dvdt) = length(data)-n
            %
            % See also DIFF
            
            if nargin==1, n=1; end
            if isscalar(obj)
                dvdt = diff(obj.getdata, n)./diff(obj.gettime);
            else
                dvdt = arrayfun(@(x) obj(x).getdvdt(n),1:numel(obj),'UniformOutput',false);
                dvdt = reshape(dvdt,size(obj));
            end
        end
        
        function obj = getdvdtts(obj,n)
            % get n-th order derivative of ap and store as its timeseries
            % Based on matlab function 'diff'. obj.getDvdt(n) returns nth difference, so n=2 would give second order
            % differential for example. When n not provided, returns 1st derivative. Note length(dvdt) = length(data)-n
            %
            % See also GETDVDT, DIFF
            
            if nargin==1, n=1; end
            if isscalar(obj)
                tmpts = obj.ts.delsample('Index',obj.ts.Length-(n-1):obj.ts.Length); % shorten timeseries
                tmpts.Data = obj.getdvdt(n);
                if strcmp(tmpts.DataInfo.Units,'mV') && n==1
                    tmpts.DataInfo.Units = 'mV/ms';
                else
                    tmpts.DataInfo.Units = '';
                end
                obj.ts = tmpts; % replace old with derived
            else
                obj = arrayfun(@(x) obj(x).getdvdtts(n),1:numel(obj),'UniformOutput',false);
                obj = [obj{:}];
            end
        end
        
        function pp = phaseplane(obj)
            % make phase plane data of an action potential.
            % If numel(obj)>1, returns cell array.
            %
            % See also GETDATA, GETDVDT.
            
            if isscalar(obj),
                d  = obj.getdata;
                pp = cat(2,d(1:end-1),obj.getdvdt);
            else
                pp = arrayfun(@(x) obj(x).phaseplane,1:numel(obj),'UniformOutput',false);
            end
        end
        
        function plot(obj,varargin)
            % plot action potential(s).
            % Supports name/value input arguments for plotting. One extra option is to choose 'superimpose' as extra input
            % argument with peak, thresh, ahp, or adp as corresponding value. In that case, APs are plotted on top of each
            % other, aligned to choice.
            %
            % See also ACTIONPOTENTIAL
            
            % check inputs
            if nargin==1, % simple plot
                hold on
                for i=1:numel(obj)
                    if ~isempty(obj(i).ts)
                        obj(i).ts.plot
                    end
                end
                hold off
            else
                if mod(numel(varargin),2)~=0
                    error('Uneven number of name/value pairs.')
                elseif ~all(cellfun(@ischar,varargin(1:2:end)))
                    error('All property names must be strings')
                elseif strcmp('superimpose',varargin{end})
                    error('please provide an AP feature to align APs. Choose one of following: start, peak, thresh, ahp, or adp');
                else
                    [found,sidx]=ismember('superimpose',varargin(1:2:end));
                end
                
                if found, % then superimpose APs
                    feature = varargin{sidx*2};
                    % check feature validity:
                    if ~ismember(feature,{'start','peak','thresh','ahp','adp'})
                        error('please provide a valid AP feature to align APs. Choose one of following: start, peak, thresh, ahp, or adp');
                    end
                    featuretime = sprintf('%s_time',feature);
                    varargin([sidx*2-1,sidx*2]) = []; % remove superimpose name/value pair from varargin (rest is for plot.m)
                    hold on
                    for i=1:numel(obj)
                        if ~isempty(obj(i).ts) && ~isempty(obj(i).(featuretime)) && ~isnan(obj(i).(featuretime))
                            obj(i).ts.Time = obj(i).ts.Time - obj(i).(featuretime);
                            obj(i).ts.plot(varargin{:})
                        end
                    end
                    hold off
                else
                    hold on
                    for i=1:numel(obj)
                        if ~isempty(obj(i).ts)
                            obj(i).ts.plot(varargin{:})
                        end
                    end
                    hold off
                end
            end
        end
        
        function phaseplot(obj,varargin)
            % plot phase plane of ap.
            % varargin are passed on to plot as name/value pair arguments.
            %
            % See also PHASEPLANE
            
            if isscalar(obj)
                pp = obj.phaseplane;
                plot(pp(:,1),pp(:,2),varargin{:})
            else
                hold on
                arrayfun(@(x) obj(x).phaseplot(varargin{:}),1:numel(obj),'UniformOutput',false);
                hold off
            end
        end
        
        function plotanalysis(obj)
            % plot ACTIONPOTENTIAL with derivative overlay and analysis details
            % See also ACTIONPOTENTIAL, TIMESERIES, TSDATA.EVENT, GETDVDTTS.
            if ~isscalar(obj), error('Object must be scalar'); end
            figure
            set(gcf,'position',[413,177,1120,801],'color',[1 1 1])
            
            % store copy for derivative
            dvdtap = obj.getdvdtts;
            
            subplot(2,2,[1 2])
            % plot ap waveform
            
            if ~isempty(obj.peak_time)   && ~isnan(obj.peak_time),   obj.ts = obj.ts.addevent('peak',obj.peak_time); end
            if ~isempty(obj.thresh_time) && ~isnan(obj.thresh_time), obj.ts = obj.ts.addevent('threshold',obj.thresh_time); end
            if ~isempty(obj.ahp_time)    && ~isnan(obj.ahp_time),    obj.ts = obj.ts.addevent('ahp',obj.ahp_time); end
            obj.ts.Name = 'AP waveform';
            yyaxis left;  obj.plot('linewidth',2)
            
            % threshold line
            try
                line(xlim,[obj.thresh obj.thresh],'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % peak time line
            try
                line([obj.peak_time obj.peak_time],[obj.thresh obj.peak],'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % ahp time line
            try
                line([obj.ahp_time  obj.ahp_time] ,[obj.ahp obj.thresh] ,'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % halfwidth line
            try
                line([obj.halfwidth_strt_time  obj.halfwidth_end_time]  ,[1 1]*obj.amp*0.5+obj.thresh,'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            
            % plot derivative
            dvdtap.ts.Name = 'AP derivative';
            if ~isempty(dvdtap.maxdvdt_time) && ~isnan(dvdtap.maxdvdt_time), dvdtap.ts = dvdtap.ts.addevent('maxdvdt',dvdtap.maxdvdt_time); end
            if ~isempty(dvdtap.mindvdt_time) && ~isnan(dvdtap.mindvdt_time), dvdtap.ts = dvdtap.ts.addevent('mindvdt',dvdtap.mindvdt_time); end
            yyaxis right; dvdtap.ts.plot('color','r')
            
            % labelling
            xlabel(sprintf('Time (%s since trace onset)',obj.ts.TimeInfo.Units))
            title(sprintf('Action potential analysis (#%d)',obj.number))
            grid on
            
            try
                subplot(2,2,3)
                obj.phaseplot
                line(xlim,[0 0],'lineStyle','--','color',[0.8 0.8 0.8])
                line([0 0],ylim,'lineStyle','--','color',[0.8 0.8 0.8])
                ylabel('AP derivative (mV/ms)')
                xlabel('AP waveform (mV)')
                title ('AP phase plot')
                rectangle('position',[obj.onsetrapvm-5, obj.apThreshRapidity-30, 10, 60],'EdgeColor','r')
            catch
            end
            try
                subplot(2,2,4)
                obj.phaseplot('o-')
                xlim([obj.onsetrapvm-5 obj.onsetrapvm+5])
                ylim([obj.apThreshRapidity-30 obj.apThreshRapidity+30])
                line(xlim,[0 0],'lineStyle','--','color',[0.8 0.8 0.8])
                line([0 0],ylim,'lineStyle','--','color',[0.8 0.8 0.8])
                ylabel('AP derivative (mV/ms)')
                xlabel('AP waveform (mV)')
                title('AP phase plot zoom-in')
                
                
                % add fitting
                line(xlim,[obj.apThreshRapidity obj.apThreshRapidity],'lineStyle','--','color','r')
                fitx = [obj.onsetrapvm-1 obj.onsetrapvm+1];
                fity = polyval(obj.onsetrapfit,fitx);
                line(fitx,fity,'linewidth',2,'color','r')
            catch
            end
        end
        
        function plotanalysis2(obj)
            % plot ACTIONPOTENTIAL with derivative overlay and analysis details
            % See also ACTIONPOTENTIAL, TIMESERIES, TSDATA.EVENT, GETDVDTTS.
            if ~isscalar(obj), error('Object must be scalar'); end
            
            % store copy for derivative
            dvdtap = obj.getdvdtts;
            % plot ap waveform
            
            if ~isempty(obj.peak_time)   && ~isnan(obj.peak_time),   obj.ts = obj.ts.addevent('peak',obj.peak_time); end
            if ~isempty(obj.thresh_time) && ~isnan(obj.thresh_time), obj.ts = obj.ts.addevent('threshold',obj.thresh_time); end
            if ~isempty(obj.ahp_time)    && ~isnan(obj.ahp_time),    obj.ts = obj.ts.addevent('ahp',obj.ahp_time); end
            obj.ts.Name = 'AP waveform';
            yyaxis left;  obj.plot('linewidth',2)
            
            % threshold line
            try
                line(xlim,[obj.thresh obj.thresh],'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % peak time line
            try
                line([obj.peak_time obj.peak_time],[obj.thresh obj.peak],'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % ahp time line
            try
                line([obj.ahp_time  obj.ahp_time] ,[obj.ahp obj.thresh] ,'linewidth',2,'lineStyle',':','color','k')
            catch
            end
            % halfwidth line
            try
                line([obj.halfwidth_strt_time  obj.halfwidth_end_time]  ,[1 1]*obj.amp*0.5+obj.thresh,'linewidth',2,'lineStyle',':','color','k')
                text(obj.halfwidth_end_time,obj.amp*0.5+obj.thresh,['  HW ' num2str(obj.halfwidth) ' ms'], 'FontSize',10)
            catch
            end
            
            % plot derivative
            dvdtap.ts.Name = 'AP derivative';
            if ~isempty(dvdtap.maxdvdt_time) && ~isnan(dvdtap.maxdvdt_time), dvdtap.ts = dvdtap.ts.addevent('maxdvdt',dvdtap.maxdvdt_time); end
            if ~isempty(dvdtap.mindvdt_time) && ~isnan(dvdtap.mindvdt_time), dvdtap.ts = dvdtap.ts.addevent('mindvdt',dvdtap.mindvdt_time); end
            yyaxis right; dvdtap.ts.plot('color','r')
            
            % labelling
            xlabel(sprintf('Time (%s since trace onset)',obj.ts.TimeInfo.Units))
            title(sprintf('Action potential analysis (#%d)',obj.number))
            grid on
            
            
        end
    end
    
end

