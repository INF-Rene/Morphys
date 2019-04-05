classdef Trace < timeseries
    %TRACE Ephys trace
    %   A trace is a piece of a recorded signal. Used as a superclass to define methods to be inherited by both the the Sweep 
    %   and Epoch classes. 
    %   
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    properties
        datetimestart
        datetimeend
        samplefreq              % sampling frequency (Hz)
        timespan                % duration of the trace. Chose timespan as name because 'duration' is already a matlab class
        aps                     % list of action potentials
        nrofaps = 0
    end
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Trace
            % contruct a Trace
            %
            % See also TIMESERIES
            
            obj = obj@timeseries;
        end

        %% ------------------------------------------ OTHER METHODS -------------------------------------------------------
        function ts = getstarttime(obj)
            ti = [obj.TimeInfo];
            ts = reshape([ti.Start],size(obj));
        end
        
        function te = getendtime(obj)
            ti = [obj.TimeInfo];
            te = reshape([ti.End],size(obj));
        end        
        
        function obj = addap(obj, varargin)
            % add an AP to the list. 
            % Option to provide an ACTIONPOTENTIAL object as input. In that case, no other optional arguments (varargin) are 
            % allowed.
            % See also ACTIONPOTENTIAL
            
            if ~isscalar(obj), error('object must be scalar'); end
            if any(cellfun(@isa,varargin,repmat({'Actionpotential'},1,numel(varargin)))) % if any input is an actionpotential
                if ~any(ismember(varargin(2:2:end), {'parent_guid'})) % added RWS 20171003 to allow guid to be passed in the function. This is to prevent sweep guids to be used as parent guids when APs are passed from sweeps to epochs
                    varargin={varargin{1:end}, 'parent_guid', obj.guid};
                end
                if numel(varargin)==1 || numel(varargin)==3
                    ap = varargin{1};
                    ap=set(ap, 'parent_guid', varargin{3});
                else
                    error('Only 1 or 3 (in case of guid passing) input arguments expected when providing Actionpotential object as input.') 
                end
            else
                varargin={varargin{1:end}, 'parent_guid', obj.guid};
                ap = Actionpotential(varargin{:});  % make the ap
            end
            
            % add it 
            if isempty(obj.aps), obj.aps        = ap; 
            else                 obj.aps(end+1) = ap;  
            end
            % update stats
            obj = obj.updateapstats;           
        end

        function obj = updateap(obj,idx,varargin)
            % update properties of an AP in list
            %
            % See also ACTIONPOTENTIAL 
            
            % updateitem can handle non-scalar objects, so no for loop necessary:
            obj = obj.updateitem('aps',idx,varargin{:}).updateapstats;
        end
        
        function obj = sortaps(obj,varargin)
            % sort APs in list
            % Returns object with list sorted according to ap FEATURE. 
            % Default feature is peak time, default order is ascending, unless 'descending' is set to 1.
            % NOTE this is peak_time relative to onset sweep, absolute datetime sorting not yet implemented, could be done by
            % adding sweep start datetime to peak time...
            % Usage
            %   obj = sortaps(OBJ,FEATURE,DESCENDING)
            %
            % Example:
            %   obj.sortaps()               --> sort aps in ascending order by time of peak
            %   obj.sortaps('thresh_time')  --> sort aps in ascending order by time of threshold
            %   obj.sortaps('amp',1)        --> sort aps in descending order by their amplitude 
            %
            % See also ACTIONPOTENTIAL, SORTITEM, SORT
            for i = 1:numel(obj); if obj.nrofaps > 0, obj(i) = obj(i).sortitems('aps',varargin{:}); end; end
        end
        
        function ap = getap(obj,varargin)
            % get an action potential from list. 
            % Returns APs as a 1xN Actionpotential object with N equal to numel(idx). If no input provided, returns all aps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also ACTIONPOTENTIAL, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'aps',0,varargin{:});
        end
        
        function obj = removeap(obj,varargin)
            % remove ap(s) from list. 
            %
            % see also ACTIONPOTENTIAL, GETITEM, REMOVEITEM, SELECTAP, SHAREDMETHODS
            for i = 1:numel(obj); obj(i) = obj(i).removeitem('aps',varargin{:}).updateapstats; end
        end
        
        function obj = selectap(obj,varargin)
            % select ap(s) from list. 
            %
            % see also ACTIONPOTENTIAL, GETITEM, SELECTITEM, REMOVEAP, SHAREDMETHODS
            for i = 1:numel(obj); obj(i) = obj(i).selectitem('aps',varargin{:}).updateapstats; end
        end 
        
        function obj = updateapstats(obj)
            % update ap counts and order by peak time. 
            %
            % see also ACTIONPOTENTIAL.
            for i = 1:numel(obj)
                obj(i).nrofaps = numel(obj(i).aps);
                if obj(i).nrofaps>0
                    [~,idx] = sort([obj(i).getap.peak_time]); 
                    obj(i).aps = obj(i).aps(idx); % sort APs by peak time 
                end
            end
        end
        
        %% ----------------------------------------- ANALYSIS METHODS -------------------------------------------------------
        function d = getdata(obj)
            % get trace data. If numel(obj)>1, length of TRACE objects must be equal
            % See also TIMESERIES
            d = ones(obj(1).Length,numel(obj)); for i = 1:numel(obj); d(:,i) = obj(i).Data; end
        end
        
        function t = gettime(obj)
            % get trace time.
            % See also TIMESERIES
            t = ones(obj(1).Length,numel(obj)); for i = 1:numel(obj); t(:,i) = obj(i).Time; end
        end
        
        function d = getdvdt(obj,n)
            % make nth derivative of trace. 
            % Based on matlab function 'diff'. obj.getDvdt(n) returns nth difference, so n=2 would give second order 
            % differential for example. When n not provided, returns 1st derivative. Note length(dvdt) = length(data)-n 
            %
            % See also DIFF, GETDATA, TIMESERIES
            if nargin==1, n=1; end
            d = diff(obj.getdata,1,n)./diff(obj.gettime,1,n);
        end  
        
        function obj = getdvdtts(obj,n)
            % get Epoch derivative as timeseries. Based on matlab function 'diff'. obj.getDvdt(n) returns nth difference, 
            % so n=2 would give second order differential for example. Note length(dvdt) = length(data)-n 
            % See also GETDVDT, TIMESERIES, DELSAMPLE.
            if nargin==1, n=1; end
            
            for i=1:numel(obj)
                dvdt   = obj(i).getdvdt(n);
                obj(i) = obj(i).delsample('Index',obj(i).Length-(n-1):obj(i).Length); % shorten timeseries
                obj(i).Data = dvdt;
                % fix units.
                if strcmp(obj(i).DataInfo.Units,'mV') && n==1
                    obj(i).DataInfo.Units = 'mV/ms';
                else
                    obj(i).DataInfo.Units = '';
                end
            end
        end
        
        function obj = getsampleusingtime(obj,varargin)
            % get section of sweep using time. Overload of TIMERSERIES' GETSAMPLEUSINGTIME method to operate on non-scalar
            % timeseries objects. 
            %
            % See also: TIMERSERIES, GETSAMPLEUSINGTIME
            for i=1:numel(obj), 
                obj(i) = getsampleusingtime@timeseries(obj(i),varargin{:}); 
            end
        end
        
        function m = meantime(obj,starttime,endtime)
            % get mean value of data between a start and end time (relative to trace onset)
            %
            % See also MEAN, TIMESERIES, GETSAMPLEUSINGTIME
            m = zeros(size(obj));
            for i=1:numel(obj); m(i) = mean(obj(i).getsampleusingtime(starttime, endtime)); end
        end
        
        function s = stdtime(obj,starttime,endtime)
            % get standard deviation of data between a start and end time (relative to trace onset)
            %
            % See also STD, TIMESERIES, GETSAMPLEUSINGTIME
            s = zeros(size(obj));
            for i=1:numel(obj); s(i) = std(obj(i).getsampleusingtime(starttime, endtime)); end
        end  
        
        function obj = runmean(obj,runwin,varargin)
            % apply a running average over timeseries data using a window of specified size (milliseconds). 
            % Common use: rm = obj.runmean(30). The varargin optional input can be used to access other features of the
            % runmean function. 
            %
            % See also RUNMEAN, GETDATA, TIMESERIES. 
            for i=1:numel(obj);  
                runwin = floor(runwin/(1e3/obj(i).samplefreq));
                obj(i) = obj(i).set('Data',runmean(obj(i).getdata, runwin, varargin{:}));
            end
        end
        
        function obj = medianfilter(obj,filtwin,varargin)
            % apply a one dimensional median filter on timeseries data.
            % This function acts as a running median filter, using a window of size filtwin (milliseconds). The varargin 
            % optional input can be used to access other features of the medfilt1 function. 
            %
            % See also MEDFILT1, GETDATA, TIMESERIES.
            for i = 1:numel(obj)
                if nargin == 1, obj(i) = obj(i).set('Data',medfilt1(obj.getdata));
                else
                    filtwin = floor(filtwin/(1e3/obj(i).samplefreq));
                    obj(i) = obj(i).set('Data',medfilt1(obj(i).getdata, filtwin, varargin{:}));
                end
            end
        end
        
        function s = avtrace(obj)
            % average of Sweep objects. Note this is different from the TIMESERIES native MEAN method, which returns a mean
            % over time, so that's 1 value. Returned sweep objects retains all properties of 1st Sweep object in array.
            s = obj(1).set('Data',mean(obj.getdata,2));
        end     
        
        function s = stdtrace(obj)
            % standard deviation of Sweep objects. Returned sweep objects retains all properties of 1st Sweep object in array.
            s = obj(1).set('Data',std(obj.getdata,0,2));
        end
        
        function s = semtrace(obj)
            % standard error of the mean of Sweep objects. Returned sweep objects retains all properties of 1st Sweep object in array.
            s = obj(1).set('Data',std(obj.getdata,0,2)/sqrt(numel(obj)));
        end
        
        %% ---------------------------------------- AP ANALYSIS METHODS -----------------------------------------------------
        function obj = analyseaps(obj)
            % find APs and all currently available features
            for i = 1:numel(obj)   
                %fprintf('Finding APs...\n')
                obj(i) = obj(i).findaps;
                %fprintf('Getting thresholds...\n')
                obj(i) = obj(i).getthresh;     % get initial estimate to determine maxdvdt a.o.
                %fprintf('maxdvdts...\n')
                obj(i) = obj(i).getmaxdvdt;   
                %fprintf('Threshold2...\n')
                obj(i) = obj(i).getthresh2;    % get final estimate based on Allen institute method
                %fprintf('Relative amps...\n')
                obj(i) = obj(i).getrelamp;
                %fprintf('HWs...\n')
                obj(i) = obj(i).gethalfwidth;
                %fprintf('Fast AHPs...\n')
                obj(i) = obj(i).getahp;
                %fprintf('Slow AHPs...\n')
                obj(i) = obj(i).getslowahp;
                %fprintf('Relative AHP...\n')
                obj(i) = obj(i).getrelahp;
                %fprintf('mindvdt...\n')
                obj(i) = obj(i).getmindvdt;
                %obj(i) = obj(i).getadp;
                %obj(i) = obj(i).getreladp;
                %fprintf('Onsetrapidity...\n')
                obj(i) = obj(i).getonsrapidity;
                %fprintf('startend...\n')
                obj(i) = obj(i).getapstartend;
                %fprintf('apnrs...\n')
                obj(i) = obj(i).findapnr;
                %fprintf('isifreqaps...\n')
                obj(i) = obj(i).isifreqaps;
                %fprintf('adding AP TS...\n')
                obj(i) = obj(i).addapts;                    
            end
        end
        
        function obj = findaps(obj)
            % find Actionpotentials in epoch. 
            % Returns a 1xN Actionpotential array. Uses findpeaks to find APs.
            %
            % See also FINDPEAKS
            if ~isscalar(obj), error('Object must be scalar'); end
            warning('off', 'signal:findpeaks:largeMinPeakHeight') % fisrt turn off the annoying findpeaks warnings that occur when no peaks are found.
            [pks, locs] = findpeaks(obj.plateaufix(obj.getdata),'MinPeakHeight',obj.apMinPeakHeightVm,'minpeakdistance',floor(obj.apMinPeakDistance*obj.samplefreq*1e-3), 'MinPeakProminence', obj.apMinPeakProminence);
            for i=1:numel(pks)
                obj = obj.addap('peak',pks(i),'peak_time',(locs(i)-1)*1e3/obj.samplefreq+obj.TimeInfo.Start);
            end
        end
     
        function obj = getthresh(obj)
            % find the threshold of action potentials using derivative of timeseries. 
            % Starts at maximum of rising phase of current AP and finds sample before AP derivative first drops below 
            % apThreshDvdt (as seen from peak, going back in time), which is usually set at 10 or 20 mV/ms. 
            % NOTE: Assumes list of ACTIONPOTENTIALS is sorted!
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;       
            for i = 1:obj.nrofaps                                                                               % for every AP event
                if i == 1, 
                     strt = obj.TimeInfo.Start;                                                                 % find start point.
                else strt = obj.getap(i-1).peak_time;
                end
                isi       = obj.getsampleusingtime(strt,obj.getap(i).peak_time);                                % get timeseries of trace before AP peak
                isidvdt   = isi.getdvdtts.medianfilter(obj.dvdtmedfilt);                                        % get derivative. Perform a modest median filter to selectively filter out capacative transients
                strtidx   = find(isidvdt.Data==max(isidvdt),1);                                                 % choose max dvdt of rising phase of current AP
                threshidx = strtidx - find(flipud(isidvdt.Data(1:strtidx))<obj.apThreshDvdt,1)+1;               % calculate threshold index.
                if ~isempty(threshidx)
                    obj   = obj.updateap(i,'thresh',isi.Data(threshidx),'thresh_time',isi.Time(threshidx));     % update AP
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no threshold.\n',cnt,obj.nrofaps); end
        end
        
        function obj = getthresh2(obj)  
            % Djai's version (based on allen institute method)
            % find the threshold of action potentials using derivative of timeseries. 
            % similar to the one used above, only difference is that the threshold for dv/dt is set as 
            % 5 percent of the max dv/dt (averaged over all APs in trace)
            % NOTE: Assumes list of ACTIONPOTENTIALS is sorted!
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;       
            maxdvdts = [] ;
            
            for i = 1:obj.nrofaps 
                if ~isempty(obj.getap(i).maxdvdt)
                maxdvdts(i) = obj.getap(i).maxdvdt ;
                end
            end
            
            if ~isempty(maxdvdts)
                Mmaxdvdt = nanmean(maxdvdts) ;
            else
                Mmaxdvdt = NaN ;
            end
            
            for i = 1:obj.nrofaps                                                                               % for every AP event
                if i == 1 
                     strt = obj.TimeInfo.Start;                                                                 % find start point.
                else strt = obj.getap(i-1).peak_time;
                end
                isi       = obj.getsampleusingtime(strt,obj.getap(i).maxdvdt_time);                                % get timeseries of trace before AP peak
                isidvdt   = isi.getdvdtts.medianfilter(obj.dvdtmedfilt);                                        % get derivative. Perform a modest median filter to selectively filter out capacative transients
                strtidx   = find(isidvdt.Data==max(isidvdt),1);                                                 % choose max dvdt of rising phase of current AP
                
                if i == 1   || obj.getap(i).peak_time-obj.getap(i-1).peak_time > 1500   %added 1500 ms condition for endurance analysis                                                                         % calculate threshold index.
                threshidx = strtidx - find(flipud(isidvdt.Data(1:strtidx)) < 0.05*obj.getap(i).maxdvdt,1)+1;    % for first AP (with usually high max dv/dt)           
                else
                threshidx = strtidx - find(flipud(isidvdt.Data(1:strtidx)) < 0.05*Mmaxdvdt,1)+1;                % average max dv/dt for other APs
                end

                threshidx = strtidx - find(flipud(isidvdt.Data(1:strtidx)) < 0.05*obj.getap(i).maxdvdt,1)+1;
                
                if ~isempty(threshidx)
                    obj   = obj.updateap(i,'thresh',isi.Data(threshidx),'thresh_time',isi.Time(threshidx));     % update AP
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no threshold.\n',cnt,obj.nrofaps); end
        end  %Djai's version               
        
        function obj = getslowahp(obj)
            % find the after-hyperpolarising potential of APs. 
            % Runs a window of few milliseconds (obj.ahpSnakeWin) from peak of AP forward in 
            % time, until it finds point where the head of this window has a larger value than the tail. Then takes minumum 
            % value within window as AHP.
            % NOTE: Assumes list of ACTIONPOTENTIALS is sorted!
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            snakelen = floor(obj.ahpSnakeWin/mean(diff(obj.Time)));                               % obtain length of ahp window
            %obj = obj.sortaps; % could add this line to be sure...
            for i = 1:obj.nrofaps                                                                   % for every AP event
                if i == obj.nrofaps, 
                    endpoint = obj.TimeInfo.End;                                                    % find end point
                else
                    endpoint = obj.getap(i+1).peak_time;
                end
                
                if ~isempty(obj.getap(i).ahp_time) && ~isnan(obj.getap(i).ahp_time)
                    startpoint = obj.getap(i).ahp_time;
                    isi        = obj.getsampleusingtime(startpoint,endpoint);                          
                    startpoint = isi.Time(find(isi.Data(1:end-snakelen*2+1)>isi.Data(snakelen*2:end),1));  % check if signal goed down again if a fast AHP happened
                else
                    startpoint = obj.getap(i).peak_time;
                end
                if ~isempty(startpoint)                                                          % do not search for slow ahp if the signal does not go down again after fast ahp
                isi        = obj.getsampleusingtime(startpoint,endpoint);                           % get timeseries of data after AP
                ahpwinstrt = find(isi.Data(1:end-snakelen*3+1)<isi.Data(snakelen*3:end),1);             % find start of ahp window
                ahpwinend  = ahpwinstrt + snakelen*3;   
                if ahpwinend > numel(isi.Data)                                                      % get end window
                    ahpwinend = numel(isi.Data);                 
                end
                if ~isempty(ahpwinstrt)
                    ahpidx = find(isi.Data==min(isi.Data((ahpwinstrt:ahpwinend))),1);             % find index of minimum
                    if ahpidx>3
                    obj    = obj.updateap(i,'ahp_slow',isi.Data(ahpidx+1),'ahp_slow_time',isi.Time(ahpidx+1));    % update AP
                    else
                    cnt=cnt+1;
                    end
                else
                    cnt=cnt+1;
                end
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, 
                %fprintf('%d out of %d action potentials: no slow AHP.\n',cnt,obj.nrofaps); 
            end
        end
 
        function obj = getahp(obj)
            % find the after-hyperpolarising potential of APs. 
            % Runs a window of few milliseconds (obj.ahpSnakeWin) from peak of AP forward in 
            % time, until it finds point where the head of this window has a larger value than the tail. Then takes minumum 
            % value within window as AHP.
            % NOTE: Assumes list of ACTIONPOTENTIALS is sorted!
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            snakelen = floor(obj.ahpSnakeWin/mean(diff(obj.Time)));                                 % obtain length of ahp window
            %obj = obj.sortaps; % could add this line to be sure...
            for i = 1:obj.nrofaps                                                                   % for every AP event
                endpoint   = obj.getap(i).peak_time+obj.getap(i).halfwidth+4.5+obj.ahpSnakeWin;     % fast ahp needs to be within 4.5+halfwidth ms after peak
                isi        = obj.getsampleusingtime(obj.getap(i).peak_time,endpoint);               % get timeseries of data after AP
                ahpwinstrt = find(isi.Data(1:end-snakelen+1)<isi.Data(snakelen:end),1);             % find start of ahp window
                ahpwinend  = ahpwinstrt + snakelen;                                                 % get end window
                if ~isempty(ahpwinstrt)
                    ahpidx = find(isi.Data==min(isi.Data((ahpwinstrt-1:ahpwinend-1))),1);             % find index of minimum
                    if isi.Time(ahpidx)-obj.getap(i).peak_time<5
                    obj    = obj.updateap(i,'ahp',isi.Data(ahpidx),'ahp_time',isi.Time(ahpidx));    % update AP (only if AHP is fast enough)
                    else
                        cnt=cnt+1;
                    end
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, 
                %fprintf('%d out of %d action potentials: no fast AHP.\n',cnt,obj.nrofaps); 
            end
        end
        
        function obj = getadp(obj)
            % find the after-depolarising potential of APs. 
            % Runs a window of few milliseconds (obj.adpSnakeWin) from ahp of AP forward in time, until it finds point where 
            % the head of this window has a lower value than the tail. Then takes maximum value within window as ADP. The 
            % length of this window may need considerable tweaking.
            % NOTE: Needs ahp time of aps to calculate. 
            % NOTE: Assumes list of ACTIONPOTENTIALS is sorted!
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            snakelen = floor(obj.adpSnakeWin/mean(diff(obj.Time)));                                         % obtain length of adp window
            %obj = obj.sortaps; % could add this line to be sure...
            for i = 1:obj.nrofaps                                                                           % for every AP event
                if ~isempty(obj.getap(i).ahp_time) && ~isnan(obj.getap(i).ahp_time)
                    if i == obj.nrofaps, 
                         endpoint = obj.TimeInfo.End;                                                       % find end point
                    else endpoint = obj.getap(i+1).peak_time;
                    end
                    isi        = obj.getsampleusingtime(obj.getap(i).ahp_time,endpoint);                    % get timeseries of data after AP
                    adpwinstrt = find(isi.Data(1:end-snakelen)>isi.Data(snakelen:end-1),1);                 % find start of adp window
                    if ~isempty(adpwinstrt) 
                        adpwinend = adpwinstrt + snakelen;                                                  % get end window
                        adpidx    = find(isi.Data==max(isi.Data(adpwinstrt:adpwinend)),1);                  % find index of maximum
                        obj       = obj.updateap(i,'adp',isi.Data(adpidx),'adp_time',isi.Time(adpidx));     % update AP 
                    else
                        cnt=cnt+1;
                    end
                else
                    cnt=cnt+1;
                end
            end
            %if cnt>0, fprintf('%d out of %d action potentials: no ADP.\n',cnt,obj.nrofaps); end
        end
        
        function obj = getrelahp(obj)
            % find the relative AHP of APs, difference between AHP and threshold (i.e. amount of undershoot).
            % NOTE: Needs threshold and ahp of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETAHP
            if ~isscalar(obj), error('Object must be scalar'); end
            
            %fast component:
            cnt=0;
            for i = 1:obj.nrofaps
                if ~isempty(obj.getap(i).ahp) && ~isempty(obj.getap(i).thresh) && ~isnan(obj.getap(i).ahp) && ~isnan(obj.getap(i).thresh)
                    obj = obj.updateap(i,'relahp',obj.getap(i).ahp-obj.getap(i).thresh); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            %if cnt>0, fprintf('%d out of %d action potentials: no relative AHP.\n',cnt,obj.nrofaps); end
            
            %slow component:
            cnt=0;
            for i = 1:obj.nrofaps
                if ~isempty(obj.getap(i).ahp_slow) && ~isempty(obj.getap(i).thresh) && ~isnan(obj.getap(i).ahp_slow) && ~isnan(obj.getap(i).thresh)
                    obj = obj.updateap(i,'relahp_slow',obj.getap(i).ahp_slow-obj.getap(i).thresh); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            %if cnt>0, fprintf('%d out of %d action potentials: no relative AHP.\n',cnt,obj.nrofaps); end
        end
        
        function obj = getreladp(obj)
            % find the relative adp of ACTIONPOTENTIALs, difference between ADP and threshold.
            % NOTE: Needs threshold and adp of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETADP
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i = 1:obj.nrofaps
                if ~isempty(obj.getap(i).adp) && ~isempty(obj.getap(i).thresh) && ~isnan(obj.getap(i).adp) && ~isnan(obj.getap(i).thresh)
                    obj = obj.updateap(i,'reladp',obj.getap(i).ahp-obj.getap(i).thresh); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            %if cnt>0, fprintf('%d out of %d action potentials: no relative ADP.\n',cnt,obj.nrofaps); end
        end   
        
        function obj = getrelamp(obj)
            % find the amplitudes of ACTIONPOTENTIALs, difference between threshold and peak 
            % NOTE: Needs peak, and threshold of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETTHRESH, FINDAPS
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i = 1:obj.nrofaps
                if ~isempty(obj.getap(i).peak) && ~isempty(obj.getap(i).thresh) && ~isnan(obj.getap(i).peak) && ~isnan(obj.getap(i).thresh)
                    obj = obj.updateap(i,'amp',obj.getap(i).peak-obj.getap(i).thresh); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no amplitude.\n',cnt,obj.nrofaps); end
        end
     
        function obj = gethalfwidth(obj)
            % find the halfwidths of ACTIONPOTENTIALs. 
            % Resamples AP to 100KHz to ensure same precision when estimating AP halfwidth, regardless of sampling frequency.
            % NOTE: Needs peak time, threshold and amplitude of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETTHRESH, FINDAPS, GETRELAMP.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i = 1:obj.nrofaps 
                if ~isempty(obj.getap(i).thresh) && ~isempty(obj.getap(i).amp) && ~isnan(obj.getap(i).thresh) && ~isnan(obj.getap(i).amp)
                    hwline  = obj.getap(i).amp*0.5+obj.getap(i).thresh;
                    endtime = obj.Time(find(obj.Data(obj.Time > obj.getap(i).peak_time & obj.Time < obj.TimeInfo.End) < hwline,1)+2) + obj.getap(i).peak_time;
                    if ~isempty(endtime)
                        hw  = obj.getsampleusingtime(obj.getap(i).thresh_time, endtime);
                        hw  = hw.resample(hw.TimeInfo.Start:obj.upsample:hw.TimeInfo.End,'linear');
                        % update AP
                        hwstrt = hw.Time(find(hw.Data>hwline,1,'first'));
                        hwend  = hw.Time(find(hw.Data>hwline,1,'last' ));
                        obj = obj.updateap(i,'halfwidth_strt_time',hwstrt,'halfwidth_end_time',hwend,'halfwidth',hwend-hwstrt);
                    else
                        cnt=cnt+1;
                    end
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no halfwidth.\n',cnt,obj.nrofaps); end
        end
        
        function obj = getmaxdvdt(obj)
            % find the maximum dvdt during rising phase of action potential.
            % NOTE: Needs peak time, and threshold time of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETTHRESH, GETMINDVDT, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i = 1:obj.nrofaps 
                if ~isempty(obj.getap(i).thresh_time) && ~isempty(obj.getap(i).peak_time) && ~isnan(obj.getap(i).thresh_time) && ~isnan(obj.getap(i).peak_time)
                    ts  = obj.getsampleusingtime(obj.getap(i).thresh_time,obj.getap(i).peak_time).getdvdtts;
                    mx  = max(ts);
                    obj = obj.updateap(i,'maxdvdt',mx,'maxdvdt_time',ts.Time(find(ts.Data==mx,1))); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no Maxdvdt.\n',cnt,obj.nrofaps); end
        end
        
        function obj = getmindvdt(obj)
            % find the minimum dvdt during repolarising phase of ACTIONPOTENTIAL.
            % NOTE: Needs peak time, and ahp time of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETTHRESH, GETAHP, GETMAXDVDT, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i = 1:obj.nrofaps 
                if ~isempty(obj.getap(i).ahp_time) && ~isnan(obj.getap(i).ahp_time)
                    ahp_time=obj.getap(i).ahp_time;
                else
                    ahp_time=obj.getap(i).ahp_slow_time;
                end
                if ~isempty(ahp_time) && ~isempty(obj.getap(i).peak_time) && ~isnan(ahp_time) && ~isnan(obj.getap(i).peak_time)
                    ts  = obj.getsampleusingtime(obj.getap(i).peak_time,ahp_time).getdvdtts;
                    mn  = min(ts);
                    obj = obj.updateap(i,'mindvdt',mn,'mindvdt_time',ts.Time(find(ts.Data==mn,1))); % update AP 
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no mindvdt.\n',cnt,obj.nrofaps); end
        end 
        
        function obj = getonsrapidity(obj)
            % find the onset rapidity of an action potential.
            % Resamples AP to 1MHz to ensure same precision when estimating AP halfwidth, regardless of sampling frequency.
            % Of course, the value obtained for onset rapidity still depends hewavily on sampling frequency!
            % NOTE: Needs maxdvdt time, and threshold time of aps to calculate. 
            % See also ACTIONPOTENTIAL, GETTHRESH, GETMAXDVDT.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt = 0;
            for i = 1:obj.nrofaps 
                if ~isempty(obj.getap(i).thresh_time) && ~isempty(obj.getap(i).maxdvdt_time) && ~isnan(obj.getap(i).thresh_time) && ~isnan(obj.getap(i).maxdvdt_time)
                    if obj.getap(i).maxdvdt > obj.apThreshRapidity
                        
                        ts   = obj.getsampleusingtime(obj.getap(i).thresh_time,obj.getap(i).maxdvdt_time);
                        vm   = ts.resample(ts.TimeInfo.Start:obj.upsample:ts.Time(end-1),'linear').getdata;
                        dvts = ts.getdvdtts;
                        dv   = dvts.resample(dvts.TimeInfo.Start:obj.upsample:dvts.TimeInfo.End,'linear').getdata;

                        % find rapidity threshold crossing and make a fitting window.
%                         strtidx = max([1, length(dv) - find(flipud(dv)<obj.apThreshRapidity,1) - obj.onsetrapfitwin + 1]);
%                         endidx  = min([strtidx + obj.onsetrapfitwin + 1, length(vm)]);

                        % temporary hack for IQ paper analysis
                        % startidx at crossing of threshold 15 mV/ms
                        strtidx = max([1, length(dv) - find(flipud(dv)<15,1)+ 1]);
%                         % endidx at 3 mV above the voltage at which
%                         % threshold was crossed
%                         endidx = find(vm > vm(strtidx)+3, 1);                  
                        

                        % linear fit
%                         x = vm(strtidx:endidx);
%                         y = dv(strtidx:endidx);
%                         if numel(vm)>= endidx
%                             x = vm([strtidx, endidx]);
%                             y = dv([strtidx, endidx]);
%                             vmcentre = min(x)+0.5*range(x);
%                             onsetrapfit = polyfit(x, y, 1); % straight fit through points enclosed in window
%                             obj = obj.updateap(i,'onsetrapidity',onsetrapfit(1),'onsetrapfit',onsetrapfit,'onsetrapvm',vmcentre); % update AP 
%                         else
%                             warning('not enough space to fit onsetrapidity window before reaching maxDVDT')
%                             obj = obj.updateap(i,'onsetrapidity',NaN,'onsetrapvm',NaN);
%                         end
                         % instead of linear fit, just look at the maximum
                         % slope of the phase plot from threshold
                        x = vm(strtidx:end);
                        y = dv(strtidx:end);
                        dvvm=diff(dv)./diff(vm);
                        [onsetrap, index]=max(dvvm);
                        if onsetrap > 150
                            onsetrap=NaN;
                        end
                        obj = obj.updateap(i,'onsetrapidity',onsetrap,'onsetrapvm',vm(index));
                    end
                else
                    cnt = cnt+1;
                end
            end
            if cnt>0, 
                fprintf('%d out of %d action potentials: no onset rapidity.\n',cnt,obj.nrofaps); 
            end
        end 
        
        function obj = getapstartend(obj)
            % find start and end time of APs. Mainly used for waveform, so takes few ms before threshold as startpoint
            % (determined by obj.preAPtime), and thresh_time of next AP or obj.postAPtime as end of AP
            % NOTE: Needs ahp time, and threshold time of aps to calculate. 
            % See also ACTIONPOTENTIAL, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            obj = obj.sortaps; % to be sure
            for i = 1:obj.nrofaps
                if ~isempty(obj.getap(i).ahp_time) && ~isnan(obj.getap(i).ahp_time)
                    ahp_time=obj.getap(i).ahp_time;
                else
                    ahp_time=obj.getap(i).ahp_slow_time;
                end
                if ~isempty(obj.getap(i).thresh_time) && ~isnan(obj.getap(i).thresh_time)
                    obj = obj.updateap(i,'start_time',obj.getap(i).thresh_time-obj.preAPtime); % update AP
                else
                    cnt=cnt+1;
                end
                if ~isempty(ahp_time) && ~isnan(ahp_time)
                    if i == obj.nrofaps
                        endtime = min([obj.TimeInfo.End,ahp_time+obj.postAPtime]);
                    else
                        endtime = min([obj.TimeInfo.End,ahp_time+obj.postAPtime,obj.getap(i+1).thresh_time]);
                    end
                    obj = obj.updateap(i,'end_time',endtime); % update AP
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no start and/or end times.\n',cnt,obj.nrofaps); end
        end
        
        function obj = findapnr(obj)
            % number ACTIONPOTENTIALs according to peak time
            % See also ACTIONPOTENTIAL, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            obj = obj.sortaps; % to be sure
            for i=1:obj.nrofaps, obj = obj.updateap(i,'number',i); end
        end
        
        function obj = isifreqaps(obj)
            % add instantaneous firing frequency to ACTIONPOTENTIALs.
            % See also ACTIONPOTENTIAL, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            obj = obj.sortaps; % to be sure
            for i=1:obj.nrofaps
                if i == 1,
                    obj = obj.updateap(i,'freq',0);
                else
                    isi = obj.getap(i).peak_time - obj.getap(i-1).peak_time;
                    obj = obj.updateap(i,'isi',isi,'freq',1e3/isi);
                end
            end
        end
        
        function obj = addapts(obj)
            % add timeseries data to ACTIONPOTENTIAL.
            % NOTE: Needs start and end time of aps to calculate. 
            % See also ACTIONPOTENTIAL, TIMESERIES.
            if ~isscalar(obj), error('Object must be scalar'); end
            cnt=0;
            for i=1:obj.nrofaps
                if ~isempty(obj.getap(i).start_time) && ~isempty(obj.getap(i).end_time) && ~isnan(obj.getap(i).start_time) && ~isnan(obj.getap(i).end_time)
                    tmp = obj.getsampleusingtime(obj.getap(i).start_time,obj.getap(i).end_time);
                    tmp = timeseries(tmp.Data,tmp.Time);
                    tmp.DataInfo.Units = obj.DataInfo.Units;
                    tmp.TimeInfo.Units = obj.TimeInfo.Units;
                    obj = obj.updateap(i,'ts',tmp);
                else
                    cnt=cnt+1;
                end
            end
            if cnt>0, fprintf('%d out of %d action potentials: no Timeseries.\n',cnt,obj.nrofaps); end
        end
        
        function obj = addapevents(obj,varargin)
            % function to add APs as tsdata.event to the timeseries
            % varargin must all be strings, corresponding to an AP feature (e.g. peak, threshold, ahp, ...). For every 
            % argument in, an event will be added for each AP. If no extra inputs are provided, only peak events are added.
            % See also ACTIONPOTENTIAL, TIMESERIES, TSDATA.EVENT.
            if ~isscalar(obj), error('Object must be scalar'); end
            if obj.nrofaps>0,
                if nargin == 1
                    obj = obj.addevent(obj.getap.eventme('peak')); 
                else
                    for i=1:numel(varargin)
                        feature = varargin{i};
                        if ischar(feature)
                            obj = obj.addevent(obj.getap.eventme(feature));
                        else
                            error('Please provide only string input')
                        end
                    end
                end
            end
        end
        
        function isits = getisits(obj)
            % get waveforms of the inter spike intervals. Returns a cell array of isi timeseries objects.
            % NOTE: Needs peak time of aps to calculate. 
            % See also ACTIONPOTENTIAL, TIMESERIES, FINDAPS.
            if ~isscalar(obj), error('Object must be scalar'); end
            if obj.nrofaps>1
                isits = arrayfun(@(x) obj.getsampleusingtime(obj.getap(x-1).peak_time,obj.getap(x).peak_time),2:obj.nrofaps,'UniformOutput',false);
            end
        end
        
        function isiahpthreshts = getisiahpthreshts(obj)
            % get waveforms of the inter spike intervals. Returns a cell array of isi timeseries objects. Assumes ahp_time
            % and thresh_time is available for all APs.
            % NOTE: Needs ahp time and threshold time of aps to calculate. 
            % See also ACTIONPOTENTIAL, TIMESERIES, FINDAPS, GETAHP, GETTHRESH.
            if ~isscalar(obj), error('Object must be scalar'); end
            if obj.nrofaps>1
                isiahpthreshts = arrayfun(@(x) obj.getsampleusingtime(obj.getap(x-1).ahp_time,obj.getap(x).thresh_time),2:obj.nrofaps,'UniformOutput',false);
            end
        end  
        
        function data = plateaufix(obj, data)
            % Peaks with plateaus (i.e. two consecutive sample points have EXACT same value) are not found by findpeaks.m.
            % Therefore, substract very small value from second value in plateau to detect peak. Trick obtained from:
            % http://stackoverflow.com/questions/4376636/unexpected-behaviour-of-function-findpeaks-in-matlabs-signal-processing-toolbox
            % In future, issue of having consecutive sample points having same exact value should be dealt with by recording with appropriately high gain...
            if ~isscalar(obj), error('Object must be scalar'); end
            prevalue = data(1);
            for i = 2:size(data,1)
                if  data(i) == prevalue
                    data(i) =  data(i) - obj.plateaucorrection; 
                else
                    prevalue = data(i);
                end
            end
        end
        
        %% ------------------------------------------ PLOTTING METHODS -------------------------------------------------------        
        function plot(obj,varargin)
            for i=1:numel(obj)
                plot@timeseries(obj(i),varargin{:})
            end
            grid on
        end
        
        function plotaps(obj,varargin)
            % plot all APs found in a trace. 
            % Supports name/value input arguments for plotting. One extra option is to choose 'superimpose' as extra input 
            % argument with peak, thresh, ahp, or adp as corresponding value. In that case, APs are plotted on top of each 
            % other, aligned to choice. 
            % See also ACTIONPOTENTIAL.
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin == 1, obj.getap.plot
            else obj.getap.plot(varargin{:})
            end
        end
        
        function plotapevents(obj,feature,varargin)
            % function to plot AP feature (e.g. peak, threshold, ahp, ...) as tsdata.event on the timeseries. 
            % Feature should be either a string or a cell array of strings.
            % See also ACTIONPOTENTIAL, TIMESERIES, TSDATA.EVENT.
            if isscalar(obj)
                if nargin == 1, obj.addapevents('peak','thresh','ahp', 'ahp_slow').plot(varargin{:})
                elseif ischar(feature)
                    obj.addapevents(feature).plot(varargin{:})
                elseif iscell(feature) && all(cellfun(@ischar,feature))
                    obj = obj.addapevents(feature{:});
                    obj.plot(varargin{:})
                else
                    error('Feature should be either a string or a cell array of strings.')
                end
            else
                if nargin == 1, feature = {'peak','thresh','ahp', 'ahp_slow'}; end
                hold on
                    arrayfun(@(x) obj(x).plotapevents(feature,varargin{:}),1:numel(obj),'UniformOutput',false)
                hold off
            end
        end
        
        function checkapplot(obj)
            % plot TRACE with derivative overlay, plus peak aligned APs and AP phase plots.
            % See also ACTIONPOTENTIAL, TIMESERIES, TSDATA.EVENT, GETDVDTTS.
            if ~isscalar(obj), error('Object must be scalar'); end
            figure
            set(gcf,'position',[100,400,1700,500],'color',[1 1 1])
            yyaxis left;  obj.plotapevents({'peak','ahp','thresh'},'linewidth',2)
            yyaxis right; obj.getdvdtts.plot('r')
        end
        
        function phaseplot(obj,varargin)
            % plot phase plane of TRACE
            % See also TIMESERIES, GETDATA, GETDVDTTS.
            if ~isscalar(obj), error('Object must be scalar'); end
            vm = obj.getdata;
            dv = obj.getdvdtts.getdata;
            plot(vm(1:end-1),dv,varargin{:})
        end
        
        function phaseplotaps(obj,fitflag,lines)
            % plot phase plots of action potentials.
            % plot the onsetrapidity fits as well by setting fitflag to 1, add lines (threshold and zeros) with lines=1.
            % See also ACTIONPOTENTIAL, GETONSETRAPIDITY.
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin == 1, fitflag=0; lines=0;
            elseif nargin ==2, lines=0;
            end           
            obj.getap.phaseplot; 
            if fitflag
                hold on
                xFit = cat(1,[obj.getap.onsetrapvm]-1,[obj.getap.onsetrapvm]+1)'; % make the x-range, 2 mV around point of onset rapidity threshold crossing
                arrayfun(@(x) plot(xFit(x,:),polyval(obj.getap(x).onsetrapfit, xFit(x,:)),'r','LineWidth',2),1:obj.nrofaps,'UniformOutput',false);             % evaluate linear fit
            end
            if lines
                line(xlim,[obj.apThreshRapidity obj.apThreshRapidity],'color','r','linestyle','--','linewidth',2); 
                line(xlim,[0 0],'color',[0.5 0.5 0.5],'linestyle','-','linewidth',2); 
                line([0 0],ylim,'color',[0.5 0.5 0.5],'linestyle','-','linewidth',2); 
                grid on
            end
        end
        
        function phasetunnel(obj,varargin)
            % plot phase 'tunnel' of TRACE. 
            % 3D plot with axes: units, time, units/time. So for current clamp recording: mV x ms x mV/ms.
            if ~isscalar(obj), error('Object must be scalar'); end
            vm = obj.Data;
            t  = obj.Time;
            dv = obj.getdvdtts.getdata;
            plot3(vm(1:end-1),t(1:end-1),dv,varargin{:})
        end
        
        function plotnormisiahpthreshts(obj,varargin)
            % plot isis normalised in time (ahp to next thresh). Just for fun. Supports name/value pairs for formatting plot.
            % See also ACTIONPOTENTIAL, TIMESERIES, GETISIAHPTHRESHTS.
            if ~isscalar(obj), error('Object must be scalar'); end
            isits = obj.getisiahpthreshts;
            hold on
                arrayfun(@(x) isits{x}.set('Time',linspace(0,100,isits{x}.Length)).plot(varargin{:}),1:numel(isits),'UniformOutput',false)
            hold off    
        end     
        
        function plotnormisi(obj,varargin)
            % plot isis normalised in time (peak to peak). 
            % Just for fun. Supports name/value pairs for formatting plot.
            % See also ACTIONPOTENTIAL, TIMESERIES, GETISITS.
            if ~isscalar(obj), error('Object must be scalar'); end
            isits = obj.getisits;
            hold on
                arrayfun(@(x) isits{x}.set('Time',linspace(0,100,isits{x}.Length)).plot(varargin{:}),1:numel(isits),'UniformOutput',false)
            hold off    
        end
        
        function plotisi(obj,varargin)
            % plot isis, superimposed. 
            % Supports name/value pairs for formatting plot.
            % See also ACTIONPOTENTIAL, TIMESERIES, GETISITS.
            if ~isscalar(obj), error('Object must be scalar'); end
            isits = obj.getisits;
            hold on
                arrayfun(@(x) isits{x}.set('Time',isits{x}.Time-isits{x}.TimeInfo.Start).plot(varargin{:}),1:numel(isits),'UniformOutput',false)
            hold off
        end
        
        function stdshade(obj,varargin)
            % plot mean and standard deviation or s.e.m. shading for a 1xN TRACE.
            % intention: apply the stdshade function found online to a set of traces 
            error('not implemented yet')
        end
                
    end  
end

