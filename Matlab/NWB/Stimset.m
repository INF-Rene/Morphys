classdef Stimset < Sharedmethods
    %Stimset object
    %   A stimset object is a group of acquisition sweeps in an NWBfile that were acquired by the same protocol/stimset
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       René Wilbers (renewilbers@gmail.com)
    %   Created:      29-03-2019     
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ########################################################## 
    properties
        guid_NWB            % Globally unique identifier of parent NWB object.
        filename            % filename of parent NWB
        filedirectory       % directory of parent NWB
        name                % Stimset name
        datetimestart
        datetimeend
        
        sweepnrs
        labbooknum
        datalocs
        stimdatalocs
        
        updatephase = 1     % phase of class; (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)   
        
        sweeps              % list of Analogins  (= Analogin  class instances)
        units
        stimwaves
        stimunits
        wavename
        samplefreq
        
        nrofsweeps  = 0     % number of Analog Input channels recorded

    end
    
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Stimset(varargin)
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
        % -------------------------------------------- OTHER METHODS --------------------------------------------------------
        function obj = loadsweepdata(obj)
            
            if ~isscalar(obj)
                for i=1:numel(obj)
                    obj(i).loadsweepdata;
                end
            else
                %read epochtable and stimwave
                fn=fullfile(obj.filedirectory, obj.filename);
                
                obj.stimunits=h5readatt(fn, [obj.stimdatalocs{1} '/data'], 'IGORWaveUnits');
                obj.units=h5readatt(fn, [obj.datalocs{1} '/data'], 'IGORWaveUnits');
                scaling=h5readatt(fn, [obj.datalocs{1} '/data'], 'IGORWaveScaling');
                
                sampleint=scaling(3);
                if strcmp(obj.units{2}, 'ms'), sampleint=sampleint./1000;end
                obj.samplefreq=1/sampleint;
                stimdata=h5read(fn, [obj.stimdatalocs{1} '/data']);
                stimtime=0:1e3/obj.samplefreq:1e3/obj.samplefreq*(numel(stimdata)-1);
                obj.stimwaves{1}=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/obj.samplefreq*(numel(stimdata)-1));
                
                epochtable = obj.makeepochtable;
                %read scalefactors belonging to sweeps
                if any(any(epochtable.deltas>0))
                    stimwavemode=2; %for every sweep the stimwave will be read from file and saved
                else
                    stimwavemode=1; %only first stimwave will be saved as the stimwave is assumed to be exactly the same for every sweep
                end
                
                
                for i=1:numel(obj.datalocs)
                    data=h5read(fn, [obj.datalocs{i} '/data']);
                    sweepLB = obj.labbooknum(obj.labbooknum.SweepNum == obj.sweepnrs(i),:);
                    scepochtable=epochtable;
                    scalefactor=unique(sweepLB.StimScaleFactor(~isnan(sweepLB.StimScaleFactor)));
                    if ~isempty(scalefactor) && ~isscalar(scalefactor), error('Multiple scalefactors found'); end
                    if ~isempty(scalefactor) && scalefactor~=1
                        scepochtable.firstlevel=epochtable.firstlevel.*scalefactor;
                        stimwavemode=2;
                    end
                    %starttime= h5read(fn, [obj.datalocs{i} '/starting_time']);
                    if stimwavemode==2 && i>1
                        stimdata=h5read(fn, [obj.stimdatalocs{i} '/data']);
                    end

                    starttime = nanmin(sweepLB.TimeStamp);
                    endtime = nanmax(sweepLB.TimeStamp);
                    nrinset = sweepLB.SetSweepCount(~isnan(sweepLB.SetSweepCount))+1;
                    sweep2add = Sweep('number', obj.sweepnrs(i), 'nrinset', nrinset, 'samplefreq', obj.samplefreq,...
                     'datetimestart', starttime,'datetimeend', endtime, ...
                    'labbooknum', sweepLB,'guid_stimset', obj.guid);
                
                    % remove fake data at end
                    removesmpl=[];
                    if data(end)==0
                        removesmpl=find(data~=0,1, 'last')+1;
                        data(removesmpl:end)=[];
                        stimdata(removesmpl:end)=[];                    
                    end
                    
                    if stimwavemode==2 && i>1
                        stimtime=0:1e3/obj.samplefreq:1e3/obj.samplefreq*(numel(stimdata)-1);
                        obj.stimwaves{i}=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/obj.samplefreq*(numel(stimdata)-1));              
                    end
                        
                
                    sweep2add=sweep2add.adddata(data, obj.units{1});
                    sweep2add=sweep2add.addtime(obj.samplefreq);
                    
                    %get epoch information to split in epochs
                    
                    sweep2add=sweep2add.addNWBepochs(epochtable);

                    if isempty(obj.sweeps) 
                        obj.sweeps        = sweep2add;
                    else
                        obj.sweeps(end+1) = sweep2add;
                    end
    %             obj = obj.updatestimstats;
                end
            end

        end        
  
        function epochtable = makeepochtable(obj)
            fn=fullfile(obj.filedirectory, obj.filename);
            typedata=h5read(fn, ['/general/stimsets/' obj.name{1} '_SegWvType']);
            typelabels=h5readatt(fn, ['/general/stimsets/' obj.name{1} '_SegWvType'], 'IGORWaveDimensionLabels');
            nrofepochs=typedata(strcmp('Total number of epochs', typelabels(2:end)));
            epochtable=table([1:nrofepochs]',[0:nrofepochs-1]',typedata(1:nrofepochs), ...
                'VariableNames', {'idx', 'number', 'type'});
            
            epochdata=h5read(fn, ['/general/stimsets/' obj.name{1} '_WP']);
            labels=h5readatt(fn, ['/general/stimsets/' obj.name{1} '_WP'], 'IGORWaveDimensionLabels');
            
            epochtable.duration=epochdata(1,epochtable.idx, strcmp(labels(1,2:end)', 'Duration'))';
            epochtable.firstlevel=epochdata(1,epochtable.idx, strcmp(labels(1,2:end), 'Amplitude'))';
            epochtable.pulsewidth=epochdata(6,epochtable.idx, strcmp(labels(1,2:end), 'Train pulse duration'))';
            epochtable.maxfrequency=epochdata(6,epochtable.idx, strcmp(labels(1,2:end), 'Sin/chirp/saw frequency'))';
            epochtable.pulseperiod=1000./epochtable.maxfrequency;
            epochtable.pulseperiod(epochtable.pulseperiod==Inf)=0;
            
            types={'Square pulse','Ramp','Noise','Sin','Saw tooth','Pulse train','PSC','Load custom wave','Combine'};
            epochtable.typestr=types(epochtable.type+1)';
            
            % check if there is any sweep delta values
            deltas={'Duration delta','Amplitude delta','Sin/chirp/saw frequency delta','Train pulse duration delta'};
            deltadata=squeeze(epochdata(1,epochtable.idx, ismember(labels(1,2:end), deltas(1:2))));
            deltadata(:,3:4)=squeeze(epochdata(1,epochtable.idx, ismember(labels(1,2:end), deltas(3:4))));
            epochtable.deltas=deltadata;
            % use this data later on the sweep level to correct the values for delta and scalefactor
            
            
            if contains(obj.name, 'Endur')
                for i=1:height(epochtable)
                    if epochtable.firstlevel(i)>0
                        midepochtime=sum(epochtable.duration(1:i-1))+0.5*epochtable.duration(i);
                        epochtable.firstlevel(i)=obj.stimwaves{1}.getsampleusingtime(midepochtime).Data(1);
                    end
                end
            end
                
            epochtable.stepdiff=[NaN; diff(epochtable.firstlevel)];
            for i=1:height(epochtable)
                if epochtable.idx(i)<=26
                    epochtable.idxstr(i)={native2unicode(65+epochtable.number(i))};
                else %in case of >26 epochs use double letter codes, supports up to 676 epochs per sweep
                    epochtable.idxstr(i)={[native2unicode(65+floor(epochtable.number(i)/26)),...
                        native2unicode(65+mod(epochtable.number(i),26))]};
                end
            end
            %now add 100 ms initial epoch as first row
            epochtable(2:end+1,:)=epochtable(1:end,:);
            epochtable(1,:)={0,-1,-1,100,0,0,0,0, 'MIES testpulse', [0 0 0 0], NaN, '<'};
                
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
        function ap = getsweep(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'sweeps',0,varargin{:});
        end
        
        function obj = analysestimset(obj)
            % perform default analysis of analog in. Only for primaries with voltage traces are analysed. 
            for i = 1:numel(obj)
                fprintf('\t\tanalysing stimset %s\n',obj(i).wavename)
                if strcmpi(obj(i).units{1},'mV')
                    obj(i).sweeps = obj(i).getsweep.analysesweep;
                end
                obj(i).updatephase = 3;
            end
        end
        

        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function plot(obj,varargin)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
            for i=1:numel(obj)
                figure
                set(gcf,'position',[413,177,1120,801],'color',[1 1 1])
                ax_handles = zeros(2,1);
                ax_handles(1) = subplot('Position', [0.1 0.3 0.85 0.65]);
                obj(i).getsweep.plot(varargin{:});
                title(strrep(obj.name{1},'_',' '))
                ax_handles(2) = subplot('Position', [0.1 0.05 0.85 0.18]);
                if numel(obj(i).stimwaves)==1
                    obj(i).stimwaves{1}.plot(varargin{:})
                    grid on
                else
                    hold on
                    for j=1:numel(obj(i).stimwaves)
                        obj(i).stimwaves{j}.plot(varargin{:})                        
                        grid on
                    end
                    hold off
                    ylabel(obj(i).stimunits{1});
                end
                xlim([xlim(ax_handles(1))]);
                linkaxes(ax_handles,'x')
%                 set(gcf,'CurrentAxes',ax_handles(1))
            end

        end       
        
        function plotanalysis(obj)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
            for i=1:numel(obj)
                figure
                set(gcf,'position',[413,177,1120,801],'color',[1 1 1])
                ax_handles = zeros(2,1);
                ax_handles(1) = subplot('Position', [0.1 0.3 0.85 0.65]);
                obj(i).getsweep.plotanalysis;
                title(strrep(obj.name,'_',' '))
                ax_handles(2) = subplot('Position', [0.1 0.05 0.85 0.18]);
                if numel(obj(i).stimwaves)==1
                    plot(obj(i).getsweep(1).Time, obj(i).stimwaves{1})
                    grid on
                else
                    hold on
                    for j=1:numel(obj(i).stimwaves)
                        obj(i).stimwaves{j}.plot
                        grid on
                    end
                    hold off
                    ylabel(obj(i).stimunits{1});
                end
                xlim([xlim(ax_handles(1))]);
                linkaxes(ax_handles,'x')
%                 set(gcf,'CurrentAxes',ax_handles(1))
            end

        end  
% 
%         function hfig = plotanalysis(obj)
%             % Plots the Channel object with Epoch boundaries and results of basic analysis. Returns figure handle.
%             % --------------------
%             if isscalar(obj)
%                 ax_handles = zeros(obj.nrofanalogins,1);
%                 for i=1:obj.nrofanalogins
%                     ax_handles(i) = subplot(obj.nrofanalogins,1,i); 
%                     analogin2plot = obj.getin(i);
%                     analogin2plot.plotanalysis
%                     if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
%                         configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.number);
%                     else
%                         configstring = '';
%                     end
%                     title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
%                     ylabel(analogin2plot.units)
%                 end
%                 xlabel('milliseconds')
%                 linkaxes(ax_handles,'x')
%                 hfig = gcf;
%                 set(hfig, 'color', [1 1 1])
%             else
%                 arrayfun(@(x) obj(x).plotanalysis,1:numel(obj),'UniformOutput',false)
%             end
%         end
    end
end

