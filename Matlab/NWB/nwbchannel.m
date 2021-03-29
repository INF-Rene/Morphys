classdef nwbchannel < Sharedmethods
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
        guid_stimset           % Globally unique identifier of parent NWB object.
        filename            % filename of parent NWB
        filedirectory       % directory of parent NWB
        name                % AD name
        number              % AD number
        
        sweeptable
        sweepnrs
        associated_stimsets  % simultaneously recorded stimsets of which sweeps will be saved here as well
        labbooknum
        
        updatephase = 1     % phase of class; (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)
        
        sweeps              % list of sweeps
        stimwaves           % stimwave objects
        
        nrofsweeps  = 0     % number of Analog Input channels recorded
        
    end
    
    %####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = nwbchannel(varargin)
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
        function obj = addsweeps(obj)
            
            if ~isscalar(obj)
                for i=1:numel(obj)
                    obj(i).addsweeps;
                end
            else
                fn=fullfile(obj.filedirectory, obj.filename);
                
                epochtables = table();
                epochtables.stimsets=obj.associated_stimsets;
                for i=1:height(epochtables)
                    epochtables.epochtable(i) = {obj.makeepochtable(obj.associated_stimsets{i})};
                    if any(any(epochtables.epochtable{i}.deltas>0))
                        epochtables.stimwavemode(i)=2; %for every sweep the stimwave will be read from file and saved
                    else
                        %epochtables.stimwavemode(i)=1; %only first stimwave will be saved as the stimwave is assumed to be exactly the same for every sweep                        
                        epochtables.stimwavemode(i)=2; % unsophisticated override of this setting
                    end
                end
                
                TPepochdur=100;
                for i=1:height(obj.sweeptable)
                    %read swp data
                    data=h5read(fn, [obj.sweeptable.dataloc{i} '/data']);
                    
                    % get sweep assocated labbook data
                    sweepLB = obj.labbooknum(obj.labbooknum.SweepNum == obj.sweeptable.sweepnr(i),:);
                    
                    tmp=strcmp(epochtables.stimsets, obj.sweeptable.protocol{i});
                    scepochtable=epochtables.epochtable{tmp};
                    scepochtable.duration(1)=TPepochdur;
                    stimwavemode=epochtables.stimwavemode(tmp);
                    
%                     scalefactor=unique(sweepLB.StimScaleFactor(~isnan(sweepLB.StimScaleFactor)));
%                     if ~isempty(scalefactor) && ~isscalar(scalefactor)
%                         warning('Multiple scalefactors found. Using the last one.');
%                         scalefactor=scalefactor(end);
%                     end
%                     if ~isempty(scalefactor) && scalefactor~=1
%                         scepochtable.firstlevel=scepochtable.firstlevel.*scalefactor;
%                         stimwavemode=2;
%                     end
                    %starttime= h5read(fn, [obj.datalocs{i} '/starting_time']);
%                     if stimwavemode==2 && i>1
%                         stimdata=h5read(fn, [obj.stimdatalocs{i} '/data']);
%                     end
                    
                    starttime = nanmin(sweepLB.TimeStamp);
                    endtime = nanmax(sweepLB.TimeStamp);
                    nrinset = sweepLB.SetSweepCount(~isnan(sweepLB.SetSweepCount))+1;
                    if isempty(nrinset), nrinset=0; end
                    
                    % samplefreq is sometimes sweep specific
                    units=h5readatt(fn, [obj.sweeptable.dataloc{i} '/data'], 'IGORWaveUnits');
                    scaling=h5readatt(fn, [obj.sweeptable.dataloc{i} '/data'], 'IGORWaveScaling');
                    sampleint=scaling(3);
                    if strcmp(units{2}, 'ms'), sampleint=sampleint./1000;end
                    samplefreq=1/sampleint;
                    
                    sweep2add = Sweep('number', obj.sweeptable.sweepnr(i), 'nrinset', nrinset, 'samplefreq', samplefreq,...
                        'datetimestart', starttime,'datetimeend', endtime, 'labbooknum', sweepLB,...
                        'stimwavename', obj.sweeptable.protocol{i}, 'stimwavemode', stimwavemode,...
                        'stimdataloc', obj.sweeptable.stimdataloc{i},'guid_AD', obj.guid);
                    
                    % remove filling zeroes at end of aborted sweep
                    removesmpl=[];
                    if data(end)==0
                        removesmpl=find(data~=0,1, 'last')+1;
                        data(removesmpl:end)=[];
%                         stimdata(removesmpl:end)=[];
                    end
                    
%                     if stimwavemode==2 && i>1
%                         stimtime=0:1e3/obj.samplefreq:1e3/obj.samplefreq*(numel(stimdata)-1);
%                         obj.stimwaves{i}=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/obj.samplefreq*(numel(stimdata)-1));
%                     end
                    
                    
                    sweep2add=sweep2add.adddata(data, units{1});
                    sweep2add=sweep2add.addtime(samplefreq);
                    sweep2add.timespan=sweep2add.Time(end);
                    
                    %find TP epoch duration from stimddata. This is really stupid but has to be done like this since the LB entries for
                    %TP duration are always empty... Hopefully will be implemented in MIES at a later time.
                    if sweep2add.timespan<100
                        stimdata=h5read(fn, [obj.sweeptable.stimdataloc{i} '/data'], 1, sweep2add.Length/10, 10); % first 100 ms (take only 1 in 10 samples for efficiency)
                    else
                        stimdata=h5read(fn, [obj.sweeptable.stimdataloc{i} '/data'], 1, samplefreq/10*100e-3, 10); % first 100 ms (take only 1 in 10 samples for efficiency)
                    end
                    TPepochdur=(find(stimdata~=0,1, 'first')-1)*10/samplefreq*1e3/0.45; % div by 0.45 bc the delay of the TP is 45% of the TP epoch (since the TP is 10% and in the middle of the epoch)
                    if isempty(TPepochdur) || round(sweep2add.Time(end)) == sum(scepochtable.duration(2:end))
                        scepochtable(1,:)=[]; %no TP inserted
                        TPepochdur=0;
                    else 
                        scepochtable.duration(1)=TPepochdur;
                    end
                    if nnz(scepochtable.type==7)==1 %type 7 is "Loaded custom wave", duration is always 0 in epochtable
                        scepochtable.duration(scepochtable.type==7)=sweep2add.Time(end)-sum(scepochtable.duration);
                    end
                    %get epoch information to split in epochs
                    if sweep2add.Time(end)>sum(scepochtable.duration) && sweep2add.Time(end)-sum(scepochtable.duration)>eps(sum(scepochtable.duration))
                        % sometimes duration of the recording is longer than the protocol
                        % because another protocol was simultaneously recorded
                        % in that case add another epoch to fill the sweep
                        dur=sweep2add.Time(end)-sum(scepochtable.duration);
                        scepochtable2=scepochtable;
                        scepochtable2(end+1,{'idx', 'number', 'duration', 'typestr', 'idxstr',...
                            'type', 'firstlevel', 'pulsewidth', 'maxfrequency', 'pulseperiod', 'deltas', 'stepdiff'})=...
                            {scepochtable.idx(end)+1, scepochtable.number(end)+1, dur, 'Square pulse',...
                            native2unicode(65+scepochtable.number(end)+1),...
                            0,0,0,0,0,[0, 0, 0, 0], 0};
                        sweep2add=sweep2add.addNWBepochs(scepochtable2,fn);
                    else
                        sweep2add=sweep2add.addNWBepochs(scepochtable,fn);
                    end
                    
                    
                    if isempty(obj.sweeps)
                        obj.sweeps        = sweep2add;
                    else
                        obj.sweeps(end+1) = sweep2add;
                    end
                end
            end
        end
        function obj = loadsweepdata(obj)
            
            if ~isscalar(obj)
                for i=1:numel(obj)
                    obj(i).loadsweepdata;
                end
            else
                %read epochtable and stimwave
                fn=fullfile(obj.filedirectory, obj.filename);
                
                obj.stimunits=h5readatt(fn, [obj.sweeptable.stimdataloc{1} '/data'], 'IGORWaveUnits');
                obj.units=h5readatt(fn, [obj.sweeptable.dataloc{1} '/data'], 'IGORWaveUnits');
                
                scaling=h5readatt(fn, [obj.sweeptable.dataloc{1} '/data'], 'IGORWaveScaling');
                sampleint=scaling(3);
                if strcmp(obj.units{2}, 'ms'), sampleint=sampleint./1000;end
                obj.samplefreq=1/sampleint;
                
                obj.stimwaves=cell(numel(obj.associated_stimsets),1);
                epochtable=cell(numel(obj.associated_stimsets),1);
                stimwavemode=[];
                for i=1:numel(obj.associated_stimsets)
                    loc = find(strcmp(obj.sweeptable.protocol, obj.associated_stimsets{i}),1,'first');
                    stimdata=h5read(fn, [obj.sweeptable.stimdataloc{loc} '/data']);
                    stimtime=0:1e3/obj.samplefreq:1e3/obj.samplefreq*(numel(stimdata)-1);
                    obj.stimwaves{i}=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/obj.samplefreq*(numel(stimdata)-1));
                    epochtable{i} = obj.makeepochtable(obj.associated_stimsets{i});
                    if any(any(epochtable{i}.deltas>0))
                        stimwavemode(i)=2; %for every sweep the stimwave will be read from file and saved
                    else
                        stimwavemode(i)=1; %only first stimwave will be saved as the stimwave is assumed to be exactly the same for every sweep
                    end
                end
                
                
                %read sweep data and correct stimwaves for scalefactors
                
                for i=1:numel(obj.sweeptable.dataloc)
                    data=h5read(fn, [obj.sweeptable.dataloc{i} '/data']);
                    
                    % get sweep assocated labbook data
                    sweepLB = obj.labbooknum(:, :, obj.labbooknum(1,strcmp(obj.labbooknum_keys, 'SweepNum'),:) == obj.sweepnrs(i));
                    scepochtable=epochtable;
                    scalefactor=unique(sweepLB.StimScaleFactor(~isnan(sweepLB.StimScaleFactor)));
                    if ~isempty(scalefactor) && ~isscalar(scalefactor)
                        warning('Multiple scalefactors found. Using the last one.');
                        scalefactor=scalefactor(end);
                    end
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
                    if isempty(nrinset), nrinset=0; end
                    % samplefreq is sometimes sweep specific
                    scaling=h5readatt(fn, [obj.datalocs{i} '/data'], 'IGORWaveScaling');
                    sampleint=scaling(3);
                    if strcmp(obj.units{2}, 'ms'), sampleint=sampleint./1000;end
                    samplefreq=1/sampleint;
                    if samplefreq~=obj.samplefreq, obj.samplefreq=NaN; end
                    
                    sweep2add = Sweep('number', obj.sweepnrs(i), 'nrinset', nrinset, 'samplefreq', samplefreq,...
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
                    sweep2add=sweep2add.addtime(samplefreq);
                    
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
        
        function epochtable = makeepochtable(obj, protocol)
            fn=fullfile(obj.filedirectory, obj.filename);
            typedata=h5read(fn, ['/general/stimsets/' protocol '_SegWvType']);
            typelabels=h5readatt(fn, ['/general/stimsets/' protocol '_SegWvType'], 'IGORWaveDimensionLabels');
            nrofepochs=typedata(strcmp('Total number of epochs', typelabels(2:end)));
            epochtable=table([1:nrofepochs]',[0:nrofepochs-1]',typedata(1:nrofepochs), ...
                'VariableNames', {'idx', 'number', 'type'});
            
            epochdata=h5read(fn, ['/general/stimsets/' protocol '_WP']);
            labels=h5readatt(fn, ['/general/stimsets/' protocol '_WP'], 'IGORWaveDimensionLabels');
            
            untypes= unique(epochtable.type);
            for i = 1:numel(untypes)
                locs=epochtable.type==untypes(i);
                epochtable.duration(locs)=epochdata(untypes(i)+1,locs, strcmp(labels(1,2:end)', 'Duration'))';
                epochtable.firstlevel(locs)=epochdata(untypes(i)+1,locs, strcmp(labels(1,2:end), 'Amplitude'))';
                epochtable.pulsewidth(locs)=epochdata(untypes(i)+1,locs, strcmp(labels(1,2:end), 'Train pulse duration'))';
                epochtable.maxfrequency(locs)=epochdata(untypes(i)+1,locs, strcmp(labels(1,2:end), 'Sin/chirp/saw frequency'))';
                epochtable.pulseperiod(locs)=1000./epochtable.maxfrequency(locs);
            end
            epochtable.pulseperiod(epochtable.pulseperiod==Inf)=0;
            types={'Square pulse','Ramp','Noise','Sin','Saw tooth','Pulse train','PSC','Load custom wave','Combine'};
            epochtable.typestr=types(epochtable.type+1)';
            % for some reason the pulsetrains have a duration of half a period longer in the log than in reality:
            epochtable.duration(epochtable.type==5) = epochtable.duration(epochtable.type==5) - 0.5*epochtable.pulseperiod(epochtable.type==5) ...
                + epochtable.pulsewidth(epochtable.type==5);
            
            % check if there is any sweep delta values
            deltas={'Duration delta','Amplitude delta','Sin/chirp/saw frequency delta','Train pulse duration delta'};
            deltadata=squeeze(epochdata(1,epochtable.idx, ismember(labels(1,2:end), deltas(1:2))));
            if size(deltadata,2)==1, deltadata=deltadata';end % in case of only 1 epoch the squeeze function makes the output in different direction
            deltadata(:,3:4)=squeeze(epochdata(1,epochtable.idx, ismember(labels(1,2:end), deltas(3:4))));
            epochtable.deltas=deltadata;
            % use this data later on the sweep level to correct the values for delta and scalefactor
            
            
            if contains(obj.name, 'Endur')
                for i=1:height(epochtable)
                    if epochtable.firstlevel(i)>0
                        midepochtime=sum(epochtable.duration(1:i-1))+0.5*epochtable.duration(i);
                        epochtable.firstlevel(i)=obj.stimwaves{1}.getsampleusingtime(midepochtime-0.1,midepochtime+0.1).mean();
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
            %now add testpulse initial epoch as first row
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
        
        function obj = addstimwaves(obj)
            fn=fullfile(obj.filedirectory, obj.filename);
            stimwavenms={obj.getsweep.stimwavename};
            stimwavemode=[obj.getsweep.stimwavemode];
            stimdatalocs={obj.getsweep.stimdataloc};
            samplefreqs=[obj.getsweep.samplefreq];
            sweepnrs2=[obj.getsweep.number];
            lengths=[obj.getsweep.Length];
            
            for j=1:numel(obj.associated_stimsets)
                loc=find(strcmp(stimwavenms, obj.associated_stimsets{j}));
                if any(stimwavemode(loc)==2) % if already 
                    for i=1:numel(loc)
                        stimdata = h5read(fn, [stimdatalocs{loc(i)} '/data']);
                        if numel(stimdata) > lengths(loc(i))
                            stimdata = stimdata(1:lengths(loc(i)));
                        end
                        units=h5readatt(fn, [stimdatalocs{loc(i)} '/data'], 'IGORWaveUnits');
                        stimtime=0:1e3/samplefreqs(loc(i)):1e3/samplefreqs(loc(i))*(numel(stimdata)-1);
                        stimwavets=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/samplefreqs(loc(i))*(numel(stimdata)-1));
                        stimwavets.DataInfo.Units=units{1};
                        stimwavets.TimeInfo.Units=units{2};
                        stimwave = Stimwave('guid_nwbchannel', obj.guid, 'filename', obj.filename, ...
                            'filedirectory', obj.filedirectory, 'name', stimwavenms{loc(i)},...
                            'sweepnrs', sweepnrs2(loc(i)), 'dataloc', stimdatalocs{loc(i)}, 'units', units );
                        stimwave=stimwave.addts(stimwavets);
                        obj = obj.addstimwave(stimwave);
                    end
                else
                    stimdata = h5read(fn, [stimdatalocs{loc(1)} '/data']);
                    units=h5readatt(fn, [stimdatalocs{loc(1)} '/data'], 'IGORWaveUnits');
                    if all(numel(stimdata) > lengths(loc))
                            stimdata = stimdata(1:lengths(loc(1)));
                    end
                    stimtime=0:1e3/samplefreqs(loc(1)):1e3/samplefreqs(loc(1))*(numel(stimdata)-1);
                    stimwavets=resample(timeseries(stimdata,stimtime), 0:1e3/obj.stimwavesampling:1e3/samplefreqs(loc(1))*(numel(stimdata)-1));
                    stimwavets.DataInfo.Units=units{1};
                    stimwavets.TimeInfo.Units=units{2};
                    stimwave = Stimwave('guid_nwbchannel', obj.guid, 'filename', obj.filename, ...
                            'filedirectory', obj.filedirectory, 'name', stimwavenms{loc(1)},...
                           'sweepnrs', sweepnrs2(loc), 'dataloc', {stimdatalocs{loc}}, 'units', units );
                    stimwave=stimwave.addts(stimwavets);
                    obj = obj.addstimwave(stimwave);
                end
            end
            
        end
        
        function obj = addstimwave(obj, stimwave)
            if isempty(obj.stimwaves)
                obj.stimwaves        = stimwave;
            else
                obj.stimwaves(end+1) = stimwave;
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
        
        function ap = getstimwave(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'stimwaves',0,varargin{:});
        end
        
        function obj = analysenwbchannel(obj)
            % perform default analysis of analog in. Only for primaries with voltage traces are analysed.
            for i = 1:numel(obj)
                obj(i).sweeps = obj(i).getsweep.analysesweep;
                obj(i).updatephase = 3;
            end
        end
        
        
        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function plot(obj,varargin)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
             for i=1:numel(obj)
                ax_handles   = zeros(2,1);
                ax_handles(1) = subplot(2,1,1);
                hold on
                title(obj(i).name)
                for ii=1:obj(i).nrofsweeps
                    obj(i).getsweep(ii).plot(varargin{:})
                end
                hold off
                ax_handles(2) = subplot(2,1,2); 
                hold on
                title(['DA' num2str(obj(i).number-1)])
                for ii=1:numel(obj(i).stimwaves)
                    obj(i).getstimwave(ii).plot(varargin{:})
                end
                hold off
                linkaxes(ax_handles,'x')
             end
             
        end
        
        function plotanalysis(obj)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
            for i=1:numel(obj)
                ax_handles   = zeros(2,1);
                ax_handles(1) = subplot(2,1,1);
                hold on
                title(obj(i).name)
                for ii=1:obj(i).nrofsweeps
                    obj(i).getsweep(ii).plotanalysis
                end
                hold off
                ax_handles(2) = subplot(2,1,2); 
                hold on
                title(['DA' num2str(obj(i).number-1)])
                for ii=1:numel(obj(i).stimwaves)
                    obj(i).getstimwave(ii).plot
                end
                hold off
                linkaxes(ax_handles,'x')
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

