%% Analysis script
% Written by D.B. Heyer
close all, clear all

%% Set path to load and save data
basedir = '/Users/elinemertens/Matlab Data/H20.29.173/2020_03_03/cell1' ; % specify path to analyzed abf objects. (.mat files)
savedir = '/Users/elinemertens/Matlab Data/Analyzed' ; % specify path for saving the data summary
savename = 'CellSummary_Par1'; % specify filename to be saved

%%
p   = ('/Users/elinemertens/Data/ephys/Hippocampus/H17.29.117.21/H17.29.117.21.IV');
fn  = '2017_03_29_0295.abf';
fp  = fullfile(p,fn);
ss  = load('/Users/elinemertens/Documents/CNCR/Data/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf = Abffile(fp,ss);
%% import CSV files %if you're loading a lot of files at the same time. not used in this version of the script
%requires specific folder- and filenames in location basedir!
% abfs = readtable(fullfile(basedir, 'Abffiles','Abffiles.txt'),'Delimiter',',') ;
% channels = readtable(fullfile(basedir, 'Channels','Channels.txt'),'Delimiter',',') ;
% ins = readtable(fullfile(basedir, 'Analogins','Analogins.txt'),'Delimiter',',') ;
% outs = readtable(fullfile(basedir, 'Analogouts','Analogout.txt'),'Delimiter',',') ;
% sweeps = readtable(fullfile(basedir, 'Sweeps','Sweeps.txt'),'Delimiter',',') ;
% epochs = readtable(fullfile(basedir, 'Epochs','Epochs.txt'),'Delimiter',',') ;
% aps = readtable(fullfile(basedir, 'Actionpotentials','Actionpotentials.txt'),'Delimiter',',') ;

% for data Thijs:
%load(fullfile(basedir,'selectedtables.mat'))
%% Locate files
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};
    
    %% from epoch
    
        % set some defaults
            obj.created = datestr(now());                           % date of object creation
            obj.guid    = upper(char(java.util.UUID.randomUUID())); % assign globally unique identifier
                      
        
        % ------------------------------------------- SHARED METHODS --------------------------------------------------------
    
 
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
      
        
       
            % Function assigns 'value' to any property matching string 'PropertyName' from object.
            % Supports Property/value pair input arguments.
            if nargin < 3 || mod(nargin-1,2)~=0, error('Property/value pairs must come in even number.'); end
            prop2val = cat(1,varargin(1:2:end),varargin(2:2:end))';
            for i=1:size(prop2val,1)
                propertyname = prop2val{i,1};
                value        = prop2val{i,2};
                for ii = 1:numel(obj), obj(ii).(propertyname) = value; end
            end
        
        
        %%
             % call to superclass constructor
           % obj = obj@Sharedmethods;
            %obj = obj@Trace;
      
            % epoch specification
            obj.type         = epochtab.type;  
            obj.typestr      = epochtab.typestr;
            obj.amplitude    = epochtab.firstlevel;
            obj.maxfrequency = epochtab.maxfrequency;
            obj.pulsewidth   = epochtab.pulsewidth; %RWS
            obj.pulseperiod  = epochtab.pulseperiod;   %RWS
            obj.jitter       = nanstd(obj.Data); %DBH
            obj.meansignal   = nanmean(obj.Data); %DBH
        
            
         
    
    %% extract data tables from obj
    abfs = struct2table(obj.getabf.metadata) ;
    channels = struct2table(obj.getabf.getchannel.metadata) ;
    ins = struct2table(obj.getabf.getchannel.getin.metadata) ;
    outs = struct2table(obj.getabf.getchannel.getout.metadata) ;
    sweeps = struct2table(obj.getabf.getchannel.getin.getsweep.metadata) ;
    epochs = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.metadata) ;
    aps = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.getap.metadata) ;
    
    % get instantaneous freq bins
    edges = 1:10:201 ;
    aps.freqbin = discretize(aps.freq, edges) ;
    aps.freqbin(isnan(aps.freqbin)) = 0 ;
    aps.currinj = aps.number*0 ;
    for ii = 1:height(aps)
        aps(ii,:).currinj = epochs(ismember(epochs.guid,aps(ii,:).parent_guid),:).amplitude ;  
    end
    aps.updownratio = aps.maxdvdt./abs(aps.mindvdt) ;
    aps.onsetrapidity(aps.onsetrapidity > 100) = NaN ;
    
    %% Make subset of data per abf file
    abf = SubsetTable2struct(abfs,channels,ins,outs,sweeps,epochs,aps) ;
    %% If abf is a stepprotocol: continue with analysis
    %[ccabf, chs] = isccstep(abf) ; %not necessary if you're not bulk-loading your abffiles
    %if ccabf == 1 
        fprintf('Retrieving analysis parameters from CC-step file %1.0f \n', index);
        %% Analyze
        %sweep = abf.channel([abf.channel.number] == chs).in(strcmp({abf.channel([abf.channel.number] == chs).in.signal},'primary')).sweep ;
        sweep = abf(end).getchannel(1).getin('signal', 'primary').getsweep 
        NrofSweeps = length(sweep) ;  
        
        % find current injection epoch and assign aps to sweep
        for step = 1:length(sweep(1).epoch(1))
            if sweep(1).epoch(step).stepdiff ~= 0 && (sweep(1).epoch(step).stepdiff + sweep(1).epoch(step+1).stepdiff) == 0 && seconds(sweep(1).epoch(step).timespan) > 0.03
                break
            end
        end
        for j = 1:NrofSweeps
            sweep(j).vmbase = sweep(j).epoch(step-1).steadystate ;
            sweep(j).jitter = sweep(j).epoch(step-1).jitter ;
            sweep(j).currinj = sweep(j).epoch(step).stepdiff ;
            sweep(j).vmresponse = sweep(j).epoch(step).vstep ;
            sweep(j).ap = sweep(j).epoch(step).ap ;
            for ap = 1:length(sweep(j).epoch(step+1).ap)
                %if sweep(j).epoch(step+1).ap(ap).start_time > (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+7 && sweep(j).epoch(step+1).ap(ap).start_time < (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+300
                    sweep(j).rbap(ap) = sweep(j).epoch(step+1).ap(ap) ;
                %end
            end
        end
  
        % find rheobase sweep
        apcrit=0;
        for frstspikeswp = 1:NrofSweeps
            if sweep(frstspikeswp).epoch(step).nrofaps > 0 && sweep(frstspikeswp).epoch(step).stepdiff > 0
                apcrit=1;
                break
            end
        end
        
        idx1 = 1 ;
        for j = frstspikeswp:NrofSweeps           
            for ii = 1:length(sweep(j).epoch(step).ap)
                apguids(idx1) = {sweep(j).epoch(step).ap(ii).guid} ;
                idx1 = idx1 + 1 ;
            end
        end             
        aps2 = aps(ismember(aps.guid,apguids),:) ;

        %if apcrit==1 %end at line 332
        % calculate variables
        vmbase = [sweep.vmbase] ;
        Freqs=[];
        StimInts=[];
        currInjections_R=[];
        voltageResponses=[];
        taus=[];
        for j = 1:NrofSweeps           
            if sweep(j,1).currinj >= -100 && sweep(j,1).currinj < 0
                voltageResponses(j,1) = sweep(j,1).vmresponse ; 
                currInjections_R(j,1) = sweep(j,1).currinj ;
                if sweep(j,1).epoch(step).tau < 100 && sweep(j,1).epoch(step).tau > 0 && sweep(j,1).epoch(step).gof > 0.95
                    taus(j,1) = sweep(j,1).epoch(step).tau ;
                else
                    taus(j,1) = NaN;
                end
            end       

            if ~isempty(sweep(j,1).vmresponse) 
                MinVmResponse(j,1) = sweep(j,1).vmresponse ;
                PkDeflect(j,1) = sweep(j,1).vmbase - sweep(j,1).vmresponse ;
            end        

            if sweep(j,1).currinj > 0
                NrofAPs(j,1) = length(sweep(j,1).ap) ;  
            else
                NrofAPs(j,1) = 0 ;    
            end
            if sweep(j,1).currinj < 0 && isfield(sweep, 'rbap')
                NrofRBAPs(j,1) = length(sweep(j,1).rbap) ;
            else
                NrofRBAPs(j,1) = 0 ;
            end

            if length(sweep(j).ap) >= 4 && sweep(j,1).currinj > 0
                Freqs(j,1) =  mean([sweep(j).ap(4:end).freq]) ;
                StimInts(j,1) = [sweep(j,1).currinj] ;
            end         
        end

        if sum(NrofRBAPs) > 0
            NrofRBAPsM = mean(NrofRBAPs(NrofRBAPs~=0)) ;
        else
            NrofRBAPsM = 0 ;
        end
        


        % find trainsweep 
        % Traincurr=rheobase+50 :
        TrainCurr = sweep(frstspikeswp).epoch(step).stepdiff +50 ;
        for j = 1:NrofSweeps
            tmp(j) = abs(sweep(j).epoch(step).stepdiff - TrainCurr) ;
        end
        
        [TrSwp TrSwp] = min(tmp) ;
        CurrAbvRheo=NaN;
        for TrainSweep = TrSwp:NrofSweeps          
            if length(sweep(TrainSweep).ap) > 3
                isis = [sweep(TrainSweep).ap(2:end).isi];
                stutterAP = [0 0 0 isis(3:end) > 3*isis(2:end-1)];
                stutterISI= [0 0 isis(3:end) > 3*isis(2:end-1)];
               if length(sweep(TrainSweep).ap(~stutterAP)) > 3
                CurrAbvRheo = sweep(TrainSweep).epoch(step).stepdiff - (TrainCurr-50) ;
                TrSweepCrit=1;
                break
               end
            elseif TrainSweep==NrofSweeps
                TrSweepCrit=0;
            end  
        end
        
 

        



        
% Traincurr= sweep with AP freq (AP4:APend) closest to 15 Hz (and reasonable number of spikes):  
%         TrSwp = find(abs(15-Freqs)<5);
%         TrainSweep=[];
%         if ~isempty(TrSwp)
%             for j = 1:length(TrSwp)
%                 nrofaps(j) = sweep(TrSwp(j)).nrofaps;
%             end
%             [~,tmp] = max(nrofaps);
%             TrainSweep=TrSwp(tmp);
%         elseif ~isempty(Freqs)
%             [~,TrainSweep] = min(abs(15-Freqs));
%         end
%         if isempty(TrainSweep)
%             TrainSweep = NrofSweeps;
%         end
%         
%         CurrAbvRheo=NaN;
%         if length(sweep(TrainSweep).ap)>=4
%             TrSweepCrit=1;
%             CurrAbvRheo = sweep(TrainSweep).epoch(step).stepdiff - (sweep(frstspikeswp).epoch(step).stepdiff) ;
%         else
%             TrSweepCrit=0;
%         end
        
        if ~isempty(sweep(TrainSweep).ap)
            for l = 1:length(sweep(TrainSweep).ap)           
                if ~isempty(sweep(TrainSweep).ap(l,1).halfwidth)
                    HWsTS(l,1) = [sweep(TrainSweep).ap(l,1).halfwidth] ;  
                end       
            end
        else
            HWsTS=NaN;
        end
        
        if TrSweepCrit==1
            TSbasetothresh = ([sweep(TrainSweep).ap.thresh]-sweep(TrainSweep).vmbase) ;
            TSpeaktoahp = ([sweep(TrainSweep).ap.ahp_time]-[sweep(TrainSweep).ap.peak_time]); 
            AmpsTSthresh = [sweep(TrainSweep).ap.amp] ;
            AHPsTS = [sweep(TrainSweep).ap.relahp] ;
            AHPslowTS = [sweep(TrainSweep).ap.relahp_slow] ;
            ISIsTS = [sweep(TrainSweep).ap(2:end).isi] ;
            FreqTrSwp = mean([sweep(TrainSweep).ap(4:end).freq]) ;
            NrOfAPsTrSwp = length(sweep(TrainSweep).ap) ; 
            OnsetTSFAP = NaN ;%sweep(TrainSweep).ap(1).thresh_time - (sum(second({sweep(TrainSweep).epoch(1:step-1).timespan}))*1000) ;
        else
            TSbasetothresh = NaN;
            TSpeaktoahp = NaN;
            AmpsTSthresh = NaN;
            AHPsTS = NaN;
            AHPslowTS = NaN;
            ISIsTS = NaN;
            FreqTrSwp = NaN;
            NrOfAPsTrSwp = NaN;
            OnsetTSFAP = NaN;
        end

        % calculate input resistance     
        f_R=fittype('R*x+b');
        tmp=length(find(currInjections_R~=0 & ~isnan(voltageResponses)));
        if tmp > 1         
            [fitR]=fit(currInjections_R(currInjections_R~=0 & voltageResponses~=0 & ~isnan(voltageResponses)),voltageResponses(currInjections_R~=0 & voltageResponses~=0 & ~isnan(voltageResponses)),f_R, 'StartPoint', [0 0]); 
            Rin=fitR.R*1e3; % In mOhm
        elseif tmp ==1  % If there is only one point, use baseline Vm to determine input resistance:
            if ~isnan(sweep([sweep.currinj] == currInjections_R(currInjections_R~=0 & ~isnan(voltageResponses))).vmbase)
                [fitR]=fit([currInjections_R(currInjections_R~=0 & voltageResponses~=0 & ~isnan(voltageResponses)); 0],[voltageResponses(currInjections_R~=0 & voltageResponses~=0 & ~isnan(voltageResponses)); sweep([sweep.currinj] == currInjections_R(currInjections_R~=0 & ~isnan(voltageResponses))).vmbase],f_R, 'StartPoint', [0 0]); 
                Rin=fitR.R*1e3; % In mOhm    
            else
                Rin=NaN;
            end
        else
            Rin=NaN;
        end
        % determine sweep for sag, (minimum voltage response closest to -100)
        tmp = abs(MinVmResponse+100) ;
        [sagswp sagswp] = min(tmp) ;

        % calculate input frequency curve
        Freqs = Freqs(Freqs~=0) ;
        StimInts = StimInts(StimInts~=0) ;
        if length(Freqs) > 1
            [fitFi]=fit(StimInts,Freqs,f_R, 'StartPoint', [0 0]); 
            FrqChngStimInt = fitFi.R ;
        else  
            FrqChngStimInt = NaN ;   
        end

        % bursting & adaptation index
        if TrSweepCrit==1
            isis=ISIsTS(~stutterISI); %exclude stuttering ISIs since that can mess with adaptation analysis
            amps=AmpsTSthresh(~stutterAP);
            hws=HWsTS(~stutterAP);
        else
            isis=ISIsTS;
            amps=AmpsTSthresh;
            hws=HWsTS;
        end
        if length(isis) > 2
            ISIRatio1toAll = mean(isis(2:end)) / mean(isis(1)) ;
            N = length(isis)-1 ;
            for n = 1:N
                ISIchanges(n,1) = (isis(n+1)-isis(n)) / (isis(n+1)+isis(n));
            end
            AdaptIdx = (1/N)*sum(ISIchanges) ;        
        else
            ISIRatio1toAll = NaN;
            AdaptIdx = NaN;
        end
        
        % Amplitude accomodation
        if TrSweepCrit==1

        end
        if length(amps) > 2           
            N = length(amps)-1 ;
            for n = 1:N
                Ampchanges(n,1) = (amps(n+1)-amps(n)) / (amps(n+1)+amps(n));
                HWchanges(n,1) = (hws(n+1)-hws(n)) / (hws(n+1)+hws(n));
            end
            AmpAccom = (1/N)*sum(Ampchanges) ;  
            HWAccom = (1/N)*sum(HWchanges) ;
        else
            AmpAccom = NaN;    
            HWAccom = NaN; 
        end        
        if ~isempty(Freqs)
            freqmax = max(Freqs) ;
        else
            freqmax = NaN ;
        end
             
        %% Create summary  
        Summary(index).File               = abf.filename ;
        Summary(index).Date               = abf.filetimestart ;
        Summary(index).UserID             = abf.userid ;
        Summary(index).guid               = abf.guid ;
        Summary(index).Channel            = NaN ;
        Summary(index).scalefactor        = abf.channel(1).out.scalefactor ;
        Summary(index).holdingcurrent     = abf.channel(1).out.holdingI ;
        Summary(index).holdingvoltage     = abf.channel(1).out.holdingV ;
        Summary(index).NrofSweeps         = NrofSweeps ;
        Summary(index).PDur               = seconds(sweep(1).epoch(step).timespan)*1000 ;
        Summary(index).FrstP              = sweep(1).currinj ;
        Summary(index).DeltaP             = sweep(2).currinj - sweep(1).currinj ;
        Summary(index).Rheobase           = sweep(frstspikeswp).currinj ;
        Summary(index).FrstSpikeSwp       = frstspikeswp ; 
        Summary(index).TrainSwp           = TrainSweep ; 
        Summary(index).CurrAbvRheo        = CurrAbvRheo ;
        Summary(index).vmbaseM            = nanmean(vmbase) ;
        Summary(index).vmbaseSD           = nanstd(vmbase) ;
        Summary(index).Jitter             = nanmean([sweep.jitter]) ;
        Summary(index).InputR             = Rin ;% in MOhm...
        Summary(index).FreqMax            = freqmax ;
        Summary(index).NrOfAPs            = sum(NrofAPs) ;
        Summary(index).NrOfAPsMax         = max(NrofAPs) ; 
        Summary(index).FreqTrSwp          = FreqTrSwp ;
        Summary(index).NrOfAPsTrSwp       = NrOfAPsTrSwp ; 
        Summary(index).FrqChngStimInt     = FrqChngStimInt ;
        Summary(index).NrofRBAPs          = sum(NrofRBAPs) ;
        Summary(index).NrofRBAPsM         = NrofRBAPsM ;
        Summary(index).Sag                = sweep(sagswp,1).epoch(step).sag / PkDeflect(sagswp,1) ;
        Summary(index).VmatSag            = MinVmResponse(sagswp,1) ;
        Summary(index).TauM               = nanmean(taus(taus~=0)) ;
        Summary(index).TauSD              = nanstd(taus(taus~=0)) ;
        Summary(index).OnsetFrstAP        = NaN ;%sweep(frstspikeswp).ap(1).thresh_time - (sum(seconds({sweep(frstspikeswp).epoch(1:step-1).timespan}))*1000) ; 
        Summary(index).ThreshFrstAP       = sweep(frstspikeswp).ap(1).thresh ; 
        Summary(index).FAPbasetothresh    = sweep(frstspikeswp).ap(1).thresh-sweep(frstspikeswp).vmbase ; 
        Summary(index).AmpFAPthresh       = sweep(frstspikeswp).ap(1).amp ;
        Summary(index).FAPpeaktoahp       = sweep(frstspikeswp).ap(1).ahp_time - sweep(frstspikeswp).ap(1).peak_time ;
        Summary(index).HalfWFrstAP        = sweep(frstspikeswp).ap(1).halfwidth ; 
        Summary(index).AHPFrstAP          = sweep(frstspikeswp).ap(1).relahp ;
        Summary(index).AHPslowFrstAP      = sweep(frstspikeswp).ap(1).relahp_slow ;
        Summary(index).UpStrkFrstAP       = sweep(frstspikeswp).ap(1).maxdvdt ;
        Summary(index).DwnStrkFrstAP      = sweep(frstspikeswp).ap(1).mindvdt ;
        Summary(index).UpDwnStrkRatio     = abs(sweep(frstspikeswp).ap(1).maxdvdt) / abs(sweep(frstspikeswp).ap(1).mindvdt) ;
        Summary(index).OnsetTSFAP         = OnsetTSFAP ;  
        Summary(index).TSbasetothreshM    = mean(TSbasetothresh) ; 
        Summary(index).TSbasetothreshSD   = std(TSbasetothresh) ; 
        Summary(index).AmpTSthreshM       = mean(AmpsTSthresh) ; 
        Summary(index).AmpTSthreshSD      = std(AmpsTSthresh) ; 
        Summary(index).TSpeaktoahpM       = nanmean(TSpeaktoahp(TSpeaktoahp~=0)) ; 
        Summary(index).TSpeaktoahpSD      = nanstd(TSpeaktoahp(TSpeaktoahp~=0)) ; 
        Summary(index).HalfWTrSwpM        = nanmean(HWsTS(HWsTS~=0)) ; 
        Summary(index).HalfWTrSwpSD       = nanstd(HWsTS(HWsTS~=0)) ; 
        Summary(index).AHPTrSwpM          = nanmean(AHPsTS(AHPsTS~=0)) ; 
        Summary(index).AHPTrSwpSD         = nanstd(AHPsTS(AHPsTS~=0)) ;
        Summary(index).AHPslowTrSwpM      = nanmean(AHPslowTS(AHPslowTS~=0)) ; 
        Summary(index).AHPslowTrSwpSD     = nanstd(AHPslowTS(AHPslowTS~=0)) ;        
        Summary(index).ISITrSwpM          = mean(ISIsTS(ISIsTS~=0)) ;
        Summary(index).ISITrSwpSD         = std(ISIsTS(ISIsTS~=0)) ;
        Summary(index).ISITrSwpCV         = (std(ISIsTS(ISIsTS~=0)) / mean(ISIsTS(ISIsTS~=0))) ; %coefficient of variation
        Summary(index).AdaptIndexTS       = AdaptIdx ;
        Summary(index).AmpAccomTS         = AmpAccom ;
        Summary(index).HWAccomTS          = HWAccom ;
        Summary(index).ISIRatio1toAll     = ISIRatio1toAll ;
        Summary(index).threshFB0     = nanmean(aps2.thresh(aps2.freqbin==0)) ;
        Summary(index).threshFB1     = nanmean(aps2.thresh(aps2.freqbin==1)) ;
        Summary(index).threshFB2     = nanmean(aps2.thresh(aps2.freqbin==2)) ;
        Summary(index).threshFB3     = nanmean(aps2.thresh(aps2.freqbin==3)) ;
        Summary(index).threshFB4     = nanmean(aps2.thresh(aps2.freqbin==4)) ;
        Summary(index).threshFB5     = nanmean(aps2.thresh(aps2.freqbin==5)) ;    
        Summary(index).HalfWFB0     = nanmean(aps2.halfwidth(aps2.freqbin==0)) ;
        Summary(index).HalfWFB1     = nanmean(aps2.halfwidth(aps2.freqbin==1)) ;
        Summary(index).HalfWFB2     = nanmean(aps2.halfwidth(aps2.freqbin==2)) ;
        Summary(index).HalfWFB3     = nanmean(aps2.halfwidth(aps2.freqbin==3)) ;
        Summary(index).HalfWFB4     = nanmean(aps2.halfwidth(aps2.freqbin==4)) ;
        Summary(index).HalfWFB5     = nanmean(aps2.halfwidth(aps2.freqbin==5)) ;
        Summary(index).upstrokeFB0     = nanmean(aps2.maxdvdt(aps2.freqbin==0)) ;
        Summary(index).upstrokeFB1     = nanmean(aps2.maxdvdt(aps2.freqbin==1)) ;
        Summary(index).upstrokeFB2     = nanmean(aps2.maxdvdt(aps2.freqbin==2)) ;
        Summary(index).upstrokeFB3     = nanmean(aps2.maxdvdt(aps2.freqbin==3)) ;
        Summary(index).upstrokeFB4     = nanmean(aps2.maxdvdt(aps2.freqbin==4)) ;
        Summary(index).upstrokeFB5     = nanmean(aps2.maxdvdt(aps2.freqbin==5)) ;      
        Summary(index).downstrokeFB0     = nanmean(aps2.mindvdt(aps2.freqbin==0)) ;
        Summary(index).downstrokeFB1     = nanmean(aps2.mindvdt(aps2.freqbin==1)) ;
        Summary(index).downstrokeFB2     = nanmean(aps2.mindvdt(aps2.freqbin==2)) ;
        Summary(index).downstrokeFB3     = nanmean(aps2.mindvdt(aps2.freqbin==3)) ;
        Summary(index).downstrokeFB4     = nanmean(aps2.mindvdt(aps2.freqbin==4)) ;
        Summary(index).downstrokeFB5     = nanmean(aps2.mindvdt(aps2.freqbin==5)) ; 
        Summary(index).onsetrapFB0     = nanmean(aps2.onsetrapidity(aps2.freqbin==0)) ;
        Summary(index).onsetrapFB1     = nanmean(aps2.onsetrapidity(aps2.freqbin==1)) ;
        Summary(index).onsetrapFB2     = nanmean(aps2.onsetrapidity(aps2.freqbin==2)) ;
        Summary(index).onsetrapFB3     = nanmean(aps2.onsetrapidity(aps2.freqbin==3)) ;
        Summary(index).onsetrapFB4     = nanmean(aps2.onsetrapidity(aps2.freqbin==4)) ;
        Summary(index).onsetrapFB5     = nanmean(aps2.onsetrapidity(aps2.freqbin==5)) ;  
        Summary(index).ampFB0     = nanmean(aps2.amp(aps2.freqbin==0)) ;
        Summary(index).ampFB1     = nanmean(aps2.amp(aps2.freqbin==1)) ;
        Summary(index).ampFB2     = nanmean(aps2.amp(aps2.freqbin==2)) ;
        Summary(index).ampFB3     = nanmean(aps2.amp(aps2.freqbin==3)) ;
        Summary(index).ampFB4     = nanmean(aps2.amp(aps2.freqbin==4)) ;
        Summary(index).ampFB5     = nanmean(aps2.amp(aps2.freqbin==5)) ;  
        Summary(index).ahpFB0     = nanmean(aps2.relahp(aps2.freqbin==0)) ;
        Summary(index).ahpFB1     = nanmean(aps2.relahp(aps2.freqbin==1)) ;
        Summary(index).ahpFB2     = nanmean(aps2.relahp(aps2.freqbin==2)) ;
        Summary(index).ahpFB3     = nanmean(aps2.relahp(aps2.freqbin==3)) ;
        Summary(index).ahpFB4     = nanmean(aps2.relahp(aps2.freqbin==4)) ;
        Summary(index).ahpFB5     = nanmean(aps2.relahp(aps2.freqbin==5)) ; 
        Summary(index).updwnratioFB0     = nanmean(aps2.updownratio(aps2.freqbin==0)) ;
        Summary(index).updwnratioFB1     = nanmean(aps2.updownratio(aps2.freqbin==1)) ;
        Summary(index).updwnratioFB2     = nanmean(aps2.updownratio(aps2.freqbin==2)) ;
        Summary(index).updwnratioFB3     = nanmean(aps2.updownratio(aps2.freqbin==3)) ;
        Summary(index).updwnratioFB4     = nanmean(aps2.updownratio(aps2.freqbin==4)) ;
        Summary(index).updwnratioFB5     = nanmean(aps2.updownratio(aps2.freqbin==5)) ;
        Summary(index).amp1to20  = nanmean(aps2.amp(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).upstroke1to20  = nanmean(aps2.maxdvdt(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).downstroke1to20  = nanmean(aps2.mindvdt(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).HalfW1to20  = nanmean(aps2.halfwidth(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).thresh1to20  = nanmean(aps2.thresh(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).onsetrap1to20  = nanmean(aps2.onsetrapidity(aps2.freqbin>=1 & aps2.freqbin<=2)) ;
        Summary(index).updwnratio1to20  = nanmean(aps2.updownratio(aps2.freqbin>=1 & aps2.freqbin<=2)) ;   
        Summary(index).amp21to40  = nanmean(aps2.amp(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).upstroke21to40  = nanmean(aps2.maxdvdt(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).downstroke21to40  = nanmean(aps2.mindvdt(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).HalfW21to40  = nanmean(aps2.halfwidth(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).thresh21to40  = nanmean(aps2.thresh(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).onsetrap21to40  = nanmean(aps2.onsetrapidity(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
        Summary(index).updwnratio21to40  = nanmean(aps2.updownratio(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
        
        %end %if at line 64
        %clear variables assigned in "sweep" For loop
        clearvars -except abf Summary i basedir savename abfs aps channels ins outs epochs sweeps index
        index = index + 1 ;
    %end

%% save
save(fullfile(basedir, savename), 'Summary') ;
clearvars -except Summary i










