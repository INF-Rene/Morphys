%% test script
close all, clear all

%% Load data from CSV tables into nested struct
% requires specific folder- and filenames in location basedir!
basedir = 'C:\Users\DBHeyer\Documents\PhD\Human Database\Test' ;
structname = 'NestedStruct2.mat' ;
summaryname = 'DataSummary.mat' ;

abf = ImportCSVstoStruct(basedir,structname) ;

%% Make selection of stepprotocols

%% Loop through abfs
  
for i = 1:length(abf)
    sweep = abf(i).channel.in.sweep ;
    NrofSweeps = length(sweep) ;  
    % find current injection epoch and assign aps to sweep
    for step = 1:length(sweep(1).epoch)
        if sweep(1).epoch(step).stepdiff < 0 && (sweep(1).epoch(step).stepdiff + sweep(1).epoch(step+1).stepdiff) == 0
            break
        end
    end
    for j = 1:NrofSweeps
        sweep(j).vmbase = sweep(j).epoch(step-1).steadystate ;
        sweep(j).currinj = sweep(j).epoch(step).stepdiff ;
        sweep(j).vmresponse = sweep(j).epoch(step).vstep ;
        sweep(j).ap = sweep(j).epoch(step).ap ;
        for ap = 1:length(sweep(j).epoch(step+1).ap)
            if sweep(j).epoch(step+1).ap(ap).start_time > (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+7 && sweep(j).epoch(step+1).ap(ap).start_time < (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+300
                sweep(j).rbap(ap) = sweep(j).epoch(step+1).ap(ap) ;
            end
        end
    end
    % find rheobase sweep
    for frstspikeswp = 1:NrofSweeps
        if sweep(frstspikeswp).epoch(step).nrofaps > 0 && sweep(frstspikeswp).epoch(step).stepdiff > 0
            break
        end
    end
    % find trainsweep
    TrainCurr = sweep(frstspikeswp).epoch(step).stepdiff +50 ;
    for j = 1:NrofSweeps
        tmp(j) = abs(sweep(j).epoch(step).stepdiff - TrainCurr) ;
    end
   
    [TrSwp TrSwp] = min(tmp) ;
 
    for TrainSweep = TrSwp:NrofSweeps      
        if length(sweep(TrainSweep).ap) > 3           
            CurrAbvRheo = sweep(TrainSweep).epoch(step).stepdiff - (TrainCurr-50) ;
            break
        end  
    end
   
    % calculate variables
    vmbase = [sweep.vmbase] ;
    
    for j = 1:NrofSweeps           
        if sweep(j,1).currinj >= -100 && sweep(j,1).currinj < 0
            voltageResponses(j,1) = sweep(j,1).vmresponse ; 
            currInjections_R(j,1) = sweep(j,1).currinj ;
            if sweep(j,1).epoch(step).tau < 100 && sweep(j,1).epoch(step).tau > 0
                taus(j,1) = sweep(j,1).epoch(step).tau ;
            end
        end       

        if ~isempty(sweep(j,1).vmresponse) 
            MinVmResponse(j,1) = sweep(j,1).vmresponse ;
            PkDeflect(j,1) = sweep(j,1).vmbase - sweep(j,1).vmresponse ;
        end        
        
        if sweep(j,1).currinj > 0
            NrofAPs(j,1) = length(sweep(j,1).ap) ;  
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
    
    for l = 1:length(sweep(TrainSweep).ap)           
        if ~isempty(sweep(TrainSweep).ap(l,1).halfwidth)
            HWsTS(l,1) = [sweep(TrainSweep).ap(l,1).halfwidth] ;  
        end       
    end
    
    TSbasetothresh = ([sweep(TrainSweep).ap.thresh]-sweep(TrainSweep).vmbase) ; 
    TSpeaktoahp = ([sweep(TrainSweep).ap.ahp_time]-[sweep(TrainSweep).ap.peak_time]) ; 
    AmpsTSthresh = [sweep(TrainSweep).ap.amp] ;               
    AHPsTS = [sweep(TrainSweep).ap.relahp] ;  
    ISIsTS = [sweep(TrainSweep).ap(2:end).isi] ;    
          

    % calculate input resistance     
    f_R=fittype('R*x+b');
    [fitR]=fit(currInjections_R(currInjections_R~=0),voltageResponses(voltageResponses~=0),f_R, 'StartPoint', [0 0]); 

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
    if length(ISIsTS) > 2
        ISIRatio1toAll = mean(ISIsTS(2:end)) / mean(ISIsTS(1)) ;
        N = length(ISIsTS)-1 ;
        for n = 1:N
            ISIchanges(n,1) = (ISIsTS(n+1)-ISIsTS(n)) / (ISIsTS(n+1)+ISIsTS(n));
        end
        AdaptIdx = (1/N)*sum(ISIchanges) ;        
    else
        ISIRatio1toAll = NaN;
        AdaptIdx = NaN;
    end

    %% Create summaries   
    Summary(i).File               = abf(i).filename ;
    Summary(i).Date               = abf(i).filetimestart ;
    Summary(i).UserID             = abf(i).userid ;
    Summary(i).guid               = abf.guid ;
    Summary(i).Channel            = abf(i).channel.in.number ;
    Summary(i).NrofSweeps         = NrofSweeps ;
    Summary(i).PDur               = second(sweep(1).epoch(step).timespan)*1000 ;
    Summary(i).FrstP              = sweep(1).currinj ;
    Summary(i).DeltaP             = sweep(2).currinj - sweep(1).currinj ;
    Summary(i).Rheobase           = sweep(frstspikeswp).currinj ;
    Summary(i).FrstSpikeSwp       = frstspikeswp ; 
    Summary(i).TrainSwp           = TrainSweep ; 
    Summary(i).CurrAbvRheo        = CurrAbvRheo ;
    Summary(i).vmbaseM            = nanmean(vmbase) ;
    Summary(i).vmbaseSD           = nanstd(vmbase) ;
    Summary(i).InputR             = fitR.R*1e3 ;% in MOhm...
    Summary(i).FreqMax            = max(Freqs) ;
    Summary(i).NrOfAPsMax         = max(NrofAPs) ; 
    Summary(i).FreqTrSwp          = mean([sweep(TrainSweep).ap(4:end).freq]) ;
    Summary(i).NrOfAPsTrSwp       = length(sweep(TrainSweep).ap) ; 
    Summary(i).FrqChngStimInt     = FrqChngStimInt ;
    Summary(i).NrofRBAPs          = sum(NrofRBAPs) ;
    Summary(i).NrofRBAPsM         = NrofRBAPsM ;
    Summary(i).Sag                = sweep(sagswp,1).epoch(step).sag / PkDeflect(sagswp,1) ;
    Summary(i).VmatSag            = MinVmResponse(sagswp,1) ;
    Summary(i).TauM               = nanmean(taus(taus~=0)) ;
    Summary(i).TauSD              = nanstd(taus(taus~=0)) ;
    Summary(i).OnsetFrstAP        = sweep(frstspikeswp).ap(1).thresh_time - (sum(second({sweep(frstspikeswp).epoch(1:step-1).timespan}))*1000) ; 
    Summary(i).ThreshFrstAP       = sweep(frstspikeswp).ap(1).thresh ; 
    Summary(i).FAPbasetothresh    = sweep(frstspikeswp).ap(1).thresh-sweep(frstspikeswp).vmbase ; 
    Summary(i).AmpFAPthresh       = sweep(frstspikeswp).ap(1).amp ;
    Summary(i).FAPpeaktoahp       = sweep(frstspikeswp).ap(1).ahp_time - sweep(frstspikeswp).ap(1).peak_time ;
    Summary(i).HalfWFrstAP        = sweep(frstspikeswp).ap(1).halfwidth ; 
    Summary(i).AHPFrstAP          = sweep(frstspikeswp).ap(1).relahp ; 
    Summary(i).UpStrkFrstAP       = sweep(frstspikeswp).ap(1).maxdvdt ;
    Summary(i).DwnStrkFrstAP      = sweep(frstspikeswp).ap(1).mindvdt ;
    Summary(i).UpDwnStrkRatio     = abs(sweep(frstspikeswp).ap(1).maxdvdt) / abs(sweep(frstspikeswp).ap(1).mindvdt) ;
    Summary(i).OnsetTSFAP         = sweep(TrainSweep).ap(1).thresh_time - (sum(second({sweep(TrainSweep).epoch(1:step-1).timespan}))*1000) ;  
    Summary(i).TSbasetothreshM    = mean(TSbasetothresh) ; 
    Summary(i).TSbasetothreshSD   = std(TSbasetothresh) ; 
    Summary(i).AmpTSthreshM       = mean(AmpsTSthresh) ; 
    Summary(i).AmpTSthreshSD      = std(AmpsTSthresh) ; 
    Summary(i).TSpeaktoahpM       = mean(TSpeaktoahp) ; 
    Summary(i).TSpeaktoahpSD      = std(TSpeaktoahp) ; 
    Summary(i).HalfWTrSwpM        = mean(HWsTS) ; 
    Summary(i).HalfWTrSwpSD       = std(HWsTS) ; 
    Summary(i).AHPTrSwpM          = mean(AHPsTS) ; 
    Summary(i).AHPTrSwpSD         = std(AHPsTS) ;
    Summary(i).ISITrSwpM          = mean(ISIsTS(ISIsTS~=0)) ;
    Summary(i).ISITrSwpSD         = std(ISIsTS(ISIsTS~=0)) ;
    Summary(i).ISITrSwpCV         = (std(ISIsTS(ISIsTS~=0)) / mean(ISIsTS(ISIsTS~=0))) ; %coefficient of variation
    Summary(i).AdaptIndexTS       = AdaptIdx ;
    Summary(i).ISIRatio1toAll     = ISIRatio1toAll ;

    %clear variables assigned in "sweep" For loop
    clearvars -except abf Summary i basedir structname summaryname
end

%% save
save(fullfile(basedir, summaryname), 'abf', 'Summary') ;











