function [Summary] = CCstepSummary(a, cellname)
    ss=load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
    ss=ss.obj;
    Summary=struct;
    
    if isa(a, 'Abffile') || iscell(a)
        if numel(a)==1
        sweep=a.getchannel.getin('signal', 'primary').getsweep;
        a2=a;
        else
            a1=a{1};
            a2=a{2};
            sweep1 = a1.getchannel.getin('signal','primary').getsweep ;
            sweep2 = a2.getchannel.getin('signal','primary').getsweep ;
            sweep=[sweep1, sweep2([2 1 3:end])];
        end
    elseif isa(a, 'Spike2Channel')
        sweep=a.getin('signal', 'primary').getsweep;
        a2=a;
    end
    
    
    %% Get CCstep summary data
    index=1;

    NrofSweeps = length(sweep) ;
    % find current injection epoch and assign aps to sweep
    for step = 1:length(sweep(1).getepoch)
        if sweep(1).getepoch(step).stepdiff < 0 && abs(sweep(1).getepoch(step).stepdiff + sweep(1).getepoch(step+1).stepdiff) < 5
            break
        end
    end
    for j = 1:NrofSweeps
        swp(j).vmbase = sweep(j).getepoch(step-1).steadystate ;
        if isnan(swp(j).vmbase)
            preAPtimes=sweep(j).getepoch(step-1).getap.peak_time;
            stepstart=sweep(j).getepoch(step).TimeInfo.Start;
            if ~isempty(preAPtimes) && preAPtimes(end)+500 < stepstart
                swp(j).vmbase = sweep(j).getepoch(step-1).getsampleusingtime(stepstart-200,stepstart-1).median;
            end
        end
        swp(j).jitter = sweep(j).getepoch(step-1).jitter ;
        if sweep(j).nrofepochs>=step
            swp(j).currinj = sweep(j).getepoch(step).stepdiff ;
            swp(j).vmresponse = sweep(j).getepoch(step).vstep ;
            %swp(j).ap = sweep(j).getepoch(step).ap ;
        else
            swp(j).currinj=0;
            swp(j).vmresponse=NaN;
        end
    end

    % find rheobase sweep
    apcrit=0;
    for frstspikeswp = 1:NrofSweeps
        if sweep(frstspikeswp).nrofepochs>=step && sweep(frstspikeswp).getepoch(step).nrofaps > 0 && sweep(frstspikeswp).getepoch(step).stepdiff > 0
            apcrit=1;
            break
        end
    end

    if apcrit==1
        % calculate variables
        vmbase = [swp.vmbase] ;
        Freqs=[];
        StimInts=[];
        for j = 1:NrofSweeps
            if swp(j).currinj >= -104 && swp(j).currinj < 0
                voltageResponses(j,1) = swp(j).vmresponse ;
                currInjections_R(j,1) = swp(j).currinj ;
                if sweep(j).getepoch(step).tau < 100 && sweep(j).getepoch(step).tau > 0 && sweep(j).getepoch(step).gof > 0.95
                    taus(j,1) = sweep(j).getepoch(step).tau ;
                else
                    taus(j,1) = NaN;
                end
            end

            if ~isempty(swp(j).vmresponse) && ~isnan(swp(j).vmresponse)
                MinVmResponse(j,1) = swp(j).vmresponse ;
                PkDeflect(j,1) = swp(j).vmbase - swp(j).vmresponse ;
            end

            if swp(j).currinj > 0
                NrofAPs(j,1) = sweep(j).nrofaps ;
            end
            if  swp(j).currinj < 0 && sweep(j).getepoch(step+1).nrofaps>0 %count RBAPs if they occur after the step within a 7-300 ms window
                NrofRBAPs(j,1)=numel(find([sweep(j).getepoch(step+1).getap.start_time] > sum(milliseconds([sweep(1).getepoch(1:step).timespan]))+7&...
                [sweep(j).getepoch(step+1).getap.start_time] < sum(milliseconds([sweep(1).getepoch(1:step).timespan]))+300));
            else
                NrofRBAPs(j,1) = 0 ;
            end

            if sweep(j).nrofepochs >= step && sweep(j).getepoch(step).nrofaps >= 4 && swp(j).currinj > 0
                Freqs(j,1) =  mean([sweep(j).getepoch(step).getap(4:end).freq]) ;
                StimInts(j,1) = [swp(j).currinj] ;
            end
        end

        if sum(NrofRBAPs) > 0
            NrofRBAPsM = mean(NrofRBAPs(NrofRBAPs~=0)) ;
        else
            NrofRBAPsM = 0 ;
        end



        % find trainsweep
        % Traincurr=rheobase+50 :
        TrainCurr = sweep(frstspikeswp).getepoch(step).stepdiff + 40 ;
        tmp = abs([swp.currinj] - TrainCurr) ;


        [~, TrSwp] = min(tmp) ;
        CurrAbvRheo=NaN;
        for TrainSweep = TrSwp:NrofSweeps
            if sweep(TrainSweep).getepoch(step).nrofaps > 3
                isis = [sweep(TrainSweep).getepoch(step).getap(2:end).isi];
                stutterAP = [0 0 0 isis(3:end) > 3*isis(2:end-1)];
                stutterISI= [0 0 isis(3:end) > 3*isis(2:end-1)];
                if numel(sweep(TrainSweep).getepoch(step).getap(~stutterAP)) > 3
                    CurrAbvRheo = sweep(TrainSweep).getepoch(step).stepdiff - (TrainCurr-50) ;
                    TrSweepCrit=1;
                    break
                elseif TrainSweep==NrofSweeps
                    TrSweepCrit=0;
                end
            elseif TrainSweep==NrofSweeps
                TrSweepCrit=0;
            end
        end

        if sweep(TrainSweep).getepoch(step).nrofaps>0
            for l = 1:sweep(TrainSweep).getepoch(step).nrofaps
                HWsTS(l,1) = sweep(TrainSweep).getepoch(step).getap(l).halfwidth ;
            end
        else
            HWsTS=NaN;
        end

        if TrSweepCrit==1
            TSbasetothresh = ([sweep(TrainSweep).getepoch(step).getap.thresh]-swp(TrainSweep).vmbase) ;
            TSpeaktoahp = ([sweep(TrainSweep).getepoch(step).getap.ahp_time]-[sweep(TrainSweep).getepoch(step).getap.peak_time]);
            AmpsTSthresh = [sweep(TrainSweep).getepoch(step).getap.amp] ;
            AHPsTS = [sweep(TrainSweep).getepoch(step).getap.relahp] ;
            AHPslowTS = [sweep(TrainSweep).getepoch(step).getap.relahp_slow] ;
            ISIsTS = [sweep(TrainSweep).getepoch(step).getap(2:end).isi] ;
            FreqTrSwp = mean([sweep(TrainSweep).getepoch(step).getap(4:end).freq]) ;
            NrOfAPsTrSwp = length(sweep(TrainSweep).getepoch(step).getap) ;
            OnsetTSFAP = sweep(TrainSweep).getepoch(step).getap(1).thresh_time - sweep(TrainSweep).getepoch(step).TimeInfo.Start ;
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
            if ~isnan([swp([swp([swp.currinj]<0 & [swp.currinj]>-104).currinj]' == currInjections_R(currInjections_R~=0 & ~isnan(voltageResponses))).vmbase])
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
        [~, sagswp] = min(tmp) ;

        % calculate input frequency curve
        Freqs = Freqs(Freqs~=0) ;
        StimInts = StimInts(StimInts~=0) ;
        if length(Freqs) > 1
            [fitFi]=fit(StimInts,Freqs,f_R, 'StartPoint', [0 0]);
            FrqChngStimInt = fitFi.R ;
        else
            FrqChngStimInt = NaN ;
        end
        if isempty(Freqs), Freqs=NaN; end

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


        %% Create summary
        Summary(index).Cellname           = cellname;
        Summary(index).File               = a2.filename;
        Summary(index).guid               = a2.guid;
        %Summary(index).Date               = a2.filetimestart ;
        if isa(a2, 'Abffile')
            Summary(index).Date               = year(a2.filetimestart)*1e4+month(a2.filetimestart)*1e2+day(a2.filetimestart) ;
            Summary(index).UserID             = a2.userid ;
            Summary(index).Channel            = a2.getchannel.number;
            Summary(index).scalefactor        = a2.getchannel.getout.scalefactor;
            Summary(index).holdingcurrent     = a2.getchannel.getout.holdingI ;
            Summary(index).holdingvoltage     = a2.getchannel.getout.holdingV ;
            Summary(index).PDur               = milliseconds(sweep(1).getepoch(step).timespan);
        elseif isa(a, 'Spike2Channel')
            Summary(index).Date               = 0;
            Summary(index).UserID             = 'RWS';
            Summary(index).Channel            = 1;
            Summary(index).scalefactor        = 1;
            Summary(index).holdingcurrent     = a2.getin('signal', 'Secondary').getsweep(1).getepoch(1).median;
            Summary(index).holdingvoltage     = NaN ;
            Summary(index).PDur               = sweep(1).getepoch(step).timespan;
        end
        Summary(index).NrofSweeps         = NrofSweeps ;  
        Summary(index).FrstP              = swp(1).currinj ;
        Summary(index).DeltaP             = mean(diff([swp.currinj])) ;
        Summary(index).Rheobase           = swp(frstspikeswp).currinj ;
        Summary(index).FrstSpikeSwp       = frstspikeswp ;
        Summary(index).TrainSwp           = TrainSweep ;
        Summary(index).CurrAbvRheo        = CurrAbvRheo ;
        Summary(index).vmbaseM            = nanmean(vmbase) ;
        Summary(index).vmbaseSD           = nanstd(vmbase) ;
        Summary(index).Jitter             = nanmean([swp.jitter]) ;
        Summary(index).InputR             = Rin ;% in MOhm...
        Summary(index).FreqMax            = max(Freqs) ;
        Summary(index).NrOfAPsMax         = max(NrofAPs) ;
        Summary(index).FreqTrSwp          = FreqTrSwp ;
        Summary(index).NrOfAPsTrSwp       = NrOfAPsTrSwp ;
        Summary(index).FrqChngStimInt     = FrqChngStimInt ;
        Summary(index).NrofRBAPs          = sum(NrofRBAPs) ;
        Summary(index).NrofRBAPsM         = NrofRBAPsM ;
        Summary(index).Sag                = sweep(sagswp).getepoch(step).sag / PkDeflect(sagswp,1) ;
        Summary(index).VmatSag            = MinVmResponse(sagswp,1) ;
        Summary(index).TauM               = nanmean(taus(taus~=0)) ;
        Summary(index).TauSD              = nanstd(taus(taus~=0)) ;
        Summary(index).OnsetFrstAP        = sweep(frstspikeswp).getepoch(step).getap(1).thresh_time - sweep(frstspikeswp).getepoch(step).TimeInfo.Start ;
        Summary(index).ThreshFrstAP       = sweep(frstspikeswp).getepoch(step).getap(1).thresh ;
        Summary(index).FAPbasetothresh    = sweep(frstspikeswp).getepoch(step).getap(1).thresh-swp(frstspikeswp).vmbase ;
        Summary(index).AmpFAPthresh       = sweep(frstspikeswp).getepoch(step).getap(1).amp ;
        Summary(index).FAPpeaktoahp       = sweep(frstspikeswp).getepoch(step).getap(1).ahp_time - sweep(frstspikeswp).getepoch(step).getap(1).peak_time ;
        Summary(index).HalfWFrstAP        = sweep(frstspikeswp).getepoch(step).getap(1).halfwidth ;
        Summary(index).AHPFrstAP          = sweep(frstspikeswp).getepoch(step).getap(1).relahp ;
        Summary(index).AHPslowFrstAP      = sweep(frstspikeswp).getepoch(step).getap(1).relahp_slow ;
        Summary(index).UpStrkFrstAP       = sweep(frstspikeswp).getepoch(step).getap(1).maxdvdt ;
        Summary(index).DwnStrkFrstAP      = sweep(frstspikeswp).getepoch(step).getap(1).mindvdt ;
        Summary(index).UpDwnStrkRatio     = abs(sweep(frstspikeswp).getepoch(step).getap(1).maxdvdt) / abs(sweep(frstspikeswp).getepoch(step).getap(1).mindvdt) ;
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

        Summary=struct2table(Summary, 'AsArray', true);
    end
    % if ~isempty(destinpath)
    %     save([destinpath,cellname, '.mat'],Summary)
    % end

end

