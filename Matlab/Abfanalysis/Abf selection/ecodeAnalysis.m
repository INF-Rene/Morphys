function [Summary, missing] = ecodeAnalysis(cellname, celldir, destinpath, cellinfo)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here





a=dir(celldir);
tmp=endsWith({a.name}, '.abf');
files=fullfile({a(tmp).folder},{a(tmp).name});
Summary=table;
Summary.Cellname(1)           = {cellname};
missing=[1,1];  
if isempty(files), return, end

if isempty(cellname)
    cellname=strsplit(celldir, {'\'});
    cellname=cellname{end};
end
%% load files
a={};

for i=1:numel(files)
    fprintf('Open file %1.0f out of %1.0f \n', i, numel(files))
    [~, ~,h]=abfload_pro(files{i});
    a{i}=h.stringSection;
end

if ~contains(h.stringSection, 'RawOutPut')
    ss=load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF700A.mat');
    ss=ss.obj;
else
    ss=load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
    ss=ss.obj;
end


tmp1=find(contains(a, 'LSCOARSE'));
tmp2=find(contains(a, 'LSFINEST'));
if ~isempty(tmp1)
    a1=Abffile(files{tmp1(end)}, ss);
    a1=a1.analyseabf;
    missing(1)=0;
else
    fprintf('LSCOARSE of cell %s is missing\n', cellname)
end
if ~isempty(tmp2)
    a2=Abffile(files{tmp2(end)}, ss);
    a2=a2.analyseabf;
    missing(2)=0;
else
    fprintf('LSFINEST of cell %s is missing\n', cellname)
end

if all(missing==1)
    Summary=table;
    Summary.Cellname           = cellname;
    return
end

    

%% make and save an overviewplot
close all
figure ('position', [0 0 1920 1000])

%table with PatchSeqInfo on top

T=cellinfo(1,ismember(cellinfo.Properties.VariableNames, {'CellName', 'MorphologyPresent_Yes_no', 'TransciptomicType','H_score',...
    'ResolutionIndex','NormMarkerSum__0_40_', 'Cluster'}));
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
TString = TString(1:end-2);
FixedWidth = get(0,'FixedWidthFontName');
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','Position',[0.04 0.975 0.90 0.02 ],...
    'FitBoxToText', 'on', 'Margin', 1);
%text(0,0, TString)



if missing(1)==0
    %passive plot
    subplot('position', [0.04 0.55 0.2 0.35 ])
    a1.getchannel.getin('signal', 'primary').getsweep.plotanalysis
    xlim([900 2800])
    xlabel 'Time (ms)';
    ylabel 'Vm (mV)';
    title('Passive properties')
end
if missing(2)==0 && ~any([a2.getchannel.getin('signal', 'primary').getsweep.getepoch('Name', 'Epoch D').nrofaps]>0), missing(2)=1; end
if missing(2)==0 
    %First AP
    subplot('position', [0.275 0.55 0.2 0.35 ])
    FrstSpikeSwp=find([a2.getchannel.getin('signal', 'primary').getsweep.getepoch('Name', 'Epoch D').nrofaps]>0,1);
    ap=a2.getchannel.getin('signal', 'primary').getsweep(FrstSpikeSwp).getepoch('Name', 'Epoch D').getap(1);
    ap.plotanalysis2;
    xlim([ap.start_time-1 ap.start_time+8]);
    title('First AP')

    %I-F plot
    steps=[a2.getchannel.getin('signal', 'primary').getsweep.getepoch('Name', 'Epoch D').amplitude];
    freqs=NaN(1,a2.nrofsweeps);
    nrofAPs=NaN(1,a2.nrofsweeps);
    for i=1:a2.nrofsweeps
        nrofAPs(i)=a2.getchannel.getin('signal', 'primary').getsweep(i).getepoch('Name', 'Epoch D').nrofaps;
        if nrofAPs(i)>1
            freqs(i)=nanmean([a2.getchannel.getin('signal', 'primary').getsweep(i).getepoch('Name', 'Epoch D').getap(4:end).freq]);
        end
    end
    subplot('position', [0.52 0.55 0.19 0.35 ])
    yyaxis left
    if any(~isnan(freqs))
        plot(steps([1 3:end])./steps(1)*100,freqs([1 3:end]), '--o')
        hold on
        if numel(steps)>1, scatter(steps(2)./steps(1)*100,freqs([2])), end
        xlabel('I/Ithresh(%)')
        ylabel('Mean Inst. Frequency of APs 4-end (Hz)')
        ylim([0 max(freqs)*1.1])
    end
    yyaxis right;
    plot(steps([1 3:end])./steps(1)*100,nrofAPs([1 3:end]), '--o')
    hold on
    if numel(steps)>1, scatter(steps(2)./steps(1)*100,nrofAPs([2])), end
    ylabel('Nr of APs')
    ylim([0 max(nrofAPs)*1.1])
    title('I vs steady-state firing rate')

    %adaptation & bursting behaviour plot
    subplot('position', [0.77 0.55 0.19 0.35 ])
    brst=NaN(1,numel(steps));
    adapt=NaN(1,numel(steps));
    for j=1:numel(steps)
        if nrofAPs(j)>1
            isis=[a2.getchannel.getin.getsweep(j).getepoch('Name', 'Epoch D').getap.isi];
            isis=isis(~isnan(isis));
            if numel(isis)>2
                brst(j)=nanmean(isis(2:end))/isis(1);
                disi=diff(isis);
                adapt(j)= nanmean(disi./(isis(1:end-1)+isis(2:end)));
            end
        end
    end
    yyaxis left;
    if any(~isnan(brst))
        plot(steps([1 3:end])./steps(1)*100,brst([1 3:end]), '--o')
        hold on
        if numel(steps)>1, scatter(steps(2)./steps(1)*100,brst([2])), end
        xlabel('I/Ithresh(%)')
        ylabel('Burst Index')
        ylim([0 max(brst)*1.1])
    end
    yyaxis right;
    if any(~isnan(adapt))
        plot(steps([1 3:end])./steps(1)*100,adapt([1 3:end]), '--o')
        hold on
        if numel(steps)>1, scatter(steps(2)./steps(1)*100,adapt([2])), end
        ylabel('Adaptation index')
        ylim([-inf max(adapt)*1.1])
    end
    title('Bursting & Adaptation')
    
    %Sweep 1, 2, 5
    subplot('position', [0.04 0.10 0.2 0.35 ])
    a2.getchannel.getin('signal', 'primary').getsweep(1).plotanalysis
    if a2.nrofsweeps>1  
        a2.getchannel.getin('signal', 'primary').getsweep(2).plot('Color', [0.6 0.6 0.6])
    end
    xlim([900 2800])
    rheo=a2.getchannel.getin('signal', 'primary').getsweep(1).getepoch('Name', 'Epoch D').amplitude;
    title(['Rheobase (' num2str(rheo) ' pA) and Rheo-10 pA'])
    
    if a2.nrofsweeps > 2
        subplot('position', [0.28 0.10 0.2 0.35 ])
        a2.getchannel.getin('signal', 'primary').getsweep(3).plotanalysis
        xlim([900 2800])
        title('Rheobase + 40')
    end
    if a2.nrofsweeps > 4
        subplot('position', [0.54 0.10 0.2 0.35 ])
        a2.getchannel.getin('signal', 'primary').getsweep(5).plotanalysis
        xlim([900 2800])
        title('Rheobase + 120')
    elseif a2.nrofsweeps == 4
        subplot('position', [0.54 0.10 0.2 0.35 ])
        a2.getchannel.getin('signal', 'primary').getsweep(4).plotanalysis
        xlim([900 2800])
        title('Rheobase + 80')
    end
     
    %at least 10 APs (or maximum)
    subplot('position', [0.78 0.10 0.2 0.35 ])
    AP10=find(nrofAPs>=10,1);
    if isempty(AP10), [~,AP10]=max(nrofAPs); end
    a2.getchannel.getin('signal', 'primary').getsweep(AP10).plotanalysis
    xlim([900 2800])
    title(['first 10 AP (' num2str(steps(AP10)) ' pA)'])
end

if ~isempty(destinpath)
    saveas(gcf, [destinpath, cellname, '.fig'])
    print([destinpath, cellname, '.jpg'],  '-djpeg', '-r300');
end

if all(missing==0) && a2.nrofsweeps>1
    Summary=table2struct(Summary);
    %% Get CCstep summary data
    index=1;
    sweep1 = a1.getchannel.getin('signal','primary').getsweep ;
    sweep2 = a2.getchannel.getin('signal','primary').getsweep ;
    sweep=[sweep1, sweep2([2 1 3:end])];
    NrofSweeps = length(sweep) ;
    % find current injection epoch and assign aps to sweep
    for step = 1:length(sweep(1).getepoch)
        if sweep(1).getepoch(step).stepdiff < 0 && abs(sweep(1).getepoch(step).stepdiff + sweep(1).getepoch(step+1).stepdiff) < 5
            break
        end
    end
    for j = 1:NrofSweeps
        swp(j).vmbase = sweep(j).getepoch(step-1).steadystate ;
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
        Summary(index).Date               = a2.filetimestart ;
        Summary(index).UserID             = a2.userid ;
        Summary(index).guid               = a2.guid;
        Summary(index).Channel            = a2.getchannel.number;
        Summary(index).scalefactor        = a2.getchannel.getout.scalefactor;
        Summary(index).holdingcurrent     = a2.getchannel.getout.holdingI ;
        Summary(index).holdingvoltage     = a2.getchannel.getout.holdingV ;
        Summary(index).NrofSweeps         = NrofSweeps ;
        Summary(index).PDur               = milliseconds(sweep(1).getepoch(step).timespan);
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
        Summary(index).OnsetFrstAP        = sweep(frstspikeswp).getepoch(step).getap(1).thresh_time - sum([sweep(frstspikeswp).getepoch(1:step-1).timespan]) ;
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

end

