clear all
clc
%% Bulk conversion for NWB2 files - bulk nwb2 analysis
% set path to find .NWB files
dir_abfs          = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9';
% set path to save analyzed files
dir_mats_analysed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/done'
% set path to save figures
dir_figs_analysed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/PatchSeq analysed/HUMAN/NOT WORKING/new nwb2/figures'
% get a list of all .nwb files in the specified path
fileinfo  = dir(fullfile(dir_abfs,'*.nwb'));
filelist  = {fileinfo.name};

%% NWB2 Convertion and analysis
% loop through filelist
errorMessages1=[];%keep a list of any errors, Anna 01-Dec
verbose=1%keep a list of any errors, Anna 01-Dec
for i=1:size(filelist,2)
    try%try to process data for a given session
        
    % add the right extensions to the filename for saving and loading
    fn    = filelist{i}(1:end-4);
    fnabf = [fn '.nwb'];
    fnmat = [fn '.mat'];
    
    % load file and convert to NWBfile object. 
    % use the final argument to specify which stimsets to include
    %a = NWBfile(fullfile(dir_abfs,fnabf), [{'LP'} {'hresh'} {'LSFINEST'} {'LSCOARSE'}]) ;
    a = NWBfile(fullfile(dir_abfs,fnabf), [{'hresh'} {'teps'} {'LSFINEST'} {'LSCOARSE'}]) ;
    %a = NWBfile(fullfile(dir_abfs,fnabf), [{'Rheo'} {'CCSteps'}]) ;
    % analyse NWBfile object
    a = a.analyseNWB ;
    % save the object to .mat file to be used for further analysis
    a.saveme(dir_mats_analysed,fnmat);
    %a.plotanalysis(dir_figs_analysed)
    
    %Anna_edition 1-Dec - handling errors
    catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end        
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages1=[errorMessages1; filelist{i} {ME}];%append the error message to a list
        end
    end    
end    

%General end error message - Anna 1-Dec
if ~isempty(errorMessages1)
    for k=1:size(errorMessages1,1)
        fprintf('\nError at session %d\n',errorMessages1{k,1});%print problematic session number to screen
        errorMessages1{k,2}.message%display error message
        errorMessages1{k,2}.stack%print function and line number at which error occurred
    end
end
save('/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/errorMessages_nwb.mat','errorMessages1');




%% Getting the nwb into summary
% Written by D.B. Heyer
close all, clear all
clc
% Set path to load and save data
basedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/done' ;
savedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9' ;
savename = 'CellSummary_NWB2';
% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};
% Loop through abfs
index = 1 ;
errorMessages2=[];%keep a list of any errors, Anna 01-Dec
verbose = 1 %keep a list of any errors, Anna 01-Dec
for i = 1:length(filelist)
    try%try to process data for a given session, Anna 01-Dec 
 % Make subset of data per abf file
    fprintf('Looking for CC-step protocols: file nr %1.0f \n', i);
    load(fullfile(basedir,filelist{i})) ;
    stimsets = struct2table(obj.getstimsets.metadata) ;
    sweeps = struct2table(obj.getstimsets.getnwbchannel.getsweep.metadata) ;
    epochs = struct2table(obj.getstimsets.getnwbchannel.getsweep.getepoch.metadata) ;
    aps = struct2table(obj.getstimsets.getnwbchannel.getsweep.getepoch.getap.metadata) ;
    % remove incomplete sweeps
    sweeps(sweeps.nrofepochs < 3,:) = [] ;
    %% get freqbins
    edges = 1:10:201 ;
    aps.freqbin = discretize(aps.freq, edges) ;
    aps.freqbin(isnan(aps.freqbin)) = 0 ;
    aps.currinj = aps.number*0 ;
    for ii = 1:height(aps)
        aps(ii,:).currinj = epochs(ismember(epochs.guid,aps(ii,:).parent_guid),:).amplitude ;  
    end
    aps.updownratio = aps.maxdvdt./abs(aps.mindvdt) ;
    aps.onsetrapidity(aps.onsetrapidity > 100) = NaN ;

    % Create nested structure
    aps2 = table2struct(aps) ;
    epochs2 = table2struct(epochs) ;
    sweeps2 = table2struct(sweeps) ;
    stimset = table2struct(stimsets) ;
    
    for ii = 1:length(epochs2)    
            epochs2(ii).ap = aps2(ismember({aps2.parent_guid},epochs2(ii).guid)) ;    
    end
    for ii = 1:length(sweeps2)    
            sweeps2(ii).epoch = epochs2(ismember({epochs2.guid_swp},sweeps2(ii).guid)) ;    
    end
    for ii = 1:length(stimset)    
            stimset(ii).sweep = sweeps2(ismember([sweeps2.number],stimset(ii).sweepnrs)) ;         
    end
    
    % If abf is a stepprotocol: continue with analysis

        fprintf('Retrieving analysis parameters from CC-step file %1.0f \n', index);
        % Analyze
        sweep = sweeps2 ;
        NrofSweeps = length(sweep) ;  
        % find current injection epoch and assign aps to sweep
%         for step = 1:length(sweep(1).epoch)
%             if sweep(1).epoch(step).amplitude ~= 0 && seconds(sweep(1).epoch(step).timespan) > 0.03
%                 break
%             end
%         end
        for j = 1:NrofSweeps
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
            sweep(j).vmbase = sweep(j).epoch(step-1).steadystate ;
            sweep(j).jitter = sweep(j).epoch(step-1).jitter ;
            sweep(j).currinj = sweep(j).epoch(step).amplitude ;
            sweep(j).vmresponse = sweep(j).epoch(step).vstep ;
            sweep(j).ap = sweep(j).epoch(step).ap ;
            %try
                for ap = 1:length(sweep(j).epoch(step+1).ap)
                    sweep(j).rbap(ap) = sweep(j).epoch(step+1).ap(ap) ;      
                end
            %end
        end
        
        sweep=sortrows(struct2table(sweep),{'currinj','number'});
        sweep=table2struct(sweep);
        % find rheobase sweep
        apcrit=0;
        for frstspikeswp = 1:NrofSweeps
            step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
            if sweep(frstspikeswp).epoch(step).nrofaps > 0 && sweep(frstspikeswp).epoch(step).amplitude > 0
                apcrit=1;
                break
            end
        end

        idx1 = 1 ;
        for j = frstspikeswp:NrofSweeps  
            step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
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
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
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
        step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
        TrainCurr = sweep(frstspikeswp).epoch(step).amplitude +50 ;
        for j = 1:NrofSweeps
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
            tmp(j) = abs(sweep(j).epoch(step).amplitude - TrainCurr) ;
        end
        
        [TrSwp TrSwp] = min(tmp) ;
        CurrAbvRheo=NaN;
        TrSweepCrit=[];
        for TrainSweep = TrSwp:NrofSweeps  
            step = find(strcmp({sweep(TrainSweep).epoch.idxstr}, 'B'));
            if length(sweep(TrainSweep).ap) > 3
                isis = [sweep(TrainSweep).ap(2:end).isi];
                stutterAP = [0 0 0 isis(3:end) > 3*isis(2:end-1)];
                stutterISI= [0 0 isis(3:end) > 3*isis(2:end-1)];
               if length(sweep(TrainSweep).ap(~stutterAP)) > 3
                CurrAbvRheo = sweep(TrainSweep).epoch(step).amplitude - (TrainCurr-50) ;
                TrSweepCrit=1;
                break
               end
            elseif TrainSweep==NrofSweeps
                TrSweepCrit=0;
            end  
        end
       
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
            OnsetTSFAP = sweep(TrainSweep).ap(1).thresh_time - (seconds(sum([sweep(TrainSweep).epoch(1:find(strcmp({sweep(TrainSweep,1).epoch.idxstr}, 'A'))).timespan]))*1000) ;
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
        Summary(index).File               = stimset(1).filename ;
        Summary(index).Date               = obj.filetimestart ;
        Summary(index).UserID             = obj.userid ;
        Summary(index).guid               = obj.guid ;
        Summary(index).Channel            = NaN ;
        Summary(index).scalefactor        = NaN ;
        Summary(index).holdingcurrent     = NaN ;
        Summary(index).holdingvoltage     = NaN ;
        Summary(index).NrofSweeps         = NrofSweeps ;
        Summary(index).PDur               = seconds(sweep(1).epoch(strcmp({sweep(1).epoch.idxstr}, 'B')).timespan)*1000 ;
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
        Summary(index).Sag                = sweep(sagswp,1).epoch(strcmp({sweep(sagswp,1).epoch.idxstr}, 'B')).sag / PkDeflect(sagswp,1) ;
        Summary(index).VmatSag            = MinVmResponse(sagswp,1) ;
        Summary(index).TauM               = nanmean(taus(taus~=0)) ;
        Summary(index).TauSD              = nanstd(taus(taus~=0)) ;
        Summary(index).OnsetFrstAP        = sweep(frstspikeswp).ap(1).thresh_time - (seconds(sum([sweep(frstspikeswp).epoch(1:find(strcmp({sweep(frstspikeswp,1).epoch.idxstr}, 'A'))).timespan]))*1000) ; 
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
        clearvars -except Summary i basedir savedir savename filelist index verbose errorMessages2
        index = index + 1 ;
        
        catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages2=[errorMessages2; filelist{i} {ME}];%append the error message to a list
        end
    end
end
% save
save(fullfile(savedir, savename), 'Summary') ;
save('/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/errorMessages_matsummary.mat','errorMessages2');
clearvars -except Summary i
% save the file into the excel summary - AG, 1-Dec-2021
Summary_table = struct2table(Summary)
writetable(Summary_table,'Summary.xlsx')
%General end error message - Anna 1-Dec
%if ~isempty(errorMessages2)
%    for k=1:size(errorMessages2,1)
    %    fprintf('\nError at session %d\n',errorMessages2{k,1});%print problematic session number to screen
     %   errorMessages2{k,2}.message%display error message
      %  errorMessages2{k,2}.stack%print function and line number at which error occurred
    %end
%end

%% sag_ratio analysis
clear all
clc
basedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/done' ;
savedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/done' ;
savename = 'Summary_sag_ratio' ;

% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};

% settings
filter_freq = 1000; % low pass filter frequency for sag traces (Hz)

% loop over files
%t_all=table();
table_initialize = false ;
errorMessages3=[];%keep a list of any errors, Anna 01-Dec-2021
verbose=1%keep a list of any errors, Anna 01-Dec-2021
for i=1:numel(filelist)
    try%try to process data for a given session, Anna 01-Dec-2021
    %load file
    fn = fullfile(basedir, filelist{i});
    nwb = load(fn);
    nwb = nwb.obj;
    
    % loop over stimsets
    for j = 1:nwb.nrofstimsets
        stimset = nwb.getstimset(j);
%         stimset_filter = nwb.getstimset.name
%         ({'Sub'} {'CCSteps_DA_0'});
        swps = stimset.getnwbchannel.getsweep();
        
        % extract from each sweep:
        % step size (pA)
        stepepoch = find([swps(1).getepoch.amplitude]~=0, 1);
        pA = [swps.getepoch(stepepoch).amplitude];
        
        %select only negative injection swps
        if all(pA>=0), continue;end
        swps = swps(pA<0);
        pA = pA(pA<0);
        
    
      %  swptime = swps.Time ;
      %  if all(swptime<300000), continue;end
        
        for k = 1:numel(swps)
            loc = find(isnan(swps(k).Data), 1); %find from where swp is NaN
            if ~isempty(loc) ,swps(k) = swps(k).getsamples(1:loc-1);end
        end
        
        % peak and ss voltage
        
%         swps_f = swps.low_pass_filter(1000);
        
        vstep = [swps.getepoch(stepepoch).vstep];
        rmp = [swps.getepoch(stepepoch-1).steadystate];
        v_delta = vstep - rmp;
        ss = [swps.getepoch(stepepoch).steadystate];
        ss_delta = ss - rmp;
        delta = ss - vstep;
        sag_ratio = delta ./ -v_delta;
        
        % peak time
        
        % time constant of sag activation
        
        % rebound amplitude
%         reb_start_time = swps(1).getepoch(stepepoch).Time(end);
%         reb_ss_time = swps(1).Time(end);
%         rebound_window=300; % in ms
%         reb_traces = swps.getsampleusingtime(reb_start_time, reb_start_time+rebound_window);
%         reb_volts = NaN(1, numel(reb_traces));
%         reb_ss=reb_volts;
%         for k = 1:numel(reb_traces)
%             reb_volts(k)=max(reb_traces(k).medianfilter(0.5,'truncate').getdata);
%             reb_ss(k) = swps(k).getsampleusingtime(reb_ss_time-200, reb_ss_time).median;
%         end
%         reb_amps = reb_volts-reb_ss;
        
        % results to table
     % Protocol_column = repmat(stimset.name, height(t), 1);
     % Filename_column = repmat(filelist{i}, height(t), 1); 
      
        t = table(pA', vstep', rmp', ss', ss_delta', delta', sag_ratio', 'VariableNames', {'pA', 'sagvolt', 'RMP', 'ss', 'ss_delta', 'delta', 'sag_ratio'});
      %  t.Protocol = repmat(stimset.name, height(t), 1);
      t.Filename = repmat({filelist{i}}, height(t), 1);
      t.Protocol = repmat({stimset.name},height(t),1) %add a ; if you dont want it printed in command
        
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
        %Anna's edition 2021-Dec-1
        writetable(t_all,'sag_ratios.xlsx');%if you add an extra argument "Sheet",
        %1/2/etc - then you can add multiplie MAtLab tables into one xlsx
    end
    catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end        
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages3=[errorMessages3; filelist{i} {ME}];%append the error message to a list
        end
    end
end
save('/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/errorMessages_matsagratio.mat','errorMessages3');
%General end error message - Anna 1-Dec
if ~isempty(errorMessages3)
    for k=1:size(errorMessages3,1)
        fprintf('\nError at session %d\n',errorMessages3{k,1});%print problematic session number to screen
        errorMessages3{k,2}.message%display error message
        errorMessages3{k,2}.stack%print function and line number at which error occurred
    end
end


%% AP upstroke analysis per cell
clear all
clc

basedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9' ;
savedir = '/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/done' ;
fileinfo  = dir(fullfile(basedir,'*.nwb'));
filelist  = {fileinfo.name};
table_initialize = false ; % Anna 2-Dec-2021
errorMessages4=[];%keep a list of any errors, Anna 01-Dec
verbose=1%keep a list of any errors, Anna 01-Dec
for i=1:numel(filelist)
    try%try to process data for a given session
    fn = fullfile(basedir, filelist{i});
    nwb = NWBfile(fn, {'X2LP_Search','teps'});
    %getting the sweeps I need
    swps = nwb.getstimset.getnwbchannel.getsweep;
    %analysis of the sweeps
    swps = swps.analysesweep;
    % getting th eAPs from the sweeps which have them
    APs = swps.getap;
    %%getting the frequency of every AP
    inst_freqs = [APs.freq];
    % dividing the APS into the bins to get the ones we need
    binedges = [0,0.00001, 10, 20, 30, 40, 500];
    bins = discretize(inst_freqs, binedges);
    binlabels = {'FirstAP', 's0to10 Hz', 's10to20 Hz','s20to30 Hz' 's30to40 Hz', 's40toINF Hz'};
   
    upstrokes = [APs.maxdvdt];
    mean_binned_upstrokes=NaN(1, numel(binlabels));
    median_binned_upstrokes=NaN(1, numel(binlabels));
    for j=1:numel(binlabels)
        mean_binned_upstrokes(j) = nanmean(upstrokes(bins==j) )  ;
        median_binned_upstrokes(j) = nanmedian(upstrokes(bins==j) )  ;
    end
    %boxplot(upstrokes, bins)
    %hold on
    %plot(swps)
    %scatter(categorical(binlabels(bins)), upstrokes)
    Mean_upstroke = array2table(mean_binned_upstrokes, 'VariableNames', {'M_FirstAP', 'M_s0to10 Hz', 'M_s10to20 Hz','M_s20to30 Hz' 'M_s30to40 Hz', 'M_s40toINF Hz'});
    Median_upstroke = array2table(median_binned_upstrokes, 'VariableNames', {'Md_FirstAP', 'Md_s0to10 Hz', 'Md_s10to20 Hz','Md_s20to30 Hz' 'Md_s30to40 Hz', 'Md_s40toINF Hz'});
    Summary_per_cell = [Mean_upstroke Median_upstroke];
    Summary_per_cell.Filename = repmat({filelist{i}}, height(Summary_per_cell), 1);
    %append results to overview
        if table_initialize 
            t_all = [t_all; Summary_per_cell];
        else 
            t_all = Summary_per_cell ; 
            table_initialize = true
        end
        %Anna_edition 1-Dec - handling errors
    catch ME%catch any errors
        if verbose
            fprintf('Error at session %s\n',filelist{i});%print session number to screen
        end        
        disp(ME);%print the error message to screen
        for k=1:length(ME.stack)
            errorMessages3=[errorMessages4; filelist{i} {ME}];%append the error message to a list
        end
    end  
end  

writetable(t_all,'APs.xlsx')

%General end error message - Anna 1-Dec
if ~isempty(errorMessages4)
    for k=1:size(errorMessages4,1)
        fprintf('\nError at session %d\n',errorMessages4{k,1});%print problematic session number to screen
        errorMessages4{k,2}.message%display error message
        errorMessages4{k,2}.stack%print function and line number at which error occurred
    end
end
save('/Users/annagalakhova/PhD INF CNCR VU/DATA/OSC/Human/H22.29.209/DIV9/errorMessages_nwbAPupstroke.mat','errorMessages4');




%% 


