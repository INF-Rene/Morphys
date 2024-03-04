%% Analysis script
close all, clear all
%%
folder = uigetdir 
cd (folder);
%make a list with all nwbsephys
list = dir();
list = struct2table(list);
list = list(list.bytes>10000,:); %only files with actual data

%% USE THIS nwb2: 
for i = 1:numel(list.name) 
    fn =cell2mat(list.name(i));
 nwb = NWBfile(fn,[{'LP'} {'hresh'} {'CC'} {'teps'} {'LSFINEST'} {'LSCOARSE'}]);
 obj =nwb.analyseNWB ;
 obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
 saveme(obj,'/Volumes/Expansion/Ephys/analyzed/MAT3', obj.savename) 
end
%% Set path to load and save data; mat data load 
basedir = '/Volumes/Expansion/Ephys/analyzed/daan test' ;
savedir = '/Volumes/Expansion/2024 summary' ;
savename = 'Summary_nwb4' ;
%load('/Users/elinemertens/Data/ephys/Masterstruct/dataStruct4.mat');
% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};
%% Loop through abfs
index = 1 ;  
% Initialize empty tables for storing  data
% dataStruct.ApsTime = {};
% dataStruct.ApsData = {};
% dataStruct.Derivatives = {};
for i = 1:length(filelist)
    %% Make subset of data per abf file
    fprintf('Looking for CC-step protocols: file nr %1.0f \n', i);
    load(fullfile(basedir,filelist{i})) ;
      % either run through all protocols, if it can't resolve this, remove %    
%      stimnms={obj.getstimsets.name};
%      CCsteploc=cellfun(@(x) contains(x, 'teps'), stimnms);
%    stimsets = struct2table(obj.getstimsets(CCsteploc).metadata, 'AsArray', true) ;
    stimsets = struct2table(obj.getstimsets.metadata) ;
    sweeps = struct2table(obj.getstimsets.getnwbchannel.getsweep.metadata) ;
    epochs = struct2table(obj.getstimsets.getnwbchannel.getsweep.getepoch.metadata) ;
    aps = struct2table(obj.getstimsets.getnwbchannel.getsweep.getepoch.getap.metadata) ;
    % remove incomplete sweeps
    sweeps(sweeps.nrofepochs < 3,:) = [] ;
    %% get freqbins
    edges = 1:20:201 ;
    aps.freqbin = discretize(aps.freq, edges) ;
    aps.freqbin(isnan(aps.freqbin)) = 0 ;
    aps.currinj = aps.number*0 ;
    for ii = 1:height(aps)
        aps(ii,:).currinj = epochs(ismember(epochs.guid,aps(ii,:).parent_guid),:).amplitude ;  
    end
    aps.updownratio = aps.maxdvdt./abs(aps.mindvdt) ;
    %aps.onsetrapidity(aps.onsetrapidity > 100) = NaN ;
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
    %% If abf is a stepprotocol: continue with analysis
        fprintf('Retrieving analysis parameters from CC-step file %1.0f \n', index);
        %% Analyze
        sweep = sweeps2 ;
        NrofSweeps = length(sweep) ; 
        NrofEpochs = length(epochs2) ;
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
            try
                for ap = 1:length(sweep(j).epoch(step+1).ap)
                    sweep(j).rbap(ap) = sweep(j).epoch(step+1).ap(ap) ;      
                end
            end
        end        
        sweep=sortrows(struct2table(sweep),{'currinj','number'});
        sweep=table2struct(sweep);
        % find rheobase sweep
        apcrit=0;
        for frstspikeswp = 1:NrofSweeps
            step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
            if sweep(frstspikeswp).epoch(step).nrofaps > 0 && sweep(frstspikeswp).epoch(step).amplitude > 0
                apcrit=1;
                firstsweepname = sweep(frstspikeswp).Name ;
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
        vmbase = [sweep.vmbase] ;
        Freqs=[];
        StimInts=[];
        currInjections_R=[];
        voltageResponses=[];
        sags = [] ;
        taus=[];
         for j = 1:NrofSweeps       
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
            if sweep(j,1).currinj >= -100 && sweep(j,1).currinj % < 0 && sweep(j,1).vmbase < -60               
                % remove the % 1 row above when there are rinputs that aren negative (means that the cell was unhealthy and got depolarized
                %throughout the subthreshold, by removing that vmbase < -60 you correct for that unsteady rmp
                voltageResponses(j,1) = sweep(j,1).vmresponse  ; 
                %sags(j,1) = sweep(j,1).sag ; 
                % for cells with a non fixed (-70mV) rmp, use the
                % difference between vmresponse and the vmbase: 
                % voltageResponses(j,1) = sweep(j,1).vmresponse - sweep(j,1).vmbase
                currInjections_R(j,1) = sweep(j,1).currinj ;
                if sweep(j,1).epoch(step).tau < 100 && sweep(j,1).epoch(step).tau > 0 && sweep(j,1).epoch(step).gof > 0.95
                    taus(j,1) = sweep(j,1).epoch(step).tau ;
                else
                    taus(j,1) = NaN;
                end
                if sweep(j,1).epoch(step).amplitude < 0
                    sags(j,1) = sweep(j,1).epoch(step).sag ; 
                    sagvoltage(j,1) = sweep(j,1).epoch(step).steadystate_diff ; 
                    sagratios(j,1) = -1 * (sags(j,1) / sagvoltage(j,1)) ;
                else 
                    sags(j,1) = NaN ;
                    sagvoltage(j,1) = NaN ;
                    sagratios(j,1) = NaN ;
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

            if length(sweep(j).ap) >= 2 && sweep(j,1).currinj > 0
                Freqs(j,1) =  mean([sweep(j).ap(2:end).freq]) ;
                StimInts(j,1) = [sweep(j,1).currinj] ;
            end         
        end

        if sum(NrofRBAPs) > 0
            NrofRBAPsM = mean(NrofRBAPs(NrofRBAPs~=0)) ;
        else
            NrofRBAPsM = 0 ;
        end
%         
%          for j = 1:NrofEpochs
%              step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
%              if epochs2(step).amplitude < -25
%                  sags(j) = epochs2(step).sag ;
%                  sagvoltage(j) = epochs2(step).steadystate_diff ; 
%                  sagratios(j) = sags / sagvoltage ; 
%                 else
%                     sags(j) = NaN;
%                     sagvoltage(j) = NaN ;
%                     sagratios(j) = NaN ;                   
%              end
%          end 

      % Find trainsweep 
step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
Rheobase = sweep(frstspikeswp).currinj; 
TrainCurr = Rheobase * 1.5;

for j = 1:NrofSweeps
    step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
    tmp(j) = abs(sweep(j).epoch(step).amplitude - TrainCurr);
end       

[TrSwp TrSwp] = min(tmp);
TrSweepCrit = 0;

for TrainSweep = TrSwp:NrofSweeps  
    step = find(strcmp({sweep(TrainSweep).epoch.idxstr}, 'B'));
    CurrAbvRheo = sweep(TrainSweep).epoch(step).amplitude - Rheobase;
    
    if length(sweep(TrainSweep).ap) <= 2                
        TrSweepCrit = 1;
        %stutterAP = NaN ; 
        %stutterISI = NaN;
    elseif length(sweep(TrainSweep).ap) >= 3
       % stutterAP = [0 0 0 isis(3:end) > 3*isis(2:end-1)];
       % stutterISI = [0 0 isis(3:end) > 3*isis(2:end-1)];
        TrSweepCrit = 1;
        break;
    elseif TrainSweep == NrofSweeps
        TrSweepCrit = 3;
    end  
end 


if ~isempty(sweep(TrainSweep).ap)
    for l = 1:length(sweep(TrainSweep).ap)           
        if ~isempty(sweep(TrainSweep).ap(l,1).halfwidth)
            HWsTS(l,1) = sweep(TrainSweep).ap(l,1).halfwidth;  
        end       
    end
else
    HWsTS = NaN;
end

if TrSweepCrit == 1
    TSbasetothresh = ([sweep(TrainSweep).ap.thresh] - sweep(TrainSweep).vmbase);
    TSpeaktoahp = ([sweep(TrainSweep).ap.ahp_time] - [sweep(TrainSweep).ap.peak_time]); 
    AmpsTSthresh = [sweep(TrainSweep).ap.amp];
    AHPsTS = [sweep(TrainSweep).ap.relahp];
    AHPslowTS = [sweep(TrainSweep).ap.relahp_slow];
    ISIsTS = [sweep(TrainSweep).ap(2:end).isi];
    FreqTrSwp = mean([sweep(TrainSweep).ap(1:end).freq]);
    NrOfAPsTrSwp = length(sweep(TrainSweep).ap); 
    OnsetTSFAP = sweep(TrainSweep).ap(1).thresh_time - ...
        (seconds(sum([sweep(TrainSweep).epoch(1:find(strcmp({sweep(TrainSweep,1).epoch.idxstr}, 'A'))).timespan])) * 1000);
    onsetrapidity = sweep(TrainSweep).ap(1).onsetrapidity ; 
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
    onsetrapidity = NaN ; 
end


 TrainCurr2 = Rheobase * 2 ;
        for j = 1:NrofSweeps
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
            tmp(j) = abs(sweep(j).epoch(step).amplitude - TrainCurr2) ;
        end
        
        [TrSwp2 TrSwp2] = min(tmp) ;
        CurrAbvRheo2=NaN;
        for TrainSweep2 = TrSwp2:NrofSweeps  
            step = find(strcmp({sweep(TrainSweep2).epoch.idxstr}, 'B'));
            if length(sweep(TrainSweep2).ap) > 1
                TrSweepCrit2=1;
                CurrAbvRheo2 = sweep(TrainSweep2).epoch(step).amplitude - Rheobase ;               
                break
               end
            if TrainSweep2==NrofSweeps
                TrSweepCrit2=0;
            end  
        end      

  for l = 1:length(sweep(TrainSweep2))           
    if ~isempty(sweep(TrainSweep2).ap(l,1))
            HWsTS2(l,1) = [sweep(TrainSweep2).ap(l,1).halfwidth] ;
            TrSweepCrit2 = 0 ; 
            TSbasetothresh2 = NaN;
            TSpeaktoahp2(l,1) = NaN;
            AmpsTSthresh2 = NaN;
            AHPsTS2 = NaN;
            AHPslowTS2 = NaN;
            ISIsTS2 = NaN;
            FreqTrSwp2 = NaN;
            NrOfAPsTrSwp2 = NaN;
            OnsetTSFAP2 = NaN;
            HWsTS2=NaN;
    end
  end
        
%             TSbasetothresh2 = NaN;
%             TSpeaktoahp2 = NaN;
%             AmpsTSthresh2 = NaN;
%             AHPsTS2 = NaN;
%             AHPslowTS2 = NaN;
%             ISIsTS2 = NaN;
%             FreqTrSwp2 = NaN;
%             NrOfAPsTrSwp2 = NaN;
%             OnsetTSFAP2 = NaN;
%             HWsTS2=NaN;
% %             for l = 1:length(sweep(TrainSweep2).ap)           
% %                 if ~isempty(sweep(TrainSweep2).ap(l,1).halfwidth)
% %                     HWsTS2(l,1) = [sweep(TrainSweep2).ap(l,1).halfwidth] ;  
% %                 end       
% %             end
%         end
        
        if TrSweepCrit2==1
            TSbasetothresh2 = ([sweep(TrainSweep2).ap.thresh]-sweep(TrainSweep2).vmbase) ;
            TSpeaktoahp2 = ([sweep(TrainSweep2).ap.ahp_time]-[sweep(TrainSweep2).ap.peak_time]); 
            AmpsTSthresh2 = [sweep(TrainSweep2).ap.amp] ;
            AHPsTS2 = [sweep(TrainSweep2).ap.relahp] ;
            AHPslowTS2 = [sweep(TrainSweep2).ap.relahp_slow] ;
            ISIsTS2 = [sweep(TrainSweep2).ap(2:end).isi] ;
            FreqTrSwp2 = mean([sweep(TrainSweep2).ap(1:end).freq]) ;
            NrOfAPsTrSwp2 = length(sweep(TrainSweep2).ap) ; 
            OnsetTSFAP2 = sweep(TrainSweep2).ap(1).thresh_time - (seconds(sum([sweep(TrainSweep2).epoch(1:find(strcmp({sweep(TrainSweep2,1).epoch.idxstr}, 'A'))).timespan]))*1000) ;
            isis_TS2 = [sweep(TrainSweep2).ap(2:end).isi];
            isis_TS2_1 = [sweep(TrainSweep2).ap(2).isi];
            NrofAPtrainSwp2 = length(sweep(TrainSweep2).ap)
            ActualTCurr2 = sweep(TrainSweep2).currinj
            onsetrapidity2 = sweep(TrainSweep2).ap(1).onsetrapidity ; 
            ThreshTSAP2 = sweep(TrainSweep2).ap(1).thresh ; 
            TSAP2basetothresh = sweep(TrainSweep2).ap(1).thresh-sweep(TrainSweep2).vmbase ;
            AmpTSAP2thresh = sweep(TrainSweep2).ap(1).amp;
            TSAP2peaktoahp = sweep(TrainSweep2).ap(1).ahp_time - sweep(TrainSweep2).ap(1).peak_time;
            HalfWTSAP2 = sweep(TrainSweep2).ap(1).halfwidth;
            AHPTSAP2 = sweep(TrainSweep2).ap(1).relahp;
            AHPslowTSAP2 = sweep(TrainSweep2).ap(1).relahp_slow;
            UpStrokeTSAP2 = sweep(TrainSweep2).ap(1).upstroke;
            DwnStrokeTSAP2 = sweep(TrainSweep2).ap(1).downstroke ;
            UpDwnStrkRatioTSAP = abs(sweep(TrainSweep2).ap(1).upstroke) / abs(sweep(TrainSweep2).ap(1).downstroke);
            MaxUpStrkTSAP2 = sweep(TrainSweep2).ap(1).maxdvdt;
            MaxDwnStrkTSAP2 = sweep(TrainSweep2).ap(1).mindvdt;
                       
        else
            TSbasetothresh2 = NaN;
            TSpeaktoahp2 = NaN;
            AmpsTSthresh2 = NaN;
            AHPsTS2 = NaN;
            AHPslowTS2 = NaN;
            ISIsTS2 = NaN;
            FreqTrSwp2 = NaN;
            NrOfAPsTrSwp2 = NaN;
            OnsetTSFAP2 = NaN;
            isis_TS2 = NaN ;
            isis_TS2_1 = NaN  ;
            NrofAPtrainSwp2 = NaN ; 
            ActualTCurr2 = NaN ; 
            onsetrapidity2 = NaN ;
            ThreshTSAP2       = NaN ;
            TSAP2basetothresh = NaN;
             AmpTSAP2thresh = NaN ;
             TSAP2peaktoahp = NaN ;
             HalfWTSAP2 = NaN ;
             AHPTSAP2 = NaN ;
             AHPslowTSAP2 = NaN ;
             UpStrokeTSAP2 = NaN ;
             DwnStrokeTSAP2 = NaN ;
             UpDwnStrkRatioTSAP2 = NaN ;
             MaxUpStrkTSAP2 = NaN ;
             MaxDwnStrkTSAP2 = NaN ;
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
        
        % determine sweep for sag, (different sweeps, with median)        
  tmp = abs(MinVmResponse+100) ;
        [sagswp sagswp] = min(tmp) ;
        tmp = abs(MinVmResponse+90) ; 
        [sagswp2 sagswp2] = min(tmp) ;
        tmp = abs(MinVmResponse+80) ;
        [sagswp3 sagswp3] = min(tmp) ;                  
        step = find(strcmp({sweep(sagswp).epoch.idxstr}, 'B'));            
        Sag                = sweep(sagswp,1).epoch(step).sag / PkDeflect(sagswp,1) ;
        Sag2               = sweep(sagswp2,1).epoch(step).sag / PkDeflect(sagswp2,1) ;
        Sag3               = sweep(sagswp3,1).epoch(step).sag / PkDeflect(sagswp3,1) ;
        SagMedian             = nanmedian([Sag Sag2 Sag3]) ; 
        
        
         
        % calculate input frequency curve
        Freqs = Freqs(Freqs~=0) ;
        StimInts = StimInts(StimInts~=0) ;
        if length(Freqs) > 1
            [fitFi]=fit(StimInts,Freqs,f_R, 'StartPoint', [0 0]); 
            FrqChngStimInt = fitFi.R ;
        else  
            FrqChngStimInt = NaN ;   
        end

       if TrSweepCrit==1
            isis=ISIsTS; %exclude stuttering ISIs since that can mess with adaptation analysis
            amps=AmpsTSthresh;
            hws=HWsTS;
        else
            isis=NaN;
            amps=NaN;
            hws=NaN;
       end
       
       if TrSweepCrit2==1
            isis2=ISIsTS; %exclude stuttering ISIs since that can mess with adaptation analysis
            amps2=AmpsTSthresh;
            hws2=HWsTS;
        else
            isis2=NaN;
            amps2=NaN;
            hws2=NaN;
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
        
        if length(isis2) > 2
            ISIRatio1toAll2 = mean(isis2(2:end)) / mean(isis2(1)) ;
            N = length(isis2)-1 ;
            for n = 1:N
                ISIchanges2(n,1) = (isis2(n+1)-isis2(n)) / (isis2(n+1)+isis2(n));
            end
            AdaptIdx2 = (1/N)*sum(ISIchanges2) ;        
        else
            ISIRatio1toAll2 = NaN;
            AdaptIdx2 = NaN;
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
        
         if TrSweepCrit2==1
         end
        
        if length(amps2) > 2           
            N = length(amps2)-1 ;
            for n = 1:N
                Ampchanges2(n,1) = (amps2(n+1)-amps2(n)) / (amps2(n+1)+amps2(n));
                HWchanges2(n,1) = (hws2(n+1)-hws2(n)) / (hws2(n+1)+hws2(n));
            end
            AmpAccom2 = (1/N)*sum(Ampchanges2) ;  
            HWAccom2 = (1/N)*sum(HWchanges2) ;
        else
            AmpAccom2 = NaN;    
            HWAccom2 = NaN; 
        end     
        
        if ~isempty(Freqs)
            freqmax = max(Freqs) ;
        else
            freqmax = NaN ;
        end        
        
         % Firstsweep and lastsweep ISIS EM edits
          % nr of APs first sweep (Eline Edit)
         NrOfAPfrstSwp = length(sweep(frstspikeswp).ap) ;
         NrofAPtrainSwp = length(sweep(TrainSweep).ap);
         NrofAPlastSwp = length(sweep(NrofSweeps).ap);
   
        if length(sweep(frstspikeswp).ap) > 1
            isis_FS = [sweep(frstspikeswp).ap(2:end).isi];
            isis_FS1 = [sweep(frstspikeswp).ap(2).isi] ; 
        else 
            isis_FS = NaN ;
            isis_FS1 = NaN ;
        end

         if length(sweep(TrainSweep).ap) > 1
            isis_TS = [sweep(TrainSweep).ap(2:end).isi];
            isis_TS1 = [sweep(TrainSweep).ap(2).isi];
          %  isis_TSend = [sweep(TrainSweep).ap(end).isi];
        else 
            isis_TS = NaN ;
            isis_TS1 = NaN  ;
         end
         
%         if length(sweep(TrainSweep2).ap) > 1
%             isis_TS2 = [sweep(TrainSweep2).ap(2:end).isi];
%             isis_TS2_1 = [sweep(TrainSweep2).ap(2).isi];
%           %  isis_TSend = [sweep(TrainSweep).ap(end).isi];
%         else 
%             isis_TS2 = NaN ;
%             isis_TS2_1 = NaN  ;
%          end
         
        if length(sweep(NrofSweeps).ap) > 1
            isis_LS = [sweep(NrofSweeps).ap(2:end).isi];
            isis_LS1 = [sweep(NrofSweeps).ap(2).isi];
        else 
            isis_LS = NaN ;
            isis_LS1 = NaN  ;
        end
     
        % tihs works for just getting the derivatives, but if you want to
        % plot them together, they need to be downsampled
%         for j = 1:length(obj.stimsets)
%             if any(ismember({obj.getstimset(j).getnwbchannel.getsweep.Name},firstsweepname))
%                 obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).plotanalysis;
%                 ApsData = obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).ts.Data ;
%                 ApsTime = obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).ts.Time ; 
%                  derivativestime = obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).getdvdt();        
%                legend(obj.filename) 
%                  ApsData = ApsData(1:2700);
%                  ApsTime = ApsTime(1:2700);
%                  derivativestime = derivatives(1:2700);
%                    ApsTimeArray = [ApsTimeArray; ApsTime'];
%                  ApsDataArray = [ApsDataArray; ApsData'];
%                     DerivativesArray = [DerivativesArray ; derivativestime'];
%              %   title('First AP')
%             %   ylabel('mV')
%             %    xlabel('ms')
%             end
%         end
%   
for j = 1:length(obj.stimsets)
    if any(ismember({obj.getstimset(j).getnwbchannel.getsweep.Name}, firstsweepname))
        % Retrieve data for the first action potential
        ap = obj.getstimset(j).getnwbchannel.getsweep('Name', firstsweepname).aps(1) ;
        obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).plotanalysis;
        sweeptest = obj.getstimset(j).getnwbchannel.getsweep('Name', firstsweepname) ;
        derivativestime = obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).aps(1).getdvdt();        
        samplefreq = sweeptest.samplefreq ; 
        desired_sampling_rate = 25000;
        Q = round(samplefreq / desired_sampling_rate);
        P = 1 ; 
        %Q = round(resample_factor) ; 
        ApsData = ap.ts.Data;
        ApsTime = ap.ts.Time;
        ApsData = double(ApsData); % Convert ApsData to double
        ApsDer = double(derivativestime) ; 
        ApsDataR = resample(ApsData, P, Q) ; 
        ApsTimeR = resample(ApsTime, P, Q) ;
        DerivativesR = resample(ApsDer, P, Q) ; 
        %AP_resampled = resample(ApsData, 0:1/resample_factor:(length(ApsData)-1)/samplefreq);
         % Append data to tables with column names from filelist
        % Append resampled data to struct
          % Preprocess file name to make it a valid field name
        fileName = strrep(obj.filename, '-', '_'); % Replace '-' with '_'
        fileName = strrep(fileName, '.', ''); % Remove '.'
        fileName = matlab.lang.makeValidName(fileName); % Make valid MATLAB identifier        
        % Add resampled data to struct with modified file name as field name
        dataStruct.ApsTime.(fileName) = ApsTimeR';
        dataStruct.ApsData.(fileName) = ApsDataR';
        dataStruct.Derivatives.(fileName) = DerivativesR';
    end
end
%         % Check if the field exists
% if ~isempty(sweep(TrainSweep).ap) && ~isempty(sweep(TrainSweep).ap(1).onsetrapidity)
%     onsetrapidityTSAP = sweep(TrainSweep).ap(1).onsetrapidity;
% else
%    onsetrapidityTSAP = NaN;
% end
        
        %get the holding current 
        lb = obj.getstimset(1).getnwbchannel(1).labbooknum ; 
        vals=lb.I0x2DClampHoldingLevel(~isnan(lb.I0x2DClampHoldingLevel)) ; 
        currentinj_avg = nanmedian(vals) ; 
        
%end if you want to make a large summary, remove the end here 

        %% Create summary  
        Summary(index).File               = stimset(1).filename ;
        Summary(index).UserID             = obj.userid ;
        Summary(index).guid               = obj.guid ;
        Summary(index).scalefactor        = NaN ;
        Summary(index).holdingcurrent     = NaN ;
        Summary(index).currentinj         = currentinj_avg ;
        Summary(index).NrofSweeps         = NrofSweeps ;
        Summary(index).samplefreq         = samplefreq ; 
        Summary(index).PDur               = seconds(sweep(1).epoch(strcmp({sweep(1).epoch.idxstr}, 'B')).timespan)*1000 ;
        Summary(index).FrstP              = sweep(1).currinj ;
        Summary(index).DeltaP             = sweep(2).currinj - sweep(1).currinj ;
        Summary(index).vmbaseM            = nanmean(vmbase) ;
        Summary(index).vmbaseSD           = nanstd(vmbase) ;
        Summary(index).Jitter             = nanmean([sweep.jitter]) ;
        Summary(index).NrofRBAPs          = sum(NrofRBAPs) ;
        Summary(index).NrofRBAPsM         = NrofRBAPsM ;
        Summary(index).TauM               = nanmean(taus(taus~=0)) ;
        Summary(index).TauSD              = nanstd(taus(taus~=0)) ;
        Summary(index).InputR             = Rin ;% in MOhm...
        Summary(index).Sag1               = Sag ;
        Summary(index).Sag2               = Sag2 ; 
        Summary(index).Sag3               = Sag3 ; 
        Summary(index).SagMedian          = SagMedian ; 
        Summary(index).SagMedianAll         = nanmedian(sagratios(sagratios~=0)) ;
        Summary(index).VmatSag1            = MinVmResponse(sagswp,1) ;
        Summary(index).VmatSag2           = MinVmResponse(sagswp2,1) ;
        Summary(index).VmatSag3           = MinVmResponse(sagswp3,1) ;
        Summary(index).FrstSpikeSwp       = frstspikeswp ;
        Summary(index).TrainSwp          = TrainSweep ;
        Summary(index).TrainSwp2          = TrainSweep2 ;
        Summary(index).NrofAPsFrstSwp     = NrOfAPfrstSwp ; 
        Summary(index).NrOfAPsTrSwp       = NrOfAPsTrSwp ;
        Summary(index).NrofAPtrainSwp2    = NrofAPtrainSwp2 ;
        Summary(index).NrofAPlastSwp      = NrofAPlastSwp ;
        Summary(index).Rheobase           = sweep(frstspikeswp).currinj ;
        Summary(index).TCurrCalculation   = TrainCurr ; 
        Summary(index).TCurrCalculation2   = TrainCurr2 ; 
        Summary(index).ActualTCurr         = sweep(TrainSweep).currinj ;
        Summary(index).ActualTCurr2         = ActualTCurr2 ;
        Summary(index).TCurrAbvRheo       = CurrAbvRheo ;
        Summary(index).CurrAbvRheo2        = CurrAbvRheo2 ;
        %Summary(index).isis_FS            = isis_FS ;
        Summary(index).isis_FS1           = isis_FS1 ;
       % Summary(index).isis_TS            = isis_TS ;
        Summary(index).isis_TS1           = isis_TS1 ;
        Summary(index).isis_TS2           = isis_TS2_1 ;
       % Summary(index).isis_LS            = isis_LS ;
        Summary(index).isis_LS1           = isis_LS1 ;
        Summary(index).FreqMax            = freqmax ;
        Summary(index).NrOfAPs            = sum(NrofAPs) ;
        Summary(index).NrOfAPsMax         = max(NrofAPs) ;
        Summary(index).FrqChngStimInt     = FrqChngStimInt ;        
        Summary(index).starttime          = sweep(frstspikeswp).ap(1).start_time ; 
        Summary(index).endtime            = sweep(frstspikeswp).ap(1).end_time ;
        Summary(index).peaktime           = sweep(frstspikeswp).ap(1).peak_time ; 
        Summary(index).OnsetFrstAP        = sweep(frstspikeswp).ap(1).thresh_time - (seconds(sum([sweep(frstspikeswp).epoch(1:find(strcmp({sweep(frstspikeswp,1).epoch.idxstr}, 'A'))).timespan]))*1000) ; 
        Summary(index).onsetrapidity      = sweep(frstspikeswp).ap(1).onsetrapidity ; 
        Summary(index).ThreshFrstAP       = sweep(frstspikeswp).ap(1).thresh ; 
        Summary(index).FAPbasetothresh    = sweep(frstspikeswp).ap(1).thresh-sweep(frstspikeswp).vmbase ; 
        Summary(index).AmpFAPthresh       = sweep(frstspikeswp).ap(1).amp ;
        Summary(index).FAPpeaktoahp       = sweep(frstspikeswp).ap(1).ahp_time - sweep(frstspikeswp).ap(1).peak_time ;
        Summary(index).HalfWFrstAP        = sweep(frstspikeswp).ap(1).halfwidth ; 
        Summary(index).AHPFrstAP          = sweep(frstspikeswp).ap(1).relahp ;
        Summary(index).AHPslowFrstAP      = sweep(frstspikeswp).ap(1).relahp_slow ;
        Summary(index).UpStrkFrstAP       = sweep(frstspikeswp).ap(1).upstroke ;
        Summary(index).DwnStrkFrstAP      = sweep(frstspikeswp).ap(1).downstroke ;
        Summary(index).UpDwnStrkRatio     = abs(sweep(frstspikeswp).ap(1).upstroke) / abs(sweep(frstspikeswp).ap(1).downstroke) ;
        Summary(index).MaxUpFrstAP        = sweep(frstspikeswp).ap(1).maxdvdt ;
        Summary(index).MaxDwnFrstAP       = sweep(frstspikeswp).ap(1).mindvdt ;
        Summary(index).TrainSwp           = TrainSweep ; 
        Summary(index).FreqTrSwp          = FreqTrSwp ;
         Summary(index).AmpAccom           = AmpAccom ; 
        Summary(index).HWAccom            = HWAccom ;  
        Summary(index).ISIRatio1toAll     = ISIRatio1toAll ;
        Summary(index).AdaptIndexTS       = AdaptIdx ;
        Summary(index).OnsetTSFAP         = OnsetTSFAP ; 
        Summary(index).onsetrapidityTSAP  = sweep(TrainSweep).ap(1).onsetrapidity ; 
        Summary(index).ThreshTSAP         = sweep(TrainSweep).ap(1).thresh ; 
        Summary(index).TSAPbasetothresh    = sweep(TrainSweep).ap(1).thresh-sweep(TrainSweep).vmbase ; 
        Summary(index).AmpTSAPthresh       = sweep(TrainSweep).ap(1).amp ;
        Summary(index).HalfWTSAP        = sweep(TrainSweep).ap(1).halfwidth ; 
        Summary(index).AHPTSAP          = nanmean(sweep(TrainSweep).ap(1).relahp) ;
        Summary(index).AHPslowTSAP      = nanmean(sweep(TrainSweep).ap(1).relahp_slow) ;
        Summary(index).UpStrokeTSAP       = sweep(TrainSweep).ap(1).upstroke ;
        Summary(index).DwnStrokeTSAP      = sweep(TrainSweep).ap(1).downstroke ;
        Summary(index).UpDwnStrkRatioTSAP     = abs(sweep(TrainSweep).ap(1).upstroke) / abs(sweep(TrainSweep).ap(1).downstroke) ;
        Summary(index).MaxUpStrkTSAP       = sweep(TrainSweep).ap(1).maxdvdt ;
        Summary(index).MaxDwnStrkTSAP      = sweep(TrainSweep).ap(1).mindvdt ;
        Summary(index).TrainSwp2           = TrainSweep2 ; 
        Summary(index).FreqTrSwp2          = FreqTrSwp2 ;
        Summary(index).AmpAccom2           = AmpAccom2 ; 
        Summary(index).HWAccom2            = HWAccom2 ; 
        Summary(index).ISIRatio1toAll2     = ISIRatio1toAll2 ;
        Summary(index).AdaptIndexTS2       = AdaptIdx2 ;
        Summary(index).OnsetTSAP2        = OnsetTSFAP2 ; 
        Summary(index).onsetrapidityTSAP2 = onsetrapidity2 ; 
        Summary(index).ThreshTSAP2       = ThreshTSAP2  ; 
        Summary(index).TSAP2basetothresh    = TSAP2basetothresh ; 
        Summary(index).AmpTSAP2thresh       = AmpTSAP2thresh ;
        Summary(index).TSAP2peaktoahp       = TSAP2peaktoahp ;
        Summary(index).HalfWTSAP2        = HalfWTSAP2 ;
        Summary(index).AHPTSAP2          = AHPTSAP2 ;
        Summary(index).AHPslowTSAP2      = AHPslowTSAP2 ; 
        Summary(index).UpStrokeTSAP2       = UpStrokeTSAP2 ; 
        Summary(index).DwnStrokeTSAP2      = DwnStrokeTSAP2 ; 
        Summary(index).UpDwnStrkRatioTSAP2     = UpDwnStrkRatioTSAP2 ;
        Summary(index).MaxUpStrkTSAP2       = MaxUpStrkTSAP2 ;
        Summary(index).MaxDwnStrkTSAP2      = MaxDwnStrkTSAP2 ; 
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
        Summary(index).threshFB0     = nanmean(aps2.thresh(aps2.freqbin==0)) ;
        Summary(index).threshFB1     = nanmean(aps2.thresh(aps2.freqbin==1)) ;
        Summary(index).threshFB2     = nanmean(aps2.thresh(aps2.freqbin==2)) ;
        Summary(index).threshFB3     = nanmean(aps2.thresh(aps2.freqbin==3)) ;
        Summary(index).threshFB4     = nanmean(aps2.thresh(aps2.freqbin==4)) ;
        Summary(index).threshFB5     = nanmean(aps2.thresh(aps2.freqbin==5)) ;  
        Summary(index).threshFB6     = nanmean(aps2.thresh(aps2.freqbin==6)) ;
        Summary(index).HalfWFB0     = nanmean(aps2.halfwidth(aps2.freqbin==0)) ;
        Summary(index).HalfWFB1     = nanmean(aps2.halfwidth(aps2.freqbin==1)) ;
        Summary(index).HalfWFB2     = nanmean(aps2.halfwidth(aps2.freqbin==2)) ;
        Summary(index).HalfWFB3     = nanmean(aps2.halfwidth(aps2.freqbin==3)) ;
        Summary(index).HalfWFB4     = nanmean(aps2.halfwidth(aps2.freqbin==4)) ;
        Summary(index).HalfWFB5     = nanmean(aps2.halfwidth(aps2.freqbin==5)) ;
        Summary(index).HalfWFB6     = nanmean(aps2.halfwidth(aps2.freqbin==6)) ;
        Summary(index).upstrokeFB0     = nanmean(aps2.upstroke(aps2.freqbin==0)) ;
        Summary(index).upstrokeFB1     = nanmean(aps2.upstroke(aps2.freqbin==1)) ;
        Summary(index).upstrokeFB2     = nanmean(aps2.upstroke(aps2.freqbin==2)) ;
        Summary(index).upstrokeFB3     = nanmean(aps2.upstroke(aps2.freqbin==3)) ;
        Summary(index).upstrokeFB4     = nanmean(aps2.upstroke(aps2.freqbin==4)) ;
        Summary(index).upstrokeFB5     = nanmean(aps2.upstroke(aps2.freqbin==5)) ; 
        Summary(index).upstrokeFB6     = nanmean(aps2.upstroke(aps2.freqbin==6)) ; 
        Summary(index).downstrokeFB0     = nanmean(aps2.downstroke(aps2.freqbin==0)) ;
        Summary(index).downstrokeFB1     = nanmean(aps2.downstroke(aps2.freqbin==1)) ;
        Summary(index).downstrokeFB2     = nanmean(aps2.downstroke(aps2.freqbin==2)) ;
        Summary(index).downstrokeFB3     = nanmean(aps2.downstroke(aps2.freqbin==3)) ;
        Summary(index).downstrokeFB4     = nanmean(aps2.downstroke(aps2.freqbin==4)) ;
        Summary(index).downstrokeFB5     = nanmean(aps2.downstroke(aps2.freqbin==5)) ; 
        Summary(index).downstrokeFB6     = nanmean(aps2.downstroke(aps2.freqbin==6)) ; 
        Summary(index).onsetrapFB0     = nanmean(aps2.onsetrapidity(aps2.freqbin==0)) ;
        Summary(index).onsetrapFB1     = nanmean(aps2.onsetrapidity(aps2.freqbin==1)) ;
        Summary(index).onsetrapFB2     = nanmean(aps2.onsetrapidity(aps2.freqbin==2)) ;
        Summary(index).onsetrapFB3     = nanmean(aps2.onsetrapidity(aps2.freqbin==3)) ;
        Summary(index).onsetrapFB4     = nanmean(aps2.onsetrapidity(aps2.freqbin==4)) ;
        Summary(index).onsetrapFB5     = nanmean(aps2.onsetrapidity(aps2.freqbin==5)) ;
        Summary(index).onsetrapFB5     = nanmean(aps2.onsetrapidity(aps2.freqbin==6)) ;
        Summary(index).ampFB0     = nanmean(aps2.amp(aps2.freqbin==0)) ;
        Summary(index).ampFB1     = nanmean(aps2.amp(aps2.freqbin==1)) ;
        Summary(index).ampFB2     = nanmean(aps2.amp(aps2.freqbin==2)) ;
        Summary(index).ampFB3     = nanmean(aps2.amp(aps2.freqbin==3)) ;
        Summary(index).ampFB4     = nanmean(aps2.amp(aps2.freqbin==4)) ;
        Summary(index).ampFB5     = nanmean(aps2.amp(aps2.freqbin==5)) ; 
        Summary(index).ampFB6     = nanmean(aps2.amp(aps2.freqbin==6)) ; 
        Summary(index).ahpFB0     = nanmean(aps2.relahp(aps2.freqbin==0)) ;
        Summary(index).ahpFB1     = nanmean(aps2.relahp(aps2.freqbin==1)) ;
        Summary(index).ahpFB2     = nanmean(aps2.relahp(aps2.freqbin==2)) ;
        Summary(index).ahpFB3     = nanmean(aps2.relahp(aps2.freqbin==3)) ;
        Summary(index).ahpFB4     = nanmean(aps2.relahp(aps2.freqbin==4)) ;
        Summary(index).ahpFB5     = nanmean(aps2.relahp(aps2.freqbin==5)) ; 
        Summary(index).ahpFB6     = nanmean(aps2.relahp(aps2.freqbin>=6)) ; 
        Summary(index).amp1to40  = nanmean(aps2.amp(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).upstroke1to40  = nanmean(aps2.maxdvdt(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).downstroke1to40  = nanmean(aps2.mindvdt(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).HalfW1to40  = nanmean(aps2.halfwidth(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).thresh1to40  = nanmean(aps2.thresh(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).onsetrap1to40  = nanmean(aps2.onsetrapidity(aps2.freqbin>=1 & aps2.freqbin<=2)) ;
        %Summary(index).updwnratio1to40  = nanmean(aps2.updownratio(aps2.freqbin>=1 & aps2.freqbin<=2)) ;   
        Summary(index).amp41to80  = nanmean(aps2.amp(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).upstroke41to80  = nanmean(aps2.maxdvdt(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).downstroke41to80  = nanmean(aps2.mindvdt(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).HalfW41to80  = nanmean(aps2.halfwidth(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).thresh41to80  = nanmean(aps2.thresh(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).onsetrap41to80  = nanmean(aps2.onsetrapidity(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
        %Summary(index).updwnratio41to60  = nanmean(aps2.updownratio(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
        Summary(index).amp100up  = nanmean(aps2.amp(aps2.freqbin>=6)) ; 
        Summary(index).upstroke100up  = nanmean(aps2.maxdvdt(aps2.freqbin>=6)) ;    
        Summary(index).downstroke100up  = nanmean(aps2.mindvdt(aps2.freqbin>=6)) ;    
        Summary(index).HalfW100up  = nanmean(aps2.halfwidth(aps2.freqbin>=6)) ;    
        Summary(index).thresh100up  = nanmean(aps2.thresh(aps2.freqbin>=6)) ; 
        Summary(index).onsetrap100up  = nanmean(aps2.onsetrapidity(aps2.freqbin>=6)) ;
        %end %if at line 64 
        %clear variables assigned in "sweep" For loop
        %clearvars -except Summary i basedir savedir savename filelist index Summary1
        index = index + 1 ;
         end
%% save
% Define the full file paths for saving
%save('/Users/elinemertens/Data/ephys/Masterstruct/dataStruct5.mat', 'dataStruct');
% Convert Summary struct to a table
Summary_T5 = struct2table(Summary);
% Construct the full path to the Excel file
full_path = fullfile(savedir, [savename '.xlsx']);
% Write the table to the Excel file
writetable(Summary_T5, full_path);
%% Clear all variables except Summary_T
clearvars -except Summary_T1 Summary_T2;
%%
    ApsTimeArray = ApsTimeArray' ; 
    ApsDataArray = ApsDataArray' ;   
    DerivativesArray = DerivativesArray' ; 
 % Calculate the average of the first 5 columns for each row
% Calculate the average of the first 5 columns for each row
ApsData_avg = mean(ApsDataArray(:, 1:end), 2);
Derivatives_avg = mean(DerivativesArray(:, 1:end),2);
% Add the average as a new column to ApsDataArray
% ApsDataArray_with_avg = [ApsData_avg, ApsDataArray];
% Derivatives_with_avg = [Derivatives_avg,DerivativesArray];

% Create a figure
 % Apply a moving average filter to smooth the data
    smoothed_ApsDataArray = smoothdata(ApsDataArray, 'movmean', 5);
    smoothed_ApsDerivativesArray = smoothdata(DerivativesArray, 'movmean', 5);
figure;
subplot(2,2,1)
% Plot individual traces in red
for i = 1:size(ApsDataArray, 2)-1
    plot(ApsDataArray(:, i), 'Color', 'red');
    hold on;
end

% Plot average trace in thick black
plot(ApsData_avg, 'Color', 'black', 'LineWidth', 2);
% Customize the plot as needed
xlabel('Time');
ylabel('ApsData');
legend('Trace 1', 'Trace 2', 'Trace 3', 'Trace 4', 'Trace 5', 'Average', 'Location', 'best');
title('ApsData Traces');

subplot(2,2,2)
for i = 1:size(DerivativesArray, 2)-1
    plot(DerivativesArray(:, i), 'Color', 'red');
    hold on;
end
plot(Derivatives_avg, 'Color', 'black', 'LineWidth', 2);
% Add any additional customization as needed
hold off;


    subplot(2,2,3);
    plot(ApsDataArray, DerivativesArray, 'b', 'LineWidth', 2);
    title('Original Phase Plane');
    xlabel('Voltage');
    ylabel('dV/dt');

    subplot(2,2,4);
    plot(smoothed_ApsDataArray, smoothed_ApsDerivativesArray, 'r', 'LineWidth', 2);
    title('Smoothed Phase Plane');
    xlabel('Voltage');
    ylabel('dV/dt');
%end
    %%
save(fullfile(savedir, savename), 'Summary') ;
clearvars -except Summary i


%%
% Convert Summary struct to a table
Summary_T = struct2table(Summ);

% Define the full path to the Excel file
full_path = '/Volumes/Overig/summary/summary_h245_test.xlsx';

% Write the table to the Excel file
writetable(Summ, full_path);


% Write the table to the Excel file
writetable(Summary_T, full_path);








