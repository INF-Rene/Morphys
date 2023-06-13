%% Analysis script
% Written by D.B. Heyer !
close all, clear all

%% Set path to load and save data; mat data load 
basedir = '/Users/elinemertens/Data/ephys/Analyzed/242 new' ;
savedir = '/Users/elinemertens/Data/ephys/Hippocampus/2022_Summary';
savename = 'Summary_198' ;  

%% load file list
fileinfo  = dir(fullfile(basedir,'*.mat'));
filelist  = {fileinfo.name};

    %% Loop through abfs
index = 1 ; 
for i = 1:length(filelist)
    %% Make subset of data per abf file
    fprintf('Looking for CC-step protocols: file nr %1.0f \n', i);
    load(fullfile(basedir,filelist{i})) ;
    
   % either run through all protocols, if it can't resolve this, remove %
   % on line 22-24 and put it at 24. It will now loop through stimsets to
   % find ccsteps 
%    
%      stimnms={obj.getstimsets.name};
%      CCsteploc=cellfun(@(x) contains(x, 'teps'), stimnms);
%    stimsets = struct2table(obj.getstimsets(CCsteploc).metadata, 'AsArray', true) ;
    stimsets = struct2table(obj.getstimsets.metadata, 'AsArray', true) ;
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

    %% If abf is a stepprotocol: continue with analysis

        fprintf('Retrieving analysis parameters from CC-step file %1.0f \n', index);
        %% Analyze
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
        % calculate variables
        vmbase = [sweep.vmbase] ;
        Freqs= NaN(NrofSweeps,1) ;
        StimInts= NaN(NrofSweeps,1) ;
        currInjections_R= NaN(NrofSweeps,1);
        voltageResponses= NaN(NrofSweeps,1);
        taus= NaN(NrofSweeps,1);
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
       

        
%         % Get the total number of stimsets
% numStimsets = numel(obj.stimsets());
% 
% % Iterate over each stimset
% for v = 1:numStimsets
%     % Get the total number of sweeps for the current stimset
%     numSweeps = numel(obj.stimsets(v).nwbchannels.sweeps());
%     
%     % Iterate over each sweep
%     for w = 1:numSweeps
%         % Check if downstroke is NaN for the current sweep
%         if obj.stimsets(v).nwbchannels.sweeps(w).aps(1, 1).downstroke == NaN
%             % Remove the entire row of sweeps for the current sweep
%             obj.stimsets(v).nwbchannels.sweeps(w) = [];
%             
%             % Decrement the sweep counter to adjust for the removed sweep
%             j = j - 1;
%             
%             % Update the total number of sweeps
%             numSweeps = numSweeps - 1;
%         end
%     end
% end
%         
%         % Create a logical index for rows containing NaN in the specified variable
% dsproblem = ismissing(aps.downstroke);
% dsproblem = ismissing(aps2.downstroke);
% 
% % Remove rows where the specified variable is NaN
% aps(dsproblem, :) = [];
% aps2(dsproblem, :) = [];
% 

% if isnan(obj.stimsets().nwbchannels.sweeps().aps(1, 1).downstroke)
%     % Remove the entire row of sweeps
%     obj.stimsets().nwbchannels.sweeps = [];
% end
% 
% for t = 1:NrofSweeps 
% if isnan(sweep(t).ap(1).downstroke)
%     sweep(t).ap(1) = [] ;
% end
% end



        % find trainsweep 
        % Traincurr=rheobase+50 :
        step = find(strcmp({sweep(frstspikeswp).epoch.idxstr}, 'B'));
       % changed this to 50 pA
        TrainCurr = sweep(frstspikeswp).epoch(step).amplitude +50 ;
        for j = 1:NrofSweeps
            step = find(strcmp({sweep(j).epoch.idxstr}, 'B'));
            tmp(j) = abs(sweep(j).epoch(step).amplitude - TrainCurr) ;
%          if isnan(sweep(j).ap(1).downstroke)
%                   sweep(j).ap(1) = [] ;
%          end
        end
        
        [TrSwp TrSwp] = min(tmp) ;
        CurrAbvRheo=NaN;
        for TrainSweep = TrSwp:NrofSweeps  
            step = find(strcmp({sweep(TrainSweep).epoch.idxstr}, 'B'));
            if length(sweep(TrainSweep).ap) > 1
                 TrSweepCrit=1;
                CurrAbvRheo = sweep(TrainSweep).epoch(step).amplitude - (TrainCurr-50) ;     
%                 if isnan(sweep(TrainSweep).ap(1).downstroke)
%                     TrainSweep = TrainSweep + 1 ;  
%                     if isnan(aps.downstroke(TrainSweep))
%                     TrainSweep = TrainSweep + 2 ; 
                break
               end
            if TrainSweep==NrofSweeps
                TrSweepCrit=0;
            end  
        end
        %    end
%         end
%         
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
%         
        % third try
         TrainCurr2 = sweep(frstspikeswp).epoch(step).amplitude +150 ;
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
                CurrAbvRheo2 = sweep(TrainSweep2).epoch(step).amplitude - (TrainCurr2-150) ;               
                break
               end
            if TrainSweep2==NrofSweeps
                TrSweepCrit2=0;
            end  
        end
       
        
        if ~isempty(sweep(TrainSweep2).ap)
            for l = 1:length(sweep(TrainSweep2).ap)           
                if ~isempty(sweep(TrainSweep2).ap(l,1).halfwidth)
                    HWsTS2(l,1) = [sweep(TrainSweep2).ap(l,1).halfwidth] ;  
                end       
            end
        else
            HWsTS2=NaN;
        end
        
        if TrSweepCrit2==1
            TSbasetothresh2 = ([sweep(TrainSweep2).ap.thresh]-sweep(TrainSweep2).vmbase) ;
            TSpeaktoahp2 = ([sweep(TrainSweep2).ap.ahp_time]-[sweep(TrainSweep2).ap.peak_time]); 
            AmpsTSthresh2 = [sweep(TrainSweep2).ap.amp] ;
            AHPsTS2 = [sweep(TrainSweep2).ap.relahp] ;
            AHPslowTS2 = [sweep(TrainSweep2).ap.relahp_slow] ;
            ISIsTS2 = [sweep(TrainSweep2).ap(2:end).isi] ;
            FreqTrSwp2 = mean([sweep(TrainSweep2).ap(4:end).freq]) ;
            NrOfAPsTrSwp2 = length(sweep(TrainSweep2).ap) ; 
            OnsetTSFAP2 = sweep(TrainSweep2).ap(1).thresh_time - (seconds(sum([sweep(TrainSweep2).epoch(1:find(strcmp({sweep(TrainSweep2,1).epoch.idxstr}, 'A'))).timespan]))*1000) ;
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
        %commented on 04-05-2023 by deduction from abf of lower part
%         tmp = abs(MinVmResponse+100) ;
%         tmp = tmp(~isnan(tmp));
%         [sagswp sagswp] = min(tmp)  ;
%  sagsweepname = sweep(sagswp).Name ;
 
  tmp = abs(MinVmResponse+100) ;
        [sagswp sagswp] = min(tmp) ;
        tmp2 = abs(MinVmResponse+90) ; 
        [sagswp2 sagswp2] = min(tmp2) ;
        tmp3 = abs(MinVmResponse+80) ;
        [sagswp3 sagswp3] = min(tmp3) ; 
         
        Sag                = sweep(sagswp,1).epoch(step).sag / PkDeflect(sagswp,1) ;
        Sag2               = sweep(sagswp2,1).epoch(step).sag / PkDeflect(sagswp2,1) ;
        Sag3               = sweep(sagswp3,1).epoch(step).sag / PkDeflect(sagswp3,1) ;
        SagMedian             = nanmedian([Sag Sag2 Sag3]) ; 
        
 
%  if you want to plot all sag sweeps, get the percentage away form
%  obj.getstimset and run this part via index 
% figure(3)
%  for j = 1:length(obj.stimsets)
%             if any(ismember({obj.getstimset(j).getnwbchannel.getsweep.Name},sagsweepname))
%                 %%figure(); obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).plot
%                 obj.getstimset(j).getnwbchannel.getsweep('Name',sagsweepname).plot;
%                 legend(filelist)
%                 xlim([0 1800])
%                 title('Sag')
%                 grid off
%                  set(gca, 'TickDir', 'out')
%                  ylabel('mV')
%               xlabel('ms')
%             end
%         end

        % calculate input frequency curve
       % this caused problems after fixing the input resistance 
%         Freqs = Freqs(Freqs~=0) ;
%         StimInts = StimInts(StimInts~=0) ;
%         if length(Freqs) > 1
%             [fitFi]=fit(StimInts,Freqs,f_R, 'StartPoint', [0 0]); 
%             FrqChngStimInt = fitFi.R ;
%         else  
%             FrqChngStimInt = NaN ;   
%         end




        % bursting & adaptation index
        if TrSweepCrit==1
            isis=ISIsTS; %exclude stuttering ISIs since that can mess with adaptation analysis
            amps=AmpsTSthresh;
            hws=HWsTS;
        else
            isis=NaN;
            amps=NaN;
            hws=NaN;
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
        
        
         % Firstsweep and lastsweep ISIS EM edits
          % nr of APs first sweep (Eline Edit)
   
        if length(sweep(frstspikeswp).ap) > 1
            isis_FS = [sweep(frstspikeswp).ap(2:end).isi];
            isis_FS1 = [sweep(frstspikeswp).ap(2).isi];
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
        
      
         if TrSweepCrit2 == 0
           isis_TS2 = NaN ;
           isis_TS2_1 = NaN ; 
         elseif length(sweep(TrainSweep2).ap) > 1
            isis_TS2 = [sweep(TrainSweep2).ap(2:end).isi];
            isis_TS2_1 = [sweep(TrainSweep2).ap(2).isi];          
         end

        if length(sweep(NrofSweeps).ap) > 1
            isis_LS = [sweep(NrofSweeps).ap(2:end).isi];
            isis_LS1 = [sweep(NrofSweeps).ap(2).isi];
        else 
            isis_LS = NaN ;
            isis_LS1 = NaN  ;
        end
     
        %get the holding current 
        lb = obj.getstimset(1).getnwbchannel(1).labbooknum ; 
        vals=lb.I0x2DClampHoldingLevel(~isnan(lb.I0x2DClampHoldingLevel)) ; 
        currentinj_avg = nanmean(vals) ; 
        
% %end if you want to make a large summary, remove the end here 
% TrainSweep1 = sweep(TrainSweep).Name ; 
if isempty(TrainSweep2)
    TrainSweep2 = NaN;
end

% figure(5)
%         for j = 1:length(obj.stimsets)
%             if any(ismember({obj.getstimset(j).getnwbchannel.getsweep.Name},TrainSweep1)) 
%                 %figure(); obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).plot
%                  obj.getstimset(j).getnwbchannel.getsweep('Name',TrainSweep1).plot ;
%                 legend(filelist)
%                 xlim([-5 2000])
%                 title('trainsweep')
%             %   ylabel('mV')
%             %    xlabel('ms')
%             end
%         end


        
%          figure(5)
%         for j = 1:length(obj.stimsets)
%             if any(ismember({obj.getstimset(j).getnwbchannel.getsweep.Name},TrainSweep)) 
%                 %figure(); obj.getstimset(j).getnwbchannel.getsweep('Name',firstsweepname).plot
%                  obj.getstimset(j).getnwbchannel.getsweep('Name',TrainSweep).getepoch(step).aps(1).plot('superimpose','peak');
%                 legend(filelist)
%                 xlim([-5 10])
%               %  title('First AP')
%             %   ylabel('mV')
%             %    xlabel('ms')
%             end
%         end
      

        % Create summary  
        Summary(index).File               = stimset(1).filename ;
        Summary(index).Date               = obj.filetimestart ;
        Summary(index).UserID             = obj.userid ;
        Summary(index).guid               = obj.guid ;
        Summary(index).Channel            = NaN ;
        Summary(index).scalefactor        = NaN ;
        Summary(index).holdingcurrent     = NaN ;
       % Summary(index).currentinj         = currentinj_avg ;
        Summary(index).holdingvoltage     = NaN ;
       % Summary(index).NrofSweeps         = NrofSweeps ;
       % Summary(index).PDur               = seconds(sweep(1).epoch(strcmp({sweep(1).epoch.idxstr}, 'B')).timespan)*1000 ;
        Summary(index).FrstP              = sweep(1).currinj ;
        Summary(index).DeltaP             = sweep(2).currinj - sweep(1).currinj ;
        Summary(index).Rheobase           = sweep(frstspikeswp).currinj ;
        Summary(index).FrstSpikeSwp       = frstspikeswp ; 
        Summary(index).TrainSwp           = TrainSweep ; 
        Summary(index).TrainSwp2          = TrainSweep2 ; 
        Summary(index).CurrAbvRheo        = CurrAbvRheo ;
        Summary(index).CurrAbvRheo2        = CurrAbvRheo2 ;
        Summary(index).NrOfAPfrstSwp = length(sweep(frstspikeswp).ap) ;
        Summary(index).NrofAPtrainSwp = length(sweep(TrainSweep).ap);
        Summary(index).NrofAPtrainSwp2 = length(sweep(TrainSweep2).ap);
        Summary(index).NrofAPlastSwp = length(sweep(NrofSweeps).ap);
        Summary(index).isis_FS            = isis_FS ;
        Summary(index).isis_FS1           = isis_FS1 ;
        Summary(index).isis_TS            = isis_TS ;
        Summary(index).isis_TS1           = isis_TS1 ;
        Summary(index).isis_TS2           = isis_TS2 ;
        Summary(index).isis_TS2_1         = isis_TS2_1 ;
        Summary(index).isis_LS            = isis_LS ;
        Summary(index).isis_LS1           = isis_LS1 ;
        Summary(index).vmbaseM            = nanmean(vmbase) ;
        Summary(index).vmbaseSD           = nanstd(vmbase) ;
        Summary(index).Jitter             = nanmean([sweep.jitter]) ;
        Summary(index).InputR             = Rin ;% in MOhm...
        Summary(index).FreqMax            = freqmax ;
        Summary(index).NrOfAPs            = sum(NrofAPs) ;
        Summary(index).NrOfAPsMax         = max(NrofAPs) ;
        Summary(index).FreqTrSwp          = FreqTrSwp ;
        %Summary(index).NrOfAPsTrSwp       = NrOfAPsTrSwp ; 
      %  Summary(index).FrqChngStimInt     = FrqChngStimInt ;
        Summary(index).FrqChngStimInt     = NaN ;
        Summary(index).NrofRBAPs          = sum(NrofRBAPs) ;
        Summary(index).NrofRBAPsM         = NrofRBAPsM ;
        Summary(index).Sag                = Sag ;
        Summary(index).Sag2               = Sag2 ; 
        Summary(index).Sag3               = Sag3 ; 
        Summary(index).SagMedian          = SagMedian ; 
        Summary(index).VmatSag            = MinVmResponse(sagswp,1) ;
        Summary(index).VmatSag2           = MinVmResponse(sagswp2,1) ;
        Summary(index).VmatSag3           = MinVmResponse(sagswp3,1) ; 
        Summary(index).TauM               = nanmean(taus(taus~=0)) ;
        Summary(index).TauSD              = nanstd(taus(taus~=0)) ;
       % Summary(index).curr_FrstAP        = sweep(frstspikeswp).ap(1).currinj ;
        Summary(index).OnsetFAP        = sweep(frstspikeswp).ap(1).thresh_time - (seconds(sum([sweep(frstspikeswp).epoch(1:find(strcmp({sweep(frstspikeswp,1).epoch.idxstr}, 'A'))).timespan]))*1000) ; 
        Summary(index).ThreshFAP       = sweep(frstspikeswp).ap(1).thresh ; 
        Summary(index).FAPbasetothresh    = sweep(frstspikeswp).ap(1).thresh-sweep(frstspikeswp).vmbase ; 
        Summary(index).AmpFAPthresh       = sweep(frstspikeswp).ap(1).amp ;
        Summary(index).PeakFrstAP         = sweep(frstspikeswp).ap(1).thresh + sweep(frstspikeswp).ap(1).amp ; 
        Summary(index).FAPpeaktoahp       = sweep(frstspikeswp).ap(1).ahp_time - sweep(frstspikeswp).ap(1).peak_time ;
        Summary(index).HalfWFAP        = sweep(frstspikeswp).ap(1).halfwidth ; 
        Summary(index).AHPFAP          = sweep(frstspikeswp).ap(1).relahp ;
        Summary(index).AHPslowFAP      = sweep(frstspikeswp).ap(1).relahp_slow ;
        Summary(index).UpStrokeFAP       = sweep(frstspikeswp).ap(1).upstroke ;
        Summary(index).DwnStrokeFAP      = sweep(frstspikeswp).ap(1).downstroke ;
        Summary(index).UpDwnStrkRatioFAP     = abs(sweep(frstspikeswp).ap(1).upstroke) / abs(sweep(frstspikeswp).ap(1).downstroke) ;
        Summary(index).MaxUpStrkFAP       = sweep(frstspikeswp).ap(1).maxdvdt ;
        Summary(index).MaxDwnStrkFAP      = sweep(frstspikeswp).ap(1).mindvdt ;
        Summary(index).OnsetTSAP        = sweep(TrainSweep).ap(1).thresh_time - (seconds(sum([sweep(TrainSweep).epoch(1:find(strcmp({sweep(TrainSweep,1).epoch.idxstr}, 'A'))).timespan]))*1000) ; 
        Summary(index).ThreshTSAP       = sweep(TrainSweep).ap(1).thresh ; 
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
        Summary(index).OnsetTSAP2        = sweep(TrainSweep2).ap(1).thresh_time - (seconds(sum([sweep(TrainSweep2).epoch(1:find(strcmp({sweep(TrainSweep2,1).epoch.idxstr}, 'A'))).timespan]))*1000) ; 
        Summary(index).ThreshTSAP2       = sweep(TrainSweep2).ap(1).thresh ; 
        Summary(index).TSAP2basetothresh    = sweep(TrainSweep2).ap(1).thresh-sweep(TrainSweep2).vmbase ; 
        Summary(index).AmpTSAP2thresh       = sweep(TrainSweep2).ap(1).amp ;
        Summary(index).TSAP2peaktoahp       = sweep(TrainSweep2).ap(1).ahp_time - sweep(TrainSweep2).ap(1).peak_time ;
        Summary(index).HalfWTSAP2        = sweep(TrainSweep2).ap(1).halfwidth ; 
        Summary(index).AHPTSAP2          = sweep(TrainSweep2).ap(1).relahp ;
        Summary(index).AHPslowTSAP2      = sweep(TrainSweep2).ap(1).relahp_slow ;
        Summary(index).UpStrokeTSAP2       = sweep(TrainSweep2).ap(1).upstroke ;
        Summary(index).DwnStrokeTSAP2      = sweep(TrainSweep2).ap(1).downstroke ;
        Summary(index).UpDwnStrkRatioTSAP2     = abs(sweep(TrainSweep2).ap(1).upstroke) / abs(sweep(TrainSweep2).ap(1).downstroke) ;
        Summary(index).MaxUpStrkTSAP2       = sweep(TrainSweep2).ap(1).maxdvdt ;
        Summary(index).MaxDwnStrkTSAP2      = sweep(TrainSweep2).ap(1).mindvdt ;
        Summary(index).AdaptIndexTS       = AdaptIdx ;
        Summary(index).AmpAccomTS         = AmpAccom ;
        Summary(index).HWAccomTS          = HWAccom ;
        Summary(index).ISIRatio1toAll     = ISIRatio1toAll ;
        Summary(index).amp1to20  = nanmean(aps2.amp(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).upstroke1to20  = nanmean(aps2.upstroke(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).downstroke1to20  = nanmean(aps2.downstroke(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).HalfW1to20  = nanmean(aps2.halfwidth(aps2.freqbin>=1 & aps2.freqbin<=2)) ;    
        Summary(index).thresh1to20  = nanmean(aps2.thresh(aps2.freqbin>=1 & aps2.freqbin<=2)) ; 
        Summary(index).onsetrap1to20  = nanmean(aps2.onsetrapidity(aps2.freqbin>=1 & aps2.freqbin<=2)) ;
        Summary(index).updwnratio1to20  = nanmean(aps2.updownratio(aps2.freqbin>=1 & aps2.freqbin<=2)) ;   
        Summary(index).amp21to40  = nanmean(aps2.amp(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).upstroke21to40  = nanmean(aps2.upstroke(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).downstroke21to40  = nanmean(aps2.downstroke(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).HalfW21to40  = nanmean(aps2.halfwidth(aps2.freqbin>=3 & aps2.freqbin<=4)) ;    
        Summary(index).thresh21to40  = nanmean(aps2.thresh(aps2.freqbin>=3 & aps2.freqbin<=4)) ; 
        Summary(index).onsetrap21to40  = nanmean(aps2.onsetrapidity(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
        Summary(index).updwnratio21to40  = nanmean(aps2.updownratio(aps2.freqbin>=3 & aps2.freqbin<=4)) ;
       
        
        
        
        %end %if at line 64 
        %clear variables assigned in "sweep" For loop
      %  clearvars -except Summary i basedir savedir savename filelist index
        index = index + 1 ;
    end
%% save
save(fullfile(savedir, savename), 'Summary') ;
clearvars -except Summary i


%%

Summary_T = struct2table(Summary) ; 
writetable(Summary_T, 'summary_h235_test.xlsx');






