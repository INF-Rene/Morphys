%% Analysis script
% Written by D.B. Heyer
close all, clear all

%% Set path to load and save data
basedir = 'C:\Users\DBHeyer\Documents\PhD\RNA analysis\Patchseq\Allen\Data\Analyzed' ; % specify path to analyzed abf objects. (.mat files)
savedir = 'C:\Users\DBHeyer\Documents\PhD\RNA analysis\Patchseq\Allen' ; % specify path for saving the data summary
savename = 'CellSummary'; % specify filename to be saved

%alle protocollen apart analyseren, dit tegelijkertijd analyseren


%% import CSV files %if you're loading a lot of files at the same time. not used in this version of the script
%requires specific folder- and filenames in location basedir!
% abfs = readtable(fullfile(basedir, 'Abffiles','Abffiles.txt'),'Delimiter',',') ;
% channels = readtable(fullfile(basedir, 'Channels','Channels.txt'),'Delimiter',',') ;
% ins = readtable(fullfile(basedir, 'Analogins','Analogins.txt'),'Delimiter',',') ;
% outs = readtable(fullfile(basedir, 'Analogouts','Analogout.txt'),'Delimiter',',') ;
% sweeps = readtable(fullfile(basedir, 'Sweeps','Sweeps.txt'),'Delimiter',',') ;
% epochs = readtable(fullfile(basedir, 'Epochs','Epochs.txt'),'Delimiter',',') ;
% aps = readtable(fullfile(basedir, 'Actionpotentials','Actionpotentials.txt'),'Delimiter',',') ;

%hierboven maak je los die tabellen 

% for data Thijs:
%load(fullfile(basedir,'selectedtables.mat'))
%% Locate files 
% je maakt van de file een .mat file, dit doe je ipv abf batch. Dit houdt
% ze apart om ze even te bekijken. Batch is wat meer bij elkaar. 
%je moet de xAbfBulkConversionAnalysis bulk a = abffile(fullfile) a =analayseabs, a.saveme 
% dit zet je erin om een bulk te maken en dan maakt ie overal mat files van
% file van abf naar mat, in files staat soms nwb, ik gebruik abf dus
% vervang dit 

%dit hieronder doet lotte niet, zij doet enkel load(en dan de juiste file)

%fileinfo  = dir(fullfile(basedir,'*.mat'));
%filelist  = {fileinfo.name};
%% Loop through abfs
index = 1 ;
for i = 1:height(filelist)
    fprintf('Looking for CC-step protocols: file nr %1.0f \n', i);
    load(fullfile(basedir,filelist{i})) ;
%hierboven begint loop wanneer je losse bestanden doet 
    
    %lotte heeft van, obj = bb er tussenstaan 
    % extract data tables from obj
    abf(i) = struct2table(obj.getabf.metadata) ; %werkt alleen met abfbatch is, of getabf eruit halen. eerst paar errors, alleen bij out kwam die error.
    
% als je een batch aanmaakt, dan kan je het laten staan 

    channels = struct2table(obj.getabf.getchannel.metadata) ;
    ins = struct2table(obj.getabf.getchannel.getin.metadata) ;
    outs = struct2table(obj.getabf.getchannel.getout.metadata) ;
    sweeps = struct2table(obj.getabf.getchannel.getin.getsweep.metadata) ;
    epochs = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.metadata) ;
    aps = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.getap.metadata) ;
    
    % get instantaneous freq bins, elke regel nieuwe abf bestand en AP
    % eruit vissen en in freq bins zetten. voor alle AP in 1x. geen
    % verschil tussen abf file. 
    edges = 1:10:201 ;
    aps.freqbin = discretize(aps.freq, edges) ;
    aps.freqbin(isnan(aps.freqbin)) = 0 ;
    aps.currinj = aps.number*0 ;
    for ii = 1:height(aps)
        aps(ii,:).currinj = epochs(ismember(epochs.guid,aps(ii,:).parent_guid),:).amplitude ;  
    end
    aps.updownratio = aps.maxdvdt./abs(aps.mindvdt) ;
    aps.onsetrapidity(aps.onsetrapidity > 100) = NaN ;
    
% onsetrapidity moet eigenlijk sampling rate van minimaal 100 kHz hebben,
% meeste zitten rond de 50 en dan krijg je snel foute schattingen. Alles
% boven 100 is standaard fout, dus maakt daar niks van. Heel onbetrouwbaar.
    
    %% Make subset of data per abf file
%
    abf = SubsetTable2struct(abfs(i),channels,ins,outs,sweeps,epochs,aps) ;
    %% If abf is a stepprotocol: continue with analysis
    %[ccabf, chs] = isccstep(abf) ; %not necessary if you're not bulk-loading your abffiles
    %if ccabf == 1 
    %checkt of protocol een step protocol is en in welk kanaal die zit. 
        fprintf('Retrieving analysis parameters from CC-step file %1.0f \n', index);
        %% Analyze
        %sweep = abf.channel([abf.channel.number] == chs).in(strcmp({abf.channel([abf.channel.number] == chs).in.signal},'primary')).sweep ;
        % alleen nodig wanener je veel verschillende channels heb, nu niet
        % nodig. alleen bij batch voeg je al die i toe. 
        sweep = abf(i).channel(1).in(strcmp({abf.channel(1).in.signal},'primary')).sweep ;
        NrofSweeps = length(sweep) ;  
        % find current injection epoch and assign aps to sweep
% hij kiest juiste epoch (kleinere versie van sweep, je hebt meerdere
% events tijdens een sweep, elk event is epoch). hij kijkt of er verschil
% is geweest in current step. als je de current steps bij elkaar optelt,
% dan moet het o zijn, dan is er een duidelijke injectie is 
% moet langer zijn dan 0.03 seconde, zodat het geen test puls is

        for step = 1:length(sweep(1).epoch)
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
              %je creert een sweep daarboven, hiermee voeg je info toe aan
              %sweeps. voeg info toe aan de 'sweep' variabele toe. 
              %if sweep(j).epoch(step+1).ap(ap).start_time > (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+7 && sweep(j).epoch(step+1).ap(ap).start_time < (sum(second({sweep(1).epoch(1:step).timespan}))*1000)+300
                    sweep(j).rbap(ap) = sweep(j).epoch(step+1).ap(ap) ;
                %end 
%niveau van epochs eruit halen, zodat je per sweep kan kijken welke info je
%wil hebben en niet in de verkeerde epochs zoekt 
            end
        end
  
        % find rheobase sweep. welke sweep is eerste met ap, gaat door
        % lijst van sweeps. wanneer er AP is > 0 en het is depolarizing,
        % dan breakt ie hem. Dit doe je zodat je geen spontane AP mee
        % neemt. 

        apcrit=0;
        for frstspikeswp = 1:NrofSweeps
            if sweep(frstspikeswp).epoch(step).nrofaps > 0 && sweep(frstspikeswp).epoch(step).stepdiff > 0d
                apcrit=1;
                break
            end
        end
        
        %je pakt alle APs in de current step zijn voorgekomen.
% dit is in het begin, alle APs van Abf file. ook die niet tijdens current
% steps zijn, of rebound spikes. die wil je niet meenemen eigenlijk, niet
% bij freq. daarom even subselectie erin gemaakt. Welke voldoen aan de
% steps. 
%aps 
        
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
% maakt eerst lege objecten aan, is beter voor de loop.
        vmbase = [sweep.vmbase] ;
        Freqs=[];
        StimInts=[];
        currInjections_R=[];
        voltageResponses=[];
        taus=[];
        for j = 1:NrofSweeps           
            if sweep(j,1).currinj >= -100 && sweep(j,1).currinj < 0 %hyperpolarizing tussen -100 en 0 current injection for input resistance. 
                %                                                   dit kan je ook berekenen van andere steps
                voltageResponses(j,1) = sweep(j,1).vmresponse ; 
                currInjections_R(j,1) = sweep(j,1).currinj ;
                if sweep(j,1).epoch(step).tau < 100 && sweep(j,1).epoch(step).tau > 0 && sweep(j,1).epoch(step).gof > 0.95 %goodness of fit = bepaalde power, tau is decay time van hyperpolerizing steps
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
            if sweep(j,1).currinj < 0 && isfield(sweep, 'rbap') %ap die hyperpolarizing current step volgen
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
        


        % find trainsweep: waar meeste AP features bekeken worden. Hoe
        % allen dit ook doet. ipv alleen eerste of alle, zeg je: ik kijk
        % naar alle AP van sweep op gestandaardiseerde current step is
        % opgenomen. Hierdoor AP vergelijkbaar met vergelijkbare current
        % step.  je kijkt naar rheobase +50 = 50pA boven rheobase zijn
        % opgenomen. maakt cellen makkelijker te vergelijken.
        
        % Traincurr=rheobase+50 :
        % stuttering cell of niet 
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
        % met traces van -100 en 0 pA current injectie. 
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
        Summary(index).File               = abf(i).filename ;
        Summary(index).Date               = abf(i).filetimestart ;
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
        %nu krijg je de freq bins, eerst keken ze naar first ap, of
        %gemiddelde AP van gestandaardiseerde sweep. Dus dat is voor
        %eerste, dan alle van train. 
        %nu bins, of naar bin 0 (gemiddelde van alle eerste AP van elke
        %sweep). Dan kan je kijken hoe features verandren wanneer frequency
        %toeneemt. 
        
        %hogere frequency bins vaak op hogere current injections. dus hoe
        %hoger de bin, hoe sneller de APs achter elkaar zitten. Bin0 is
        %altijd alleen de eerste AP, dus hier zit 1 AP per cel in. 
        
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
        %nu kan je de samenvatting van de cel in 1 keer opslaan. 
        clearvars -except abf Summary i basedir savename abfs aps channels ins outs epochs sweeps index
        index = index + 1 ;
    %end
end
%% save
save(fullfile(basedir, savename), 'Summary') ;
clearvars -except Summary i










