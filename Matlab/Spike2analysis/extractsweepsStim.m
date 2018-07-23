%chdata=load('C:\Users\DBHeyer\Documents\PhD\China\Djai\Recordings\mats\stim\180621_012.mat');
chdata=load('C:\Users\DBHeyer\Documents\PhD\China\Djai\Recordings\mats\stim\180705_002.mat');
%settings
ch='Ch1';
secch='Ch2';
stimch='Ch3';
StringCh='Ch31';
CCStepkey='B';
stopkey='Z';

%sweeplength=10001.3425; %length of each sweep in ms
sweeplength=10001.621;


%find beginning and end of CC step protocol
CCStepkey=unicode2native(CCStepkey);
stopkey=unicode2native(stopkey);
starttime=chdata.(StringCh).times(chdata.(StringCh).codes == CCStepkey);
stoptime=chdata.(StringCh).times(chdata.(StringCh).codes == stopkey);
starttime=starttime(end);
stoptime=stoptime(end);
%stoptime= 184;
secoffset=chdata.(stimch).start-chdata.(ch).start;

% create timeseries
ts=timeseries(chdata.(ch).values,[0:(chdata.(ch).length-1)].*chdata.(ch).interval);
ts=ts.getsampleusingtime(starttime,stoptime);
sects=timeseries(chdata.(stimch).values,[0:(chdata.(stimch).length-1)].*chdata.(stimch).interval+secoffset);
sects=sects.getsampleusingtime(starttime,stoptime);

%find epoch edges and holding values in secondary timeseries
%baseline=ts.getsamples(1:1000).mean;    %warning: assumption that the first 1000 samples are always without current injection
%sects=sects-baseline;
Wn = [(10000/(25000/2))];
[B,A] = butter(2,Wn);
sectscorr=sects.filter(B,A);
[~,edgeidx]=findpeaks(abs(diff(sectscorr.data)), 'MinPeakHeight',3, 'MinPeakDistance', 10);

%add sweep edges as epoch edges (assuming that sweeps always start and end
%with the same holding value!!)
keydelay=1; %approximate delay of actually starting the protocol and registering the keystroke in Ch31. May differ from computer to computer due to speed
numsweeps=floor((stoptime-starttime-keydelay/1000)/((sweeplength-1)/1000));
sweepidx1 = round([1 [1:numsweeps]*(sweeplength-1)/chdata.(ch).interval/1000]+keydelay/1000/chdata.(ch).interval)';
edgeidx1=edgeidx*(chdata.(stimch).interval/chdata.(ch).interval);
edgeidx1(edgeidx1>max(sweepidx1))=[];
edgeidx1=sort([edgeidx1 ; sweepidx1]);

sweepidx = round([1 [1:numsweeps]*(sweeplength-1)/chdata.(stimch).interval/1000]+keydelay/1000/chdata.(stimch).interval)';
sweepedges=zeros(numsweeps, 2);
edgeidx(edgeidx>max(sweepidx))=[];
edgeidx=sort([edgeidx ; sweepidx]);
%make analogin with sweeps and epochs
a=Spike2Analogin;
a.units=chdata.(ch).units;
a.signal='primary';
a2=Spike2Analogin;
a2.units=chdata.(stimch).units;
a2.signal='Secondary';
for i=1:numsweeps
    sweepedges(i,:) = [sweepidx(i) sweepidx(i+1)-1];
    s = Spike2Sweep('number',i,'guid_IN',a.guid);
    s = s.adddata(get(ts.getsamples(sweepidx1(i):sweepidx1(i+1)-1), 'Data'),chdata.(ch).units);
    s = s.addtime(1/chdata.(ch).interval);
    s2 = Spike2Sweep('number',i,'guid_IN',a.guid);
    s2 = s2.adddata(get(sects.getsamples(sweepidx(i):sweepidx(i+1)-1), 'Data'),chdata.(stimch).units);
    s2 = s2.addtime(1/chdata.(stimch).interval);
    edges=edgeidx(edgeidx>=sweepidx(i) & edgeidx<=sweepidx(i+1));
    for j=1:numel(edges)-1
        epochtab=struct;
        epochtab.idx=j;
        epochtab.number=j-1;
        alphabet={'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
        epochtab.idxstr=alphabet{j};
        epochtab.type=1; %Only steps work for now!
        epochtab.typestr='step';
        epochtab.firstlevel=sects.getsamples(edges(j):edges(j+1)).mean;
        epochtab.pulseperiod=0;
        epochtab.pulsewidth=0;
        epochtab.maxfrequency=0;
        epochtab.datetimestart = '';%epochdatetimestart; %no datetme information                          % add start time
        strttime = (sects.getsamples(edges(j)).Time-sects.getsamples(sweepedges(i,1)).Time)*1000;
        endtime  = (sects.getsamples(edges(j+1)).Time-sects.getsamples(sweepedges(i,1)).Time)*1000;
        epochtab.timespan      = endtime-strttime;
        epochtab.datetimeend   = '';%epochtab.datetimestart + epochtab.timespan;                 % add end time
        
        swpsection          = s.getsampleusingtime(strttime,endtime);
        epochtab.data       = swpsection.Data';
        epochtab.units      = s.DataInfo.Units;
        epochtab.samplefreq = s.samplefreq;
        epochtab.strttime   = strttime;
        
        s = s.addepoch(epochtab); % add epoch object to list
        epochtab2=epochtab;
        swpsection2          = s2.getsampleusingtime(strttime,endtime);
        epochtab2.data       = swpsection.Data';
        epochtab2.units      = s2.DataInfo.Units;
        epochtab2.samplefreq = s2.samplefreq;
        epochtab2.strttime   = strttime;
        s2 = s2.addepoch(epochtab2);
    end
    if isempty(a.sweeps), a.sweeps        = s;
    else                  a.sweeps(end+1) = s;
    end
    a = updatesweepstats(a);
    if isempty(a2.sweeps), a2.sweeps        = s2;
    else                  a2.sweeps(end+1) = s2;
    end
    a2 = updatesweepstats(a2);
end
c=Spike2Channel;
c.analogins=a;
c.analogins(2)=a2;
c.updatephase=2;
c=c.analysechannel;

c.getin(1).plotanalysis

%% get average trace

for i = 1:numsweeps
    x(:,i) = c.analogins(1,1).sweeps(1,i).Data(1:8000) ;
end

x = mean(x,2) ;

figure; plot(x)







