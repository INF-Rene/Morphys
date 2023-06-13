%% Make an Abffile object.
% To make an Abffile object, we need a file location of a ".abf" file, and a Setupsettings object describing the setup
% settings that apply to the file. 

%We dont use any BB files here!
% if the batch isn't needed, we dont need the bb. We won't add this to any
% of the upcoming codes either. We will remove bb.getabf(1).plot and replace it
% with abf.plot 

p1   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn1  = '2020_03_04_0002.abf';
fp1  = fullfile(p1,fn1);
ss  = load('/Users/elinemertens/Downloads/Morphys-master/Morphys/Data/Electrophysiology/SetupSettings/Setupsettings_INF.mat');
ss  = ss.obj;
abf1 = Abffile(fp1,ss);

%%
    
p2   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn2  = '2020_03_04_0021.abf';
fp2  = fullfile(p2,fn2);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf2 = Abffile(fp2,ss);

p3   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn3  = '2020_03_04_0045.abf';
fp3  = fullfile(p3,fn3);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf3 = Abffile(fp3,ss);

p4   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn4  = '2020_03_04_0064.abf';
fp4  = fullfile(p4,fn4);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf4 = Abffile(fp4,ss);

p5   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn5  = '2020_03_05_0001.abf';
fp5  = fullfile(p5,fn5);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf5 = Abffile(fp5,ss);

p6   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn6  = '2020_03_05_0021.abf';
fp6  = fullfile(p6,fn6);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf6 = Abffile(fp6,ss);

p7   = ('/Users/elinemertens/Matlab Data/H20.29.174');
fn7  = '2020_03_05_0041.abf';
fp7  = fullfile(p7,fn7);
ss  = load('/Users/elinemertens/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
abf7 = Abffile(fp7,ss);

%%

%abf.channels(1).analogins(1).units = 'mV';

  abf1.channels(1).analogins(1).units = 'mV';
    swps1 = abf1.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps1)
        swps1(j).Data = swps1(j).Data./20;
        swps1(j).DataInfo.Units = 'mV';
        epochs = swps1(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps1(j).epochs = epochs;
    end
    abf1.channels.analogins(1).sweeps=swps1;
    
    abf2.channels(1).analogins(1).units = 'mV';
    swps2 = abf2.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps2)
        swps2(j).Data = swps2(j).Data./20;
        swps2(j).DataInfo.Units = 'mV';
        epochs = swps2(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps2(j).epochs = epochs;
    end
    abf2.channels.analogins(1).sweeps=swps2;
   

    abf3.channels(1).analogins(1).units = 'mV';
    swps3 = abf3.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps3)
        swps3(j).Data = swps3(j).Data./20;
        swps3(j).DataInfo.Units = 'mV';
        epochs = swps3(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps3(j).epochs = epochs;
    end
    abf3.channels.analogins(1).sweeps=swps3;
   
    
    abf4.channels(1).analogins(1).units = 'mV';
    swps4 = abf4.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps4)
        swps4(j).Data = swps4(j).Data./20;
        swps4(j).DataInfo.Units = 'mV';
        epochs = swps4(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps4(j).epochs = epochs;
    end
    abf4.channels.analogins(1).sweeps=swps4;
   
      abf5.channels(1).analogins(1).units = 'mV';
    swps5 = abf5.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps5)
        swps5(j).Data = swps5(j).Data./20;
        swps5(j).DataInfo.Units = 'mV';
        epochs = swps5(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps5(j).epochs = epochs;
    end
    abf5.channels.analogins(1).sweeps=swps5;
 
  
      abf6.channels(1).analogins(1).units = 'mV';
    swps6 = abf6.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps6)
        swps6(j).Data = swps6(j).Data./20;
        swps6(j).DataInfo.Units = 'mV';
        epochs = swps6(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps6(j).epochs = epochs;
    end
    abf6.channels.analogins(1).sweeps=swps6;
 
      abf7.channels(1).analogins(1).units = 'mV';
    swps7 = abf7.getchannel.getin('signal', 'primary').getsweep;
    for j=1:numel(swps7)
        swps7(j).Data = swps7(j).Data./20;
        swps7(j).DataInfo.Units = 'mV';
        epochs = swps7(j).getepoch;
        for k=1:numel(epochs)
            epochs(k).Data = epochs(k).Data./20;
            epochs(k).DataInfo.Units = 'mV';
           end
        swps7(j).epochs = epochs;
    end
    abf7.channels.analogins(1).sweeps=swps7;
    
  

% hierboven zet je normaal het script neer voor de omrekening pA mV


%% 5a. Plotting objects
% plot the Abffile object
close all
abf.plot


%% 5b. Plot a Channel
close all
%bb.getabf(1).getchannel(1).plot
abf.getchannel(1).plot

%% 5c. Plot a selection of Sweeps
close all
% bb.getabf(end).getchannel(1).selectsweep(1:3).plot
abf(end).getchannel(1).selectsweep(12).plot
abf2(end).getchannel(1).selectsweep(13).plot

%% 5d. Use plot just as you would use the native Matlab plot function
s = abf(end).getchannel(1).getin('signal','primary').getsweep(11);
%upper is just in this case, i want to see sweep 10, EM edit
close all
s.plot('linewidth',2,'linestyle','-.','color','g')

%% 8a. Analysis
% EM! super handige functie voor analyze!
% A set of basic analysis functions are available. These are applied only to primary Analogin signals, with mV units. 
% In short, the analysis procedures will analyse any Actionpotential in Trace, fit any voltage deflection that is the result
% of a step current injection, and calculate sag currents etc when necessary. 
s1 = abf1(end).getchannel(1).getin('signal','primary').getsweep(12);
s1 = s1.analysesweep;
s1.plotanalysis

s2 = abf2(end).getchannel(1).getin('signal','primary').getsweep(13); 
s2 = s2.analysesweep;
s2.plotanalysis

s3 = abf3(end).getchannel(1).getin('signal','primary').getsweep(12); 
s3 = s3.analysesweep;
s3.plotanalysis

s4 = abf4(end).getchannel(1).getin('signal','primary').getsweep(13); 
s4 = s4.analysesweep;
s4.plotanalysis

s5 = abf5(end).getchannel(1).getin('signal','primary').getsweep(12); 
s5 = s5.analysesweep;
s5.plotanalysis

s6 = abf6(end).getchannel(1).getin('signal','primary').getsweep(12); 
s6 = s6.analysesweep;
s6.plotanalysis

s7= abf7(end).getchannel(1).getin('signal','primary').getsweep(13); 
s7 = s7.analysesweep;
s7.plotanalysis



%% 8b. The events in the plot are timeseries 'events'. 
help tsdata.event

%% 8c. Inspect the Actionpotentials details
close all
% all of the AP = aps = s.getap;
% only analyze the first AP:
aps1 = s1.getap(1); 
aps2 = s2.getap(1);
aps3 = s3.getap(1);
aps4 = s4.getap(1);
aps5 = s5.getap(1);
aps6 = s6.getap(1);
aps7 = s7.getap(1);

%%
amplitudes = aps2.get('amp');

%% threshold = aps.get('thresh');

aps1.plot('superimpose','peak');
hold on
aps2.plot('superimpose','peak');
aps3.plot('superimpose','peak');
aps4.plot('superimpose','peak');
aps5.plot('superimpose','peak');
aps6.plot('superimpose','peak');
aps7.plot('superimpose','peak');
hold off

legend show
figure(1)
title('AP traces Temp')
xlabel('time (ms)')
ylabel('mV')
xlim([-2 6])
%ylim([-0.4 0.8])
%line([obj.peak_time obj.peak_time],[obj.thresh obj.peak],'linewidth',2,'lineStyle',':','color','k')

           
% NOTE: although an Actionpotential really is just a short section of Trace, Actionpotentials do not inherit from Trace. This
% would probably have been a lot nicer and could be done in the future, but for now left as is. The Trace is now stored as a
% property called "ts", which is a timeseries object.

%%
epoch1 = s1.getepoch.getsteadystate
epoch2 = s2.getepoch.getsteadystate
epoch3 = s3.getepoch.getsteadystate
epoch4 = s4.getepoch.getsteadystate
epoch5 = s5.getepoch.getsteadystate
epoch6 = s6.getepoch.getsteadystate
epoch7 = s7.getepoch.getsteadystate

%%

epoch1b = s1.getepoch(1).getsteadystate
epoch2b = s2.getepoch(1).getsteadystate
epoch3b = s3.getepoch(1).getsteadystate
epoch4b = s4.getepoch(1).getsteadystate
epoch5b = s5.getepoch(1).getsteadystate
epoch6b = s6.getepoch(1).getsteadystate
epoch7b = s7.getepoch(1).getsteadystate

%%
% this way it can directly be set into a table
T = table(epoch1, epoch2, epoch3, epoch4, epoch5, epoch6, epoch7)

%%
hw1 = aps1.get('halfwidth');
hw2 = aps2.get('halfwidth');
hw3 = aps3.get('halfwidth');
hw4 = aps4.get('halfwidth');
hw5 = aps5.get('halfwidth');
hw6 = aps6.get('halfwidth');
hw7 = aps7.get('halfwidth');

Thw = table(hw1, hw2, hw3, hw4, hw5, hw6, hw7)

%% if you want to get the max / min dvdt 
%(also known as upstroke & downstroke) you can get these measures from the first AP as well 

us1 = aps1.get('maxdvdt');
us2 = aps2.get('maxdvdt');
us3 = aps3.get('maxdvdt');
us4 = aps4.get('maxdvdt');
us5 = aps5.get('maxdvdt');
us6 = aps6.get('maxdvdt');
us7 = aps7.get('maxdvdt');


Tus = table(us1, us2, us3, us4, us5, us6, us7)

%% if you want to know the risetime

rt1 = aps1.get('maxdvdt_time');
rt2 = aps2.get('maxdvdt_time');
rt3 = aps3.get('maxdvdt_time');
rt4 = aps4.get('maxdvdt_time');
rt5 = aps5.get('maxdvdt_time');
rt6 = aps6.get('maxdvdt_time');
rt7 = aps7.get('maxdvdt_time');

Trt = table(rt1, rt2, rt3, rt4, rt5, rt6, rt7)


%% downstroke
ds1 = aps1.get('mindvdt');
ds2 = aps2.get('mindvdt');
ds3 = aps3.get('mindvdt');
ds4 = aps4.get('mindvdt');
ds5 = aps5.get('mindvdt');
ds6 = aps6.get('mindvdt');
ds7 = aps7.get('mindvdt');

Tds = table(ds1, ds2, ds3, ds4, ds5, ds6, ds7)

%%
dt1 = aps1.get('mindvdt_time');
dt2 = aps2.get('mindvdt_time');
dt3 = aps3.get('mindvdt_time');
dt4 = aps4.get('mindvdt_time');
dt5 = aps5.get('mindvdt_time');
dt6 = aps6.get('mindvdt_time');
dt7 = aps7.get('mindvdt_time');

Tdt = table(dt1, dt2, dt3, dt4, dt5, dt6, dt7)


%% 8d. Analyse whole batch
bb = bb.analysebatch;

%% 8e. output tables with (meta) data
bb.abftable
bb.channeltable
bb.analogintable
% this one is very relaxed 
bb.sweeptable
bb.epochtable
bb.aptable


%%
%Take average
% AVG = ((aps1 + aps2 + aps3 + aps4 + aps5 + aps6)/6)
%apss = [aps1 aps2 aps3 aps4 aps5 aps6];
apsmean = mean([aps1 aps2 aps3 aps4 aps5 aps6]);


figure(1)
plot(apsmean)
title('average APs')
xlabel('tijd')
ylabel('amplitude')

%% 9. Save abfs as mat objects
path2save = fullfile('/Users/elinemertens/Morphys/Data/Electrophysiology/Abffiles/MBV/mat_tests');
if ~isdir(path2save), mkdir(path2save); end
for i=1:bb.nrofabfs
    bb.getabf(i).saveme(path2save,bb.getabf(i).filename(1:end-4))
    %bb.getabf(i).saveme(path2save,bb.getabf(i).filename('fn'))
end
 
%%
figure(1)
abf.getchannel(1).selectsweep([1 3 5 7 9 11 13 15]).plot  ;

