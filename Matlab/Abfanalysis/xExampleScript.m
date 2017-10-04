% example script
close all, clear all

%% 1a. Set base paths
% specify which folder the Morphys directories are placed
% MAC USERS WILL NEED TO CHANGE SLASHES!!!
dir_base = 'C:\Users\Thijs\Documents';

%% 1b. Send Matlab to correct folder
dir_oop = fullfile(dir_base,'Morphys','Code','Matlab','Abfanalysis');
cd(dir_oop)

%% 2a. Make a setup settings object
% make an empty object
open Setupsettings
ss = Setupsettings;
disp(ss)
help Setupsettings

%% 2b. Inspect setupsettings properties 
properties(ss)

%% 2c. Set a property: give it a name
ss.setupsettingsname = 'Setupsettings_INF';      % directly
ss.set('setupsettingsname','Setupsettings_INF'); % via the 'set' method (safer, more controlled)

% get a property
n = ss.get('setupsettingsname');

%% 2d. Add channels using name/value pairs
ss = ss.addchannel('number',1,'dacnum',0,'primary',0,'secondary',4); % 1st amplifier channel
ss = ss.addchannel('number',2,'dacnum',1,'primary',1,'secondary',5); % 2nd amplifier channel
ss = ss.addchannel('number',3,'dacnum',2,'primary',2,'secondary',6); % 3rd amplifier channel
ss = ss.addchannel('number',4,'dacnum',3,'primary',3,'secondary',7); % 4th amplifier channel
disp(ss)

%% 2e. Remove a channel
ss = ss.removechannel('number',1); % remove amplifier channel with number 1
ss = ss.removechannel(3);          % remove amplifier channel at index 3
disp(ss)

%% 2f. Add complete channels to setup settings object
c1 = Channel('number',1,'dacnum',0,'primary',0,'secondary',4);
c2 = Channel('number',4,'dacnum',3,'primary',3,'secondary',7);
ss = ss.addchannel(c1);
ss = ss.addchannel(c2);
disp(ss)

%% 2g. See how correct inputs are enforced.
% Following lines will give errors:
ss = ss.addchannel('number',4,'dacnum',0,'primary',3,'secondary',7);  % amplifier channel number already specified
ss = ss.addchannel('number',5,'dacnum',4,'primary',8,'secondary',9);  % dacnumber outside expected range (0-3)

%% 2h. Save Setupsettings object.
% The saveme method is inherited from Sharedmethods. All objects inheriting from Sharedmethods will share that function
destination = fullfile(dir_base,'Morphys','Data','Electrophysiology','SetupSettings');
fn = ss.setupsettingsname;
ss.saveme(destination,fn);

%% Make an Abffile object.
% To make an Abffile object, we need a file location of a ".abf" file, and a Setupsettings object describing the setup
% settings that apply to the file. 
p   = fullfile(dir_base,'Morphys','Data','Electrophysiology','Abffiles','MBV','abf');
fn  = '2009_10_01Cel4_0000.abf';
fp  = fullfile(p,fn);
ss  = load('C:\Users\Thijs\Documents\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_MBV.mat');
ss  = ss.obj;
abf = Abffile(fp,ss);

%% 3a. Make a batch with gui
% The Abfbatch object allows you to handle a collection of Abffiles. 
% Use GUI to select abf/mat files to add to batch, and select an setup settings object
bb = Abfbatch('gui');

%% for example, select: '2009_10_01Cel4_0000.abf'    '2009_10_01Cel5_0002.abf'    '2012_09_26Cel03_0013.abf'    '2017_03_30_0169.abf'    '2017_03_30_0170.abf'
% testfiles = {'2009_10_01Cel4_0000.abf'    '2009_10_01Cel5_0002.abf'    '2012_09_26Cel03_0013.abf'    '2017_03_30_0169.abf'    '2017_03_30_0170.abf'};
% testpaths = arrayfun(@(x) fullfile(dir_base,'Morphys','Data','Electrophysiology','Abffiles','MBV',x),testfiles);
% bb = Abfbatch(testfiles,ss);

%% 3b. Inspect the batch, remove nonsense Channels and unwanted Abffiles
bb = bb.inspectbatch;

%% 4a. Navigate through the structures.
% get an Abffile from batch using the "getabf" method
help getabf
a = bb.getabf(1);
disp(a)

%% 4b. Get vs select methods
% to 'get' means obtain an object from list, to 'select' means filter the list for objects selected
bb.getabf(1);                       % returns Abffile object at index 1 in list of abfs, as a 1x1 Abffile object
bb.selectabf(1);                    % returns batch with only the first abf in the list
bb.getabf;                          % returns whole list of abfs
bb.getabf.selectchannel('number',1) % returns 1xN Abffile object only containing data from amplifier Channel number 1.
bb.getabf.getchannel('number',1)    % returns array of all Channels with number equal to 1.
bb.getabf(1).getchannel(2).getin('signal','primary').getsweep(end).getepoch(2) % returns 1 epoch.
bb.getabf.getchannel.getin('signal','primary').getsweep.getepoch % returns every single epoch in batch

%% 5a. Plotting objects
% plot the Abffile object
close all
bb.getabf(1).plot

%% 5b. Plot a Channel
close all
bb.getabf(1).getchannel(2).plot

%% 5c. Plot a selection of Sweeps
close all
bb.getabf(end).getchannel(2).selectsweep(1:3).plot

%% 5d. Use plot just as you would use the native Matlab plot function
s = bb.getabf(end).getchannel(2).getin('signal','primary').getsweep(end-1);
close all
s.plot('linewidth',2,'linestyle','-.','color','g')

%% 6a. Timeseries objects
% Sweeps and epochs are timeseries objects because they inherit from Trace, which in turn inherits from timeseries.
% This means timeseries methods and properties are available for Sweep objects!
t = timeseries; % click "More properties" and "Methods" to see all

%% 6b. View some inherited properties
s.TimeInfo
s.DataInfo
% The 'plot' function called up to now in fact always referred to the timeseries plot function. For this reason, time units
% and data units did not have to be specified when plotting, as the timeseries plot method reads these values from the
% Timeinfo and Datainfo properties.

%% 6c. Use some inherited timeseries methods:
close all; 
s = s.getsampleusingtime(690,720);
s.plot('bo-');

% NOTE 1. no adjustment or updating is done to the Sweep properties when using timeseries methods. For example, when using
% getsampleusingtime, all timeseries inherited properties are updated by Matlab, but the number of Epochs in the
% Sweep epochs list are still the same and have not been adjusted for the shorter sweep length. Keep this in mind. Also
% Sweep start date, end date, and duration remain unchanged.
%
% NOTE 2. tscollection objects cannot be made from Trace/Sweep/Epoch objects. Results in unpredictable behaviour, bad idea.

%% 6d. Some more handy timeseries methods: easily upsample/downsample traces
upsample = 1e-2;
s.resample(s.TimeInfo.Start:upsample:s.TimeInfo.End,'linear').plot('rx-')

%% 7a. Protcol information
% In fixed-wavelength episodic recording mode (so if your protocol has "sweeps", and is not Gap-free), there usually is a
% protocol of steps/ramps/.../sinewaves, with each sweep sectioned into "Epochs". Protocols can be specified in 2 ways:
% 1. Via the pClamp analog waveform tab: a table specifying what to inject (step/ramp/...) where (DAC channel / Analogout)
%    and when (which Epoch) and how long (duration)
% 2. Via the use of a stimulus file, which should be an axon text file (".atf")
% This protocol information is required for any analysis, so that we can section the sweep into its different epochs.
% In case 1, the epoch information is in the Abffile (h.EpochSec, returned by abfload_pro.m), can be entirely obtained
% programmatically and is automatically added when converting an abffile into an Abffile object. The protocol injected into 
% the cell via a Channel/Analogout is stored as an Analogwaveformtab object under the property name "analogwaveformtable" 
% in the Analogout of a channel:
out = bb.getabf(1).getchannel(1).getout;
pro = bb.getabf(1).getchannel(1).getout.analogwaveformtable;
help Analogwaveformtab

%% In case a stimulus file was used: since the epoch information cannot be extracted (yet?) from the atf file alone, it has 
% to be provided. Therefore, a pre-made Analogwaveformtable object is searched in a fixed, predesignated folder that has the 
% exact same name as the stimulus file name listed for the Analogout object.
out = bb.getabf('filename','2017_03_30_0170.abf').getchannel(2).getout;
pro = bb.getabf('filename','2017_03_30_0170.abf').getchannel(2).getout.analogwaveformtable;

% NOTE if the stimulus file is not found, the Analogout will not be analysed with the auto-analysis routine!
out = bb.getabf('filename','2012_09_26Cel03_0013.abf').getchannel(1).getout;
pro = bb.getabf('filename','2012_09_26Cel03_0013.abf').getchannel(1).getout.analogwaveformtable;

%% Making Analogwaveformtab objects
% Now done using an external excell file:
path2ecode = fullfile(dir_base,'Morphys','Data','Electrophysiology','Protocols','AnalogWaveforms','MBV');
fn = 'eCodeEpochList.xlsx';
t  = readtable(fullfile(path2ecode,fn));
open xMakeAnalogwaveformtabs

%% Leading and lagging epochs
% Now we've got our epochs specified, you may think we're done, but pClamp wouldn't be pClamp if that were true. Enter:
% the misterious pClamp 'waiting time'. Specify a step or whatever to occur after 100 ms from sweep onset, and you'll see it
% doesn't start at 100ms. pClamp waits and records (but never mentions this anywhere!), precisely 1/64th of the total number 
% of samples in the sweep before actually starting the protocol you designed in the analog waveform tab. 
% Strangely, pClamp does not do this for protocols where stimfiles are used. So, depending on how the protocol is specified, 
% we may need to add a "leading epoch" (duration = (1/64)*nrofsamples/samplefrequency) in front of the protocol. Also, if the 
% sweep is set to be longer than the total sum of epoch durations, a lagging epoch should also be added to assign all of the 
% Trace to an epoch. Thankfully, it's taken care of now programmatically, but it is good to be aware of this when indexing 
% over epochs. Thank the legendary Tim Kroon for figuring that 1/64 rule out.
e = bb.getabf('filename','2009_10_01Cel4_0000.abf').getchannel(1).getin('signal','primary').getsweep(1);

%% 8a. Analysis
% A set of basic analysis functions are available. These are applied only to primary Analogin signals, with mV units. 
% In short, the analysis procedures will analyse any Actionpotential in Trace, fit any voltage deflection that is the result
% of a step current injection, and calculate sag currents etc when necessary. 
s = bb.getabf(end).getchannel(2).getin('signal','primary').getsweep(1);
s = s.analysesweep;
close all;
s.plotanalysis

%% 8b. The events in the plot are timeseries 'events'. 
help tsdata.event

%% 8c. Inspect the Actionpotentials details
close all
aps = s.getap;
amplitudes = aps.get('amp');
aps.plot('superimpose','peak')

% NOTE: although an Actionpotential really is just a short section of Trace, Actionpotentials do not inherit from Trace. This
% would probably have been a lot nicer and could be done in the future, but for now left as is. The Trace is now stored as a
% property called "ts", which is a timeseries object.

%% 8d. Analyse whole batch
bb = bb.analysebatch;

%% 8e. output tables with (meta) data
bb.abftable
bb.channeltable
bb.analogintable
bb.sweeptable
bb.epochtable
bb.aptable

%% 9. Save abfs as mat objects
path2save = fullfile(dir_base,'Morphys','Data','Electrophysiology','Abffiles','MBV','mat_tests');
if ~isdir(path2save), mkdir(path2save); end
for i=1:bb.nrofabfs
    bb.getabf(i).saveme(path2save,bb.getabf(i).filename(1:end-4))
end

%% 10. Reload example
close all, clear all
% --- load example manually by dragging file into workspace ---
obj.getchannel(2).getin('signal','primary').getsweep(1).plotanalysis



