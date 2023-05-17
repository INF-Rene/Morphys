%select folder to analyze 
folder = uigetdir 
cd (folder);
%make a list with all nwbsephys

list = dir();
list = struct2table(list);
list = list(list.bytes>10000,:); %only files with actual data

%% USE THIS nwb2: 
for i = 1:numel(list.name) 
    fn =cell2mat(list.name(i));
 nwb = NWBfile(fn,[{'LP'} {'hresh'} {'Steps'} {'CC'} {'LSFINEST'} {'LSCOARSE'}]);
 obj =nwb.analyseNWB ;
 obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
 saveme(obj,'/Users/elinemertens/Data/ephys/Analyzed/242', obj.savename) 
end

    
    %%
fn= '/Users/elinemertens/Data/Projects/NEW_Hippocampus_2022_07/Data/mouse/test/M23.29.45995.21.01.02.nwb'

%% let op of je wel of geen chirp wilt 
nwb = NWBfile(fn,[{'LP'} {'hresh'} {'CC'} {'teps'} {'LSFINEST'} {'LSCOARSE'}]);

%% ONLY CC
nwb = NWBfile(fn,[{'CC'} {'teps'}]);

%% everything
nwb = NWBfile(fn,[{'TRIPLE'} {'X6S'}]); 
%%
obj=nwb.analyseNWB
%%
obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
saveme(obj,'/Users/elinemertens/Data/Projects/NEW_Hippocampus_2022_07/Data/mouse/mat analysed/new', obj.savename) 

%%
%%obj.saveme('/Users/elinemertens/Data/ephys/Human/H20.29.185.21.01/nwb analyzed','185_cell1.mat');

%%
%%save('myfile.mat', 'nwb')


%% 01092021 
fn= '/Users/elinemertens/Data/Projects/NEW_Hippocampus_2022_07/Data/mouse/resonance/to be analysed/M23.29.45994.21.03.02_M-compressed.nwb'
nwb = NWBfile(fn,[{'LP'} {'CCSteps_DA_0'}]);

%%
%%ccstep = nwb.getstimset('name', 'CCSteps_DA_0');
%%ccstep_channel = ccstep.getnwbchannel(1);

%% nwb 1

%%fn='/Users/elinemertens/Data/Analysis/Hippocampus/Hipp patient 3/Kopie van H19.29.165.21.41.01.nwb';
%%nwb = NWBfile(fn) 

%% make a selection of stimsets and/or sweeps
%%stimsetfilter={'Ramp', 'TRIPLE', 'CC' }
%%sweepselect=[1:58 60:100];
%%nwb = NWBfile(fn, {'CC'}, sweepselect) 