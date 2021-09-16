%select folder to analyze 
folder = uigetdir 
cd (folder);
%make a list with all nwbs
list = dir();
list = struct2table(list);
list = list(list.bytes>1000,:); %only files with actual data

%% USE THIS nwb2: 
for i = 1:numel(list.name)
    fn =cell2mat(list.name(i));
 nwb = NWBfile(fn,[{'LP'} {'hresh'} {'CCSteps_DA_0'} {'LSFINEST'} {'LSCOARSE'}]);
 obj =nwb.analyseNWB ;
 obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
 saveme(obj,'/Users/elinemertens/Data/ephys/nwb2 analyzed/198_2', obj.savename) 
end

    
    %%
fn= '/Users/elinemertens/Data/ephys/Hippocampus/165/H19.29.165.21.41_ephys_converted/H19.29.165.21.41.02-nwb2.nwb'


%%
nwb = NWBfile(fn,[{'LP'} {'hresh'} {'CCSteps_DA_0'} {'LSFINEST'} {'LSCOARSE'}])
%%
obj=nwb.analyseNWB
%%
obj.saveme('/Users/elinemertens/Data/ephys/Hippocampus/165/165_nwb2_analysis', '165_rec2');

%%
%%obj.saveme('/Users/elinemertens/Data/ephys/Human/H20.29.185.21.01/nwb analyzed','185_cell1.mat');

%%
%%save('myfile.mat', 'nwb')


%% 01092021 
fn= '/Users/elinemertens/Data/ephys/Human_nwb2/H21.29.198.21/H21.29.198.21.11.02.nwb'
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