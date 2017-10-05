%% test script
close all, clear all

%% load bb
load('C:\Users\DBHeyer\Documents\PhD\Human Database\bbtest.mat')

%% get tables
abfs     = bb.abftable ;
channels = bb.channeltable ;
ins      = bb.analogintable ;
sweeps   = bb.sweeptable ;
epochs   = bb.epochtable ;
aps      = bb.aptable ;

%% get new batch
cc = Abfbatch('gui');
%% analyze
tic
cc = cc.analysebatch ;
toc
%% test
bb.getabf(1).getchannel(1).getin('signal','primary').getsweep(14).plotanalysis

cc.getabf(1).getchannel(1).getin('signal','primary').getsweep(8).plotanalysis

%%


for i = 1:cc.getabf(1).getchannel(1).getin('signal','primary').getsweep(21).nrofaps 
    if ~isempty(cc.getabf(1).getchannel(1).getin('signal','primary').getsweep(21).getap(i).thresh)
    maxdvdts(i) = cc.getabf(1).getchannel(1).getin('signal','primary').getsweep(21).getap(i).thresh ;
    end
end






