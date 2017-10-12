function [ccabf] = getccstep(abf)
%UNTITLED4 Select CC steps protocols from Abf struct created by
%ImportCSVstoStruct function
%   Detailed explanation goes here



%% Abf level
% Only abfs with numberofsweeps higher than 5
ccabf=abf([abf.nrofsweeps]>5);
% Only abfs with samplefreq of at least 10000
ccabf=ccabf([ccabf.samplefreq]>=10000);

%% Analogout level
% might be inserted later if we would want to exclude E-Codes

%% Channel&Sweep level

% delete primary channels that are not in current clamp
% delete channels that have sweeps that:
% consist of less then 4 epochs
% do not last between 0.5-6 seconds
% contain no APs in any of the sweeps

for i=1:size(ccabf,1)
    for j=size(ccabf(i).channel,1):-1:1
        if ~strcmp(ccabf(i).channel(j).in(strcmp({ccabf(i).channel(j).in.signal},'primary')).units, 'mV')
            ccabf(i).channel(j)=[];
        elseif any([ccabf(i).channel(j).in(strcmp({ccabf(i).channel(j).in.signal},'primary')).sweep.nrofepochs]<4)
            ccabf(i).channel(j)=[];
        elseif ~all(0.5<cellfun(@(x) seconds(duration(str2double(strsplit(x, ':')))),{ccabf(i).channel(j).in(strcmp({ccabf(i).channel(j).in.signal},'primary')).sweep.timespan})<6)
            ccabf(i).channel(j)=[];
        elseif all([ccabf(i).channel(j).in(strcmp({ccabf(i).channel(j).in.signal},'primary')).sweep.nrofaps]==0)
            ccabf(i).channel(j)=[];
        end
    end
end
deleteemptyabfs()

%% Epochs

% all epoch A, B and C need to be steps
% epoch B needs to be 500-1500 ms
% epoch A has a stepdiff of 0
% the stepdiff of C is always equal to -1*(stepdiff of b)
% epoch B needs to be at least 4 different stepdiffs
for i=1:size(ccabf,1)
    for j=size(ccabf(i).channel,1):-1:1
        epochABC=[ccabf(i).channel(j).in(strcmp({ccabf(i).channel(j).in.signal},'primary')).sweep(:).epoch];
        epochABC=epochABC(ismember({epochABC.Name},{'Epoch A', 'Epoch B', 'Epoch C'} ));
        if ~all(strcmp({epochABC.typestr}, 'step'))
            ccabf(i).channel(j)=[];
        elseif ~all(0.5<cellfun(@(x) seconds(duration(str2double(strsplit(x,':')))),{epochABC([epochABC.number]==1).timespan})<1.5)
            ccabf(i).channel(j)=[];
        elseif ~all([epochABC([epochABC.number]==0).stepdiff]==0)
            ccabf(i).channel(j)=[];
        elseif ~all([epochABC([epochABC.number]==1).stepdiff]== -[epochABC([epochABC.number]==2).stepdiff])
            ccabf(i).channel(j)=[];
        elseif numel(unique([epochABC([epochABC.number]==1).stepdiff]))<4
            ccabf(i).channel(j)=[];
        elseif ~(any([epochABC([epochABC.number]==1).stepdiff]<0) && any([epochABC([epochABC.number]==1).stepdiff]>0))
            ccabf(i).channel(j)=[];
        end
    end
end
deleteemptyabfs()



%%
    function deleteemptyabfs()
       ccabf=ccabf(~cellfun(@isempty,{ccabf.channel}));
    end


end

