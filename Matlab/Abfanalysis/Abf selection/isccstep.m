function [ccabf, chs] = isccstep(abf)
%ISCCSTEP determine if abf file was created with a current step protocol
%   Returns logical value ccabf
%   Return chs which gives the channelnumbers for which it is true that a
%   ccstep protocol was applied

chs=[abf.channel(:).number];
ccabf=logical(false);

%% Abf level
% Numberofsweeps higher than 5
if abf.nrofsweeps<5, chs=[]; return; end
% Only abfs with samplefreq of at least 10000
if abf.samplefreq<10000, chs=[]; return; end

%% Analogout level
% might be inserted later if we would want to exclude E-Codes

%% Channel&Sweep level

% delete primary channels that are not in current clamp
% delete channels that have sweeps that:
% consist of less then 4 epochs
% do not last between 0.5-6 seconds
% contain no APs in any of the sweeps
for j=chs
        if ~strcmp(abf.channel([abf.channel.number]==j).in(strcmp({abf.channel([abf.channel.number]==j).in.signal},'primary')).units, 'mV')
            chs(chs==j)=[];
        elseif any([abf.channel([abf.channel.number]==j).in(strcmp({abf.channel([abf.channel.number]==j).in.signal},'primary')).sweep.nrofepochs]<4)
            chs(chs==j)=[];
        elseif ~all(0.5<cellfun(@(x) seconds(duration(str2double(strsplit(x, ':')))),{abf.channel([abf.channel.number]==j).in(strcmp({abf.channel([abf.channel.number]==j).in.signal},'primary')).sweep.timespan})<6)
            chs(chs==j)=[];
        elseif all([abf.channel([abf.channel.number]==j).in(strcmp({abf.channel([abf.channel.number]==j).in.signal},'primary')).sweep.nrofaps]==0)
            chs(chs==j)=[];
        end
end
if isempty(chs), chs=[]; return; end

%% Epochs

% all epoch A, B and C need to be steps
% epoch B needs to be 500-1500 ms
% epoch A has a stepdiff of 0
% the stepdiff of C is always equal to -1*(stepdiff of b)
% epoch B needs to be at least 4 different stepdiffs
for j=chs
    epochABC=[abf.channel([abf.channel.number]==j).in(strcmp({abf.channel([abf.channel.number]==j).in.signal},'primary')).sweep(:).epoch];
    epochABC=epochABC(ismember({epochABC.Name},{'Epoch A', 'Epoch B', 'Epoch C'} ));
    if ~all(strcmp({epochABC.typestr}, 'step'))
        chs(chs==j)=[];
    elseif ~all(0.5<cellfun(@(x) seconds(duration(str2double(strsplit(x,':')))),{epochABC([epochABC.number]==1).timespan})<1.5)
        chs(chs==j)=[];
    elseif ~all([epochABC([epochABC.number]==0).stepdiff]==0)
        chs(chs==j)=[];
    elseif ~all([epochABC([epochABC.number]==1).stepdiff]== -[epochABC([epochABC.number]==2).stepdiff])
        chs(chs==j)=[];
    elseif numel(unique([epochABC([epochABC.number]==1).stepdiff]))<4
        chs(chs==j)=[];
    elseif ~(any([epochABC([epochABC.number]==1).stepdiff]<0) && any([epochABC([epochABC.number]==1).stepdiff]>0))
        chs(chs==j)=[];
    end
end
if isempty(chs), chs=[]; return; end
%All conditions for CCstep are met:
ccabf=logical(true);

end

