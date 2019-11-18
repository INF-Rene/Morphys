fn='D:\Allen Data\Ephys\Human\H19.29.144\H19.29.144.11.41\H19.29.144.11.41.01\H19_29_144_11_41_01.nwb';
info=h5info(fn);

swps={info.Groups(1).Groups(2).Groups.Name};
stimswps={info.Groups(6).Groups(1).Groups.Name};
swpstimnames={};
for i=1:numel(swps)
    swpstimnames(i)=h5read(fn, [swps{i} '/stimulus_description']);
end

LBN=h5read(fn, '/general/labnotebook/Dev1/numericalValues');
LBN=squeeze(LBN(1,:,:))';
LBN_keys=h5read(fn, '/general/labnotebook/Dev1/numericalKeys');
LBN=cell2table(num2cell(LBN));
LBN.Properties.VariableNames=genvarname(LBN_keys(:,1));
LBN.Properties.VariableUnits=LBN_keys(:,2);
LBT=h5read(fn, '/general/labnotebook/Dev1/textualValues');
LBT=squeeze(LBT(1,:,:))';
LBT_keys=h5read(fn, '/general/labnotebook/Dev1/textualKeys');
LBT=cell2table(LBT);
LBT.Properties.VariableNames=genvarname(LBT_keys(:,1));

QCfails=any(LBN{:,contains(LBN.Properties.VariableNames, 'QC')}==0,2);
QCfails=unique(LBN.SweepNum(QCfails))+1;

plotswps=find(contains(swpstimnames,'X4PS_SupraThresh'));
%plotswps(ismember(plotswps,QCfails))=[];
data=[];
stimwave=[];
for i=1:numel(plotswps)
    n=h5read(fn, [swps{plotswps(i)} '/num_samples']);
    data(1:n,i)=h5read(fn, [swps{plotswps(i)} '/data']);
    stimwave(1:n,i)=h5read(fn, [stimswps{plotswps(i)} '/data']);
end
plot(data)
max(stimwave)
