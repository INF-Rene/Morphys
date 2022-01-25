
data=readtable('/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/All_L1_matlab_human_nmsscore.csv');
%%
%scatter plot sag_ratio VS Ap upstroke

t_types={,'Inh L1 SST NMBR (ADARB2+)',...
    'Inh L1-2 GAD1 MC4R (ADARB2+)', 'Inh L1-2 LAMP5 DBP', ...
    'Inh L1-2 PAX6 CDH12', ... }
    'Inh L1-2 VIP PCDH20','Inh L1-2 VIP TSPAN12', ...
    'Inh L1-3 PAX6 SYT6 (Sncg)',...
    'Inh L1-4 LAMP5 LCP2 (rosehip)'};

% missing t type = 'Inh L1-3 VIP ADAMTSL1', 'Inh L1-3 SST CALB1', 'Inh L1-2 VIP LBH'
% t types with NaN p-value (small N?) 'Inh L1-2 PAX6 TNFAIP8L3', 'Inh L1-2
% SST BAGE2 (ADARB2+)', 'Inh L1 SST CHRNA4 (ADARB2+)', n = 3 , 7, 4

%assigning colors to t-type

%for i=1:length(t_types)
 % rgb_colors_t_type(i,1:3)=uisetcolor;
%end
    
%%
%scatter plot for AP vs sag

figure('Position', [100 100 1000 700]);
%k=1;
%for i=1:size(data,1)
%    sag(i)=cell2mat(data.Sag_ratio_100mV(i));
%   upstroke(i,1)=cell2mat(data.M_s20to30Hz(i));
%end

legends = {};
property1='Sag_ratio_100mV';
property2='M_s20to30Hz';
for i=1:size(t_types,2)
    idx=[];
    idx= strcmp(t_types{i},data.CellType);
    %M_s20to30Hz
    %hold on
    scatter(data.(property1)(idx),data.(property2)(idx),100,rgb_colors_t_type(i,:),'filled')
    hold on
    
    [r, p]=corrcoef(data.(property1)(idx),data.(property2)(idx),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', t_types{i}, string(r(1,2)),string(p(1,2)))
end
legend(legends);
ylabel(property2)
xlabel(property1)

%%
%scatter plot AP vs HCN expression
figure
legends = {};
gene='CPM_HCN3';
%property='Sag_ratio_100mV';
property='Sag_ratio_100mV';
for i=1:size(t_types,2)
    idx2=[];
    idx2= strcmp(t_types{i},data.CellType);
    %M_s20to30Hz
    %hold on
    scatter(data.(gene)(idx2), data.(property)(idx2),100,rgb_colors_t_type(i,:),'filled')

    hold on
    [r, p]=corrcoef(data.(gene)(idx2),data.(property)(idx2),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', t_types{i}, string(r(1,2)),string(p(1,2)));
end
legend(legends)
ylabel(property)
labelxaxis=sprintf('%s expression',gene);
xlabel(labelxaxis);
%%
%scatter plot AP vs SCN expression
figure
legends = {};
gene='CPM_SCN8A';
property='M_s20to30Hz';
for i=1:size(t_types,2)
    idx2=[];
    idx2= strcmp(t_types{i},data.CellType);
    %M_s20to30Hz
    %hold on
    scatter(data.(gene)(idx2), data.(property)(idx2),100,rgb_colors_t_type(i,:),'filled')

    hold on
    [r, p]=corrcoef(data.(gene)(idx2),data.(property)(idx2),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', t_types{i}, string(r(1,2)),string(p(1,2)))
end
legend(legends)
ylabel(property)
labelxaxis=sprintf('%s expression',gene);
xlabel(labelxaxis);