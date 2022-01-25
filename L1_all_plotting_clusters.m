fulldata=readtable('/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_PatchSeq/All_L1_matlab_human_nmsscore.csv');
data=fulldata(fulldata.Norm_Marker_Sum04_label==1,:);
%%
%scatter plot sag_ratio VS Ap upstroke
cluster_name={'LAMP5','PAX','Rxrg', 'Vip/Sncg', 'other human'};%, 'NA'};

%%assigning colors to cluster
% for i=1:length(cluster_name)
%     rgb_colors_cluster(i,1:3)=uisetcolor;
% end
%     
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
for i=1:size(cluster_name,2)
    idx=[];
    idx= strcmp(cluster_name{i},data.Cell_cluster);
    %M_s20to30Hz
    %hold on
    scatter(data.(property1)(idx),data.(property2)(idx),100,rgb_colors_cluster(i,:),'filled')
    hold on
    
    [r, p]=corrcoef(data.(property1)(idx),data.(property2)(idx),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', cluster_name{i}, string(r(1,2)),string(p(1,2)))
end
legend(legends)
ylabel(property2)
xlabel(property1)



%%
%scatter plot AP vs HCN expression
figure
legends = {};
gene='CPM_HCN4';
property='Sag_ratio_100mV';
%data.sumHCN=sum(data{:,6:9},2,'omitnan');
for i=1:size(cluster_name,2)
    idx2=[];
    idx2= strcmp(cluster_name{i},data.Cell_cluster);
    %M_s20to30Hz
    %hold on
    scatter(data.(gene)(idx2), data.(property)(idx2),100,rgb_colors_cluster(i,:),'filled')

    hold on
    [r, p]=corrcoef(data.(gene)(idx2),data.(property)(idx2),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', cluster_name{i}, string(r(1,2)),string(p(1,2)));
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
for i=1:size(cluster_name,2)
    idx2=[];
    idx2= strcmp(cluster_name{i},data.Cell_cluster);
 
    scatter(data.(gene)(idx2), data.(property)(idx2),100,rgb_colors_cluster(i,:),'filled')

    hold on
    [r, p]=corrcoef(data.(gene)(idx2),data.(property)(idx2),'Rows','complete');
    r(1,2);
    p(1,2);
    legends{end+1}=sprintf('%s (r=%0.2f, p=%0.2f)', cluster_name{i}, string(r(1,2)),string(p(1,2)))
end
legend(legends)
ylabel(property)
labelxaxis=sprintf('%s expression',gene);
xlabel(labelxaxis);

%%
%plotting HCN1 expression per cluster in boxplot
figure;    
boxplot(data.CPM_HCN4, data.Cell_cluster,'Symbol', 'o', 'Widths',0.5, 'ColorGroup', data.Cell_cluster);
    %'GroupOrder',data.Cell_cluster);
    %'Colors', rgb_colors_cluster(i,:), 'Symbol', 'o', 'Widths',0.5)%,
title('HCN4 expression per clusters')

figure;
for i=1:size(cluster_name,2)
    idx2=[];
    idx2= strcmp(cluster_name{i},data.Cell_cluster);
    boxplot(data.CPM_HCN4, data.Cell_cluster,'Symbol', 'o', 'Widths',0.5);
    hold on
end

