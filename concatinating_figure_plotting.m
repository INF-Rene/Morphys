data1.marker_sum_norm_label=nan(size(data1,1),1);

for i=1:size(data2,1)
    idx=[];
    idx= strcmp(data2.cell_name(i),data1.Cell_id);
    data1.marker_sum_norm_label(idx)=data2.marker_sum_norm_label(i);
    data1=data1;
end 


data1.Norm_Marker_Sum04_label=nan(size(data1,1),1);

for i=1:size(data2,1)
    idx=[];
    idx= strcmp(data2.cell_name(i),data1.Cell_id);
    data1.Norm_Marker_Sum04_label(idx)=data2.Norm_Marker_Sum_04_label(i);
    data1=data1;
end 

%%
figure
for i=1:size(t_types,2)
    idx2=[];
    idx2= strcmp(t_types{i},data.CellType);
    %M_s20to30Hz
    hold on
    histogram(data.CPM_HCN1(idx2),'FaceColor',rgb_colors_t_type(i,:));
end


%%
figure
for i=1:size(cluster_name,2)
    idx2=[];
    idx2= strcmp(cluster_name{i},data.Cell_cluster);
    %M_s20to30Hz
    hold on
    histogram(data.CPM_HCN1(idx2),'FaceColor',rgb_colors_cluster(i,:));
end
