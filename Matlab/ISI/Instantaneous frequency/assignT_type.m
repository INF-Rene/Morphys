% gets t-types from 'human_shiny_273.csv' file
file_ttypes='/Users/natalia/Documents/Research/Projects_active/Transcriptomics/Shiny_output_data/human_shiny_273.csv';
A=readtable(file_ttypes);
ar_id=num2str(A.ar_id);
for i=1:size(A.ar_id,1)
    cell_id_number(i,1)=A.cell_id(i);
end
for i=1:size(upstroke,1)
    idx(i,1)=find(cell_id_number==str2num(cell2mat(upstroke{i,1})));
    t_type(i,1)=A.topLeaf_label(idx(i,1));
end
upstroke.t_type=t_type;
downstroke.t_type=t_type;
peak_v.t_type=t_type;
width.t_type=t_type;
threshold.t_type=t_type;
upstroke_downstroke_ratio.t_type=t_type;
% find how many ttypes there are
t_types=unique(t_type);
% make an empty table
upstroke_per_t_type=table;
downstroke_per_t_type=table;
peak_v_per_t_type=table;
width_per_t_type=table;
threshold_per_t_type=table;
upstroke_downstroke_ratio_per_t_type=table;

% calculates mean feature per t-type
for i=1:size(t_types,1)
    idx2=[];
    idx2=strcmp(upstroke.t_type,t_types(i));
    newupstroke=table;
    newupstroke=upstroke(idx2,:);
    mean_upstroke(i,:)=nanmean(table2array(upstroke(idx2,2:12)),1);
    stdev_upstroke(i,:)=std(table2array(upstroke(idx2,2:12)),1,'omitnan');
    upstroke_per_t_type=[upstroke_per_t_type;newupstroke];
    
    newdownstroke=table;
    newdownstroke=downstroke(idx2,:);
    downstroke_per_t_type=[downstroke_per_t_type;newdownstroke];
    
    newpeak_v=table;
    newpeak_v=peak_v(idx2,:);
    peak_v_per_t_type=[peak_v_per_t_type;newpeak_v];
    
    newwidth=table;
    newwidth=width(idx2,:);
    width_per_t_type=[width_per_t_type;newwidth];
    
    newthreshold=table;
    newthreshold=threshold(idx2,:);
    threshold_per_t_type=[threshold_per_t_type;newthreshold];
    
    newupstroke_downstroke_ratio=table;
    newupstroke_downstroke_ratio=upstroke_downstroke_ratio(idx2,:);
    upstroke_downstroke_ratio_per_t_type=[upstroke_downstroke_ratio_per_t_type;newupstroke_downstroke_ratio];
end