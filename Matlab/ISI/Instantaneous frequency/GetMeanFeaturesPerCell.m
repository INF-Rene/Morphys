% calculates mean feature for all spikes per cell
function GetMeanFeaturesPerCell(cellnames, path2saveData, path2saveSpikes)
binnames={'firstAP','s0to10Hz','s10to20Hz','s20to30Hz','s30to40Hz','s40to50Hz',...
     's50to60Hz','s60to70Hz','s70to80Hz', 's80to90Hz','s90to100Hz'}
upstroke=table;
downstroke=table;
peak_v=table;
width=table;
threshold=table;
upstroke_downstroke_ratio=table;
upstroke1=table;
downstroke1=table;
peak_v1=table;
width1=table;
threshold1=table;
upstroke_downstroke_ratio1=table;


for m=1:size(cellnames,1)
    name=sprintf('%s',cellnames{m},'_spikes.mat');
    file2load=fullfile(path2saveSpikes,name);
    load(file2load,'data_all');
    if isempty(data_all)
        continue
    end
    x=data_all.inst_freq; 
    edges=[0:10:100]; 
    vars={'upstroke','downstroke','peak_v','width', 'threshold','upstroke_downstroke_ratio'};
    for k=1:size(vars,2)   
        y=data_all.(vars{k}); 
        binNr =data_all.freq_bin ;
       %get the data for each bin separately
        tempData=[];
        for i=1:length(edges)
        tempData=y(binNr==i);  
        mean.(vars{k})(1)=y(1);
        stdev.(vars{k})(1)=0;
     
            if size(tempData,1)>1
            mean.(vars{k})(i+1)=nanmean(tempData);
            stdev.(vars{k})(i+1)=std(tempData, 'omitnan');
            else
                mean.(vars{k})(i+1)=NaN;
                stdev.(vars{k})(i+1)=NaN;   
            end
         end   
    end
% display the results
%     figure
%         xas=[0:10:110];
%         errorbar(xas, mean.(vars{1}),stdev.(vars{1}),'o','MarkerFaceColor','b');
%         xlabel('instantaneous frequency, Hz');
%         ylabel(vars{1});
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize= 14;
        
  % put all the data in a table (only mean) with rows for each cell  

     upstroke1.cell_id=cellnames(m);
     downstroke1.cell_id=cellnames(m);
     peak_v1.cell_id=cellnames(m);
     width1.cell_id=cellnames(m);
     threshold1.cell_id=cellnames(m);
     upstroke_downstroke_ratio1.cell_id=cellnames(m);
     for n=1:length(binnames) 
         upstroke1.(binnames{n})=mean.upstroke(n);
         downstroke1.(binnames{n})=mean.downstroke(n);
         peak_v1.(binnames{n})=mean.peak_v(n);
         width1.(binnames{n})=mean.width(n);
         threshold1.(binnames{n})=mean.threshold(n);
         upstroke_downstroke_ratio1.(binnames{n})=mean.upstroke_downstroke_ratio(n);
     end
    % add new data for each cell to the big table
    upstroke=[upstroke;upstroke1];
    downstroke=[downstroke;downstroke1];
    peak_v=[peak_v;peak_v1];
    width=[width;width1];
    threshold=[threshold;threshold1];
    upstroke_downstroke_ratio=[upstroke_downstroke_ratio;upstroke_downstroke_ratio1];

    
     % save the table
    file2save=fullfile(path2saveData,'data_per_cell_binned.mat');
    save(file2save,'upstroke','downstroke', 'peak_v','width','threshold',...
    'upstroke_downstroke_ratio')


end


end