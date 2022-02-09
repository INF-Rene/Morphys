%select folder to analyze 
folder = uigetdir 
cd (folder);
%make a list with all nwbs
list = dir();
list = struct2table(list);
list = list(list.bytes>10000,:); %only files with actual data

%% USE THIS nwb2: 
for i = 1:numel(list.name) 
    fn =cell2mat(list.name(i));
 nwb = NWBfile(fn,[ {'hresh'} {'Steps'} {'CC'}]);
 obj =nwb.analyseNWB ;
 obj.savename = sprintf('NWB_%s.mat',obj.filename(1:end-4));
 saveme(obj,'/Users/elinemertens/Data/ephys/Hippocampus/nwb2_analyzed/2022_NEW/198', obj.savename) 
end
