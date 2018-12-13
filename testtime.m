% load all analyzed ccstep files and create analysis overviews for
% inspection
load('D:\Morphys\Data\Labbooks\NAG\MetadataABF\DataSummaryNAG2.mat')

filelist=squeeze(struct2cell(Summary));
filelist=filelist(1,:)';
filelist=cellfun(@(x) strsplit(x, '.'), filelist, 'UniformOutput', false);
filelist=cellfun(@(x) x{1}, filelist, 'UniformOutput', false);
filelist=strcat(filelist, '.mat');
filelist=strcat('D:\Morphys\Data\Electrophysiology\Abffiles\NAG\abf\Converted\',filelist);

 load(filelist{81});
 
 aa=obj.analyseabf;