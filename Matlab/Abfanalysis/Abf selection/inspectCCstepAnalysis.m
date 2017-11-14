
% load all analyzed ccstep files and create analysis overviews for
% inspection
load('D:\Morphys\Data\Labbooks\NAG\MetadataABF2\DataSummaryNAG3.mat')

filelist=squeeze(struct2cell(Summary));
filelist=filelist(1,:)';
filelist=cellfun(@(x) strsplit(x, '.'), filelist, 'UniformOutput', false);
filelist=cellfun(@(x) x{1}, filelist, 'UniformOutput', false);
filelist=strcat(filelist, '.mat');
filelist=strcat('D:\Morphys\Data\Electrophysiology\Abffiles\NAG\abf\Analyzed\',filelist);

for i=247:numel(filelist)
   fprintf('Open file %1.0f out of %1.0f \n', i, numel(filelist))
   load(filelist{i});
   
   close all
   figure ('position', [0 0 1920 1000])
   subplot(2,3,1);
   [~, zeroswp]=min(abs([obj.getchannel.getin('signal', 'primary').getsweep.getepoch('Name', 'Epoch B').amplitude]));
   obj.getchannel.getin('signal', 'primary').getsweep(1:zeroswp-1).plotanalysis
   title('Passive properties')
   
   subplot(2,3,2);
   [~, swp100]=min(abs(-100-[obj.getchannel.getin('signal', 'primary').getsweep.getepoch('Name', 'Epoch B').amplitude]));
   obj.getchannel.getin('signal', 'primary').getsweep(swp100:zeroswp-1).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch('Name', 'Epoch B').Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch('Name', 'Epoch B').Time(1)+100])
   title('Tau')
   
   subplot(2,3,3);
   obj.getchannel.getin('signal', 'primary').getsweep(zeroswp).plotanalysis
   if ~isnan(Summary(i).vmbaseM)
   ylim([Summary(i).vmbaseM-4 Summary(i).vmbaseM+4])
   end
   title('Zero current injection')
   
   subplot(2,3,4);
   obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).getepoch('Name', 'Epoch B').getap(1).start_time-10 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).getepoch('Name', 'Epoch B').getap(1).start_time+20])
   title('First AP')
   
   subplot(2,3,5);
   obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch('Name', 'Epoch B').Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch('Name', 'Epoch B').Time(end)+300])
   title('Trainsweep')
   
   subplot(2,3,6);
   freqs=[];
   for j=1:obj.nrofsweeps
       if numel(obj.getchannel.getin('signal', 'primary').getsweep(j).getepoch('Name', 'Epoch B').getap)>=4
           freqs(j)=mean([obj.getchannel.getin('signal', 'primary').getsweep(j).getepoch('Name', 'Epoch B').getap(4:end).freq]);
       end
   end
   [~, maxfreqswp]=max(freqs);
   if ~isempty(maxfreqswp)
   obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).getepoch('Name', 'Epoch B').Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).getepoch('Name', 'Epoch B').Time(end)+300])
   end
   title('MaxFreqSweep')
   
   filename=strsplit(filelist{i}, {'.', '\'});
   saveas(gcf, ['D:\Connectivity analysis\AnalysisRene\Clustering Morphys\CCStepAnalysisSnapshots2\MatlabFigs\', filename{end-1}, '.fig'])
   print(['D:\Connectivity analysis\AnalysisRene\Clustering Morphys\CCStepAnalysisSnapshots2\' filename{end-1} '.jpg'],  '-djpeg', '-r300');
      
end