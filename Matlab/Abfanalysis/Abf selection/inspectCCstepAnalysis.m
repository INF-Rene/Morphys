
% load all analyzed ccstep files and create analysis overviews for
% inspection
load('D:\Morphys\Data\Labbooks\NAG\MetadataABF\DataSummaryNAG2.mat')

filelist=squeeze(struct2cell(Summary));
filelist=filelist(1,:)';
filelist=cellfun(@(x) strsplit(x, '.'), filelist, 'UniformOutput', false);
filelist=cellfun(@(x) x{1}, filelist, 'UniformOutput', false);
filelist=strcat(filelist, '.mat');
filelist=strcat('C:\Users\DBHeyer\Documents\PhD\Human Database\hippocampus\Steps\analyzed\',filelist);

for i=1:numel(filelist)
   fprintf('Open file %1.0f out of %1.0f \n', i, numel(filelist))
   load(filelist{i});
   
   close all
   figure ('position', [15 45 1500 720])
   subplot(2,3,1);
   zeroswp=-Summary(i).FrstP/Summary(i).DeltaP+1;
   obj.getchannel.getin('signal', 'primary').getsweep(1:zeroswp-1).plotanalysis
   title('Passive properties')
   
   subplot(2,3,2);
   swp100=round((-Summary(i).FrstP-100)/Summary(i).DeltaP+1);
   obj.getchannel.getin('signal', 'primary').getsweep(swp100:zeroswp-1).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch(2).Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch(2).Time(1)+100])
   title('Tau')
   
   subplot(2,3,3);
   obj.getchannel.getin('signal', 'primary').getsweep(zeroswp).plotanalysis
   if ~isnan(Summary(i).vmbaseM)
   ylim([Summary(i).vmbaseM-4 Summary(i).vmbaseM+4])
   end
   title('Zero current injection')
   
   subplot(2,3,4);
   obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).getap(1).start_time-10 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).FrstSpikeSwp).getap(1).start_time+20])
   title('First AP')
   
   subplot(2,3,5);
   obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch(3).Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(Summary(i).TrainSwp).getepoch(3).Time(end)+300])
   title('Trainsweep')
   
   subplot(2,3,6);
   freqs=[];
   for j=1:obj.nrofsweeps
       if numel(obj.getchannel.getin('signal', 'primary').getsweep(j).getepoch(3).getap)>=4
           freqs(j)=mean([obj.getchannel.getin('signal', 'primary').getsweep(j).getepoch(3).getap(4:end).freq]);
       end
   end
   [~, maxfreqswp]=max(freqs);
   if ~isempty(maxfreqswp)
   obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).plotanalysis
   xlim([obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).getepoch(3).Time(1)-50 obj.getchannel.getin('signal', 'primary').getsweep(maxfreqswp).getepoch(3).Time(end)+300])
   end
   title('MaxFreqSweep')
   
   filename=strsplit(filelist{i}, {'.', '\'});
   saveas(gcf, ['C:\Users\DBHeyer\Documents\PhD\Human Database\hippocampus\Steps\figures\', filename{end-1}, '.fig'])
   print(['C:\Users\DBHeyer\Documents\PhD\Human Database\hippocampus\Steps\figures\' filename{end-1} '.jpg'],  '-djpeg', '-r300');
      
end