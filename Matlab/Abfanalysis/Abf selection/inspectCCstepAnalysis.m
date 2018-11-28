
% load all analyzed ccstep files and create analysis overviews for
% inspection
load('D:\Human inhibition\Cell classification\AllenClusterMethod\superclusterdata.mat')
files={SuperCellSummary.File};
files=ChangeFilenameExts(files, 'mat');
users={SuperCellSummary.UserID};

for i=1:numel(files)
   fprintf('Open file %1.0f out of %1.0f \n', i, numel(files))
    if strcmp(users{i}, 'NGA') || contains(files{i}, '2017_08_30')
        load(['D:\Morphys\Data\Electrophysiology\Abffiles\NAG\abf\Analyzed\' files{i}]);
    elseif strcmp(users{i}, 'RWS')
        load(['D:\Morphys\Data\Electrophysiology\Abffiles\RWS\abf\Analyzed\' files{i}]);
    else
        error('User not found');
    end
   
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
   
   filename=strsplit(files{i}, {'.', '\'});
   saveas(gcf, ['D:\Human inhibition\Cell classification\ReneCCstepSnapshots\MatlabFigs\', filename{end-1}, '.fig'])
   print(['D:\Human inhibition\Cell classification\ReneCCstepSnapshots\' filename{end-1} '.jpg'],  '-djpeg', '-r300');
      
end