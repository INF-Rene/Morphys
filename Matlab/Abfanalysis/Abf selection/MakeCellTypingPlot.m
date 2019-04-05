
% load all analyzed ccstep files and create analysis overviews for
% inspection
% load('D:\Human inhibition\Cell classification\AllenClusterMethod\superclusterdata.mat')
% files={SuperCellSummary.File};
% files=ChangeFilenameExts(files, 'mat');
% users={SuperCellSummary.UserID};
destinpath='D:\Human inhibition\Cell classification\AllenClusterMethod\Snaps\';
%files={'2018_02_07_0000.abf';'2018_02_07_0001.abf';'2018_02_07_0004.abf';'2018_02_07_0005.abf'};
%files=fullfile('D:\ConnectivityABF\Batch 2\2018-02-07', files);
files={'2015_11_04_0095.abf';'2016_03_16_0007.abf';'2016_06_22_0002.abf';'2016_06_22_0001.abf';'2016_09_14_0060.abf';...
    '2016_10_05_0003.abf';'2016_10_05_0001.abf';'2016_10_05_0108.abf';'2016_10_05_0104.abf';...
    '2017_01_25_0003.abf';'2017_01_25_0004.abf';'2017_02_01_0002.abf';'2017_08_30_0001.abf';...
    '2017_11_01_0055.abf'};
files=fullfile('D:\Morphys\Data\Electrophysiology\Abffiles\NAG\abf\', files);
users={'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA', 'NGA'};
ss=load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
ss=ss.obj;

for i=14:numel(files)
    fprintf('Open file %1.0f out of %1.0f \n', i, numel(files))
%     if strcmp(users{i}, 'NGA') || contains(files{i}, '2017_08_30')
%         a=load(['D:\Morphys\Data\Electrophysiology\Abffiles\NAG\abf\Analyzed\' files{i}]);
%         a=a.obj;
%     elseif strcmp(users{i}, 'RWS')
%         a=load(['D:\Morphys\Data\Electrophysiology\Abffiles\RWS\abf\Analyzed\' files{i}]);
%         a=a.obj;
%     else
%         error('User not found');
%     end
    a = Abffile(files{i}, ss);
    if strcmp(a.getchannel.getin('signal', 'primary').units, 'pA')==1
        error('File has wrong units. If it is a current clamp recording, you can adjust the units in clampfit')
    end
    
    a = a.analyseabf;
    
    delta=[a.getchannel.getin('signal','primary').getsweep(2).getepoch.amplitude]-[a.getchannel.getin('signal','primary').getsweep(1).getepoch.amplitude];
    idxstr=a.getchannel.getin('signal','primary').getsweep(2).getepoch(delta>0).idxstr;
    
    steps=[a.getchannel.getin('signal','primary').getsweep.getepoch('idxstr',idxstr).amplitude];
    frstap=find([a.getchannel.getin('signal','primary').getsweep(steps>0).getepoch('idxstr',idxstr).nrofaps]>0,1)+find(steps>=0, 1);
    rheo=a.getchannel.getin('signal','primary').getsweep(frstap).getepoch('idxstr',idxstr).amplitude;
    epstart=a.getchannel.getin('signal','primary').getsweep(1).getepoch('idxstr',idxstr).TimeInfo.Start;
    epend=a.getchannel.getin('signal','primary').getsweep(1).getepoch('idxstr',idxstr).TimeInfo.End;
    
    close all
    figure ('position', [0 0 1920 1000])
    
    %passive plot
%     subplot(2,4,1);
    subplot('position', [0.04 0.6 0.2 0.35 ])
    a.getchannel.getin('signal', 'primary').getsweep(steps<0).plotanalysis
    xlim([epstart-100 epstart+1200])
    title('Passive properties')
    
    %First AP
    subplot('position', [0.28 0.6 0.2 0.35 ])
    ap=a.getchannel.getin('signal', 'primary').getsweep(frstap).getepoch('idxstr',idxstr).getap(1);
    ap.plotanalysis2;
    xlim([ap.start_time-1 ap.start_time+8])
    title('First AP')
    
    %I-F plot
    subplot('position', [0.54 0.6 0.19 0.35 ])
    freqs=NaN(1,a.nrofsweeps);
    nrofAPs=NaN(1,a.nrofsweeps);
    for j=frstap:a.nrofsweeps
        freqs(j)=nanmean([a.getchannel.getin('signal', 'primary').getsweep(j).getepoch('idxstr',idxstr).getap(4:end).freq]);
        nrofAPs(j)=a.getchannel.getin('signal', 'primary').getsweep(j).getepoch('idxstr',idxstr).nrofaps;
    end
    yyaxis left;
    plot(steps./steps(frstap)*100,freqs, '--o')
    xlabel('I/Ithresh(%)')
    ylabel('Mean Inst. Frequency of APs 4-end (Hz)')
    ylim([0 max(freqs)*1.1])
    yyaxis right;
    plot(steps./steps(frstap)*100,nrofAPs, '--o')
    ylabel('Nr of APs')
    ylim([0 max(nrofAPs)*1.1])
    title('I vs steady-state firing rate')
    
    %adaptation & bursting behaviour plot
    subplot('position', [0.77 0.6 0.19 0.35 ])
    tmp=steps(steps>rheo);
    brst=NaN(1,numel(tmp));
    adapt=NaN(1,numel(tmp));
    for j=1:numel(tmp)
        isis=[a.getchannel.getin.getsweep(steps==tmp(j)).getepoch('idxstr',idxstr).getap.isi];
        isis=isis(~isnan(isis));
        if numel(isis)>2
            brst(j)=nanmean(isis(2:end))/isis(1);
            disi=diff(isis);
            adapt(j)= nanmean(disi./(isis(1:end-1)+isis(2:end)));
        end
    end
    yyaxis left;
    plot(tmp./tmp(1)*100,brst, '--o')
    xlabel('I/Ithresh(%)')
    ylabel('Burst Index')
    ylim([0 max(brst)*1.1])
    yyaxis right;
    plot(tmp./tmp(1)*100,adapt, '--o')
    ylabel('Adaptation index')
    ylim(sort([0 max(adapt)*1.1]))
    title('Bursting & Adaptation')
    
    
    %Sweep 1, 2, 5
    subplot('position', [0.04 0.14 0.2 0.35 ])
    a.getchannel.getin('signal', 'primary').getsweep(frstap).plotanalysis
    xlim([epstart-100 epstart+1100])
     title(['Rheobase (' num2str(rheo) ' pA)'])
    
    [~,rheo50]=min(abs(steps-rheo-50));
    subplot('position', [0.28 0.14 0.2 0.35 ])
    a.getchannel.getin('signal', 'primary').getsweep(rheo50).plotanalysis
    xlim([epstart-100 epstart+1100])
    title('Rheobase + 50')
    
    [~,rheo100]=min(abs(steps-rheo-100));
    subplot('position', [0.54 0.14 0.2 0.35 ])
    a.getchannel.getin('signal', 'primary').getsweep(rheo100).plotanalysis
    xlim([epstart-100 epstart+1100])
    title('Rheobase + 100')
    
    
    %at least 10 APs (or maximum)
    subplot('position', [0.78 0.14 0.2 0.35 ])
    AP10=find(nrofAPs>=10,1);
    if isempty(AP10), [~,AP10]=max(nrofAPs); end
    a.getchannel.getin('signal', 'primary').getsweep(AP10).plotanalysis
    xlim([epstart-100 epstart+1100])
    title(['first 10 AP (' num2str(steps(AP10)) ' pA)'])
    
    filename=strsplit(files{i}, {'.', '\'});
    saveas(gcf, [destinpath, users{i},'_', filename{end-1}, '.fig'])
    print([destinpath, users{i},'_' filename{end-1} '.jpg'],  '-djpeg', '-r300');
    
end