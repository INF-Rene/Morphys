%select folder to analyze 
folder = uigetdir 
cd (folder);
%make a list with all nwbs
list = dir();
list = struct2table(list);
list = list(list.bytes>1000,:); %only files with actual data

%% USE THIS nwb2: 
for i = 1:numel(list.name)
 nwb = NWBfile([{'X8'}]);
 nwb.savename = sprintf('res_%s.nwb',nwb.filename(1:end-4));
 saveme(nwb,'/Users/elinemertens/Data/ephys/Hippocampus/Resonance/198', nwb.savename) 
end


%% Resonance analysis  single files 

fn= '/Users/elinemertens/Data/ephys/Hippocampus/165/H19.29.165.21.41_ephys_converted/H19.29.165.21.41.03-nwb2.nwb'

nwb = NWBfile(fn,[{'X8_C'}])


%% Zorg dat hiervoor de traces individueel al bekeken zijn of ze vergelijkbaar zijn (indien niet, geen average trace nemen)

%% know the list (only for list)
basedir = '/Users/elinemertens/Data/ephys/Hippocampus/Resonance' ;
fileinfo  = dir(fullfile(basedir,'*.nwb'));
filelist  = {fileinfo.name};

%% Bij errors 'expected input to be finite' moet je 
% even de losse sweeps bekijken en dan de juiste sweeps selecteren
% signal = nwb.getstimset.getnwbchannel.getsweep(1:3).avtrace
  
signal = nwb.getstimset.getnwbchannel.getsweep(2:4).avtrace;

plot(signal)
%%
signal = signal.low_pass_filter(1000);
stim = nwb.getstimset.getnwbchannel.getstimwave(2).Data ;
% if the first sweep is bad or incomplete, take getstimwave(2)
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;


%%
%   signal=normalize(mean(squeeze(data(:,1,:))'), 'range');
%     stim=normalize(mean(squeeze(data(:,2,:))'), 'range');

    %%
    Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;

    plot(signal)
    grid off
    
    %%
    signal_fft=fft(signal);
    stim_fft=fft(stim);

    %plot

%     %raw
    plot(f, abs((signal_fft)./abs(stim_fft)*1000));
    xlim([0.5 35])
    hold on
        %smoothed
    smoothed=smooth(abs((signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 35])
    
    %% plot only smoothed 
    
     smoothed=smooth((abs(signal_fft)./abs(stim_fft)*1000),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 35])
    grid off   
            
                 ylabel('MOhm')
              xlabel('Frequency (Hz)')
    
   %% resonance freq
    [~, loc] = max(smoothed(f>0.5 & f<35));
    f2 = f(f>0.5 & f<35);
    res_freq = f2(loc);
    

