% -------------------------------------------------------------------
%  Generated by MATLAB on 10-Sep-2021 16:59:21
%  MATLAB version: 9.10.0.1602886 (R2021a)
% -------------------------------------------------------------------

p   = ('/Volumes/Ephys2/Ephys hipp/healthy/cluster1/chirp');
fn  = 'hipp4_2017_03_29_0312.abf';
fp  = fullfile(p,fn);
ss  = load('/Users/elinemertens/Documents/CNCR/Data/Matlab Data/Morphys-master/Matlab/Abfanalysis/Setupsettings_INF.mat');
ss  = ss.obj;
zabf = Abffile(fp,ss);
%%
signal =  zabf.channels(2).analogins(1, 1).sweeps(1:3).avtrace ;
plot(signal)

%% Thijs werkend
signal = signal.low_pass_filter(1000);
% data = abfload_pro(fp);
% stim = squeeze(data(:,2,1));
signal=signal.resample(0:0.1:signal.TimeInfo.End).Data;
stim = zabf.channels(1).analogins(2).sweeps(1).Data   ;
% signal = signal.Data  ;
% time = zabf.channels(2).analogins(1, 1).sweeps(1) ;
% time = time.low_pass_filter(1000) ; 
% time=time.resample(0:0.1:time.TimeInfo.End).Time;
%%
 Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;
    plot(signal)
    grid off
                 ylabel('Voltage (mV)')
              xlabel('Time (sec)')  
    signal_fft=fft(signal);
    stim_fft=fft(stim);

%% raw
    figure(1)
    plot(f, abs((signal_fft)./abs(stim_fft)));
    xlim([1.0 25])
    hold on
 %% smoothed
   figure(2)
   smoothed=smooth(abs((signal_fft)./abs(stim_fft)),23,  'lowess');
    plot(f, smoothed);
    xlim([1.0 25])
    % plot 3dB 
   [max_imp, loc] = max(smoothed(f>1.2 & f<25));
    f2 = f(f>1.2 & f<25);
    res_freq = f2(loc);
    max_imp_mOhm = (1000 * max_imp) ; 
    
   % 3db cutoff
   normalized = smoothed(f>1.2 & f<25)./max_imp;
   % when does normalized frequency response go below square root of 0.5?
   % This is the definition of 3 dB cutoff!
   cutoff_3db=f2(find(normalized>=sqrt(0.5), 1, 'last'));
   figure (3)
   plot(f2, normalized)
   hold on
   line(xlim,[sqrt(0.5), sqrt(0.5)],'linewidth',2,'lineStyle',':','color','k')
   %line([cutoff_3db, cutoff_3db],ylim,'linewidth',2,'lineStyle',':','color','k')
 ylim([0.1 1.2])
 set(gca, 'TickDir', 'out')


%%

for i=1:numel(fn)
    data=abfload_pro(fn(i));
 
   
    signal=mean(squeeze(data(:,1,:))');
    stim=mean(squeeze(data(:,2,:))');

     Fs = 10000;           % Sampling frequency
    L = length(signal);      % Signal length
    f = (0:L-1)*Fs/L;

    signal_fft=fft(signal);
    stim_fft=fft(stim);

    %plot

%     %raw
    plot(f, abs(signal_fft)./abs(stim_fft));
    xlim([0.5 35])
    hold on
        %smoothed
    smoothed=smooth(abs(signal_fft)./abs(stim_fft),23,  'lowess');
    plot(f, smoothed);
    xlim([0.5 35])
    
    %resonance freq
    [~, loc] = max(smoothed(f>0.5 & f<35));
    f2 = f(f>0.5 & f<35);
    res_freq = f2(loc);

end