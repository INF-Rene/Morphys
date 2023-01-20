%%Load in analyzed mat files 
basedir = '/Volumes/Ephys2/Ephys hipp/healthy/cluster3' ;
savedir = '/Users/elinemertens/Data/ephys/Summary/Human' ;
savename = 'Summary_209' ;

%% load file list
fileinfo  = dir(fullfile(basedir,'*.nwb'));
filelist  = {fileinfo.name};

%% settings
filter_freq = 1000; % low pass filter frequency for sag traces (Hz)

%% loop over files
%t_all=table();
table_initialize = false ;

for i=1:numel(filelist)
    %load file
    fn = fullfile(basedir, filelist{i});
    nwb = NWBfile(fn,[{'CHIRP'}]) ;
   
     signal = nwb.getstimset.getnwbchannel.getsweep(1:3);
 plot(signal) ; 
 
   signal = nwb.getstimset.getnwbchannel.getsweep(1:3).avtrace;
plot(signal) ; 
        
%check if you take correct stimwave(#) (after wash in this is diff)
stim = nwb.getstimset.getnwbchannel.getstimwave(2).Data ;
signal = signal.low_pass_filter(1000);
% if the first sweep is bad or incomplete, take getstimwave(2)
signal= signal.resample(0:0.1:signal.TimeInfo.End).Data;
time = nwb.getstimset.getnwbchannel.getsweep(2).Time  ; 
time= time.low_pass_filter(1000);
time= time.resample(0:0.1:time.TimeInfo.End).Data;

%signal = lowpass(signal,1000,10000);

  t = table(stim', signal', 'VariableNames', {'stim', 'signal'});
      %  t.Protocol = repmat(stimset.name, height(t), 1);
      t.Filename = repmat(filelist{i}, height(t), 1);
     
        
        %append results to overview
        if table_initialize 
            t_all = [t_all; t];
        else 
            t_all = t ; 
            table_initialize = true
        end
           
    end
    