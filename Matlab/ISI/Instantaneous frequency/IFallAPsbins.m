%BULK loading
% set path to find .NWB files
dir_abfs          = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_nAChRs currents/Human_ephys/test_IF';
% set path to save analyzed files
dir_mats_analysed = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_nAChRs currents/Human_ephys/test_IF'
% get a list of all .nwb files in the specified path
fileinfo  = dir(fullfile(dir_abfs,'*.nwb'));
filelist  = {fileinfo.name};

%% 
%%how to load specific stimsets of the file - stimsets containing APs (and
%%supposed to contain APs)

fn = '/Users/annagalakhova/PhD INF CNCR VU/DATA/L1_nAChRs currents/Human_ephys/test_IF/H21.29.201.11.11.01.nwb'

nwb = NWBfile(fn,[{'Search'} {'teps'}])
%%

    % nr of APs first sweep (Eline Edit)
         NrofAPfrstSwp = length(sweep.ap);
         NrofAPtrainSwp = length(sweep(TrainSweep).ap);
         NrofAPlastSwp = length(sweep(NrofSweeps).ap);
   
        if length(sweep(frstspikeswp).ap) > 1
            isis_FS = [sweep(frstspikeswp).ap(2:end).isi];
        else 
            isis_FS = NaN ;
        end

        if length(sweep(NrofSweeps).ap) > 1
            isis_LS = [sweep(NrofSweeps).ap(2:end).isi];
            isis_LS1 = [sweep(NrofSweeps).ap(2).isi];
        else 
            isis_LS = NaN ;
            isis_LS1 = NaN 
        end
    %%
    
    