function protocolset = xMakeAnalogwaveformtabs(destinationpath)
    % Make Analogwaveform objects using info from an excel file.
    % See excel file specified below (path_in) for how to setup such a file
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
    %% get eCode info from excel doc
    path_in  = 'D:\Morphys\Data\Electrophysiology\Protocols\AnalogWaveforms\RWS\AllenEpochList.xlsx';
    if nargin == 0, destinationpath=''; end 
    %% load it and create struct
    ecodetable   = readtable(path_in);
    procolumns   = {'protocolName','protocolGuid','nrOfSweeps','sampleFreq'};
    epochcolumns = ecodetable.Properties.VariableNames(~ismember(ecodetable.Properties.VariableNames,[procolumns,'comments','fmt']));
    [ecodepros,~,IC] = unique(ecodetable(:,procolumns),'rows','stable');
    ecodestruct      = table2struct(ecodepros);
    
    %% make epochs
    uniqueidx = unique(IC);
    for i = 1:numel(uniqueidx)
        ecodestruct(i).epochs = ecodetable(IC==uniqueidx(i),epochcolumns);
    end
    
    %% make protocolObjects
    protocolset = arrayfun(@(x) Analogwaveformtab(ecodestruct(x).epochs,'name',        ecodestruct(x).protocolName,...
                                                                        'guid',        ecodestruct(x).protocolGuid,...
                                                                        'samplefreq',  ecodestruct(x).sampleFreq,  ...
                                                                        'nrofsweeps',  ecodestruct(x).nrOfSweeps), ...
                                                                        1:numel(ecodestruct),'UniformOutput',false);
    
    protocolset = [protocolset{:}];    
    
    %% save protocolObjects
    if ~isempty(destinationpath) && ischar(destinationpath), 
        arrayfun(@(x) protocolset(x).saveme(destinationpath,1),1:numel(protocolset),'UniformOutput',false); 
    else
        disp('no destination path provided, nothing saved')
    end
        
      
end