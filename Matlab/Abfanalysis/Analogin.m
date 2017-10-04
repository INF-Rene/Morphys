classdef Analogin < Sharedmethods
    
    %   ANALOGIN Analog-to-Digital conversion channel
    %   An Analogin is a recording channel. For ABF files from the INF, this usually is the output of a Primary, Secondary, 
    %   Scope or Mode BNC of an amplifier. Master 8 triggers are also often recorded. 
    %
    %   Still thinking about whether to let this channel inherit from timeseriescollection class... Conceptually, an Analogin
    %   is a timeseries collection (every sweep is one timeseries in the collection...). But what extra methods would this
    %   offer? The addts and removets methods could be extended for addSweep/removeSweep 
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################

    properties (Access = public)
        guid_channel          % guid of 'parant' channel
        number                % ADC number (corresponds to number of digidata analog input)(h.ADCSec.nADCNum)
        adcusername           % ADC name as specified in pClamp lab bench and selected as input on protocol inputs tab)
        units                 % units of signal recorded in this ADC channel
        telegraphenabled      % boolean indication of whether telegraphed instrument was enabled (h.ADCSec.nTelegraphEnable) 
        telegraphinstrument   % number of telegraphed instrument (h.ADCSec.nTelegraphInstrument)
        gain                  % gain setting (h.ADCSec.fTelegraphAdditGain)
        lowpassfilter         % low-pass filter (h.ADCSec.fTelegraphFilter)
        instrumentscalefactor % instrument scalefactor (Note this is different from the "eCode"-related scalefactor!) (h.ADCSec.fInstrumentScaleFactor)
        protocollag           % lag in protocol. 
        samplefreq            % sampling frequency (Hz)
        nrofsweeps = 0        % number of sweeps in analog IN
        signal                % Type of signal recorded by this in (e.g. Primary, Secondary, ... )
        updatephase = 1       % phase of class; (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)   
        sweeps                % list of Sweeps (= Sweep class instances)
        
    end
    
%####################################################### METHODS ############################################################
    methods
        
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Analogin(infostruct,data,guid_channel)
            % This function adds all properties to Analog IN channel. Most derive from the "h.ADCSec" section from the "h" 
            % struct in output of abfload_pro.
            % ininfo is a struct with the following fields:
            % number                  
            % adcusername            
            % units                  
            % samplefreq             
            % telegraphenabled       
            % telegraphinstrument    
            % gain                   
            % lowpassfilter          
            % instrumentscalefactor
            % --------------------
            
            % call to superclass constructor
            obj = obj@Sharedmethods;
            
            % return empty if no input arguments
            if nargin == 0, return; end
            
            obj.guid_channel = guid_channel;
            if isstruct(infostruct)
                flds = fieldnames(infostruct);
                flds = flds(~ismember(flds,{'h'}));
                for i=1:numel(flds)
                    obj = obj.set(flds{i},infostruct.(flds{i}));
                end
            else
                error('Please provide an info struct')
            end
            obj.nrofsweeps              = size(data,2);
            
            obj.updatephase             = 1; % default
            obj.protocollag             = 0; % default
            
            % Populate sweepList with Sweep objects
            for i=1:obj.nrofsweeps
                obj = obj.addsweep(i,infostruct.h,data(:,i));
            end
            
        end
        
        % ------------------------------------------- HELPER METHODS --------------------------------------------------------
        function identify(obj)
            % This function returns a small summary of the Analogin object.
            if ~isscalar(obj), error('Analogin object must be scalar.'); end
            if isempty(inputname(1)), myname = 'This'; else myname = sprintf('"%s"',inputname(1)); end
            idmessage = sprintf(['%s is an Analogin object:\n'...
                                 'IN name:      %s\n',...
                                 'IN number:    %d\n',...
                                 'Units:        %s\n'...
                                 'Sweeps:       %d\n',...
                                 'Sample freq:  %d kHz\n',...
                                 'Lowpass filt: %d kHz\n'],...
                                 myname, obj.adcusername, obj.number, obj.units, ...
                                 obj.nrofsweeps, floor(obj.samplefreq/1000),floor(obj.lowpassfilter/1000));
             disp(idmessage)
        end
        
        function obj = addsweep(obj,idx,h,swpdata)
            % Instantiates one Sweep object using info in the initialisation struct and appends
            % this to the Analogin objects sweepList.
            % --------------------
            if ~isscalar(obj), error('Analogin object must be scalar.'); end
            if ~isscalar(idx), error('Sweep index must be scalar.'); end
            s = Sweep('number',idx,'guid_IN',obj.guid);
            s = s.adddata(swpdata,obj.units);
            s = s.addtime(obj.samplefreq);
            s = s.adddates(h);
            if isempty(obj.sweeps), obj.sweeps        = s;
            else                    obj.sweeps(end+1) = s;
            end
            obj = updatesweepstats(obj);
        end

        function obj = sortsweeps(obj,varargin)
            % sort SWEEPs in list
            % Returns object with list sorted according to sweep FEATURE. 
            % Default feature is number, default order is ascending, unless DESCENDING is set to 1.
            % Usage
            %   obj = sortsweeps(OBJ,FEATURE,DESCENDING)
            %
            % Example:
            %   obj.sortsweeps()               --> sort aps in ascending order by number
            %   obj.sortsweeps('nrofaps')      --> sort aps in ascending order by amplitude
            %   obj.sortsweeps('guid',1)       --> sort aps in descending order by their guid 
            %
            % See also SWEEP, SORTITEM, SORT
            for i = 1:numel(obj), obj(i) = obj(i).sortitems('sweeps',varargin{:}); end
        end
        
        function ap = getsweep(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'sweeps',0,varargin{:});
        end
        
        function obj = removesweep(obj,varargin)
            % remove SWEEP(s) from list. 
            %
            % see also SWEEP, GETITEM, REMOVEITEM, SELECTEPOCH, SHAREDMETHODS
            for i = 1:numel(obj), obj(i) = obj(i).removeitem('sweeps',varargin{:}).updatesweepstats; end
        end
        
        function obj = selectsweep(obj,varargin)
            % select SWEEP(s) from list. 
            %
            % see also SWEEP, GETITEM, SELECTITEM, REMOVEEPOCH, SHAREDMETHODS
            for i = 1:numel(obj), obj(i) = obj(i).selectitem('sweeps',varargin{:}).updatesweepstats; end
        end
        
        function data = getdata(obj)
            % function to retrieve all channel data (now stored downstream in sweeps). Data returned as 2D matrix (samples x
            % sweeps)
            if isscalar(obj)
                data = [obj.getsweep.getdata];
                if iscell(data), data = cell2mat(data); end
            else
                data = arrayfun(@(x) obj(x).getdata,1:numel(obj),'UniformOutput',false);
                data = cell2mat(data);
                data = reshape(data,size(data,1),numel(obj),obj(1).nrofsweeps);
            end
        end
        
        % ------------------------------------------ UPDATING METHODS -------------------------------------------------------
        function obj = analysein(obj)
            % perform default analysis of analog in. Only for primaries with voltage traces are analysed. 
            for i = 1:numel(obj)
                fprintf('\t\tanalysing analogin %d\n',obj(i).number)
                if strcmp(obj(i).signal,'primary') && strcmpi(obj(i).units,'mV')
                    obj(i).sweeps = obj(i).getsweep.analysesweep;
                end
                obj(i).updatephase = 3;
            end
        end
        
        function obj = updatesweepstats(obj)
            if ~isscalar(obj), error('Analogin object must be scalar.'); end
            obj.nrofsweeps = numel(obj.sweeps);
        end
        
        function obj = makeepochs(obj,analogwaveformtable,addlag)
            % updates the sweeps of ANALOGIN with epoch information contained in Analogwaveformtab object
            % (analogwaveformtable). If addlag ==1 , epochs will be delayed with the annoying pClamp delay (1/64th of trace 
            % length delay)
            if ~isscalar(obj), error('Analogin object must be scalar.'); end
            if addlag == 1, 
                obj.protocollag = floor(size(obj.getdata,1)/64); 
            end                    

            for i = 1:obj.nrofsweeps
                swptab = analogwaveformtable.epochs4sweep(i);                       % make table of epochs for this sweep              
                obj.sweeps(i) = obj.getsweep(i).makeepochs(swptab,obj.protocollag); % replace sweep in list with updated one with populated epoch list
            end
            obj.updatephase = 2;
        end

        % ------------------------------------------ PLOTTING METHODS -------------------------------------------------------
        function plot(obj,varargin)
            % plot the ANALOGIN object. Note that when ANALOGIN object is not scalar, all will be plot in same figure panel!
            if strcmp(obj.signal,'secondary'), varargin = cat(2,varargin,{'linewidth',2}); end
            for i=1:numel(obj)
                hold on
                for ii=1:obj(i).nrofsweeps
                    obj(i).getsweep(ii).plot(varargin{:})
                end
                hold off
            end
        end
        
        function plotanalysis(obj)
            % plot ANALOGIN with analysis details (e.g. baselines, ap peak events, tau fits and sags etc...). Note to plot
            % these, the ANALOGIN must be analysed first with the ANALYSEIN method.
            % Note that when ANALOGIN object is not scalar, all will be plot in same figure panel!
            % See also ANALYSEIN
            for i = 1:numel(obj)
                hold on
                for ii=1:obj(i).nrofsweeps
                    obj(i).getsweep(ii).plotanalysis
                end
                hold off
            end
        end
    end

end

