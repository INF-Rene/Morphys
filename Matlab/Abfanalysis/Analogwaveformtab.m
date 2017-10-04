classdef Analogwaveformtab < Sharedmethods
    %ANALOGWAVEFORMTAB Specification of an analog waveform
    %   Analog waveform object is an object specifying the analog waveform signal sent out via an OUT (DAC) channel. It 
    %   contains the same information as the analog waveform tab seen when editing a pClamp protocol. AnalogWaveform objects
    %   cannot inherit from the Matlab table class (Sealed unfortunately). 
    %   Analogwaveformtab objects can be saved as atf files and used as stimulus files in pClamp.
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
     
%###################################################### PROPERTIES ##########################################################
    properties
        name        = '';    % name of protocol/atf
        samplefreq  = [];    % default sampling frequency for this protocol.
        nrofsweeps  = [];    % default number of sweeps for this protocol
        tab         = table; % table with epoch specification (analogous to pClamp Analog Waveform tab)
        scalefactor = 1;     % scalefactor
        units       = 'pA'   % important!! default is current clamp protocols!       
    end

%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Analogwaveformtab(t,varargin)
            % create an ANALOGWAVEFORMTAB object. 
            % First input must be table, rest optional Name/Value pairs.            
            obj = obj@Sharedmethods;
            
            if nargin == 0, return;
            elseif ~istable(t), 
                error('First input to AnalogWaveform must be a table.')
            else
                % fix variable names
                t.Properties.VariableNames = lower(t.Properties.VariableNames);
                obj.tab = t;
            end
            
            % check name/value pairs
            if mod(numel(varargin),2)~=0
                error('Uneven number of name/value pairs.')
            elseif ~all(cellfun(@ischar,varargin(1:2:end))) % names should all be strings of characters
                error('All property names must be strings')
            else
                for i=1:2:numel(varargin)
                    obj = obj.set(varargin{i},varargin{i+1});
                end
            end
        end
        
        % ------------------------------------------- OTHER METHODS ---------------------------------------------------------    
        function identify(obj)
            idmssg = sprintf('I am the Analog Waveform object describing %s.',obj.name);
            disp(idmssg)
        end
        
        function swptab = epochs4sweep(obj,idx)
            % function to make a table of epochs for a particular sweep number (sweepNo) using obj.tab. obj and idx must be 
            % scalar. So this function uses the index of the sweep provided to calculate all the (incremented) amplitudes and 
            % durations of epochs.
            if ~isscalar(obj),     error('Object must be scalar')
            elseif ~isscalar(idx), error('Index must be scalar.')
            end
            swptab = obj.tab;
            swptab.firstlevel = swptab.firstlevel + swptab.deltalevel*(idx-1);
            
            epdur = cellfun(@(x) eval(x),swptab.timespan,'UniformOutput',false);
            if all(cellfun('length',epdur)==1)
                swptab.timespan = cell2mat(epdur) + swptab.deltaduration *(idx-1);
            else
                swptab.tmp = zeros(height(swptab),1);
                for i=1:height(swptab)
                    if numel(epdur{i})==1, swptab{i,'tmp'} = cell2mat(epdur(i)) + swptab{i,'deltaduration'} *(idx-1);
                    else
                        swptab{i,'tmp'} = epdur{i}(idx);
                    end
                end
                swptab.timespan = []; % ditch old
                swptab.Properties.VariableNames{'tmp'} = 'timespan'; % replace with new
            end
            
            % remove obsolete 
            swptab.deltalevel    = [];
            swptab.deltaduration = [];
        end
        
        function wf = waveform(obj,idx,lag)
            % function generates a waveform using tab and sweep index. If index (idx) is scalar, wf is a timeseries object. 
            % If idx is a vector, wf is a cell array of timeseries objects. lag is optional, specifies number of points delay 
            % before protocol start (1/64th of sweep length in samples). This way, the waveform can be synchronised with the
            % recorded signals.
            error('The Analogwaveformtab waveform function is not yet fully operational, seems to be a problem with setting proper time units!')
            if nargin == 2, lag = 0; end

            if isscalar(idx) 
                s  = Sweep('samplefreq',obj.samplefreq);    % make temporary Sweep object
                e  = obj.epochs4sweep(idx);                 % obtain epoch table for that 
                e.firstlevel = e.firstlevel*obj.scalefactor;% scale amplitudes with scalefactor
                s  = s.updatesweep(e,lag);                  % update sweep with epochs, including lag
                wf = s.waveform;                            % get waveform by using Sweep waveform method  
                wf.DataInfo.Units = obj.units;              % add units
                wf.Name = sprintf('sweep_%d',idx);          % add name
            elseif isvector(idx)
                wf = arrayfun(@(x) waveform(obj,x,lag),idx,'UniformOutput',false); 
            else
                error('idx must be scalar or vector')
            end
        end
        
        % ------------------------------------------- OTHER METHODS ---------------------------------------------------------
        function atfme(obj,destinationpath)
            % function to save an AnalogWaveform object as an axon text file ('.atf'), which can be injected into a cell by 
            % using it as a 'stimulus file' in your pClamp protocol. 
            % Set properties of Analogwaveformtab object to tailor atf file to your needs:
            % obj.set('scalefactor',2,'samplefreq',20000,'nrofsweeps',5).atfme will print an atf file with 5 sweeps, sampled 
            % at 20kHz and a scalefactor of 2. Filename and folder can be changed by setting the 'name' and 'atfdestination'
            % properties to desired values, respectively.
            
            if ~isscalar(obj), error('object must be scalar'); end
            
            % make a matrix of time and sample points, scaled accoding to scalefactor
            swpstartms = round((0.000:obj.nrofsweeps-1)*obj.waveform(1).TimeInfo.End,3);
            wf_array   = obj.waveform(1:obj.nrofsweeps);
            mtx        = cat(2,wf_array{1}.Time, cell2mat(arrayfun(@(x) wf_array{x}.Data,1:numel(wf_array),'UniformOutput', false))*obj.scalefactor);

            % To let pClamp understand the protocol in the axon text file, let's make the atf header lines (the 'special' 
            % pClampy bit that makes an atf) as shown in following example:

            % ATF	1.0
            % 7	13     
            % "AcquisitionMode=Fixed-Length Event-Driven"
            % "Comment="
            % "YTop=200.004"
            % "YBottom=-300.004"
            % "SweepStartTimesMS=0.000,1400.500,2801.000,4201.500,5602.000,7002.500,8403.000,9803.500,11204.000,12604.500,14005.000,15405.500"
            % "SignalsExported=Signal 00"
            % "Signals="	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"	"Signal 00"
            % "Time (s)"	"Trace #1 (pA)"	"Trace #2 (pA)"	"Trace #3 (pA)"	"Trace #4 (pA)"	"Trace #5 (pA)"	"Trace #6 (pA)"	"Trace #7 (pA)"	"Trace #8 (pA)"	"Trace #9 (pA)"	"Trace #10 (pA)"	"Trace #11 (pA)"	"Trace #12 (pA)"
            
            cellstr_starttimes = strjoin(cellstr(num2str(swpstartms(2:end)','%0.3f'))',',');
            cellstr_starttimes(ismember(cellstr_starttimes,' '))=[]; % remove empties
            signals_exported    = 'Signal 00';

            % zooming of yaxis (nonsense, but pClamp needs it so whatever)
            centre = max(mtx(size(mtx,1)+1:end))-range(mtx(size(mtx,1)+1:end))/2;
            Ytop   = centre+range(mtx(size(mtx,1)+1:end))/2*1.2;
            Ybot   = centre-range(mtx(size(mtx,1)+1:end))/2*1.2;

            % now collect all to make full header section
            headers = {sprintf('ATF\t1.0');
                       sprintf('7\t%d',obj.nrofsweeps+1);
                       '"AcquisitionMode=Fixed-Length Event-Driven"';
                       '"Comment="';
                       sprintf('"YTop=%d"',   Ytop);
                       sprintf('"YBottom=%d"',Ybot);
                       sprintf('"SweepStartTimesMS=0.000,%s"',cellstr_starttimes);
                       '"SignalsExported=Signal 00"';
                       ['"Signals="' repmat(sprintf('\t"%s"',signals_exported),1,obj.nrofsweeps)]
                       sprintf('"Time (s)"%s',strjoin(arrayfun(@(x) sprintf('\t"Trace #%d (pA)"',x),1:obj.nrofsweeps,'UniformOutput', false)))
                       };

            % store header in a new text file
            fullpath = fullfile(destinationpath,obj.name);
            fid = fopen([fullpath '.txt'], 'wt');      
            for i = 1:size(headers,1)
                h = headers{i,:};
                fprintf(fid, '%s\n', h);
            end
            fclose(fid);

            % add signals to txt file and change extension to atf
            save([fullpath '.txt'],'mtx','-ascii', '-tabs','-append') % save as text
            fprintf('\nSaving as: %s... ',     [fullpath '.atf']);
            movefile([fullpath '.txt'],[fullpath '.atf']);    % rename to atf
            fprintf('Done.\n')
        end
                           
        function saveme(obj,destinationpath,overwrite)
            % save the AnalogWaveform object as a .mat file. Default setting is not to overwrite any existing files. Set
            % overwrite to 1 to do so anyway.
            if nargin < 2, overwrite = 0; end
            if ~isscalar(obj), error('AnalogWaveform object must be scalar.'); end
            if exist(fullfile(destinationpath,[obj.name '.mat']),'file') && ~overwrite
                warning('Filename already exists, saving terminated. To overwrite, set overwrite flag to 1.')
                return
            end
            save(fullfile(destinationpath,obj.name),'obj')
        end
    end
    
end

