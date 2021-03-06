classdef Analogout < Sharedmethods
    
    %   ANALOGOUT Digital-to-Analog conversion channel
    %   An Analogout is a channel through which an analog waveform is sent out to the "command" BNC of an amplifier channel. 
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    properties (Access = public)
        guid_channel        % guid of 'parant' channel
        number              % (nDACNum) DAC number (corresponds to number of digidata analog input)
        dacusername         % OUT name specified in pClamp lab bench and selected as input on protocol inputs tab)
        units               % units of signal outputted by this DAC channel
        epochinfosource     % source of epoch info (epoch section or stringsection)
        analogwaveformtable % table specifying epoch information
        path2pro            % path to pClamp protocol file ('.pro') used
        protocolfile        % name of the pClamp protocol used
        path2stim           % path to stimulus file ('.atf') used
        stimulusfile        % name of the stimulus file used for analog waveform output (is empty if no atf file was used)
        scalefactor = 1     % scalefactor of protocol injected waveform. Applies to analog waveforms defined by ATF stimfiles (else scalefactor=1)
        holdingI = nan      % holding current determined from secondary analog in
        holdingV = nan      % holding voltage determined from secondary analog in
    end
       
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Analogout(infostruct,guid_channel)
            % This function adds all properties to OUT channel. 
            % outinfo is a struct containing following fields:
            % 'name'        = name of the OUT channel (e.g. "OUT 1" or "Cmd 0")
            % 'number'      = number of OUT (DAC) channel 
            % 'units'       = units, usually pA or mV
            % 'path2pro'    = path to pClamp protocol file ('.pro') used
            % 'protocolfile'= name of pClamp protocol file ('.pro') used
            % 'path2stim'   = path to stimulus file ('.atf') used
            % 'stimulusfile'= name of stimulus file ('.atf') used
            % 'analogwaveformtable' = Analogwaveformtab object containing info of analog waveform injected into cell, per epoch.
            % --------------------
            
            % call to superclass constructor
            obj = obj@Sharedmethods;
            
            % return empty if no input arguments
            if nargin == 0, return; end
            
            obj.guid_channel = guid_channel;
            
            obj.dacusername  = infostruct.name;
            obj.number       = infostruct.number;
            obj.units        = infostruct.units;
            obj.path2pro     = infostruct.path2pro;
            obj.protocolfile = infostruct.profilename;
            obj.path2stim    = infostruct.path2stim;
            obj.stimulusfile = infostruct.stimfilename;            
            if isempty(obj.stimulusfile), obj.epochinfosource = 'pClampProFile';
            else                          obj.epochinfosource = 'pClampAtfFile';
            end
            obj.analogwaveformtable = infostruct.analogwaveformtable;          
        end
        
        % ------------------------------------------- HELPER METHODS --------------------------------------------------------
        function identify(obj)
            % This function returns a small summary of the Analogout object. 
            if ~isscalar(obj), error('Analogout object must be scalar.'); end
            if isempty(inputname(1)), myname = 'This'; else myname = sprintf('"%s"',inputname(1)); end
            IDmessage = sprintf(['%s is an Analogout object:\n'...
                                 'Analogout name:   %s\n',...
                                 'Analogout number: %d\n',...
                                 'Units:            %s\n',...
                                 'EpochInfo:        %s\n',...
                                 'Protocol file:    %s.pro\n',...
                                 'Simulus file:     %s.atf\n'],...
                                 myname, obj.dacusername, obj.number, obj.units, obj.epochinfosource, obj.protocolfile, obj.stimulusfile);
             disp(IDmessage)
        end
        
        % --------------------------------------- GET/SET/SELECT METHODS ----------------------------------------------------
        
        function sf = getscalefactor(obj,INobj)
            % function extracts scalefactor by comparing the epochTable default values for the analog waveform to that 
            % recorded in the analog IN object INobj. This channel is assumed to be the analog IN recording the corresponding
            % secondary output of this channel. For different protocols, the method of scalefactor determination may differ.
            % Switch statement is used to determine how. 
            % --------------------
            
            if ~isscalar(obj), error('Analogout object must be scalar.'); end
            % 2. Calculate scalefactor
            switch obj.stimulusfile 
                case {'eCode_1_IV'} % all protocols with steps
                    epo1dur = fix((str2double(obj.analogwaveformtable.tab(1,:).timespan)/1000)*obj.analogwaveformtable.samplefreq);
                    epo2dur = fix(((str2double(obj.analogwaveformtable.tab(2,:).timespan)/1000)*obj.analogwaveformtable.samplefreq)+epo1dur);
                    epo1 = INobj.getsweep(1).Data(1:epo1dur);
                    epo2 = INobj.getsweep(1).Data(epo1dur:epo2dur);
                    sf = (nanmedian(epo2)-nanmedian(epo1))/-140;
                case {'eCode_2_Noise'}      % noise protocols
                    sf = 1;
                case {'eCode_2_Resonance'}  % resonance...
                    sf = 1; 
                case {'LSFINEST'}
                    epo1dur = fix((sum(obj.analogwaveformtable.tab(1:3,:).timespan)/1000)*obj.analogwaveformtable.samplefreq);
                    epo2dur = fix((sum(obj.analogwaveformtable.tab(1:4,:).timespan)/1000)*obj.analogwaveformtable.samplefreq);
                    epo1 = INobj.getsweep(1).Data(1:epo1dur);
                    epo2 = INobj.getsweep(1).Data(epo1dur:epo2dur);
                    if strcmp(INobj.units, 'nA')
                        sf = (nanmedian(epo2)-nanmedian(epo1))*10;
                    else
                        sf = (nanmedian(epo2)-nanmedian(epo1))/100;
                    end
                case {'LSCOARSE'}
                    epo1dur = fix((sum(obj.analogwaveformtable.tab(1:3,:).timespan)/1000)*obj.analogwaveformtable.samplefreq);
                    epo2dur = fix((sum(obj.analogwaveformtable.tab(1:4,:).timespan)/1000)*obj.analogwaveformtable.samplefreq);
                    epo1 = INobj.getsweep(1).Data(1:epo1dur);
                    epo2 = INobj.getsweep(1).Data(epo1dur:epo2dur);
                    if strcmp(INobj.units, 'nA')
                        sf = (nanmedian(epo2)-nanmedian(epo1))/-150*1000;
                    else
                        sf = (nanmedian(epo2)-nanmedian(epo1))/-150;
                    end
                otherwise
                    sf = 1;
            end
            % 3. Standardise scalefactor. Now scalefactor is determined, it has to be standardised, or rather, rounded to the 
            % nearest 2.5% increment. So for example, 1.269 scalefactor becomes 1.275)
            resolution = 2.5;
            sf = round(sf*100/resolution,0)*resolution/100; % round to nearest X%
        end
        
        function [holdI, holdV] = getholdingIorV(obj,CHobj)
            % function extracts holding current or voltage by looking at the intitial current injected in first step epoch 
            % recorded in the secondary analog IN object in CHobj. If epoch 'A' is not a step, holding current can not be determined 
            % (unless manually specified otherwise for certain protocols). 
            % --------------------      
            if ~isscalar(obj), error('Analogout object must be scalar.'); end
            
            holdI = NaN;
            holdV = NaN;
            % if current clamp; get holding current
            if strcmp(CHobj.getin('signal','primary').units, 'mV') && strcmp(CHobj.getin('signal','secondary').units, 'pA')                
                epoA = CHobj.getin('signal','secondary').getsweep(1).getepoch('Name', 'Epoch A');
                if strcmp(epoA.typestr, 'step') && epoA.amplitude == 0
                    holdI = nanmedian(epoA.Data);
                end
            % elseif voltage clamp; get holding voltage
            elseif strcmp(CHobj.getin('signal','primary').units, 'pA') && strcmp(CHobj.getin('signal','secondary').units, 'mV')                
                epoA = CHobj.getin('signal','secondary').getsweep(1).getepoch('Name', 'Epoch A');
                if strcmp(epoA.typestr, 'step') && epoA.amplitude == 0
                    holdV = nanmedian(epoA.Data);
                end            
            end        
        end
                
        

        % ------------------------------------------ PLOTTING METHODS -------------------------------------------------------
        function plot(obj,varargin)
            % Plot the analog output waveform for sweep(s) of this channel. Analogous to the "waveform preview" option in pClamp.
            % Note, if numel(obj)>1, then all waveforms are plotted in same graph panel.
            for i = 1:numel(obj)
                hold on
                    arrayfun(@(x) obj(i).waveform(x).plot(varargin{:}),1:obj.analogwaveformtable.nrofsweeps,'UniformOutput',false)
                hold off
            end
        end  
            
        function wf = waveform(obj,sweepidx)
            % generates the analog output waveform for sweep number 'sweepidx'.
            wf = cell(size(obj));
            for i = 1:numel(obj)
                wf{i} = obj(i).analogwaveformtable.waveform(sweepidx)*obj(i).scalefactor;
            end
            if numel(wf)==1, wf=wf{:}; end
        end
        
    end
    
end

