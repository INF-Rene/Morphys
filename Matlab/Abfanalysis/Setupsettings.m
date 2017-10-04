classdef Setupsettings < Sharedmethods
    %SETUPSETTINGS Setup configuration information
    %   This object contains all the information about the set-up configuration that is required to know what Analogout 
    %   (DAC) channels are associated with what Analogin (ADC) channels. This object is used when creating Abffile objects, 
    %   to allow the proper matching of inputs and outputs, such that for every analog input in an Abffile, the appropriate 
    %   Epoch information is available.
    %
    %   Example, create a setupsettings object with 2 amplifier channels:
    %   ss = Setupsettings;
    %   ss = ss.addchannel('number',1,'dacnum',0,'primary',0,'secondary',4)
    %   ss = ss.addchannel('number',2,'dacnum',1,'primary',1,'secondary',5)
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    %   See also SHAREDMETHODS
    
%###################################################### PROPERTIES ##########################################################    
    properties
        setupsettingsname = ''; % name of setup settings object. To be used as file name when saving object.
        nrofchannels = 0        % number of channels in Channel list
        channelnrs              % list of numbers of the channels. Usually 1-4, in case of 2 amplifiers
        analoginputnrs          % list of numbers (in range of 0-15, in case of digidata 1440A/B)
        analogoutputnrs         % list of numbers (in range of 0-3,  in case of digidata 1440A/B)
        channels                % list of channels defined for this setup
    end
    
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Setupsettings(setuptemplate)
            % create a Setupsettings object.
            % Start either with empty SETUPSETTINGS object and add channels manually, or with a Setupsettings template file.
            % 
            % FUTURE: Add option to construct object from xls/xml file?
            
            
            % call to superclass constructor
            obj = obj@Sharedmethods; 
            if nargin == 0, return
            elseif nargin == 1,
                if isa(setuptemplate,'Setupsettings')
                    % then a Setupsettings template has been provided, copy name and channels. For channels, since this is a
                    % template, a new GUID must be assigned to the channel. 
                    obj.setupsettingsname = setuptemplate.setupsettingsname;
                    for i=1:setuptemplate.nrofchannels
                        chnnl2add = setuptemplate.getchannel(i).set('guid',upper(char(java.util.UUID.randomUUID())));
                        obj = obj.addchannel(chnnl2add).updatechannelstats;
                    end
                else
                    error('Please provide a setup settings object.')
                end
            else
                error('too many inputs')
            end
            
        end
        
        % ------------------------------------------- GENERAL METHODS -------------------------------------------------------       
       
        function obj = addchannel(obj, varargin)
            % add a Channel to the list using either name/value pair arguments or an existing Channel object.
            % Example:
            % ss = ss.addchannel('number',1,'dacnum',0,'primary',0,'secondary',4)
            if ~isscalar(obj), error('object must be scalar'); end
            
            % in case a fully assembled Channel is provided
            if numel(varargin)==1 && isa(varargin{1},'Channel')
                chnnl2add = varargin{1};
            else
                chnnl2add = Channel(varargin{:});% make the channel
            end
            obj.checkchannel(chnnl2add)          % check it

            % add it
            if isempty(obj.channels), obj.channels        = chnnl2add;     
            else                      obj.channels(end+1) = chnnl2add;    
            end
            obj = sortchannels(obj);        % sort 'em
            obj = obj.updatechannelstats;   % update stats
        end
        
        function obj = updatechannel(obj,idx,varargin)
            % update properties of a channel in list
            %
            % See also CHANNEL, UPDATEITEM 
            if ~isscalar(obj), error('object must be scalar'); end
            chnnl = obj.getchannel(idx);    % retrieve the channel
            obj   = obj.removechannel(idx); % remove from list
            chnnl = chnnl.set(varargin{:}); % update it
            obj.checkchannel(chnnl)         % check it
            obj.channels(end+1) = chnnl;    % add it 
            obj = sortchannels(obj);        % sort 'em
            obj = obj.updatechannelstats;   % update stats
        end

        function obj = sortchannels(obj,varargin)
            % sort CHANNELs in list
            % Returns object with list sorted according to CHANNEL PROPERTY. 
            % Default PROPERTY is NUMBER, default order is ascending, unless DESCENDING is set to 1.
            % Usage:
            %   obj = sortaps(OBJ,PROPERTY,DESCENDING)
            %
            % See also ACTIONPOTENTIAL, SORTITEM, SORT
            
            if ~isscalar(obj), error('object must be scalar'); end
            if obj.nrofchannels > 0, 
                 obj = obj.sortitems('channels',varargin{:});
            end
        end
        
        function chnnl = getchannel(obj,varargin)
            % get CHANNEL(s) from list. 
            % Returns CHANNEL as a 1xN CHANNEL object with N equal to numel(idx). If no input provided, returns all CHANNELS. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % varargin for selecting items. 
            %
            % see also CHANNEL, GETITEM, SHAREDMETHODS
            chnnl = getitem(obj,'channels',0,varargin{:});
        end
        
        function obj = removechannel(obj,varargin)
            % remove CHANNEL(s) from list. Logical or numeric indexing can be done, but end keyword cannot be used. 
            %
            % see also CHANNEL, GETITEM, REMOVEITEM, SELECTAP, SHAREDMETHODS
            for i = 1:numel(obj); obj(i) = obj(i).removeitem('channels',varargin{:}).updatechannelstats; end
        end
        
        function obj = selectchannel(obj,varargin)
            % select CHANNEL(s) from list. Logical or numeric indexing can be done, but end keyword cannot be used.
            %
            % see also CHANNEL, GETITEM, SELECTITEM, REMOVEAP, SHAREDMETHODS
            for i = 1:numel(obj); obj(i) = obj(i).selectitem('channels',varargin{:}).updatechannelstats; end
        end 
        
        function obj = updatechannelstats(obj)
            % update CHANNEL numbers and counts
            % see also CHANNEL, GETCHANNELNRS, GETANALOGINPUTNRS, GETANALOGOUTPUTNRS.
            for i = 1:numel(obj)
                obj(i).nrofchannels = numel(obj(i).channels);
                if obj(i).nrofchannels == 0,
                    obj(i).channelnrs      = [];
                    obj(i).analoginputnrs  = [];
                    obj(i).analogoutputnrs = [];
                else
                    obj(i).channelnrs      = obj(i).getchannelnrs;
                    obj(i).analoginputnrs  = obj(i).getanaloginputnrs;
                    obj(i).analogoutputnrs = obj(i).getanalogoutputnrs;
                end
            end
        end
        
        function n = getchannelnrs(obj)
            % get all channel numbers of CHANNEL(s) in list (sorted).
            if ~isscalar(obj), error('object must be scalar'); end
            n = obj.getchannel.get('number');
            if iscell(n), n = cell2mat(n); end
            n = sort(n);
        end
        
        function n = getanalogoutputnrs(obj)
            % get all POSSIBLE/DEFAULT channel analog output numbers of CHANNEL(s) in list (sorted).
            %
            % See also CHANNEL, ANALOGOUT.
            if ~isscalar(obj), error('object must be scalar'); end
            n = obj.getchannel.get('dacnum');
            if iscell(n), n = cell2mat(n); end
            n = sort(n);
        end
        
        function n = getanalogoutputnrs_current(obj)
            analogouts = obj.getchannel.getout;
            if ~isempty(analogouts), 
                n = [analogouts.number];
            else
                n = [];
            end
        end
        
        function n = getanaloginputnrs(obj)
            % get all POSSIBLE/DEFAULT analog input numbers of CHANNEL(s) in list
            %
            % See also CHANNEL, ANALOGIN.
            if ~isscalar(obj), error('object must be scalar'); end
            n = obj.getchannel.getanaloginputnrs;
            n = sort(n);
        end
        
        function n = getanaloginputnrs_current(obj)
            n = [obj.getchannel.getin.number];
        end
        
        function checkchannel(obj,chnnl)
            % check if inputs/outputs of CHANNEL are not already in use or outside expected ranges
            %
            % See also CHANNEL, ANALOGIN, ANALOGOUT.
            if ~isscalar(obj),                                  error('object must be scalar'); end
            if     isempty (chnnl.number),                      error('All amplifier Channels require a number.');
            elseif ismember(chnnl.number,obj.channelnrs),       error('Channel number already in use.');
            elseif ismember(chnnl.dacnum,obj.analogoutputnrs),  error('DAC number already in use.');
            elseif ismember([chnnl.primary,chnnl.secondary,chnnl.scope,chnnl.mode],obj.analoginputnrs), 
                                                                error('Analog input number already in use.')
            elseif isempty([chnnl.number,chnnl.dacnum,chnnl.primary,chnnl.secondary,chnnl.scope,chnnl.mode]); 
                                                                error('Please do not add empty channels')
            end
            
            % check ranges
            if chnnl.dacnum < 0 || chnnl.dacnum > 3, error('dacnumbers expected to be in range of 0 to 3. ')
            elseif any([chnnl.primary,chnnl.secondary,chnnl.scope,chnnl.mode]<0) || any([chnnl.primary,chnnl.secondary,chnnl.scope,chnnl.mode]>15)
                error('AnalogInput numbers expected to be in range of 0 to 15')
            end
        end
        
        function [s, n] = getinsource(obj,idx)
            % get the source of input to an analog input channel. 
            % 1st output is a string of property name, 2nd is number of channel.
            s = ''; n = [];
            if ~isscalar(idx), error('idx must be scalar')
            elseif idx<0 || idx>15, error('idx must be in range of 0 to 15');
            elseif ~ismember(idx, obj.analoginputnrs), return; 
            end
            flds = {'primary','secondary','scope','mode'};
            for i = 1:obj.nrofchannels
                if ismember(idx,obj.getchannel(i).getanaloginputnrs),
                    for ii = 1:numel(flds)
                        if obj.getchannel(i).(flds{ii}) == idx, 
                            s = flds{ii}; 
                            n = obj.getchannel(i).number;
                        end
                    end
                end
            end
                
        end
    end
    
end

