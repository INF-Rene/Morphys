classdef Sharedpaths
    %SHAREDPATHS Class to hold all relevant paths
    %   Class to hold all relevant paths. Not sure whether this is common practice, but seems useful up to now. Used so that 
    %   ABFBATCH and ABFFILE objects, which inherit from SHAREDPATHS know where to retrieve / store data.
    %   
    %   Note: Use of a global variable to store username. 
    %   Upon creation of a SHAREDPATH object or any class that inherits from SHAREDPATH (now only ABFBATCH & ABFFILE objects)
    %   a check is made of whether a current user name is defined in the global variable USERNAME. This is to make sure that
    %   all paths point to user-specific storage folders (see use of USERID in SHAREDPATH constructor function).
    %   
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    properties (Hidden = false)
        userid              = 'EJM'; % pulled from global variable 'USERNAME'. 
    end
    
    properties (Hidden = true)
        listofuserids       = {'AKS','DBH','DRU','EJM', 'GTS','IKS','JDZ','JOR','JSR','MBV','NAG','RBP','RWS','SHT','THK','TKN'};
        dir_base            = '/Users/elinemertens/Downloads/Morphys-master';
        dir_ephys           = fullfile('Morphys','Data','Electrophysiology');
        dir_abfs            = fullfile('Morphys','Data','Electrophysiology','Abffiles');
        dir_setupsettings   = fullfile('Morphys','Data','Electrophysiology','SetupSettings');
        dir_analogwaveform  = fullfile('Morphys','Data','Electrophysiology','Protocols','AnalogWaveforms');
        dir_atfs            = fullfile('Morphys','Data','Electrophysiology','Protocols','StimulusFiles');
        dir_noise           = fullfile('Morphys','Data','Electrophysiology','Protocols','StimulusFiles','INF','eCodeNoise');
        path2abfs           = ''; % path to place where abf files are stored
        path2setupsettings  = ''; % path to setup settings object
        path2analogwaveforms= ''; % path to store/retrieve AnalogWaveform objects
        path2atfs           = ''; % path to store/retrieve axon text files
        path2noise          = ''; % place where noise files are stored (Originals derive from Rodrigo Perin at EPFL, these have been interpolated to accommodate an 8kHz sampling frequency. 8kHz chosen because it is compatible with the 4kHz lowpass Bessel filter used for rest of eCode_2.)
        defaultsetupsettings= 'setupSettings_INF.mat'; % default INF setup settings
    end
    
%####################################################### METHODS ############################################################
    methods (Access = public)
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Sharedpaths(username)
                        
            % set user id
            obj.setglobalusername('EJM')
%             Commented to allow parfor loop to execute
%             if nargin == 1 && ~isempty(username)
%                 obj.userid = username;
%             else
%                 obj.userid = obj.getglobalusername;
%                 if isempty(obj.userid)
%                     idx = listdlg('PromptString','Select a userid:',...
%                                 'SelectionMode','single',...
%                                 'ListString',obj.listofuserids);
%                     if isempty(idx)
%                         disp('No user name selected, returning empty Sharedpaths object')
%                         return
%                     else
%                         obj.setglobalusername(obj.listofuserids{idx})
%                         obj.userid = obj.getglobalusername;
%                     end
%                 end
%             end
            
            % complete stored paths
            if ~isempty(obj.dir_base)
                obj = obj.updatesharedpaths;
            end
        end
        
        function obj = updatesharedpaths(obj)
            % update paths depending on base directory. When base directory is changed, call this function
            obj.path2abfs=fullfile(obj.dir_base,obj.dir_abfs,obj.userid);
            obj.path2setupsettings=fullfile(obj.dir_base,obj.dir_setupsettings);
            obj.path2analogwaveforms=fullfile(obj.dir_base,obj.dir_analogwaveform,obj.userid);
            obj.path2atfs=fullfile(obj.dir_base,obj.dir_atfs);
            obj.path2noise=fullfile(obj.dir_base,obj.dir_noise);
        end
    end
    
%     methods (Access = {?Sharedpaths}) %(Access = private)
%         function prop = get(obj,propertyname)
%             % Get object properties. 
%             % V = get(obj,'PropertyName') returns the value of the specified property for the object 'obj' (scalar or non-scalar).
%             if nargin < 2, error('Too little input')
%             elseif ~ischar(propertyname) %|| ~ismember(propertyname,properties(obj))
%                 error('Please provide a valid property name (as string).')
%             end
%             if isscalar(obj)
%                 prop = obj.(propertyname);
%             else
%                 prop = arrayfun(@(x) obj(x).get(propertyname),1:numel(obj),'UniformOutput',false);
%                 prop = reshape(prop,size(obj));
%             end
%         end
%         
%         function obj = set(obj,varargin)
%             % Function assigns 'value' to any property matching string 'PropertyName' from object.
%             % Supports Property/value pair input arguments.
%             if nargin < 3 || mod(nargin-1,2)~=0, error('Property/value pairs must come in even number.'); end
%             prop2val = cat(1,varargin(1:2:end),varargin(2:2:end))';
%             for i=1:size(prop2val,1)
%                 propertyname = prop2val{i,1};
%                 value        = prop2val{i,2};
%                 for ii = 1:numel(obj), obj(ii).(propertyname) = value; end
%             end
%             
%             % if base dir was changed, update all dependent paths.
%             if ismember('dir_base',varargin(1:2:end))
%                 obj = obj.updatesharedpaths;
%             end
%         end
%         
%         function obj = updatesharedpaths(obj)
%             % update paths depending on base directory. When base directory is changed, call this function
%             obj = obj.set(  'path2abfs',            fullfile(obj.dir_base,obj.dir_abfs,obj.userid), ...
%                             'path2setupsettings',   fullfile(obj.dir_base,obj.dir_setupsettings), ...
%                             'path2analogwaveforms', fullfile(obj.dir_base,obj.dir_analogwaveform,obj.userid),...
%                             'path2atfs',            fullfile(obj.dir_base,obj.dir_atfs), ...
%                             'path2noise',           fullfile(obj.dir_base,obj.dir_noise));
%         end
%         
%     end
    
    methods (Static)
        
        function setglobalusername(id)
            % Set the current user id to be a global variable
            global USERNAME
            USERNAME = id;
        end
        
        function id = getglobalusername
            % Get current user id (from global variable USERNAME). If not defined id = empty.
            global USERNAME
            id = USERNAME;
        end
    end
    
end
