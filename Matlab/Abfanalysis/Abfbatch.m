classdef Abfbatch < Sharedmethods & Sharedpaths
    %ABFBATCH List or batch of Abffile objects. 
    %   A class to oversee the handling of abf files. Includes methods to:
    %   - Convert .abf files into matlab abffile objects
    %   - Save abffile objects in a designated folder
    %   - Return tables of metadata from all relevant levels, ready for uploading to database using guids as match fields.
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ##########################################################
    
    properties (Access = public)
        nrofabfs = 0    % number of abfs in list
        abfs            % list of Abffile objects 
        savename        % auto-generated batch name (batch-date-guid)
    end

%####################################################### METHODS ############################################################
    
    methods
        
        % ----------------------------------------- CLASS CONSTRUCTOR--------------------------------------------------------
        function obj = Abfbatch(varargin)
            % create an ABFBATCH object.
            % provide information on where to find abf files to include in batch, and the location of a SETUPSETTINGS object
            % that describes the set-up configuration, which is required to match inputs to outputs. 
            % 
            % Usage:
            % [...] =   Abfbatch() makes an empty ABFBATCH object 
            % [...] =   Abfbatch('gui') makes an ABFBATCH object using a gui to select abffiles and setup settings object. 
            % [...] =   Abfbatch(PATHS2MATS) makes an ABFBATCH object using PATHS2MATS. PATHS2MATS must be a scalar or
            %           vector cell array of strings, each pointing to an existing Abffile object saved as a ".mat" file.
            % [...] =   Abfbatch(PATHS2ABFS,SETUPSETTINGS) makes an ABFBATCH object using PATHS2ABFS. PATHS2MATS must be a 
            %           scalar or vector cell array of strings, each pointing to an existing abffile (".abf").
            %           SETUPSETTINGS must be a SETUPSETTINGS object, or a string of full path to a saved setup settings
            %           object.
            % See also ABFFILE, SETUPSETTINGS, SHAREDMETHODS            
            
            % call superclass constructors
            obj = obj@Sharedmethods;
            obj = obj@Sharedpaths;

            % check inputs... this has become a monster. Think I checked all conditions, common uses are ok. leave for now.
            if isempty(varargin), 
                 disp('No inputs provoded, returning empty');
                 return
            elseif nargin == 1    % can be empty, can be the string "gui", or cell array containing paths to mat files
                
                % empty input
                if isempty(varargin{1}),
                    disp('Empty input provoded, returning empty');
                    return
                
                % string input, must "gui")
                elseif ischar(varargin{1}) && strcmpi(varargin{1},'gui')
                    
                    % first choose whether abf or mat
                    filetype = questdlg('Please choose the file type you wish to work with:','Choose file type','.abf','.mat','.abf');
                    if isempty(filetype); 
                        disp('No file type chosen, returning empty'); 
                        return 
                    end
                    
                    [filenames, path2abf] = uigetfile(fullfile(obj.get('path2abfs'),filetype), sprintf('Pick the "%s" files you wish to include in batch',filetype),'MultiSelect', 'on');
                    if ischar(filenames), 
                        filenames = cellstr(filenames); 
                    end
                    
                    if ~iscell(filenames), % not a string, not a cell, then uigetfile did not return anything, list empty
                        disp('No files chosen, returning empty'); 
                        return 
                    end
                    if isempty(filenames), 
                        disp('No files chosen, returning empty batch')
                        return
                    end
                    
                    % make full paths to selected files
                    filepathlist = cellfun(@(x) fullfile(path2abf,x),filenames,'UniformOutput',false);                     
                    switch filetype
                        case '.abf'
                            % abfs need a setup settings object
                            [setupfile, setuppath] = uigetfile(obj.path2setupsettings, 'Pick the Setup Settings object for batch analysis','MultiSelect', 'off'); 
                            if ~ischar(setupfile), % not a string, then uigetfile did not return anything
                                disp('No setup settings chosen, returning empty batch')
                                return
                            else
                                % feed into itself to jump to nargin == 2 condition
                                obj = Abfbatch(filepathlist,fullfile(setuppath,setupfile));
                            end
                        case '.mat'
                            % mat files are assumed to already have setup settings icorporated.
                            % feed into itself to jump to nargin = 1 & iscell condition
                            obj = Abfbatch(filepathlist);
                    end
                    
                % cell input, in which case a cell array of existing ".mat" files should be provided
                elseif iscell(varargin{1})
                    if ~isscalar(varargin{1}) && ~isvector(varargin{1})
                        error('Input must be a scalar or vector cell array')
                    elseif ~all(cellfun(@(x) exist(x,'file'),varargin{1})==2),
                        error('Input contains file paths that do not exist')
                    elseif ~all(cellfun(@strncmp,cellfun(@fliplr,reshape(varargin{1},length(varargin{1}),1),'UniformOutput',false),...
                                                     repmat({'tam.'},1,numel(varargin{1}))',...
                                                     repmat({4},1,numel(varargin{1}))')); % test if they all end with '.tam'
                        error('Input must only contain paths to mat files.')
                    end
                    % all ok, so start adding mats using paths provided
                    for i = 1:numel(varargin{1})
                        fprintf('Adding mat file #%05d to batch:',i)
                        obj = obj.addabf(varargin{1}{i});
                    end                    
                else
                    error('Only the string "gui" or cell input allowed, provided: %s.',class(varargin{1}))
                end
                
            elseif nargin == 2
                % test first input
                if isempty(varargin{1}),
                    disp('Empty input provided, returning empty')
                    return
                elseif ~iscell(varargin{1}),
                    error('First input argument must be a cell array or strings describing full path to abf file(s)')
                elseif ~isscalar(varargin{1}) && ~isvector(varargin{1})
                    error('First input argument must be a scalar or vector cell array')
                elseif ~all(cellfun(@(x) exist(x,'file'),varargin{1})==2),
                    error('First input argument contains file paths that do not exist')
                elseif ~all(cellfun(@strncmp,cellfun(@fliplr,reshape(varargin{1},length(varargin{1}),1),'UniformOutput',false),...
                                                 repmat({'fba.'},1,numel(varargin{1}))',...
                                                 repmat({4},1,numel(varargin{1}))')); % test if they all end with '.abf'
                    error('First input argument must only contain paths to abf files.')
                end
                
                % test second input
                if isa(varargin{2},'Setupsettings')
                    setup = varargin{2};
                else
                    if ~ischar(varargin{2})
                        error('Second input argument should be a character array describing full path to a saved setup settings object or a setup settings object.')
                    elseif exist(varargin{2},'file')~=2
                        error('Path to setup settings object provided as second input argument does not exist.')
                    else
                        try
                            setup = load(varargin{2});
                            flds  = fieldnames(setup);
                            setup = setup.(flds{1}); % When loading setup settings, the Setupsettings object is named 'obj' and becomes a field in struct setup. Hence setup.obj.
                        catch err
                            assignin('base','setupsettingsloadingerror',err)
                            error('Setup settings object provided in second input argument doesn''t load properly, please check correctness.')
                        end
                        if ~isa(setup,'Setupsettings')
                            error('Loaded file is not a setup settings object.')
                        end
                    end
                end
                % if tests passed, we now have a non-empty cell array vector of existing paths to .abf files, and we have a valid setup settings object
                
                % add abfs using path2abf and setupsettings object
                for i = 1:numel(varargin{1})
                    fprintf('\nAdding abf file #%05d to batch:',i)
                    obj = obj.addabf(varargin{1}{i},setup);
                end
            end
            
            % give it a name
            obj.savename = sprintf('Batch_%s_%s.mat',sprintf('%d-%0.2d-%0.2d',year(obj.created),month(obj.created),day(obj.created)),obj.guid);
        end
        
        % ------------------------------------------ HELPER METHODS ---------------------------------------------------------
        function obj = addabf(obj,varargin)
            % Add ABFFILE objects to list of abfs.
            % Use:
            % ... = ADDABF(OBJ, ABFFILE); where ABFFILE is an ABFFILE object or a path to a saved ABFFILE object.
            % ... = ADDABF(OBJ, PATH2FILE, SETUPSETTINGS); where PATH2FILE is a full path to an ABF file in Axon Binary Format, or
            %       or MAT file of a saved ABFFILE object. SETUPSETTINGS is a SETUPSETTINGS object describing setup settings
            %       used during recording.
            
            if ~isscalar(obj), error('Object must be scalar.'); end 
            if isempty(varargin), disp('no abf file specified, no action taken.'), return
            elseif numel(varargin) == 1, % add abffile object directly or from path
                if isa(varargin{1},'Abffile'),
                    abf2add = varargin{1};
                else
                    path2file  = varargin{1};
                    [~,fn,ext] = fileparts(path2file);
                    switch ext
                        case '.mat', 
                            if exist(path2file, 'file') == 2
                                fprintf(' >> %s%s...',fn,ext)
                                abf2add = load(path2file);
                                flds    = fieldnames(abf2add);
                                abf2add = abf2add.(flds{1}); 
                                if ~isa(abf2add,'Abffile'), 
                                    warning('\tEncountered ".mat" file that did not contain an Abffile object: %s.mat\n\t\t\tFile not appended',fn); 
                                    return % so nothing is appended to list in this case
                                end
                                fprintf('done\n')
                            else error('Following file does not exist: %s',path2file)
                            end
                        otherwise, error('Please provide either (1) an Abffile object, (2) a path to an Abffile object or (3) a path to an Abffile (*.abf) and a setup settings object.')
                    end
                end
            elseif numel(varargin) == 2,
                path2file     = varargin{1};
                setupsettings = varargin{2};
                
                [~,fn,ext] = fileparts(path2file);
                switch ext
                    case '.mat', 
                        if exist(path2file, 'file') == 2
                            fprintf(' >> %s%s...',fn,ext)
                            abf2add = load(path2file);
                            flds    = fieldnames(abf2add);
                            abf2add = abf2add.(flds{1}); 
                            if ~isa(abf2add,'Abffile'), 
                                warning('\tEncountered ".mat" file that did not contain an Abffile object: %s.mat\n\t\t\tFile not appended',fn); 
                                return % so nothing is appended to list in this case
                            end
                            fprintf('done\n')
                        else error('Following file does not exist: %s',path2file)
                        end
                    case '.abf', abf2add = Abffile(path2file,setupsettings);
                    otherwise, error('Unsupported extension: %s',ext)
                end
            else
                error('too many inputs provided')
            end
            
            % add guid to Abffile object
            abf2add = abf2add.set('guid_batch',obj.guid);
            if isempty(obj.abfs), 
                 obj.abfs        = abf2add; 
            else obj.abfs(end+1) = abf2add; 
            end
            obj = obj.updateabfstats;
        end
        
        function e = isempty(obj)
            % overload of the isempty function
            if isempty(obj.abfs), e=true; else e=false; end
        end

        function obj = updateabf(obj,idx,varargin)
            % update properties of an abffile in list
            %
            % See also ABFFILE, UPDATEITEM. 
            if ~isscalar(obj), error('Object must be scalar'); end
            obj = obj.updateitem('abfs',idx,varargin{:}).updateabfstats; % update stats
        end
        
        function obj = sortabfs(obj,varargin)
            % sort ABFFILEs in list
            % Returns object with list sorted according to value of a property. Default is to sort on Abffile property
            % filename.
            % Usage:
            %   obj = sortabfs(OBJ,PROPERTY,DESCENDING)
            %
            % Example:
            %   obj.sortabfs()                  --> sort abfs in ascending order by file name
            %   obj.sortabfs('fileSystemDate')  --> sort abfs in ascending order by time of saving file to disk
            %   obj.sortabfs('proName',1)       --> sort abfs in descending order by the protocol name
            %
            % See also ABFFILE, SORTITEM, SORT
            for i = 1:numel(obj); obj(i) = obj(i).sortitems('abfs',varargin{:}); end; 
        end
        
        function abf = getabf(obj,varargin)
            % get an ABFFILE from list. Returns ABFFILE(s) object specified by varargin. If no input provided, returns all. 
            % Usage example:
            % ABFFILE(s) can be retrieved using indexing with real positive integers, logical arrays, or by using 
            % name/value pairs:
            %
            % ABFFILE = ABFBATCH.GETABF(1)                          % returns first ABFFILE in list
            % ABFFILE = ABFBATCH.GETABF(1:3:end)                    % returns every 3rd ABFFILE in list
            % ABFFILE = ABFBATCH.GETABF([true true false])          % returns 1st two ABFFILES in list
            % ABFFILE = ABFBATCH.GETABF('filename','cell01_0000')   % returns all ABFFILE objects with matching filename
            %
            % Uses GETITEM from SHAREDMETHODS superclass to get items; see function description of GETITEM for more detailed
            % information on how to use varargin for selecting items. 
            %
            % see also ABFFILE, GETITEM, SHAREDMETHODS
            abf = obj.getitem('abfs',0,varargin{:});
        end
        
        function obj = removeabf(obj,varargin)
            % remove ABFFILE(s) from list. 
            %
            % see also ABFFILE, GETITEM, REMOVEITEM, SELECTAP, SHAREDMETHODS
            for i=1:numel(obj), obj = obj.removeitem('abfs',varargin{:}).updateabfstats; end
        end
        
        function obj = inspectbatch(obj)
            % browse through plots of ABFFILES in using the BATCHVIEWER gui
            obj = batchviewer(obj);
        end
        
        function obj = plus(obj1, obj2)
            % Add two ABFBATCH objects. Note the creation date and guid are kept from 1st ABFBATCH object 
            % provided. Future: add duplicate abffiles warning?
            obj1.abfs = [obj1.abfs,obj2.abfs];
            obj = obj1.updateabfstats;
        end 
      
        % ----------------------------------------- SELECTING METHODS -------------------------------------------------------
        function obj = selectabf(obj,varargin)
            % select ABFFILE(s) from list. 
            %
            % see also ABFFILE, GETITEM, SELECTITEM, REMOVEAP, SHAREDMETHODS
            for i=1:numel(obj), obj = obj.selectitem('abfs',varargin{:}).updateabfstats; end
        end    
        
        function obj = selectchannel(obj,varargin)
            % make ABFBATCH object including only ANALOGINS/OUTS associated with CHANNEL (specified by channel number).
            % See also ABFFILE, ANALOGIN, CHANNEL, SELECTITEM.
            if ~isscalar(obj), error('Object must be scalar.'); end
            if nargin == 1, return
            else
                for i=1:obj.nrofabfs
                    obj.abfs(i)=obj.getabf(i).selectchannel(varargin{:});
                end
            end
        end
        
        function obj = selectin(obj,varargin)
            % make ABFBATCH object including only specified ANALOGINs.
            % See also ABFFILE, ANALOGIN, SELECTITEM.
            if ~isscalar(obj), error('Object must be scalar.'); end
            if nargin == 1, return
            else
                for i=1:obj.nrofabfs
                    for ii=1:obj.getabf.nrofchannels
                        obj.abfs(i).channels(ii)=obj.getabf(i).getchannel(ii).selectin(varargin{:});
                    end
                end
            end
        end
        
        function obj = selectout(obj,varargin)
            % make ABFBATCH object including only specified ANALOGOUTs.
            % See also ABFFILE, ANALOGOUT, SELECTITEM.
            if ~isscalar(obj), error('Object must be scalar.'); end 
            if nargin == 1, return
            else
                for i=1:obj.nrofabfs
                    for ii=1:obj.getabf.nrofchannels
                        obj.abfs(i).channels(ii)=obj.getabf(i).getchannel(ii).selectout(varargin{:});
                    end
                end
            end
        end
        
        % ------------------------------------------ UPDATING METHODS -------------------------------------------------------
        function obj = updatenrofabfs(obj)
            if ~isscalar(obj), error('Object must be scalar.'); end
            obj.nrofabfs = numel(obj.abfs);
        end
        
        function obj = updateabfstats(obj)
            % update all stats relating to the list of abfs
            % currently only nrofabfs, could be extended in future...
            obj = obj.updatenrofabfs;
        end
        
        function obj = analysebatch(obj)
            % perform default analysis on abfs in batch
            disp('Analysing batch:')
            obj.abfs = obj.getabf.analyseabf;
        end
        
        % ----------------------------------------- EXTRACTING METHODS ------------------------------------------------------

        function t = abftable(obj,metadataproperties)
            % make a table listing all the properties of ABFFILES in list. METADATAPROPERTIES is a cell array of
            % strings corresponding to property names of ABFFILE objects. If no METADATAPROPERTIES provided, chooses default.
            % See also ABFFILE
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.metadata);
            else          t = struct2table(obj.getabf.metadata(metadataproperties));
            end
        end
        
        function t = channeltable(obj,metadataproperties)
            % make a table listing all the properties of CHANNELS of ABFFILES in list. METADATAPROPERTIES is a cell array of
            % strings corresponding to property names of ABFFILE objects. If no METADATAPROPERTIES provided, chooses default.
            % See also CHANNEL
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.getchannel.metadata);
            else          t = struct2table(obj.getabf.getchannel.metadata(metadataproperties));
            end
        end        
        
        function t = analogintable(obj,metadataproperties)
            % make a table listing all the properties of ANALOGINS of ABFFILES in list. METADATAPROPERTIES is a cell 
            % array of strings corresponding to property names of ANALOGIN objects. If no METADATAPROPERTIES provided, 
            % chooses default.
            % See also ANALOGINPUT
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.getchannel.getin.metadata);
            else          t = struct2table(obj.getabf.getchannel.getin.metadata(metadataproperties));
            end
        end        
        
        function t = analogouttable(obj,metadataproperties)
            % make a table listing all the properties of ANALOGOUTS of ABFFILES in list. METADATAPROPERTIES is a cell 
            % array of strings corresponding to property names of ANALOGOUT objects. If no METADATAPROPERTIES provided, 
            % chooses default.
            % See also ANALOGOUTPUT
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.getchannel.getout.metadata);
            else          t = struct2table(obj.getabf.getchannel.getout.metadata(metadataproperties));
            end
        end

        function t = sweeptable(obj,metadataproperties)
            % make a table listing all the properties of every SWEEP of ABFFILES in list. METADATAPROPERTIES is a cell 
            % array of strings corresponding to property names of SWEEP objects. If no METADATAPROPERTIES provided, chooses 
            % default.
            % See also SWEEP
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.getchannel.getin.getsweep.metadata);
            else          t = struct2table(obj.getabf.getchannel.getin.getsweep.metadata(metadataproperties));
            end
        end
        
        function t = epochtable(obj,metadataproperties)
            % make a table listing all the properties of every EPOCH of ABFFILES in list. METADATAPROPERTIES is a cell 
            % array of strings corresponding to property names of EPOCH objects. If no METADATAPROPERTIES provided, chooses 
            % default.
            % See also EPOCH
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, t = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.metadata);
            else          t = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.metadata(metadataproperties));
            end
        end
        
        function t = aptable(obj,level,metadataproperties)
            % make a table listing all the properties of every ACTIONPOTENTIAL of ABFFILES in list. METADATAPROPERTIES 
            % is a cell array of strings corresponding to property names of ACTIONPOTENTIAL objects. If no METADATAPROPERTIES 
            % provided, chooses default. LEVEL is the level in hierarchy from which ACTIONPOTENTIALS should be retrieved,
            % that of SWEEPS or EPOCHS. 
            % See also ACTIONPOTENTIAL
            if ~isscalar(obj), error('Object must be scalar'); end
            if nargin==1, level='epoch'; end
            if nargin<3,
                switch level
                    case 'epoch', t = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.getap.metadata);
                    case 'sweep', t = struct2table(obj.getabf.getchannel.getin.getsweep.getap.metadata);
                    otherwise, error('unknown level specified. choose epoch or sweep.')
                end
            else
                switch level
                    case 'epoch', t = struct2table(obj.getabf.getchannel.getin.getsweep.getepoch.getap.metadata(metadataproperties));
                    case 'sweep', t = struct2table(obj.getabf.getchannel.getin.getsweep.getap.metadata(metadataproperties));
                    otherwise, error('unknown level specified. choose epoch or sweep.')
                end
            end
        end   
        
        % ------------------------------------------ PLOTTING METHODS -------------------------------------------------------
        function plot(obj,varargin)
            % function plots all Abffile objects. 
            if ~isscalar(obj), error('Object must be scalar.'); end
            for i=1:obj.nrofabfs
                figure()
                obj.getabf(i).plot(varargin{:});
            end
        end    

        % ------------------------------------------ EXPORTING METHODS ------------------------------------------------------
        function exporttables4db(obj,destinationdir,extension)
            % function to create tables ready for database import. 
            % save with batch guid in different folders: abfs, analogouts, analogins, swps, ... 
            
            if ~isscalar(obj), error('Abfbatch object must be scalar.'); end
            % check inputs
            supportedext = {'xls', 'xlsm', 'xlsx', 'csv', 'txt'};
            if isempty(extension) || ~ismember(extension,{'xls', 'xlsm', 'xlsx', 'csv', 'txt'})
                error('Only ".%s" and ".%s" extensions are supported.',strjoin(supportedext(1:end-1),'", ".'),supportedext{end})
            end
            % save 'em
            tabletypes = {'abfs','channels','analogouts','analogins','sweeps','epochs','aps'};
            cellfun(@(x) obj.savetableas(x,destinationdir,extension),tabletypes,'UniformOutput',false)
        end
        
        function savetableas(obj,whichTable,destinationdir,extension)
            % Function to print tables to file. 
            % destinationdir specifies the directory to store
            % Supported filetypes/extensions are: 'xls', 'xlsm', 'xlsx', 'csv' and 'txt'.
            % Uses Batch GUID as file name: Abfbatch_date_GUID
            % See also ABFTABLE, INTABLE, OUTTABLE, SWEEPTABLE, EPOCHTABLE, APTABLE, WRITETABLE.
            
            if ~isscalar(obj), error('Object must be scalar.'); end
            %check inputs
            supportedext = {'xls', 'xlsm', 'xlsx', 'csv', 'txt'};
            if isempty(extension) || ~ismember(extension,{'xls', 'xlsm', 'xlsx', 'csv', 'txt'})
                error('Only ".%s" and ".%s" extensions are supported.',strjoin(supportedext(1:end-1),'", ".'),supportedext{end})
            end
            
            switch whichTable
                case 'abfs',        table2print = abftable(obj);
                case 'channels',    table2print = channeltable(obj);
                case 'analogins',   table2print = analogintable(obj);  
                case 'analogouts',  table2print = analogouttable(obj);
                case 'sweeps',      table2print = sweeptable(obj);
                case 'epochs',      table2print = epochtable(obj);
                case 'aps',         table2print = aptable(obj,'epoch');
                otherwise,          error('Requested table not available, choose one of following: abfs,channels,analogouts,analogins,sweeps,epochs,aps')
            end
            
            % fix datetime and duration formats by converting to string
            for i = table2print.Properties.VariableNames, 
                switch i{:}
                    case {'fileduration','filetimestart','filetimeend','prostopwatchtime','timespan','datetimestart','datetimeend','sweeplagtime'}
                        table2print.(i{:}) = cellstr(char(table2print.(i{:})));
                end
            end
            
            % make a name and save it
            datestring = datetime(year(now),month(now),day(now),hour(now),minute(now),second(now),'Format','yyyy-MM-dd_HHmmss');
            path2file  = fullfile(destinationdir,sprintf('Batch_%s_%s_%s.%s',datestring,whichTable,obj.guid,extension));
            writetable(table2print,path2file,'WriteVariableNames',false)
        end
    end       
end

