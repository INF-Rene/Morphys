classdef Stimset < Sharedmethods
    %Stimset object
    %   A stimset object is a group of acquisition sweeps in an NWBfile that were acquired by the same protocol/stimset
    %
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       René Wilbers (renewilbers@gmail.com)
    %   Created:      29-03-2019     
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
%###################################################### PROPERTIES ########################################################## 
    properties
        guid_NWB            % Globally unique identifier of parent NWB object.
        filename            % filename of parent NWB
        filedirectory       % directory of parent NWB
        name                % Stimset name
        datetimestart
        datetimeend
        
        sweeptable
        sweepnrs
        activeHS
        associated_stimsets  % simultaneously recorded stimsets of which sweeps will be saved here as well
        labbooknum
        labbooknum_keys
        labbook_timestamps
        datalocs
        stimdatalocs
        
        updatephase = 1     % phase of class; (phase 1 = no setup settings info, phase 2 = settings added, phase 3 = settings added and analysed)   
        
        nwbchannels              % list of nwbchannel objects
        wavename
        samplefreq
        
        nrofsweeps  = 0     % number of Analog Input channels recorded

    end
    
%####################################################### METHODS ############################################################
    methods
        % ----------------------------------------- CLASS CONSTRUCTOR -------------------------------------------------------
        function obj = Stimset(varargin)
            % call to superclass constructor. Takes a struct or name/value pair arguments as inputs.
            obj = obj@Sharedmethods;
                        
            % check/process inputs
            if     numel(varargin) == 0, return % returns empty Channel object   
            elseif numel(varargin) == 1
                s = varargin{:};
                if isstruct(s)
                    flds = fieldnames(s);
                    for i=1:numel(flds)
                        obj = obj.set(flds{i},s.(flds{i}));
                    end
                else
                    error('Please provide a struct or name/value pairs as inputs')
                end
            elseif mod(numel(varargin),2)~=0
                error('Uneven number of name/value pairs.')
            elseif ~all(cellfun(@ischar,varargin(1:2:end)))
                error('All property names must be strings')
            else
                for i=1:2:numel(varargin)
                    obj = obj.set(varargin{i},varargin{i+1});
                end
            end
%             if ~isempty(obj.number),
%                 obj.name = sprintf('Channel %d',obj.number);
%             end

        end
        % -------------------------------------------- OTHER METHODS --------------------------------------------------------
        function obj = addnwbchannels(obj)
            %create AD objects and add sweeps
            [ADnames, ~, ic] = unique(obj.sweeptable.ADname);
            for i=1:numel(ADnames)
                swpt = obj.sweeptable(ic==i,:);
                adnr=str2num(ADnames{i}(end))+1;
                
                adLB = squeeze(obj.labbooknum(adnr,:,:))';
                adLB=cell2table(num2cell(adLB));
                adLB.Properties.VariableNames=genvarname(obj.labbooknum_keys(:,1));
                adLB.Properties.VariableUnits=obj.labbooknum_keys(:,2);
                adLB.TimeStamp = datetime(adLB.TimeStamp,'ConvertFrom', 'epochtime', 'Epoch', '1904-01-01');
                
                ad2add = nwbchannel('guid_stimset', obj.guid, 'filename', obj.filename, 'filedirectory', obj.filedirectory,...
                    'name', ADnames{i}, 'number', adnr, 'sweeptable', swpt, 'sweepnrs', unique(swpt.sweepnr),...
                    'associated_stimsets', unique(swpt.protocol),'labbooknum', adLB, 'nrofsweeps', numel(unique(swpt.sweepnr)));
                
                ad2add = ad2add.addsweeps;
                ad2add = ad2add.addstimwaves; % add stim data (either single sweep or multiple depending on protocol)
                if isempty(obj.nwbchannels)
                    obj.nwbchannels        = ad2add;
                else
                    obj.nwbchannels(end+1) = ad2add;
                end
            end
            
        end        
        
        function ap = getnwbchannel(obj,varargin)
            % get SWEEP(s) from list. 
            % Returns SWEEP as a 1xN SWEEP object with N equal to numel(idx). If no input provided, returns all sweeps. 
            % Uses GETITEM from SHAREDMETHODS superclass. See function description of GETITEM for information on how to use 
            % variable arguments in for selecting items. 
            %
            % see also SWEEP, GETITEM, SHAREDMETHODS
            ap = getitem(obj,'nwbchannels',0,varargin{:});
        end
        
        function obj = analysestimset(obj)
            % perform default analysis of analog in. Only for primaries with voltage traces are analysed. 
            for i = 1:numel(obj)
                fprintf('\t\tanalysing stimset %s\n',obj(i).wavename)
                obj(i).nwbchannels = obj(i).getnwbchannel.analysenwbchannel;
                obj(i).updatephase = 3;
            end
        end
        

        % ----------------------------------------- PLOTTING METHODS --------------------------------------------------------
        
        function plot(obj,varargin)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
            for i=1:numel(obj)
                figure %create figure for each stimset
                set(gcf,'position',[2100,-100,1120,801],'color',[1 1 1])             
                
                nrofsubplots = numel(obj.getnwbchannel)*2;
                ax_handles   = zeros(nrofsubplots,1);
                subplotcount = 1;
                for j=1:numel(obj.getnwbchannel)
                    ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                    subplotcount  = subplotcount+1;

                    obj.getnwbchannel(j).getsweep.plot(varargin{:});
                    title([strrep(obj.filename,'_','\_') ' : ' obj.getnwbchannel(j).name])

                    ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                    subplotcount  = subplotcount+1;
                    
                    obj.getnwbchannel(j).getstimwave.plot(varargin{:});
                    title([strrep(obj.name,'_','\_') ' : DA' num2str(obj.getnwbchannel(j).number-1)])

                    
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
            end

        end       
        
        function plotanalysis(obj)
            % Plots all sweeps in the Stimset object
            % Add here option to filter out sweeps that did not pass QC
            % --------------------
            for i=1:numel(obj)
                figure %create figure for each stimset
                set(gcf,'position',[2100,-100,1120,801],'color',[1 1 1])
                title(strrep(obj.name,'_',' '))
                
                nrofsubplots = numel(obj.getnwbchannel)*2;
                ax_handles   = zeros(nrofsubplots,1);
                subplotcount = 1;
                for j=1:numel(obj.getnwbchannel)
                    ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                    subplotcount  = subplotcount+1;

                    obj.getnwbchannel(j).getsweep.plotanalysis;
                    title([strrep(obj.filename,'_','\_') ' : ' obj.getnwbchannel(j).name])

                    ax_handles(subplotcount) = subplot(nrofsubplots,1,subplotcount); 
                    subplotcount  = subplotcount+1;
                    
                    obj.getnwbchannel(j).getstimwave.plot;
                    title([strrep(obj.name,'_','\_') ' : DA' num2str(obj.getnwbchannel(j).number-1)])

                    
                end
                xlabel('milliseconds')
                linkaxes(ax_handles,'x')
            end

        end  
% 
%         function hfig = plotanalysis(obj)
%             % Plots the Channel object with Epoch boundaries and results of basic analysis. Returns figure handle.
%             % --------------------
%             if isscalar(obj)
%                 ax_handles = zeros(obj.nrofanalogins,1);
%                 for i=1:obj.nrofanalogins
%                     ax_handles(i) = subplot(obj.nrofanalogins,1,i); 
%                     analogin2plot = obj.getin(i);
%                     analogin2plot.plotanalysis
%                     if analogin2plot.updatephase == 2 || analogin2plot.updatephase == 3
%                         configstring = sprintf('- %s (Channel %d)',analogin2plot.signal,obj.number);
%                     else
%                         configstring = '';
%                     end
%                     title(sprintf('Analog Input %d %s',analogin2plot.number,configstring))
%                     ylabel(analogin2plot.units)
%                 end
%                 xlabel('milliseconds')
%                 linkaxes(ax_handles,'x')
%                 hfig = gcf;
%                 set(hfig, 'color', [1 1 1])
%             else
%                 arrayfun(@(x) obj(x).plotanalysis,1:numel(obj),'UniformOutput',false)
%             end
%         end
    end
end

