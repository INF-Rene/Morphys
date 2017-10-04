function varargout = batchviewer(varargin)
    % BATCHVIEWER MATLAB code for batchviewer.fig
    %      BATCHVIEWER, by itself, creates a new BATCHVIEWER or raises the existing
    %      singleton*.
    %
    %      H = BATCHVIEWER returns the handle to a new BATCHVIEWER or the handle to
    %      the existing singleton*.
    %
    %      BATCHVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in BATCHVIEWER.M with the given input arguments.
    %
    %      BATCHVIEWER('Property','Value',...) creates a new BATCHVIEWER or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before batchviewer_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to batchviewer_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help batchviewer
    %   ---------------------------------------------------------------------------------------------------------------------
    %   Author:       Thijs Verhoog (thijs.verhoog@gmail.com)
    %   Created:      2017-08-18       
    %   Modifications (Date/Name/Description):
    %   ---------------------------------------------------------------------------------------------------------------------
    
    % Last Modified by GUIDE v2.5 13-Sep-2017 20:28:37

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @batchviewer_OpeningFcn, ...
                       'gui_OutputFcn',  @batchviewer_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT

% --- Executes just before batchviewer is made visible.
function batchviewer_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to batchviewer (see VARARGIN)

    % Choose default command line output for batchviewer
    handles.output = hObject;

    % Make a batch
    if isempty(varargin), 
         handles.mybatch = dialogbox;
    else handles.mybatch = varargin{1};
    end

    % Update handles structure
    guidata(hObject, handles);

    % Check if batch is valid
    if isempty(handles.mybatch), disp('Empty batch provided, nothing to view in batch viewer'), return, end
    
    % set some default values
    handles.defaultchannelnames = {'channel1','channel2','channel3','channel4'};
    handles.defaultnrofchannels = numel(handles.defaultchannelnames);
    handles.plotpripos_small    = [86 30 130 18];
    handles.plotsecpos_small    = [86 7.5 130 18];
    handles.plotpripos_large    = [86 15 130 30]; 
    handles.activechannelnr     = [];
    handles.activechannelstr    = '';
    handles.buttongrp_channels.Children = handles.buttongrp_channels.Children(end:-1:1); % fix children to be numbered ascending
    
    % Initialise gui
    initialize_gui(hObject, handles)

    % wait for user response
    uiwait(handles.figure1);
    
% --- Outputs from this function are returned to the command line.
function varargout = batchviewer_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if isfield(handles, 'mybatch') && ~isempty(handles.mybatch)
         varargout{1} = filtermybatch(handles);
    else varargout{1} = [];
    end
    close(handles.figure1)

% --- Executes on selection change in listbox_abflist.
function listbox_abflist_Callback(hObject, eventdata, handles)
    % hObject    handle to listbox_abflist (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns listbox_abflist contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from listbox_abflist
    handles.abfidx = get(hObject,'Value');
    guidata(handles.figure1, handles);
    update_gui(handles)

% --- Executes during object creation, after setting all properties.
function listbox_abflist_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to listbox_abflist (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function initialize_gui(fig_handle, handles)
    % What to do upon initialisation.
    % populate list box with list of abf file names
    handles.abfnames = handles.mybatch.getabf.get('filename');
    if ischar(handles.abfnames), handles.abfnames = {handles.abfnames}; end
    set(handles.listbox_abflist,'String', handles.abfnames)
    set(handles.txt_username,   'String', sprintf('User: %s',handles.mybatch.userid));

    % starting value for current abf index in list
    handles.abfidx = 1;
    handles.abf    = handles.mybatch.getabf(handles.abfidx);
    
    % set starting value of channel radio button group to lowest available channel
    handles.activechannelnr = min(handles.abf.channelnrs);
    if isempty(handles.activechannelnr) || handles.activechannelnr > handles.defaultnrofchannels, 
        error('no channels present in abf "%s" or no valid channel numbers',handles.abf.filename); 
    end
    handles.activechannelstr = handles.defaultchannelnames{handles.activechannelnr};

    % set button group selected object to the active channel and store handles
    set(handles.buttongrp_channels,'SelectedObject',handles.buttongrp_channels.Children(handles.activechannelnr));
    
    % struct storing include/exclude status
    handles.abfcheck = false(handles.mybatch.nrofabfs,handles.defaultnrofchannels); % so thats nr of abfs by max nr of channels   
    for i=1:handles.mybatch.nrofabfs
        handles.abfcheck(i,handles.mybatch.getabf(i).channelnrs) = true;
    end
    
    % turn on toolbars for zooming, panning
    set(handles.figure1,'toolbar','figure');

    % Update handles structure and gui controls
    guidata(handles.figure1, handles);
    update_gui(handles)
    
function update_gui(handles)
    % update UI controls to current abffile
    handles.abf = handles.mybatch.getabf(handles.abfidx);
    set(handles.txt_header_abfinfo, 'String', handles.abf.filename);
    set(handles.txt_abfinfo,        'String', handles.abf.idmessage);
    set(handles.listbox_abflist,    'Value',  handles.abfidx);

    % include abf check box
    set(handles.checkbox_includeabf,'Value',iscurrentabfincluded(handles))

    % store present amplifier channels
    handles.presentchannelnames = handles.defaultchannelnames(handles.abf.channelnrs);
    
    % set check boxes according to selection
    % first hide all
    for i=1:numel(handles.defaultchannelnames)
        cstr = handles.defaultchannelnames{i};
        set(handles.(sprintf('radio_%s',cstr)),'Visible','off')
    end
    % now reveal applicable
    for i=1:numel(handles.abf.channelnrs)
        cstr = handles.presentchannelnames{i};
        set(handles.(sprintf('radio_%s',cstr)),'Visible','on')
    end
    
    % reset channel number if current active channel nr is empty or not available
    if isempty(handles.activechannelnr) || ~ismember(handles.activechannelnr,handles.abf.channelnrs)
        % set starting value of channel radio button group to lowest available channel
        handles.activechannelnr = min(handles.abf.channelnrs);
        if isempty(handles.activechannelnr) || handles.activechannelnr > handles.defaultnrofchannels, 
            error('no channels present in abf "%s" or no valid channel numbers',handles.abf.filename); 
        end
    else
        % find current channel from radiobutton group (dodgy stuff...)
        handles.activechannelstr = sscanf(get(handles.buttongrp_channels.SelectedObject,'Tag'),'radio_%s');
        handles.activechannelnr  = sscanf(handles.activechannelstr,'channel%d');
    end    
    
    % set button group selected object to the active channel
    handles.activechannelstr = handles.defaultchannelnames{handles.activechannelnr};
    set(handles.buttongrp_channels,'SelectedObject',handles.buttongrp_channels.Children(handles.activechannelnr));

    % update the check status of include channel checkbox
    set(handles.checkbox_includechannel,'Value',handles.abfcheck(handles.abfidx,handles.activechannelnr))

    % colour code excluded abfs
    if ~iscurrentabfincluded(handles)
        set(handles.txt_header_abfinfo,'ForegroundColor',[1 0 0])
    else
        set(handles.txt_header_abfinfo,'ForegroundColor',[0 0.4980 0])
    end
    
    % colour code excluded channels
    for i=1:numel(handles.abf.channelnrs)
        cstr = handles.presentchannelnames{i};
        if ~handles.abfcheck(handles.abfidx,handles.abf.channelnrs(i))
            set(handles.(sprintf('radio_%s',cstr)),'ForegroundColor',[1 0 0])
        else
            set(handles.(sprintf('radio_%s',cstr)),'ForegroundColor',[0 0.4980 0])
        end
    end

    % update handles structure
    guidata(handles.figure1, handles);

    % plot abf
    plotchannels(handles)

function plotchannels(handles)

    % first reset all to remove any lingering traces from previous display round
    cla(handles.axis_primary,'reset')
    cla(handles.axis_secondary,'reset')

    % plot primary
    adc = handles.abf.getchannel('number',handles.activechannelnr).getin('signal','primary');
    if isa(adc,'Analogin') && ~isempty(adc)        
        tsdata = adc.getdata;
        tstime = repmat(adc.getsweep(1).gettime,1,adc.nrofsweeps);                    
        plot(handles.axis_primary,tstime,tsdata)
        if adc.updatephase > 1
            configstring = sprintf('- %s (Channel %d)',adc.signal,handles.activechannelnr);
        else
            configstring = '';
        end
        title (handles.axis_primary,sprintf('Analog Input %d %s',adc.number,configstring))
        ylabel(handles.axis_primary,adc.units)
        xlabel(handles.axis_primary,'milliseconds')
        xlim  (handles.axis_primary,[adc.getsweep(1).TimeInfo.Start adc.getsweep(1).TimeInfo.End])
    end

    % plot secondary
    adc = handles.abf.getchannel('number',handles.activechannelnr).getin('signal','secondary');
    if isa(adc,'Analogin') && ~isempty(adc)
        
        % set axis positions and visibility
        set(handles.txt_SecondaryAbsent,'Visible','off')
        set(handles.axis_secondary,     'Visible','on')
        set(handles.axis_primary,  'Position',handles.plotpripos_small)
        set(handles.axis_secondary,'Position',handles.plotsecpos_small)
        
        tsdata = adc.getdata;
        tstime = repmat(adc.getsweep(1).gettime,1,adc.nrofsweeps);                    
        plot(handles.axis_secondary,tstime,tsdata)
        if adc.updatephase > 1
            configstring = sprintf('- %s (Channel %d)',adc.signal,handles.activechannelnr);
        else
            configstring = '';
        end
        title (handles.axis_secondary,sprintf('Analog Input %d %s',adc.number,configstring))
        ylabel(handles.axis_secondary,adc.units)
        xlabel(handles.axis_secondary,'milliseconds')
        xlim  (handles.axis_secondary,[adc.getsweep(1).TimeInfo.Start adc.getsweep(1).TimeInfo.End])
    else
        % set axis positions and visibility
        set(handles.txt_SecondaryAbsent,'Visible','on')
        set(handles.axis_secondary,'Visible','off')
        set(handles.axis_primary,'Position',handles.plotpripos_large)
    end
    
    % link axes
    linkaxes([handles.axis_primary,handles.axis_secondary],'x')


% Callbacks ###################################################################

% --- Executes on button press in checkbox_includeabf.
function checkbox_includeabf_Callback(hObject, eventdata, handles)
    % hObject    handle to checkbox_includeabf (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of checkbox_includeabf
    
    % So for the row of the abf in the truth table, assign the columns corresponding to present channels in current abf to
    % the logical value of the includeabf checkbox status:
    handles.abfcheck(handles.abfidx,handles.abf.channelnrs) = repmat(get(hObject,'Value'),1,numel(handles.abf.channelnrs));
    guidata(handles.figure1, handles);
    update_gui(handles)

% --- Executes on button press in checkbox_includechannel.
function checkbox_includechannel_Callback(hObject, eventdata, handles)
    % hObject    handle to checkbox_includechannel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of checkbox_includechannel
    
    % So for the row of the abf in the truth table, assign the column corresponding to active channel number to
    % the logical value of the includechannel checkbox status:
    handles.abfcheck(handles.abfidx,handles.activechannelnr) = get(hObject,'Value');
    guidata(handles.figure1, handles);
    update_gui(handles)

% --- Executes on button press in button_previous.
function button_previous_Callback(hObject, eventdata, handles)
    % hObject    handle to button_previous (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if handles.abfidx>1,handles.abfidx = handles.abfidx-1; end
    guidata(hObject,handles)
    update_gui(handles)

% --- Executes on button press in button_next.
function button_next_Callback(hObject, eventdata, handles)
    % hObject    handle to button_next (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if handles.abfidx<handles.mybatch.nrofabfs, handles.abfidx = handles.abfidx+1; end
    update_gui(handles)
    guidata(hObject,handles)

% --- Executes on button press in button_viewallchannels.
function button_viewallchannels_Callback(hObject, eventdata, handles)
    % hObject    handle to button_viewallchannels (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    figure()
    fh_abf = handles.abf.plot;
    set(fh_abf,'position',[ 680   160   560   820])

% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)
    % hObject    handle to button_done (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    uiresume(handles.figure1)

% --- Executes when selected object is changed in buttongrp_channels.
function buttongrp_channels_SelectionChangedFcn(hObject, eventdata, handles)
    % hObject    handle to the selected object in buttongrp_channels 
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    update_gui(handles)

% ########### Helpers ##############
function b = dialogbox
    % Construct a questdlg with three options
    choice = questdlg('No input provided. What would you like to do?', ...
        'No input', ...
        'New batch','Load batch','Quit','New batch');
    % Handle response
    switch choice
        case 'New batch'
            b = newbatch;
        case 'Load batch'
            b = loadbatch;
        case 'Quit'
            b = Abfbatch; % creates an empty ABFBATCH
    end

function b = newbatch
    % make a new ABFBATCH object using gui
    b = Abfbatch('gui');

function b = loadbatch
    % navigate to and load an existing ABFBATCH object
    [batchfile, path2batch] = uigetfile('*.mat','Select Abfbatch object','MultiSelect','off');  
    if ischar(batchfile) && ~isempty(batchfile)
        b = load(fullfile(path2batch,batchfile));
        b = b.obj;
        if ~isa(b,'Abfbatch'), error('Chosen file is not an Abfbatch object'), end
    else
        b = Abfbatch; % creates an empty ABFBATCH
    end

function b = filtermybatch(handles)
    % remove all abfs and abf channels that were not inluded
    b        = handles.mybatch;
    idx2keep = any(handles.abfcheck,2);
    b        = b.selectabf(idx2keep);
    handles.abfcheck = handles.abfcheck(idx2keep,:); 
    
    % remove unwanted abfs and channels
    chnnls = 1:handles.defaultnrofchannels;
    for i  = 1:size(handles.abfcheck,1)
        abf  = b.getabf(i);
        lose = abf.channelnrs(ismember(abf.channelnrs,chnnls(~handles.abfcheck(i,:))));
        for ii = 1:numel(lose)
            b.abfs(i) = abf.removechannel('number',lose(ii));
        end
    end   

function i = iscurrentabfincluded(handles)
    % check if current abf is included. The include status of entire abf depends on there being at least one channel included 
    % (that is, having a logical value true in the abfcheck truth table).
    i = any(handles.abfcheck(handles.abfidx,:));    
    
% --- Executes on button press in button_savebatch.
function button_savebatch_Callback(hObject, eventdata, handles)
    % hObject    handle to button_savebatch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    mybatchfilt = filtermybatch(handles); % allow selection to take effect
    path2file   = fullfile(mybatchfilt.path2abfs,'batch',mybatchfilt.savename);
    [filename,path2file] = uiputfile(path2file,'Save Abfbatch object');
    if path2file == 0,
        disp('No path selected, nothing saved')
        return
    end
    mybatchfilt.saveme(path2file,filename);

% --- Executes on button press in button_saveabfs.
function button_saveabfs_Callback(hObject, eventdata, handles)
    % hObject    handle to button_saveabfs (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    mybatchfilt = filtermybatch(handles); % allow selection to take effect
    path2file   = uigetdir(mybatchfilt.path2abfs,'Select folder to store Abffile objects');
    if path2file == 0,
        disp('No path selected, nothing saved')
        return
    end
    for i=1:mybatchfilt.nrofabfs
        mybatchfilt.getabf(i).saveme(path2file,mybatchfilt.getabf(i).savename);
    end

% --- Executes on button press in button_savecurrentabf.
function button_savecurrentabf_Callback(hObject, eventdata, handles)
    % hObject    handle to button_savecurrentabf (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    mybatchfilt = filtermybatch(handles); % allow selection to take effect
    path2file   = fullfile(mybatchfilt.path2abfs,'mat',mybatchfilt.getabf(handles.abfidx).savename);
    [filename,path2file] = uiputfile(path2file,'Save current Abffile object');
    if path2file == 0,
        disp('No path selected, nothing saved')
        return
    end
    mybatchfilt.getabf(handles.abfidx).saveme(path2file,filename);
