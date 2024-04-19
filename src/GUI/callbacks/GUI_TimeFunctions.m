%% Dynamic Nodal Loads Dialog Callback Functions
% This file contains the callback functions associated with the "Dynamic
% Nodal Loads" dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
%#ok<*DEFNU>
function varargout = GUI_TimeFunctions(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_TimeFunctions_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_TimeFunctions_OutputFcn, ...
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


% --- Executes just before GUI_TimeFunctions is made visible.
function GUI_TimeFunctions_OpeningFcn(hObject, ~, handles, varargin)
include_constants;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Set callback new pushbutton
set(handles.pushbutton_New,'Callback',@pushbutton_New_Callback);

% Set callback delete pushbutton
set(handles.pushbutton_Delete,'Callback',@pushbutton_Delete_Callback);

% Set callback create pushbutton
set(handles.pushbutton_CreateLFcn,'Callback',@pushbutton_CreateLFcn_Callback);

% Set callback Apply pushbutton
set(handles.pushbutton_Apply,'Callback',@pushbutton_Apply_Callback);
            
% Get model object from root
model = getappdata(0,'model');

if isempty(model.strTimeFcns)
    cEdit = [ false(1,2) true(1,6) ];
    cFormat = {[] [] [] [] [] [] [] 'logical'};
    set(handles.uitable_LoadFcn,'Data',{},...
        'CellEditCallback',@uitable_LoadFcn_CellEditCallback,...
        'CellSelectionCallback',@uitable_LoadFcn_CellSelectionCallback,...
        'ColumnEditable',cEdit,'ColumnFormat',cFormat,...
        'Enable','off');
    % Choose default command line output for GUI_TimeFunctions
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    return
end
    
% Set popupmenu_FcnList properties
set(handles.popupmenu_FcnList,'Enable','on','value',1,'string',model.strTimeFcns,...
    'max',length(model.strTimeFcns));

% Enable delete pushbutton
set(handles.pushbutton_Delete,'Enable','On');

% Set popupmenu_LoadFunction properties
set(handles.popupmenu_LoadFunction,'Enable','on','value',1);

% Factor
set(handles.edit_Amplitude,'Enable','on','String',num2str(1));

% Initial time
set(handles.edit_Ti,'Enable','on','String',num2str(0));

% Final time
set(handles.edit_dT,'Enable','on','String',num2str(0));

% Frequency
set(handles.edit_Freq,'Enable','on','String',num2str(0));

% Phase
set(handles.edit_Phase,'Enable','on','String',num2str(0));

% Enable create pushbutton
set(handles.pushbutton_CreateLFcn,'Enable','On');

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(handles.axes_LFcn);
step = model.t/model.n_steps;

if ~isempty(model.timeFcns{1})
    x = 0:step:model.t;
    y = model.timeFcns{1}.evalAll;
else
    x = [0, model.t];
    y = [0, 0];  
end

ymax  = max(y);
ymin  = min(y);
ymean = (ymax + ymin) / 2;
ylen  = ymax - ymin;
if ~ylen
    ylen = 1 / 0.6;
end
plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
ylim([ymean - 0.6 * ylen, ymean + 0.6 * ylen]);
grid on

% Set table data
if ~isempty(model.timeFcns{1})
    enable_table = 'On';
    
    % Get number of time functions
    [~,nFcn] = model.timeFcns{1}.goThrough();
    
    % Predimension cell to hold time function info
    lfcn = cell(nFcn,8);
    
    ptr =  model.timeFcns{1};
    count = 1;
    while ~isempty(ptr) && count <= nFcn
        switch ptr.type
            case PERIODIC
                lfcn(count,:) = {ptr.id, 'Harmonic', ptr.weightFctr, ptr.w, ptr.phi, ptr.ti, ptr.tf,false};
            case STATIC
                lfcn(count,:) = {ptr.id, 'Constant', ptr.weightFctr, ' ', ' ', ptr.ti, ptr.tf,false};
            case SLOPE
                lfcn(count,:) = {ptr.id, 'Slope', ptr.weightFctr, ' ', ' ', ptr.ti, ptr.tf,false};
            case TABLE
                lfcn(count,:) = {ptr.id, 'Table', ptr.weightFctr, ' ', ' ',ptr.ti, ptr.tf,false};
        end
        ptr = ptr.next;
        count = count + 1;
    end
else
    enable_table = 'Off';
    lfcn = {};
end
cEdit = [ false(1,2) true(1,6) ];
cFormat = {[] [] [] [] [] [] [] 'logical'};
set(handles.uitable_LoadFcn,'Data',lfcn,...
    'CellEditCallback',@uitable_LoadFcn_CellEditCallback,...
    'CellSelectionCallback',@uitable_LoadFcn_CellSelectionCallback,...
    'ColumnEditable',cEdit,'ColumnFormat',cFormat,...
    'Enable',enable_table);

% Disable Apply pushbutton
set(handles.pushbutton_Apply,'Enable','Off');

% Choose default command line output for GUI_TimeFunctions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_TimeFunctions_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_New.
function pushbutton_New_Callback(hObject, ~, ~)
include_constants;
model = getappdata(0,'model');
set(hObject,'enable','off');

gui = guidata(findobj('Tag','GUI_TimeFunctions'));

% Define time function list name
name = char(inputdlg('Provide a name for the new time function','New time function',[1,50]));

if isempty(name)
    set(hObject,'enable','on');
    return
end

new_id = model.createFcnList(name);
setappdata(0,'model',model);

% Set popupmenu_LoadFunction properties
set(gui.popupmenu_FcnList,'Enable','on','value',new_id,'max',new_id,'string',model.strTimeFcns);

% Enable delete pushbutton
set(gui.pushbutton_Delete,'Enable','On');

% Set popupmenu_LoadFunction properties
set(gui.popupmenu_LoadFunction,'Enable','on','value',1);

% Factor
set(gui.edit_Amplitude,'Enable','on','String',num2str(1));

% Initial time
set(gui.edit_Ti,'Enable','on','String',num2str(0));

% Final time
set(gui.edit_dT,'Enable','on','String',num2str(0));

% Frequency
set(gui.edit_Freq,'Enable','on','String',num2str(0));

% Phase
set(gui.edit_Phase,'Enable','on','String',num2str(0));

% Enable create pushbutton
set(gui.pushbutton_CreateLFcn,'Enable','On');

popupmenu_FcnList_Callback(gui.popupmenu_FcnList,[],gui)

set(hObject,'enable','on');

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Delete.
function pushbutton_Delete_Callback(hObject, ~, ~)
include_constants;

set(hObject,'enable','off');
mdata = guidata(findobj('Tag','GUI_Main'));
gui  = guidata(findobj('Tag','GUI_TimeFunctions'));

model = getappdata(0,'model');
fcn_id = get(gui.popupmenu_FcnList,'value');

%Deleting fcn list
model.deleteFcnList(fcn_id);

% Return model to root
setappdata(0,'model',model);

shouldReenableButton = true;
if get(gui.popupmenu_FcnList,'max') == 1  
    shouldReenableButton = false;
    
    % Set popupmenu_FcnList properties
    set(gui.popupmenu_FcnList,'Enable','off','value',1,'string',{' '},'max',1);

    % Enable delete pushbutton
    set(gui.pushbutton_Delete,'Enable','off');

    % Set popupmenu_LoadFunction properties
    set(gui.popupmenu_LoadFunction,'Enable','off','value',1);

    % Factor
    set(gui.edit_Amplitude,'Enable','off','String',' ');

    % Initial time
    set(gui.edit_Ti,'Enable','off','String',' ');

    % Final time
    set(gui.edit_dT,'Enable','off','String',' ');

    % Frequency
    set(gui.edit_Freq,'Enable','off','String',' ');

    % Phase
    set(gui.edit_Phase,'Enable','off','String',' ');

    % Enable create pushbutton
    set(gui.pushbutton_CreateLFcn,'Enable','off');

else
    
    % Set popupmenu_FcnList properties
    val = get(gui.popupmenu_FcnList,'value');
    val = val - (val==get(gui.popupmenu_FcnList,'max'));

    set(gui.popupmenu_FcnList,'value',val,'string',model.strTimeFcns,'max',length(model.strTimeFcns));
end

popupmenu_FcnList_Callback(gui.popupmenu_FcnList,[],gui)

% Enable "Process Data" button in main GUI
set(mdata.pushbutton_ProcessData,'Enable','on');

% Disable result buttons
allLoadsNeedToBeDrawn = false;
if get(mdata.popupmenu_Results,'Value') ~= 1
    allLoadsNeedToBeDrawn = true;
end
set(mdata.popupmenu_Results,'Enable','off','value',1);
set(mdata.pushbutton_Textual,'Enable','off');
set(mdata.text_Element,'string','Elements');
set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(mdata.edit_Scale,'enable','off','visible','off');
set(mdata.pushbutton_DynamicResults,'enable','off');
set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);

% Disable/Enable visualization options
set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);

% Make GUI a normal window
this_fig = findobj('Tag','GUI_TimeFunctions');
set(this_fig,'WindowStyle','normal');

% Draw updated model
if allLoadsNeedToBeDrawn, redraw(mdata,'Loads'); end

% Make GUI a modal window
set(this_fig,'WindowStyle','modal');

if shouldReenableButton, set(hObject,'enable','on'); end

% --- Executes on selection change in popupmenu_LoadFunction.
function popupmenu_FcnList_Callback(hObject, ~, handles)
include_constants;

model=getappdata(0,'model');

if strcmp(get(hObject,'Enable'),'off')
    % Update axes
    delete(findobj('Tag','DrawLfcn'))
    axes(handles.axes_LFcn);
    
    x = [0, model.t];
    y = [0, 0];
    
    plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
    ylim([- 0.6, 0.6]);
    grid on
    
    set(handles.uitable_LoadFcn,'Data',{},'Enable','off');
    return
end

fcnId = get(hObject,'value');

% Set table data
if ~isempty(model.timeFcns{fcnId})
    enable_table = 'On';
    
    % Get number of time functions
    [~,nFcn] = model.timeFcns{fcnId}.goThrough();
    
    % Predimension cell to hold time function info
    lfcn = cell(nFcn,8);
    
    ptr =  model.timeFcns{fcnId};
    count = 1;
    while ~isempty(ptr) && count <= nFcn
        switch ptr.type
            case PERIODIC
                lfcn(count,:) = {ptr.id,'Harmonic', ptr.weightFctr, ptr.w, ptr.phi,  ptr.ti, ptr.tf,false};
            case STATIC
                lfcn(count,:) = {ptr.id,'Constant', ptr.weightFctr, ' ', ' ', ptr.ti, ptr.tf,false};
            case SLOPE
                lfcn(count,:) = {ptr.id,'Slope', ptr.weightFctr, ' ', ' ', ptr.ti, ptr.tf,false};
            case TABLE
                lfcn(count,:) = {ptr.id, 'Table', ptr.weightFctr, ' ', ' ', ptr.ti, ptr.tf,false};
        end
        ptr = ptr.next;
        count = count + 1;
    end
else
    enable_table = 'Off';
    lfcn = {};
end
set(handles.uitable_LoadFcn,'Data',lfcn,'Enable',enable_table);

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(handles.axes_LFcn);
step = model.t/model.n_steps;

if ~isempty(model.timeFcns{fcnId})
    x = 0:step:model.t;
    y = model.timeFcns{fcnId}.evalAll;
else
    x = [0, model.t];
    y = [0, 0];  
end

ymax  = max(y);
ymin  = min(y);
ymean = (ymax + ymin) / 2;
ylen  = ymax - ymin;

if ~ylen
    ylen = 1 / 0.6;
end

plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
ylim([ymean - 0.6 * ylen, ymean + 0.6 * ylen]);
grid on

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_CreateLFcn.
function pushbutton_CreateLFcn_Callback(hObject, ~, ~)
include_constants;

mdata = guidata(findobj('Tag','GUI_Main'));
gui   = guidata(findobj('Tag','GUI_TimeFunctions'));

model = getappdata(0,'model');

% Disable button while input data is being set
set(hObject,'enable','off')

% Get current function list id
fcnId = get(gui.popupmenu_FcnList,'Value');

% Get load factor
fctr = str2double(get(gui.edit_Amplitude,'String'));
if isnan(fctr)
    fctr = 0;
end

% Get load frequency values
freq = str2double(get(gui.edit_Freq,'String'));
if isnan(freq)
    freq = 0;
end

% Get load phase angle values
phi = str2double(get(gui.edit_Phase,'String'));
if isnan(phi)
    phi = 0;
end

% Get load initial time
ti = str2double(get(gui.edit_Ti,'String'));
if isnan(ti) || ti < 0
    ti = 0;
elseif ti > model.t
    % Enable button for futre use
    set(hObject,'enable','on')

    msgbox('Invalid input data! Initial time must be within analysis interval.', 'Error','error');
    return
end

% Get load final time
tf = str2double(get(gui.edit_dT,'String'));
if isnan(tf) 
    tf = ti;
end

% Get load duration
dt = tf-ti;
if  dt < 0
    dt = 0;
elseif ti + dt > model.t
    % Enable button for futre use
    set(hObject,'enable','on')
    msgbox('Invalid input data! Time function duration must be within analysis interval.', 'Error','error');
    return
end

% Get model fcn ids
if isempty(model.timeFcns)
    new_fcn_id = 1;
elseif isempty(model.timeFcns{fcnId})
    new_fcn_id = 1;
else
    id = model.timeFcns{fcnId}.getIds();
    if ~isempty(id)
        new_fcn_id = max(id) + 1;
    else
        new_fcn_id = 1;
    end
    clear id
end

% Set new time function to dynamic load
if  get(gui.popupmenu_LoadFunction,'value') == 1      % Periodic
    model.addFcn(fcnId,PERIODIC,{new_fcn_id,fctr,freq,phi,ti,ti+dt});
    type = 'Harmonic';
elseif get(gui.popupmenu_LoadFunction,'value') == 2  % Static
    model.addFcn(fcnId,STATIC,{new_fcn_id,fctr,ti,ti+dt});
    type = 'Constant';
    freq = ' ';
    phi  = ' ';
elseif get(gui.popupmenu_LoadFunction,'value') == 3  % Slope
    model.addFcn(fcnId,SLOPE,{new_fcn_id,fctr,ti,ti+dt});
    type = 'Slope';
    freq = ' ';
    phi  = ' ';
elseif get(gui.popupmenu_LoadFunction,'value') == 4  % Table
    
    filterspec = {'*.txt'};
    [filename,pathname] = uigetfile(filterspec,'LESM - Get Time Table');
    if isequal(filename,0)
        set(hObject,'enable','on')
        return;
    end
    [x,F] = read_time_table(strcat(pathname,filename));
    if ~all(x>=0)
        F(x<0) = [];
        x(x<0) = [];
        uiwait(msgbox('Table has negative time domain values. Only positive time values were considered.','Warning','warn','modal'));
    end    
     
    model.addFcn(fcnId,TABLE,{new_fcn_id,fctr,x,F});
    [ptr,~] = model.timeFcns{fcnId}.goThrough();
    ptr.src_file = strcat(pathname,filename);
    ti = x(1);
    if ti < 0
        ti = 0;
    end
    tf = x(end);
    if tf > model.t
        tf = model.t;
    end
    type = 'Table';
    freq = ' ';
    phi  = ' ';
end

% Set model object properties
setappdata(0,'model',model);

% Update uitable data
tableData = get(gui.uitable_LoadFcn,'Data');
if isempty(tableData)    
    set(gui.uitable_LoadFcn,'Enable','On');
    tableData = {new_fcn_id,type,fctr,freq,phi,ti,tf,false};
else
    tableData(end+1,:) = {new_fcn_id,type,fctr,freq,phi,ti,tf,false};
end
set(gui.uitable_LoadFcn,'Data',tableData);

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(gui.axes_LFcn);
step = model.t/model.n_steps;

if ~isempty(model.timeFcns{fcnId})
    x = 0:step:model.t;
    y = model.timeFcns{fcnId}.evalAll;
else
    x = [0, model.t];
    y = [0, 0];  
end

ymax  = max(y);
ymin  = min(y);
ymean = (ymax + ymin) / 2;
ylen  = ymax - ymin;

if ~ylen
    ylen = 1 / 0.6;
end

plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
ylim([ymean - 0.6 * ylen, ymean + 0.6 * ylen]);
grid on

% Enable "Process Data" button in main GUI
set(mdata.pushbutton_ProcessData,'Enable','on');

% Disable model type option
set(mdata.popupmenu_Anm,'Enable','off');

% Disable result buttons
allLoadsNeedToBeDrawn = false;
if get(mdata.popupmenu_Results,'Value') ~= 1
    allLoadsNeedToBeDrawn = true;
end
set(mdata.popupmenu_Results,'Enable','off','value',1);
set(mdata.pushbutton_Textual,'Enable','off');
set(mdata.text_Element,'string','Elements');
set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(mdata.edit_Scale,'enable','off','visible','off');
set(mdata.pushbutton_DynamicResults,'enable','off');
set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);

% Disable/Enable visualization options
set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);

% Return variables to root
setappdata(0,'resultType',0);
setappdata(0,'model',model);

% Make GUI a normal window
gui_dialog = findobj('Tag','GUI_TimeFunctions');
set(gui_dialog,'WindowStyle','normal');

% Draw updated model
if allLoadsNeedToBeDrawn
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui_dialog,'WindowStyle','modal');

% Enable button for futre use
set(hObject,'enable','on')

%--------------------------------------------------------------------------
function uitable_LoadFcn_CellEditCallback(~, eventdata, ~)
% Get selected cell on uitable
if isempty(eventdata.Indices)
    return
elseif size(eventdata.Indices,1) > 1
    id = [eventdata.Indices(end,1),1];
else
    id = eventdata.Indices;
end

% Check if edited cell was on the delete column
if id(2) == 8
    return
end

% Get handle to this figure
gui = guidata(findobj('tag','GUI_TimeFunctions'));

% Enable Apply pushbutton
set(gui.pushbutton_Apply,'Enable','On');

%--------------------------------------------------------------------------
function uitable_LoadFcn_CellSelectionCallback(~,eventdata,~)
% Get selected cell on uitable
if isempty(eventdata.Indices)
    return
elseif size(eventdata.Indices,1) > 1
    id = [eventdata.Indices(end,1),1];
else
    id = eventdata.Indices;
end

% Check if click was on the delete column
if id(2) ~= 8
    return
end

% Get handle to this figure
gui = guidata(findobj('tag','GUI_TimeFunctions'));

% Get table data
tableData = get(gui.uitable_LoadFcn,'Data');

if id(1) == 1 && size(tableData,1) > 1
    choice = questdlg('Deleting the first function will erase the whole function table. Continue?',...
        'Warning!','Yes','No','No');
    switch choice
        case 'No'
            tableData{1,8} = false;
            set(gui.uitable_LoadFcn,'Data',tableData);
            return
    end
end

% Disable uitable while function is being deleted
set(gui.uitable_LoadFcn,'Enable','Off');

% Get handle to model
model = getappdata(0,'model');

% Get current function list id
n = get(gui.popupmenu_FcnList,'value');

% Get id of fcn to be deleted
fcn_id = tableData{id(1),1};

% Remove selected function from list
model.removeFcn(n,fcn_id);

% Remove function info from table data
% If the deleted function is the first, the whole pattern is deleted
% (avoids errors)
if id(1) == 1 
    tableData = [];
else 
    tableData(id(1),:) = [];
end

% Set updated table data to uitable
set(gui.uitable_LoadFcn,'Data',tableData);

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(gui.axes_LFcn);
step = model.t/model.n_steps;

if ~isempty(model.timeFcns{n})
    x = 0:step:model.t;
    y = model.timeFcns{n}.evalAll;
else
    x = [0, model.t];
    y = [0, 0];  
end

ymax  = max(y);
ymin  = min(y);
ymean = (ymax + ymin) / 2;
ylen  = ymax - ymin;
if ~ylen
    ylen = 1 / 0.6;
end
plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
ylim([ymean - 0.6 * ylen, ymean + 0.6 * ylen]);
grid on

% Check if there still are other time functions associated to node.
% If so, enable uitable for future use.
if ~isempty(tableData)
    set(gui.uitable_LoadFcn,'Enable','On');
else
    % Disable Apply pushbutton
    set(gui.pushbutton_Apply,'Enable','Off');
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(hObject, ~, ~)
% Disable pushbutton
set(hObject,'Enable','Off');

% Get handle to this figure
gui = guidata(findobj('tag','GUI_TimeFunctions'));

% Get handle to model
model = getappdata(0,'model');

% Get current function list id
n = get(gui.popupmenu_FcnList,'value');

% Get uitable data
tableData = get(gui.uitable_LoadFcn,'Data');

% Loop through all rows of tableData
for i = 1:size(tableData,1)
    % Get handle to Lfcn object to be changed
    fcn_id = tableData{i,1};
    ptr = model.timeFcns{n}.getById(fcn_id);
    
    % Set new info to Lfcn object, if input is valid
    % Function weight factor
    fctr = tableData{i,3};
    if ~isnan(fctr)
        ptr.weightFctr = fctr;
    else
        fctr = ptr.weightFctr;
    end
    
    % Initial instant
    ti = tableData{i,6};
    if isnan(ti)
        ti = ptr.ti;
    elseif ti < 0 || ti > model.t
        ti = ptr.ti;
    else
        ptr.ti = ti;
    end
    
    % Final Instant
    tf = tableData{i,7};
    if isnan(tf)
        tf = ptr.tf;
    elseif tf < 0
        tf = ptr.tf;
    elseif tf > model.t
        ptr.tf = model.t;
        tf = ptr.tf;
    else
        ptr.tf = tf;
    end
    
    % Check for Lfcn type
    type   = tableData{i,2};
    if strcmp(type,'Harmonic')
        % Angular speed
        w   = tableData{i,4};
        if ~isnan(w)
            ptr.w = w;
        else
            w = ptr.w;
        end
        
        % Phase angle
        phi = tableData{i,5};
        if ~isnan(phi)
            ptr.phi = phi;
        else
            phi = ptr.phi;
        end
    elseif strcmp(type,'Table')
        f_totalT = ptr.x(end) - ptr.x(1);
        if tf > ti + f_totalT
            tf = ti + f_totalT;
            ptr.tf = tf;
        end
        w   = ' ';
        phi = ' ';
    else
        w   = ' ';
        phi = ' ';
    end
    
    % Update tableData
    tableData(i,:) = {fcn_id,type,fctr,w,phi,ti,tf,false};
    
    % Compute new time function values
    ptr.evalFcn();
end

% Reset uitable data
set(gui.uitable_LoadFcn,'Data',tableData);

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(gui.axes_LFcn);
step = model.t/model.n_steps;

if ~isempty(model.timeFcns{n})
    x = 0:step:model.t;
    y = model.timeFcns{n}.evalAll;
else
    x = [0, model.t];
    y = [0, 0];  
end

ymax  = max(y);
ymin  = min(y);
ymean = (ymax + ymin) / 2;
ylen  = ymax - ymin;
if ~ylen
    ylen = 1 / 0.6;
end
plot(x,y,'Tag','DrawLfcn','LineWidth',1.1);
ylim([ymean - 0.6 * ylen, ymean + 0.6 * ylen]);
grid on

% Get handle to GUI_main
mdata = guidata(findobj('Tag','GUI_Main'));

% Enable "Process Data" button in main GUI
set(mdata.pushbutton_ProcessData,'Enable','on');

% Disable result buttons
allLoadsNeedToBeDrawn = false;
if get(mdata.popupmenu_Results,'Value') ~= 1
    allLoadsNeedToBeDrawn = true;
end
set(mdata.popupmenu_Results,'Enable','off','value',1);
set(mdata.pushbutton_Textual,'Enable','off');
set(mdata.text_Element,'string','Elements');
set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(mdata.edit_Scale,'enable','off','visible','off');
set(mdata.pushbutton_DynamicResults,'enable','off');
set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);

% Disable/Enable visualization options
set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);

% Return model to root
setappdata(0,'model',model);

% Make GUI a normal window
gui = findobj('Tag','GUI_TimeFunctions');
set(gui,'WindowStyle','normal');

% Draw updated model
if allLoadsNeedToBeDrawn
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');


%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
% function popupmenu_FcnList_CreateFcn(hObject, ~, ~)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% --- Executes on selection change in popupmenu_LoadFunction.
function popupmenu_LoadFunction_Callback(~, ~, handles)
gui   = guidata(findobj('Tag','GUI_TimeFunctions'));

% Factor
    set(handles.edit_Amplitude,'Enable','on','String',num2str(1));

if get(gui.popupmenu_LoadFunction,'value') == 1
    
      
    % Initial time
    set(handles.edit_Ti,'Enable','on','String',num2str(0));
    
    % Duration
    set(handles.edit_dT,'Enable','on','String',num2str(0));
    
    % Frequency
    set(handles.edit_Freq,'Enable','on','String',num2str(0));
    
    % Phase
    set(handles.edit_Phase,'Enable','on','String',num2str(0));
    
    % Enable create pushbutton
    set(handles.pushbutton_CreateLFcn,'Enable','On','Callback',@pushbutton_CreateLFcn_Callback);
    
elseif get(gui.popupmenu_LoadFunction,'value') == 2
    
     
    % Initial time
    set(handles.edit_Ti,'Enable','on','String',num2str(0));
    
    % Duration
    set(handles.edit_dT,'Enable','on','String',num2str(0));
    
    % Frequency
    set(handles.edit_Freq,'Enable','off','String','');
    
    % Phase
    set(handles.edit_Phase,'Enable','off','String','');
    
    % Enable create pushbutton
    set(handles.pushbutton_CreateLFcn,'Enable','On','Callback',@pushbutton_CreateLFcn_Callback);
    
elseif get(gui.popupmenu_LoadFunction,'value') == 3
      
    % Initial time
    set(handles.edit_Ti,'Enable','on','String',num2str(0));
    
    % Duration
    set(handles.edit_dT,'Enable','on','String',num2str(0));
    
    % Frequency
    set(handles.edit_Freq,'Enable','off','String','');
    
    % Phase
    set(handles.edit_Phase,'Enable','off','String','');
    
    % Enable create pushbutton
    set(handles.pushbutton_CreateLFcn,'Enable','On','Callback',@pushbutton_CreateLFcn_Callback);
    
elseif get(gui.popupmenu_LoadFunction,'value') == 4
    
   
    % Initial time
    set(handles.edit_Ti,'Enable','off','String','');
    
    % Duration
    set(handles.edit_dT,'Enable','off','String','');
    
    % Frequency
    set(handles.edit_Freq,'Enable','off','String','');
    
    % Phase
    set(handles.edit_Phase,'Enable','off','String','');
    
    % Enable create pushbutton
    set(handles.pushbutton_CreateLFcn,'Enable','On','Callback',@pushbutton_CreateLFcn_Callback);
end


% --- Executes during object creation, after setting all properties.
function popupmenu_LoadFunction_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_dT_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_dT_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Ti_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Ti_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Phase_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Phase_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Freq_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Freq_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


