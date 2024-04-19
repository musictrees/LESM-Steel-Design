%% Nodal Load Dialog Callback Functions
% This file contains the callback functions associated with the "Nodal Loads"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_NodalLoads(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_NodalLoads_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_NodalLoads_OutputFcn, ...
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

%--------------------------------------------------------------------------
% Executes just before NodalLoads GUI is made visible.
% Sets GUI initial properties.
function GUI_NodalLoads_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_NodalLoads
handles.output = hObject;

% Move GUI to the center of the screen
if getappdata(0,'move') == 1
    movegui(hObject,'center');
end

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_NodalLoads_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes during popupmenu_Node creation, after setting all properties.
function popupmenu_Node_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Create list of nodes
nnp = getappdata(0,'nnp');
n = zeros(1,nnp);
for i = 1:nnp
    n(i) = i;
end
n = num2str(n,'%d\n');
set(hObject,'string',n)

%--------------------------------------------------------------------------
% Executes on popupmenu_Node selection change.
function popupmenu_Node_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));

% Get ID of selected load case
lc = get(mdata.popupmenu_LoadCase,'Value');

% Get ID of selected node
n = get(hObject,'Value');

% Get load components
nodes = getappdata(0,'nodes');
if lc > size(nodes(n).nodalLoadCase,2)
    fx = num2str(0);
    fy = num2str(0);
    fz = num2str(0);
    mx = num2str(0);
    my = num2str(0);
    mz = num2str(0);
else    
    fx = num2str(nodes(n).nodalLoadCase(1,lc));
    fy = num2str(nodes(n).nodalLoadCase(2,lc));
    fz = num2str(nodes(n).nodalLoadCase(3,lc));
    mx = num2str(nodes(n).nodalLoadCase(4,lc));
    my = num2str(nodes(n).nodalLoadCase(5,lc));
    mz = num2str(nodes(n).nodalLoadCase(6,lc));
end    

% Show load components of selected node
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 1
    set(handles.edit_Fx,'String',fx)
    set(handles.edit_Fy,'String',fy)
    set(handles.edit_Fz,'String','')
    set(handles.edit_Mx,'String','')
    set(handles.edit_My,'String','')
    set(handles.edit_Mz,'String','')
elseif anm == 2
    set(handles.edit_Fx,'String',fx)
    set(handles.edit_Fy,'String',fy)
    set(handles.edit_Fz,'String','')
    set(handles.edit_Mx,'String','')
    set(handles.edit_My,'String','')
    set(handles.edit_Mz,'String',mz)
elseif anm == 3
    set(handles.edit_Fx,'String','')
    set(handles.edit_Fy,'String','')
    set(handles.edit_Fz,'String',fz)
    set(handles.edit_Mx,'String',mx)
    set(handles.edit_My,'String',my)
    set(handles.edit_Mz,'String','')
elseif anm == 4
    set(handles.edit_Fx,'String',fx)
    set(handles.edit_Fy,'String',fy)
    set(handles.edit_Fz,'String',fz)
    set(handles.edit_Mx,'String','')
    set(handles.edit_My,'String','')
    set(handles.edit_Mz,'String','')
elseif anm == 5
    set(handles.edit_Fx,'String',fx)
    set(handles.edit_Fy,'String',fy)
    set(handles.edit_Fz,'String',fz)
    set(handles.edit_Mx,'String',mx)
    set(handles.edit_My,'String',my)
    set(handles.edit_Mz,'String',mz)
end

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Set.
% Sets nodal load properties of a Node object.
function pushbutton_Set_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
nodes = getappdata(0,'nodes');

% Disable button while input data is being set
set(hObject,'enable','off')

% Get ID of selected load case
lc = get(mdata.popupmenu_LoadCase,'Value');

% Get ID of selected node
if ~get(handles.checkbox_MultiNodes,'value')
    n_ID = get(handles.popupmenu_Node,'Value');
else
    n_str = get(handles.edit_MultiNodes,'string');
    if strcmp(n_str,'all')
        n_ID = 1:getappdata(0,'nnp');
    elseif strcmp(n_str,'ALL')
        n_ID = 1:getappdata(0,'nnp');
    elseif strcmp(n_str,'All')
        n_ID = 1:getappdata(0,'nnp');
    else
        [flag,n_ID] = readStr(n_str,getappdata(0,'nnp'));
        if ~flag
            % Enable button for futre use
            set(hObject,'enable','on')
            msgbox('Invalid input data!', 'Error','error');
            return
        end
    end
end

% Get load component values
fx = str2double(get(handles.edit_Fx,'String'));
if isnan(fx)
    fx = 0;
end
fy = str2double(get(handles.edit_Fy,'String'));
if isnan(fy)
    fy = 0;
end
fz = str2double(get(handles.edit_Fz,'String'));
if isnan(fz)
    fz = 0;
end
mx = str2double(get(handles.edit_Mx,'String'));
if isnan(mx)
    mx = 0;
end
my = str2double(get(handles.edit_My,'String'));
if isnan(my)
    my = 0;
end
mz = str2double(get(handles.edit_Mz,'String'));
if isnan(mz)
    mz = 0;
end

if (isnan(fx)) || (isnan(fy)) || (isnan(fz)) || (isnan(mx)) || (isnan(my)) || (isnan(mz))
    msgbox('Invalid input data!', 'Error','error');
    return
end

for n = n_ID
    % Set load cases of Node object
    if isempty(nodes(n).nodalLoadCase) == 1
        nodes(n).nodalLoadCase = zeros(12,lc);
    end    
    nodes(n).nodalLoadCase(1,lc) = fx;
    nodes(n).nodalLoadCase(2,lc) = fy;
    nodes(n).nodalLoadCase(3,lc) = fz;
    nodes(n).nodalLoadCase(4,lc) = mx;
    nodes(n).nodalLoadCase(5,lc) = my;
    nodes(n).nodalLoadCase(6,lc) = mz;

    % Set load property of Node object
    if (fx ~= 0) || (fy ~= 0) || (fz ~= 0) || (mx ~= 0) || (my ~= 0) || (mz ~= 0)
        nodes(n).load.static(1) = fx;
        nodes(n).load.static(2) = fy;
        nodes(n).load.static(3) = fz;
        nodes(n).load.static(4) = mx;
        nodes(n).load.static(5) = my;
        nodes(n).load.static(6) = mz;
    else
        nodes(n).load.static = [];
    end
end

% Set model object properties
model.nodes = nodes;

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
setappdata(0,'nodes',nodes);
setappdata(0,'model',model);

% Make GUI a normal window
gui = findobj('Tag','GUI_NodalLoads');
set(gui,'WindowStyle','normal');

% Draw updated model
if allLoadsNeedToBeDrawn == false
    redraw(mdata,'Nodal Loads')
else
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');

% Enable button for futre use
set(hObject,'enable','on')

% Relaunch GUI
setappdata(0,'move',0);
GUI_NodalLoads

%--------------------------------------------------------------------------
% Executes during edit_Fx creation, after setting all properties.
function edit_Fx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        fx = num2str(nodes(1).nodalLoadCase(1,lc));
    else
        fx = num2str(0);
    end
    set(hObject,'Enable','on','String',fx)
end

function edit_Fx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Fy creation, after setting all properties.
function edit_Fy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        fy = num2str(nodes(1).nodalLoadCase(2,lc));
    else
        fy = num2str(0);
    end
    set(hObject,'Enable','on','String',fy)
end

function edit_Fy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Fz creation, after setting all properties.
function edit_Fz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        fz = num2str(nodes(1).nodalLoadCase(3,lc));
    else
        fz = num2str(0);
    end
    set(hObject,'Enable','on','String',fz)
end

function edit_Fz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Mx creation, after setting all properties.
function edit_Mx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        mx = num2str(nodes(1).nodalLoadCase(4,lc));
    else
        mx = num2str(0);
    end
    set(hObject,'Enable','on','String',mx)
end

function edit_Mx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_My creation, after setting all properties.
function edit_My_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        my = num2str(nodes(1).nodalLoadCase(5,lc));
    else
        my = num2str(0);
    end
    set(hObject,'Enable','on','String',my)
end

function edit_My_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Mz creation, after setting all properties.
function edit_Mz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');
nodes = getappdata(0,'nodes');

if (anm == 1) || (anm == 3) || (anm == 4)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(nodes(1).nodalLoadCase,2) && isempty(nodes(1).nodalLoadCase) == 0
        mz = num2str(nodes(1).nodalLoadCase(6,lc));
    else
        mz = num2str(0);
    end
    set(hObject,'Enable','on','String',mz)
end

function edit_Mz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
function checkbox_MultiNodes_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.popupmenu_Node,'enable','off','string',' ','value',1)
    set(handles.edit_MultiNodes,'enable','on','string','1')
else
    set(handles.popupmenu_Node,'enable','on','string',num2str(1:getappdata(0,'nnp'),'%d\n'),'value',1)
    popupmenu_Node_Callback(handles.popupmenu_Node, eventdata, handles)
    set(handles.edit_MultiNodes,'enable','off','string','')
end

%--------------------------------------------------------------------------
% Auxiliary function
% Reads string to get numbers or vectors
% Input:
% * str -> string to be read
% * max -> maximum value to be read (pre-dimensions output)
% Output:
% * flag -> flag for erros (0 = error; 1 = success)
% * output -> vector of integer numbers read from string
function [flag,output] = readStr(str,max)
output = zeros(1,max);
count = 0;           % counter for output index
numFlag = false;     % flag for number being read
vctrFlag = false;    % flag for vector being read
errorFlag = false;   % flag for errors on string input

for aux = 1:size(str,2)
    if strcmp(str(aux),' ')
        numFlag = false;
    elseif ~numFlag
        if ~isnan(str2double(str(aux)))
            numFlag = true;
            count = count + 1;
            output(count) = str2double(str(aux));
            if vctrFlag && aux == size(str,2)
                vctr = linspace(output(count-1),output(count),abs(output(count)-output(count-1))+1);
                if ~all(vctr <= max)
                    errorFlag = true;
                    break
                end
                output(count-1:count-2+size(vctr,2)) = vctr;
            end
        else
            errorFlag = true;
            break
        end
    elseif numFlag
        if ~isnan(str2double(str(aux)))
            numFlag = true;
            output(count) = output(count)*10 + str2double(str(aux));
            if vctrFlag && aux == size(str,2)
                vctr = linspace(output(count-1),output(count),abs(output(count)-output(count-1))+1);
                if ~all(vctr <= max)
                    errorFlag = true;
                    break
                end
                output(count-1:count-2+size(vctr,2)) = vctr;
            end
        elseif strcmp(str(aux),';')
            numFlag = false;
            if vctrFlag
                vctr = linspace(output(count-1),output(count),abs(output(count)-output(count-1))+1);
                if ~all(vctr <= max)
                    errorFlag = true;
                    break
                end
                output(count-1:count-2+size(vctr,2)) = vctr;
                count = count-1+size(vctr,2);
                vctrFlag = false;
            end
        elseif strcmp(str(aux),'-')
            if vctrFlag || aux == size(str,2)
                errorFlag = true;
                break
            end
            numFlag = false;
            vctrFlag = true;
        else
            errorFlag = true;
            break
        end
    end
end
output = nonzeros(output)';

if errorFlag || ~all(output <= max) || isempty(output)
    flag = false;
    output = [];
else
    flag = true;
end
