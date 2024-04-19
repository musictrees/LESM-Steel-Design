%% Element Load Dialog Callback Functions
% This file contains the callback functions associated with the "Element Loads"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_ElementLoads(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ElementLoads_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ElementLoads_OutputFcn, ...
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
% Executes just before ElementLoads GUI is made visible.
% Sets GUI initial properties.
function GUI_ElementLoads_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_ElementLoads
handles.output = hObject;

% Move GUI to the center of the screen
if getappdata(0,'move') == 1
    movegui(hObject,'center')
end

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_ElementLoads_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes during popupmenu_Element creation, after setting all properties.
function popupmenu_Element_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Create list of elements
nel = getappdata(0,'nel');
e = zeros(1,nel);
for i = 1:nel
    e(i) = i;
end
e = num2str(e,'%d\n');
set(hObject,'string',e)

%--------------------------------------------------------------------------
% Executes on popupmenu_Element selection change.
function popupmenu_Element_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));

% Get ID of selected load case
lc = get(mdata.popupmenu_LoadCase,'Value');

% Get ID of selected node
e = get(hObject,'Value');

% Get uniform load components
elems = getappdata(0,'elems');
if lc <= size(elems(e).load.elemLoadCase,2) && ~isempty(elems(e).load.elemLoadCase)
    unifdir = elems(e).load.elemLoadCase(1,lc) + 1;
    qx = num2str(elems(e).load.elemLoadCase(2,lc));
    qy = num2str(elems(e).load.elemLoadCase(3,lc));
    qz = num2str(elems(e).load.elemLoadCase(4,lc));
else   
    unifdir = 1;
    qx = num2str(0);
    qy = num2str(0);
    qz = num2str(0);
end

% Get linear load components
if lc <= size(elems(e).load.elemLoadCase,2) && ~isempty(elems(e).load.elemLoadCase) &&...
   size(elems(e).load.elemLoadCase,1) > 4
    lineardir = elems(e).load.elemLoadCase(5,lc) + 1;
    qx1 = num2str(elems(e).load.elemLoadCase(6,lc));
    qy1 = num2str(elems(e).load.elemLoadCase(7,lc));
    qz1 = num2str(elems(e).load.elemLoadCase(8,lc));
    qx2 = num2str(elems(e).load.elemLoadCase(9,lc));
    qy2 = num2str(elems(e).load.elemLoadCase(10,lc));
    qz2 = num2str(elems(e).load.elemLoadCase(11,lc));
else
    lineardir = 1;
    qx1 = num2str(0);
    qy1 = num2str(0);
    qz1 = num2str(0);
    qx2 = num2str(0);
    qy2 = num2str(0);
    qz2 = num2str(0);
end

% Get thermal load components
if lc <= size(elems(e).load.elemLoadCase,2) && ~isempty(elems(e).load.elemLoadCase) &&...
   size(elems(e).load.elemLoadCase,1) > 11
    dtx = num2str(elems(e).load.elemLoadCase(12,lc));
    dty = num2str(elems(e).load.elemLoadCase(13,lc));
    dtz = num2str(elems(e).load.elemLoadCase(14,lc));
else
    dtx = num2str(0);
    dty = num2str(0);
    dtz = num2str(0);
end

% Show load components of selected element
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
set(handles.popupmenu_UnifDirection,'Value',unifdir)
set(handles.popupmenu_LinearDirection,'Value',lineardir)
if (anm == 1) || (anm == 2)
    set(handles.edit_Qx,'String',qx)
    set(handles.edit_Qy,'String',qy)
    set(handles.edit_Qz,'String','')
    set(handles.edit_Qx1,'String',qx1)
    set(handles.edit_Qx2,'String',qx2)
    set(handles.edit_Qy1,'String',qy1)
    set(handles.edit_Qy2,'String',qy2)
    set(handles.edit_Qz1,'String','')
    set(handles.edit_Qz2,'String','')
    set(handles.edit_dtx,'String',dtx)
    set(handles.edit_dty,'String',dty)
    set(handles.edit_dtz,'String','')

elseif anm == 3
    set(handles.edit_Qx,'String','')
    set(handles.edit_Qy,'String','')
    set(handles.edit_Qz,'String',qz)
    set(handles.edit_Qx1,'String','')
    set(handles.edit_Qx2,'String','')
    set(handles.edit_Qy1,'String','')
    set(handles.edit_Qy2,'String','')
    set(handles.edit_Qz1,'String',qz1)
    set(handles.edit_Qz2,'String',qz2)
    set(handles.edit_dtx,'String','')
    set(handles.edit_dty,'String','')
    set(handles.edit_dtz,'String',dtz)
	
elseif (anm == 4) || (anm == 5)
    set(handles.edit_Qx,'String',qx)
    set(handles.edit_Qy,'String',qy)
    set(handles.edit_Qz,'String',qz)
    set(handles.edit_Qx1,'String',qx1)
    set(handles.edit_Qx2,'String',qx2)
    set(handles.edit_Qy1,'String',qy1)
    set(handles.edit_Qy2,'String',qy2)
    set(handles.edit_Qz1,'String',qz1)
    set(handles.edit_Qz2,'String',qz2)
    set(handles.edit_dtx,'String',dtx)
    set(handles.edit_dty,'String',dty)
    set(handles.edit_dtz,'String',dtz)
end

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Set.
% Sets nodal load properties of a Node object.
function pushbutton_Set_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
elems = getappdata(0,'elems');

% Disable button while input data is being set
set(hObject,'enable','off')

% Get ID of selected load case
lc = get(mdata.popupmenu_LoadCase,'Value');

% Get ID of selected element
if ~get(handles.checkbox_MultiElems,'value')
    e_ID = get(handles.popupmenu_Element,'Value');
else
    e_str = get(handles.edit_MultiElems,'string');
    if strcmp(e_str,'all')
        e_ID = 1:getappdata(0,'nel');
    elseif strcmp(e_str,'ALL')
        e_ID = 1:getappdata(0,'nel');
    elseif strcmp(e_str,'All')
        e_ID = 1:getappdata(0,'nel');
    else
        [flag,e_ID] = readStr(e_str,getappdata(0,'nel'));
        if ~flag
            % Enable button for futre use
            set(hObject,'enable','on')
    
            msgbox('Invalid input data!', 'Error','error');
            return
        end
    end
end

% Get uniform load component values
qx = str2double(get(handles.edit_Qx,'String'));
if isnan(qx)
  qx = 0;  
end
qy = str2double(get(handles.edit_Qy,'String'));
if isnan(qy)
  qy = 0;  
end
qz = str2double(get(handles.edit_Qz,'String'));
if isnan(qz)
  qz = 0;  
end

if (isnan(qx)) || (isnan(qy)) || (isnan(qz))
    % Enable button for futre use
    set(hObject,'enable','on')
    
    msgbox('Invalid input data!', 'Error','error');
    return
end

for e = e_ID

    % Set uniform load to load case
    elems(e).load.elemLoadCase(1,lc) = get(handles.popupmenu_UnifDirection,'Value') - 1;
    elems(e).load.elemLoadCase(2:4,lc) = [ qx ;
                                           qy ;
                                           qz ];
    % Clear previous uniform loads
    elems(e).load.uniformDir = elems(e).load.elemLoadCase(1,lc);
    elems(e).load.uniformGbl = [];
    elems(e).load.uniformLcl = [];
    
    % Set uniform load property of Elem object
    if (qx ~= 0) || (qy ~= 0) || (qz ~= 0)
        qu = [qx,qy,qz];
        elems(e).load.setUnifLoad(qu,elems(e).load.elemLoadCase(1,lc));
    end
    
    % Avoids numeric problems with new loads
    if elems(e).load.uniformDir == 0
        for i = 1:size(elems(e).load.uniformLcl,1)
            if abs(elems(e).load.uniformLcl(i)) <= 10^-10
                elems(e).load.uniformLcl(i) = 0;
            end
        end
    elseif elems(e).load.uniformDir == 1
        for i = 1:size(elems(e).load.uniformGbl,1)
            if abs(elems(e).load.uniformGbl(i)) <= 10^-10
                elems(e).load.uniformGbl(i) = 0;
            end
        end
    end
    
    % Initialize linear loads load case (if necessary)
    if size(elems(e).load.elemLoadCase,1) < 5
        elems(e).load.elemLoadCase(5:11,lc) = zeros(7,1);
    end
    
    % Get linear load component values
    % In addtion to checking if input is not a number, check if input is
    % different from current value, considering string precision. This is to
    % avoid unwanted changes in values entered via mouse modeling, where all
    % loads are set as linearly distributed.
    qx1 = str2double(get(handles.edit_Qx1,'String'));
    if isnan(qx1)
        qx1 = 0;
    elseif abs(qx1 - elems(e).load.elemLoadCase(6,lc)) <= 0.5*10^(floor(log10(abs(qx1)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qx1 = elems(e).load.elemLoadCase(6,lc);
    end
    qy1 = str2double(get(handles.edit_Qy1,'String'));
    if isnan(qy1)
        qy1 = 0;
    elseif abs(qy1 - elems(e).load.elemLoadCase(7,lc)) <= 0.5*10^(floor(log10(abs(qy1)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qy1 = elems(e).load.elemLoadCase(7,lc);
    end
    qz1 = str2double(get(handles.edit_Qz1,'String'));
    if isnan(qz1)
        qz1 = 0;
    elseif abs(qz1 - elems(e).load.elemLoadCase(8,lc)) <= 0.5*10^(floor(log10(abs(qz1)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qz1 = elems(e).load.elemLoadCase(8,lc);
    end
    qx2 = str2double(get(handles.edit_Qx2,'String'));
    if isnan(qx2)
        qx2 = 0;
    elseif abs(qx2 - elems(e).load.elemLoadCase(9,lc)) <= 0.5*10^(floor(log10(abs(qx2)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qx2 = elems(e).load.elemLoadCase(9,lc);
    end
    qy2 = str2double(get(handles.edit_Qy2,'String'));
    if isnan(qy2)
        qy2 = 0;
    elseif abs(qy2 - elems(e).load.elemLoadCase(10,lc)) <= 0.5*10^(floor(log10(abs(qy2)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qy2 = elems(e).load.elemLoadCase(10,lc);
    end
    qz2 = str2double(get(handles.edit_Qz2,'String'));
    if isnan(qz2)
        qz2 = 0;
    elseif abs(qz2 - elems(e).load.elemLoadCase(11,lc)) <= 0.5*10^(floor(log10(abs(qz2)))-4) &&...
            elems(e).load.linearDir == get(handles.popupmenu_LinearDirection,'Value') - 1
        qz2 = elems(e).load.elemLoadCase(11,lc);
    end
    
    if (isnan(qx1)) || (isnan(qy1)) || (isnan(qz1)) || (isnan(qx2)) || (isnan(qy2)) || (isnan(qz2))
        % Enable button for futre use
        set(hObject,'enable','on')
    
        msgbox('Invalid input data!', 'Error','error');
        return
    end
    
    % Set linear load to load case
    elems(e).load.elemLoadCase(5,lc) = get(handles.popupmenu_LinearDirection,'Value') - 1;
    elems(e).load.elemLoadCase(6:11,lc) = [ qx1 ;
                                            qy1 ;
                                            qz1 ;
                                            qx2 ;
                                            qy2 ;
                                            qz2 ];
    % Clear previous linear loads
    elems(e).load.linearDir = elems(e).load.elemLoadCase(5,lc);
    elems(e).load.linearGbl = [];
    elems(e).load.linearLcl = [];
    
    % Set linear load property of Elem object
    if (qx1 ~= 0) || (qy1 ~= 0) || (qz1 ~= 0) || (qx2 ~= 0) || (qy2 ~= 0) || (qz2 ~= 0)
        ql = [qx1,qy1,qz1,qx2,qy2,qz2];
        elems(e).load.setLinearLoad(ql,elems(e).load.elemLoadCase(5,lc));
    end
    
    % Avoids numeric problems with new loads
    if elems(e).load.linearDir == 0
        for i = 1:size(elems(e).load.linearLcl,1)
            if abs(elems(e).load.linearLcl(i)) <= 10^-10
                elems(e).load.linearLcl(i) = 0;
            end
        end
    elseif elems(e).load.linearDir == 1
        for i = 1:size(elems(e).load.linearGbl,1)
            if abs(elems(e).load.linearGbl(i)) <= 10^-10
                elems(e).load.linearGbl(i) = 0;
            end
        end
    end
    
    % Get thermal load values
    dtx = str2double(get(handles.edit_dtx,'String'));
    if isnan(dtx)
        dtx = 0;
    end
    dty = str2double(get(handles.edit_dty,'String'));
    if isnan(dty)
        dty = 0;
    end
    dtz = str2double(get(handles.edit_dtz,'String'));
    if isnan(dtz)
        dtz = 0;
    end
    
    if (isnan(dtz)) || (isnan(dty)) || (isnan(dtz))
        % Enable button for futre use
        set(hObject,'enable','on')
    
        msgbox('Invalid input data!', 'Error','error');
        return
    end
    
    % Set thermal load to load case
    elems(e).load.elemLoadCase(12:14,lc) = [ dtx ;
                                             dty ;
                                             dtz ];
    
    
    % Set thermal load property of Elem object
    elems(e).load.tempVar_X = dtx;
    elems(e).load.tempVar_Y = dty;
    elems(e).load.tempVar_Z = dtz;

end

% Set model object properties
model.elems = elems;

% Enable "Process Data" button in main GUI
set(mdata.pushbutton_ProcessData,'Enable','on');

% Disable model type option
set(mdata.popupmenu_Anm,'Enable','off');

% Disable result buttons
allLoadsNeedToBeRedrawn = false;
if get(mdata.popupmenu_Results,'Value') ~= 1
    allLoadsNeedToBeRedrawn = true;
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
setappdata(0,'elems',elems);
setappdata(0,'model',model);

% Make GUI a normal window
gui = findobj('Tag','GUI_ElementLoads');
set(gui,'WindowStyle','normal');

% Draw updated model
if allLoadsNeedToBeRedrawn == false
    redraw(mdata,'Element Loads')
else
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');

% Enable button for futre use
set(hObject,'enable','on')

% Relaunch GUI
setappdata(0,'move',0);
GUI_ElementLoads

%--------------------------------------------------------------------------
% Executes during popupmenu_UnifDirection creation, after setting all properties.
function popupmenu_UnifDirection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
lc = get(mdata.popupmenu_LoadCase,'Value');
elems = getappdata(0,'elems');
if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
    unifdir = elems(1).load.elemLoadCase(1,lc) + 1;
else
    unifdir = 1;
end    
set(hObject,'Value',unifdir,'enable','on')

anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'value',1,'Enable','off');
end

% Executes on selection change in popupmenu_UnifDirection.
function popupmenu_UnifDirection_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qx creation, after setting all properties.
function edit_Qx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qx = num2str(elems(1).load.elemLoadCase(2,lc));
    else
        qx = num2str(0);
    end
    set(hObject,'String',qx,'enable','on')
end

function edit_Qx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qy creation, after setting all properties.
function edit_Qy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qy = num2str(elems(1).load.elemLoadCase(3,lc));
    else
        qy = num2str(0);
    end
    set(hObject,'String',qy,'enable','on')
end

function edit_Qy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qz creation, after setting all properties.
function edit_Qz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qz = num2str(elems(1).load.elemLoadCase(4,lc));
    else
        qz = num2str(0);
    end
    set(hObject,'String',qz,'enable','on')
end

function edit_Qz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_LinearDirection creation, after setting all properties.
function popupmenu_LinearDirection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
lc = get(mdata.popupmenu_LoadCase,'Value');
elems = getappdata(0,'elems');
if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
    lineardir = elems(1).load.elemLoadCase(5,lc) + 1;
else
    lineardir = 1;
end    
set(hObject,'Value',lineardir,'enable','on')

anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'value',1,'Enable','off');
end

% Executes on selection change in popupmenu_LinearDirection.
function popupmenu_LinearDirection_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qx1 creation, after setting all properties.
function edit_Qx1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qx1 = num2str(elems(1).load.elemLoadCase(6,lc));
    else
        qx1 = num2str(0);
    end
    set(hObject,'String',qx1,'enable','on')
end

function edit_Qx1_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qx2 creation, after setting all properties.
function edit_Qx2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qx2 = num2str(elems(1).load.elemLoadCase(9,lc));
    else
        qx2 = num2str(0);
    end
    set(hObject,'String',qx2,'enable','on')
end

function edit_Qx2_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qy1 creation, after setting all properties.
function edit_Qy1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qy1 = num2str(elems(1).load.elemLoadCase(7,lc));
    else
        qy1 = num2str(0);
    end
    set(hObject,'String',qy1,'enable','on')
end

function edit_Qy1_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qy2 creation, after setting all properties.
function edit_Qy2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if anm == 3
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
       qy2 = num2str(elems(1).load.elemLoadCase(10,lc));
    else
        qy2 = num2str(0);
    end
    set(hObject,'String',qy2,'enable','on')
end

function edit_Qy2_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qz1 creation, after setting all properties.
function edit_Qz1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qz1 = num2str(elems(1).load.elemLoadCase(8,lc));
    else
        qz1 = num2str(0);
    end
    set(hObject,'String',qz1,'enable','on')
end

function edit_Qz1_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Qz2 creation, after setting all properties.
function edit_Qz2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

elems = getappdata(0,'elems');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
lc = get(mdata.popupmenu_LoadCase,'Value');

if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
else
    if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
        qz2 = num2str(elems(1).load.elemLoadCase(11,lc));
    else
        qz2 = num2str(0);
    end
    set(hObject,'String',qz2,'enable','on')
end

function edit_Qz2_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_dtx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
lc = get(mdata.popupmenu_LoadCase,'Value');
elems = getappdata(0,'elems');
if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
    dtx = elems(1).load.elemLoadCase(12,lc);
else
    dtx = 0;
end    
set(hObject,'String',dtx,'enable','on')

anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','');
end

function edit_dtx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_dty_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
lc = get(mdata.popupmenu_LoadCase,'Value');
elems = getappdata(0,'elems');
if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
    dty = elems(1).load.elemLoadCase(13,lc);
else
    dty = 0;
end
set(hObject,'String',dty,'enable','on')

anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','');
end

function edit_dty_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_dtz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
lc = get(mdata.popupmenu_LoadCase,'Value');
elems = getappdata(0,'elems');
if lc <= size(elems(1).load.elemLoadCase,2) && isempty(elems(1).load.elemLoadCase) == 0
    dtz = elems(1).load.elemLoadCase(14,lc);
else
    dtz = 0;
end
set(hObject,'String',dtz,'enable','on')

anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
end

function edit_dtz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on selection change in popupmenu_TempDirection.
function popupmenu_TempDirection_Callback(hObject, eventdata, handles)

% Executes during object creation, after setting all properties.
function popupmenu_TempDirection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
function checkbox_MultiElems_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.popupmenu_Element,'enable','off','string',' ','value',1)
    set(handles.edit_MultiElems,'enable','on','string','1')
else
    set(handles.popupmenu_Element,'enable','on','string',num2str(1:getappdata(0,'nel'),'%d\n'),'value',1)
    popupmenu_Element_Callback(handles.popupmenu_Element, eventdata, handles)
    set(handles.edit_MultiElems,'enable','off','string','')
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
