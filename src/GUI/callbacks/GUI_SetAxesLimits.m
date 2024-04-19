%% Set Axes Limits Dialog Callback Functions
% This file contains the callback functions associated with the Set Axes 
% Limits dialog of the graphical version of the LESM program.
% Called by Grillage models to set axes where user can interact with mouse.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_SetAxesLimits(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SetAxesLimits_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SetAxesLimits_OutputFcn, ...
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

% --- Executes just before GUI_SetAxesLimits is made visible.
function GUI_SetAxesLimits_OpeningFcn(hObject, ~, handles, varargin)
% Choose default command line output for GUI_SetAxesLimits
handles.output = hObject;

% Move GUI to the center of the screen
if getappdata(0,'move') == 1
    movegui(hObject,'center');
end

% Get handle to GUI_Main
mdata = guidata(findobj('Tag','GUI_Main'));

% Get axes limits
origXLim = get(mdata.axes_Canvas,'XLim');
origYLim = get(mdata.axes_Canvas,'YLim');

% Set current axes limits to editable text strings
set(handles.editText_XLim_From,'Enable','on','String',num2str(origXLim(1)));
set(handles.editText_XLim_To,  'Enable','on','String',num2str(origXLim(2)));
set(handles.editText_YLim_From,'Enable','on','String',num2str(origYLim(1)));
set(handles.editText_YLim_To,  'Enable','on','String',num2str(origYLim(2)));
set(handles.editText_ZLim_From,'Enable','off','String','0');
set(handles.editText_ZLim_To,  'Enable','off','String','0');

% Update handles structure
guidata(hObject, handles);
%--------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = GUI_SetAxesLimits_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(~, ~, handles) %#ok<DEFNU>
% Get handle to GUI_Main
mdata = guidata(findobj('Tag','GUI_Main'));

% Get new limits
newXLim_From = str2double(get(handles.editText_XLim_From,'String'));
newXLim_To = str2double(get(handles.editText_XLim_To,'String'));
newYLim_From = str2double(get(handles.editText_YLim_From,'String'));
newYLim_To = str2double(get(handles.editText_YLim_To,'String'));

% Get original axes limits
origXLim = get(mdata.axes_Canvas,'XLim');
origYLim = get(mdata.axes_Canvas,'YLim');

% Check for validity of input XLim
if isnan(newXLim_From) || isnan(newXLim_To)
    newXLim_From = origXLim(1);
    newXLim_To   = origXLim(2);
elseif newXLim_From == newXLim_To
    newXLim_From = origXLim(1);
    newXLim_To   = origXLim(2);
elseif newXLim_From > newXLim_To
    aux = newXLim_From;
    newXLim_From = newXLim_To;
    newXLim_To   = aux;
end

if isnan(newYLim_From) || isnan(newYLim_To)
    newYLim_From = origYLim(1);
    newYLim_To   = origYLim(2);
elseif newYLim_From == newYLim_To
    newYLim_From = origYLim(1);
    newYLim_To   = origYLim(2);
elseif newYLim_From > newYLim_To
    aux = newYLim_From;
    newYLim_From = newYLim_To;
    newYLim_To   = aux;
end

% Make GUI a normal window to get axes without error sound
gui = findobj('Tag','GUI_SetAxesLimits');
set(gui,'WindowStyle','normal');
axes(mdata.axes_Canvas)
set(gui,'WindowStyle','modal');

% Set new axes limits
axis equal
set(mdata.axes_Canvas,'XLim',[newXLim_From, newXLim_To]);
set(mdata.axes_Canvas,'YLim',[newYLim_From, newYLim_To]);
set(mdata.axes_Canvas,'ZLim',[0, ((max(xlim)-min(xlim))+(max(ylim)-min(ylim)))/1000]);
xlabel('X');
ylabel('Y');
zlabel(' ');
zticks(0);
zticklabels({' '});

% Set ruler as visible
mdata.axes_Canvas.XAxis.Visible = 'on';
mdata.axes_Canvas.YAxis.Visible = 'on';
mdata.axes_Canvas.ZAxis.Visible = 'on';
set(mdata.rulerButton,'Checked','on')

% Turn grid on
grid on
set(mdata.gridButton,'Checked','on')

% Compute limtis difference
dx = (newXLim_To - newXLim_From)/diff(origXLim);
dy = (newYLim_To - newYLim_From)/diff(origYLim);
d_max = max([dx dy]);

% Update size property of draw object
draw = getappdata(0,'draw');
draw.size = draw.size * d_max;
setappdata(0,'draw',draw)

% Reset mouse zoom and original limits properties
mouse = getappdata(0,'mouse');
mouse.currentZoom = 1;
mouse.originalXLim = get(mdata.axes_Canvas, 'XLim');
mouse.originalYLim = get(mdata.axes_Canvas, 'YLim');
mouse.originalZLim = get(mdata.axes_Canvas, 'ZLim');
setappdata(0,'mouse',mouse)

% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(~, ~, ~) %#ok<DEFNU>
delete(gcf)

%--------------------------------------------------------------------------
function editText_XLim_From_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_XLim_From_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editText_XLim_To_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_XLim_To_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editText_YLim_From_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_YLim_From_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editText_YLim_To_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_YLim_To_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editText_ZLim_From_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_ZLim_From_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editText_ZLim_To_Callback(~, ~, ~) %#ok<DEFNU>

% --- Executes during object creation, after setting all properties.
function editText_ZLim_To_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
