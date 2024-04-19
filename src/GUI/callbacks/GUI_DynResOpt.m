%% Main Dialog Callback Functions
% This file contains the callback functions associated with the auxliary
% dialog for dynamic results visualization options of the graphical version
% of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
function varargout = GUI_DynResOpt(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_DynResOpt_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_DynResOpt_OutputFcn, ...
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


% --- Executes just before GUI_DynResOpt is made visible.
function GUI_DynResOpt_OpeningFcn(hObject, ~, handles, varargin)
% Move GUI to the center of the screen
movegui(hObject,'center');

main_dlg = findobj('Tag','GUI_Main');
main_handles = guidata(main_dlg);
val = main_handles.dynResSpeed;
switch val
    case 0.25
        set(handles.popupmenu_Speed,'Value',1);
    case 0.5
        set(handles.popupmenu_Speed,'Value',2);
    case 1
        set(handles.popupmenu_Speed,'Value',3);
    case 2
        set(handles.popupmenu_Speed,'Value',4);
    case 4
        set(handles.popupmenu_Speed,'Value',5);
    case 8
        set(handles.popupmenu_Speed,'Value',6);
    case 16
        set(handles.popupmenu_Speed,'Value',7);
    case 32
        set(handles.popupmenu_Speed,'Value',8);
    case -1
        set(handles.popupmenu_Speed,'Value',9);
end

% Choose default command line output for GUI_DynResOpt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_DynResOpt_OutputFcn(~, ~, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(~, ~, handles) %#ok<DEFNU>
val = get(handles.popupmenu_Speed,'Value');
main_dlg = findobj('Tag','GUI_Main');
main_handles = guidata(main_dlg);
switch val
    case 1  % 0.25x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 0.25;
    case 2  % 0.5x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 0.5;
    case 3  % 1.0x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 1;
    case 4  % 2.0x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 2;
    case 5  % 4.0x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 4;
    case 6  % 8.0x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 8;
    case 7  % 16x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 16;
    case 8  % 32x
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = 32;
    case 9  % Real Time
        % Dynamic results reproduction speed
        main_handles.dynResSpeed = -1;
end
% Update handles structure
guidata(main_dlg, main_handles);

% Close this dialog
delete(findobj('Tag','GUI_DynResOpt'));

% --- Executes on selection change in popupmenu_Speed.
function popupmenu_Speed_Callback(~, ~, ~) %#ok<DEFNU>


% --- Executes during object creation, after setting all properties.
function popupmenu_Speed_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
