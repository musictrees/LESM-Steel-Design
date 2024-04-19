%% Support Dialog Callback Functions
% This file contains the callback functions associated with the "Supports"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Supports(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Supports_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Supports_OutputFcn, ...
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
% Executes just before Supports GUI is made visible.
% Sets GUI initial properties.
function GUI_Supports_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Supports
handles.output = hObject;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Get handle to GUI_Main
mdata = guidata(findobj('Tag','GUI_Main'));

% Get flag of current analysis model
anm = get(mdata.popupmenu_Anm,'Value');
anl = get(mdata.popupmenu_AnalysisType,'Value');

% Get ID of selected node
nodes = getappdata(0,'nodes');
n = get(handles.popupmenu_Node,'Value');

% Check if inclined node option is on or off
inclSupp = nodes(n).isInclinedSupp;

% Set value of inclined support check box
set(handles.checkbox_InclinedSupps,'value',inclSupp);

% Call checkbox callback
if inclSupp
    checkbox_InclinedSupps_Callback(handles.checkbox_InclinedSupps, eventdata, handles)
end

% Set direction vector values to editable box strings
if ~isempty(nodes(n).inclSuppDir)
    vx = nodes(n).inclSuppDir(1);
    if abs(vx) < 10^-6
        vx = 0;
    end
    vy = nodes(n).inclSuppDir(2);
    if abs(vy) < 10^-6
        vy = 0;
    end
    vz = nodes(n).inclSuppDir(3);
    if abs(vz) < 10^-6
        vz = 0;
    end
    set(handles.edit_InclSupp_1,'string',num2str(vx,3))
    set(handles.edit_InclSupp_2,'string',num2str(vy,3))
    set(handles.edit_InclSupp_3,'string',num2str(vz,3))
    
    if anm == 4 || anm == 5
        if ~isempty(nodes(n).inclSupp_vy)
            vyx = nodes(n).inclSupp_vy(1);
            if abs(vyx) < 10^-6
                vyx = 0;
            end
            vyy = nodes(n).inclSupp_vy(2);
            if abs(vyy) < 10^-6
                vyy = 0;
            end
            vyz = nodes(n).inclSupp_vy(3);
            if abs(vyz) < 10^-6
                vyz = 0;
            end
            set(handles.edit_InclSupp_4,'string',num2str(vyx,3),'enable','on')
            set(handles.edit_InclSupp_5,'string',num2str(vyy,3),'enable','on')
            set(handles.edit_InclSupp_6,'string',num2str(vyz,3),'enable','on')
        end
    end
else
    set(handles.edit_InclSupp_1,'string','1')
    set(handles.edit_InclSupp_2,'string','0')
    set(handles.edit_InclSupp_3,'string','0')
end

% Set selected node support conditions
if nodes(n).ebc(1) == 1
    set(handles.checkbox_Dx,'Value',1);
    set(handles.edit_Dx,'Enable','on');
    set(handles.checkbox_Kx,'Value',0);
    set(handles.edit_Kx,'Enable','off');
elseif nodes(n).ebc(1) == 2
    set(handles.checkbox_Dx,'Value',0);
    set(handles.edit_Dx,'Enable','off');
    set(handles.checkbox_Kx,'Value',1);
    set(handles.edit_Kx,'Enable','on');
else
    set(handles.checkbox_Dx,'Value',0);
    set(handles.edit_Dx,'Enable','off');
    set(handles.checkbox_Kx,'Value',0);
    set(handles.edit_Kx,'Enable','off');
end

if nodes(n).ebc(2) == 1
    set(handles.checkbox_Dy,'Value',1);
    set(handles.edit_Dy,'Enable','on');
    set(handles.checkbox_Ky,'Value',0);
    set(handles.edit_Ky,'Enable','off');
elseif nodes(n).ebc(2) == 2
    set(handles.checkbox_Dy,'Value',0);
    set(handles.edit_Dy,'Enable','off');
    set(handles.checkbox_Ky,'Value',1);
    set(handles.edit_Ky,'Enable','on');
else
    set(handles.checkbox_Dy,'Value',0);
    set(handles.edit_Dy,'Enable','off');
    set(handles.checkbox_Ky,'Value',0);
    set(handles.edit_Ky,'Enable','off');
end

if nodes(n).ebc(3) == 1
    set(handles.checkbox_Dz,'Value',1);
    set(handles.edit_Dz,'Enable','on');
    set(handles.checkbox_Kz,'Value',0);
    set(handles.edit_Kz,'Enable','off');
elseif nodes(n).ebc(3) == 2
    set(handles.checkbox_Dz,'Value',0);
    set(handles.edit_Dz,'Enable','off');
    set(handles.checkbox_Kz,'Value',1);
    set(handles.edit_Kz,'Enable','on');
else
    set(handles.checkbox_Dz,'Value',0);
    set(handles.edit_Dz,'Enable','off');
    set(handles.checkbox_Kz,'Value',0);
    set(handles.edit_Kz,'Enable','off');
end

if nodes(n).ebc(4) == 1
    set(handles.checkbox_Rx,'Value',1);
    set(handles.edit_Rx,'Enable','on');
    set(handles.checkbox_Krx,'Value',0);
    set(handles.edit_Krx,'Enable','off');
elseif nodes(n).ebc(4) == 2
    set(handles.checkbox_Rx,'Value',0);
    set(handles.edit_Rx,'Enable','off');
    set(handles.checkbox_Krx,'Value',1);
    set(handles.edit_Krx,'Enable','on');
else
    set(handles.checkbox_Rx,'Value',0);
    set(handles.edit_Rx,'Enable','off');
    set(handles.checkbox_Krx,'Value',0);
    set(handles.edit_Krx,'Enable','off');
end

if nodes(n).ebc(5) == 1
    set(handles.checkbox_Ry,'Value',1);
    set(handles.edit_Ry,'Enable','on');
    set(handles.checkbox_Kry,'Value',0);
    set(handles.edit_Kry,'Enable','off');
elseif nodes(n).ebc(5) == 2
    set(handles.checkbox_Ry,'Value',0);
    set(handles.edit_Ry,'Enable','off');
    set(handles.checkbox_Kry,'Value',1);
    set(handles.edit_Kry,'Enable','on');
else
    set(handles.checkbox_Ry,'Value',0);
    set(handles.edit_Ry,'Enable','off');
    set(handles.checkbox_Kry,'Value',0);
    set(handles.edit_Kry,'Enable','off');
end

if nodes(n).ebc(6) == 1
    set(handles.checkbox_Rz,'Value',1);
    set(handles.edit_Rz,'Enable','on');
    set(handles.checkbox_Krz,'Value',0);
    set(handles.edit_Krz,'Enable','off');
elseif nodes(n).ebc(6) == 2
    set(handles.checkbox_Rz,'Value',0);
    set(handles.edit_Rz,'Enable','off');
    set(handles.checkbox_Krz,'Value',1);
    set(handles.edit_Krz,'Enable','on');
else
    set(handles.checkbox_Rz,'Value',0);
    set(handles.edit_Rz,'Enable','off');
    set(handles.checkbox_Krz,'Value',0);
    set(handles.edit_Krz,'Enable','off');
end

% Get selected load case ID
lc = get(mdata.popupmenu_LoadCase,'Value');

% Set selected node presc. displ.
if lc > size(nodes(n).nodalLoadCase,2)
    dx = num2str(0);
    dy = num2str(0);
    dz = num2str(0);
    rx = num2str(0);
    ry = num2str(0);
    rz = num2str(0);
else    
    dx = num2str(1e3*nodes(n).nodalLoadCase(7,lc));
    dy = num2str(1e3*nodes(n).nodalLoadCase(8,lc));
    dz = num2str(1e3*nodes(n).nodalLoadCase(9,lc));
    rx = num2str(nodes(n).nodalLoadCase(10,lc));
    ry = num2str(nodes(n).nodalLoadCase(11,lc));
    rz = num2str(nodes(n).nodalLoadCase(12,lc));
end

% Set selected node spring stiffness
if isempty (nodes(n).springStiff) == 0
    kx = num2str(nodes(n).springStiff(1));
    ky = num2str(nodes(n).springStiff(2));
    kz = num2str(nodes(n).springStiff(3));
    krx = num2str(nodes(n).springStiff(4));
    kry = num2str(nodes(n).springStiff(5));
    krz = num2str(nodes(n).springStiff(6));
else
    kx = num2str(0);
    ky = num2str(0);
    kz = num2str(0);
    krx = num2str(0);
    kry = num2str(0);
    krz = num2str(0);
end

% Disable constraints not used in the current analysis model and
% set node presc. displ. and spring stiffness
if anm == 1
    set(handles.checkbox_Dz,'Enable','off');
    set(handles.checkbox_Rx,'Enable','off');
    set(handles.checkbox_Ry,'Enable','off');
    set(handles.checkbox_Rz,'Enable','off');
    
    set(handles.edit_Dx,'String',dx)
    set(handles.edit_Dy,'String',dy)
    set(handles.edit_Dz,'String','')
    set(handles.edit_Rx,'String','')
    set(handles.edit_Ry,'String','')
    set(handles.edit_Rz,'String','')
    
    set(handles.checkbox_Kz,'Enable','off');
    set(handles.checkbox_Krx,'Enable','off');
    set(handles.checkbox_Kry,'Enable','off');
    set(handles.checkbox_Krz,'Enable','off');
    
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String','')
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String','')
    
elseif anm == 2
    set(handles.checkbox_Dz,'Enable','off');
    set(handles.checkbox_Rx,'Enable','off');
    set(handles.checkbox_Ry,'Enable','off');
    
    set(handles.edit_Dx,'String',dx)
    set(handles.edit_Dy,'String',dy)
    set(handles.edit_Dz,'String','')
    set(handles.edit_Rx,'String','')
    set(handles.edit_Ry,'String','')
    set(handles.edit_Rz,'String',rz)
    
    set(handles.checkbox_Kz,'Enable','off');
    set(handles.checkbox_Krx,'Enable','off');
    set(handles.checkbox_Kry,'Enable','off');
    
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String','')
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String',krz)
    
elseif anm == 3
    set(handles.checkbox_Dx,'Enable','off');
    set(handles.checkbox_Dy,'Enable','off');
    set(handles.checkbox_Rz,'Enable','off');
    
    set(handles.edit_Dx,'String','')
    set(handles.edit_Dy,'String','')
    set(handles.edit_Dz,'String',dz)
    set(handles.edit_Rx,'String',rx)
    set(handles.edit_Ry,'String',ry)
    set(handles.edit_Rz,'String','')
    
    set(handles.checkbox_Kx,'Enable','off');
    set(handles.checkbox_Ky,'Enable','off');
    set(handles.checkbox_Krz,'Enable','off');
    
    set(handles.edit_Kx,'String','')
    set(handles.edit_Ky,'String','')
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String',krx)
    set(handles.edit_Kry,'String',kry)
    set(handles.edit_Krz,'String','')
    
    set(handles.checkbox_InclinedSupps,'enable','off') % NO INCLINED SUPPORTS
                                                       % ON GRILLAGE MODELS!
elseif anm == 4
    set(handles.checkbox_Rx,'Enable','off');
    set(handles.checkbox_Ry,'Enable','off');
    set(handles.checkbox_Rz,'Enable','off');
    
    set(handles.edit_Dx,'String',dx)
    set(handles.edit_Dy,'String',dy)
    set(handles.edit_Dz,'String',dz)
    set(handles.edit_Rx,'String','')
    set(handles.edit_Ry,'String','')
    set(handles.edit_Rz,'String','')
    
    set(handles.checkbox_Krx,'Enable','off');
    set(handles.checkbox_Kry,'Enable','off');
    set(handles.checkbox_Krz,'Enable','off');
    
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String','')
    
elseif anm == 5
    set(handles.edit_Dx,'String',dx)
    set(handles.edit_Dy,'String',dy)
    set(handles.edit_Dz,'String',dz)
    set(handles.edit_Rx,'String',rx)
    set(handles.edit_Ry,'String',ry)
    set(handles.edit_Rz,'String',rz)
    
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String',krx)
    set(handles.edit_Kry,'String',kry)
    set(handles.edit_Krz,'String',krz)
end

% Disable presc displ in dynamic analysis
if anl == 2
    set(handles.edit_Dx,'Enable','off');
    set(handles.edit_Dy,'Enable','off');
    set(handles.edit_Dz,'Enable','off');
    set(handles.edit_Rx,'Enable','off');
    set(handles.edit_Ry,'Enable','off');
    set(handles.edit_Rz,'Enable','off');
    if anm == 1
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 2
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','0')
    elseif anm == 3
        set(handles.edit_Dx,'String','')
        set(handles.edit_Dy,'String','')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','0')
        set(handles.edit_Ry,'String','0')
        set(handles.edit_Rz,'String','')
    elseif anm == 4
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 5
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','0')
        set(handles.edit_Ry,'String','0')
        set(handles.edit_Rz,'String','0')
    end
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Supports_OutputFcn(hObject, eventdata, handles)
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
% Executes on node selection change.
function popupmenu_Node_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
anl = get(mdata.popupmenu_AnalysisType,'Value');
nodes = getappdata(0,'nodes');

% Get ID of selected node
n = get(hObject,'Value');

% Set selected node support conditions
if nodes(n).ebc(1) == 1
    set(handles.checkbox_Dx,'Value',1);
    if anl == 1
        set(handles.edit_Dx,'Enable','on');
    elseif anl == 2
        set(handles.edit_Dx,'Enable','off');
    end
    set(handles.checkbox_Kx,'Value',0);
    set(handles.edit_Kx,'Enable','off');
elseif nodes(n).ebc(1) == 2
    set(handles.checkbox_Dx,'Value',0);
    set(handles.edit_Dx,'Enable','off');
    set(handles.checkbox_Kx,'Value',1);
    set(handles.edit_Kx,'Enable','on');
else
    set(handles.checkbox_Dx,'Value',0);
    set(handles.edit_Dx,'Enable','off');
    set(handles.checkbox_Kx,'Value',0);
    set(handles.edit_Kx,'Enable','off');
end

if nodes(n).ebc(2) == 1
    set(handles.checkbox_Dy,'Value',1);
    if anl == 1
        set(handles.edit_Dy,'Enable','on');
    elseif anl == 2
        set(handles.edit_Dy,'Enable','off');
    end
    set(handles.checkbox_Ky,'Value',0);
    set(handles.edit_Ky,'Enable','off');
elseif nodes(n).ebc(2) == 2
    set(handles.checkbox_Dy,'Value',0);
    set(handles.edit_Dy,'Enable','off');
    set(handles.checkbox_Ky,'Value',1);
    set(handles.edit_Ky,'Enable','on');
else
    set(handles.checkbox_Dy,'Value',0);
    set(handles.edit_Dy,'Enable','off');
    set(handles.checkbox_Ky,'Value',0);
    set(handles.edit_Ky,'Enable','off');
end

if nodes(n).ebc(3) == 1
    set(handles.checkbox_Dz,'Value',1);
    if anl == 1
        set(handles.edit_Dz,'Enable','on');
    elseif anl == 2
        set(handles.edit_Dz,'Enable','off');
    end
    set(handles.checkbox_Kz,'Value',0);
    set(handles.edit_Kz,'Enable','off');
elseif nodes(n).ebc(3) == 2
    set(handles.checkbox_Dz,'Value',0);
    set(handles.edit_Dz,'Enable','off');
    set(handles.checkbox_Kz,'Value',1);
    set(handles.edit_Kz,'Enable','on');
else
    set(handles.checkbox_Dz,'Value',0);
    set(handles.edit_Dz,'Enable','off');
    set(handles.checkbox_Kz,'Value',0);
    set(handles.edit_Kz,'Enable','off');
end

if nodes(n).ebc(4) == 1
    set(handles.checkbox_Rx,'Value',1);
    if anl == 1
        set(handles.edit_Rx,'Enable','on');
    elseif anl == 2
        set(handles.edit_Rx,'Enable','off');
    end
    set(handles.checkbox_Krx,'Value',0);
    set(handles.edit_Krx,'Enable','off');
elseif nodes(n).ebc(4) == 2
    set(handles.checkbox_Rx,'Value',0);
    set(handles.edit_Rx,'Enable','off');
    set(handles.checkbox_Krx,'Value',1);
    set(handles.edit_Krx,'Enable','on');
else
    set(handles.checkbox_Rx,'Value',0);
    set(handles.edit_Rx,'Enable','off');
    set(handles.checkbox_Krx,'Value',0);
    set(handles.edit_Krx,'Enable','off');
end

if nodes(n).ebc(5) == 1
    set(handles.checkbox_Ry,'Value',1);
    if anl == 1
        set(handles.edit_Ry,'Enable','on');
    elseif anl == 2
        set(handles.edit_Ry,'Enable','off');
    end
    set(handles.checkbox_Kry,'Value',0);
    set(handles.edit_Kry,'Enable','off');
elseif nodes(n).ebc(5) == 2
    set(handles.checkbox_Ry,'Value',0);
    set(handles.edit_Ry,'Enable','off');
    set(handles.checkbox_Kry,'Value',1);
    set(handles.edit_Kry,'Enable','on');
else
    set(handles.checkbox_Ry,'Value',0);
    set(handles.edit_Ry,'Enable','off');
    set(handles.checkbox_Kry,'Value',0);
    set(handles.edit_Kry,'Enable','off');
end

if nodes(n).ebc(6) == 1
    set(handles.checkbox_Rz,'Value',1);
    if anl == 1
        set(handles.edit_Rz,'Enable','on');
    elseif anl == 2
        set(handles.edit_Rz,'Enable','off');
    end
    set(handles.checkbox_Krz,'Value',0);
    set(handles.edit_Krz,'Enable','off');
elseif nodes(n).ebc(6) == 2
    set(handles.checkbox_Rz,'Value',0);
    set(handles.edit_Rz,'Enable','off');
    set(handles.checkbox_Krz,'Value',1);
    set(handles.edit_Krz,'Enable','on');
else
    set(handles.checkbox_Rz,'Value',0);
    set(handles.edit_Rz,'Enable','off');
    set(handles.checkbox_Krz,'Value',0);
    set(handles.edit_Krz,'Enable','off');
end

if anl == 1 % static
    % Get selected load case ID
    lc = get(mdata.popupmenu_LoadCase,'Value');

    % Set selected node presc. displ.
    if lc > size(nodes(n).nodalLoadCase,2)
        dx = num2str(0);
        dy = num2str(0);
        dz = num2str(0);
        rx = num2str(0);
        ry = num2str(0);
        rz = num2str(0);
    else
        dx = num2str(1e3*nodes(n).nodalLoadCase(7,lc));
        dy = num2str(1e3*nodes(n).nodalLoadCase(8,lc));
        dz = num2str(1e3*nodes(n).nodalLoadCase(9,lc));
        rx = num2str(nodes(n).nodalLoadCase(10,lc));
        ry = num2str(nodes(n).nodalLoadCase(11,lc));
        rz = num2str(nodes(n).nodalLoadCase(12,lc));
    end

    % Set node presc. displ.
    if anm == 1
        set(handles.edit_Dx,'String',dx)
        set(handles.edit_Dy,'String',dy)
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 2
        set(handles.edit_Dx,'String',dx)
        set(handles.edit_Dy,'String',dy)
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String',rz)
    elseif anm == 3
        set(handles.edit_Dx,'String','')
        set(handles.edit_Dy,'String','')
        set(handles.edit_Dz,'String',dz)
        set(handles.edit_Rx,'String',rx)
        set(handles.edit_Ry,'String',ry)
        set(handles.edit_Rz,'String','')
    elseif anm == 4
        set(handles.edit_Dx,'String',dx)
        set(handles.edit_Dy,'String',dy)
        set(handles.edit_Dz,'String',dz)
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 5
        set(handles.edit_Dx,'String',dx)
        set(handles.edit_Dy,'String',dy)
        set(handles.edit_Dz,'String',dz)
        set(handles.edit_Rx,'String',rx)
        set(handles.edit_Ry,'String',ry)
        set(handles.edit_Rz,'String',rz)
    end
elseif anl == 2 % dynamic
    if anm == 1
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 2
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','0')
    elseif anm == 3
        set(handles.edit_Dx,'String','')
        set(handles.edit_Dy,'String','')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','0')
        set(handles.edit_Ry,'String','0')
        set(handles.edit_Rz,'String','')
    elseif anm == 4
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','')
        set(handles.edit_Ry,'String','')
        set(handles.edit_Rz,'String','')
    elseif anm == 5
        set(handles.edit_Dx,'String','0')
        set(handles.edit_Dy,'String','0')
        set(handles.edit_Dz,'String','0')
        set(handles.edit_Rx,'String','0')
        set(handles.edit_Ry,'String','0')
        set(handles.edit_Rz,'String','0')
    end
end

% Set selected node spring stiffness
if isempty(nodes(n).springStiff) == 0
    kx = num2str(nodes(n).springStiff(1));
    ky = num2str(nodes(n).springStiff(2));
    kz = num2str(nodes(n).springStiff(3));
    krx = num2str(nodes(n).springStiff(4));
    kry = num2str(nodes(n).springStiff(5));
    krz = num2str(nodes(n).springStiff(6));
else
    kx = num2str(0);
    ky = num2str(0);
    kz = num2str(0);
    krx = num2str(0);
    kry = num2str(0);
    krz = num2str(0);
end

% Set node spring stiffness
if anm == 1
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String','')
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String','')
elseif anm == 2
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String','')
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String',krz)
elseif anm == 3
    set(handles.edit_Kx,'String','')
    set(handles.edit_Ky,'String','')
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String',krx)
    set(handles.edit_Kry,'String',kry)
    set(handles.edit_Krz,'String','')
elseif anm == 4
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String','')
    set(handles.edit_Kry,'String','')
    set(handles.edit_Krz,'String','')
elseif anm == 5
    set(handles.edit_Kx,'String',kx)
    set(handles.edit_Ky,'String',ky)
    set(handles.edit_Kz,'String',kz)
    set(handles.edit_Krx,'String',krx)
    set(handles.edit_Kry,'String',kry)
    set(handles.edit_Krz,'String',krz)
end

% Check if inclined node option is on or off
inclSupp = nodes(n).isInclinedSupp;

% Set value of inclined support check box
set(handles.checkbox_InclinedSupps,'value',inclSupp);

% Hide vy input options
if (anm == 4 || anm == 5) &&...
   (get(handles.radiobutton_DirectionVector,'value') && inclSupp ||...
    get(handles.radiobutton_Angles,'value') && ~inclSupp)
    setVy(handles,[],false);
end

% Call checkbox callback
checkbox_InclinedSupps_Callback(handles.checkbox_InclinedSupps, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on button press in "Set" pushbutton.
% Sets nodal supports and presc. displ. of a Node object.
function pushbutton_Set_Callback(hObject, eventdata, handles)
% Get handle to GUI_Main
mdata = guidata(findobj('Tag','GUI_Main'));

% Get variables from root
nodes = getappdata(0,'nodes');
model = getappdata(0,'model');

% Get analysis model id
anm = model.anm.analysis_type + 1;
        
% Disable button while input data is being set
set(hObject,'enable','off')

% Initialize dof counters
countNewFree = 0;
countNewFixed = 0;
countNewSpring = 0;
countRemovedFree = 0;
countRemovedFixed = 0;
countRemovedSpring = 0;

% Get selected load case ID 
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

% Check if new support is inclined
if get(handles.checkbox_InclinedSupps,'value')
    inclSupp = true;
    
    % Check if user entered a direction vector or angle
    if get(handles.radiobutton_DirectionVector,'value')
        % Get user input from editable text boxes
        dir_x = str2double(get(handles.edit_InclSupp_1,'string'));
        if isnan(dir_x)
            dir_x = 0;
        end
        dir_y = str2double(get(handles.edit_InclSupp_2,'string'));
        if isnan(dir_y)
            dir_y = 0;
        end
        dir_z = str2double(get(handles.edit_InclSupp_3,'string'));
        if isnan(dir_z)
            dir_z = 0;
        end
        
        % Assemble direction vector
        dir = [dir_x, dir_y, dir_z];
        
        % Check for input errors
        if all(dir == 0)
            % Enable button for futre use
            set(hObject,'enable','on')
            
            msgbox('Support local axes X cannot be a null vector. Please, change vector vx.', 'Error','error');
            return
        end
        
        % Normalize direction vector
        dir = dir / norm(dir);
        
        % Check if a vy direction versor must be collected
        if strcmp(get(handles.edit_InclSupp_4,'enable'),'on') &&...
           strcmp(get(handles.edit_InclSupp_5,'enable'),'on') &&...
           strcmp(get(handles.edit_InclSupp_6,'enable'),'on')
            % Get user input from editable text boxes
            diry_x = str2double(get(handles.edit_InclSupp_4,'string'));
            if isnan(diry_x)
                diry_x = 0;
            end
            diry_y = str2double(get(handles.edit_InclSupp_5,'string'));
            if isnan(diry_y)
                diry_y = 0;
            end
            diry_z = str2double(get(handles.edit_InclSupp_6,'string'));
            if isnan(diry_z)
                diry_z = 0;
            end

            % Assemble direction vector
            dir_y = [diry_x, diry_y, diry_z];
            
            % Check for input errors
            if all(dir_y == 0)
                % Enable button for futre use
                set(hObject,'enable','on')
                
                msgbox('Support local axes Y cannot be a null vector. Please, change vector vy.', 'Error','error');
                return
            end
            
            % Compare vectors 'dir' and 'dir_y'
            w = cross(dir,dir_y);
            if (abs(w(1)) < 1e-10) && (abs(w(2)) < 1e-10) && (abs(w(3)) < 1e-10)
                % Enable button for futre use
                set(hObject,'enable','on')
                msgbox('Support local axes X and Y are parallels. Please, change vector vx or vy.', 'Error','error');
                return
            end

            % Normalize direction vector
            dir_y = dir_y / norm(dir_y);
        else
            dir_y = [];
        end

    else % User entered angle
        % Get user input from editable text boxes
        if anm == 1 || anm == 2
            thetaX = 0;
            thetaY = 0;
            thetaZ = str2double(get(handles.edit_InclSupp_2,'string'));
            if isnan(thetaZ)
                thetaZ = 0;
            end
        elseif anm == 4 || anm == 5
            thetaX = str2double(get(handles.edit_InclSupp_1,'string'));
            if isnan(thetaX)
                thetaX = 0;
            end
            thetaY = str2double(get(handles.edit_InclSupp_2,'string'));
            if isnan(thetaY)
                thetaY = 0;
            elseif abs(thetaY) > 90
                % Enable button for futre use
                set(hObject,'enable','on')
                msgbox('Invalid input data - absolute rot. Y cannot be greater then 90°.', 'Error','error');
                return
            end
            thetaZ = str2double(get(handles.edit_InclSupp_3,'string'));
            if isnan(thetaZ)
                thetaZ = 0;
            elseif abs(thetaY) == 90 && thetaX ~= 0 && thetaZ ~= 0
                choice = questdlg('When absolute Rot. Y is 90°, Rot. X and Rot. Z are on the same direction. For that reason, only the specified value for Rot. X will be set.',...
                                  'Attention','OK','Set Rot. Z instead','Cancel','OK');
                if strcmp(choice,'Cancel')
                    % Enable button for futre use
                    set(hObject,'enable','on')
                    return
                elseif strcmp(choice,'Set Rot. Z instead')
                    thetaX = 0;
                    set(handles.edit_InclSupp_1,'string','0')
                else
                    thetaZ = 0;
                    set(handles.edit_InclSupp_3,'string','0')
                end
            end
        end
        
        % Convert from degrees to rad
        thetaX = thetaX * pi / 180;
        if abs(thetaX) < 10^-6
            thetaX = 0;
        end
        
        thetaY = thetaY * pi / 180;
        if abs(thetaY) < 10^-6
            thetaY = 0;
        end
        
        thetaZ = thetaZ * pi / 180;
        if abs(thetaZ) < 10^-6
            thetaZ = 0;
        end
        
        % Assemble direction vector vx
        dir = [cos(thetaZ)*cos(thetaY), sin(thetaZ)*cos(thetaY), sin(thetaY)];
        
        % Normalize direction vector vx
        dir = dir / norm(dir);
        
        % Compute defalut local axis y direction
        if anm == 4 || anm == 5
            % Auxiliar vz node local axis z orientation versor
            vz = [sin(thetaX)*sin(thetaZ) - cos(thetaX)*cos(thetaZ)*sin(thetaY), - sin(thetaX)*cos(thetaZ) - cos(thetaX)*sin(thetaY)*sin(thetaZ), cos(thetaX)*cos(thetaY)];
            % Compute direction vector vy
            y = cross(vz,dir);
            dir_y = y / norm(y);
        else
            dir_y = [];
        end
    end
else
    inclSupp = false;
end

% Get constraint flags
Dx = get(handles.checkbox_Dx,'Value');
Dy = get(handles.checkbox_Dy,'Value');
Dz = get(handles.checkbox_Dz,'Value');
Rx = get(handles.checkbox_Rx,'Value');
Ry = get(handles.checkbox_Ry,'Value');
Rz = get(handles.checkbox_Rz,'Value');
ebc = [Dx, Dy, Dz, Rx, Ry, Rz];

% Get presc. displ.
dx = 1e-3*str2double(get(handles.edit_Dx,'String'));
if isnan(dx)
    dx = 0;
end
dy = 1e-3*str2double(get(handles.edit_Dy,'String'));
if isnan(dy)
    dy = 0;
end
dz = 1e-3*str2double(get(handles.edit_Dz,'String'));
if isnan(dz)
    dz = 0;
end
rx = str2double(get(handles.edit_Rx,'String'));
if isnan(rx)
    rx = 0;
end
ry = str2double(get(handles.edit_Ry,'String'));
if isnan(ry)
    ry = 0;
end
rz = str2double(get(handles.edit_Rz,'String'));
if isnan(rz)
    rz = 0;
end

if (isnan(dx)) || (isnan(dy)) || (isnan(dz)) || (isnan(rx)) || (isnan(ry)) || (isnan(rz))
    % Enable button for futre use
    set(hObject,'enable','on')
    msgbox('Invalid input data!', 'Error','error');
    return
end

if (dx ~= 0) || (dy ~= 0) || (dz ~= 0) || (rx ~= 0) || (ry ~= 0) || (rz ~= 0)
    prescDispl = [dx, dy, dz, rx, ry, rz];
else
    prescDispl = [];
end

% Get spring flags
springFlag = [get(handles.checkbox_Kx,'Value')
              get(handles.checkbox_Ky,'Value')
              get(handles.checkbox_Kz,'Value')
              get(handles.checkbox_Krx,'Value')
              get(handles.checkbox_Kry,'Value')
              get(handles.checkbox_Krz,'Value')];  


% Get spring stiffness
kx = str2double(get(handles.edit_Kx,'String'));
if isnan(kx)
    kx = 0;    
end
ky = str2double(get(handles.edit_Ky,'String'));
if isnan(ky)
    ky = 0;
end
kz = str2double(get(handles.edit_Kz,'String'));
if isnan(kz)
    kz = 0;
end
krx = str2double(get(handles.edit_Krx,'String'));
if isnan(krx)
    krx = 0;
end
kry = str2double(get(handles.edit_Kry,'String'));
if isnan(kry)
    kry = 0;
end
krz = str2double(get(handles.edit_Krz,'String'));
if isnan(krz)
    krz = 0;
end

if (kx < 0) || (ky < 0) || (kz < 0) || (krx < 0) || (kry < 0) || (krz <0)
    % Enable button for futre use
    set(hObject,'enable','on')
    
    msgbox('Invalid input data! Spring stiffness must be a non-negative.','Error','error');
    return
end  

% Assemble spring flags in ebc vector
if (springFlag(1) == 1) && (kx > 0)
    ebc(1) = 2;
end    
if (springFlag(2) == 1) && (ky > 0)
    ebc(2) = 2;
end
if (springFlag(3) == 1) && (kz > 0)
    ebc(3) = 2;
end  
if (springFlag(4) == 1) && (krx > 0)
    ebc(4) = 2;
end  
if (springFlag(5) == 1) && (kry > 0)
    ebc(5) = 2;
end 
if (springFlag(6) == 1) && (krz > 0)
    ebc(6) = 2;
end

if (isnan(kx)) || (isnan(ky)) || (isnan(kz)) || (isnan(krx)) || (isnan(kry)) || (isnan(krz))
    % Enable button for futre use
    set(hObject,'enable','on')
    
    msgbox('Invalid input data!', 'Error','error');
    return
end  

if (kx > 0) || (ky > 0) || (kz > 0) || (krx > 0) || (kry > 0) || (krz > 0)
    springStiff = [kx, ky, kz, krx, kry, krz];
else
    springStiff = [];
end

% Check if user restrained any displacements. If not, there can't be
% inclined supports, so the inclSupp flag is updated to false.
if all(ebc(1:3) == 0)
    inclSupp = false;
    dir = [];
    dir_y = [];
end

% Initialize flag
nodalLoadsNeedToBeRedrawn = false;

for n = n_ID
    % Set load case of Node object
    if ~isempty(prescDispl)
        if isempty(nodes(n).nodalLoadCase)
            nodes(n).nodalLoadCase = zeros(12,lc);
        end
        nodes(n).nodalLoadCase(7:12,lc) = prescDispl';
    elseif ~isempty(nodes(n).nodalLoadCase)
        nodes(n).nodalLoadCase(7:12,lc) = zeros(6,1);
    end

    % Check if new presc. displ. were set
    if isempty(nodes(n).prescDispl) && isempty(prescDispl)
    elseif (isempty(nodes(n).prescDispl) && ~isempty(prescDispl)) || (~isempty(nodes(n).prescDispl) && isempty(prescDispl))
        nodalLoadsNeedToBeRedrawn = true;
    elseif ~all((nodes(n).prescDispl == prescDispl) == true)
        nodalLoadsNeedToBeRedrawn = true;
    elseif ((inclSupp && ~nodes(n).isInclinedSupp) || (~inclSupp && nodes(n).isInclinedSupp)) && ~isempty(nodes(n).prescDispl)
       nodalLoadsNeedToBeRedrawn = true;
    elseif ~isempty(nodes(n).inclSuppDir)
        if ~all((dir == nodes(n).inclSuppDir) == true)
            nodalLoadsNeedToBeRedrawn = true;
        elseif ~isempty(dir_y)
            if ~all((dir_y == nodes(n).inclSupp_vy) == true)
                nodalLoadsNeedToBeRedrawn = true;
            end
        end
    end
    
    % If there are prescribed displacements, and a rotational support can
    % be added/removed, prescribed displacements must be redraw
    if anm == 2
        if ~isempty(prescDispl) && (prescDispl(1) ~= 0 || prescDispl(2) ~= 0)
            nodalLoadsNeedToBeRedrawn = true;
        end
    elseif anm == 3
        if ~isempty(prescDispl) && (prescDispl(3) ~= 0)
            nodalLoadsNeedToBeRedrawn = true;
        end
    elseif anm == 5
        if ~isempty(prescDispl) && (prescDispl(1) ~= 0 || prescDispl(2) ~= 0 || prescDispl(3) ~= 0 || prescDispl(4) ~= 0 || prescDispl(5) ~= 0 || prescDispl(6) ~= 0)
            nodalLoadsNeedToBeRedrawn = true;
        end
    end
    
    % Get number of supports added to node and update information panel
    for i = 1:6
        if nodes(n).ebc(i) ~= 0 && ebc(i) == 0
            countNewFree = countNewFree + 1;
        elseif nodes(n).ebc(i) ~= 1 && ebc(i) == 1
            countNewFixed = countNewFixed + 1;
        elseif nodes(n).ebc(i) ~= 2 && ebc(i) == 2
            countNewSpring = countNewSpring + 1;
        end    
        if nodes(n).ebc(i) == 0 && ebc(i) ~= 0
            countRemovedFree = countRemovedFree + 1;
        elseif nodes(n).ebc(i) == 1 && ebc(i) ~= 1
            countRemovedFixed = countRemovedFixed + 1;
        elseif nodes(n).ebc(i) == 2 && ebc(i) ~= 2
            countRemovedSpring = countRemovedSpring + 1;    
        end
    end

    % Set Node object support conditions and presc. displ.
    nodes(n).ebc = ebc;
    nodes(n).springStiff = springStiff;
    nodes(n).prescDispl = prescDispl;
    
    if inclSupp && isempty(dir_y)
        nodes(n).setInclinedSupp(dir);
    elseif inclSupp && ~isempty(dir_y)
        nodes(n).setInclinedSupp(dir,dir_y);
    elseif ~inclSupp && nodes(n).isInclinedSupp
        nodes(n).removeInclinedSupp();
    end
end

% Update uitable data (on GUI_Main)
infoPanelData = get(mdata.uitable_infoPanel,'Data');
nfreedof = cell2mat(infoPanelData(6,2)) + countNewFree - countRemovedFree;
nfixeddof = cell2mat(infoPanelData(7,2)) + countNewFixed - countRemovedFixed;
nspringdof = cell2mat(infoPanelData(8,2)) + countNewSpring - countRemovedSpring;
infoPanelData(6,:) = {'Free DOFs',nfreedof};
infoPanelData(7,:) = {'Fixed DOFs',nfixeddof};
infoPanelData(8,:) = {'Springs',nspringdof};
set(mdata.uitable_infoPanel,'Data',infoPanelData)

% Update mouse property
mouse = getappdata(0,'mouse');
if ~isempty(mouse.originalData)
    mouse.originalData = infoPanelData;
end
setappdata(0,'mouse',mouse)

% Set model object properties
model.nodes = nodes;

% Enable "Process Data" button in main GUI
mdata = guidata(findobj('Tag','GUI_Main'));
set(mdata.pushbutton_ProcessData,'Enable','on');

% Disable model type option
set(mdata.popupmenu_Anm,'Enable','off');

% Disable result buttons
allLoadsNeedToBeRedrawn = false; % Initialize flag
if get(mdata.popupmenu_Results,'value') ~= 1
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
setappdata(0,'nodes',nodes);
setappdata(0,'model',model);

% Make GUI a normal window
gui = findobj('Tag','GUI_Supports');
set(gui,'WindowStyle','normal');

% Draw updated model
redraw(mdata,'Nodes',false)
if allLoadsNeedToBeRedrawn == true
    redraw(mdata,'Loads')
elseif nodalLoadsNeedToBeRedrawn == true
    redraw(mdata,'Nodal Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');

% Enable button for futre use
set(hObject,'enable','on')

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Dx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','0');
end

function edit_Dx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Dy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','0');
end

function edit_Dy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Dz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 1 || anm == 2
    set(hObject,'Enable','off','String','0');
end

function edit_Dz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Rx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 1 || anm == 2 || anm == 4
    set(hObject,'Enable','off','String','0');
end

function edit_Rx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Ry_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 1 || anm == 2 || anm == 4
    set(hObject,'Enable','off','String','0');
end

function edit_Ry_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Rz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 1 || anm == 3 || anm == 4
    set(hObject,'Enable','off','String','0');
end

function edit_Rz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Dx.
function checkbox_Dx_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
if get(hObject,'Value') == 1
    if get(mdata.popupmenu_AnalysisType,'Value') == 1
       set(handles.edit_Dx,'Enable','on'); 
    end
    set(handles.checkbox_Kx,'Value',0);
    set(handles.edit_Kx,'Enable','off','String','0');
else
    set(handles.edit_Dx,'Enable','off','String','0');  
end

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Dy.
function checkbox_Dy_Callback(hObject, eventdata, handles)
% model = getappdata(0,'model');
% anm = model.anm.analysis_type + 1;
% clear model
mdata = guidata(findobj('Tag','GUI_Main'));
if get(hObject,'Value') == 1
    if get(mdata.popupmenu_AnalysisType,'Value') == 1
       set(handles.edit_Dy,'Enable','on');
    end
    set(handles.checkbox_Ky,'Value',0);
    set(handles.edit_Ky,'Enable','off','String','0');
else
    set(handles.edit_Dy,'Enable','off','String','0');
end
%     if (anm == 4 || anm == 5) %&& strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         %strcmp(get(handles.edit_InclSupp_4,'enable'),'off')
%         dx = str2double(get(handles.edit_InclSupp_1,'string'));
%         if isnan(dx)
%             dx = 0;
%         end
%         dy = str2double(get(handles.edit_InclSupp_2,'string'));
%         if isnan(dy)
%             dy = 0;
%         end
%         dz = str2double(get(handles.edit_InclSupp_3,'string'));
%         if isnan(dz)
%             dz = 0;
%         end
%         
%         if dx == 0 && dy == 0 && dz == 0
%             dir_1 = [1, 0];
%             dir_2 = [1, 0];
%             vx = [1, 0, 0];
%         else
%             dir_1 = [dx, dy];
%             dir_2 = [norm(dir_1), dz];
%             dir_1 = dir_1 / norm(dir_1);
%             dir_2 = dir_2 / norm(dir_2);
%             vx = [dx, dy, dz] / norm([dx, dy, dz]);
%         end
%         
%         % Compute thetaZ
%         thetaZ = acos(dir_1(1));
%         if dir_1(2) < 0
%             thetaZ = - thetaZ;
%         end
%         
%         % Compute thetaY
%         thetaY = asin(dir_2(2));
%         
%         vz = [- cos(thetaZ)*sin(thetaY), - sin(thetaY)*sin(thetaZ), cos(thetaY)]; % thetaX = 0
%         vy = cross(vz,vx);
%         vy = round(vy / norm(vy),3);
%         set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vy(1),3));
%         set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vy(2),3));
%         set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vy(3),3));

%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%             strcmp(get(handles.edit_InclSupp_1,'enable'),'off')
       
%         set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','0');
%    end
% else
%     set(handles.edit_Dy,'Enable','off','String','0');
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Ky,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_4,'Enable','off','string','');
%         set(handles.edit_InclSupp_5,'Enable','off','string','');
%         set(handles.edit_InclSupp_6,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%         ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Ky,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_1,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     end
% end

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Dz.
function checkbox_Dz_Callback(hObject, eventdata, handles)
% model = getappdata(0,'model');
% anm = model.anm.analysis_type + 1;
% clear model
mdata = guidata(findobj('Tag','GUI_Main'));
if get(hObject,'Value') == 1
    if get(mdata.popupmenu_AnalysisType,'Value') == 1
       set(handles.edit_Dz,'Enable','on');
    end
    set(handles.checkbox_Kz,'Value',0);
    set(handles.edit_Kz,'Enable','off','String','0');
else
    set(handles.edit_Dz,'Enable','off','String','0');
end
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         strcmp(get(handles.edit_InclSupp_4,'enable'),'off')
%         dx = str2double(get(handles.edit_InclSupp_1,'string'));
%         if isnan(dx)
%             dx = 0;
%         end
%         dy = str2double(get(handles.edit_InclSupp_2,'string'));
%         if isnan(dy)
%             dy = 0;
%         end
%         dz = str2double(get(handles.edit_InclSupp_3,'string'));
%         if isnan(dz)
%             dz = 0;
%         end
%         
%         if dx == 0 && dy == 0 && dz == 0
%             dir_1 = [1, 0];
%             dir_2 = [1, 0];
%             vx = [1, 0, 0];
%         else
%             dir_1 = [dx, dy];
%             dir_2 = [norm(dir_1), dz];
%             dir_1 = dir_1 / norm(dir_1);
%             dir_2 = dir_2 / norm(dir_2);
%             vx = [dx, dy, dz] / norm([dx, dy, dz]);
%         end
%         
%         % Compute thetaZ
%         thetaZ = acos(dir_1(1));
%         if dir_1(2) < 0
%             thetaZ = - thetaZ;
%         end
%         
%         % Compute thetaY
%         thetaY = asin(dir_2(2));
%         
%         vz = [- cos(thetaZ)*sin(thetaY), - sin(thetaY)*sin(thetaZ), cos(thetaY)]; % thetaX = 0
%         vy = cross(vz,vx);
%         vy = round(vy / norm(vy),3);
%         set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vy(1),3));
%         set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vy(2),3));
%         set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vy(3),3));
% 
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%             strcmp(get(handles.edit_InclSupp_1,'enable'),'off')
%        
%         set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','0');
%     end
% else
%     set(handles.edit_Dz,'Enable','off','String','0');
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Ky,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_4,'Enable','off','string','');
%         set(handles.edit_InclSupp_5,'Enable','off','string','');
%         set(handles.edit_InclSupp_6,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Ky,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_1,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     end
% end

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Rx.
function checkbox_Rx_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
anl = get(mdata.popupmenu_AnalysisType,'Value');

if anm == 3
    set(handles.checkbox_Ry,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rx,'Enable','on');
            set(handles.edit_Ry,'Enable','on');
        end
        set(handles.checkbox_Krx,'Value',0);
        set(handles.checkbox_Kry,'Value',0);
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
    else
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Ry,'Value',get(hObject,'Value'));
    set(handles.checkbox_Rz,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rx,'Enable','on');
            set(handles.edit_Ry,'Enable','on');
            set(handles.edit_Rz,'Enable','on');
        end
        set(handles.checkbox_Krx,'Value',0);
        set(handles.checkbox_Kry,'Value',0);
        set(handles.checkbox_Krz,'Value',0);
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    else
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    end
end

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Ry.
function checkbox_Ry_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
anl = get(mdata.popupmenu_AnalysisType,'Value');

if anm == 3
    set(handles.checkbox_Rx,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rx,'Enable','on');
            set(handles.edit_Ry,'Enable','on');
        end
        set(handles.checkbox_Krx,'Value',0);
        set(handles.checkbox_Kry,'Value',0);
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
    else
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Rx,'Value',get(hObject,'Value'));
    set(handles.checkbox_Rz,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rx,'Enable','on');
            set(handles.edit_Ry,'Enable','on');
            set(handles.edit_Rz,'Enable','on');
        end
        set(handles.checkbox_Krx,'Value',0);
        set(handles.checkbox_Kry,'Value',0);
        set(handles.checkbox_Krz,'Value',0);
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    else
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    end
end

%--------------------------------------------------------------------------
% Executes on button press in checkbox_Rz.
function checkbox_Rz_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
anl = get(mdata.popupmenu_AnalysisType,'Value');

if anm == 2
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rz,'Enable','on');
        end
        set(handles.checkbox_Krz,'Value',0);
        set(handles.edit_Krz,'Enable','off','String','0');
    else
        set(handles.edit_Rz,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Rx,'Value',get(hObject,'Value'));
    set(handles.checkbox_Ry,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        if anl == 1
            set(handles.edit_Rx,'Enable','on');
            set(handles.edit_Ry,'Enable','on');
            set(handles.edit_Rz,'Enable','on');
        end
        set(handles.checkbox_Krx,'Value',0);
        set(handles.checkbox_Kry,'Value',0);
        set(handles.checkbox_Krz,'Value',0);
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    else
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    end
end

%--------------------------------------------------------------------------
function edit_Kx_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_Kx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','0');
end

function edit_Ky_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_Ky_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if anm == 3
    set(hObject,'Enable','off','String','0');
end

function edit_Kz_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_Kz_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','0');
end

function edit_Krx_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

set(handles.edit_Kry,'String',get(hObject,'String'));
if anm == 5
    set(handles.edit_Krz,'String',get(hObject,'String'));
end    

% --- Executes during object creation, after setting all properties.
function edit_Krx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','0');
end

function edit_Kry_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

set(handles.edit_Krx,'String',get(hObject,'String'));
if anm == 5
    set(handles.edit_Krz,'String',get(hObject,'String'));
end    

% --- Executes during object creation, after setting all properties.
function edit_Kry_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','0');
end

function edit_Krz_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

if anm == 5
    set(handles.edit_Krx,'String',get(hObject,'String'));
    set(handles.edit_Kry,'String',get(hObject,'String'));
end    

% --- Executes during object creation, after setting all properties.
function edit_Krz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 3) || (anm == 4)
    set(hObject,'Enable','off','String','0');
end

% --- Executes on button press in checkbox_Kx.
function checkbox_Kx_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    set(handles.edit_Kx,'Enable','on');
    set(handles.checkbox_Dx,'Value',0);
    set(handles.edit_Dx,'Enable','off','String','0');
else
    set(handles.edit_Kx,'Enable','off','String','0');
end

% --- Executes on button press in checkbox_Ky.
function checkbox_Ky_Callback(hObject, eventdata, handles)
% model = getappdata(0,'model');
% anm = model.anm.analysis_type + 1;
% clear model

if get(hObject,'Value') == 1
    set(handles.edit_Ky,'Enable','on');
    set(handles.checkbox_Dy,'Value',0);
    set(handles.edit_Dy,'Enable','off','String','0');
else
    set(handles.edit_Ky,'Enable','off','String','0');
end
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         strcmp(get(handles.edit_InclSupp_4,'enable'),'off')
%         dx = str2double(get(handles.edit_InclSupp_1,'string'));
%         if isnan(dx)
%             dx = 0;
%         end
%         dy = str2double(get(handles.edit_InclSupp_2,'string'));
%         if isnan(dy)
%             dy = 0;
%         end
%         dz = str2double(get(handles.edit_InclSupp_3,'string'));
%         if isnan(dz)
%             dz = 0;
%         end
%         
%         if dx == 0 && dy == 0 && dz == 0
%             dir_1 = [1, 0];
%             dir_2 = [1, 0];
%             vx = [1, 0, 0];
%         else
%             dir_1 = [dx, dy];
%             dir_2 = [norm(dir_1), dz];
%             dir_1 = dir_1 / norm(dir_1);
%             dir_2 = dir_2 / norm(dir_2);
%             vx = [dx, dy, dz] / norm([dx, dy, dz]);
%         end
%         
%         % Compute thetaZ
%         thetaZ = acos(dir_1(1));
%         if dir_1(2) < 0
%             thetaZ = - thetaZ;
%         end
%         
%         % Compute thetaY
%         thetaY = asin(dir_2(2));
%         
%         vz = [- cos(thetaZ)*sin(thetaY), - sin(thetaY)*sin(thetaZ), cos(thetaY)]; % thetaX = 0
%         vy = cross(vz,vx);
%         vy = round(vy / norm(vy),3);
%         set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vy(1),3));
%         set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vy(2),3));
%         set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vy(3),3));
% 
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%             strcmp(get(handles.edit_InclSupp_1,'enable'),'off')
%        
%         set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','0');
%     end
% else
%     set(handles.edit_Ky,'Enable','off','String','0');
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_4,'Enable','off','string','');
%         set(handles.edit_InclSupp_5,'Enable','off','string','');
%         set(handles.edit_InclSupp_6,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Kz,'value')
%         set(handles.edit_InclSupp_1,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     end
% end

% --- Executes on button press in checkbox_Kz.
function checkbox_Kz_Callback(hObject, eventdata, handles)
% model = getappdata(0,'model');
% anm = model.anm.analysis_type + 1;
% clear model

if get(hObject,'Value') == 1
    set(handles.edit_Kz,'Enable','on');
    set(handles.checkbox_Dz,'Value',0);
    set(handles.edit_Dz,'Enable','off','String','0');
else
    set(handles.edit_Kz,'Enable','off','String','0');
end
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         strcmp(get(handles.edit_InclSupp_4,'enable'),'off')
%         dx = str2double(get(handles.edit_InclSupp_1,'string'));
%         if isnan(dx)
%             dx = 0;
%         end
%         dy = str2double(get(handles.edit_InclSupp_2,'string'));
%         if isnan(dy)
%             dy = 0;
%         end
%         dz = str2double(get(handles.edit_InclSupp_3,'string'));
%         if isnan(dz)
%             dz = 0;
%         end
%         
%         if dx == 0 && dy == 0 && dz == 0
%             dir_1 = [1, 0];
%             dir_2 = [1, 0];
%             vx = [1, 0, 0];
%         else
%             dir_1 = [dx, dy];
%             dir_2 = [norm(dir_1), dz];
%             dir_1 = dir_1 / norm(dir_1);
%             dir_2 = dir_2 / norm(dir_2);
%             vx = [dx, dy, dz] / norm([dx, dy, dz]);
%         end
%         
%         % Compute thetaZ
%         thetaZ = acos(dir_1(1));
%         if dir_1(2) < 0
%             thetaZ = - thetaZ;
%         end
%         
%         % Compute thetaY
%         thetaY = asin(dir_2(2));
%         
%         vz = [- cos(thetaZ)*sin(thetaY), - sin(thetaY)*sin(thetaZ), cos(thetaY)]; % thetaX = 0
%         vy = cross(vz,vx);
%         vy = round(vy / norm(vy),3);
%         set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vy(1),3));
%         set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vy(2),3));
%         set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vy(3),3));
% 
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%             strcmp(get(handles.edit_InclSupp_1,'enable'),'off')
%        
%         set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','0');
%     end
% else
%     set(handles.edit_Kz,'Enable','off','String','0');
%     if (anm == 4 || anm == 5) && strcmp(get(handles.text_InclSupp_4,'Visible'),'on') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Ky,'value')
%         set(handles.edit_InclSupp_4,'Enable','off','string','');
%         set(handles.edit_InclSupp_5,'Enable','off','string','');
%         set(handles.edit_InclSupp_6,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     elseif (anm == 4 || anm == 5) && get(handles.radiobutton_Angles,'value') &&...
%         ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Ky,'value')
%         set(handles.edit_InclSupp_1,'Enable','off','string','');
%         if get(handles.checkbox_AxesInclSupp3D,'value')
%             updateInclSupp3DCnvs(handles,true);
%         end
%     end
% end

% --- Executes on button press in checkbox_Krx.
function checkbox_Krx_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

if anm == 3
    set(handles.checkbox_Kry,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        set(handles.edit_Krx,'Enable','on');
        set(handles.edit_Kry,'Enable','on');
        set(handles.checkbox_Rx,'Value',0);
        set(handles.checkbox_Ry,'Value',0);
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
    else
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Kry,'Value',get(hObject,'Value'));
    set(handles.checkbox_Krz,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        set(handles.edit_Krx,'Enable','on');
        set(handles.edit_Kry,'Enable','on');
        set(handles.edit_Krz,'Enable','on');
        set(handles.checkbox_Rx,'Value',0);
        set(handles.checkbox_Ry,'Value',0);
        set(handles.checkbox_Rz,'Value',0);
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    else
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    end
end

% --- Executes on button press in checkbox_Kry.
function checkbox_Kry_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

if anm == 3
    set(handles.checkbox_Krx,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        set(handles.edit_Krx,'Enable','on');
        set(handles.edit_Kry,'Enable','on');
        set(handles.checkbox_Rx,'Value',0);
        set(handles.checkbox_Ry,'Value',0);
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
    else
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Krx,'Value',get(hObject,'Value'));
    set(handles.checkbox_Krz,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        set(handles.edit_Krx,'Enable','on');
        set(handles.edit_Kry,'Enable','on');
        set(handles.edit_Krz,'Enable','on');
        set(handles.checkbox_Rx,'Value',0);
        set(handles.checkbox_Ry,'Value',0);
        set(handles.checkbox_Rz,'Value',0);
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    else
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    end
end

% --- Executes on button press in checkbox_Krz.
function checkbox_Krz_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');

if anm == 2
    if get(hObject,'Value') == 1
        set(handles.edit_Krz,'Enable','on');
        set(handles.checkbox_Rz,'Value',0);
        set(handles.edit_Rz,'Enable','off','String','0');
    else
        set(handles.edit_Krz,'Enable','off','String','0');
    end
elseif anm == 5
    set(handles.checkbox_Krx,'Value',get(hObject,'Value'));
    set(handles.checkbox_Kry,'Value',get(hObject,'Value'));
    if get(hObject,'Value') == 1
        set(handles.edit_Krx,'Enable','on');
        set(handles.edit_Kry,'Enable','on');
        set(handles.edit_Krz,'Enable','on');
        set(handles.checkbox_Rx,'Value',0);
        set(handles.checkbox_Ry,'Value',0);
        set(handles.checkbox_Rz,'Value',0);
        set(handles.edit_Rx,'Enable','off','String','0');
        set(handles.edit_Ry,'Enable','off','String','0');
        set(handles.edit_Rz,'Enable','off','String','0');
    else
        set(handles.edit_Krx,'Enable','off','String','0');
        set(handles.edit_Kry,'Enable','off','String','0');
        set(handles.edit_Krz,'Enable','off','String','0');
    end
end

%--------------------------------------------------------------------------
function checkbox_MultiNodes_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.popupmenu_Node,'enable','off','string',' ','value',1)
    set(handles.edit_MultiNodes,'enable','on','string','1')
    if get(handles.checkbox_AxesInclSupp3D,'value')
        updateInclSupp3DCnvs(handles,true);
    end
else
    set(handles.popupmenu_Node,'enable','on','string',num2str(1:getappdata(0,'nnp'),'%d\n'),'value',1)
    popupmenu_Node_Callback(handles.popupmenu_Node, eventdata, handles)
    set(handles.edit_MultiNodes,'enable','off','string','')
end

%--------------------------------------------------------------------------
function checkbox_InclinedSupps_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    % Enable radio buttons
    set(handles.radiobutton_DirectionVector,'Enable','on','value',1);
    set(handles.radiobutton_Angles,'Enable','on','value',0);
    
    % Get current node
    nodes = getappdata(0,'nodes');
    node = nodes(get(handles.popupmenu_Node,'value'));
    
    % Get analysis model id
    model = getappdata(0,'model');
    anm = model.anm.analysis_type + 1;
    clear model
    
    % Check if current node has inclined support
    if node.isInclinedSupp && ~get(handles.checkbox_MultiNodes,'value')
        % Enable edtable text boxes
        set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string',num2str(round(node.inclSuppDir(1),3),3));
        set(handles.edit_InclSupp_2,'Enable','on','Visible','on','string',num2str(round(node.inclSuppDir(2),3),3));
        if anm == 1 || anm == 2
            set(handles.edit_InclSupp_3,'Enable','off','Visible','on','string','0');
        elseif anm == 4 || anm == 5
            set(handles.edit_InclSupp_3,'Enable','on','Visible','on','string',num2str(round(node.inclSuppDir(3),3),3));
            setVy(handles,node,true);
        end
        set(handles.text_InclSupp_1,'Visible','on','string','vx_X');
        set(handles.text_InclSupp_2,'Visible','on','string','vx_Y');
        set(handles.text_InclSupp_3,'Visible','on','string','vx_Z');
        set(handles.text_Rot_1,'Visible','off');
        set(handles.text_Rot_2,'Visible','off');
        set(handles.text_Rot_3,'Visible','off');
    else
        % Enable edtable text boxes
        set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','1');
        set(handles.edit_InclSupp_2,'Enable','on','Visible','on','string','0');
        if anm == 1 || anm == 2
            set(handles.edit_InclSupp_3,'Enable','off','Visible','on','string','0');
        elseif anm == 4 || anm == 5
            set(handles.edit_InclSupp_3,'Enable','on','Visible','on','string','0');
            setVy(handles,[],true);
        end
        set(handles.text_InclSupp_1,'Visible','on','string','vx_X');
        set(handles.text_InclSupp_2,'Visible','on','string','vx_Y');
        set(handles.text_InclSupp_3,'Visible','on','string','vx_Z');
        set(handles.text_Rot_1,'Visible','off');
        set(handles.text_Rot_2,'Visible','off');
        set(handles.text_Rot_3,'Visible','off');
    end
    
    if ~checkVyPos(handles) && (anm == 4 || anm == 5)
        changeVyPos(handles)
    elseif checkVyPos(handles) && (anm == 1 || anm == 2)
        changeVyPos(handles)
    end
    
    if anm == 4 || anm == 5
        set(handles.checkbox_AxesInclSupp3D,'enable','on','visible','on')
        if get(handles.checkbox_AxesInclSupp3D,'value')
            checkbox_AxesInclSupp3D_Callback(handles.checkbox_AxesInclSupp3D, eventdata, handles,false)
        end
    end

else
    % Get analysis model id
    model = getappdata(0,'model');
    anm = model.anm.analysis_type + 1;
    clear model
    
    % Disable radio buttons
    set(handles.radiobutton_DirectionVector,'Enable','off','value',0);
    set(handles.radiobutton_Angles,'Enable','off','value',0);
    
    % Disable edtable text boxes
    set(handles.edit_InclSupp_1,'Enable','off','Visible','off');
    set(handles.edit_InclSupp_2,'Enable','off','Visible','off');
    set(handles.edit_InclSupp_3,'Enable','off','Visible','off');
    set(handles.text_InclSupp_1,'Visible','off');
    set(handles.text_InclSupp_2,'Visible','off');
    set(handles.text_InclSupp_3,'Visible','off');
    set(handles.text_Rot_1,'Visible','off');
    set(handles.text_Rot_2,'Visible','off');
    set(handles.text_Rot_3,'Visible','off');
    if anm == 4 || anm == 5
        setVy(handles,[],false);
    end
    
    if checkVyPos(handles)
        changeVyPos(handles)
    end
    
    if anm == 4 || anm == 5
        if get(handles.checkbox_AxesInclSupp3D,'value')
            flag = true;
        else
            flag = false;
        end
        set(handles.checkbox_AxesInclSupp3D,'enable','off','visible','off','value',0)
        checkbox_AxesInclSupp3D_Callback(handles.checkbox_AxesInclSupp3D, eventdata, handles,flag)
    end
end

%--------------------------------------------------------------------------
function checkbox_AxesInclSupp3D_Callback(hObject, eventdata, handles,szChngFlag)
% Check if window size must be changed
if nargin < 4
    szChngFlag = true;
end

% Disable checkbox while window is being changed (avoid user errors)
set(hObject,'enable','off')

% Get handle to GUI_Supports figure
gui = findobj('Tag','GUI_Supports');

% Show / hide inclined supports canvas, by changing window size
%-------------------------------------------------------------
if szChngFlag
    % Adjust GUI_Supports size
    pos = get(gui,'Position');
    if get(hObject,'value')
        addedSpace = 300;
    else
        addedSpace = -300;
    end
    pos = [pos(1), pos(2), pos(3)+addedSpace, pos(4)];
    set(gui,'Position',pos);
end
%-------------------------------------------------------------

% Adjust inclined supports canvas and its drawings
%-------------------------------------------------------------
if get(hObject,'value')
    % Reset inclined supports canvas
    axes(handles.axes_InclSupp3D);
    cla reset
    
    % Make inclined supports canvas visible
    set(handles.axes_InclSupp3D,'Visible','on')
    
    % Show static text with instructions on how to use canvas
    set(handles.text_HowToUse_3DInclSuppAxis,'Visible','on')
    
    % Set inclined supports canvas view to 3D
    view(3);
    
    % Set inclined supports canvas properties
    axis equal
    grid on
    hold on
    set(handles.axes_InclSupp3D,'Clipping','off')
    
    % Set inclined supports canvas axis limits
    xlim([-1,1])
    ylim([-1,1])
    zlim([-1,1])
    
    % Hide inclined supports canvas axis rulers
    handles.axes_InclSupp3D.XAxis.Visible = 'off';
    handles.axes_InclSupp3D.YAxis.Visible = 'off';
    handles.axes_InclSupp3D.ZAxis.Visible = 'off';
    
    % Create Emouse object to handle inclined supports canvas
    mouse_inclSupp = Emouse_3D_InclSupp(gui,handles.axes_InclSupp3D);
    set(handles.axes_InclSupp3D,'UserData',mouse_inclSupp)
    
    % Update inclined supports canvas
    updateInclSupp3DCnvs(handles,true)
else
    % Clean inclined supports canvas
    delete(findobj('tag','inclSupp_gblAx'))
    delete(findobj('tag','inclSupp_lclAx'))

    % Make inclined supports canvas invisible and clean UserData
    set(handles.axes_InclSupp3D,'Visible','off','UserData',[])
    
    % Hide static text with instructions on how to use canvas
    set(handles.text_HowToUse_3DInclSuppAxis,'Visible','off')
end

% Enable checkbox for future use
set(hObject,'enable','on')

%--------------------------------------------------------------------------
function updateInclSupp3DCnvs(handles,getDirFromTextBxs)
if nargin < 2
    getDirFromTextBxs = false;
end

% Initialize auxiliary flag
thetaY_isGT90 = false;

% Make sure objects are drawn on the inclined supports canvas
axes(handles.axes_InclSupp3D);

% Clean inclined supports canvas
delete(findobj('tag','inclSupp_gblAx'))
delete(findobj('tag','inclSupp_lclAx'))

% Get handle to Draw object
drawIncl = getappdata(0,'draw');

% Draw global axis symbol
clr = [1 0 0;
       0 0.8 0;
       0 0 1];
drawIncl.axis3D(drawIncl,-0.8,-0.8,-1,0.5,clr,[],'inclSupp_gblAx');

% Get local axis current direction
if ~getDirFromTextBxs
    if strcmp(get(handles.popupmenu_Node,'enable'),'on')
        nodes = getappdata(0,'nodes');
        node = nodes(get(handles.popupmenu_Node,'value'));
        if node.isInclinedSupp
            [dx,dy,dz] = node.getInclinedSuppLocAxis;
            dir = [dx;dy;dz];
        else
            dir = eye(3);
        end
    else
        dir = eye(3);
    end
elseif get(handles.radiobutton_DirectionVector,'value')
    vxx = str2double(get(handles.edit_InclSupp_1,'string'));
    vxy = str2double(get(handles.edit_InclSupp_2,'string'));
    vxz = str2double(get(handles.edit_InclSupp_3,'string'));
    if isnan(vxx) || isnan(vxy) || isnan(vxz)
        dx = [1,0,0];
        set(handles.edit_InclSupp_1,'string','1')
        set(handles.edit_InclSupp_2,'string','0')
        set(handles.edit_InclSupp_3,'string','0')
    elseif vxx == 0 && vxy == 0 && vxz == 0
        dx = [1,0,0];
        set(handles.edit_InclSupp_1,'string','1')
        set(handles.edit_InclSupp_2,'string','0')
        set(handles.edit_InclSupp_3,'string','0')
    else
        dx = [vxx,vxy,vxz]/norm([vxx,vxy,vxz]);
    end
    if strcmp(get(handles.edit_InclSupp_4,'enable'),'on')
        vyx = str2double(get(handles.edit_InclSupp_4,'string'));
        vyy = str2double(get(handles.edit_InclSupp_5,'string'));
        vyz = str2double(get(handles.edit_InclSupp_6,'string'));
        if isnan(vyx) || isnan(vyy) || isnan(vyz)
            vy = [0,1,0];
            set(handles.edit_InclSupp_4,'string','0')
            set(handles.edit_InclSupp_5,'string','1')
            set(handles.edit_InclSupp_6,'string','0')
        elseif vyx == 0 && vyy == 0 && vyz == 0
            vy = [0,1,0];
            set(handles.edit_InclSupp_4,'string','0')
            set(handles.edit_InclSupp_5,'string','1')
            set(handles.edit_InclSupp_6,'string','0')
        else
            vy = [vyx,vyy,vyz]/norm([vyx,vyy,vyz]);
        end
        dz = cross(dx,vy);
        if (abs(dz(1)) < 1e-10) && (abs(dz(2)) < 1e-10) && (abs(dz(3)) < 1e-10)
            dir = eye(3);
            set(handles.edit_InclSupp_4,'string','0')
            set(handles.edit_InclSupp_5,'string','1')
            set(handles.edit_InclSupp_6,'string','0')
        else
            dz = dz / norm(dz);
            dy = cross(dz,dx);
            dir = [dx;dy;dz];
        end
    else
        if abs(dx(1)) <= (10^-10) && abs(dx(2)) <= (10^-10) && abs(dx(3)) > (10^-10)
            vz = [-(dx(3)/abs(dx(3))), 0, 0];
        else
            vz = [0, 0, 1];
        end
        dy = cross(vz,dx);
        dy = dy / norm(dy);
        dz = cross(dx,dy);
        dir = [dx;dy;dz];
    end
else % angles
    thetaY_is90 = false;   % flag
    
    thetaX = str2double(get(handles.edit_InclSupp_1,'string'));
    if isnan(thetaX)
        thetaX = 0;
    elseif abs(thetaX) < 10^-6
        thetaX = 0;
    else
        thetaX = thetaX * pi / 180;
    end
    thetaY = str2double(get(handles.edit_InclSupp_2,'string'));
    if isnan(thetaY)
        thetaY = 0;
    elseif abs(thetaY) < 10^-6
        thetaY = 0;
    else
        if abs(abs(thetaY) - 90) < 10^-6
            thetaY_is90 = true;
        elseif abs(thetaY) > 90
            thetaY_isGT90 = true;
            thetaY = 0;
            set(handles.edit_InclSupp_2,'string','0')
        end
        thetaY = thetaY * pi / 180;
    end
    thetaZ = str2double(get(handles.edit_InclSupp_3,'string'));
    if isnan(thetaZ)
        thetaZ = 0;
    elseif abs(thetaZ) < 10^-6
        thetaZ = 0;
    else
        thetaZ = thetaZ * pi / 180;
    end
    if thetaY_is90 && thetaX ~= 0 && thetaZ ~= 0
        choice = questdlg('When absolute Rot. Y is 90°, Rot. X and Rot. Z are on the same direction. For that reason, only the specified value for Rot. X will be set.',...
                          'Attention','OK','Set Rot. Z instead','OK');
        if strcmp(choice,'Set Rot. Z instead')
          thetaX = 0;
          set(handles.edit_InclSupp_1,'string','0')
        else
          thetaZ = 0;
          set(handles.edit_InclSupp_3,'string','0')
        end
    end
        
    dx = [cos(thetaZ)*cos(thetaY), sin(thetaZ)*cos(thetaY), sin(thetaY)];
    dz = [sin(thetaX)*sin(thetaZ) - cos(thetaX)*cos(thetaZ)*sin(thetaY), - sin(thetaX)*cos(thetaZ) - cos(thetaX)*sin(thetaY)*sin(thetaZ), cos(thetaX)*cos(thetaY)];
    dy = cross(dz,dx);
    dy = dy / norm(dy);
    dir = [dx;dy;dz];
end

% Draw local axis symbol
drawIncl.axis3D(drawIncl,0,0,0,1,clr,dir,'inclSupp_lclAx',false);

% Draw local inclined support plane
drawIncl.plane3D(0,0,0,2,[0.75,0.75,0.75],dir,'inclSupp_lclAx');

% Display warning message regarding Rot Y, if input is greater than 90°
% This happens after all drawings are made, so that there is no
% misinformation as to which is the current canvas
if thetaY_isGT90
    msgbox('Absolute rot. Y cannot be greater than 90°.', 'Attention','warn');
end

%--------------------------------------------------------------------------
function radiobutton_DirectionVector_Callback(hObject, eventdata, handles)
if ~get(hObject,'value')
    if ~get(handles.radiobutton_Angles,'value')
        set(hObject,'value',1)
        return
    end
else
    % Unselect angles radio button
    set(handles.radiobutton_Angles,'Enable','on','value',0);
    
    % Get analysis model id
    model = getappdata(0,'model');
    anm = model.anm.analysis_type + 1;
    
    clear model

    % Get user input from editable text boxes
    if anm == 1 || anm == 2
        thetaX = 0;
        thetaY = 0;
        thetaZ = str2double(get(handles.edit_InclSupp_2,'string'));
        if isnan(thetaZ)
            thetaZ = 0;
        end
    elseif anm == 4 || anm == 5
        thetaX = str2double(get(handles.edit_InclSupp_1,'string'));
        if isnan(thetaX)
            thetaX = 0;
        end
        thetaY = str2double(get(handles.edit_InclSupp_2,'string'));
        if isnan(thetaY)
            thetaY = 0;
        elseif abs(thetaY) > 90
            thetaY = 0;
        end
        thetaZ = str2double(get(handles.edit_InclSupp_3,'string'));
        if isnan(thetaZ)
            thetaZ = 0;
        elseif abs(thetaY) == 90 && thetaX ~= 0 && thetaZ ~= 0
            choice = questdlg('When absolute Rot. Y is 90°, Rot. X and Rot. Z are on the same direction. For that reason, only the specified value for Rot. X will be set.',...
                'Attention','OK','Set Rot. Z instead','OK');
            if strcmp(choice,'Set Rot. Z instead')
                thetaX = 0;
            else
                thetaZ = 0;
            end
        end
    end
    
    % Convert from degrees to rad
    thetaX = thetaX * pi / 180;
    if abs(thetaX) < 10^-6
        thetaX = 0;
    end
    
    thetaY = thetaY * pi / 180;
    if abs(thetaY) < 10^-6
        thetaY = 0;
    end
    
    thetaZ = thetaZ * pi / 180;
    if abs(thetaZ) < 10^-6
        thetaZ = 0;
    end
    
    % Assemble direction vector vx
    vx = [cos(thetaZ)*cos(thetaY), sin(thetaZ)*cos(thetaY), sin(thetaY)];
    id = find((abs(vx) < 10^-6));
    if ~isempty(id)
        vx(id) = zeros(1,length(id));
    end
    
    % Compute defalut local axis y direction
    if anm == 4 || anm == 5
        % Auxiliar vz node local axis z orientation versor
        vz = [sin(thetaX)*sin(thetaZ) - cos(thetaX)*cos(thetaZ)*sin(thetaY), - sin(thetaX)*cos(thetaZ) - cos(thetaX)*sin(thetaY)*sin(thetaZ), cos(thetaX)*cos(thetaY)];
        % Compute direction vector vy
        vy = cross(vz,vx);
        vy = vy / norm(vy);
        id = find((abs(vy) < 10^-6));
        if ~isempty(id)
            vy(id) = zeros(1,length(id));
        end
    else
        vy = [];
    end

    % Enable edtable text boxes
    set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string',num2str(vx(1),3));
    set(handles.edit_InclSupp_2,'Enable','on','Visible','on','string',num2str(vx(2),3));
    if anm == 1 || anm == 2
        set(handles.edit_InclSupp_3,'Enable','off','Visible','on','string','0');
    elseif anm == 4 || anm == 5
        set(handles.edit_InclSupp_3,'Enable','on','Visible','on','string',num2str(vx(3),3));
        setVy(handles,[],true,vy);
    end
    set(handles.text_InclSupp_1,'Visible','on','string','vx_X');
    set(handles.text_InclSupp_2,'Visible','on','string','vx_Y');
    set(handles.text_InclSupp_3,'Visible','on','string','vx_Z');
    set(handles.text_Rot_1,'Visible','off');
    set(handles.text_Rot_2,'Visible','off');
    set(handles.text_Rot_3,'Visible','off');

    if get(handles.checkbox_AxesInclSupp3D,'value')
        updateInclSupp3DCnvs(handles,true);
    end
end

%--------------------------------------------------------------------------
function radiobutton_Angles_Callback(hObject, eventdata, handles)
if ~get(hObject,'value')
    if ~get(handles.radiobutton_DirectionVector,'value')
        set(hObject,'value',1)
        return
    end
else
    % Unselect direction vector radio button
    set(handles.radiobutton_DirectionVector,'Enable','on','value',0);
    
    % Get analysis model id
    model = getappdata(0,'model');
    anm = model.anm.analysis_type + 1;
    clear model
    
    % Check if current analysis mode is 2D or 3D
    if anm ~= 4 && anm ~= 5
        % Compute direction vector as angle
        dx = str2double(get(handles.edit_InclSupp_1,'string'));
        if isnan(dx)
            dx = 0;
        end
        dy = str2double(get(handles.edit_InclSupp_2,'string'));
        if isnan(dy)
            dy = 0;
        end
        
        if dx == 0 && dy == 0
            dir = [1, 0];
        else
            dir = [dx, dy] / norm([dx, dy]);
        end
        
        thetaZ = acos(dir(1));
        if dir(2) < 0
            thetaZ = - thetaZ;
        end

        % Convert from rad to degrees
        thetaZ = thetaZ * 180 / pi;
        if abs(thetaZ) < 10^-10
            thetaZ = 0;
        end

    else
        % Compute direction vector as angle
        dx = str2double(get(handles.edit_InclSupp_1,'string'));
        if isnan(dx)
            dx = 0;
        end
        dy = str2double(get(handles.edit_InclSupp_2,'string'));
        if isnan(dy)
            dy = 0;
        end
        dz = str2double(get(handles.edit_InclSupp_3,'string'));
        if isnan(dz)
            dz = 0;
        end
        
        if dx == 0 && dy == 0 && dz == 0
            dir_1 = [1, 0];
            dir_2 = [1, 0];
            dir = [1, 0, 0];
        else
            dir_1 = [dx, dy];
            dir_2 = [norm(dir_1), dz];
            dir_1 = dir_1 / norm(dir_1);
            dir_2 = dir_2 / norm(dir_2);
            dir = [dx, dy, dz] / norm([dx, dy, dz]);
        end

        % Compute thetaZ
        thetaZ = acos(dir_1(1));
        if dir_1(2) < 0
            thetaZ = - thetaZ;
        end

        % Convert from rad to degrees
        thetaZ = thetaZ * 180 / pi;
        if abs(thetaZ) < 10^-10
            thetaZ = 0;
        end

        % Compute thetaY
        thetaY = asin(dir_2(2));

        % Convert from rad to degrees
        thetaY = thetaY * 180 / pi;
        if abs(thetaY) < 10^-10
            thetaY = 0;
        end
        
        % Round angle values, to avoid numeric errors from vector to angle
        % conversion
        if abs(thetaZ) >= 10
            thetaZ = round(thetaZ,1);
        end
        if abs(thetaY) >= 10
            thetaY = round(thetaY,1);
        end
    end

    % Enable edtable text boxes
    if anm == 1 || anm == 2
        set(handles.edit_InclSupp_1,'Enable','off','Visible','off','string','0');
        set(handles.edit_InclSupp_2,'Enable','on','Visible','on','string',num2str(thetaZ,4));
        set(handles.edit_InclSupp_3,'Enable','off','Visible','off','string','0');
        set(handles.text_Rot_1,'Visible','off');
        set(handles.text_Rot_2,'Visible','on','String','Rotation Z');
        set(handles.text_Rot_3,'Visible','off');
        set(handles.text_InclSupp_1,'Visible','off');
        set(handles.text_InclSupp_3,'Visible','off');
    elseif anm == 4 || anm == 5
%         if ~get(handles.checkbox_Dy,'value') && ~get(handles.checkbox_Ky,'value') &&...
%            ~get(handles.checkbox_Dz,'value') && ~get(handles.checkbox_Kz,'value')
%             set(handles.edit_InclSupp_1,'Enable','off','Visible','on','string','');
%         else
            if abs(abs(thetaY)-90) < 10^-6 && abs(thetaZ) > 10^-6
                set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string','0');
            else
                x = dir;
                % Auxiliar vz node local axis z orientation versor
                if abs(abs(thetaY)-90) < 10^-6
                    if thetaY > 0
                        z1 = [-1, 0, 0];
                    else
                        z1 = [1, 0, 0];
                    end
                else
                    z1 = [0, 0, 1];
                end

                % Compute local axis y and z directions
                y1 = cross(z1,x);
                y1 = y1 / norm(y1);
                
                dyx = str2double(get(handles.edit_InclSupp_4,'string'));
                if isnan(dyx)
                    dyx = 0;
                end
                dyy = str2double(get(handles.edit_InclSupp_5,'string'));
                if isnan(dyy)
                    dyy = 0;
                end
                dyz = str2double(get(handles.edit_InclSupp_6,'string'));
                if isnan(dyz)
                    dyz = 0;
                end
                if dyx == 0 && dyy == 0 && dyz == 0
                    y2 = y1;
                else
                    vy = [dyx, dyy, dyz] / norm([dyx, dyy, dyz]);
                    z2 = cross(x,vy);
                    z2 = z2 / norm(z2);
                    y2 = cross(z2,x);
                end
                
                w = cross(y1,y2);
                aux = x * w';
                thetaX = acos((y1*y2')/(norm(y1)*norm(y2))) * 180 / pi;
                if abs(thetaX) < 10^-10
                    thetaX = 0;
                elseif aux < 0
                    thetaX = - thetaX;
                end
                if abs(thetaX) >= 10
                    thetaX = round(thetaX,1);
                else
                    thetaX = round(thetaX,2);
                end
                set(handles.edit_InclSupp_1,'Enable','on','Visible','on','string',num2str(thetaX,4));
            end
%        end
        set(handles.edit_InclSupp_2,'Enable','on','Visible','on','string',num2str(thetaY,4));
        if isnan(thetaZ)
            thetaZ = 0;
        end
        set(handles.edit_InclSupp_3,'Enable','on','Visible','on','string',num2str(thetaZ,4));
        set(handles.text_Rot_1,'Visible','on');
        set(handles.text_Rot_2,'Visible','on','String','Rotation Y');
        set(handles.text_Rot_3,'Visible','on');
        set(handles.text_InclSupp_1,'Visible','on','string','degrees');
        set(handles.text_InclSupp_3,'Visible','on','string','degrees');
    end
    set(handles.text_InclSupp_2,'Visible','on','string','degrees');
    
    % Hide vy input options
    if anm == 4 || anm == 5
        setVy(handles,[],false);
    end

    if get(handles.checkbox_AxesInclSupp3D,'value')
        updateInclSupp3DCnvs(handles,true);
    end
end

%--------------------------------------------------------------------------
function edit_InclSupp_Callback(hObject, eventdata, handles)
if get(handles.checkbox_AxesInclSupp3D,'value')
    updateInclSupp3DCnvs(handles,true);
end

%--------------------------------------------------------------------------
% Callback for pressed keys on the keyboard
function GUI_Supports_KeyPressFcn(hObject, eventdata, handles)
if ~get(handles.checkbox_AxesInclSupp3D,'value')
    return
end

key = get(gcf, 'CurrentKey');
switch key
    case 'escape'
        mouse_InclSupp = get(handles.axes_InclSupp3D,'UserData');
        mouse_InclSupp.dir_x_beingCollected = false;
        mouse_InclSupp.dir_x_1stPt = [];
        mouse_InclSupp.dir_x_2ndPt = [];
        delete(findobj('tag','dynamicLine_InclSupp_dir_x'))
        delete(findobj('tag','text_InclSupp_dir_x'))
        set(handles.axes_InclSupp3D,'UserData',mouse_InclSupp);
end

%--------------------------------------------------------------------------
% Auxiliary function
% Sets buttons and texts properties to fit vy input components inside the
% inclined supports uipanel.
% Input:
% * handles -> handle to struct containing gui information
% * node -> handle to node object according to selection on pop-up menu
% * flag -> indicates if operation is do (true) or undo (false)
function setVy(handles,node,flag,vy)
if nargin < 4
    vy = [];
end
if flag % Vy input is being shown
    % Set string, visualization and enabled status
%     if get(handles.checkbox_Dy,'value') || get(handles.checkbox_Dz,'value') ||...
%        get(handles.checkbox_Ky,'value') || get(handles.checkbox_Kz,'value')
        if ~isempty(vy)
            set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vy(1),3));
            set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vy(2),3));
            set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vy(3),3));
        elseif isempty(node)
            set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string','0');
            set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string','1');
            set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string','0');
        elseif ~isempty(node.inclSupp_vy)
            vyx = node.inclSupp_vy(1);
            if abs(vyx) < 10^-10
                vyx = 0;
            end
            vyy = node.inclSupp_vy(2);
            if abs(vyy) < 10^-10
                vyy = 0;
            end
            vyz = node.inclSupp_vy(3);
            if abs(vyz) < 10^-10
                vyz = 0;
            end
            set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string',num2str(vyx,3));
            set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string',num2str(vyy,3));
            set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string',num2str(vyz,3));
        else
            set(handles.edit_InclSupp_4,'Enable','on','Visible','on','string','0');
            set(handles.edit_InclSupp_5,'Enable','on','Visible','on','string','1');
            set(handles.edit_InclSupp_6,'Enable','on','Visible','on','string','0');
        end
%     else
%         set(handles.edit_InclSupp_4,'Enable','off','Visible','on','string','');
%         set(handles.edit_InclSupp_5,'Enable','off','Visible','on','string','');
%         set(handles.edit_InclSupp_6,'Enable','off','Visible','on','string','');
%     end
    set(handles.text_InclSupp_4,'Visible','on');
    set(handles.text_InclSupp_5,'Visible','on');
    set(handles.text_InclSupp_6,'Visible','on');
    
else % Vy input is being hidden
    set(handles.edit_InclSupp_4,'Enable','off','Visible','off','string','');
    set(handles.edit_InclSupp_5,'Enable','off','Visible','off','string','');
    set(handles.edit_InclSupp_6,'Enable','off','Visible','off','string','');
    set(handles.text_InclSupp_4,'Visible','off');
    set(handles.text_InclSupp_5,'Visible','off');
    set(handles.text_InclSupp_6,'Visible','off');
end

changeVyPos(handles);

%--------------------------------------------------------------------------
% Auxiliary function
% Sets buttons and texts position to fit vy input components inside the
% inclined supports uipanel.
% Input:
% * handles -> handle to struct containing gui information
function changeVyPos(handles)
% Change components positions
old_pos = get(handles.edit_InclSupp_1,'Position');
new_pos = get(handles.edit_InclSupp_1,'UserData');
set(handles.edit_InclSupp_1,'Position',new_pos);
set(handles.edit_InclSupp_1,'UserData',old_pos);

old_pos = get(handles.edit_InclSupp_2,'Position');
new_pos = get(handles.edit_InclSupp_2,'UserData');
set(handles.edit_InclSupp_2,'Position',new_pos);
set(handles.edit_InclSupp_2,'UserData',old_pos);

old_pos = get(handles.edit_InclSupp_3,'Position');
new_pos = get(handles.edit_InclSupp_3,'UserData');
set(handles.edit_InclSupp_3,'Position',new_pos);
set(handles.edit_InclSupp_3,'UserData',old_pos);

old_pos = get(handles.text_InclSupp_1,'Position');
new_pos = get(handles.text_InclSupp_1,'UserData');
set(handles.text_InclSupp_1,'Position',new_pos);
set(handles.text_InclSupp_1,'UserData',old_pos);

old_pos = get(handles.text_InclSupp_2,'Position');
new_pos = get(handles.text_InclSupp_2,'UserData');
set(handles.text_InclSupp_2,'Position',new_pos);
set(handles.text_InclSupp_2,'UserData',old_pos);

old_pos = get(handles.text_InclSupp_3,'Position');
new_pos = get(handles.text_InclSupp_3,'UserData');
set(handles.text_InclSupp_3,'Position',new_pos);
set(handles.text_InclSupp_3,'UserData',old_pos);

old_pos = get(handles.edit_InclSupp_4,'Position');
new_pos = get(handles.edit_InclSupp_4,'UserData');
set(handles.edit_InclSupp_4,'Position',new_pos);
set(handles.edit_InclSupp_4,'UserData',old_pos);

old_pos = get(handles.edit_InclSupp_5,'Position');
new_pos = get(handles.edit_InclSupp_5,'UserData');
set(handles.edit_InclSupp_5,'Position',new_pos);
set(handles.edit_InclSupp_5,'UserData',old_pos);

old_pos = get(handles.edit_InclSupp_6,'Position');
new_pos = get(handles.edit_InclSupp_6,'UserData');
set(handles.edit_InclSupp_6,'Position',new_pos);
set(handles.edit_InclSupp_6,'UserData',old_pos);

%--------------------------------------------------------------------------
% Auxiliary function
% Checks if inclined support buttons and texts postion is up or down.
% Input:
% * handles -> handle to struct containing gui information
% Output:
% * bttnsAreUp -> flag (true or false)
function bttnsAreUp = checkVyPos(handles)
current_pos = get(handles.edit_InclSupp_1,'Position');
stored_pos = get(handles.edit_InclSupp_1,'UserData');
if current_pos(2) > stored_pos(2)
    bttnsAreUp = true;
else
    bttnsAreUp = false;
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
