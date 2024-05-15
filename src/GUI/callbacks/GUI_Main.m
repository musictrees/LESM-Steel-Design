%% Main Dialog Callback Functions
% This file contains the callback functions associated with the main
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Main(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Main_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Main_OutputFcn, ...
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
% Executes just before main GUI is made visible.
% Sets GUI initial properties.
function GUI_Main_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Main
handles.output = hObject;

% Create progress bar for opening interface
prog = waitbar(.0,'Please wait...','Name','Starting LESM');

% Move GUI to the center of the screen
movegui(gcf,'center')

% Disable camera rotation option
set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')

% Set stamp figure
set(hObject,'CurrentAxes',handles.axes_Stamp);
imshow('logo_lesm_stamp.jpg')
axis image
set(handles.axes_Stamp,'HandleVisibility','off') % make it inacessible to gcf or gca

waitbar(.2,prog);

% Change the default configurations for pan and zoom
z = zoom;
p = pan;
z.ActionPreCallback = @zoomprecallback;
p.ActionPreCallback = @panprecallback;

% Disable axes default interactivity
disableDefaultInteractivity(handles.axes_Canvas);

% Set callbacks for pushbuttons that will handle tabs
%set(handles.pushbutton_TabModeling,'Callback',@pushbutton_TabModeling_Callback);
%set(handles.pushbutton_TabResults,'Callback',@pushbutton_TabResults_Callback);

set(handles.uipanel_Modeling,'Visible','on');
%set(handles.uipanel_Results,'Visible','off','Position',get(handles.uipanel_Modeling,'Position'));

% Set axes limits and ratio
set(hObject,'CurrentAxes',handles.axes_Canvas);
view(2)
axis equal
xlim([-5,5])

% Turn on grid/ruler
grid on
handles.axes_Canvas.XAxis.Visible = 'on';
handles.axes_Canvas.YAxis.Visible = 'on';
handles.axes_Canvas.ZAxis.Visible = 'on';

waitbar(.4,prog);

% Include global constants
include_constants;

% Set initial analysis type as linear elastic static
set(handles.popupmenu_AnalysisType,'value',1);

% Set initial analysis model to frame 2D
set(handles.popupmenu_Anm,'Value',2);

% Material (100 GPa by default)
initialMat = Material(1,100*10^6,0.3,10^-5,1.0);
%RONALD INICIADO O LESM COM MATERIAL DE Aï¿½O
initialMat = Steel(1,200000,0.32,12*1e-6,7.85,250,400);

% Section
initialSec = Section(1,0.01,0.01,0.01,0.00001,0.00001,0.00001,0.1,0.1);
initialSec = W_Beam(1,30,5,10,60); % section do Ronald
%Ronald add mais doi paramatros defaul type="wbeam" s = Section(nsec,Ax,Ay,Az,Ix,Iy,Iz,Hy,Hz);

% Initialize variables and save them in root
% fullname = path to file
% nmat = number of materials
% nsec = number of cross-sections
% nnp = number of nodal points
% nel = number of elements
% currentLc = current load case
% decPrec = decimal precision
setappdata(0,'fullname',{})
setappdata(0,'nnp',0);
setappdata(0,'nel',0);
setappdata(0,'currentLc',1);
setappdata(0,'decPrec',1);

waitbar(.6,prog);

% Initialize objects
model = Model();
print = Print_Frame2D();
draw = Draw_Frame2D();
drv = Drv_LES(true,model);
mouse = Emouse_2D(gcf,gca);

% Set object properties
model.drv = drv;
model.anm = Anm_Frame2D();
model.nmat = 1;
model.materials = initialMat;
model.nsec = 1;
model.sections = initialSec;
model.nlc = 1;
model.strLc = {'CASE 01'};
model.whichSolver = STATIC_LINEAR;
print.model = model;
draw.mdl = model;

% Save objects in root (make them accessible to all GUIs)
setappdata(0,'model',model)
setappdata(0,'nmat',1)
setappdata(0,'materials',initialMat)
setappdata(0,'nsec',1)
setappdata(0,'sections',initialSec)
setappdata(0,'print',print)
setappdata(0,'draw',draw)
setappdata(0,'mouse',mouse)

% Update information panel
infoPanelData(1,:) = {'Materials',1};
infoPanelData(2,:) = {'Cross-Sections',1};
infoPanelData(3,:) = {'Nodes',0};
infoPanelData(4,:) = {'Elements',0};
infoPanelData(5,:) = {'DOFs',0};
infoPanelData(6,:) = {'Free DOFs',0};
infoPanelData(7,:) = {'Fixed DOFs',0};
infoPanelData(8,:) = {'Springs',0};
infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',0};
set(handles.uitable_infoPanel,'Data',infoPanelData)
set(handles.uitable_infoPanelEditable,'Data',{},'CellEditCallback',@uitable_infoPanelEditable_CellEditCallback,'enable','off')

% Initialize flag for allowing visualization options
setappdata(0,'vis',0);

% Adjust canvas position
dfltUnits = get(handles.axes_Canvas, 'Units');
set(handles.axes_Canvas, 'Units', 'normalized');
set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
set(handles.axes_Canvas, 'Units', dfltUnits);

% Dynamic results reproduction speed
handles.dynResSpeed = 1;

waitbar(.8,prog);

% Create buttons in toolbar for cursor coordinates indication
hToolbar   = findall(hObject,'tag','toolbar');
jToolbar   = hToolbar.JavaContainer.getComponentPeer;
size_axis  = java.awt.Dimension(20,23);
size_value = java.awt.Dimension(60,23);

axisx = uipushtool(hToolbar);
set(axisx,'Enable','off','Separator','on');
drawnow;
jaxisx = jToolbar.getComponent(jToolbar.getComponentCount-1);
jaxisx.setText('X:');
jaxisx.setMaximumSize(size_axis);
jaxisx.setPreferredSize(size_axis);
jaxisx.setSize(size_axis);

coordx = uipushtool(hToolbar);
set(coordx,'Tag','text_coordx','Enable','off');
drawnow;
jcoordx = jToolbar.getComponent(jToolbar.getComponentCount-1);
jcoordx.setText(' ');
jcoordx.setMaximumSize(size_value);
jcoordx.setPreferredSize(size_value);
jcoordx.setSize(size_value);

axisy = uipushtool(hToolbar);
set(axisy,'Enable','off');
drawnow;
jaxisy = jToolbar.getComponent(jToolbar.getComponentCount-1);
jaxisy.setText('Y:');
jaxisy.setMaximumSize(size_axis);
jaxisy.setPreferredSize(size_axis);
jaxisy.setSize(size_axis);

coordy = uipushtool(hToolbar);
set(coordy,'Tag','text_coordy','Enable','off');
drawnow;
jcoordy = jToolbar.getComponent(jToolbar.getComponentCount-1);
jcoordy.setText(' ');
jcoordy.setMaximumSize(size_value);
jcoordy.setPreferredSize(size_value);
jcoordy.setSize(size_value);

waitbar(1,prog);
close(prog);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Main_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes on button press in "Pan" pushbutton of toolbar.
% Change cursor for pan
function panprecallback(varargin)
setAxes3DPanAndZoomStyle(pan,gca,'limits');

%--------------------------------------------------------------------------
% Executes on button press in "Zoom In" and "Zoom out" pushbuttons of toolbar.
% Change cursor for zoom
function zoomprecallback(varargin)
setAxes3DPanAndZoomStyle(pan,gca,'camera');

% --------------------------------------------------------------------
function fileMenu_Callback(~, ~, ~) %#ok<DEFNU>

% --------------------------------------------------------------------
function newButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

% Open dialog box to confirm user action
choice = questdlg('All unsaved data will be lost. Do you want to continue?','New Model','OK','Cancel','Cancel');

if strcmp(choice,'OK')
    % Disable save button
    setappdata(0,'fullname',{})
    setappdata(0,'filename',{})
    %set(handles.saveButton,'Enable','off')
    
    % Turn off 2D view togglebutton
    if strcmp(get(handles.togglebutton_2DView,'state'),'on')
        setappdata(0,'reset3DViewFlag',1)
        set(handles.togglebutton_2DView,'state','off')
    end
    
    % Change window title
    set(gcf,'Name','LESM - Linear Elements Structure Model')
    
    % Clean canvas
    axes(handles.axes_Canvas);
    cla reset
    
    % Reset toggle buttons on toolbar
    set(handles.togglebutton_Node,'enable','on','state','off')
    set(handles.togglebutton_Element,'enable','on','state','off')
    set(handles.togglebutton_CrossElements,'state','off')
    set(handles.togglebutton_Polyline,'state','off')
    set(handles.togglebutton_Ortho,'state','off')
    set(handles.pushbutton_SolveIntSects,'enable','off')
    
    % Reset model visualization options
    set(handles.viewNodalLoadsButton,'Checked','on')
    set(handles.viewNodalMassButton,'Checked','on')
    set(handles.viewPrescDisplButton,'Checked','on')
    set(handles.viewInitialConditionsButton,'Checked','on')
    set(handles.viewDistribLoadsButton,'Checked','on')
    set(handles.viewThermalLoadsButton,'Checked','on')
    set(handles.viewSupportsButton,'Checked','on')
    set(handles.viewSemiRigidButton,'Checked','on')
    
    % Clear coordinates
    hToolbar = findall(groot,'tag','toolbar');
    jToolbar = hToolbar.JavaContainer.getComponentPeer;
    jaxisx   = jToolbar.getComponent(jToolbar.getComponentCount-3);
    jaxisy   = jToolbar.getComponent(jToolbar.getComponentCount-1);
    jaxisx.setText(' ');
    jaxisy.setText(' ');
    
    % Set properties that depend on the selected analysis model
    anm = get(handles.popupmenu_Anm,'value');
    if anm == 1
        anm = Anm_Truss2D();
        print = Print_Truss2D();
        draw = Draw_Truss2D();
        view(2);
%         axis equal
%         xlim([0,10])
        set(handles.planeButton,'Enable','off');
        set(handles.togglebutton_2DView,'enable','off')
        set(handles.pushbutton_AxesLimits,'enable','off')
    elseif anm == 2
        anm = Anm_Frame2D();
        print = Print_Frame2D();
        draw = Draw_Frame2D();
        view(2);
%         axis equal
%         xlim([0,10])
        set(handles.planeButton,'Enable','off');
        set(handles.togglebutton_2DView,'enable','off')
        set(handles.pushbutton_AxesLimits,'enable','off')
    elseif anm == 3
        anm = Anm_Grillage();
        print = Print_Grillage();
        draw = Draw_Grillage();
        view(3);
        axis equal
        xlim([-5,5])
        ylim([-5,5])
        zlim([0,0.02])
        xlabel('X');
        ylabel('Y');
        zlabel(' ');
        zticks(0);
        zticklabels({' '});
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
        set(handles.togglebutton_2DView,'enable','on')
        set(handles.pushbutton_AxesLimits,'enable','on')
    elseif anm == 4
        anm = Anm_Truss3D();
        print = Print_Truss3D();
        draw = Draw_Truss3D();
        view(3);
        axis equal
        xlim([0,10])
        ylim([0,10])
        zlim([-1,1])
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
        set(handles.togglebutton_2DView,'enable','off')
        set(handles.pushbutton_AxesLimits,'enable','off')
    elseif anm == 5
        anm = Anm_Frame3D();
        print = Print_Frame3D();
        draw = Draw_Frame3D();
        view(3);
        axis equal
        xlim([0,10])
        ylim([0,10])
        zlim([-1,1])
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
        set(handles.togglebutton_2DView,'enable','off')
        set(handles.pushbutton_AxesLimits,'enable','off')
    end
    
    % Material 
    initialMat = Material(1,100*10^6,0.3,10^-5,1.0);
     initialMat = Steel(1,200000,0.32,12*1e-6,7.85,250,400); % Material Ronald
    
    % Section
    initialSec = Section(1,0.01,0.01,0.01,0.00001,0.00001,0.00001,0.1,0.1);
    initialSec = W_Beam(1,30,5,10,60); % Initial Sec do Ronald -> Possível fonte de problemas por atribuir um valor diferente da linha de código anterior.
    
    % Clean model object and set its new properties
    model = getappdata(0,'model');
    model.clean();
    model.anm = anm;
    model.nmat = 1;
    model.materials = initialMat;
    model.nsec = 1;
    model.sections = initialSec;
    model.nlc = 1;
    model.strLc = {'CASE 01'};
    
    % Update information panel
    infoPanelData(1,:) = {'Materials',1};
    infoPanelData(2,:) = {'Cross-Sections',1};
    infoPanelData(3,:) = {'Nodes',0};
    infoPanelData(4,:) = {'Elements',0};
    infoPanelData(5,:) = {'DOFs',0};
    infoPanelData(6,:) = {'Free DOFs',0};
    infoPanelData(7,:) = {'Fixed DOFs',0};
    infoPanelData(8,:) = {'Springs',0};
    infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',0};
    set(handles.uitable_infoPanel,'Data',infoPanelData)
    set(handles.uitable_infoPanelEditable,'Data',{},'CellEditCallback',@uitable_infoPanelEditable_CellEditCallback,'enable','off')
    set(handles.pushbutton_ApplyInfoPanel,'enable','off')
    
    % Deactivate visualizing tools
    pan off
    rotate3d off
    zoom off
    
    % Enable modeling buttons 
    set(handles.popupmenu_Anm,'Enable','on');
    set(handles.popupmenu_AnalysisType,'Enable','on');
    set(handles.pushbutton_Materials,'Enable','on');
    set(handles.pushbutton_Sections,'Enable','on');
    set(handles.pushbutton_Nodes,'Enable','on');
    set(handles.pushbutton_Elements,'Enable','on');
    set(handles.pushbutton_NodalLoads,'Enable','on');
    set(handles.pushbutton_ElementLoads,'Enable','on');
    set(handles.pushbutton_Supports,'Enable','on');
    
    % Disable "Process Data" button
    set(handles.pushbutton_ProcessData,'Enable','off');
    
    % Disable result options and adjust colors
    set(handles.popupmenu_Results,'Enable','off','value',1);
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.edit_Scale,'Visible','off','Enable','off','String',' ');
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    
    % Create drv object and enable according modeling options
    if get(handles.popupmenu_AnalysisType,'Value') == 1
        drv = Drv_LES(true,model);
        staticAnalysisCallback(handles);
        set(handles.popupmenu_LoadCase,'enable','on','string','CASE 01','value',1,'Max',1);
    elseif get(handles.popupmenu_AnalysisType,'Value') == 2
        drv = Drv_LED(1,true,model);
        dynamicAnalysisCallback(handles);
        set(handles.popupmenu_LoadCase,'enable','off','value',1,'string',' ');
    end
    model.drv = drv;
    print.model = model;
    draw.mdl = model;
    
    % Update variables in root
    setappdata(0,'model',model);
    setappdata(0,'print',print);
    setappdata(0,'draw',draw);
    setappdata(0,'nmat',1);
    setappdata(0,'materials',initialMat);
    setappdata(0,'nsec',1);
    setappdata(0,'sections',initialSec);
    setappdata(0,'nnp',0);
    setappdata(0,'nel',0);
    setappdata(0,'currentLc',1);
    setappdata(0,'intersections',[]);
    
    % Turn grid on/off
    if strcmp(get(handles.gridButton,'Checked'),'on') == 1
        grid on
    else
        grid off
    end
    
    % Turn ruler on/off, if checkbox is clicked
    if strcmp(get(handles.rulerButton,'Checked'),'on') == 1
        handles.axes_Canvas.XAxis.Visible = 'on';
        handles.axes_Canvas.YAxis.Visible = 'on';
        handles.axes_Canvas.ZAxis.Visible = 'on';
    else
        handles.axes_Canvas.XAxis.Visible = 'off';
        handles.axes_Canvas.YAxis.Visible = 'off';
        handles.axes_Canvas.ZAxis.Visible = 'off';
    end
    
    % Reset flag for allowing visualization options
    setappdata(0,'vis',0);
    
    % Adjust canvas position
    dfltUnits = get(handles.axes_Canvas, 'Units');
    set(handles.axes_Canvas, 'Units', 'normalized');
    if (anm.analysis_type == 0) || (anm.analysis_type == 1)
        set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
    else
        set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
    end
    set(handles.axes_Canvas, 'Units', dfltUnits);
    
    % Reinitialize object for mouse events and save it in root
    if (anm.analysis_type == 0) || (anm.analysis_type == 1)
        mouse = Emouse_2D(gcf,handles.axes_Canvas);
    else
        mouse = Emouse_3D(gcf,handles.axes_Canvas);
    end
    setappdata(0,'mouse',mouse);
    
    % Refresh draw to avoid issues
    % (it corrects a bug of distorted draws after opening new blank model)
    redraw(handles);
end

% Save changes to handles structure
guidata(hObject,handles)

% --------------------------------------------------------------------
function openButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

include_constants;
model = getappdata(0,'model');

if model.nnp ~= 0
    choice = questdlg('All unsaved data will be lost. Do you want to continue?',...
        'Open file','Save and open','Open without saving','Cancel','Cancel');

    switch choice
        case 'Save and open'
            filename = saveButton_Callback([],eventdata,handles);
            if filename == 0
                return
            end
        case 'Cancel'
            return
    end
else
    choice = 'Open without saving';
end

if ~strcmp(choice,'Save and open') && ~strcmp(choice,'Open without saving')
    return
end

% Open neutral-format file and get file id
filterspec = {'*.lsm'};
[filename,pathname] = uigetfile(filterspec,'LESM - Input file');
if isequal(filename,0)
    fid = 0;
else
    fullname = strcat(pathname,filename);
    fid = fopen(fullname,'rt');
end

% Create window title
filename = fopen(fid);
[~,name,~] = fileparts(filename);
setappdata(0,'filename',{name});
name = strcat('LESM - Linear Elements Structure Model -',name);

% Check if a valid file was opened
if fid ~= 0
    % Enable save button
    setappdata(0,'fullname',{fullname})
    %set(handles.saveButton,'Enable','on')
    
    % Reset toggle buttons on toobar
    set(handles.togglebutton_Node,'enable','on','state','off')
    set(handles.togglebutton_Element,'enable','on','state','off')
    set(handles.togglebutton_CrossElements,'state','off')
    set(handles.togglebutton_Polyline,'state','off')
    set(handles.togglebutton_Ortho,'state','off')
    
    % Reset model visualization options
    set(handles.viewNodalLoadsButton,'Checked','on')
    set(handles.viewNodalMassButton,'Checked','on')
    set(handles.viewPrescDisplButton,'Checked','on')
    set(handles.viewInitialConditionsButton,'Checked','on')
    set(handles.viewDistribLoadsButton,'Checked','on')
    set(handles.viewThermalLoadsButton,'Checked','on')
    set(handles.viewSupportsButton,'Checked','on')
    set(handles.viewSemiRigidButton,'Checked','on')
    
    % Clear coordinates
    hToolbar = findall(groot,'tag','toolbar');
    jToolbar = hToolbar.JavaContainer.getComponentPeer;
    jaxisx   = jToolbar.getComponent(jToolbar.getComponentCount-3);
    jaxisy   = jToolbar.getComponent(jToolbar.getComponentCount-1);
    jaxisx.setText(' ');
    jaxisy.setText(' ');
    
    % Create new model object
    model = Model();
    
    % Check if there were unsolved intersections on previous model
    ints = getappdata(0,'intersections');
    if ~isempty(ints)
        thereWereIntersections = true;
    else
        thereWereIntersections = false;
    end
    setappdata(0,'intersections',[]);
    
    % Read input file data
    [vs,print,draw,nclc] = readFile(fid,model,true,pathname);
    
    % Check input file version compatibility
    if vs == 1
        % Set window title
        set(gcf,'Name',name)

        % Turn off 2D view togglebutton
        if strcmp(get(handles.togglebutton_2DView,'state'),'on')
            setappdata(0,'reset3DViewFlag',1)
            set(handles.togglebutton_2DView,'state','off')
        end
    
        % Clean canvas
        axes(handles.axes_Canvas)
        cla reset
        
        % Update variables in root
        setappdata(0,'nmat',model.nmat);
        setappdata(0,'materials',model.materials);
        setappdata(0,'nsec',model.nsec);
        setappdata(0,'sections',model.sections);
        setappdata(0,'nnp',model.nnp);
        setappdata(0,'nodes',model.nodes);
        setappdata(0,'nel',model.nel);
        setappdata(0,'elems',model.elems);
        setappdata(0,'currentLc',nclc);
        
        % Update combination properties in model object
        if isempty(model.loadComb)
            model.ncomb = 0;
            model.strComb = ' ';
        end
        
        % Update information panel
        infoPanelData(1,:) = {'Materials',getappdata(0,'nmat')};
        infoPanelData(2,:) = {'Cross-Sections',getappdata(0,'nsec')};
        infoPanelData(3,:) = {'Nodes',getappdata(0,'nnp')};
        infoPanelData(4,:) = {'Elements',getappdata(0,'nel')};
        infoPanelData(5,:) = {'DOFs',model.neq};
        infoPanelData(6,:) = {'Free DOFs',(model.neq - model.neqfixed - model.neqspring)};
        infoPanelData(7,:) = {'Fixed DOFs',model.neqfixed};
        infoPanelData(8,:) = {'Springs',model.neqspring};
        infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',model.njoints * model.anm.nrdof};
        set(handles.uitable_infoPanel,'Data',infoPanelData)
        set(handles.uitable_infoPanelEditable,'Data',{},'CellEditCallback',@uitable_infoPanelEditable_CellEditCallback,'enable','off')
        set(handles.pushbutton_ApplyInfoPanel,'enable','off')
        
        % Deactivate visualizing tools
        pan off
        rotate3d off
        zoom off

        % Set properties of interface objects that depend on analysis model
        anm = model.anm.analysis_type;
        if anm == 0
            set(handles.popupmenu_Anm,'Value',1)
            view(2)
            set(handles.planeButton,'Enable','off')
            set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')
            set(handles.togglebutton_2DView,'enable','off')
            set(handles.pushbutton_AxesLimits,'enable','off')
            ax = handles.axes_Canvas;
            ax.Clipping = 'on';
        elseif anm == 1
            set(handles.popupmenu_Anm,'Value',2)
            view(2)
            set(handles.planeButton,'Enable','off')
            set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')
            set(handles.togglebutton_2DView,'enable','off')
            set(handles.pushbutton_AxesLimits,'enable','off')
            ax = handles.axes_Canvas;
            ax.Clipping = 'on';
        elseif anm == 2
            set(handles.popupmenu_Anm,'Value',3)
            view(3)
            set(handles.planeButton,'Enable','on');
            set(handles.plane3D,'Checked','on');
            set(handles.planeXY,'Checked','off');
            set(handles.planeXZ,'Checked','off');
            set(handles.planeYZ,'Checked','off');
            set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
            set(handles.togglebutton_2DView,'enable','on')
            set(handles.pushbutton_AxesLimits,'enable','on')
            ax = handles.axes_Canvas;
            ax.Clipping = 'off';
        elseif anm == 3
            set(handles.popupmenu_Anm,'Value',4)
            view(3)
            set(handles.planeButton,'Enable','on');
            set(handles.plane3D,'Checked','on');
            set(handles.planeXY,'Checked','off');
            set(handles.planeXZ,'Checked','off');
            set(handles.planeYZ,'Checked','off');
            set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
            set(handles.togglebutton_2DView,'enable','off')
            set(handles.pushbutton_AxesLimits,'enable','off')
            ax = handles.axes_Canvas;
            ax.Clipping = 'off';
        elseif anm == 4
            set(handles.popupmenu_Anm,'Value',5)
            view(3)
            set(handles.planeButton,'Enable','on');
            set(handles.plane3D,'Checked','on');
            set(handles.planeXY,'Checked','off');
            set(handles.planeXZ,'Checked','off');
            set(handles.planeYZ,'Checked','off');
            set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
            set(handles.togglebutton_2DView,'enable','off')
            set(handles.pushbutton_AxesLimits,'enable','off')
            ax = handles.axes_Canvas;
            ax.Clipping = 'off';
        end
        
        % Disable analysis model type option if there is at least one node
        if model.nnp > 0
            set(handles.popupmenu_Anm,'Enable','off');
        else
            set(handles.popupmenu_Anm,'Enable','on');
        end
        
        % Enable modeling buttons
        if nclc <= model.nlc
            set(handles.popupmenu_AnalysisType,'Enable','on');
            set(handles.pushbutton_Materials,'Enable','on');
            set(handles.pushbutton_Sections,'Enable','on');
            set(handles.pushbutton_Nodes,'Enable','on');
            set(handles.pushbutton_Elements,'Enable','on');
            set(handles.pushbutton_NodalLoads,'Enable','on');
            set(handles.pushbutton_ElementLoads,'Enable','on');
            set(handles.pushbutton_Supports,'Enable','on');
        else
            set(handles.popupmenu_AnalysisType,'Enable','off');
            set(handles.pushbutton_Materials,'Enable','off');
            set(handles.pushbutton_Sections,'Enable','off');
            set(handles.pushbutton_Nodes,'Enable','off');
            set(handles.pushbutton_Elements,'Enable','off');
            set(handles.pushbutton_NodalLoads,'Enable','off');
            set(handles.pushbutton_ElementLoads,'Enable','off');
            set(handles.pushbutton_Supports,'Enable','off');
        end
        set(handles.pushbutton_ProcessData,'Enable','on');
        
        % Check which analysis option was read
        switch model.whichSolver
            case STATIC_LINEAR
                staticAnalysisCallback(handles);
                drv = Drv_LES(true,model);
            case DYNAMIC_NEWMARK_LINEAR
                dynamicAnalysisCallback(handles);
                drv = Drv_LED(DYNAMIC_NEWMARK_LINEAR,true,model);
            case DYNAMIC_MODALSUP_LINEAR
                dynamicAnalysisCallback(handles);
                drv = Drv_LED(DYNAMIC_MODALSUP_LINEAR,true,model);
            case DYNAMIC_RK4_LINEAR
                dynamicAnalysisCallback(handles);
                drv = Drv_LED(DYNAMIC_RK4_LINEAR,true,model);
            case DYNAMIC_AM3_LINEAR    
                dynamicAnalysisCallback(handles);
                drv = Drv_LED(DYNAMIC_AM3_LINEAR,true,model);
            case DYNAMIC_WILSON_LINEAR
                dynamicAnalysisCallback(handles);
                drv = Drv_LED(DYNAMIC_WILSON_LINEAR,true,model);
        end
        model.drv = drv;
        
        % Manage Load Case, Edit, and TimeFcn buttons in case of dynamic analysis
        if get(handles.popupmenu_AnalysisType,'Value') == 1 % static
            set(handles.pushbutton_TimeFcn,'Enable','off');
            set(handles.dynResOptButton,'Enable','off');
            set(handles.pushbutton_EditLoadCase,'Enable','on');
            set(handles.popupmenu_LoadCase,'Enable','on');
            if model.ncomb ~=0
                set(handles.popupmenu_LoadCase,'String',{char(model.strLc),char(model.strComb)},'Value',nclc,'Max',model.nlc+model.ncomb)
            else
                set(handles.popupmenu_LoadCase,'String',char(model.strLc),'Value',nclc,'Max',model.nlc)
            end
        elseif get(handles.popupmenu_AnalysisType,'Value') == 2 % dynamic
            set(handles.pushbutton_TimeFcn,'Enable','on');
            set(handles.dynResOptButton,'Enable','on');
            set(handles.pushbutton_EditLoadCase,'Enable','off');
            set(handles.popupmenu_LoadCase,'Enable','off','Value',1,'String',' ');
        end
        
        % Disable result options and adjust colors
        set(handles.popupmenu_Results,'Enable','off','value',1);
        set(handles.pushbutton_Textual,'Enable','off');
        set(handles.checkbox_Reactions,'Enable','off','Value',0);
        set(handles.edit_Scale,'Visible','off','Enable','off','String',' ');
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(handles.pushbutton_DynamicResults,'Enable','off','Visible',model.drv.analysis ~= STATIC_LINEAR);
        set(handles.pushbutton_PlayDynamicResults,'enable','off','visible',model.drv.analysis ~= STATIC_LINEAR);

        % Set flag for allowing visualization options
        setappdata(0,'vis',1);
        
        % Draw model (nodes, elements, supports, hinges, loads and presc. displ.)
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 1 % static
            draw.setSize();
            draw.elemLoadsScaleFactor();
            draw.model();
            draw.nodalLoads();
            draw.elemLoads();
            draw.thermalLoads();
            draw.nodalPrescDispl();
        elseif anl == 2 % dynamic
            draw.setSize();
            draw.model();
            draw.dynamicNodalLoads();
            draw.nodalMass();
            draw.nodalInitialConditions();
        end
        axis equal
        draw.setLimits();
        
        % Turn nodes ID on, if checkbox is clicked
        if strcmp(get(handles.nodeIDButton,'Checked'),'on') == 1
            draw.nodeID();
        end

        % Turn elements ID on, if checkbox is clicked
        if strcmp(get(handles.elemIDButton,'Checked'),'on') == 1
            draw.elementID();
        end

        % Turn elements orientation on, if checkbox is clicked
        if strcmp(get(handles.orientationButton,'Checked'),'on') == 1
            draw.elementOrientation();
        end
        
        % Turn grid on/off
        if strcmp(get(handles.gridButton,'Checked'),'on') == 1
            grid on
        else
            grid off
        end

        % Turn ruler on/off, if checkbox is clicked
        if strcmp(get(handles.rulerButton,'Checked'),'on') == 1
            handles.axes_Canvas.XAxis.Visible = 'on';
            handles.axes_Canvas.YAxis.Visible = 'on';
            handles.axes_Canvas.ZAxis.Visible = 'on';
        else
            handles.axes_Canvas.XAxis.Visible = 'off';
            handles.axes_Canvas.YAxis.Visible = 'off';
            handles.axes_Canvas.ZAxis.Visible = 'off';
        end

        % Adjust canvas position
        dfltUnits = get(handles.axes_Canvas, 'Units');
        set(handles.axes_Canvas, 'Units', 'normalized');
        if (anm == 0) || (anm == 1)
            set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
        else
            set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
        end
        set(handles.axes_Canvas, 'Units', dfltUnits);
        
        % Get canvas borders
        dfltUnits = get(gca,'units');
        set(gca,'units','normalized');
        limits = get(gca,'Position');
        set(gca,'units',dfltUnits);
        axisWidth = limits(3);
        
        % Update unsolved intersections
        if ~isempty(getappdata(0,'intersections'))
            set(handles.pushbutton_SolveIntSects,'enable','on');
        else
            for e = 1:model.nel
                elemCoords = [model.elems(e).nodes(1).coord(1) model.elems(e).nodes(1).coord(2) model.elems(e).nodes(1).coord(3);
                              model.elems(e).nodes(2).coord(1) model.elems(e).nodes(2).coord(2) model.elems(e).nodes(2).coord(3)];
                crEPtsOut = auxModelFctn('getCrossElemPoints',{elemCoords,[],[]});
                crossPoints = crEPtsOut{1};
                for nn = 1:size(crossPoints,1)
                    % Get crossing point coordinates
                    x = crossPoints(nn,2);
                    y = crossPoints(nn,3);
                    z = crossPoints(nn,4);
                    
                    % Get id of elements to be divided
                    whichElems = auxModelFctn('isPointInElem',[x y z]);
                    
                    % Update vector of handles to intersections
                    if ~isempty(whichElems)
                        intersections = getappdata(0,'intersections');
                        existingIntSect = 0;
                        for nis = 1:size(intersections,2)
                            if norm(intersections(nis).coord - [x y z]) <= axisWidth/25
                                existingIntSect = nis;
                                break
                            end
                        end
                        if existingIntSect == 0
                            newIntSect.coord = [x y z];
                            newIntSect.elems = whichElems;
                            if ~isempty(intersections)
                                intSects = horzcat(intersections,newIntSect);
                            else
                                intSects = newIntSect;
                            end
                        else
                            intersections(existingIntSect).elems = whichElems;
                            intSects = intersections;
                        end
                        setappdata(0,'intersections',intSects);
                    end
                end
            end
            % Enable/disable solve intersections pushbutton (toolbar)
            if size(getappdata(0,'intersections'),2) >= 1
                set(handles.pushbutton_SolveIntSects,'enable','on')
            else
                set(handles.pushbutton_SolveIntSects,'enable','off')
            end
        end

        % Reinitialize object for mouse events and save it in root
        if (anm == 0) || (anm == 1)
            mouse = Emouse_2D(gcf,handles.axes_Canvas);
            mouse.sizeFlag = draw.size;
        else
            mouse = Emouse_3D(gcf,handles.axes_Canvas);
        end
        
        % Return objects to root
        setappdata(0,'mouse',mouse);
        setappdata(0,'model',model);
        setappdata(0,'print',print);
        setappdata(0,'draw',draw);
        
        % Rafael Rangel - During verifications before releasing 3.0:
        % I put the next call here to avoid a strange behavior of the canvas
        % when using the snap-to-grid while inserting nodes
        % (the axes limits were changing, but it was fixed if fit-to-view was previously called)
        pushbutton_FitWorld_ClickedCallback(handles.pushbutton_FitWorld,eventdata,handles);
    else
        if thereWereIntersections
            setappdata(0,'intersections',ints);
        end
        msgbox('This file version is not compatible with the program version!', 'Error','error');
    end
end

% --------------------------------------------------------------------
function filename = saveButton_Callback(~,eventdata,handles)
unselectEntities(handles);

gui_Name = get(gcf,'Name');
model = getappdata(0,'model');

if strcmp(gui_Name,'LESM - Linear Elements Structure Model')
    filename = saveAsButton_Callback(handles.saveAsButton, eventdata, handles);
else
    set(handles.popupmenu_LoadCase,'Value',1);
    popupmenu_LoadCase_Callback(handles.popupmenu_LoadCase, eventdata, handles);
    % Get file name
    fullname = char(getappdata(0,'fullname'));
    filename = gui_Name;
    % Save file
    fid = fopen(fullname,'wt');
    saveFile(model,fid);
    fclose(fid);
end

% --------------------------------------------------------------------
function filename = saveAsButton_Callback(hObject, eventdata, handles)
unselectEntities(handles);

filterspec = '*.lsm';
DialogTitle = 'LESM - Save file';
DefaultName = 'untitled';
[filename,pathname] = uiputfile(filterspec,DialogTitle,DefaultName);

if filename ~= 0
    model = getappdata(0,'model');
    fullname = strcat(pathname,filename);
    save_time_tables(model,pathname);
    fid = fopen(fullname,'wt');
    set(handles.popupmenu_LoadCase,'Value',1);
    popupmenu_LoadCase_Callback(handles.popupmenu_LoadCase, eventdata, handles);    
    saveFile(model,fid);
    fclose(fid);
    
    % Enable save button
    setappdata(0,'fullname',{fullname})
    setappdata(0,'filename',{filename})
    %set(handles.saveButton,'Enable','on')
    
    % Change window title
    name = strcat('LESM - Linear Elements Structure Model -',filename);
    set(handles.GUI_Main,'Name',name)
    setappdata(0,'model',model)
    
    % Reinitialize object for mouse events and save it in root
    anm = get(handles.popupmenu_Anm,'value');
    if (anm == 1) || (anm == 2)
        mouse = Emouse_2D(gcf,handles.axes_Canvas);
        draw = getappdata(0,'draw');
        mouse.sizeFlag = draw.size;
    else
        mouse = Emouse_3D(gcf,handles.axes_Canvas);
    end
    setappdata(0,'mouse',mouse);
end

% --------------------------------------------------------------------
function saveFigureButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
unselectEntities(handles);

name = get(handles.GUI_Main,'Name');
if strcmp(name,'LESM - Linear Elements Structure Model')
    name = 'untitled';
else
    name = name(length('LESM - Linear Elements Structure Model -')+1:end-4);
end

filterspec = {'*.jpeg';'*.png';'*.eps';'*.emf';'*.svg'};
DialogTitle = 'LESM - Save figure';
[filename,pathname,indx] = uiputfile(filterspec,DialogTitle,name);

if filename ~= 0
    fullname = strcat(pathname,filename);
    ax_old = handles.axes_Canvas;
    screensize = get(groot,'Screensize');
    f_new = figure;
    set(f_new, 'Visible', 'off');
    set(f_new,'Position',screensize);
    ax_new = copyobj(ax_old,f_new);
    
    if indx == 1
        saveas(ax_new,fullname,'jpeg')
	elseif indx == 2
        saveas(ax_new,fullname,'png')
	elseif indx == 3
        saveas(ax_new,fullname,'eps')
	elseif indx == 4
        saveas(ax_new,fullname,'meta')
	elseif indx == 5
        saveas(ax_new,fullname,'svg')
    end
    
    close(f_new)
end

% --------------------------------------------------------------------
function aboutButton_Callback(~,~,handles) %#ok<DEFNU>
unselectEntities(handles);
GUI_About

%--------------------------------------------------------------------------
function togglebutton_Cursor_OnCallback(hObject, eventdata, handles)
toggleButtonsOff(hObject,handles)
unselectEntities(handles);
delete(findobj('tag','snapGrid'));
set(handles.GUI_Main,'Pointer','arrow');
mouse = getappdata(0,'mouse');
mouse.mouseCursor = get(hObject,'state');
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function togglebutton_Cursor_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
flag = areAllToggleOff(handles);

if flag == 1
    set(hObject,'state','on')
else
    mouse = getappdata(0,'mouse');
    mouse.mouseCursor = get(hObject,'state');
    setappdata(0,'mouse',mouse)
end

%--------------------------------------------------------------------------
function toolbar_pan_OnCallback(hObject, eventdata, handles) %#ok<DEFNU>
% ax = handles.axes_Canvas;
% if get(handles.popupmenu_Anm,'value') >= 4
%     ax.Clipping = 'off';
% elseif get(handles.popupmenu_Anm,'value') == 3 &&...
%        strcmp(get(handles.togglebutton_2DView,'state'),'off')
%     ax.Clipping = 'off';
% end
toggleButtonsOff(hObject,handles)
unselectEntities(handles)
pan off

%--------------------------------------------------------------------------
function toolbar_pan_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
pan off

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

%--------------------------------------------------------------------------
function toolbar_zoomIn_OnCallback(hObject, eventdata, handles) %#ok<DEFNU>
toggleButtonsOff(hObject,handles)
unselectEntities(handles)

%--------------------------------------------------------------------------
function toolbar_zoomIn_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
if ~strcmp(get(handles.toolbar_zoomOut,'state'),'on')
    zoom off
end

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

%--------------------------------------------------------------------------
function toolbar_zoomOut_OnCallback(hObject, eventdata, handles) %#ok<DEFNU>
toggleButtonsOff(hObject,handles)
unselectEntities(handles)

%--------------------------------------------------------------------------
function toolbar_zoomOut_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
if ~strcmp(get(handles.toolbar_zoomIn,'state'),'on')
    zoom off
end

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

%--------------------------------------------------------------------------
function toolbar_rotate3d_OnCallback(hObject, eventdata, handles) %#ok<DEFNU>
toggleButtonsOff(hObject,handles)
unselectEntities(handles)

%--------------------------------------------------------------------------
function toolbar_rotate3d_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
rotate3d off

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

%--------------------------------------------------------------------------
function pushbutton_FitWorld_ClickedCallback(~, ~, ~)
mouse = getappdata(0,'mouse');
mouse.doubleClick();
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function pushbutton_RefreshModel_ClickedCallback(~, ~, handles)
redraw(handles)

%--------------------------------------------------------------------------
function togglebutton_2DView_OnCallback(hObject, ~, handles) %#ok<DEFNU>       
% Get analysis model
anm = get(handles.popupmenu_Anm,'value');
if anm ~= 3
    setappdata(0,'reset3DViewFlag',true);
    set(hObject,'state','off');
    return
end

model = getappdata(0,'model');
nclc  = getappdata(0,'currentLc');
if nclc > model.nlc
    msgbox('2D view is habilitated for modeling only, and it is not accessible for load combinations. Select a load case to access the 2D view.','Attention','warn');
    set(hObject,'state','off');
    setappdata(0,'reset3DViewFlag',true);
    return;
end

if strcmp(get(handles.pushbutton_ProcessData,'Enable'),'off')
    choice = questdlg('2D view is habilitated for modeling only, not for viewing results. Are you sure you want to continue? Proceeding will reset analysis results.',...
                      '2D view','Continue','Cancel','Cancel');
    if strcmp(choice,'Cancel')
        setappdata(0,'reset3DViewFlag',true);
        set(hObject,'state','off');
        return
    end
end

unselectEntities(handles)
if strcmp(get(handles.pushbutton_ProcessData,'Enable'),'off')
    % Enable process data button and disable results
    if getappdata(0,'nnp') ~= 0
        set(handles.pushbutton_ProcessData,'Enable','on');
    else
        set(handles.pushbutton_ProcessData,'Enable','off');
    end
    set(handles.popupmenu_Results,'Enable','off','Value',1);
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.edit_Scale,'Enable','off','Visible','off');
    set(handles.pushbutton_DynamicResults,'enable','off');
    set(handles.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    set(handles.pushbutton_Materials,'Enable','on');
    set(handles.pushbutton_Sections,'Enable','on');
    set(handles.pushbutton_Nodes,'Enable','on');
    set(handles.pushbutton_Elements,'Enable','on');
    set(handles.pushbutton_NodalLoads,'Enable','on');
    set(handles.pushbutton_ElementLoads,'Enable','on');
    set(handles.pushbutton_Supports,'Enable','on');
end

% Disable set axes limits pushbutton
set(handles.pushbutton_AxesLimits,'Enable','off')

% Reset canvas to 2D
axes(handles.axes_Canvas);
cla reset

% Reset canvas position
dfltUnits = get(handles.axes_Canvas, 'Units');
set(handles.axes_Canvas, 'Units', 'normalized');
set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
set(handles.axes_Canvas, 'Units', dfltUnits);

% Reset mouse properties
prevMouse = getappdata(0,'mouse');
mouse = Emouse_2D(gcf,handles.axes_Canvas);
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.originalData = prevMouse.originalData;
clear prevMouse
setappdata(0,'mouse',mouse)

% Reset draw object and redraw model
draw = Draw_Grillage2D(getappdata(0,'model'));
setappdata(0,'draw',draw);
redraw(handles);

% Reset visualization options
set(handles.plane3D,'Checked','off')
set(handles.planeXY,'Checked','on')
set(handles.planeXZ,'Checked','off','enable','off')
set(handles.planeYZ,'Checked','off','enable','off')
set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')

% Update zoom variables
mouse = getappdata(0,'mouse');
mouse.originalXLim = get(handles.axes_Canvas,'XLim');
mouse.originalYLim = get(handles.axes_Canvas,'YLim');
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function togglebutton_2DView_OffCallback(~, ~, handles) %#ok<DEFNU>
anm = get(handles.popupmenu_Anm,'value');
if anm == 3 && isempty(getappdata(0,'reset3DViewFlag'))
    unselectEntities(handles)
    
    % Enable set axes limits pushbutton
    set(handles.pushbutton_AxesLimits,'Enable','on')
    
    % Reset modeling options
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    % Reset canvas to 3D
    axes(handles.axes_Canvas);
    cla reset
    
    % Avoid showing ruler while canvas is being modified
    handles.axes_Canvas.XAxis.Visible = 'off';
    handles.axes_Canvas.YAxis.Visible = 'off';
    handles.axes_Canvas.ZAxis.Visible = 'off';
    
    % Reset canvas position
    dfltUnits = get(handles.axes_Canvas, 'Units');
    set(handles.axes_Canvas, 'Units', 'normalized');
    set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
    set(handles.axes_Canvas, 'Units', dfltUnits);
    view(3)
    
    % Reset mouse properties
    prevMouse = getappdata(0,'mouse');
    mouse = Emouse_3D(gcf,handles.axes_Canvas);
    mouse.drawNode = get(handles.togglebutton_Node,'state');
    mouse.drawElem = get(handles.togglebutton_Element,'state');
    mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
    mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
    mouse.polyline = get(handles.togglebutton_Polyline,'state');
    mouse.ortho = get(handles.togglebutton_Ortho,'state');
    mouse.originalData = prevMouse.originalData;
    clear prevMouse
    setappdata(0,'mouse',mouse)
    
    % Reset draw object and redraw model
    draw = Draw_Grillage(getappdata(0,'model'));
    setappdata(0,'draw',draw)
    redraw(handles)
    axis equal
    draw.setLimits();
    
    % Reset visualization options
    set(handles.plane3D,'Checked','on')
    set(handles.planeXY,'Checked','off')
    set(handles.planeXZ,'Checked','off','enable','on')
    set(handles.planeYZ,'Checked','off','enable','on')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
    
    % Update zoom variables
    mouse = getappdata(0,'mouse');
    mouse.originalXLim = get(handles.axes_Canvas,'XLim');
    mouse.originalYLim = get(handles.axes_Canvas,'YLim');
    mouse.originalZLim = get(handles.axes_Canvas,'ZLim');
    setappdata(0,'mouse',mouse)
    
elseif ~isempty(getappdata(0,'reset3DViewFlag'))
    rmappdata(0,'reset3DViewFlag');
end

%--------------------------------------------------------------------------
function pushbutton_AxesLimits_ClickedCallback(hObject, eventdata, handles) %#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);

setappdata(0,'move',1);
GUI_SetAxesLimits

%--------------------------------------------------------------------------
function togglebutton_Node_OnCallback(hObject, eventdata, handles)
% Get analysis model
anm = get(handles.popupmenu_Anm,'value');
% if anm == 3 && strcmp(get(handles.togglebutton_2DView,'state'),'off')
%     set(handles.togglebutton_2DView,'state','on')
% end

toggleButtonsOff(hObject,handles)

if anm < 4
    set(handles.togglebutton_SnapToGrid,'enable','on')
end

if get(handles.popupmenu_Results,'value') ~= 1
    set(handles.popupmenu_Results,'value',1);
    popupmenu_Results_Callback(handles.popupmenu_Results, [], handles);
end

unselectEntities(handles);

% Set mouse properties
mouse = getappdata(0,'mouse');
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

% --------------------------------------------------------------------
function togglebutton_Node_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
if strcmp(get(handles.togglebutton_Element,'state'),'off')
    set(handles.togglebutton_SnapToGrid,'state','off','enable','off')
end

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

unselectEntities(handles);

% Set mouse properties
mouse = getappdata(0,'mouse');
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

% --------------------------------------------------------------------
function togglebutton_Element_OnCallback(hObject, eventdata, handles)
% Get analysis model
anm = get(handles.popupmenu_Anm,'value');

% get number of materials and cross-sections
nmat = getappdata(0,'nmat');
nsec = getappdata(0,'nsec');

if nmat == 0
    msgbox('No material has been created.','Error','error');
    set(hObject,'state','off')
    return
elseif nsec == 0
    msgbox('No cross-section has been created.','Error','error');
    set(hObject,'state','off')
    return
end

toggleButtonsOff(hObject,handles)

if anm < 4
    set(handles.togglebutton_SnapToGrid,'enable','on')
    if strcmp(get(handles.togglebutton_SnapToGrid,'state'),'off')
        set(handles.togglebutton_Ortho,'enable','on')
    end
end
set(handles.togglebutton_CrossElements,'enable','on')
set(handles.togglebutton_Polyline,'enable','on')

unselectEntities(handles);

if get(handles.popupmenu_Results,'value') ~= 1
    set(handles.popupmenu_Results,'value',1);
    popupmenu_Results_Callback(handles.popupmenu_Results, [], handles);
end

% Set mouse properties
mouse = getappdata(0,'mouse');
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

% --------------------------------------------------------------------
function togglebutton_Element_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
if strcmp(get(handles.togglebutton_Node,'state'),'off')
    set(handles.togglebutton_SnapToGrid,'state','off','enable','off')
end
set(handles.togglebutton_CrossElements,'enable','off')
set(handles.togglebutton_Polyline,'enable','off')
set(handles.togglebutton_Ortho,'enable','off')

flag = areAllToggleOff(handles);

if flag == 1
    set(handles.togglebutton_Cursor,'state','on')
end

unselectEntities(handles);

% Set mouse properties
mouse = getappdata(0,'mouse');
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.elemNode = 0;
mouse.selectedNode = 0;
mouse.elemNodeID = [];
mouse.elemCoords = [];
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.orthoCoords = [];
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function pushbutton_SolveIntSects_ClickedCallback(hObject, ~, ~) %#ok<DEFNU>
choice = questdlg('A node will be created on every intersection. Continue?','Continue','OK','Cancel','Cancel');
if strcmp(choice,'OK')
    intersections = getappdata(0,'intersections');
    mouse = getappdata(0,'mouse');
    for nis = 1:size(intersections,2)
        if nis >= 2
            intersections = getappdata(0,'intersections');
        end
        [~] = auxMouseFctn('drawNodes',mouse,{intersections(1).coord,intersections(1).elems});
    end
    setappdata(0,'intersections',[])

    set(hObject,'enable','off')
end

%--------------------------------------------------------------------------
function pushbutton_DeleteEntities_ClickedCallback(hObject, eventdata, handles)
mouse = getappdata(0,'mouse');

if strcmp(get(handles.togglebutton_Element,'state'),'on')
    mouse.elemNode = 0;
    mouse.elemNodeID = [];
    mouse.elemCoords = [];
    mouse.moveAction();
end

if mouse.selectedNode ~= 0
  n = mouse.selectedNode;
  [~] = auxModelFctn('deleteNodes',n);
  unselectEntities(handles)
elseif mouse.selectedElem ~= 0
  e = mouse.selectedElem;
  [~] = auxModelFctn('deleteElems',e);
  unselectEntities(handles)
end

% --------------------------------------------------------------------
function togglebutton_SnapToGrid_OnCallback(hObject, eventdata, handles) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
dfltAns = num2str(mouse.snapPrecision);
snapPrecision = char(inputdlg('Enter grid spacing','Grid',[1 50],{dfltAns}));
if ~isnan(str2double(snapPrecision)) && str2double(snapPrecision) > 0
    mouse.snapPrecision = str2double(snapPrecision);
elseif isempty(snapPrecision)
    set(hObject,'state','off')
    return
end

unselectEntities(handles);
set(handles.togglebutton_Ortho,'state','off','enable','off')

% Set mouse properties
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

% --------------------------------------------------------------------
function togglebutton_SnapToGrid_OffCallback(hObject, eventdata, handles) %#ok<DEFNU>
unselectEntities(handles);
if strcmp(get(handles.togglebutton_Element,'state'),'on')
    set(handles.togglebutton_Ortho,'enable','on')
end

% Set mouse properties
mouse = getappdata(0,'mouse');
mouse.drawNode = get(handles.togglebutton_Node,'state');
mouse.drawElem = get(handles.togglebutton_Element,'state');
mouse.snapToGrid = get(handles.togglebutton_SnapToGrid,'state');
mouse.polyline = get(handles.togglebutton_Polyline,'state');
mouse.ortho = get(handles.togglebutton_Ortho,'state');
mouse.intersectElem = get(handles.togglebutton_CrossElements,'state');
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function togglebutton_CrossElements_OnCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.intersectElem = 'on';
setappdata(0,'mouse',mouse);

%--------------------------------------------------------------------------
function togglebutton_CrossElements_OffCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.intersectElem = 'off';
setappdata(0,'mouse',mouse);

%--------------------------------------------------------------------------
function togglebutton_Polyline_OnCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.polyline = 'on';
setappdata(0,'mouse',mouse);

%--------------------------------------------------------------------------
function togglebutton_Polyline_OffCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.polyline = 'off';
setappdata(0,'mouse',mouse);

% --------------------------------------------------------------------
function togglebutton_Ortho_OnCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.ortho = 'on';
setappdata(0,'mouse',mouse);

% --------------------------------------------------------------------
function togglebutton_Ortho_OffCallback(~, ~, ~) %#ok<DEFNU>
mouse = getappdata(0,'mouse');
mouse.ortho = 'off';
setappdata(0,'mouse',mouse);

%--------------------------------------------------------------------------
% Executes during "popupmenu_AnalysisType" creation.
function popupmenu_AnalysisType_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes on selection change in "Analysis Type" popupmenu.
% Selects the intended analysis type. (Static, Dynamic)
function popupmenu_AnalysisType_Callback(hObject, ~, handles) %#ok<DEFNU>
include_constants;

% Reset modeling/result options
setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);

% Get popupmenu selected value
val = get(hObject,'value');

% Get handle to model object
model = getappdata(0,'model');

% Check if selection is the current analysis type
if (val == LE_STATIC_ANALYSIS  && model.whichSolver == STATIC_LINEAR) ||...
   (val == LE_DYNAMIC_ANALYSIS &&...
    (model.whichSolver == DYNAMIC_NEWMARK_LINEAR ||...
     model.whichSolver == DYNAMIC_MODALSUP_LINEAR))
    return
end

% Disable popupmenu while tasks are being performed
set(hObject,'enable','off');

switch val
    case LE_STATIC_ANALYSIS
        % Disable dynamic parameters and results buttons
        %set(handles.pushbutton_DynamicParameters,'Enable','off');
        set(handles.pushbutton_DynamicResults,'Enable','off','Visible','off');
        set(handles.pushbutton_PlayDynamicResults,'Enable','off','Visible','off');
        
        % Reset solver as static
        model.whichSolver = STATIC_LINEAR;
        model.drv = [];
        model.drv = Drv_LES(true,model);
        setappdata(0,'model',model);
              
        % Rearrange buttons positions on the modeling panel
        set(handles.pushbutton_DynamicAnalysis,'Enable','off');
        set(handles.pushbutton_TimeFcn,'Enable','off');
        set(handles.dynResOptButton,'Enable','off');
        set(handles.pushbutton_EditLoadCase,'Enable','on');
        set(handles.popupmenu_LoadCase,'Enable','on');
        
        % Update load case popupmenu
        nclc = getappdata(0,'currentLc');
        if model.ncomb ~=0
            set(handles.popupmenu_LoadCase,'String',{char(model.strLc),char(model.strComb)},'Value',nclc,'Max',model.nlc+model.ncomb)
        else
            set(handles.popupmenu_LoadCase,'String',char(model.strLc),'Value',nclc,'Max',model.nlc)
        end
        
        % Check if selection is a load case or a load combination
        if nclc <= model.nlc
            set(handles.pushbutton_ElementLoads,'Enable','on');
        else
            set(handles.pushbutton_ElementLoads,'Enable','off');
        end
        
        % Reset graphical objects and update model drawing
        if strcmp(get(handles.popupmenu_Anm,'Enable'),'off')
            set(handles.pushbutton_ProcessData,'Enable','on');
        end
       
        set(handles.popupmenu_Results,'Enable','off','value',1);
        set(handles.pushbutton_Textual,'Enable','off');
        set(handles.checkbox_Reactions,'Enable','off', 'Value', 0);
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(handles.edit_Scale,'enable','off','visible','off');
        redraw(handles,'Loads')
        
    case LE_DYNAMIC_ANALYSIS
        % Enable dynamic parameters and results buttons
        %set(handles.pushbutton_DynamicParameters,'Enable','on');
        set(handles.pushbutton_DynamicResults,'Enable','off','Visible','on');
        set(handles.pushbutton_PlayDynamicResults,'Enable','off','Visible','on');
        
        % Reset solver as dynamic (numerical)
        model.whichSolver = DYNAMIC_NEWMARK_LINEAR;
        model.drv = [];
        model.drv = Drv_LED(DYNAMIC_NEWMARK_LINEAR,true,model);
        setappdata(0,'model',model);
                
        % Rearrange buttons positions on the modeling panel
        set(handles.pushbutton_DynamicAnalysis,'Enable','on');
        set(handles.dynResOptButton,'Enable','on');
        set(handles.pushbutton_ElementLoads,'Enable','off');
        set(handles.popupmenu_LoadCase,'Enable','off','Value',1,'String',' ');
        set(handles.pushbutton_EditLoadCase,'Enable','off');
        
        % Check if selection is a load case or a load combination
        if getappdata(0,'currentLc') <= model.nlc
            set(handles.pushbutton_TimeFcn,'Enable','on');
        else
            set(handles.pushbutton_TimeFcn,'Enable','off');
        end
        
        % Reset graphical objects and update model drawing
        if strcmp(get(handles.popupmenu_Anm,'Enable'),'off')
            set(handles.pushbutton_ProcessData,'Enable','on');
        end
        
        set(handles.popupmenu_Results,'Enable','off','value',1);
        set(handles.pushbutton_Textual,'Enable','off');
        set(handles.checkbox_Reactions,'Enable','off', 'Value', 0);
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(handles.edit_Scale,'enable','off','visible','off');
        redraw(handles,'Loads');
end
% Enable popupmenu for future use
set(hObject,'enable','on');

%--------------------------------------------------------------------------
% Executes during "popupmenu_Anm" creation.
function popupmenu_Anm_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes on selection change in "Analysis Model" popupmenu.
% Selects the intended analysis model.
function popupmenu_Anm_Callback(hObject, eventdata, handles) %#ok<DEFNU>
model = getappdata(0,'model');

% Get current analysis model option
anm = get(hObject,'value');

% Check if analysis model has been altered. If not, do nothing
if anm == model.anm.analysis_type + 1
    return
end

% Turn off 2D view togglebutton
if strcmp(get(handles.togglebutton_2DView,'state'),'on')
    set(handles.togglebutton_2DView,'state','off')
end

% Clean canvas
axes(handles.axes_Canvas);
cla reset

toggleButtonsOff([],handles)

% Set properties that depend on analysis model
if anm == 1
    anm = Anm_Truss2D();
    print = Print_Truss2D();
    draw = Draw_Truss2D();
    view(2);
    axis equal
    xlim([0,10])
    set(handles.planeButton,'Enable','off')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')
    set(handles.togglebutton_Node,'Enable','on')
    set(handles.togglebutton_Element,'Enable','on')
    set(handles.togglebutton_2DView,'enable','off')
    set(handles.pushbutton_AxesLimits,'enable','off')
elseif anm == 2
    anm = Anm_Frame2D();
    print = Print_Frame2D();
    draw = Draw_Frame2D();
    view(2);
    axis equal
    xlim([0,10])
    set(handles.planeButton,'Enable','off')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','off')
    set(handles.togglebutton_Node,'Enable','on')
    set(handles.togglebutton_Element,'Enable','on')
    set(handles.togglebutton_2DView,'enable','off')
    set(handles.pushbutton_AxesLimits,'enable','off')
elseif anm == 3
    anm = Anm_Grillage();
    print = Print_Grillage();
    draw = Draw_Grillage();
    view(3);
    axis equal
    xlim([-5,5])
    ylim([-5,5])
    zlim([0,0.02])
    xlabel('X');
    ylabel('Y');
    zlabel(' ');
    zticks(0);
    zticklabels({' '});
    set(handles.planeButton,'Enable','on')
    set(handles.plane3D,'Checked','on')
    set(handles.planeXY,'Checked','off')
    set(handles.planeXZ,'Checked','off')
    set(handles.planeYZ,'Checked','off')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
    set(handles.togglebutton_Node,'Enable','on')
    set(handles.togglebutton_Element,'Enable','on')
    set(handles.togglebutton_2DView,'enable','on')
    set(handles.pushbutton_AxesLimits,'enable','on')
elseif anm == 4
    anm = Anm_Truss3D();
    print = Print_Truss3D();
    draw = Draw_Truss3D();
    view(3);
    axis equal
    xlim([0,10])
    ylim([0,10])
    zlim([-1,1])
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(handles.planeButton,'Enable','on')
    set(handles.plane3D,'Checked','on')
    set(handles.planeXY,'Checked','off')
    set(handles.planeXZ,'Checked','off')
    set(handles.planeYZ,'Checked','off')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
    set(handles.togglebutton_Node,'Enable','on')
    set(handles.togglebutton_Element,'Enable','on')
    set(handles.togglebutton_2DView,'enable','off')
    set(handles.pushbutton_AxesLimits,'enable','off')
elseif anm == 5
    anm = Anm_Frame3D();
    print = Print_Frame3D();
    draw = Draw_Frame3D();
    view(3);
    axis equal
    xlim([0,10])
    ylim([0,10])
    zlim([-1,1])
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(handles.planeButton,'Enable','on')
    set(handles.plane3D,'Checked','on')
    set(handles.planeXY,'Checked','off')
    set(handles.planeXZ,'Checked','off')
    set(handles.planeYZ,'Checked','off')
    set(findall(findall(gcf),'Tag','toolbar_rotate3d'),'Enable','on')
    set(handles.togglebutton_Node,'Enable','on')
    set(handles.togglebutton_Element,'Enable','on')
    set(handles.togglebutton_2DView,'enable','off')
    set(handles.pushbutton_AxesLimits,'enable','off')
end

% Deactivate visualizing tools
pan off
rotate3d off
zoom off

% Turn grid on/off
if strcmp(get(handles.gridButton,'Checked'),'on') == 1
    grid on
else
    grid off
end

% Turn ruler on/off
if strcmp(get(handles.rulerButton,'Checked'),'on') == 1
    handles.axes_Canvas.XAxis.Visible = 'on';
    handles.axes_Canvas.YAxis.Visible = 'on';
    handles.axes_Canvas.ZAxis.Visible = 'on';
else
    handles.axes_Canvas.XAxis.Visible = 'off';
    handles.axes_Canvas.YAxis.Visible = 'off';
    handles.axes_Canvas.ZAxis.Visible = 'off';
end

% Delete cross-sections and set default
model = getappdata(0,'model');
nsec = getappdata(0,'nsec');
if nsec ~= 0
    initialSec = Section(1,0.01,0.01,0.01,0.00001,0.00001,0.00001,0.1,0.1);
    model.nsec = 1;
    model.sections = initialSec;
    setappdata(0,'nsec',1);
    setappdata(0,'sections',initialSec);
end

% Update information panel
infoPanelData = get(handles.uitable_infoPanel,'Data');
infoPanelData(2,:) = {'Cross-Sections',1};
set(handles.uitable_infoPanel,'Data',infoPanelData)

% Update objects properties
model.anm = anm;
print.model = model;
draw.mdl = model;

% Adjust canvas position
dfltUnits = get(handles.axes_Canvas, 'Units');
set(handles.axes_Canvas, 'Units', 'normalized');
if (model.anm.analysis_type == 0) || (model.anm.analysis_type == 1)
    set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
else
    set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
end
set(handles.axes_Canvas, 'Units', dfltUnits);

% Reinitialize object for mouse events and save it in root
if (model.anm.analysis_type == 0) || (model.anm.analysis_type == 1)
    mouse = Emouse_2D(gcf,handles.axes_Canvas);
else
    mouse = Emouse_3D(gcf,handles.axes_Canvas);
end

% Return objects to root
setappdata(0,'mouse',mouse);
setappdata(0,'model',model);
setappdata(0,'print',print);
setappdata(0,'draw',draw);

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_EditLoadCase.
function pushbutton_EditLoadCase_Callback(~,~,handles)
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

setappdata(0,'move',1);
GUI_EditLoadCase

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function popupmenu_LoadCase_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes on selection change in popupmenu_LoadCase.
function popupmenu_LoadCase_Callback(hObject, eventdata, handles)
model = getappdata(0,'model');
nodes = getappdata(0,'nodes');
elems = getappdata(0,'elems');

% Get current load case (or combination)
lc = get(hObject,'Value');

if lc == getappdata(0,'currentLc')
    return
else
    setappdata(0,'currentLc',lc)
end    

unselectEntities(handles);
toggleButtonsOff([],handles);

% Get number of load cases
nlc = model.nlc;

% Check if selection is a load case or a load combination
if lc <= nlc
    
    set(handles.popupmenu_AnalysisType,'Enable','on');
    set(handles.togglebutton_Node,'enable','on')
    set(handles.togglebutton_Element,'enable','on')
    set(handles.togglebutton_SnapToGrid,'enable','off')
    
    % Allocate nodal loads and prescribed displacements in each Node object
    % obs.: nodalLoadCase(:,i) = [ fx        
    %                              fy
    %                              fz     (nodalLoadCase is a matrix, each
    %                              mx      column refers to one specific 
    %                              my      load case)
    %                              mz
    %                              dx
    %                              dy
    %                              dz
    %                              rx
    %                              ry
    %                              rz ]
    for n = 1:model.nnp
        if isempty(nodes(n).nodalLoadCase) == 0 && lc <= size(nodes(n).nodalLoadCase,2)
            if all(nodes(n).nodalLoadCase(1:6,lc) == 0)
                nodes(n).load.static = []; 
            else %if there are any nodal loads, set them to load.static
                nodes(n).load.static = nodes(n).nodalLoadCase(1:6,lc);
            end
            if size(nodes(n).nodalLoadCase,1) <= 6
                nodes(n).prescDispl = [];
            elseif all(nodes(n).nodalLoadCase(7:12,lc) == 0)
                nodes(n).prescDispl = [];
            else %if there are any prescribed displacements, set them to prescDispl
                nodes(n).prescDispl = nodes(n).nodalLoadCase(7:12,lc);
            end
        else
            nodes(n).load.static = []; 
            nodes(n).prescDispl = [];
        end
    end    
    
    % Allocate element loads (distributed loads and thermal loads) in each
    % Lelem object
    % obs.: elemLoadCase(:,i) = [  unifDir
    %                                qx
    %                                qy
    %                                qz        (elemLoadCase is a matrix, 
    %                             linearDir     each column refers to one  
    %                                qx1        specific load case)
    %                                qy1
    %                                qz1
    %                                qx2
    %                                qy2
    %                                qz2
    %                                dtx
    %                                dty
    %                                dtz   ]
    for e = 1:model.nel
        if isempty(elems(e).load.elemLoadCase) == 0 && lc <= size(elems(e).load.elemLoadCase,2)
            % Clear previous uniform loads
            elems(e).load.uniformGbl = [];    
            elems(e).load.uniformLcl = [];
            if all(elems(e).load.elemLoadCase(2:4,lc) == 0)
                elems(e).load.uniformDir = 0;  
            else  % if there are uniform loads, set their local and global components
                elems(e).load.uniformDir = elems(e).load.elemLoadCase(1,lc);
                elems(e).load.setUnifLoad((elems(e).load.elemLoadCase(2:4,lc))',elems(e).load.elemLoadCase(1,lc));
            end
            % Clear previous linear loads
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            if size(elems(e).load.elemLoadCase,1) <= 5
                elems(e).load.linearDir = 0;
            elseif all(elems(e).load.elemLoadCase(6:11,lc) == 0)
                elems(e).load.linearDir = 0;
            else  % if there are linear loads, set their local and global components
                elems(e).load.linearDir = elems(e).load.elemLoadCase(5,lc);
                elems(e).load.setLinearLoad((elems(e).load.elemLoadCase(6:11,lc))',elems(e).load.elemLoadCase(5,lc));
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
            end
            % Set thermal loads
            if size(elems(e).load.elemLoadCase,1) <= 11
                elems(e).load.tempVar_X = 0;
                elems(e).load.tempVar_Y = 0;  
                elems(e).load.tempVar_Z = 0;
            else
                elems(e).load.tempVar_X = elems(e).load.elemLoadCase(12,lc);
                elems(e).load.tempVar_Y = elems(e).load.elemLoadCase(13,lc);
                elems(e).load.tempVar_Z = elems(e).load.elemLoadCase(14,lc);
            end
        else
            elems(e).load.uniformDir = 0;      
            elems(e).load.uniformGbl = [];    
            elems(e).load.uniformLcl = [];    
            elems(e).load.linearDir = 0;       
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            elems(e).load.tempVar_X = 0;
            elems(e).load.tempVar_Y = 0;  
            elems(e).load.tempVar_Z = 0; 
        end
    end    

    % Enable process data button and disable results
    if model.nnp ~= 0
        set(handles.pushbutton_ProcessData,'Enable','on');
    else
        set(handles.pushbutton_ProcessData,'Enable','off');
    end
    set(handles.popupmenu_Results,'Enable','off','Value',1);
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.pushbutton_DynamicResults,'enable','off');
    set(handles.pushbutton_PlayDynamicResults,'Enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    set(handles.pushbutton_Materials,'Enable','on');
    set(handles.pushbutton_Sections,'Enable','on');
    set(handles.pushbutton_Nodes,'Enable','on');
    set(handles.pushbutton_Elements,'Enable','on');
    set(handles.pushbutton_NodalLoads,'Enable','on');
    set(handles.pushbutton_ElementLoads,'Enable','on');
    set(handles.pushbutton_Supports,'Enable','on');
    
else  % if lc > nlc, the current load case is a combination
    
    set(handles.popupmenu_AnalysisType,'Enable','off');
    set(handles.togglebutton_Node,'enable','off')
    set(handles.togglebutton_Element,'enable','off')
    set(handles.togglebutton_SnapToGrid,'enable','off')
    
    % Consider load cases factors for the selected combination and allocate
    % the resulting nodal loads and prescribed displacements in each Node
    % object.
    % obs.: nodalLoadCase(:,i) = [ fx        
    %                              fy
    %                              fz     (nodalLoadCase is a matrix, each
    %                              mx      column refers to a specific load
    %                              my      case)
    %                              mz
    %                              dx
    %                              dy
    %                              dz
    %                              rx
    %                              ry
    %                              rz ]
    %
    % obs.2: loadComb is a matrix that holds combination factors for load
    % case combinations, where each column refers to a specific load comb.
    
    for n = 1:model.nnp
        if isempty(nodes(n).nodalLoadCase) == 0
            allLogic = all(nodes(n).nodalLoadCase(1:6,:) == 0);
            if all(allLogic == 1)
                nodes(n).load.static = []; 
            else  % if there are nodal loads, ,multiply them by their 
                  % proper comb factors, sum the results and set them to 
                  % load.static.
                nodes(n).load.static(1) = nodes(n).nodalLoadCase(1,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).load.static(2) = nodes(n).nodalLoadCase(2,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).load.static(3) = nodes(n).nodalLoadCase(3,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).load.static(4) = nodes(n).nodalLoadCase(4,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).load.static(5) = nodes(n).nodalLoadCase(5,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).load.static(6) = nodes(n).nodalLoadCase(6,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                if all(nodes(n).load.static == 0)
                    nodes(n).load.static = [];
                end    
            end
            if size(nodes(n).nodalLoadCase,1) <= 6
                allLogic = 1;
            else
                allLogic = all(nodes(n).nodalLoadCase(7:12,:) == 0);
            end
            if all(allLogic == 1)
                nodes(n).prescDispl = [];
            else  % if there are prescribed displacements, ,multiply them 
                  % by their proper comb factors, sum the results and set  
                  % them to prescDispl.
                nodes(n).prescDispl(1) = nodes(n).nodalLoadCase(7,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).prescDispl(2) = nodes(n).nodalLoadCase(8,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).prescDispl(3) = nodes(n).nodalLoadCase(9,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).prescDispl(4) = nodes(n).nodalLoadCase(10,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).prescDispl(5) = nodes(n).nodalLoadCase(11,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                nodes(n).prescDispl(6) = nodes(n).nodalLoadCase(12,:) * model.loadComb(1:size(nodes(n).nodalLoadCase,2),lc-nlc);
                if all(nodes(n).prescDispl == 0)
                    nodes(n).prescDispl = [];
                end
            end
        else
            nodes(n).load.static = []; 
            nodes(n).prescDispl = [];
        end
    end    

    % Consider load cases factors for the selected combination and allocate
    % the resulting element loads and prescribed displacements in each
    % Lelem object.
    % obs.: elemLoadCase(:,i) = [  unifDir
    %                                qx
    %                                qy
    %                                qz        (elemLoadCase is a matrix, 
    %                             linearDir     each column refers to a  
    %                                qx1        specific load case)
    %                                qy1
    %                                qz1
    %                                qx2
    %                                qy2
    %                                qz2
    %                                dtx
    %                                dty
    %                                dtz   ]
    %
    % obs.2: loadComb is a matrix that holds combination factors for load
    % case combinations, where each column refers to a specific load comb.
    
    for e = 1:model.nel
        if isempty(elems(e).load.elemLoadCase) == 0
            % Clear previous uniform loads
            elems(e).load.uniformGbl = [];    
            elems(e).load.uniformLcl = [];
            allLogic = all(elems(e).load.elemLoadCase(2:4,:) == 0);
            if all(allLogic == 1)
                elems(e).load.uniformDir = 0;
            else  % if there are uniform loads, set their local and global
                  % components and consider their comb factors
                unifLoad = zeros(size(elems(e).load.elemLoadCase,2),3);
                for i = 1:size(elems(e).load.elemLoadCase,2)
                    unifLoad(i,:) = model.loadComb(i,lc-nlc) * elems(e).load.elemLoadCase(2:4,i);
                    if ~all(unifLoad(i,:) == 0)
                        elems(e).load.setUnifLoad(unifLoad(i,:),elems(e).load.elemLoadCase(1,i));
                    end    
                end
            end
            
            % Clear previous linear loads
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            if size(elems(e).load.elemLoadCase,1) <= 5
                allLogic = 1;
            else
                allLogic = all(elems(e).load.elemLoadCase(6:11,:) == 0);
            end
            if all(allLogic == 1)
                elems(e).load.linearDir = 0;
            else  % if there are linear loads, set their local and global
                  % components and consider their comb factors    
                linLoad = zeros(size(elems(e).load.elemLoadCase,2),6);
                for i = 1:size(elems(e).load.elemLoadCase,2)
                    linLoad(i,:) = model.loadComb(i,lc-nlc) * elems(e).load.elemLoadCase(6:11,i);
                    if ~all(linLoad(i,:) == 0)
                        elems(e).load.setLinearLoad(linLoad(i,:),elems(e).load.elemLoadCase(5,i));
                    end    
                end
                % Avoids numeric problems with new loads
                for i = 1:size(elems(e).load.linearLcl,1)
                    if abs(elems(e).load.linearLcl(i)) <= 10^-10
                        elems(e).load.linearLcl(i) = 0;
                    end
                     if abs(elems(e).load.linearGbl(i)) <= 10^-10
                        elems(e).load.linearGbl(i) = 0;
                    end
                end
            end
            % Set thermal loads
            if size(elems(e).load.elemLoadCase,1) <= 11
                elems(e).load.tempVar_X = 0;
                elems(e).load.tempVar_Y = 0;  
                elems(e).load.tempVar_Z = 0; 
            else
                elems(e).load.tempVar_X = elems(e).load.elemLoadCase(12,:) * model.loadComb(1:size(elems(e).load.elemLoadCase,2),lc-nlc);
                elems(e).load.tempVar_Y = elems(e).load.elemLoadCase(13,:) * model.loadComb(1:size(elems(e).load.elemLoadCase,2),lc-nlc);  
                elems(e).load.tempVar_Z = elems(e).load.elemLoadCase(14,:) * model.loadComb(1:size(elems(e).load.elemLoadCase,2),lc-nlc);
            end
        else
            elems(e).load.uniformDir = 0;      
            elems(e).load.uniformGbl = [];    
            elems(e).load.uniformLcl = [];    
            elems(e).load.linearDir = 0;       
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            elems(e).load.tempVar_X = 0;
            elems(e).load.tempVar_Y = 0;  
            elems(e).load.tempVar_Z = 0; 
        end
    end    

    % Enable process data button and disable results
    if model.nnp ~= 0
        set(handles.pushbutton_ProcessData,'Enable','on');
    else
        set(handles.pushbutton_ProcessData,'Enable','off');
    end
    set(handles.popupmenu_Results,'Enable','off','Value',1);
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.pushbutton_DynamicResults,'enable','off');
    set(handles.pushbutton_PlayDynamicResults,'Enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    set(handles.pushbutton_Materials,'Enable','off','Value',0);
    set(handles.pushbutton_Sections,'Enable','off','Value',0);
    set(handles.pushbutton_Nodes,'Enable','off','Value',0);
    set(handles.pushbutton_Elements,'Enable','off','Value',0);
    set(handles.pushbutton_NodalLoads,'Enable','off','Value',0);
    set(handles.pushbutton_ElementLoads,'Enable','off','Value',0);
    set(handles.pushbutton_Supports,'Enable','off','Value',0);
end    
    
model.nodes = nodes;
model.elems = elems;
setappdata(0,'model',model);
setappdata(0,'nodes',nodes);
setappdata(0,'elems',elems);    

% Draw updated model
set(handles.popupmenu_Results,'Value',1)
redraw(handles,'Loads');
delete(findobj('tag','drawReactions'));
delete(findobj('tag','textForceReactions'));
delete(findobj('tag','textMomentReactions'));

%--------------------------------------------------------------------------
% Executes on button press in "Dynamic Analysis Parameters" pushbutton.
function pushbutton_DynamicAnalysis_Callback(~, ~, handles)%#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

setappdata(0,'move',1);
GUI_Dynamic

%--------------------------------------------------------------------------
% Executes on button press in "Materials" pushbutton.
% Opens a new window to create new materials or delete existing ones.
function pushbutton_Materials_Callback(~,~,handles)
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

setappdata(0,'move',1);
GUI_Materials

%--------------------------------------------------------------------------
% Executes on button press in "Sections" pushbutton.
% Opens a new window to create new cross-sections or delete existing ones.
function pushbutton_Sections_Callback(~,~,handles)
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

setappdata(0,'move',1);
GUI_Sections

%--------------------------------------------------------------------------
% Executes on button press in "Nodes" pushbutton.
% Opens a new window to create new nodal points or delete existing ones.
function pushbutton_Nodes_Callback(~,~,handles)%#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

setappdata(0,'move',1);
GUI_Nodes

%--------------------------------------------------------------------------
% Executes on button press in "Elements" pushbutton.
% Opens a new window to create new elements or delete existing ones.
function pushbutton_Elements_Callback(~,~,handles)%#ok<DEFNU>
setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

if getappdata(0,'nmat') == 0
    msgbox('No material has been created.', 'Error','error');
elseif getappdata(0,'nsec') == 0
    msgbox('No cross-section has been created.', 'Error','error');
elseif getappdata(0,'nnp') < 2
    msgbox('At least two nodes must be created.', 'Error','error');
else
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    unselectEntities(handles);
    
    GUI_Elements
end

%--------------------------------------------------------------------------
% Executes on button press in "Supports" pushbutton.
% Opens a new window to set support conditions and precr. displ. of nodes.
function pushbutton_Supports_Callback(~,~,handles)
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

if getappdata(0,'nnp') == 0
    msgbox('No node has been created.', 'Error','error');
else
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    unselectEntities(handles);
    
    GUI_Supports
end

%--------------------------------------------------------------------------
% Executes on button press in "Nodal Loads" pushbutton.
% Opens a new window to set nodal load values.
function pushbutton_NodalLoads_Callback(~,~,handles)
include_constants;

setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

% Get popupmenu selected value
val = get(handles.popupmenu_AnalysisType,'value');

if getappdata(0,'nnp') == 0
    msgbox('No node has been created.', 'Error','error');
    
elseif val == LE_STATIC_ANALYSIS
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    unselectEntities(handles);
    
    GUI_NodalLoads
else
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    unselectEntities(handles);
    
%     model = getappdata(0,'model');
%     if model.t <= 0 || model.n_steps <= 0
%         msgbox('Number of time steps and time interval for dynamic analysis must be set prior to setting dynamic loads.', 'Attention','warn');
%         return
%     end
    
    GUI_NodalLoads_Dynamic
end

%--------------------------------------------------------------------------
% Executes on button press in "Nodal Loads" pushbutton.
% Opens a new window to set nodal load values.
function pushbutton_TimeFcn_Callback(~,~,handles)
include_constants;

setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

% Get popupmenu selected value
val = get(handles.popupmenu_AnalysisType,'value');

if  val == LE_STATIC_ANALYSIS
%     set(handles.togglebutton_Node,'state','off')
%     set(handles.togglebutton_Element,'state','off')
%     
%     unselectEntities(handles);
%     
%     GUI_NodalLoads
return
else
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    unselectEntities(handles);
    
    model = getappdata(0,'model');
    if model.t <= 0 || model.n_steps <= 0
        msgbox('Total analysis time must be set before creating time functions.','Error','error');
        return
    end
    
    GUI_TimeFunctions
end

%--------------------------------------------------------------------------
% Executes on button press in "Element Loads" pushbutton.
% Opens a new window to set element load values.
function pushbutton_ElementLoads_Callback(~,~,handles)
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

if getappdata(0,'nel') == 0
    msgbox('No element has been created.', 'Error','error');
else
    unselectEntities(handles);
    
    GUI_ElementLoads
end

%--------------------------------------------------------------------------
% Executes on button press in "Model Information" pushbutton.
% Creates a text file with model information.
function pushbutton_ModelInfo_Callback(hObject, eventdata, handles)%#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

% Open file
filterspec = '*.txt';
DialogTitle = 'LESM - Write model information';
DefaultName = 'model_info';
[filename,pathname] = uiputfile(filterspec,DialogTitle,DefaultName);
if ~filename
    return
end
fullname = strcat(pathname,filename);
fid = fopen(fullname,'wt');

% Get handle to print and model objects
print = getappdata(0,'print');
model = getappdata(0,'model');
print.model = model;

% Create file and write results
print.txt = fid;

% Get info about current load case/comb
lc = get(handles.popupmenu_LoadCase,'Value');
strLc = get(handles.popupmenu_LoadCase,'string');
currentLc = strLc(lc,:);

% Print information
print.header();
fprintf( fid, '\n\n\n____________ M O D E L  I N F O R M A T I O N ____________\n' );
print.modelLabel();
[n_srj,n_ns,n_nl,n_pd,n_ul,n_ll,n_tv] = print.modelDescrip(lc,currentLc);
print.material();
print.section();
print.nodalCoords();
print.nodalSupport();
print.spring(n_ns);
print.nodalLoads(n_nl);
print.nodalPrescDisp(n_pd);
print.elements();
print.srjoints(n_srj);
print.unifElementLoads(n_ul);
print.linearElementLoads(n_ll);
print.temperatureVariation(n_tv);

% Close file
fclose(fid);

% Return object to root
setappdata(0,'print',print)

% Show up text file
winopen(fullname);

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function pushbutton_DynamicParameters_ClickedCallback(~, ~, handles) %#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

setappdata(0,'move',1);

GUI_Dynamic

%--------------------------------------------------------------------------
function pushbutton_DynamicResults_ClickedCallback(~, ~, handles) %#ok<DEFNU>
include_constants;
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

model = getappdata(0,'model');

if strcmp(get(handles.pushbutton_ProcessData,'Enable'),'on') || model.nnp < 1
    msgbox('Data must be processed to obtain results.','Attention','warn')
    return
elseif model.drv.whichResponse == MODAL_ANALYSIS
    msgbox('Only modal analysis was performed. No transient response to show.','Attention','warn')
    return
elseif  model.n_steps < 1 || model.t <= 0
    msgbox('No time interval was provided. No transient response to show.','Attention','warn')
    return
end

setappdata(0,'move',1);
GUI_DynamicResults

%--------------------------------------------------------------------------
% Executes on button press in "Process Data" pushbutton.
% Processes structural model information data to get the analysis results.
% This step is needed so the results can be stored to be used in the
% post-processing phase.
% This step is executed when the data is not up to date, according to the
% following criteria:
%  Update data: No change was made to the model after its data have been
%               processed
%  Out-of-date data: Changes were made to the model after its data has been
%                    processed
function pushbutton_ProcessData_Callback(hObject, eventdata, handles)
include_constants;

% Turn off "Process Data" button
set(hObject,'Enable','off');

unselectEntities(handles);

if strcmp(get(handles.togglebutton_2DView,'state'),'on')
    set(handles.togglebutton_2DView,'state','off')
end

% Check if analysis is static or dynamic
if get(handles.popupmenu_AnalysisType,'Value') == 1
    analysis = 'Static';
elseif get(handles.popupmenu_AnalysisType,'Value') == 2
    analysis = 'Dynamic';
else
    set(hObject,'Enable','on');
    return
end

% Get objects from root
model = getappdata(0,'model');
print = getappdata(0,'print');
elems = getappdata(0,'elems');

% Process data
if (model.nnp > 0)
    % Process data
    data = model.drv.process();
else
    msgbox('This model is empty.', 'Error','error');
    data = 0;
end

% Check if data was successfully processed
if data && strcmp(analysis,'Static')
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    % Reset popupmenu_Results string
    resStr(1) = {'Model'};
    resStr(2) = {'Deformation'};
    resStr(3) = {'Axial Force'};
    resStr(4) = {'Torsion Moment'};
    resStr(5) = {'Shear Force Y'};
    resStr(6) = {'Shear Force Z'};
    resStr(7) = {'Bending Moment Y'};
    resStr(8) = {'Bending Moment Z'};
    
    % Set "model" property of Print object
    print.model = model;
    
    % Uncheck "Reactions" option
    set(handles.checkbox_Reactions,'Value',0);
    
    % Enable results according to analysis model
    set(handles.pushbutton_Textual,'Enable','on');
    set(handles.checkbox_Reactions,'Enable','on');
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    if model.anm.analysis_type == 0
        resStr(4:8) = [];
        view(2)
        set(handles.planeButton,'Enable','off');
    elseif model.anm.analysis_type == 1
        resStr(4) = [];
        resStr(5) = [];
        resStr(5) = [];
        view(2)
        set(handles.planeButton,'Enable','off');
    elseif model.anm.analysis_type == 2
        resStr(3) = [];
        resStr(4) = [];
        resStr(6) = [];
        view(3)
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
    elseif model.anm.analysis_type == 3
        resStr(4:8) = [];
        view(3)
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
    elseif model.anm.analysis_type == 4
        view(3)
        set(handles.planeButton,'Enable','on');
        set(handles.plane3D,'Checked','on');
        set(handles.planeXY,'Checked','off');
        set(handles.planeXZ,'Checked','off');
        set(handles.planeYZ,'Checked','off');
    end
    
    % Enable and update popupmenu_Results based on analysis model
    set(handles.popupmenu_Results,'string',resStr,'Value',1,'Enable','on');
    
    % Disable graphs and play buttons for dynamic results
    set(handles.pushbutton_DynamicResults,'Enable','off','Visible','off');
    set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','off','Visible','off');
    
    % Calculate scale factors
    draw = getappdata(0,'draw');
    draw.mdl = model;
    draw.setSize();
    draw.deformScaleFactor();
    draw.axialScaleFactor();
    draw.shearScaleFactor_XY();
    draw.shearScaleFactor_XZ();
    draw.bendingMomentScaleFactor_XY();
    draw.bendingMomentScaleFactor_XZ();
    dsf = getappdata(0,'deform_sf');
    set(handles.edit_Scale,'Visible','off','Enable','off','String',num2str(dsf*get(handles.slider_Scale,'Value')))
    
    % Recenter model
    if model.anm.analysis_type ~= 2 % grillage
        set(handles.popupmenu_Results,'Value',1)
        mouse = getappdata(0,'mouse');
        mouse.doubleClick();
        setappdata(0,'mouse',mouse);
    end
    
    % Adjust canvas position
    dfltUnits = get(handles.axes_Canvas, 'Units');
    set(handles.axes_Canvas, 'Units', 'normalized');
    if (model.anm.analysis_type == 0) || (model.anm.analysis_type == 1)
        set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
    else
        set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
    end
    set(handles.axes_Canvas, 'Units', dfltUnits);
    
elseif data && strcmp(analysis,'Dynamic')
    set(handles.togglebutton_Node,'state','off')
    set(handles.togglebutton_Element,'state','off')
    
    % Reset popupmenu_Results string
    resStr(1) = {'Model'};
    resStr(2) = {'Vibration Modes'};
    resStr(3) = {'Motion'};
    resStr(4) = {'Axial Force'};
    resStr(5) = {'Torsion Moment'};
    resStr(6) = {'Shear Force Y'};
    resStr(7) = {'Shear Force Z'};
    resStr(8) = {'Bending Moment Y'};
    resStr(9) = {'Bending Moment Z'};
    
    % Set "model" property of Print object
    print.model = model;
    
    if model.anm.analysis_type == 0
        resStr(5:9) = [];
    elseif model.anm.analysis_type == 1
        resStr(5) = [];
        resStr(6) = [];
        resStr(6) = [];
    elseif model.anm.analysis_type == 2
        resStr(4) = [];
        resStr(5) = [];
        resStr(7) = [];
    elseif model.anm.analysis_type == 3
        resStr(5:9) = [];
    end
    
    if model.drv.whichResponse ~= MODAL_ANALYSIS && model.n_steps > 0
        % Enable and update popupmenu_Results based on analysis model
        set(handles.popupmenu_Results,'string',resStr,'Value',1,'Enable','on')
    else
        set(handles.popupmenu_Results,'string',resStr(1:2),'Value',1,'Enable','on')
    end
    
    % Calculate scale factors
    if (model.nel > 0)
    draw = getappdata(0,'draw');
    draw.mdl = model;
    draw.setSize();
    if draw.mdl.drv.whichResponse ~= MODAL_ANALYSIS && draw.mdl.n_steps > 0
        draw.dynamicDeformScaleFactor();
        draw.axialScaleFactor();
        draw.torsionScaleFactor();
        draw.shearScaleFactor_XY();
        draw.shearScaleFactor_XZ();
        draw.bendingMomentScaleFactor_XY();
        draw.bendingMomentScaleFactor_XZ();
    end
    end
        
    % Uncheck "Reactions" option
    set(handles.checkbox_Reactions,'Value',0);
    
    % Enable results according to analysis model
    set(handles.pushbutton_Textual,'Enable','on');
    set(handles.checkbox_Reactions,'Enable','off');
    set(handles.pushbutton_DynamicResults,'Enable','on','Visible','on');
    set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','off','Visible','on');
    mouse = getappdata(0,'mouse');
    mouse.doubleClick();
else
    % Turn on "Process Data" button for future use
    set(hObject,'Enable','on');
end

% Return objects to root
setappdata(0,'model',model);
setappdata(0,'print',print);
setappdata(0,'elems',elems);

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Results.
function popupmenu_Results_Callback(hObject, ~, handles)
unselectEntities(handles);
result = get(hObject,'Value');
setappdata(0,'isDrawingDynamic',false)

if get(handles.popupmenu_AnalysisType,'Value') == 1
    if result == 1 % Model
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        redraw(handles,'Loads');
    elseif result == 2 % Deformation
        set(handles.togglebutton_Node,'state','off')
        set(handles.togglebutton_Element,'state','off')
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        redraw(handles);
    else
        set(handles.togglebutton_Node,'state','off')
        set(handles.togglebutton_Element,'state','off')
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','on');
        draw = getappdata(0,'draw');
        draw.axialScaleFactor();
        draw.torsionScaleFactor();
        draw.shearScaleFactor_XY();
        draw.shearScaleFactor_XZ();
        draw.bendingMomentScaleFactor_XY();
        draw.bendingMomentScaleFactor_XZ();
        redraw(handles);
    end
else
    delete(findobj('tag','drawVibrationMode'))
    delete(findobj('tag','drawDynamicDeform'))
    model = getappdata(0,'model');
    if result == 1 % Model
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','off','Visible','on')
        redraw(handles,'Nodal Loads');
    elseif result == 2 % Vibration Modes
        set(handles.text_Element,'string','Mode');
        set(handles.popupmenu_ElementResults,'Visible','on','Enable','on','string',1:model.n_modes,'value',1,'max',model.n_modes);
        set(handles.edit_ElementResults,'Visible','off','Enable','off','string','All');
        set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','on','Visible','on')
        redraw(handles);
    elseif result == 3 % Motion (Dynamic Deformation)
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','on','Visible','on')
        redraw(handles);
    else % Envelopes
        set(handles.text_Element,'string','Elements');
        set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(handles.edit_ElementResults,'Visible','on','Enable','on');
        set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2],'Enable','off','Visible','on')
        draw = getappdata(0,'draw');
        draw.axialScaleFactor();
        draw.torsionScaleFactor();
        draw.shearScaleFactor_XY();
        draw.shearScaleFactor_XZ();
        draw.bendingMomentScaleFactor_XY();
        draw.bendingMomentScaleFactor_XZ();
        redraw(handles);
    end
end

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function pushbutton_PlayDynamicResults_Callback(hObject, ~, handles) %#ok<DEFNU>
% Check if button was pressed when it should not have been enabled
% Temporary, avoids errors.
if ~strcmp(get(handles.popupmenu_Results,'Enable'),'on')
    return
end

% Switch between play being on or off
switch get(hObject,'UserData')
    case true % Was on when user clicked
        set(hObject,'UserData',false,'ForegroundColor',[0,0.7,0.2])
        setappdata(0,'isDrawingDynamic',false)
    case false % Was off when user clicked
        set(hObject,'UserData',true,'ForegroundColor','Red')
        redraw(handles);
end

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function popupmenu_Results_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes on button press in "Textual" pushbutton from the Options menu.
% Creates a text file with the analysis results.
function pushbutton_Textual_Callback(hObject, eventdata, handles)%#ok<DEFNU>
% Disable button while txt is being printed
set(hObject,'Enable','off')

set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);

% Check if analysis is static or dynamic
if get(handles.popupmenu_AnalysisType,'Value') == 1
    analysis = 'Static';
elseif get(handles.popupmenu_AnalysisType,'Value') == 2
    analysis = 'Dynamic';
else
    set(hObject,'Enable','on');
    return
end

% Set file name
filterspec = '*.txt';
DialogTitle = 'LESM - Write analysis report';
gui_name = get(gcf,'Name');
if strcmp(gui_name,'LESM - Linear Elements Structure Model')
    DefaultName = 'untitled';
elseif strcmp(gui_name(end-3:end),'.lsm')
    DefaultName = gui_name(41:end-4);
else
    DefaultName = gui_name(41:end);
end
[filename,pathname] = uiputfile(filterspec,DialogTitle,DefaultName);
if ~filename
    % Enable button for future use
    set(hObject,'Enable','On')
    return
end
fullname = strcat(pathname,filename);

% Get handle to print and model objects
print = getappdata(0,'print');
model = getappdata(0,'model');
print.model = model;

switch analysis
    case 'Static'
        % Get current load case
        lc = get(handles.popupmenu_LoadCase,'Value');
        str = get(handles.popupmenu_LoadCase,'String');
        currentLc = str(lc,:);

        % Create file and write results
        print.txt = fopen(fullname,'wt');
        print.results(lc,currentLc,true);

    case 'Dynamic'
        % Create file and write results
        print.txt = fopen(fullname,'wt');
        print.results_modalDynamic();
end

% Close file
fclose(print.txt);

% Return object to root
setappdata(0,'print',print)

% Enable button for future use
set(hObject,'Enable','on')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% Executes on button press in "Reactions" checkbox.
% Displays reactions values next to nodal supports.
function checkbox_Reactions_Callback(hObject, ~, handles) %#ok<DEFNU>
draw = getappdata(0,'draw');
anm = get(handles.popupmenu_Anm,'value');

unselectEntities(handles);

if get(hObject,'Value') == 1
      draw.reactions();
      if anm == 1 || anm == 2
          ax = handles.axes_Canvas;
          ax.Clipping = 'on';
      else
          ax = handles.axes_Canvas;
          ax.Clipping = 'off';
      end
else
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

% Reinitialize object for mouse events and save it in root
mouse = getappdata(0,'mouse');
mouse.originalXLim = get(handles.axes_Canvas,'XLim');
mouse.originalYLim = get(handles.axes_Canvas,'YLim');
if (anm == 3) || (anm == 4) || (anm == 5)
    mouse.originalZLim = get(handles.axes_Canvas,'ZLim');
end
setappdata(0,'mouse',mouse);

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% Executes during object creation, after setting all properties.
function edit_Scale_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes when scale is moved or by user textual input.
function edit_Scale_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% Get calculated scale factor
dsf = getappdata(0,'deform_sf');

% Get scale input and slider position
editScale = str2double(get(handles.edit_Scale,'String'));
sliderScale = get(handles.slider_Scale,'Value');

if isnan(editScale) % if the input isnan, the slider sets the scale
     scale = dsf * sliderScale;
     set(handles.edit_Scale,'Enable','on','String',num2str(scale,3))
elseif (editScale/dsf) < get(handles.slider_Scale,'Min') % if the input is less than the Min value of slider, the Min value of the slider sets the scale
     scale = dsf * get(handles.slider_Scale,'Min');
     set(handles.slider_Scale,'Value',get(handles.slider_Scale,'Min'))
     set(handles.edit_Scale,'Enable','on','String',num2str(scale,3))
elseif (editScale/dsf) > get(handles.slider_Scale,'Max') % if the input is greater than the Max value of slider, the Max value of the slider sets the scale
     scale = dsf * get(handles.slider_Scale,'Max');
     set(handles.slider_Scale,'Value',get(handles.slider_Scale,'Max'))
     set(handles.edit_Scale,'Enable','on','String',num2str(scale,3))
else % if there is no problem with the input, it sets the scale
     scale = editScale;
     set(handles.slider_Scale,'Value',scale/dsf)
end     

unselectEntities(handles);

% Draw updated model
if get(handles.popupmenu_Results,'Value') ~=1
    redraw(handles);
end

%--------------------------------------------------------------------------
% Executes during slider creation, after setting all properties.
function slider_Scale_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%--------------------------------------------------------------------------
% Executes on slider movement.
% Adjusts deformed configuration to a new scale set by the slider.
function slider_Scale_Callback(hObject, eventdata, handles) %#ok<DEFNU>
if get(handles.popupmenu_AnalysisType,'Value') == 1
    dsf = getappdata(0,'deform_sf');
    set(handles.edit_Scale,'String',num2str(dsf*get(hObject,'Value')))
    if get(handles.popupmenu_Results,'Value') == 1
        return
    end    
    unselectEntities(handles);
    redraw(handles);
else
    if get(handles.popupmenu_Results,'Value') == 1
        return
    end
    setappdata(0,'isDrawingDynamic',false)
    set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])
    unselectEntities(handles);
    redraw(handles);
end

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% Executes on selection change in popupmenu_ElementResults.
% Selects which mode will have its results shown.
function popupmenu_ElementResults_Callback(hObject, eventdata, handles) %#ok<DEFNU>
anl = get(handles.popupmenu_AnalysisType,'Value');
res = get(handles.popupmenu_Results,'Value');
if anl == 1 || res ~= 2
    return;
end
setappdata(0,'isDrawingDynamic',false);
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2]);
unselectEntities(handles);
redraw(handles);

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_ElementResults creation.
function popupmenu_ElementResults_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes when enter key is pressed on edit_ElementResults.
function edit_ElementResults_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% Get selected nodes
n_str = get(hObject,'string');
if strcmp(n_str,'all') || strcmp(n_str,'ALL') || strcmp(n_str,'All')
    flag = 1;
else
    model = getappdata(0,'model');
    [flag,~] = readStr(n_str,model.nel);
end
if ~flag
    waitfor(msgbox('Invalid elements!','Error','error'));
    set(hObject,'string','All');
end

% Redraw diagram
draw = getappdata(0,'draw');
draw.axialScaleFactor();
draw.torsionScaleFactor();
draw.shearScaleFactor_XY();
draw.shearScaleFactor_XZ();
draw.bendingMomentScaleFactor_XY();
draw.bendingMomentScaleFactor_XZ();
if get(handles.popupmenu_AnalysisType,'Value') == 2
    setappdata(0,'isDrawingDynamic',false);
    set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2]);
end
unselectEntities(handles);
redraw(handles);

function [flag,output] = readStr(str,max)
output    = zeros(1,max);
count     = 0;       % counter for output index
numFlag   = false;   % flag for number being read
vctrFlag  = false;   % flag for vector being read
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

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_infoPanelEditable.
function uitable_infoPanelEditable_CellEditCallback(~, ~, ~)
mdata = guidata(findobj('Tag','GUI_Main'));
set(mdata.pushbutton_ApplyInfoPanel,'enable','on');

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_ApplyInfoPanel.
function pushbutton_ApplyInfoPanel_Callback(hObject, eventdata, handles) %#ok<DEFNU>
panelData = get(handles.uitable_infoPanelEditable,'Data');
mouse     = getappdata(0,'mouse');
anm       = get(handles.popupmenu_Anm,'value');
anl       = get(handles.popupmenu_AnalysisType,'Value');

% Check if there selected entity is a node or an element
if mouse.selectedNode ~= 0
    % Get node id, handle to node objects and current load case
    n     = mouse.selectedNode;
    nodes = getappdata(0,'nodes');
    lc    = get(handles.popupmenu_LoadCase,'value');
    
    % Get node ebc vector before changing it
    prevEbc = nodes(n).ebc;
    
    if anm <= 2 % TRUSS_2D or FRAME_2D
        i = 0; % initialize frame 2D index flag

        % Get new info provided by user
        newRestrX = char(panelData(1,2));
        switch newRestrX
            case 'y'
              nodes(n).ebc(1) = 1;
            case 'Y'
              nodes(n).ebc(1) = 1;
            case 'yes'
              nodes(n).ebc(1) = 1;
            case 'Yes'
              nodes(n).ebc(1) = 1;
            case 'YES'
              nodes(n).ebc(1) = 1;
            case 'n'
              nodes(n).ebc(1) = 0;
            case 'N'
              nodes(n).ebc(1) = 0;
            case 'no'
              nodes(n).ebc(1) = 0;
            case 'No'
              nodes(n).ebc(1) = 0;
            case 'NO'
              nodes(n).ebc(1) = 0;
        end
        newRestrY = char(panelData(2,2));
        switch newRestrY
            case 'y'
              nodes(n).ebc(2) = 1;
            case 'Y'
              nodes(n).ebc(2) = 1;
            case 'yes'
              nodes(n).ebc(2) = 1;
            case 'Yes'
              nodes(n).ebc(2) = 1;
            case 'YES'
              nodes(n).ebc(2) = 1;
            case 'n'
              nodes(n).ebc(2) = 0;
            case 'N'
              nodes(n).ebc(2) = 0;
            case 'no'
              nodes(n).ebc(2) = 0;
            case 'No'
              nodes(n).ebc(2) = 0;
            case 'NO'
              nodes(n).ebc(2) = 0;
        end
        if anm == 2 % FRAME 2D
            newRestrRotZ = char(panelData(3,2));
            i = 1;
            switch newRestrRotZ
                case 'y'
                  nodes(n).ebc(6) = 1;
                case 'Y'
                  nodes(n).ebc(6) = 1;
                case 'yes'
                  nodes(n).ebc(6) = 1;
                case 'Yes'
                  nodes(n).ebc(6) = 1;
                case 'YES'
                  nodes(n).ebc(6) = 1;
                case 'n'
                  nodes(n).ebc(6) = 0;
                case 'N'
                  nodes(n).ebc(6) = 0;
                case 'no'
                  nodes(n).ebc(6) = 0;
                case 'No'
                  nodes(n).ebc(6) = 0;
                case 'NO'
                  nodes(n).ebc(6) = 0;
            end
        end
        
        if anl == 1
            newFx = cell2mat(panelData(3+i,2));
            if isnan(newFx)
                newFx = 0;
            end
            newFy = cell2mat(panelData(4+i,2));
            if isnan(newFy)
                newFy = 0;
            end
            if anm == 2 % FRAME 2D
                newMz = cell2mat(panelData(5+i,2));
                if isnan(newMz)
                    newMz = 0;
                end
            else
                newMz = 0;
            end
            
            % Set new loads to current load case
            nodes(n).load.static(1) = newFx;
            nodes(n).load.static(2) = newFy;
            nodes(n).load.static(6) = newMz;
            if isempty(nodes(n).nodalLoadCase)
                nodes(n).nodalLoadCase = zeros(12,lc);
            end
            nodes(n).nodalLoadCase(1,lc) = newFx;
            nodes(n).nodalLoadCase(2,lc) = newFy;
            nodes(n).nodalLoadCase(6,lc) = newMz;
            
        elseif anl == 2
            newM = cell2mat(panelData(3+i,2));
            if isnan(newM)
                newM = 0;
            end
            
            % Set new mass
            nodes(n).displMass = newM*0.001;
        end
        
    elseif anm == 3 % GRILLAGE
        
        % Get new info provided by user
        newRestrZ = char(panelData(1,2));
        switch newRestrZ
            case 'y'
              nodes(n).ebc(3) = 1;
            case 'Y'
              nodes(n).ebc(3) = 1;
            case 'yes'
              nodes(n).ebc(3) = 1;
            case 'Yes'
              nodes(n).ebc(3) = 1;
            case 'YES'
              nodes(n).ebc(3) = 1;
            case 'n'
              nodes(n).ebc(3) = 0;
            case 'N'
              nodes(n).ebc(3) = 0;
            case 'no'
              nodes(n).ebc(3) = 0;
            case 'No'
              nodes(n).ebc(3) = 0;
            case 'NO'
              nodes(n).ebc(3) = 0;
        end
        newRestrRotXY = char(panelData(2,2));
        switch newRestrRotXY
            case 'y'
                nodes(n).ebc(4) = 1;
                nodes(n).ebc(5) = 1;
            case 'Y'
                nodes(n).ebc(4) = 1;
                nodes(n).ebc(5) = 1;
            case 'yes'
                nodes(n).ebc(4) = 1;
                nodes(n).ebc(5) = 1;
            case 'Yes'
                nodes(n).ebc(4) = 1;
                nodes(n).ebc(5) = 1;
            case 'YES'
                nodes(n).ebc(4) = 1;
                nodes(n).ebc(5) = 1;
            case 'n'
                nodes(n).ebc(4) = 0;
                nodes(n).ebc(5) = 0;
            case 'N'
                nodes(n).ebc(4) = 0;
                nodes(n).ebc(5) = 0;
            case 'no'
                nodes(n).ebc(4) = 0;
                nodes(n).ebc(5) = 0;
            case 'No'
                nodes(n).ebc(4) = 0;
                nodes(n).ebc(5) = 0;
            case 'NO'
                nodes(n).ebc(4) = 0;
                nodes(n).ebc(5) = 0;
        end
        
        if anl == 1
            newFz = cell2mat(panelData(3,2));
            if isnan(newFz)
                newFz = 0;
            end
            newMx = cell2mat(panelData(4,2));
            if isnan(newMx)
                newMx = 0;
            end
            newMy = cell2mat(panelData(5,2));
            if isnan(newMy)
                newMy = 0;
            end
            
            % Set new loads to current load case
            if isempty(nodes(n).load.static)
                nodes(n).load.static = zeros(1,6);
            end
            nodes(n).load.static(3) = newFz;
            nodes(n).load.static(4) = newMx;
            nodes(n).load.static(5) = newMy;
            if isempty(nodes(n).nodalLoadCase)
                nodes(n).nodalLoadCase = zeros(12,lc);
            end
            nodes(n).nodalLoadCase(3,lc) = newFz;
            nodes(n).nodalLoadCase(4,lc) = newMx;
            nodes(n).nodalLoadCase(5,lc) = newMy;
            
        elseif anl == 2
            newM = cell2mat(panelData(3,2));
            if isnan(newM)
                newM = 0;
            end
            
            % Set new mass
            nodes(n).displMass = newM*0.001;
        end
        
    elseif anm >= 4 % TRUSS_3D or FRAME_3D
        i = 0; % initialize frame 3D index flag

        % Get new info provided by user
        newRestrX = char(panelData(1,2));
        switch newRestrX
            case 'y'
              nodes(n).ebc(1) = 1;
            case 'Y'
              nodes(n).ebc(1) = 1;
            case 'yes'
              nodes(n).ebc(1) = 1;
            case 'Yes'
              nodes(n).ebc(1) = 1;
            case 'YES'
              nodes(n).ebc(1) = 1;
            case 'n'
              nodes(n).ebc(1) = 0;
            case 'N'
              nodes(n).ebc(1) = 0;
            case 'no'
              nodes(n).ebc(1) = 0;
            case 'No'
              nodes(n).ebc(1) = 0;
            case 'NO'
              nodes(n).ebc(1) = 0;
        end
        newRestrY = char(panelData(2,2));
        switch newRestrY
            case 'y'
              nodes(n).ebc(2) = 1;
            case 'Y'
              nodes(n).ebc(2) = 1;
            case 'yes'
              nodes(n).ebc(2) = 1;
            case 'Yes'
              nodes(n).ebc(2) = 1;
            case 'YES'
              nodes(n).ebc(2) = 1;
            case 'n'
              nodes(n).ebc(2) = 0;
            case 'N'
              nodes(n).ebc(2) = 0;
            case 'no'
              nodes(n).ebc(2) = 0;
            case 'No'
              nodes(n).ebc(2) = 0;
            case 'NO'
              nodes(n).ebc(2) = 0;
        end
        newRestrZ = char(panelData(3,2));
        switch newRestrZ
            case 'y'
              nodes(n).ebc(3) = 1;
            case 'Y'
              nodes(n).ebc(3) = 1;
            case 'yes'
              nodes(n).ebc(3) = 1;
            case 'Yes'
              nodes(n).ebc(3) = 1;
            case 'YES'
              nodes(n).ebc(3) = 1;
            case 'n'
              nodes(n).ebc(3) = 0;
            case 'N'
              nodes(n).ebc(3) = 0;
            case 'no'
              nodes(n).ebc(3) = 0;
            case 'No'
              nodes(n).ebc(3) = 0;
            case 'NO'
              nodes(n).ebc(3) = 0;
        end
        if anm == 5 % FRAME 3D
            newRestrRot = char(panelData(4,2));
            i = 1;
            switch newRestrRot
                case 'y'
                  nodes(n).ebc(4) = 1;
                  nodes(n).ebc(5) = 1;
                  nodes(n).ebc(6) = 1;
                case 'Y'
                  nodes(n).ebc(4) = 1;
                  nodes(n).ebc(5) = 1;
                  nodes(n).ebc(6) = 1;
                case 'yes'
                  nodes(n).ebc(4) = 1;
                  nodes(n).ebc(5) = 1;
                  nodes(n).ebc(6) = 1;
                case 'Yes'
                  nodes(n).ebc(4) = 1;
                  nodes(n).ebc(5) = 1;
                  nodes(n).ebc(6) = 1;
                case 'YES'
                  nodes(n).ebc(4) = 1;
                  nodes(n).ebc(5) = 1;
                  nodes(n).ebc(6) = 1;
                case 'n'
                  nodes(n).ebc(4) = 0;
                  nodes(n).ebc(5) = 0;
                  nodes(n).ebc(6) = 0;
                case 'N'
                  nodes(n).ebc(4) = 0;
                  nodes(n).ebc(5) = 0;
                  nodes(n).ebc(6) = 0;
                case 'no'
                  nodes(n).ebc(4) = 0;
                  nodes(n).ebc(5) = 0;
                  nodes(n).ebc(6) = 0;
                case 'No'
                  nodes(n).ebc(4) = 0;
                  nodes(n).ebc(5) = 0;
                  nodes(n).ebc(6) = 0;
                case 'NO'
                  nodes(n).ebc(4) = 0;
                  nodes(n).ebc(5) = 0;
                  nodes(n).ebc(6) = 0;
            end
        end
        
        if anl == 1
            newFx = cell2mat(panelData(4+i,2));
            if isnan(newFx)
                newFx = 0;
            end
            newFy = cell2mat(panelData(5+i,2));
            if isnan(newFy)
                newFy = 0;
            end
            newFz = cell2mat(panelData(6+i,2));
            if isnan(newFz)
                newFz = 0;
            end
            if anm == 5 % FRAME 3D
                newMx = cell2mat(panelData(7+i,2));
                if isnan(newMx)
                    newMx = 0;
                end
                newMy = cell2mat(panelData(8+i,2));
                if isnan(newMy)
                    newMy = 0;
                end
                newMz = cell2mat(panelData(9+i,2));
                if isnan(newMz)
                    newMz = 0;
                end
            else
                newMx = 0;
                newMy = 0;
                newMz = 0;
            end
            
            % Set new loads to current load case
            nodes(n).load.static(1) = newFx;
            nodes(n).load.static(2) = newFy;
            nodes(n).load.static(3) = newFz;
            nodes(n).load.static(4) = newMx;
            nodes(n).load.static(5) = newMy;
            nodes(n).load.static(6) = newMz;
            if isempty(nodes(n).nodalLoadCase)
                nodes(n).nodalLoadCase = zeros(12,lc);
            end
            nodes(n).nodalLoadCase(1,lc) = newFx;
            nodes(n).nodalLoadCase(2,lc) = newFy;
            nodes(n).nodalLoadCase(3,lc) = newFz;
            nodes(n).nodalLoadCase(4,lc) = newMx;
            nodes(n).nodalLoadCase(5,lc) = newMy;
            nodes(n).nodalLoadCase(6,lc) = newMz;
            
        elseif anl == 2
            newM = cell2mat(panelData(4+i,2));
            if isnan(newM)
                newM = 0;
            end
            
            % Set new mass
            nodes(n).displMass = newM*0.001;
        end
    end
    
    % Check if inclined supports were changed
    if ~all((prevEbc(1:3) == nodes(n).ebc(1:3))) && nodes(n).isInclinedSupp
        nodes(n).removeInclinedSupp;
    end
    
    % Check if any support was removed
    for j = 1:size(prevEbc,2)
        if nodes(n).ebc(j) ~= prevEbc(j)
            if isempty(nodes(n).prescDispl)
                nodes(n).prescDispl = zeros(1,size(prevEbc,2));
            else
                nodes(n).prescDispl(j) = 0;
            end
            nodes(n).nodalLoadCase(7:12,lc) = nodes(n).prescDispl';
            if nodes(n).ebc(j) == 0 && prevEbc(j) == 1
               nodes(n).nodalLoadCase(j+6,1:end) = zeros(1,size(nodes(n).nodalLoadCase,2));
            end
        end
    end
    
    % Get number of supports added to node and update information panel
    countNewFree = 0;
    countNewFixed = 0;
    countNewSpring = 0;
    countRemovedFree = 0;
    countRemovedFixed = 0;
    countRemovedSpring = 0;
    for i = 1:size(prevEbc,2)
        if nodes(n).ebc(i) == 0 && prevEbc(i) ~= 0
            countNewFree = countNewFree + 1;
        elseif nodes(n).ebc(i) == 1 && prevEbc(i) ~= 1
            countNewFixed = countNewFixed + 1;
        elseif nodes(n).ebc(i) == 2 && prevEbc(i) ~= 2
            countNewSpring = countNewSpring + 1;
        end
        if nodes(n).ebc(i) ~= 0 && prevEbc(i) == 0
            countRemovedFree = countRemovedFree + 1;
        elseif nodes(n).ebc(i) ~= 1 && prevEbc(i) == 1
            countRemovedFixed = countRemovedFixed + 1;
        elseif nodes(n).ebc(i) ~= 2 && prevEbc(i) == 2
            countRemovedSpring = countRemovedSpring + 1;
            % Delete spring stiffness coefficient from removed spring
            nodes(n).springStiff(i) = 0;
        end
    end
    
    % Update info panel data
    infoPanelData = mouse.originalData;
    nfreedof = cell2mat(infoPanelData(6,2)) + countNewFree - countRemovedFree;
    nfixeddof = cell2mat(infoPanelData(7,2)) + countNewFixed - countRemovedFixed;
    nspringdof = cell2mat(infoPanelData(8,2)) + countNewSpring - countRemovedSpring;
    infoPanelData(6,:) = {'Free DOFs',nfreedof};
    infoPanelData(7,:) = {'Fixed DOFs',nfixeddof};
    infoPanelData(8,:) = {'Springs',nspringdof};
    mouse.originalData = infoPanelData;
    
    % Update nodes in root
    setappdata(0,'nodes',nodes)
    
    % Disable results, enable process data button and modelling options
    set(handles.pushbutton_ProcessData,'Enable','on');
    if get(handles.popupmenu_Results,'Value') ~= 1
        allLoadsNeedToBeRedrawn = true;
    else
        allLoadsNeedToBeRedrawn = false;
    end
    set(handles.popupmenu_Results,'Enable','off','Value',1);
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.edit_Scale,'enable','off','visible','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.pushbutton_DynamicResults,'enable','off');
    set(handles.pushbutton_PlayDynamicResults,'Enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    set(handles.pushbutton_Materials,'Enable','on');
    set(handles.pushbutton_Sections,'Enable','on');
    set(handles.pushbutton_Nodes,'Enable','on');
    set(handles.pushbutton_Elements,'Enable','on');
    set(handles.pushbutton_NodalLoads,'Enable','on');
    set(handles.pushbutton_Supports,'Enable','on');
    if anl == 1
        set(handles.pushbutton_ElementLoads,'Enable','on');
    elseif anl == 2
        set(handles.pushbutton_ElementLoads,'Enable','off');
    end
    
    % Redraw model
    redraw(handles,'Nodes',false)
    if allLoadsNeedToBeRedrawn == false
        redraw(handles,'Nodal Loads')
    else
        redraw(handles,'Loads')
    end
    
elseif mouse.selectedElem ~= 0
    
    % Get element id, handle to elem objects, model and current load case
    e = mouse.selectedElem;
    elems = getappdata(0,'elems');
    model = getappdata(0,'model');
    lc = get(handles.popupmenu_LoadCase,'value');
    
    % Initialize flag
    elemsNeedToBeRedrawn = false;
    
    % Get new info provided by user
    if anm == 2 || anm == 3 || anm == 5 % FRAME 2D or GRILLAGE or FRAME 3D
        newType = char(panelData(1,2));
        switch newType
            case 'y'
              elems(e).type = 1;
            case 'Y'
              elems(e).type = 1;
            case 'yes'
              elems(e).type = 1;
            case 'Yes'
              elems(e).type = 1;
            case 'YES'
              elems(e).type = 1;
            case 'n'
              elems(e).type = 0;
            case 'N'
              elems(e).type = 0;
            case 'no'
              elems(e).type = 0;
            case 'No'
              elems(e).type = 0;
            case 'NO'
              elems(e).type = 0;
        end
        newHi = char(panelData(2,2));
        hingei = elems(e).hingei;
        switch newHi
            case 'y'
              elems(e).hingei = 0;
              elems(e).kri = [];
            case 'Y'
              elems(e).hingei = 0;
              elems(e).kri = [];
            case 'yes'
              elems(e).hingei = 0;
              elems(e).kri = [];
            case 'Yes'
              elems(e).hingei = 0;
              elems(e).kri = [];
            case 'YES'
              elems(e).hingei = 0;
              elems(e).kri = [];
            case 'n'
              elems(e).hingei = 1;
              elems(e).kri = [];
            case 'N'
              elems(e).hingei = 1;
              elems(e).kri = [];
            case 'no'
              elems(e).hingei = 1;
              elems(e).kri = [];
            case 'No'
              elems(e).hingei = 1;
              elems(e).kri = [];
            case 'NO'
              elems(e).hingei = 1;
              elems(e).kri = [];
        end
        if hingei ~= elems(e).hingei
            elemsNeedToBeRedrawn = true;
        end
        newHf = char(panelData(3,2));
        hingef = elems(e).hingef;
        switch newHf
            case 'y'
              elems(e).hingef = 0;
              elems(e).krf = [];
            case 'Y'
              elems(e).hingef = 0;
              elems(e).krf = [];
            case 'yes'
              elems(e).hingef = 0;
              elems(e).krf = [];
            case 'Yes'
              elems(e).hingef = 0;
              elems(e).krf = [];
            case 'YES'
              elems(e).hingef = 0;
              elems(e).krf = [];
            case 'n'
              elems(e).hingef = 1;
              elems(e).krf = [];
            case 'N'
              elems(e).hingef = 1;
              elems(e).krf = [];
            case 'no'
              elems(e).hingef = 1;
              elems(e).krf = [];
            case 'No'
              elems(e).hingef = 1;
              elems(e).krf = [];
            case 'NO'
              elems(e).hingef = 1;
              elems(e).krf = [];
        end
        if hingef ~= elems(e).hingef
            elemsNeedToBeRedrawn = true;
        end
        i = 3;
    else
        i = 0;
    end
    
    newMat = cell2mat(panelData(1+i,2));
    if ~isnan(newMat) && newMat <= getappdata(0,'nmat') && newMat > 0 && mod(newMat,1) == 0
        elems(e).material = model.materials(newMat);
    end
    newSec = cell2mat(panelData(2+i,2));
    if ~isnan(newSec) && newSec <= getappdata(0,'nsec') && newSec > 0 && mod(newSec,1) == 0
        elems(e).section = model.sections(newSec);
    end
    
    % Loads
    if anl == 1
        if anm <= 2 % TRUSS 2D or FRAME 2D
            newQx1 = cell2mat(panelData(3+i,2));
            if isnan(newQx1)
                newQx1 = 0;
            end
            newQy1 = cell2mat(panelData(4+i,2));
            if isnan(newQy1)
                newQy1 = 0;
            end
            newQx2 = cell2mat(panelData(5+i,2));
            if isnan(newQx2)
                newQx2 = 0;
            end
            newQy2 = cell2mat(panelData(6+i,2));
            if isnan(newQy2)
                newQy2 = 0;
            end
            newDtx = cell2mat(panelData(7+i,2));
            if isnan(newDtx)
                newDtx = 0;
            end
            newDty = cell2mat(panelData(8+i,2));
            if isnan(newDty)
                newDty = 0;
            end

            % Set new loads to current load case
            elems(e).load.elemLoadCase(:,lc) = zeros(14,1);
            elems(e).load.elemLoadCase(5,lc) = 0;
            elems(e).load.elemLoadCase(6,lc) = newQx1;
            elems(e).load.elemLoadCase(7,lc) = newQy1;
            elems(e).load.elemLoadCase(9,lc) = newQx2;
            elems(e).load.elemLoadCase(10,lc) = newQy2;
            elems(e).load.elemLoadCase(12,lc) = newDtx;
            elems(e).load.elemLoadCase(13,lc) = newDty;
            elems(e).load.uniformGbl = []; % clear previous distributed loads
            elems(e).load.uniformLcl = [];
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            elems(e).load.setLinearLoad((elems(e).load.elemLoadCase(6:11,lc))',elems(e).load.elemLoadCase(5,lc));
            elems(e).load.tempVar_X = newDtx;
            elems(e).load.tempVar_Y = newDty;

        elseif anm == 3 % GRILLAGE
            newQz1 = cell2mat(panelData(3+i,2));
            if isnan(newQz1)
                newQz1 = 0;
            end
            newQz2 = cell2mat(panelData(4+i,2));
            if isnan(newQz2)
                newQz2 = 0;
            end
            newDtz = cell2mat(panelData(5+i,2));
            if isnan(newDtz)
                newDtz = 0;
            end

            % Set new loads to current load case
            elems(e).load.elemLoadCase(:,lc) = zeros(14,1);
            elems(e).load.elemLoadCase(5,lc) = 0;
            elems(e).load.elemLoadCase(8,lc) = newQz1;
            elems(e).load.elemLoadCase(11,lc) = newQz2;
            elems(e).load.elemLoadCase(14,lc) = newDtz;
            elems(e).load.uniformGbl = []; % clear previous distributed loads
            elems(e).load.uniformLcl = [];
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            elems(e).load.setLinearLoad((elems(e).load.elemLoadCase(6:11,lc))',elems(e).load.elemLoadCase(5,lc));
            elems(e).load.tempVar_Z = newDtz;

        elseif anm == 4 || anm == 5 % TRUSS 3D or FRAME 3D

            newQx1 = cell2mat(panelData(3+i,2));
            if isnan(newQx1)
                newQx1 = 0;
            end
            newQy1 = cell2mat(panelData(4+i,2));
            if isnan(newQy1)
                newQy1 = 0;
            end
            newQz1 = cell2mat(panelData(5+i,2));
            if isnan(newQz1)
                newQz1 = 0;
            end
            newQx2 = cell2mat(panelData(6+i,2));
            if isnan(newQx2)
                newQx2 = 0;
            end
            newQy2 = cell2mat(panelData(7+i,2));
            if isnan(newQy2)
                newQy2 = 0;
            end
            newQz2 = cell2mat(panelData(8+i,2));
            if isnan(newQz2)
                newQz2 = 0;
            end
            newDtx = cell2mat(panelData(9+i,2));
            if isnan(newDtx)
                newDtx = 0;
            end
            newDty = cell2mat(panelData(10+i,2));
            if isnan(newDty)
                newDty = 0;
            end
            newDtz = cell2mat(panelData(11+i,2));
            if isnan(newDtz)
                newDtz = 0;
            end

            % Set new loads to current load case
            elems(e).load.elemLoadCase(:,lc) = zeros(14,1);
            elems(e).load.elemLoadCase(5,lc) = 0;
            elems(e).load.elemLoadCase(6,lc) = newQx1;
            elems(e).load.elemLoadCase(7,lc) = newQy1;
            elems(e).load.elemLoadCase(8,lc) = newQz1;
            elems(e).load.elemLoadCase(9,lc) = newQx2;
            elems(e).load.elemLoadCase(10,lc) = newQy2;
            elems(e).load.elemLoadCase(11,lc) = newQz2;
            elems(e).load.elemLoadCase(12,lc) = newDtx;
            elems(e).load.elemLoadCase(13,lc) = newDty;
            elems(e).load.elemLoadCase(14,lc) = newDtz;
            elems(e).load.uniformGbl = []; % clear previous distributed loads
            elems(e).load.uniformLcl = [];
            elems(e).load.linearGbl = [];
            elems(e).load.linearLcl = [];
            elems(e).load.setLinearLoad((elems(e).load.elemLoadCase(6:11,lc))',elems(e).load.elemLoadCase(5,lc));
            elems(e).load.tempVar_X = newDtx;
            elems(e).load.tempVar_Y = newDty;
            elems(e).load.tempVar_Z = newDtz;
        end

        % Avoids numeric problems with new loads
        for i = 1:size(elems(e).load.linearLcl,1)
            if abs(elems(e).load.linearLcl(i)) <= 10^-10
                elems(e).load.linearLcl(i) = 0;
            end
        end
    end
    
    % Update elems in root
    setappdata(0,'elems',elems)
    
    % Disable results, enable process data button and modelling options
    set(handles.pushbutton_ProcessData,'Enable','on');
    if get(handles.popupmenu_Results,'Value') ~= 1
        allLoadsNeedToBeRedrawn = true;
    else
        allLoadsNeedToBeRedrawn = false;
    end
    set(handles.popupmenu_Results,'Enable','off','Value',1);
    set(handles.text_Element,'string','Elements');
    set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(handles.pushbutton_Textual,'Enable','off');
    set(handles.edit_Scale,'enable','off','visible','off');
    set(handles.checkbox_Reactions,'Enable','off','Value',0);
    set(handles.pushbutton_DynamicResults,'enable','off');
    set(handles.pushbutton_PlayDynamicResults,'Enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    set(handles.pushbutton_Materials,'Enable','on');
    set(handles.pushbutton_Sections,'Enable','on');
    set(handles.pushbutton_Nodes,'Enable','on');
    set(handles.pushbutton_Elements,'Enable','on');
    set(handles.pushbutton_NodalLoads,'Enable','on');
    set(handles.pushbutton_Supports,'Enable','on');
    if anl == 1
        set(handles.pushbutton_ElementLoads,'Enable','on');
    elseif anl == 2
        set(handles.pushbutton_ElementLoads,'Enable','off');
    end
    
    % Redraw model
    if elemsNeedToBeRedrawn == true
        redraw(handles,'Elements')
    end
    if allLoadsNeedToBeRedrawn == false
        redraw(handles,'Element Loads')
    else
        redraw(handles,'Loads')
    end
end

% Disable button and reset info panels
set(hObject,'enable','off')
set(handles.uitable_infoPanelEditable,'enable','off','Data',{})
set(handles.uitable_infoPanel,'Data',mouse.originalData)

% Disable delete entities pushbutton (toolbar)
set(handles.pushbutton_DeleteEntities,'enable','off');

% Update mouse
mouse.selectedNode = 0;
mouse.selectedElem = 0;
delete(findobj('tag', 'snapNode'));
delete(findobj('tag', 'snapNode2'));
delete(findobj('tag', 'snapElem'));
delete(findobj('tag', 'snapElem2'));
delete(findobj('tag', 'selectedNode'));
delete(findobj('tag', 'selectedElem'));
setappdata(0,'mouse',mouse)

%--------------------------------------------------------------------------
function visualizationMenu_Callback(~,~,~) %#ok<DEFNU>

%--------------------------------------------------------------------------
function rulerButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        handles.axes_Canvas.XAxis.Visible = 'off';
        handles.axes_Canvas.YAxis.Visible = 'off';
        handles.axes_Canvas.ZAxis.Visible = 'off';
    case 'off'
        set(hObject,'Checked','on')
        handles.axes_Canvas.XAxis.Visible = 'on';
        handles.axes_Canvas.YAxis.Visible = 'on';
        handles.axes_Canvas.ZAxis.Visible = 'on';
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function gridButton_Callback(hObject, ~, ~) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        grid off
    case 'off'    
        set(hObject,'Checked','on')
        grid on
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function unitsButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
    case 'off'    
        set(hObject,'Checked','on')
end

% unselectEntities(handles)

redraw(handles,'Units')

if getappdata(0,'vis') == 1
    % Save changes to handles structure
    guidata(hObject,handles)
end

%--------------------------------------------------------------------------
function nodeIDButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% unselectEntities(handles);

check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        value = 0;
    case 'off'    
        set(hObject,'Checked','on')
        value = 1;
end

if getappdata(0,'vis') == 1
    draw = getappdata(0,'draw');
    
    if value == 1
        draw.nodeID(); 
    else
        delete(findobj('tag','textNodeID'))
    end
end

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function elemIDButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% unselectEntities(handles);

check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        value = 0;
    case 'off'    
        set(hObject,'Checked','on')
        value = 1;
end

if getappdata(0,'vis') == 1
    draw = getappdata(0,'draw');
    
    if value == 1
        draw.elementID(); 
    else
        delete(findobj('tag','textElemID'))
    end
end

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function orientationButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% unselectEntities(handles);

check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        value = 0;
    case 'off'    
        set(hObject,'Checked','on')
        value = 1;
end

if getappdata(0,'vis') == 1
    draw = getappdata(0,'draw');
    
    if value == 1
        draw.elementOrientation(); 
    else
        delete(findobj('tag','drawElemOrient'))
    end
end

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function viewNodalLoadsButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawNodalLoads'))
        delete(findobj('tag','textNodalLoads'))
        delete(findobj('tag','textNodalMoments'))
    case 'off'
        set(hObject,'Checked','on')
        if get(handles.popupmenu_Results,'value') == 1  % MODEL
            draw = getappdata(0,'draw');
            anl = get(handles.popupmenu_AnalysisType,'Value');
            if anl == 1
                draw.nodalLoads();
            elseif anl == 2
                draw.dynamicNodalLoads();
            end
        end
end

%--------------------------------------------------------------------------
function viewNodalMassButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawNodalMass'))
        delete(findobj('tag','textNodalMass'))
    case 'off'    
        set(hObject,'Checked','on')
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 2 && get(handles.popupmenu_Results,'value') == 1 % MODEL
            draw = getappdata(0,'draw');
            draw.nodalMass();
        end
end

%--------------------------------------------------------------------------
function viewPrescDisplButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawPrescDispl'))
        delete(findobj('tag','textPrescDispl'))
        delete(findobj('tag','textPrescRot'))
    case 'off'    
        set(hObject,'Checked','on')
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 1 && get(handles.popupmenu_Results,'value') == 1  % MODEL
            draw = getappdata(0,'draw');
            draw.nodalPrescDispl();
        end
end

%--------------------------------------------------------------------------
function viewInitialConditionsButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','textInitialConditions'))
    case 'off'    
        set(hObject,'Checked','on')
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 2 && get(handles.popupmenu_Results,'value') == 1  % MODEL
            draw = getappdata(0,'draw');
            draw.nodalInitialConditions();
        end
end

%--------------------------------------------------------------------------
function viewDistribLoadsButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawElemLoads'))
        delete(findobj('tag','textElemLoads'))
    case 'off'    
        set(hObject,'Checked','on')
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 1 && get(handles.popupmenu_Results,'value') == 1  % MODEL
            draw = getappdata(0,'draw');
            draw.elemLoadsScaleFactor();
            draw.elemLoads();
        end
end

%--------------------------------------------------------------------------
function viewThermalLoadsButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawThermalLoads'))
        delete(findobj('tag','textThermalLoads'))
    case 'off'    
        set(hObject,'Checked','on')
        anl = get(handles.popupmenu_AnalysisType,'Value');
        if anl == 1 && get(handles.popupmenu_Results,'value') == 1  % MODEL
            draw = getappdata(0,'draw');
            draw.thermalLoads();
        end
end

%--------------------------------------------------------------------------
function viewSupportsButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawSupports'))
        delete(findobj('tag','textSprings'))
        delete(findobj('tag','textRotSprings'))
    case 'off'    
        set(hObject,'Checked','on')
        redraw(handles,'Nodes',false)
end

%--------------------------------------------------------------------------
function viewSemiRigidButton_Callback(hObject, ~, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'on'
        set(hObject,'Checked','off')
        delete(findobj('tag','drawSemiRigid'))
        delete(findobj('tag','textSemiRigid'))
    case 'off'    
        set(hObject,'Checked','on')
        delete(findobj('tag','drawSemiRigidTemp'))
end
redraw(handles,'Elements')

%--------------------------------------------------------------------------
function planeButton_Callback(~,~,~) %#ok<DEFNU>

%--------------------------------------------------------------------------
function plane3D_Callback(hObject, eventdata, handles) %#ok<DEFNU>
view2D = get(handles.togglebutton_2DView,'state');
switch view2D
    case 'off'
        check = get(hObject,'Checked');
        view(3);
        switch check
            case 'off'    
                set(hObject,'Checked','on')
                set(handles.planeXY,'Checked','off')
                set(handles.planeXZ,'Checked','off')
                set(handles.planeYZ,'Checked','off')
        end
    case 'on'
        set(handles.togglebutton_2DView,'state','off')
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function planeXY_Callback(hObject, ~, handles) %#ok<DEFNU>
if get(handles.popupmenu_Anm,'Value') - 1 == 2 && strcmp(get(handles.togglebutton_2DView,'state'),'off')
    view(0,90+1e-10); % at exactly 90 degrees, zooming behavior gets strange for grillage
else
    view(0,90)
end
switch get(hObject,'Checked')
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.plane3D,'Checked','off')
        set(handles.planeXZ,'Checked','off')
        set(handles.planeYZ,'Checked','off')
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function planeXZ_Callback(hObject, eventdata, handles) %#ok<DEFNU>
view(0,0);
switch get(hObject,'Checked')
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.planeXY,'Checked','off')
        set(handles.plane3D,'Checked','off')
        set(handles.planeYZ,'Checked','off')
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function planeYZ_Callback(hObject, eventdata, handles) %#ok<DEFNU>
view(90,0);
switch get(hObject,'Checked')
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.planeXY,'Checked','off')
        set(handles.planeXZ,'Checked','off')
        set(handles.plane3D,'Checked','off')
end

% unselectEntities(handles);

%--------------------------------------------------------------------------
function formatButton_Callback(~,~,~) %#ok<DEFNU>

%--------------------------------------------------------------------------
function precision0_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',0)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision1_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision0,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',1)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision2_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',2)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision3_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',3)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision4_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',4)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision5_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',5)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision6_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',6)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision7_Callback(hObject, eventdata, handles)  %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',7)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision8_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision0,'Checked','off')
        set(handles.precision9,'Checked','off')
        setappdata(0,'decPrec',8)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function precision9_Callback(hObject, eventdata, handles) %#ok<DEFNU>
check = get(hObject,'Checked');
switch check
    case 'off'    
        set(hObject,'Checked','on')
        set(handles.precision1,'Checked','off')
        set(handles.precision2,'Checked','off')
        set(handles.precision3,'Checked','off')
        set(handles.precision4,'Checked','off')
        set(handles.precision5,'Checked','off')
        set(handles.precision6,'Checked','off')
        set(handles.precision7,'Checked','off')
        set(handles.precision8,'Checked','off')
        set(handles.precision0,'Checked','off')
        setappdata(0,'decPrec',9)
end

% unselectEntities(handles);
redraw(handles,'Decimal Precision')

% Save changes to handles structure
guidata(hObject,handles)

%--------------------------------------------------------------------------
function dynResOptButton_Callback(~, ~, handles) %#ok<DEFNU>
setappdata(0,'move',1);
setappdata(0,'isDrawingDynamic',false)
set(handles.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])

set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

unselectEntities(handles);

GUI_DynResOpt();

%--------------------------------------------------------------------------
function staticAnalysisCallback(handles)
setappdata(0,'move',1);
set(handles.popupmenu_AnalysisType,'Value',1);
set(handles.pushbutton_DynamicAnalysis,'Enable','off');
set(handles.pushbutton_TimeFcn,'Enable','off');
set(handles.dynResOptButton,'Enable','off');
set(handles.pushbutton_ElementLoads,'Enable','on');
set(handles.text_Element,'string','Elements');
set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(handles.pushbutton_DynamicResults,'Enable','off','Visible','off');
set(handles.pushbutton_PlayDynamicResults,'Enable','off','Visible','off');

% Reset solver as static
model = getappdata(0,'model');
model.whichSolver = 0;
model.drv = [];
model.drv = Drv_LES(true,model);
setappdata(0,'model',model);

%--------------------------------------------------------------------------
function dynamicAnalysisCallback(handles)
setappdata(0,'move',1);
set(handles.popupmenu_AnalysisType,'value',2);
set(handles.pushbutton_DynamicAnalysis,'Enable','on');
set(handles.pushbutton_TimeFcn,'Enable','on');
set(handles.dynResOptButton,'Enable','on');
set(handles.pushbutton_ElementLoads,'Enable','off');
set(handles.text_Element,'string','Elements');
set(handles.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(handles.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(handles.pushbutton_DynamicResults,'Enable','off','Visible','on');
set(handles.pushbutton_PlayDynamicResults,'Enable','off','Visible','on');

% Reset solver as dynamic (numerical)
model = getappdata(0,'model');
model.whichSolver = 1;
model.drv = [];
model.drv = Drv_LED(1,true,model);
setappdata(0,'model',model);

%--------------------------------------------------------------------------
% Auxiliary function.
% Unselect entities (nodes or elements) and resets info panel.
function unselectEntities(handles)

delete(findobj('tag','selectedNode'));
delete(findobj('tag','selectedElem'));
delete(findobj('tag','snapNode'));
delete(findobj('tag','snapNode2'));
delete(findobj('tag','snapElem'));
delete(findobj('tag','snapElem2'));
delete(findobj('tag','snapGrid'));
delete(findobj('tag','snapIntSect'))
delete(findobj('tag','dynamicLine'));

mouse = getappdata(0,'mouse');
if ~isempty(mouse.originalData)
    if mouse.selectedNode ~= 0
        n = mouse.selectedNode;
        mouse.selectedNode = 0;
        mouse.whichNodeSnap = n;
    elseif mouse.selectedElem ~= 0
        e = mouse.selectedElem;
        mouse.selectedElem = 0;
        mouse.whichElemSnap = e;
    end
    mouse.moveAction();

    set(handles.uitable_infoPanelEditable,'enable','off','Data',{})
    set(handles.uitable_infoPanel,'Data',mouse.originalData)
end
setappdata(0,'mouse',mouse)
set(handles.pushbutton_ApplyInfoPanel,'enable','off')
set(handles.pushbutton_DeleteEntities,'enable','off')

%--------------------------------------------------------------------------
% Auxiliary function.
% Set all togglebuttons' state to 'off', except the one that was just
% toggled on.
function toggleButtonsOff(hObject,handles)
% Get tag of the toggle button that was toggled on.
if ~isempty(hObject)
    buttonOnName = get(hObject,'Tag');
else
    set(handles.togglebutton_Cursor,'state','on')
    buttonOnName = 'togglebutton_Cursor';
end

if ~strcmp(buttonOnName,'togglebutton_Cursor')
    set(handles.togglebutton_Cursor,'state','off')
end

if ~strcmp(buttonOnName,'toolbar_pan')
    set(handles.toolbar_pan,'state','off')
end

if ~strcmp(buttonOnName,'toolbar_zoomIn')
    set(handles.toolbar_zoomIn,'state','off')
end

if ~strcmp(buttonOnName,'toolbar_zoomOut')
    set(handles.toolbar_zoomOut,'state','off')
end

if ~strcmp(buttonOnName,'toolbar_rotate3d')
    set(handles.toolbar_rotate3d,'state','off')
end

if ~strcmp(buttonOnName,'togglebutton_Node')
    set(handles.togglebutton_Node,'state','off')
end

if ~strcmp(buttonOnName,'togglebutton_Element')
    set(handles.togglebutton_Element,'state','off')
end

if ~strcmp(buttonOnName,'togglebutton_SnapToGrid')
    if ~strcmp(buttonOnName,'togglebutton_Node') && ~strcmp(buttonOnName,'togglebutton_Element')
        set(handles.togglebutton_SnapToGrid,'state','off')
    end
end

if ~strcmp(buttonOnName,'togglebutton_CrossElements')
    if ~strcmp(buttonOnName,'togglebutton_Element')
        set(handles.togglebutton_CrossElements,'enable','off')
    end
end

if ~strcmp(buttonOnName,'togglebutton_Polyline')
    if ~strcmp(buttonOnName,'togglebutton_Element')
        set(handles.togglebutton_Polyline,'enable','off')
    end
end

%--------------------------------------------------------------------------
% Auxiliary function.
% Check if all toggle buttons are off.
function flag = areAllToggleOff(handles)
flag = 0;

t1 = char(get(handles.toolbar_pan,'state'));
t2 = char(get(handles.toolbar_zoomIn,'state'));
t3 = char(get(handles.toolbar_zoomOut,'state'));
t4 = char(get(handles.toolbar_rotate3d,'state'));
t5 = char(get(handles.togglebutton_Node,'state'));
t6 = char(get(handles.togglebutton_Element,'state'));
t7 = char(get(handles.togglebutton_SnapToGrid,'state'));

toggleButtonState = char(t1,t2,t3,t4,t5,t6,t7);

flagVector = zeros(1,size(toggleButtonState,1));

for i = 1:size(toggleButtonState,1)
    if strcmp(toggleButtonState(i,:),'off')
        flagVector(i) = 1;
    end
end

if all(flagVector == 1)
    flag = 1;
end

% Equivalent to:
%
% if all(toggleButtonState == 'off')
%     flag = 1;
% end

%--------------------------------------------------------------------------
% Executes when GUI_Main is resized.
function GUI_Main_SizeChangedFcn(hObject, ~, handles) 
anm = get(handles.popupmenu_Anm,'Value') - 1;

% Make sure that modifications will be made on canvas
axes(handles.axes_Canvas);

% Adjust 2D Canvas
if anm == 0 || anm == 1 || (anm == 2 && strcmp(get(handles.togglebutton_2DView,'state'),'on'))
    % Adjust canvas position
    dfltUnits = get(handles.axes_Canvas, 'Units');
    set(handles.axes_Canvas, 'Units', 'normalized');
    set(handles.axes_Canvas, 'Position', [0.207,0.049,0.776,0.926]);
    set(handles.axes_Canvas, 'Units', dfltUnits);
    
     % Get draw object
     draw = getappdata(0,'draw');
     axis equal;
     if ~isempty(draw)
         draw.setLimits();
     else
         xlim([-5,5]);
     end
    
    % Update mouse properties
    mouse = getappdata(0,'mouse');
    mouse.originalXLim = get(handles.axes_Canvas,'XLim');
    mouse.originalYLim = get(handles.axes_Canvas,'YLim');
    setappdata(0,'mouse',mouse);
    
else
    % Adjust canvas position
    dfltUnits = get(handles.axes_Canvas, 'Units');
    set(handles.axes_Canvas, 'Units', 'normalized');
    set(handles.axes_Canvas, 'Position', [0.23,0.1266,0.72,0.8576]);
    set(handles.axes_Canvas, 'Units', dfltUnits);
    
    % Get draw object
     draw = getappdata(0,'draw');
     axis equal;
     draw.setLimits();
    
    % Update mouse properties
    mouse = getappdata(0,'mouse');
    mouse.originalXLim = get(handles.axes_Canvas,'XLim');
    mouse.originalYLim = get(handles.axes_Canvas,'YLim');
    
    set(handles.axes_Canvas, 'Units', 'pixels');
    mouse.originalAxesPos = get(handles.axes_Canvas,'Position');
    set(handles.axes_Canvas, 'Units', dfltUnits);

    setappdata(0,'mouse',mouse);
end

% Adjust uitable_infoPanel position
dfltUnitsInfoTable = get(handles.uitable_infoPanel,'Units');
set(handles.uitable_infoPanel,'Units','pixels')
infoTablePosition = get(handles.uitable_infoPanel,'Position');
set(handles.uitable_infoPanel, 'ColumnWidth', {95*infoTablePosition(3)/166.3  52.5*infoTablePosition(3)/166.3})
set(handles.uitable_infoPanel,'Units',dfltUnitsInfoTable)

% Adjust uitable_infoPanelEditable position
dfltUnitsInfoTable = get(handles.uitable_infoPanelEditable,'Units');
set(handles.uitable_infoPanelEditable,'Units','pixels')
infoTablePosition = get(handles.uitable_infoPanelEditable,'Position');
set(handles.uitable_infoPanelEditable, 'ColumnWidth', {95*infoTablePosition(3)/166.3  52.5*infoTablePosition(3)/166.3})
set(handles.uitable_infoPanelEditable,'Units',dfltUnitsInfoTable)

%--------------------------------------------------------------------------
% Executes when specific key is pressed.
function GUI_Main_KeyPressFcn(~, eventdata, handles) %#ok<DEFNU>
key = get(gcf, 'CurrentKey');

resetMouse = true;

switch key
    case 'return'
        if strcmp(get(handles.pushbutton_ProcessData,'Enable'),'on')         
            pushbutton_ProcessData_Callback(handles.pushbutton_ProcessData, eventdata, handles);
        end
        resetMouse = false;
        
    case 'escape'
        togglebutton_Cursor_OnCallback(handles.togglebutton_Cursor, eventdata, handles);
        resetMouse = false;
        
    case 'delete'
        if strcmp(get(handles.togglebutton_Cursor,'state'),'on')
            pushbutton_DeleteEntities_ClickedCallback(handles.pushbutton_DeleteEntities, eventdata, handles)
        else
            resetMouse = false;
        end
        
    case 'n'
        if strcmp(get(handles.togglebutton_Node,'state'),'off')
            togglebutton_Node_OnCallback(handles.togglebutton_Node, eventdata, handles);
            set(handles.togglebutton_Node,'state','on');
        end
        resetMouse = false;
        
    case 'e'
        if strcmp(get(handles.togglebutton_Element,'state'),'off')
            togglebutton_Element_OnCallback(handles.togglebutton_Element, eventdata, handles);
            set(handles.togglebutton_Element,'state','on');
        end
        resetMouse = false;
        
    case 'f'
        pushbutton_FitWorld_ClickedCallback(handles.pushbutton_FitWorld,eventdata,handles);
        resetMouse = false;
        
    case 'r'
        pushbutton_RefreshModel_ClickedCallback(handles.pushbutton_RefreshModel,eventdata,handles);
        resetMouse = false;
        
    case 'i'
        nodeIDButton_Callback(handles.nodeIDButton,eventdata,handles);
        resetMouse = false;
        
    case 'j'
        elemIDButton_Callback(handles.elemIDButton,eventdata,handles);
        resetMouse = false;
        
    case 'o'
        orientationButton_Callback(handles.orientationButton,eventdata,handles);
        resetMouse = false;
end

if resetMouse == true
    % Delete dynamic plots
    delete(findobj('tag','dynamicLine'));
    delete(findobj('tag','snapGrid'));
    delete(findobj('tag','snapElem'));
    delete(findobj('tag','snapElem2'));
    delete(findobj('tag','snapNode'));
    delete(findobj('tag','snapNode2'));
    delete(findobj('tag','snapIntSect'))
    delete(findobj('tag','selectedNode'));
    delete(findobj('tag','selectedElem'));
    
    % Reset modeling variables
    mouse = getappdata(0,'mouse');
    mouse.elemNode = 0;
    mouse.elemNodeID = [];
    mouse.elemCoords = [];
    mouse.elemResults = [];
    mouse.selectedNode = 0;
    setappdata(0,'mouse',mouse)
end

%--------------------------------------------------------------------------
% Executes when user attempts to close GUI_Main.
function GUI_Main_CloseRequestFcn(hObject, eventdata, handles) %#ok<DEFNU>
set(handles.togglebutton_Node,'state','off')
set(handles.togglebutton_Element,'state','off')

model = getappdata(0,'model');

if model.nnp ~= 0
    choice = questdlg('All unsaved data will be lost. Do you want to continue?',...
        'Exit','Save and exit','Exit','Cancel','Cancel');

    switch choice
        case 'Save and exit'
            filename = saveButton_Callback([],eventdata,handles);
            if filename ~= 0
                appdata = get(0,'ApplicationData');
                fns = fieldnames(appdata);
                for i = 1:numel(fns)
                  rmappdata(0,fns{i});
                end
                delete(hObject);
            end
        case 'Exit'
            appdata = get(0,'ApplicationData');
            fns = fieldnames(appdata);
            for i = 1:numel(fns)
              rmappdata(0,fns{i});
            end
            delete(hObject);
        case 'Cancel'
            return
    end
else
    appdata = get(0,'ApplicationData');
    fns = fieldnames(appdata);
    for i = 1:numel(fns)
      rmappdata(0,fns{i});
    end
    delete(hObject);
end

function pushbutton_MetalCalculate_ClickedCallback(~, ~, handles) %#ok<DEFNU> %ronald

designReport=DesignReport();
model = getappdata(0,'model');
fid = fopen('teste_site.html','w');
fprintf(fid,'%s',designReport.criarHtml());
fprintf(fid,'%s',designReport.createSideMenu(size(model.elems,2)));


 wrapper=["<div id='page-content-wrapper'>"
          "<div class='container'>"];
 
      
 fprintf(fid,'%s',wrapper);     
%cabeï¿½alho

header=[     "<div class='text-center' style='width: 100%;'>"
             "<img src='lesm_logo.png' />"
             "<p style='margin-top: 40px;' class='font-weight-bold'>Pontifical Catholic University of Rio de Janeiro</p>"
             "</div>"
             "<div class='row' style='margin-top: 50px;'>"
             "<div class='col-md-12 text-center'>"
             "<h1 class='mt-4'>Memorial de Calculo</h1>"
             "<hr class='my-4'>"
             "</div>"
];
fprintf(fid,'%s',header);
%% 

answer = questdlg('which one standard would you like to use?', ...
	'Standars', ...
	'EUROCODE3','NBR8800','NBR8800');
% Handle response
switch answer
    case 'EUROCODE3'
      
    case 'NBR8800'
        norma=NBR8800(fid);
        model.designSolver(norma,fid);
        
    case 'No thank you'
       
end









endWrapper=["</div></div>"];
  
       
fprintf(fid,'%s',endWrapper);   
fprintf(fid,'%s',designReport.fecharHtml());

fclose(fid);
