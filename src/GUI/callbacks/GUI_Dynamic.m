%% Dynamic Parameters Dialog Callback Functions
% This file contains the callback functions associated with the "Dynamic
% Parameters" dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Dynamic(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Dynamic_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Dynamic_OutputFcn, ...
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
%--------------------------------------------------------------------------
% --- Executes just before GUI_Dynamic is made visible.
function GUI_Dynamic_OpeningFcn(hObject, ~, handles, varargin)
include_constants;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Get model object from root and handle to GUI_Main
model = getappdata(0,'model');
mdata = guidata(findobj('Tag','GUI_Main'));

% Set options according to response type
response = model.drv.whichResponse;
if response == TRANS_MODAL_ANALYSIS
    turnOnTransient(handles);
else
    turnOffTransient(handles);
end

% Set number of modes
infoPanelData = get(mdata.uitable_infoPanel,'Data');
nf = infoPanelData{6,2} + infoPanelData{8,2};
if model.n_modes > 0
    current_nmodes = model.n_modes;
elseif nf > 0
    current_nmodes = nf;
else
    current_nmodes = 0;
end
set(handles.edit_Nmodes,'string',num2str(current_nmodes));

% Set mass matrix
if model.mass_type == CONSISTENT_MASS
    set(handles.radiobutton_ConsistentMass,'value',1);
    set(handles.radiobutton_LumpedMass,'value',0);
    set(handles.radiobutton_MixedMass,'value',0);
    set(handles.edit_MiMassCoeff,'enable','off','string','1');
    set(handles.text_MiMassEqn,'enable','off','visible','on');
elseif model.mass_type == LUMPED_MASS
    set(handles.radiobutton_ConsistentMass,'value',0);
    set(handles.radiobutton_LumpedMass,'value',1);
    set(handles.radiobutton_MixedMass,'value',0);
    set(handles.edit_MiMassCoeff,'enable','off','string','0');
    set(handles.text_MiMassEqn,'enable','off','visible','on');
else
    set(handles.radiobutton_ConsistentMass,'value',0);
    set(handles.radiobutton_LumpedMass,'value',0);
    set(handles.radiobutton_MixedMass,'value',1);
    set(handles.edit_MiMassCoeff,'enable','on','string',num2str(model.mass_mi));
    set(handles.text_MiMassEqn,'enable','on','visible','on');
end

% Choose default command line output for GUI_Dynamic
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = GUI_Dynamic_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in radiobutton_Transient.
function radiobutton_Transient_Callback(~, ~, handles) %#ok<*DEFNU>
turnOnTransient(handles);

%--------------------------------------------------------------------------
% --- Executes on button press in radiobutton_Modal.
function radiobutton_Modal_Callback(~, ~, handles) %#ok<*DEFNU>
turnOffTransient(handles);

%--------------------------------------------------------------------------
function edit_Nmodes_Callback(hObject, ~, ~)
% Get handle to model object from root and handle to GUI_Main
model = getappdata(0,'model');
mdata = guidata(findobj('Tag','GUI_Main'));

% Get user input
nmodes = str2double(get(hObject,'string'));

% Get number of free dofs from panel in GUI_Main
infoPanelData = get(mdata.uitable_infoPanel,'Data');
nf = infoPanelData{6,2} + infoPanelData{8,2};

% Get current number of modes
if model.n_modes > 0
    current_nmodes = model.n_modes;
elseif nf > 0
    current_nmodes = nf;
else
    current_nmodes = 0;
end

% Check if input is valid
if isnan(nmodes) || nmodes <= 0 || nmodes > nf || rem(nmodes,1) ~= 0
    set(hObject,'string',num2str(current_nmodes));
    return
end

%--------------------------------------------------------------------------
% --- Executes on button press in radiobutton_ConsistentMass.
function radiobutton_ConsistentMass_Callback(~, ~, handles) %#ok<*DEFNU>
set(handles.edit_MiMassCoeff,'enable','off','string','1');
set(handles.text_MiMassEqn,'enable','off','visible','on');

%--------------------------------------------------------------------------
% --- Executes on button press in radiobutton_LumpedMass.
function radiobutton_LumpedMass_Callback(~, ~, handles) %#ok<*DEFNU>
set(handles.edit_MiMassCoeff,'enable','off','string','0');
set(handles.text_MiMassEqn,'enable','off','visible','on');

%--------------------------------------------------------------------------
% --- Executes on button press in radiobutton_MixedMass.
function radiobutton_MixedMass_Callback(~, ~, handles) %#ok<*DEFNU>
set(handles.edit_MiMassCoeff,'enable','on');
set(handles.text_MiMassEqn,'enable','on','visible','on');

%--------------------------------------------------------------------------
% --- Executes on edited text in edit_MiMassCoeff.
function edit_MiMassCoeff_Callback(hObject, ~, ~) %#ok<*DEFNU>
input = str2double(get(hObject,'string'));
if isnan(input)
    set(hObject,'string','1');
elseif input > 1
    set(hObject,'string','1');
elseif input < 0
    set(hObject,'string','0');
end

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Solver.
function popupmenu_Solver_Callback(~, ~, ~)

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Damping.
function popupmenu_Damping_Callback(hObject, ~, handles)
include_constants;

% Get model object from root
model = getappdata(0,'model');

% Get user selection on popupmenu
set(handles.edit_Damp1,'enable','on','visible','on');

switch get(hObject,'Value')
    case 1 % Critical damping ratio of 1st vibration mode
        set(handles.text_Damp1,'enable','on','visible','on');
        set(handles.text_Damp2,'enable','off','visible','off');
        set(handles.text_DampAlpha,'enable','off','visible','off');
        set(handles.text_DampBeta,'enable','off','visible','off');
        
        if isempty(model.xi)
            xi = 0;
        else
            xi = model.xi(1);
        end
        set(handles.edit_Damp1,'string',num2str(xi));
        set(handles.edit_Damp2,'enable','off','visible','off');
        set(handles.text_Coefficients,'enable','off','visible','off');
        
    case 2 % Critical damping ratio of 1st and 2nd vibration mode
        set(handles.text_Damp1,'enable','on','visible','on');
        set(handles.text_Damp2,'enable','on','visible','on');
        set(handles.text_DampAlpha,'enable','off','visible','off');
        set(handles.text_DampBeta,'enable','off','visible','off');
        
        if isempty(model.xi)
            xi = [0 0];
        else
            xi = model.xi;
        end
        set(handles.edit_Damp1,'string',num2str(xi(1)));
        set(handles.edit_Damp2,'enable','on','visible','on','string',num2str(xi(2)));
        set(handles.text_Coefficients,'enable','off','visible','off');
        
    case 3 % Rayleigh coefficients
        set(handles.text_Damp1,'enable','off','visible','off');
        set(handles.text_Damp2,'enable','off','visible','off');
        set(handles.text_DampAlpha,'enable','on','visible','on');
        set(handles.text_DampBeta,'enable','on','visible','on');
        
        set(handles.edit_Damp1,'string',num2str(model.massDampCoeff));
        set(handles.edit_Damp2,'enable','on','visible','on','string',num2str(model.stiffDampCoeff));
        set(handles.text_Coefficients,'enable','on','visible','on');
end

% Return model object to root
setappdata(0,'model',model);

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Nodes.
function popupmenu_Nodes_Callback(hObject, ~, handles)
% Get selected node and dof
n          = get(hObject,'value');
DOF_string = get(handles.popupmenu_DOF,'string');
DOF_value  = get(handles.popupmenu_DOF,'value');
DOF        = DOF_string{DOF_value};
switch DOF
    case 'Dx'
        dof = 1;
    case 'Dy'
        dof = 2;
    case 'Dz'
        dof = 3;
    case 'Rx'
        dof = 4;
    case 'Ry'
        dof = 5;
    case 'Rz'
        dof = 6;
end

% Get model object from root
model = getappdata(0,'model');

% Get selected dof initial conditions
c0 = model.nodes(n).initCond(dof,:);

% Set initial conditions of dof as strings of editable textboxes
set(handles.edit_d0,'string',num2str(c0(1),4))
set(handles.edit_v0,'string',num2str(c0(2),4))

% Set string of static texts (units)
if dof <= 3
    set(handles.text_d0,'string','[m]')
    set(handles.text_v0,'string','[m/s]')    
else % dof > 3
    set(handles.text_d0,'string','[rad]')
    set(handles.text_v0,'string','[rad/s]')    
end

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_DOF.
function popupmenu_DOF_Callback(hObject, ~, handles)
% Get selected node and dof
n          = get(handles.popupmenu_Nodes,'value');
DOF_string = get(hObject,'string');
DOF_value  = get(hObject,'value');
DOF        = DOF_string{DOF_value};
switch DOF
    case 'Dx'
        dof = 1;
    case 'Dy'
        dof = 2;
    case 'Dz'
        dof = 3;
    case 'Rx'
        dof = 4;
    case 'Ry'
        dof = 5;
    case 'Rz'
        dof = 6;
end

% Get model object from root
model = getappdata(0,'model');

% Get selected dof initial conditions
c0 = model.nodes(n).initCond(dof,:);

% Set initial conditions of dof as strings of editable textboxes
set(handles.edit_d0,'string',num2str(c0(1),4))
set(handles.edit_v0,'string',num2str(c0(2),4))

% Set string of static texts (units)
if dof <= 3
    set(handles.text_d0,'string','[m]')
    set(handles.text_v0,'string','[m/s]')    
else % dof > 3
    set(handles.text_d0,'string','[rad]')
    set(handles.text_v0,'string','[rad/s]')    
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(~,~, handles)
% Get model object from root and handle to GUI_Main
include_constants;
model = getappdata(0,'model');
mdata = guidata(findobj('Tag','GUI_Main'));

% Response type
if get(handles.radiobutton_Transient,'Value')
    response = TRANS_MODAL_ANALYSIS;
else
    response = MODAL_ANALYSIS;
end
loadsNeedToBeRedrawn    = model.drv.whichResponse ~= response;
model.drv.whichResponse = response;

% Number of modes
edit_Nmodes_Callback(handles.edit_Nmodes,[],[]);
model.n_modes = str2double(get(handles.edit_Nmodes,'string'));

% Mass matrix
if get(handles.radiobutton_ConsistentMass,'Value')
    model.mass_type = CONSISTENT_MASS;
elseif get(handles.radiobutton_LumpedMass,'Value')
    model.mass_type = LUMPED_MASS;
else
    model.mass_type = MIXED_MASS;
    model.mass_mi = str2double(get(handles.edit_MiMassCoeff,'string'));
end

% Solver, damping options, and initial conditions
if model.drv.whichResponse == TRANS_MODAL_ANALYSIS
    % Algorithm
    switch get(handles.popupmenu_Solver,'Value')
        case 1
            model.whichSolver = DYNAMIC_NEWMARK_LINEAR;
            model.drv = [];
            model.drv = Drv_LED(DYNAMIC_NEWMARK_LINEAR,true,model);
        case 2
            model.whichSolver = DYNAMIC_RK4_LINEAR;
            model.drv = [];
            model.drv = Drv_LED(DYNAMIC_RK4_LINEAR,true,model);
        case 3
            model.whichSolver = DYNAMIC_AM3_LINEAR;
            model.drv = [];
            model.drv = Drv_LED(DYNAMIC_AM3_LINEAR,true,model);
        case 4
            model.whichSolver = DYNAMIC_WILSON_LINEAR;
            model.drv = [];
            model.drv = Drv_LED(DYNAMIC_WILSON_LINEAR,true,model);
        case 5
            model.whichSolver = DYNAMIC_MODALSUP_LINEAR;
            model.drv = [];
            model.drv = Drv_LED(DYNAMIC_MODALSUP_LINEAR,true,model);
    end
    
    % Get time parameters
    dt = str2double(get(handles.edit_dt,'string'));
    ns = str2double(get(handles.edit_nSteps,'string'));
    t  = ns*dt;
    if isnan(dt) || isnan(ns)
        msgbox('Invalid solver parameters!','Error','error');
        return
    elseif dt <= 0
        msgbox('Invalid solver parameters! Time interval must be a positive value.','Error','error');
        return
    elseif rem(ns,1) ~= 0 || ns <= 0
        msgbox('Invalid solver parameters! Number of steps must be a positive integer.','Error','error');
        return
    elseif ns > 10^6
        msgbox('Invalid solver parameters! Maximum number of steps is 10^6.','Error','error');
        return
    end
    
    % Update time functions
    if (model.t ~= t) && (model.n_steps ~= ns)
        for fcn = model.timeFcns
            if ~isempty(fcn{1})
                fcn{1}.update_time_nsteps(t,ns);
            end
        end
    elseif (model.t ~= t)
        for fcn = model.timeFcns
            if ~isempty(fcn{1})
                fcn{1}.update_time(t);
            end
        end
    elseif (model.n_steps ~= ns)
        for fcn = model.timeFcns
            if ~isempty(fcn{1})
                fcn{1}.update_nsteps(ns);
            end
        end
    end
    
    % Set time parameters
    model.t       = t;
    model.n_steps = ns;
    
    % Damping parameters
    switch get(handles.popupmenu_Damping,'value')
        case 1
            model.damping = XI_1ST_MODE;
            xi = str2double(get(handles.edit_Damp1,'string'));
            if isnan(xi)
                xi = model.xi(1);
                set(handles.edit_Damp1,'string',num2str(xi));
            elseif xi < 0
                xi = model.xi(1);
                set(handles.edit_Damp1,'string',num2str(xi));
            else
                model.xi(1) = xi;
                model.xi(2) = xi;
            end
        case 2
            model.damping = XI_1ST_2ND_MODE;
            xi_1 = str2double(get(handles.edit_Damp1,'string'));
            xi_2 = str2double(get(handles.edit_Damp2,'string'));
            if isnan(xi_1)
                xi_1 = model.xi(1);
                set(handles.edit_Damp1,'string',num2str(xi_1));
            elseif xi_1 < 0
                xi_1 = model.xi(1);
                set(handles.edit_Damp1,'string',num2str(xi_1));
            else
                model.xi(1) = xi_1;
            end
            if isnan(xi_2)
                xi_2 = model.xi(2);
                set(handles.edit_Damp2,'string',num2str(xi_2));
            elseif xi_2 < 0
                xi_2 = model.xi(2);
                set(handles.edit_Damp2,'string',num2str(xi_2));
            else
                model.xi(2) = xi_2;
            end
        case 3
            model.damping = RAYLEIGH_COEFFS;
            alpha = str2double(get(handles.edit_Damp1,'string'));
            beta  = str2double(get(handles.edit_Damp2,'string'));
            if isnan(alpha)
                alpha = model.massDampCoeff;
                set(handles.edit_Damp1,'string',num2str(alpha));
            elseif alpha < 0
                alpha = model.massDampCoeff;
                set(handles.edit_Damp1,'string',num2str(alpha));
            else
                model.massDampCoeff = alpha;
            end
            if isnan(beta)
                beta = model.stiffDampCoeff;
                set(handles.edit_Damp2,'string',num2str(beta));
            elseif beta < 0
                beta = model.stiffDampCoeff;
                set(handles.edit_Damp2,'string',num2str(beta));
            else
                model.stiffDampCoeff = beta;
            end
    end
    
    % Initial conditions
    if model.nnp ~= 0
        nodes = getappdata(0,'nodes');
        
        % Get ICs
        d0 = str2double(get(handles.edit_d0,'string'));
        v0 = str2double(get(handles.edit_v0,'string'));
        if isnan(d0) || isnan(v0)
            msgbox('Invalid initial condition values!','Error','error');
            return
        end
        
        % Get node and dof IDs
        n          = get(handles.popupmenu_Nodes,'value');
        DOF_string = get(handles.popupmenu_DOF,'string');
        DOF_value  = get(handles.popupmenu_DOF,'value');
        DOF        = DOF_string{DOF_value};
        switch DOF
            case 'Dx'
                dof = 1;
            case 'Dy'
                dof = 2;
            case 'Dz'
                dof = 3;
            case 'Rx'
                dof = 4;
            case 'Ry'
                dof = 5;
            case 'Rz'
                dof = 6;
        end
        
        % Set redraw flag
        if ~isequal(full(nodes(n).initCond(dof,:)),[d0,v0])
            loadsNeedToBeRedrawn = true;
        end
        
        % Set IC to dof and return node objects to root
        nodes(n).initCond(dof,:) = [d0,v0];
        model.nodes = nodes;
        setappdata(0,'nodes',nodes);
        
        % Update uitable data
        tableData = get(handles.uitable_InitCond,'Data');
        if isempty(tableData) && ~all(nodes(n).initCond(dof,:)==0)
            tableData(1,:) = {n,DOF,...
                              nodes(n).initCond(dof,1),...
                              nodes(n).initCond(dof,2)};
        else
            nodesData = zeros(size(tableData,1),1);
            dofData = char(zeros(size(tableData,1),2));
            for aux = 1:size(tableData,1)
                nodesData(aux,1) = tableData{aux,1};
                dofData(aux,:)   = tableData{aux,2};
            end
            
            equalFlag = false;
            stop = false;
            i = 1;
            max_i = length(nodesData);
            while i <= max_i && ~stop
                if n < nodesData(i)
                    stop = true;
                elseif n == nodesData(i)
                    while n == nodesData(i) && ~stop
                        if all(DOF == dofData(i,:))
                            equalFlag = true;
                            stop = true;
                        else
                            logicVctr = (DOF < dofData(i,:));
                            if logicVctr(1) || (~logicVctr(1) && logicVctr(2))
                                stop = true;
                            elseif i >= max_i
                                i = i + 1;
                                break
                            else
                                i = i + 1;
                            end
                        end
                    end
                else
                    i = i + 1;
                end
            end
            
            if i <= max_i && ~equalFlag && (nodes(n).initCond(dof,1) ~= 0 || nodes(n).initCond(dof,2) ~= 0)
                tableData(i+1:end+1,:) = tableData(i:end,:);
            end
            if equalFlag && (nodes(n).initCond(dof,1) == 0 && nodes(n).initCond(dof,2) == 0)
                tableData(i,:) = [];
            elseif nodes(n).initCond(dof,1) ~= 0 || nodes(n).initCond(dof,2) ~= 0
                tableData(i,:) = {n,DOF,...
                                  nodes(n).initCond(dof,1),...
                                  nodes(n).initCond(dof,2)};
            end
        end
        set(handles.uitable_InitCond,'enable','on','Data',tableData);
    end
end

% Return model object to root
setappdata(0,'model',model);

% Enable "Process Data" button in main GUI
if model.nnp > 0
    set(mdata.pushbutton_ProcessData,'Enable','on');
end

% Disable result buttons
if get(mdata.popupmenu_Results,'value') ~= 1
    loadsNeedToBeRedrawn = true;
end
set(mdata.popupmenu_Results,'Enable','off','value',1);
set(mdata.pushbutton_Textual,'Enable','off');
set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
set(mdata.text_Element,'string','Elements');
set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
set(mdata.edit_Scale,'enable','off','visible','off');
set(mdata.pushbutton_DynamicResults,'enable','off');
set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);

% Make GUI a normal window
gui = findobj('Tag','GUI_Dynamic');
set(gui,'WindowStyle','normal');

% Update model drawing, if necessary
if loadsNeedToBeRedrawn
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function edit_Nmodes_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Damp1_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Damp1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Damp2_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Damp2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_d0_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_d0_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_v0_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_v0_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_nSteps_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Nodes_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_DOF_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Solver_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Damping_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Callback edit_dt
function edit_dt_Callback(~,~,~)

% --- Callback edit_nSteps
function edit_nSteps_Callback(~,~,~)

%--------------------------------------------------------------------------
function turnOnTransient(handles)
% Get model object from root
include_constants;
model = getappdata(0,'model');

% Response type
set(handles.radiobutton_Transient,'Value',1);
set(handles.radiobutton_Modal,'Value',0);

% Solver
solver_opt = {'Newmark';'Runge-Kutta';'Adams-Moulton';'Wilson-Theta';'Modal superposition'};
set(handles.popupmenu_Solver,'enable','on','string',solver_opt);
switch model.whichSolver
    case DYNAMIC_NEWMARK_LINEAR
        set(handles.popupmenu_Solver,'value',1);
    case DYNAMIC_RK4_LINEAR
        set(handles.popupmenu_Solver,'value',2);
    case DYNAMIC_AM3_LINEAR
        set(handles.popupmenu_Solver,'value',3);
    case DYNAMIC_WILSON_LINEAR
        set(handles.popupmenu_Solver,'value',4);
    case DYNAMIC_MODALSUP_LINEAR
        set(handles.popupmenu_Solver,'value',5);
end
if model.n_steps ~= 0
    set(handles.edit_dt,'enable','on','string',num2str(model.t/model.n_steps,4));
else
    set(handles.edit_dt,'enable','on','string',num2str(0));
end
set(handles.edit_nSteps,'enable','on','string',num2str(model.n_steps,4));

% Damping
damp_opt = {'Critical ratio (1st mode)';'Critical ratio (1st & 2nd modes)';'Rayleigh coefficients'};
set(handles.popupmenu_Damping,'enable','on','string',damp_opt);
switch model.damping
    case XI_1ST_MODE
        set(handles.popupmenu_Damping,'value',1);
    case XI_1ST_2ND_MODE
        set(handles.popupmenu_Damping,'value',2);
    case RAYLEIGH_COEFFS
        set(handles.popupmenu_Damping,'value',3);
end
popupmenu_Damping_Callback(handles.popupmenu_Damping,[],handles);

% Initial conditions
if model.nnp > 0
    % Set nodes popupmenu properties
    string_nodes = num2str(1:model.nnp,'%d\n');
    set(handles.popupmenu_Nodes,'string',string_nodes,'value',1,'Max',model.nnp,'enable','on')
    
    % Determine string and maximum value of the dof popupmenu
    switch model.anm.analysis_type
        case TRUSS2D_ANALYSIS
            string_dof = {'Dx';'Dy'};
            max_dof    = 2;
            id         = 1;
        case FRAME2D_ANALYSIS
            string_dof = {'Dx';'Dy';'Rz'};
            max_dof    = 3;
            id         = 1;
        case GRILLAGE_ANALYSIS
            string_dof = {'Rx';'Ry';'Dz'};
            max_dof    = 3;
            id         = 4;
        case TRUSS3D_ANALYSIS
            string_dof = {'Dx';'Dy';'Dz'};
            max_dof    = 3;
            id         = 1;
        case FRAME3D_ANALYSIS
            string_dof = {'Dx';'Dy';'Dz';'Rx';'Ry';'Rz'};
            max_dof    = 6;
            id         = 1;
    end
    
    % Set dof popupmenu properties
    set(handles.popupmenu_DOF,'string',string_dof,'value',1,'Max',max_dof,'enable','on')
    
    % Set string of static texts (units)
    if id == 1
        set(handles.text_d0,'string','[m]')
        set(handles.text_v0,'string','[m/s]')
    else % id == 4
        set(handles.text_d0,'string','[rad]')
        set(handles.text_v0,'string','[rad/s]')
    end
    
    % Enable initial conditions editable textboxes
    set(handles.edit_d0,'enable','on','string',num2str(model.nodes(1).initCond(id,1),4))
    set(handles.edit_v0,'enable','on','string',num2str(model.nodes(1).initCond(id,2),4))
    
    % Assemble table data to be displayed at uitable_InitCond
    DOF_string = {'Dx';'Dy';'Dz';'Rx';'Ry';'Rz'};
    tableData = cell(model.nnp*6,4);
    counter = 0;
    for n = 1:model.nnp
        for i = 1:6
            counter = counter + 1;
            if isempty(model.nodes(n).initCond)
                tableData(counter,:) = [];
                counter = counter - 1;
            elseif all(model.nodes(n).initCond(i,:) == zeros(1,2))
                tableData(counter,:) = [];
                counter = counter - 1;
            else
                DOF = DOF_string{i};
                tableData(counter,:) = {n,DOF,...
                    model.nodes(n).initCond(i,1),...
                    model.nodes(n).initCond(i,2)};
            end
        end
    end
    set(handles.uitable_InitCond,'enable','on','Data',tableData);
else
    set(handles.popupmenu_Nodes,'string',' ','value',1,'Max',1,'enable','off')
    set(handles.popupmenu_DOF,'string',' ','value',1,'Max',1,'enable','off')
    set(handles.text_d0,'string','[m]')
    set(handles.text_v0,'string','[m/s]')
    set(handles.edit_d0,'enable','off','string','')
    set(handles.edit_v0,'enable','off','string','')
    set(handles.uitable_InitCond,'enable','off','Data',{});
end

%--------------------------------------------------------------------------
function turnOffTransient(handles)
set(handles.radiobutton_Transient,'Value',0);
set(handles.radiobutton_Modal,'Value',1);
set(handles.popupmenu_Solver,'value',1,'enable','off','string',' ');
set(handles.edit_dt,'string','','enable','off');
set(handles.edit_nSteps,'string','','enable','off');
set(handles.popupmenu_Damping,'value',1,'enable','off','string',' ');
set(handles.text_Damp1,'enable','off','visible','off');
set(handles.edit_Damp1,'enable','off','visible','off');
set(handles.text_Damp2,'enable','off','visible','off');
set(handles.edit_Damp2,'enable','off','visible','off');
set(handles.text_DampAlpha,'enable','off','visible','off');
set(handles.text_DampBeta,'enable','off','visible','off');
set(handles.text_Coefficients,'enable','off','visible','off');
set(handles.popupmenu_Nodes,'string',' ','value',1,'Max',1,'enable','off')
set(handles.popupmenu_DOF,'string',' ','value',1,'Max',1,'enable','off')
set(handles.text_d0,'string','[m]')
set(handles.text_v0,'string','[m/s]')
set(handles.edit_d0,'enable','off','string','')
set(handles.edit_v0,'enable','off','string','')
set(handles.uitable_InitCond,'enable','off','Data',{});
