%% Dynamic Nodal Loads Dialog Callback Functions
% This file contains the callback functions associated with the "Dynamic
% Nodal Loads" dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
function varargout = GUI_NodalLoads_Dynamic(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_NodalLoads_Dynamic_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_NodalLoads_Dynamic_OutputFcn, ...
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


% --- Executes just before GUI_NodalLoads_Dynamic is made visible.
function GUI_NodalLoads_Dynamic_OpeningFcn(hObject, ~, handles, varargin)
include_constants;

% Move GUI to the center of the screen
movegui(hObject,'center');
            
% Get model object from root
model = getappdata(0,'model');

% Compute number of nodes
nnp = length(model.nodes);

% Set popupmenu_Node properties
set(handles.popupmenu_Node,'string',num2str(1:nnp,'%d\n'),'value',1,'Max',nnp);

% Set popupmenu_TimeFcns properties
fcnId = model.findFcnListByHandle(model.nodes(1).load.getFcn());
if isempty(model.strTimeFcns)
    set(handles.popupmenu_TimeFcns,'value',1,'max',1,'enable','off','string','None');
else
    str = model.strTimeFcns; str(2:end+1) = str; str{1} = 'None';
    set(handles.popupmenu_TimeFcns,'value',fcnId+1,'max',length(model.strTimeFcns)+1,'enable','on','string',str);
end

% Check if there are dynamic loads on node 1
if ~isempty(model.nodes(1).load.dynamic)
    f  = model.nodes(1).load.dynamic;
else
    f  = zeros(1,6);
end

% Switch actions according to anm
switch model.anm.analysis_type
    case TRUSS2D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(f(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(f(2)));
        set(handles.edit_Fz,'Enable','off');
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','off');
        
    case FRAME2D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(f(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(f(2)));
        set(handles.edit_Fz,'Enable','off');
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','on','String',num2str(f(6)));
        
    case GRILLAGE_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','off');
        set(handles.edit_Fy,'Enable','off');
        set(handles.edit_Fz,'Enable','on','String',num2str(f(3)));
        set(handles.edit_Mx,'Enable','on','String',num2str(f(4)));
        set(handles.edit_My,'Enable','on','String',num2str(f(5)));
        set(handles.edit_Mz,'Enable','off');
        
    case TRUSS3D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(f(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(f(2)));
        set(handles.edit_Fz,'Enable','on','String',num2str(f(3)));
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','off');
        
    case FRAME3D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(f(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(f(2)));
        set(handles.edit_Fz,'Enable','on','String',num2str(f(3)));
        set(handles.edit_Mx,'Enable','on','String',num2str(f(4)));
        set(handles.edit_My,'Enable','on','String',num2str(f(5)));
        set(handles.edit_Mz,'Enable','on','String',num2str(f(6)));
end

% Concentrated mass values
set(handles.edit_DisplMass,'Enable','on','String',num2str(model.nodes(1).displMass*1000));
% if model.anm.analysis_type == TRUSS2D_ANALYSIS || model.anm.analysis_type == TRUSS3D_ANALYSIS
%     set(handles.edit_RotMass,'Enable','off');
% else
%     set(handles.edit_RotMass,'Enable','on','String',num2str(model.nodes(1).rotMass));
% end

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(handles.axes_LFcn);
if ~fcnId
    x = [0, (model.t<=0)+model.t];
    y = [0, 0];
elseif ~isempty(model.timeFcns{fcnId})
    step = model.t/model.n_steps;
    x = 0:step:model.t;
    y = model.timeFcns{fcnId}.evalAll;
else
    x = [0, (model.t<=0)+model.t];
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

% Choose default command line output for GUI_NodalLoads_Dynamic
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_NodalLoads_Dynamic_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Node.
function popupmenu_Node_Callback(hObject, ~, handles) %#ok<*DEFNU>
include_constants;
            
% Get model object from root
model = getappdata(0,'model');

% Get popupmenu value
n = get(hObject,'Value');

% Set popupmenu_TimeFcns properties
fcnId = model.findFcnListByHandle(model.nodes(n).load.getFcn());
if isempty(model.strTimeFcns)
    set(handles.popupmenu_TimeFcns,'value',1,'max',1,'enable','off','string','None ');
else
    str = model.strTimeFcns; str(2:end+1) = str; str{1} = 'None';
    set(handles.popupmenu_TimeFcns,'value',fcnId+1,'max',length(model.strTimeFcns)+1,'enable','on','string',str);
end

% Check if there are dynamic loads on node n
if ~isempty(model.nodes(n).load.dynamic)
    ampl  = model.nodes(n).load.dynamic;
else
    ampl  = zeros(1,6);
end

% Switch actions according to anm
switch model.anm.analysis_type
    case TRUSS2D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(ampl(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(ampl(2)));
        set(handles.edit_Fz,'Enable','off');
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','off');
        
    case FRAME2D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(ampl(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(ampl(2)));
        set(handles.edit_Fz,'Enable','off');
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','on','String',num2str(ampl(6)));
        
    case GRILLAGE_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','off');
        set(handles.edit_Fy,'Enable','off');
        set(handles.edit_Fz,'Enable','on','String',num2str(ampl(3)));
        set(handles.edit_Mx,'Enable','on','String',num2str(ampl(4)));
        set(handles.edit_My,'Enable','on','String',num2str(ampl(5)));
        set(handles.edit_Mz,'Enable','off');
        
    case TRUSS3D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(ampl(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(ampl(2)));
        set(handles.edit_Fz,'Enable','on','String',num2str(ampl(3)));
        set(handles.edit_Mx,'Enable','off');
        set(handles.edit_My,'Enable','off');
        set(handles.edit_Mz,'Enable','off');
        
    case FRAME3D_ANALYSIS
        % Amplitude
        set(handles.edit_Fx,'Enable','on','String',num2str(ampl(1)));
        set(handles.edit_Fy,'Enable','on','String',num2str(ampl(2)));
        set(handles.edit_Fz,'Enable','on','String',num2str(ampl(3)));
        set(handles.edit_Mx,'Enable','on','String',num2str(ampl(4)));
        set(handles.edit_My,'Enable','on','String',num2str(ampl(5)));
        set(handles.edit_Mz,'Enable','on','String',num2str(ampl(6)));
end

%Concentrated mass values
set(handles.edit_DisplMass,'Enable','on','String',num2str(model.nodes(n).displMass*1000));
% if model.anm.analysis_type == TRUSS2D_ANALYSIS || model.anm.analysis_type == TRUSS3D_ANALYSIS
%     set(handles.edit_RotMass,'Enable','off');
% else
%     set(handles.edit_RotMass,'Enable','on','String',num2str(model.nodes(n).rotMass));
% end

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(handles.axes_LFcn);
if ~fcnId
    x = [0, (model.t<=0)+model.t];
    y = [0, 0];
elseif ~isempty(model.timeFcns{fcnId})
    step = model.t/model.n_steps;
    x = 0:step:model.t;
    y = model.timeFcns{fcnId}.evalAll;
else
    x = [0, (model.t<=0)+model.t];
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
% --- Executes on selection change in popupmenu_TimeFcns.
function popupmenu_TimeFcns_Callback(hObject, ~, handles) %#ok<*DEFNU>
            
% Get model object from root
model = getappdata(0,'model');

% Get fcn id
fcnId = get(hObject,'value')-1;

% Update axes
delete(findobj('Tag','DrawLfcn'))
axes(handles.axes_LFcn);
if ~fcnId
    x = [0, (model.t<=0)+model.t];
    y = [0, 0];
elseif ~isempty(model.timeFcns{fcnId})
    step = model.t/model.n_steps;
    x = 0:step:model.t;
    y = model.timeFcns{fcnId}.evalAll;
else
    x = [0, (model.t<=0)+model.t];
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
% --- Executes on button press in checkbox_MultiNodes.
function checkbox_MultiNodes_Callback(hObject, ~, handles)
if get(hObject,'value')
    set(hObject,'UserData',get(handles.popupmenu_Node,'value'));
    set(handles.popupmenu_Node,'enable','off','string',' ','value',1)
    set(handles.edit_MultiNodes,'enable','on','string','1')
    set(handles.popupmenu_MultiApply,'Enable','on','value',1);
else
    set(hObject,'UserData',[]);
    set(handles.popupmenu_Node,'enable','on','string',num2str(1:getappdata(0,'nnp'),'%d\n'),'value',1)
    popupmenu_Node_Callback(handles.popupmenu_Node, [], handles)
    set(handles.edit_MultiNodes,'enable','off','string','')
    set(handles.popupmenu_MultiApply,'Enable','off','value',1);
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Set.
function pushbutton_Set_Callback(hObject, ~, handles)
include_constants;
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
nodes = getappdata(0,'nodes');

% Disable button while input data is being set
set(hObject,'enable','off')

% Get ID of selected node
if ~get(handles.checkbox_MultiNodes,'value')
    n_ID = get(handles.popupmenu_Node,'Value');
    multiPopup = 0;
else
    n_str = get(handles.edit_MultiNodes,'string');
    multiPopup = get(handles.popupmenu_MultiApply,'value');
    if strcmp(n_str,'all')
        n_ID = 1:getappdata(0,'nnp');
    elseif strcmp(n_str,'ALL')
        n_ID = 1:getappdata(0,'nnp');
    elseif strcmp(n_str,'All')
        n_ID = 1:getappdata(0,'nnp');
    else
        [flag,n_ID] = readStr(n_str,getappdata(0,'nnp'));
        if ~flag
            % Enable button for future use
            set(hObject,'enable','on');
    
            msgbox('Invalid input data!', 'Error','error');
            return
        end
    end
end

if ~multiPopup
    % Get load amplitude values
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
        % Enable button for future use
        set(hObject,'enable','on')

        msgbox('Invalid input data!', 'Error','error');
        return
    end

    % Set load components
    nodes(n_ID).load.dynamic = [fx  , fy  , fz  , mx  , my  , mz];
    
    % Set concentrated displacement mass (kg to Ton)
    mass = str2double(get(handles.edit_DisplMass,'String'));
    if isnan(mass), mass=0; end
    if mass >= 0, nodes(n_ID).displMass = mass* 0.001; end
    
    % Set time fcn
    fcnId = get(handles.popupmenu_TimeFcns,'value')-1;
    fcn = [];
    if fcnId, fcn = model.timeFcns{fcnId}; end
    nodes(n_ID).load.setFcn(fcn);

elseif multiPopup == 2  % Set only load components
    % Get load amplitude values
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
        % Enable button for future use
        set(hObject,'enable','on')

        msgbox('Invalid input data!', 'Error','error');
        return
    end

    % Loop through nodes that will have new assigned loads
    for n = n_ID
        % Set load components
        nodes(n).load.dynamic = [fx  , fy  , fz  , mx  , my  , mz];
    end
    
    % Get concentrated displacement mass (kg to Ton)
    mass = str2double(get(handles.edit_DisplMass,'String'));
    if isnan(mass), mass=0; end
    mass = mass*0.001;
    for n = n_ID
        if mass >= 0, nodes(n).displMass = mass; end
    end
    
elseif multiPopup == 1 % Set load components and functions
    % Get load amplitude values
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
        % Enable button for futre use
        set(hObject,'enable','on')

        msgbox('Invalid input data!', 'Error','error');
        return
    end

    % Get fcn id
    fcnId = get(handles.popupmenu_TimeFcns,'value')-1;
    fcn = [];
    if fcnId, fcn = model.timeFcns{fcnId}; end
    
    % Loop through nodes that will have new assigned loads
    for n = n_ID
        % Set load components
        nodes(n).load.dynamic = [fx  , fy  , fz  , mx  , my  , mz];
        % Set load time functions
        nodes(n).load.setFcn(fcn);
    end
    
    % Get concentrated displacement mass (kg to Ton)
    mass = str2double(get(handles.edit_DisplMass,'String'));
    if isnan(mass), mass=0; end
    mass = mass*0.001;
    for n = n_ID
        if mass >= 0, nodes(n).displMass = mass; end
    end

%     % Get concentrated rotational inertia (kg to Ton)
%     mass = str2double(get(handles.edit_RotMass,'String')) * 0.001;
%     for n = n_ID
%         if mass > 0, nodes(n).rotMass = mass; end
%     end

else % Set only time functions
    
    % Get fcn id
    fcnId = get(handles.popupmenu_TimeFcns,'value')-1;
    fcn = [];
    if fcnId, fcn = model.timeFcns{fcnId}; end
    
    % Loop through nodes that will have new assigned loads
    for n = n_ID
        % Set load time functions
        nodes(n).load.setFcn(fcn);
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
gui = findobj('Tag','GUI_NodalLoads_Dynamic');
set(gui,'WindowStyle','normal');

% % Draw updated model
% if allLoadsNeedToBeDrawn
%     redraw(mdata,'Loads')
% end

% Draw updated model
if allLoadsNeedToBeDrawn == false
    redraw(mdata,'Nodal Loads')
else
    redraw(mdata,'Loads')
end

% Make GUI a modal window
set(gui,'WindowStyle','modal');

% Reset checkbox (RLR: Why? I will remove it)
% if get(handles.checkbox_MultiNodes,'value')
%     set(handles.checkbox_MultiNodes,'Value',0);
%     checkbox_MultiNodes_Callback(handles.checkbox_MultiNodes,[],handles);
% end

% Enable button for futre use
set(hObject,'enable','on')

%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function popupmenu_Node_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_TimeFcns_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_MultiNodes_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_MultiNodes_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Fx_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Fx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Fy_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Fy_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Fz_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Fz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Mx_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Mx_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_My_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_My_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Mz_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function edit_Mz_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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