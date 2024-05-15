%% Material Dialog Callback Functions
% This file contains the callback functions associated with the "Materials"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Materials(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Materials_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Materials_OutputFcn, ...
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
% Executes just before Materials GUI is made visible.
% Sets GUI initial properties.
function GUI_Materials_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Materials
handles.output = hObject;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Set table data
nmat = getappdata(0,'nmat');
if nmat > 0
    materials = getappdata(0,'materials');
    mat = cell(nmat,6);
    for m = 1:nmat
        id = materials(m).id;
        E = 1e-3*materials(m).elasticity;
        v = materials(m).poisson;
        G = 1e-3*materials(m).shear;
        alpha = materials(m).thermExp;
        rho = 1e3*materials(m).density;
        mat(m,:) = {id,E,v,G,alpha,rho};
    end
    set(handles.uitable_Materials,'Data',mat,'CellEditCallback',@uitable_Materials_CellEditCallback);
else
    set(handles.uitable_Materials,'Data',{'','','','','',''},'CellEditCallback',@uitable_Materials_CellEditCallback);
end

set(handles.pushbutton_Apply,'Enable','off')

% Create list of materials to be deleted
if nmat > 0
    matid =zeros(1,nmat);
    for m = 1:nmat
        matid(m) = m;
    end
    matid = num2str(matid,'%d\n');
    set(handles.popupmenu_Delete,'Value',nmat,'string',matid)
else
    set(handles.popupmenu_Delete,'Value',1,'string','No itens')
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Materials_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes on button press in "Add" pushbutton.
% Adds a Material object with input properties to the list of materials
function pushbutton_Add_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
e = 1e3*str2double(get(handles.edit_Elasticity,'String'));
v = str2double(get(handles.edit_Poisson,'String'));
te = str2double(get(handles.edit_Thermal,'String'));
rho = 1e-3*str2double(get(handles.edit_Rho,'String'));

if (e > 0) && (v < 0.5) && (v > 0) && ~isnan(te) && ~isnan(rho)
    % Increment number of materials and create a Material object
    nmat = getappdata(0,'nmat');
    nmat = nmat + 1;
    m = Material(nmat,e,v,te,rho);
    
    % Insert created Material object in a vector of materials
    if nmat ~= 1
        materials = getappdata(0,'materials');
    end
    materials(nmat) = m;
    
    % Set model object properties
    model = getappdata(0,'model');
    model.materials = materials;
    model.nmat = nmat;
    
    % Update materials uitable
    id = m.id;
    E = 1e-3*m.elasticity;
    v = m.poisson;
    G = 1e-3*m.shear;
    alpha = m.thermExp;
    rho = 1e3*m.density;
    newMat = {id,E,v,G,alpha,rho};
    tableData = get(handles.uitable_Materials,'Data');
    if isempty(tableData)
        set(handles.uitable_Materials,'Data',newMat)
    elseif strcmp(tableData{1},'')
        set(handles.uitable_Materials,'Data',newMat)
    else
        set(handles.uitable_Materials,'Data',vertcat(tableData,newMat))
    end
    clear m
    
    % Update list of materials to be deleted
    mat_id = zeros(1,nmat);
    for m = 1:nmat
        mat_id(m) = m;
    end
    mat_id = num2str(mat_id,'%d\n');
    set(handles.popupmenu_Delete,'Value',nmat,'string',mat_id,'Max',nmat)
    
    % Update information panel in GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(1,:) = {'Materials',nmat};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Update mouse property
    mouse = getappdata(0,'mouse');
    if ~isempty(mouse.originalData)
        mouse.originalData = infoPanelData;
    end
    setappdata(0,'mouse',mouse)
    
    % Return variables to root
    setappdata(0,'model',model);
    setappdata(0,'materials',materials);
    setappdata(0,'nmat',nmat);
    setappdata(0,'move',0);
else
    msgbox('Invalid input data!', 'Error','error');
end

%--------------------------------------------------------------------------
% Executes on button press in "Cancel" pushbutton.
% Returns to main GUI without creating any Material object.
function pushbutton_Cancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
delete(gcf)

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Delete.
% Deletes a Material object from the list of materials.
function pushbutton_Delete_Callback(hObject, eventdata, handles)
model = getappdata(0,'model');
materials = getappdata(0,'materials');
nmat = getappdata(0,'nmat');
elems = getappdata(0,'elems');
nel = getappdata(0,'nel');

if ~strcmp(get(handles.popupmenu_Delete,'String'),'No itens')
    % Get ID of deleted material
    del_m = get(handles.popupmenu_Delete,'Value');
    
    % Check if deleted material is being used by any element
    mat_use = false;
    for e = 1:nel
        if elems(e).material.id == materials(del_m).id
            mat_use = true;
            break
        end
    end
    
    if mat_use == true
        msgbox('This material cannot be deleted because it is being used by elements.', 'Error','error');
        return
    else
        % Remove deleted material from vector of materials
        materials(del_m) = [];
        
        % Update number of materials
        nmat = nmat - 1;
        
        % Update list of materials to be deleted and uitable data
        if nmat ~= 0
            tableData = get(handles.uitable_Materials,'Data');
            tableData(del_m,:) = [];
            mat_id = zeros(1,nmat);
            for m = 1:nmat
                tableData(m,1) = {m};
                materials(m).id = m;
                mat_id(m) = m;
            end
            mat_id = num2str(mat_id,'%d\n');
            set(handles.popupmenu_Delete,'value',nmat,'string',mat_id,...
                'Max',nmat)
        else
            set(handles.popupmenu_Delete,'Value',1,'string','No itens',...
                'Max',1)
        end
        
        % Set updated materials uitable data
        if nmat ~= 0
            set(handles.uitable_Materials,'Data',tableData)
        else
            set(handles.uitable_Materials,'Data',{'','','','','',''})
        end
        
        % Set model object properties
        model.materials = materials;
        model.nmat = nmat;
        
        % Update information panel in GUI_Main
        mdata = guidata(findobj('Tag','GUI_Main'));
        infoPanelData = get(mdata.uitable_infoPanel,'Data');
        infoPanelData(1,:) = {'Materials',nmat};
        set(mdata.uitable_infoPanel,'Data',infoPanelData)
        
        % Update mouse property
        mouse = getappdata(0,'mouse');
        if ~isempty(mouse.originalData)
            mouse.originalData = infoPanelData;
        end
        setappdata(0,'mouse',mouse)
        
        % Return variables to root
        setappdata(0,'nmat',nmat);
        setappdata(0,'materials',materials);
        setappdata(0,'elems',elems);
        setappdata(0,'nel',nel);
        setappdata(0,'model',model);
        setappdata(0,'move',0);
    end
end

%--------------------------------------------------------------------------
% Executes during edit_Elasticity creation, after setting all properties.
function edit_Elasticity_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Elasticity_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Poisson creation, after setting all properties.
function edit_Poisson_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Poisson_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Thermal creation, after setting all properties.
function edit_Thermal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Thermal_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Delete creation, after setting all properties.
function popupmenu_Delete_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_Delete_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_Materials
function uitable_Materials_CellEditCallback(hObject, eventdata, handles)
mdata = guidata(findobj('tag','GUI_Materials'));
set(mdata.pushbutton_Apply,'Enable','on')

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(hObject, eventdata, handles)
% Get materials info and uitable data
nmat = getappdata(0,'nmat');
materials = getappdata(0,'materials');
nel = getappdata(0,'nel');
elems = getappdata(0,'elems');
tableData = get(handles.uitable_Materials,'Data');

% Check which materials are in use. (0 = not used, 1 = in use)
matsInUse = zeros(1,nmat);
for e = 1:nel
    id = elems(e).material.id;
    matsInUse(id) = 1;
end

% Get new info and set it to materials
materialChangedFlag = 0; % flag for changes in material properties
for nm = 1:nmat
    shearFlag = 0; % flag for changes in E or v (new G must be calculated)
    
    % E
    newE = 1e3 * tableData{nm,2};
    if ~isnan(newE) && newE > 0 && newE ~= materials(nm).elasticity 
        materials(nm).elasticity = newE;
        shearFlag = 1;
        if matsInUse(nm) == 1
            materialChangedFlag = 1;
        end
    end
    tableData(nm,2) = {1e-3 * materials(nm).elasticity};
    
    % v
    newPoisson = tableData{nm,3};
    if ~isnan(newPoisson) && newPoisson < 0.5 && newPoisson > 0 && newPoisson ~= materials(nm).poisson
        materials(nm).poisson = newPoisson;
        shearFlag = 1;
        if matsInUse(nm) == 1
            materialChangedFlag = 1;
        end
    end
    tableData(nm,3) = {materials(nm).poisson};
    
    % G
    if shearFlag ~= 0
        newShear = materials(nm).elasticity / (2 * (1 + materials(nm).poisson));
        materials(nm).shear = newShear;
        tableData(nm,4) = {1e-3 * materials(nm).shear};
    end
    
    % Alpha
    newAlpha = tableData{nm,5};
    if ~isnan(newAlpha) && newAlpha ~= materials(nm).thermExp
        materials(nm).thermExp = newAlpha;
        if matsInUse(nm) == 1
            materialChangedFlag = 1;
        end
    end
    tableData(nm,5) = {materials(nm).thermExp};
    
    % Rho
    newRho = 1e-3 * tableData{nm,6};
    if ~isnan(newRho) && newRho ~= materials(nm).density
        materials(nm).density = newRho;
        if matsInUse(nm) == 1
            materialChangedFlag = 1;
        end
    end
    tableData(nm,6) = {1e3 * materials(nm).density};
end

% Update uitable
set(handles.uitable_Materials,'Data',tableData)

% Update model object
model = getappdata(0,'model');
model.materials = materials;

% Update variables in root
setappdata(0,'model',model)
setappdata(0,'materials',materials)

% Disable apply pushbutton
set(hObject,'enable','off')

% Disable results and enable process data pushbutton
if materialChangedFlag ~= 0
    mdata = guidata(findobj('tag','GUI_Main'));
    resType = get(mdata.popupmenu_Results,'value');
    
    % Disable result buttons
    set(mdata.popupmenu_Results,'Enable','off','value',1);
    set(mdata.pushbutton_Textual,'Enable','off');
    set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
    set(mdata.text_Element,'string','Elements');
    set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(mdata.edit_Scale,'enable','off','visible','off');
    set(mdata.pushbutton_DynamicResults,'enable','off');
    set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    % Enable process data
    set(mdata.pushbutton_ProcessData,'Enable','on');
    
    % Redraw
    if resType ~= 1
        redraw(mdata,'Loads')
    end
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

%ronald
function pushbutton_MaterialList_Callback(hObject, eventdata, handles)
list=[];


for i=1:1:length(Steel.steelList)
   
list = [list,(Steel.steelList{i}{1,1})]; %function for load a steel list

end


[indx,tf] = listdlg('SelectionMode','single','ListString',list,'Name','Material List','ListSize',[240 300]);%put a steel list in list and receive a select index response
if isempty(indx)
    
else
    % Increment number of materials and create a Material object
    nmat = getappdata(0,'nmat');
    nmat = nmat + 1;
    e=(Steel.steelList{indx}{1,2});
    v=(Steel.steelList{indx}{1,3});
    te=(Steel.steelList{indx}{1,5});
    rho=(Steel.steelList{indx}{1,6});
    leak=(Steel.steelList{indx}{1,7});
    yld=(Steel.steelList{indx}{1,8});
    m = Steel(nmat,e,v,te,rho,leak,yld);
    
    % Insert created Material object in a vector of materials
    if nmat ~= 1
        materials = getappdata(0,'materials');
    end
    materials(nmat) = m;
    
    % Set model object properties
    model = getappdata(0,'model');
    model.materials = materials;
    model.nmat = nmat;
    
    % Update materials uitable
    id = m.id;
    E = 1e-3*m.elasticity;
    v = m.poisson;
    G = 1e-3*m.shear;
    alpha = m.thermExp;
    rho = 1e3*m.density;
    newMat = {id,E,v,G,alpha,rho};
    tableData = get(handles.uitable_Materials,'Data');
    if isempty(tableData)
        set(handles.uitable_Materials,'Data',newMat)
    elseif strcmp(tableData{1},'')
        set(handles.uitable_Materials,'Data',newMat)
    else
        set(handles.uitable_Materials,'Data',vertcat(tableData,newMat))
    end
    clear m
    
    % Update list of materials to be deleted
    mat_id = zeros(1,nmat);
    for m = 1:nmat
        mat_id(m) = m;
    end
    mat_id = num2str(mat_id,'%d\n');
    set(handles.popupmenu_Delete,'Value',nmat,'string',mat_id,'Max',nmat)
    
    % Update information panel in GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(1,:) = {'Materials',nmat};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Update mouse property
    mouse = getappdata(0,'mouse');
    if ~isempty(mouse.originalData)
        mouse.originalData = infoPanelData;
    end
    setappdata(0,'mouse',mouse)
    
    % Return variables to root
    setappdata(0,'model',model);
    setappdata(0,'materials',materials);
    setappdata(0,'nmat',nmat);
    setappdata(0,'move',0);
end