%% Cross-Section Dialog Callback Functions
% This file contains the callback functions associated with the "Sections"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Sections(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Sections_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Sections_OutputFcn, ...
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
% Executes just before Sections GUI is made visible.
% Sets GUI initial properties.
function GUI_Sections_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Sections
handles.output = hObject;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Set default editable text boxes strings and table data
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
nsec = getappdata(0,'nsec');
if nsec > 0
    sections = getappdata(0,'sections');
    sec = cell(nsec,9);
    if anm ==1
        for s = 1:nsec
            id = sections(s).id;
            ax = 1e4*sections(s).area_x;
            hy = 1e2*sections(s).height_y;
            sec(s,:) = {id,ax,'          -','          -','          -','          -','          -',hy,'          -'};
            cEdit = [false(1,1) true(1,1) false(1,5) true(1,1) false(1,2)];
        end
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Hy,'Enable','on','String',10)
    elseif anm ==2
        for s = 1:nsec
            id = sections(s).id;
            ax = 1e4*sections(s).area_x;
            ay = 1e4*sections(s).area_y;
            iz = 1e8*sections(s).inertia_z;
            hy = 1e2*sections(s).height_y;
            sec(s,:) = {id,ax,ay,'          -','          -','          -',iz,hy,'          -'};
            cEdit = [false(1,1) true(1,2) false(1,3) true(1,2) false(1,2)];
        end
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Ay,'Enable','on','String',100)
        set(handles.edit_Iz,'Enable','on','String',1000)
        set(handles.edit_Hy,'Enable','on','String',10)
    elseif anm == 3
        for s = 1:nsec
            id = sections(s).id;
            ax = 1e4*sections(s).area_x;
            az = 1e4*sections(s).area_z;
            ix = 1e8*sections(s).inertia_x;
            iy = 1e8*sections(s).inertia_y;
            hz = 1e2*sections(s).height_z;
            sec(s,:) = {id,ax,'          -',az,ix,iy,'          -','          -',hz};
            cEdit = [false(1,1) true(1,1) false(1,1) true(1,3) false(1,2) true(1,1) false(1,1)];
        end
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Az,'Enable','on','String',100)
        set(handles.edit_Ix,'Enable','on','String',1000)
        set(handles.edit_Iy,'Enable','on','String',1000)
        set(handles.edit_Hz,'Enable','on','String',10)
    elseif anm == 4
        for s = 1:nsec
            id = sections(s).id;
            ax = 1e4*sections(s).area_x;
            hy = 1e2*sections(s).height_y;
            hz = 1e2*sections(s).height_z;
            sec(s,:) = {id,ax,'          -','          -','          -','          -','          -',hy,hz};
            cEdit = [false(1,1) true(1,1) false(1,5) true(1,2) false(1,1)];
        end
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Hy,'Enable','on','String',10)
        set(handles.edit_Hz,'Enable','on','String',10)
    elseif anm == 5
        for s = 1:nsec
            id = sections(s).id;
            ax = 1e4*sections(s).area_x;
            ay = 1e4*sections(s).area_y;
            az = 1e4*sections(s).area_z;
            ix = 1e8*sections(s).inertia_x;
            iy = 1e8*sections(s).inertia_y;
            iz = 1e8*sections(s).inertia_z;
            hy = 1e2*sections(s).height_y;
            hz = 1e2*sections(s).height_z;
            sec(s,:) = {id,ax,ay,az,ix,iy,iz,hy,hz};
            cEdit = [false(1,1) true(1,8) false(1,1)];
        end
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Ay,'Enable','on','String',100)
        set(handles.edit_Az,'Enable','on','String',100)
        set(handles.edit_Ix,'Enable','on','String',1000)
        set(handles.edit_Iy,'Enable','on','String',1000)
        set(handles.edit_Iz,'Enable','on','String',1000)
        set(handles.edit_Hy,'Enable','on','String',10)
        set(handles.edit_Hz,'Enable','on','String',10)
    end
    set(handles.uitable_Sections,'Data',sec,'CellEditCallback',@uitable_Sections_CellEditCallback,'ColumnEditable',cEdit);
else
    if anm ==1
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Hy,'Enable','on','String',10)
    elseif anm ==2
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Ay,'Enable','on','String',100)
        set(handles.edit_Iz,'Enable','on','String',1000)
        set(handles.edit_Hy,'Enable','on','String',10)
    elseif anm == 3
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Az,'Enable','on','String',100)
        set(handles.edit_Ix,'Enable','on','String',1000)
        set(handles.edit_Iy,'Enable','on','String',1000)
        set(handles.edit_Hz,'Enable','on','String',10)
    elseif anm == 4
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Hy,'Enable','on','String',10)
        set(handles.edit_Hz,'Enable','on','String',10)
    elseif anm == 5
        set(handles.edit_Ax,'Enable','on','String',100)
        set(handles.edit_Ay,'Enable','on','String',100)
        set(handles.edit_Az,'Enable','on','String',100)
        set(handles.edit_Ix,'Enable','on','String',1000)
        set(handles.edit_Iy,'Enable','on','String',1000)
        set(handles.edit_Iz,'Enable','on','String',1000)
        set(handles.edit_Hy,'Enable','on','String',10)
        set(handles.edit_Hz,'Enable','on','String',10)
    end
    set(handles.uitable_Sections,'Data',{'','','','','','','','',''},'CellEditCallback',@uitable_Sections_CellEditCallback);
end

% Create list of cross-sections to be deleted
if nsec > 0
    secid =zeros(1,nsec);
    for s = 1:nsec
        secid(s) = s;
    end
    secid = num2str(secid,'%d\n');
    set(handles.popupmenu_Delete,'Value',nsec,'string',secid)
else
    set(handles.popupmenu_Delete,'Value',1,'string','No itens')
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Sections_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes on button press in "Add" pushbutton.
% Adds a Section object with input properties to the list of cross-sections
function pushbutton_Add_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
Ax = 1e-4*str2double(get(handles.edit_Ax,'String'));
Ay = 1e-4*str2double(get(handles.edit_Ay,'String'));
if isnan(Ay)
    Ay = 1;
end
Az = 1e-4*str2double(get(handles.edit_Az,'String'));
if isnan(Az)
    Az = 1;
end
Ix = 1e-8*str2double(get(handles.edit_Ix,'String'));
if isnan(Ix)
    Ix = 1;
end
Iy = 1e-8*str2double(get(handles.edit_Iy,'String'));
if isnan(Iy)
    Iy = 1;
end
Iz = 1e-8*str2double(get(handles.edit_Iz,'String'));
if isnan(Iz)
    Iz = 1;
end
Hy = 1e-2*str2double(get(handles.edit_Hy,'String'));
if isnan(Hy)
    Hy = 1;
end
Hz = 1e-2*str2double(get(handles.edit_Hz,'String'));
if isnan(Hz)
    Hz = 1;
end

if (Ax > 0) && (Ay > 0) && (Az > 0) &&...
   (Ix > 0) && (Iy > 0) && (Iz > 0) &&...
   (Hy > 0) && (Hz > 0)
    % Increment number of cross-sections and create a Section object
    nsec = getappdata(0,'nsec');
    nsec = nsec + 1;
    s = Section(nsec,Ax,Ay,Az,Ix,Iy,Iz,Hy,Hz);
    
    % Insert created Section object in a vector of cross-sections
    if nsec ~= 1
        sections = getappdata(0,'sections');
    end
    sections(nsec) = s;
    
    % Update sections uitable
    model = getappdata(0,'model'); % get handle to model object
    include_constants
    if model.anm.analysis_type == TRUSS2D_ANALYSIS
        id = s.id;
        ax = 1e4*s.area_x;
        hy = 1e2*s.height_y;
        newSec = {id,ax,'          -','          -','          -','          -','          -',hy,'          -'};
        cEdit = [false(1,1) true(1,1) false(1,5) true(1,1) false(1,2)];
    elseif model.anm.analysis_type == FRAME2D_ANALYSIS
        id = s.id;
        ax = 1e4*s.area_x;
        ay = 1e4*s.area_y;
        iz = 1e8*s.inertia_z;
        hy = 1e2*s.height_y;
        newSec = {id,ax,ay,'          -','          -','          -',iz,hy,'          -'};
        cEdit = [false(1,1) true(1,2) false(1,3) true(1,2) false(1,2)];
    elseif model.anm.analysis_type == GRILLAGE_ANALYSIS
        id = s.id;
        ax = 1e4*s.area_x;
        az = 1e4*s.area_z;
        ix = 1e8*s.inertia_x;
        iy = 1e8*s.inertia_y;
        hz = 1e2*s.height_z;
        newSec = {id,ax,'          -',az,ix,iy,'          -','          -',hz};
        cEdit = [false(1,1) true(1,1) false(1,1) true(1,3) false(1,2) true(1,1) false(1,1)];
    elseif model.anm.analysis_type == TRUSS3D_ANALYSIS
        id = s.id;
        ax = 1e4*s.area_x;
        hy = 1e2*s.height_y;
        hz = 1e2*s.height_z;
        newSec = {id,ax,'          -','          -','          -','          -','          -',hy,hz};
        cEdit = [false(1,1) true(1,1) false(1,5) true(1,2) false(1,1)];
    elseif model.anm.analysis_type == FRAME3D_ANALYSIS
        id = s.id;
        ax = 1e4*s.area_x;
        ay = 1e4*s.area_y;
        az = 1e4*s.area_z;
        ix = 1e8*s.inertia_x;
        iy = 1e8*s.inertia_y;
        iz = 1e8*s.inertia_z;
        hy = 1e2*s.height_y;
        hz = 1e2*s.height_z;
        newSec = {id,ax,ay,az,ix,iy,iz,hy,hz};
        cEdit = [false(1,1) true(1,8) false(1,1)];
    end
    tableData = get(handles.uitable_Sections,'Data');
    if isempty(tableData)
        set(handles.uitable_Sections,'Data',newSec,'ColumnEditable',cEdit)
    elseif strcmp(tableData{1},'')
        set(handles.uitable_Sections,'Data',newSec,'ColumnEditable',cEdit)
    else
        set(handles.uitable_Sections,'Data',vertcat(tableData,newSec),'ColumnEditable',cEdit)
    end
    clear s
    
    % Update list of sections to be deleted
    sec_id = zeros(1,nsec);
    for s = 1:nsec
        sec_id(s) = s;
    end
    sec_id = num2str(sec_id,'%d\n');
    set(handles.popupmenu_Delete,'Value',nsec,'string',sec_id,'Max',nsec)
    
    % Set model object properties
    model.sections = sections;
    model.nsec = nsec;
    
    % Update information panel in GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(2,:) = {'Cross-Sections',nsec};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Update mouse property
    mouse = getappdata(0,'mouse');
    if ~isempty(mouse.originalData)
        mouse.originalData = infoPanelData;
    end
    setappdata(0,'mouse',mouse)
    
    % Return variables to root
    setappdata(0,'model',model);
    setappdata(0,'sections',sections);
    setappdata(0,'nsec',nsec);
    setappdata(0,'move',0);
else
    msgbox('Invalid input data!', 'Error','error');
end

%--------------------------------------------------------------------------
% Executes on button press in "Cancel" pushbutton.
% Returns to main GUI without creating any Section object.
function pushbutton_Cancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
delete(gcf)

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Delete.
% Deletes a Section object from the list of cross-sections.
function pushbutton_Delete_Callback(hObject, eventdata, handles)
model = getappdata(0,'model');
sections = getappdata(0,'sections');
nsec = getappdata(0,'nsec');
elems = getappdata(0,'elems');
nel = getappdata(0,'nel');

if ~strcmp(get(handles.popupmenu_Delete,'String'),'No itens')
    % Get ID of deleted cross-section
    del_s = get(handles.popupmenu_Delete,'Value');
    
    % Check if deleted cross-section is being used by any element
    sec_use = 0;
    for e = 1:nel
        if elems(e).section.id == sections(del_s).id
            sec_use = sec_use + 1;
        end
    end
    
    if sec_use ~= 0
        msgbox('This cross-section cannot be deleted because it is being used by elements.', 'Error','error');
        return
    else
        % Remove deleted cross-section from vector of cross-sections
        sections(del_s) = [];
        
        % Update number of cross-sections
        nsec = nsec - 1;
        
        % Update list of sections to be deleted and uitable data
        if nsec ~= 0
            tableData = get(handles.uitable_Sections,'Data');
            tableData(del_s,:) = [];
            sec_id = zeros(1,nsec);
            for s = 1:nsec
                tableData(s,1) = {s};
                sections(s).id = s;
                sec_id(s) = s;
            end
            sec_id = num2str(sec_id,'%d\n');
            set(handles.popupmenu_Delete,'value',nsec,'string',sec_id,...
                'Max',nsec)
        else
            set(handles.popupmenu_Delete,'Value',1,'string','No itens',...
                'Max',1)
        end
        
        % Set updated sections uitable data
        if nsec ~= 0
            set(handles.uitable_Sections,'Data',tableData)
        else
            set(handles.uitable_Sections,'Data',{'','','','','','','','',''})
        end
        
        % Set model object properties
        model.sections = sections;
        model.nsec = nsec;
        
        % Update information panel in GUI_Main
        mdata = guidata(findobj('Tag','GUI_Main'));
        infoPanelData = get(mdata.uitable_infoPanel,'Data');
        infoPanelData(2,:) = {'Cross-Sections',nsec};
        set(mdata.uitable_infoPanel,'Data',infoPanelData)
        
        % Update mouse property
        mouse = getappdata(0,'mouse');
        if ~isempty(mouse.originalData)
            mouse.originalData = infoPanelData;
        end
        setappdata(0,'mouse',mouse)
        
        % Return variables to root
        setappdata(0,'nsec',nsec);
        setappdata(0,'sections',sections);
        setappdata(0,'elems',elems);
        setappdata(0,'nel',nel);
        setappdata(0,'model',model);
        setappdata(0,'move',0);
    end
end

%--------------------------------------------------------------------------
% Executes during edit_Ax creation, after setting all properties.
function edit_Ax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Ax_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Ay creation, after setting all properties.
function edit_Ay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 3) || (anm == 4)
    set(hObject,'Enable','off','String','');
end

function edit_Ay_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Az creation, after setting all properties.
function edit_Az_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','');
end

function edit_Az_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Ix creation, after setting all properties.
function edit_Ix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','');
end

function edit_Ix_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Iy creation, after setting all properties.
function edit_Iy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 4)
    set(hObject,'Enable','off','String','');
end

function edit_Iy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Iz creation, after setting all properties.
function edit_Iz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 3) || (anm == 4)
    set(hObject,'Enable','off','String','');
end

function edit_Iz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Hy creation, after setting all properties.
function edit_Hy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 3)
    set(hObject,'Enable','off','String','');
end

function edit_Hy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Hz creation, after setting all properties.
function edit_Hz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2)
    set(hObject,'Enable','off','String','');
end

function edit_Hz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Delete creation, after setting all properties.
function popupmenu_Delete_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_Delete_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_Sections
function uitable_Sections_CellEditCallback(hObject, eventdata, handles)
mdata = guidata(findobj('tag','GUI_Sections'));
set(mdata.pushbutton_Apply,'Enable','on')

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(hObject, eventdata, handles)
% Get sections info and uitable data
nsec = getappdata(0,'nsec');
sections = getappdata(0,'sections');
nel = getappdata(0,'nel');
elems = getappdata(0,'elems');
tableData = get(handles.uitable_Sections,'Data');

% Get handle to gui_main
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'value');

% Check which sections are in use. (0 = not used, 1 = in use)
secsInUse = zeros(1,nsec);
for e = 1:nel
    id = elems(e).section.id;
    secsInUse(id) = 1;
end

% Get new info and set it to materials
sectionChangedFlag = 0; % flag for changes in section properties
for ns = 1:nsec
    % Ax
    newAx = 1e-4 * cell2mat(tableData(ns,2));
    if ~isnan(newAx) && newAx > 0 && newAx ~= sections(ns).area_x
        sections(ns).area_x = newAx;
        if secsInUse(ns) == 1
            sectionChangedFlag = 1;
        end
    end
    tableData(ns,2) = {1e4 * sections(ns).area_x};
    
    % Ay
    if anm == 2 || anm == 5
        newAy = 1e-4 * cell2mat(tableData(ns,3));
        if ~isnan(newAy) && newAy > 0 && newAy ~= sections(ns).area_y
            sections(ns).area_y = newAy;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,3) = {1e4 * sections(ns).area_y};
    end
    
    % Az
    if anm == 3 || anm == 5
        newAz = 1e-4 * cell2mat(tableData(ns,4));
        if ~isnan(newAz) && newAz > 0 && newAz ~= sections(ns).area_z
            sections(ns).area_z = newAz;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,4) = {1e4 * sections(ns).area_z};
    end
    
    % Ix
    if anm == 3 || anm == 5
        newIx = 1e-8 * cell2mat(tableData(ns,5));
        if ~isnan(newIx) && newIx > 0 && newIx ~= sections(ns).inertia_x
            sections(ns).inertia_x = newIx;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,5) = {1e8 * sections(ns).inertia_x};
    end
    
    % Iy
    if anm ==3 || anm == 5
        newIy = 1e-8 * cell2mat(tableData(ns,6));
        if ~isnan(newIy) && newIy > 0 && newIy ~= sections(ns).inertia_y
            sections(ns).inertia_y = newIy;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,6) = {1e8 * sections(ns).inertia_y};
    end
    
    % Iz
    if anm == 2 || anm ==5
        newIz = 1e-8 * cell2mat(tableData(ns,7));
        if ~isnan(newIz) && newIz > 0 && newIz ~= sections(ns).inertia_z
            sections(ns).inertia_z = newIz;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,7) = {1e8 * sections(ns).inertia_z};
    end
    
    % hy
    if anm ~= 3
        newHy = 1e-2 * cell2mat(tableData(ns,8));
        if ~isnan(newHy) && newHy > 0 && newHy ~= sections(ns).height_y
            sections(ns).height_y = newHy;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,8) = {1e2 * sections(ns).height_y};
    end
    
    % hz
    if anm == 3 || anm == 4 || anm == 5
        newHz = 1e-2 * cell2mat(tableData(ns,9));
        if ~isnan(newHz) && newHz > 0 && newHz ~= sections(ns).height_z
            sections(ns).height_z = newHz;
            if secsInUse(ns) == 1
                sectionChangedFlag = 1;
            end
        end
        tableData(ns,9) = {1e2 * sections(ns).height_z};
    end
end

% Update uitable
set(handles.uitable_Sections,'Data',tableData)

% Update model object
model = getappdata(0,'model');
model.sections = sections;

% Update variables in root
setappdata(0,'model',model)
setappdata(0,'sections',sections)

% Disable apply pushbutton
set(hObject,'enable','off')

% Disable results and enable process data pushbutton
if sectionChangedFlag ~= 0
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
