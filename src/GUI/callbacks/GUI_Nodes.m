%% Node Dialog Callback Functions
% This file contains the callback functions associated with the "Nodes"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Nodes(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Nodes_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Nodes_OutputFcn, ...
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
% Executes just before Nodes GUI is made visible.
% Sets GUI initial properties.
function GUI_Nodes_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Nodes
handles.output = hObject;

% Move GUI to the center of the screen
if getappdata(0,'move') == 1
    movegui(hObject,'center');
end

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Set table data
nnp = getappdata(0,'nnp');
if nnp > 0
    nodes = getappdata(0,'nodes');
    node = cell(nnp,4);
    for n = 1:nnp
        id = nodes(n).id;
        x = nodes(n).coord(1);
        y = nodes(n).coord(2);
        z = nodes(n).coord(3);
        node(n,:) = {id,x,y,z};
    end
    set(handles.uitable_Nodes,'Data',node);
else
    set(handles.uitable_Nodes,'Data',{'','','',''});
end

% Enable z coordinates for 3D models
model = getappdata(0,'model');
if model.anm.analysis_type >= 3
    set(handles.edit_Z,'enable','on')
end

% Create list of nodes to be deleted
if nnp ~= 0
    nd = zeros(1,nnp);
    for n = 1:nnp
        nd(n) = n;
    end
    nd = num2str(nd,'%d\n');
    set(handles.popupmenu_Delete,'Value',nnp,'string',nd)
else
    set(handles.popupmenu_Delete,'Value',1,'string','No itens')
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Nodes_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes on button press in "Add" pushbutton.
% Adds a Node object with input properties to the list of nodes
function pushbutton_Add_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
include_constants;

% Check if input coordinates are not the same of a previously created node
equal = 0;
nnp = getappdata(0,'nnp');
if nnp ~= 0
    nodes = getappdata(0,'nodes');
    for i = 1:nnp
        xi = str2double(get(handles.edit_X,'String'));
        yi = str2double(get(handles.edit_Y,'String'));
        zi = str2double(get(handles.edit_Z,'String'));
        xj = nodes(i).coord(1);
        yj = nodes(i).coord(2);
        zj = nodes(i).coord(3);
        if (xi == xj) && (yi == yj) && (zi == zj)
            equal = 1;
        end
    end
end

if equal == 0
    % Get coordinates
    x = str2double(get(handles.edit_X,'String'));
    y = str2double(get(handles.edit_Y,'String'));
    z = str2double(get(handles.edit_Z,'String'));
    
    if (isnan(x)) || (isnan(y)) || (isnan(z))
        msgbox('Invalid input data!', 'Error','error');
    else
        % Get graphical tolerance
        draw = getappdata(0,'draw');
        if draw.size == 0
            draw.setSize();
        end
        if draw.mdl.anm.analysis_type == TRUSS2D_ANALYSIS || draw.mdl.anm.analysis_type == TRUSS3D_ANALYSIS % TRUSS
            tol = draw.size/125;
        else
            tol = draw.size/400;
        end
        
        % Check if input coordinates are inside an element
        inWhichElems = auxModelFctn('isPointInElem',{[x y z],tol});

        % Initialize support condition
        ebc = [0,0,0,0,0,0];
        
        % Initialize load case
        nodalLoadCase = zeros(12,1);
        
        % Increment number of nodes and create a Node object
        nnp = nnp + 1;
        n = Node(nnp,[x y z],ebc,nodalLoadCase,[],[],[]);
        
        % Insert created Node object in a vector of nodes
        nodes(nnp) = n;
        
        % Update node coordinates uitable
        id = n.id;
        x = n.coord(1);
        y = n.coord(2);
        z = n.coord(3);
        newNode = {id,x,y,z};
        tableData = get(handles.uitable_Nodes,'Data');
        if isempty(tableData)
            set(handles.uitable_Nodes,'Data',newNode)
        elseif strcmp(tableData{1},'')
            set(handles.uitable_Nodes,'Data',newNode)
        else
            set(handles.uitable_Nodes,'Data',vertcat(tableData,newNode))
        end
        clear n
        
        % Update list of nodes to be deleted
        if nnp ~= 0
            nd = zeros(1,nnp);
            for n = 1:nnp
                nd(n) = n;
            end
            nd = num2str(nd,'%d\n');
            set(handles.popupmenu_Delete,'Value',nnp,'string',nd,'Max',nnp)
        else
            set(handles.popupmenu_Delete,'Value',1,'string','No itens',...
                'Max',1)
        end
        
        % Set model object properties
        model = getappdata(0,'model');
        model.nodes = nodes;
        model.nnp = nnp;
        
        % Return draw and model to root
        draw.mdl = model;
        setappdata(0,'draw',draw);
        setappdata(0,'model',model);
        
        % Enable "Process Data" button in main GUI
        mdata = guidata(findobj('Tag','GUI_Main'));
        set(mdata.pushbutton_ProcessData,'Enable','on');
        
        % Disable model type option
        set(mdata.popupmenu_Anm,'Enable','off');
        
        % Disable result buttons
        loadsNeedToBeRedrawn = false; % initialize flag
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
        
        % Update flags
        setappdata(0,'resultType',0);
        setappdata(0,'vis',1);
        
        % Update information panel in GUI_Main
        infoPanelData = get(mdata.uitable_infoPanel,'Data');
        infoPanelData(3,:) = {'Nodes',nnp};
        
        anm = get(mdata.popupmenu_Anm,'Value');
        ndof = cell2mat(infoPanelData(5,2));
        nfreedof = cell2mat(infoPanelData(6,2));
        if anm == 1
            ndof = ndof + 2;
            nfreedof = nfreedof + 2;
        elseif anm == 2 || anm == 3 || anm == 4
            ndof = ndof + 3;
            nfreedof = nfreedof + 3;
        else
            ndof = ndof + 6;
            nfreedof = nfreedof + 6;
        end
        infoPanelData(5,:) = {'DOFs',ndof};
        infoPanelData(6,:) = {'Free DOFs',nfreedof};
        set(mdata.uitable_infoPanel,'Data',infoPanelData)
        
        % Update mouse property
        mouse = getappdata(0,'mouse');
        if ~isempty(mouse.originalData)
            mouse.originalData = infoPanelData;
        end
        setappdata(0,'mouse',mouse)
        
        % Return variables to root
        setappdata(0,'nodes',nodes);
        setappdata(0,'nnp',nnp);
        
        % Check if node is in an intersection
        intersections = getappdata(0,'intersections');
        if ~isempty(intersections)
            for nis = 1:size(intersections,2)
                if norm(intersections(nis).coord - [x y z]) <= 1e-10
                    intersections(nis) = [];
                    break
                end
            end
            if isempty(intersections)
                set(mdata.pushbutton_SolveIntSects,'enable','off')
            end
            setappdata(0,'intersections',intersections)
        end
        
        % Check if node is in an element, if so, divide it in two.
        if ~isempty(inWhichElems)
            elems = getappdata(0,'elems');
            elemsNeedToBeRedrawn = true;
            for i = 1:size(inWhichElems,2)
               % Get elem object to be deleted
               e = inWhichElems(i);
               
               % Check if deleted element had distributed or thermal loads
               if ~isempty(elems(e).load.uniformGbl) || ~isempty(elems(e).load.uniformLcl) ||...
                  ~isempty(elems(e).load.linearGbl) || ~isempty(elems(e).load.linearLcl) ||...
                  elems(e).load.tempVar_X ~= 0 || elems(e).load.tempVar_Y ~= 0 ||...
                  elems(e).load.tempVar_Z ~= 0
                   loadsNeedToBeRedrawn = true;
               end
               
               [~] = auxModelFctn('divideElement',{e,nnp});
            end
        else
            elemsNeedToBeRedrawn = false;
        end
        
        % Enable/disable solve intersections pushbutton (toolbar)
        if size(getappdata(0,'intersections'),2) >= 1
            set(mdata.pushbutton_SolveIntSects,'enable','on')
        else
            set(mdata.pushbutton_SolveIntSects,'enable','off')
        end
        
        % Make GUI a normal window
        gui = findobj('Tag','GUI_Nodes');
        set(gui,'WindowStyle','normal');
        
        % Draw updated model
        redraw(mdata,'Nodes')
        if elemsNeedToBeRedrawn == true
            redraw(mdata,'Elements')
        end
        if loadsNeedToBeRedrawn == true
            redraw(mdata,'Loads')
        end
        
        % Make GUI a modal window
        set(gui,'WindowStyle','modal');
    end
else
    msgbox('There is already a node with these coordinates.', 'Error','error');
end

%--------------------------------------------------------------------------
% Executes on button press in "Delete" pushbutton.
% Deletes a Node object from the list of nodes.
function pushbutton_Delete_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
draw = getappdata(0,'draw');
nodes = getappdata(0,'nodes');
nnp = getappdata(0,'nnp');
elems = getappdata(0,'elems');
nel = getappdata(0,'nel');

if ~strcmp(get(handles.popupmenu_Delete,'String'),'No itens')
    % Get ID of deleted node
    del_n = get(handles.popupmenu_Delete,'Value');
    
    % Get number of elements connected to the deleted node
    n_elem = 0;
    del_elem = zeros(1,model.nel);
    for e = 1:nel
        if (elems(e).nodes(1).id == nodes(del_n).id) || (elems(e).nodes(2).id == nodes(del_n).id)
            n_elem = n_elem + 1;
            del_elem(n_elem) = e;
        end
    end
    del_elem = nonzeros(del_elem)';
    
    if n_elem ~= 0
        choice = questdlg('Deleting this node will automatically delete all elements connected to it.','Delete Node','OK','Cancel','Cancel');
    else
        choice = 'OK';
    end
    
    if strcmp(choice,'OK')
        % Initialize flags
        elemsNeedToBeRedrawn = false;
        nodalLoadsNeedToBeRedrawn = false;
        elemLoadsNeedToBeRedrawn = false;
        
        % Delete all elements connected to the deleted node
        for e = 1:n_elem
            % Get ID of deleted element
            elem_id = del_elem(e) - (e-1);
            
            % Check if deleted element had distributed or thermal loads
            if ~isempty(elems(elem_id).load.uniformGbl) || ~isempty(elems(elem_id).load.uniformLcl) ||...
               ~isempty(elems(elem_id).load.linearGbl) || ~isempty(elems(elem_id).load.linearLcl) ||...
               elems(elem_id).load.tempVar_X ~= 0 || elems(elem_id).load.tempVar_Y ~= 0 ||...
               elems(elem_id).load.tempVar_Z ~= 0
                elemLoadsNeedToBeRedrawn = true;
            end
            
            % Remove deleted element from vector of elements
            elems(elem_id) = [];
            
            % Update number of elements
            nel = nel - 1;
            
            % Set model object properties
            model.elems = elems;
            model.nel = nel;
            
            % Check if element had any unresolved intersections
            intersections = getappdata(0,'intersections');
            if ~isempty(intersections)
                delIntSec = zeros(1,size(intersections,2));
                for nis = 1:size(intersections,2)
                    if ~all(intersections(nis).elems ~= del_elem(e))
                        intersections(nis).elems = nonzeros(intersections(nis).elems - del_elem(e))' + del_elem(e);
                        if size(intersections(nis).elems,2) <= 1
                            delIntSec(nis) = nis;
                        end
                    end
                end
                if ~all(delIntSec == 0)
                    intersections(nonzeros(delIntSec)') = [];
                end
                setappdata(0,'intersections',intersections)
            end
        end
        
        % Update element ids on intersections struct
        intersections = getappdata(0,'intersections');
        if ~isempty(intersections)
            for nis = 1:size(intersections,2)
                for nis_e = 1:size(intersections(nis).elems,2)
                    intersections(nis).elems(nis_e) = intersections(nis).elems(nis_e) - size(find(del_elem < intersections(nis).elems(nis_e)),2);
                end
            end
            setappdata(0,'intersections',intersections)
        end
        
        % Check if elements were deleted
        if n_elem ~= 0
            elemsNeedToBeRedrawn = true;
        end
        
        % Check if deleted node had nodal loads or precribed displacements
        if ~isempty(nodes(del_n).load.static) || ~isempty(nodes(del_n).prescDispl)
            if ~all(nodes(del_n).load.static == 0) || ~all(nodes(del_n).prescDispl == 0)
                nodalLoadsNeedToBeRedrawn = true;
            end
        end
        
        % Check if deleted node had dynamic nodal loads, precribed
        % displacements or concentrated nodal mass
        if ~isempty(nodes(del_n).load.dynamic) || ~isempty(nodes(del_n).prescDispl) || ~isempty(nodes(del_n).displMass)
            if ~all(nodes(del_n).load.dynamic == 0) || ~all(nodes(del_n).prescDispl == 0) || ~all(nodes(del_n).displMass == 0)
                nodalLoadsNeedToBeRedrawn = true;
            end
        end              
        
        
        % Get number of fixed and spring dofs related to the deleted node
        countFixed = 0;
        countSpring = 0;
        for i = 1:6
            if nodes(del_n).ebc(i) == 1
                countFixed = countFixed + 1;
            elseif nodes(del_n).ebc(i) == 2
                countSpring = countSpring + 1;
            end    
        end
        
        % Remove deleted node from vector of nodes
        nodes(del_n) = [];
        
        % Update number of nodes
        nnp = nnp - 1;
        
        % Update list of nodes to be deleted and uitable data
        if nnp ~= 0
            tableData = get(handles.uitable_Nodes,'Data');
            tableData(del_n,:) = [];
            node_id = zeros(1,nnp);
            for n = 1:nnp
                tableData(n,1) = {n};
                nodes(n).id = n;
                node_id(n) = n;
            end
            node_id = num2str(node_id,'%d\n');
            set(handles.popupmenu_Delete,'value',nnp,'string',node_id,...
                'Max',nnp)
        else
            set(handles.popupmenu_Delete,'Value',1,'string','No itens',...
                'Max',1)
        end
        
        % Set updated nodes uitable data
        if nnp ~= 0
            set(handles.uitable_Nodes,'Data',tableData)
        else
            set(handles.uitable_Nodes,'Data',{'','','',''})
        end
        
        % Update information panel in GUI_Main
        infoPanelData = get(mdata.uitable_infoPanel,'Data');
        infoPanelData(3,:) = {'Nodes',nnp};
        infoPanelData(4,:) = {'Elements',nel};
        anm = get(mdata.popupmenu_Anm,'Value');
        ndof = cell2mat(infoPanelData(5,2));
        nfreedof = cell2mat(infoPanelData(6,2));
        nfixeddof = cell2mat(infoPanelData(7,2));
        nspringdof = cell2mat(infoPanelData(8,2));
        if anm == 1
            ndof = ndof - 2;
            nfreedof = nfreedof - 2 + countFixed + countSpring;
        elseif anm == 2 || anm == 3 || anm == 4
            ndof = ndof - 3;
            nfreedof = nfreedof - 3 + countFixed + countSpring;
        else
            ndof = ndof - 6;
            nfreedof = nfreedof - 6 + countFixed + countSpring;
        end
        nfixeddof = nfixeddof - countFixed;
        nspringdof = nspringdof - countSpring;
        infoPanelData(5,:) = {'DOFs',ndof};
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
        model.nnp = nnp;
        
        % Enable/disable solve intersections pushbutton (toolbar)
        if size(getappdata(0,'intersections'),2) >= 1
            set(mdata.pushbutton_SolveIntSects,'enable','on')
        else
            set(mdata.pushbutton_SolveIntSects,'enable','off')
        end
        
        % Enable "Process Data" button in main GUI
        set(mdata.pushbutton_ProcessData,'Enable','on');
        
        % Enable model type option
        if nnp == 0
            set(mdata.popupmenu_Anm,'Enable','on');
            setappdata(0,'vis',0);
        end
        
        % Disable result buttons
        if get(mdata.popupmenu_Results,'value') ~= 1
            nodalLoadsNeedToBeRedrawn = true;
            elemLoadsNeedToBeRedrawn = true;
        end
        set(mdata.popupmenu_Results,'Enable','off','value',1)
        set(mdata.pushbutton_Textual,'Enable','off');
        set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
        set(mdata.text_Element,'string','Elements');
        set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(mdata.edit_Scale,'enable','off','visible','off');
        set(mdata.pushbutton_DynamicResults,'enable','off');
        set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
        
        % Return variables to root
        setappdata(0,'resultType',0);
        setappdata(0,'elems',elems);
        setappdata(0,'nodes',nodes);
        setappdata(0,'nel',nel);
        setappdata(0,'nnp',nnp);
        setappdata(0,'model',model);
        draw.mdl = model;
        setappdata(0,'draw',draw);
        
        % Make GUI a normal window
        gui = findobj('Tag','GUI_Nodes');
        set(gui,'WindowStyle','normal');
        
        % Draw updated model
        redraw(mdata,'Nodes')
        if elemsNeedToBeRedrawn == true
            redraw(mdata,'Elements')
        end
        if nodalLoadsNeedToBeRedrawn == true && elemLoadsNeedToBeRedrawn == false
            redraw(mdata,'Nodal Loads')
        elseif nodalLoadsNeedToBeRedrawn == false && elemLoadsNeedToBeRedrawn == true
            redraw(mdata,'Element Loads')
        elseif nodalLoadsNeedToBeRedrawn == true && elemLoadsNeedToBeRedrawn == true
            redraw(mdata,'Loads')
        end
        
        % Make GUI a modal window
        set(gui,'WindowStyle','modal');
    end
end

%--------------------------------------------------------------------------
% Executes on button press in "Cancel" pushbutton.
% Returns to main GUI without creating any Node object.
function pushbutton_Cancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
delete(gcf)

%--------------------------------------------------------------------------
% Executes during edit_X creation, after setting all properties.
function edit_X_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_X_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Y creation, after setting all properties.
function edit_Y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Y_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Z creation, after setting all properties.
function edit_Z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 2) || (anm == 3)
    set(hObject,'Enable','off');
end

function edit_Z_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Delete creation, after setting all properties.
function popupmenu_Delete_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Delete.
function popupmenu_Delete_Callback(hObject, eventdata, handles)
