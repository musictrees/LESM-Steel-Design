%% Element Dialog Callback Functions
% This file contains the callback functions associated with the "Elements"
% dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_Elements(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Elements_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Elements_OutputFcn, ...
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
% Executes just before Elements GUI is made visible.
% Sets GUI initial properties.
function GUI_Elements_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% Choose default command line output for GUI_Elements
handles.output = hObject;

% Move GUI to the center of the screen
movegui(hObject,'center');

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Set table data
nel = getappdata(0,'nel');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if nel > 0
    elems = getappdata(0,'elems');
    elem = cell(nel,18);
    for e = 1:nel
        id = e;
        etype = elems(e).type;
        if anm == 1 || anm == 4
            type = 'Truss';
        else
            if etype == 1
                type = 'Timoshenko';
            elseif etype == 0
                type = 'Navier';
            end
        end
        
        if (elems(e).mass_consideration)
            mc = 'Yes';
        else
            mc = 'No';
        end
        
        mat = elems(e).material.id;
        sec = elems(e).section.id;
        ni = elems(e).nodes(1).id;
        nf = elems(e).nodes(2).id;
        
        srj_i = zeros(1,3);
        hi = elems(e).hingei;
        if hi == 1
            hi= 'Continuous';
        elseif hi == 0
            hi= 'Hinged';
        elseif hi == 2
            hi= 'Semi-Rigid';
            %Get the rotational stiffness at initial node
            srj_i = [elems(e).kri(1), elems(e).kri(2), elems(e).kri(3)];
        end
        
        srj_f = zeros(1,3);
        hf = elems(e).hingef;
        if hf == 1
            hf = 'Continuous';
        elseif hf == 0
            hf= 'Hinged';
        elseif hf == 2
            hf= 'Semi-Rigid';
            %Get the rotational stiffness at final node
            srj_f = [elems(e).krf(1), elems(e).krf(2), elems(e).krf(3)];
        end
        
        aux_str = '         -';
        if anm == 1  % TRUSS 2D
            k1 = {aux_str,aux_str,aux_str};
            k2 = {aux_str,aux_str,aux_str};
        elseif anm == 2  % FRAME 2D
            k1 = {aux_str,aux_str,srj_i(3)};
            k2 = {aux_str,aux_str,srj_f(3)};
        elseif anm == 3  % GRILLAGE
            k1 = {srj_i(1),srj_i(2),aux_str};
            k2 = {srj_f(1),srj_f(2),aux_str};
        elseif anm == 4  % TRUSS 3D
            k1 = {aux_str,aux_str,aux_str};
            k2 = {aux_str,aux_str,aux_str};
        else  % FRAME 3D
            k1 = {srj_i(1),srj_i(2),srj_i(3)};
            k2 = {srj_f(1),srj_f(2),srj_f(3)};
        end

        v = [elems(e).vz(1), elems(e).vz(2), elems(e).vz(3)];
        v = v / norm(v);

        elem(e,:) = {id,type,mc,mat,sec,ni,nf,hi,hf,k1{1},...
                     k1{2},k1{3},k2{1},k2{2},k2{3},v(1),v(2),v(3)};
    end
    
    nmat = getappdata(0,'nmat');
    mats = cell(1,nmat);
    for nm = 1:nmat
        mats(nm) = {num2str(nm)};
    end
    
    nsec = getappdata(0,'nsec');
    secs = cell(1,nsec);
    for ns = 1:nsec
        secs(ns) = {num2str(ns)};
    end
    
    if anm == 1 % TRUSS 2D
        cEdit = [false(1,2) true(1,3) false(1,2) true(1,2) false(1,9)];
        cFormat = {[] [] {'Yes' 'No'} mats secs [] [] [] [] [] [] [] [] [] [] [] [] []};
    elseif anm == 2  % FRAME 2D
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,2) false(1,2) true(1,1) false(1,2) true(1,1) false(1,3)];
        cFormat = {[] {'Navier' 'Timoshenko'} {'Yes' 'No'} mats secs [] [] {'Hinged' 'Continuous' 'Semi-Rigid'} {'Hinged' 'Continuous' 'Semi-Rigid'} [] [] [] [] [] [] [] [] []};
    elseif anm == 3  % GRILLAGE
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,4) false(1,1) true(1,2) false(1,4)];
        cFormat = {[] {'Navier' 'Timoshenko'} {'Yes' 'No'} mats secs [] [] {'Hinged' 'Continuous' 'Semi-Rigid'} {'Hinged' 'Continuous' 'Semi-Rigid'} [] [] [] [] [] [] [] [] []};
    elseif anm == 4  % TRUSS 3D
        cEdit = [false(1,2) true(1,3) false(1,10) true(1,3)];
        cFormat = {[] [] {'Yes' 'No'} mats secs [] [] [] [] [] [] []};
    else  % FRAME 3D
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,11)];
        cFormat = {[] {'Navier' 'Timoshenko'} {'Yes' 'No'} mats secs [] [] {'Hinged' 'Continuous' 'Semi-Rigid'} {'Hinged' 'Continuous' 'Semi-Rigid'} [] [] [] [] [] [] [] [] []};
    end
    set(handles.uitable_Elements,'Data',elem,...
        'CellEditCallback',@uitable_Elements_CellEditCallback,...
        'CellSelectionCallback',@uitable_Elements_CellSelectionCallback,...
        'ColumnEditable',cEdit,'ColumnFormat',cFormat);
else
    set(handles.uitable_Elements,'Data',{'','','','','','','','','','','','',...
        '','','','','',''},...
        'CellEditCallback',@uitable_Elements_CellEditCallback,...
        'CellSelectionCallback',@uitable_Elements_CellSelectionCallback);
end

% Create list of materials
nmat = getappdata(0,'nmat');
m = zeros(1,nmat);
for i = 1:nmat
    m(i) = i;
end
m = num2str(m,'%d\n');
set(handles.popupmenu_Material,'string',m)

% Create list of cross-sections
nsec = getappdata(0,'nsec');
s = zeros(1,nsec);
for i = 1:nsec
    s(i) = i;
end
s = num2str(s,'%d\n');
set(handles.popupmenu_Section,'string',s)

% Create list of nodes
nnp = getappdata(0,'nnp');
e = zeros(1,nnp);
for i = 1:nnp
    e(i) = i;
end
e = num2str(e,'%d\n');
set(handles.popupmenu_Node1,'string',e)
set(handles.popupmenu_Node2,'string',e)

% Create list of elements to be deleted
set(handles.popupmenu_Delete,'Value',nel)
if nel ~= 0
    elm = 1:nel;
    elm = num2str(elm,'%d\n');
    set(handles.popupmenu_Delete,'Value',nel,'string',elm)
else
    set(handles.popupmenu_Delete,'Value',1,'string','No itens')
end

% If the analysis model is a 2D truss or a 3D truss, set default option to 
% hinged ends.
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 4)
    set(handles.popupmenu_Hinge1,'Value',2,'Enable','off')
    set(handles.popupmenu_Hinge2,'Value',2,'Enable','off')
else    
    set(handles.popupmenu_Hinge1,'Value',1,'Enable','on')
    set(handles.popupmenu_Hinge2,'Value',1,'Enable','on')
end

% If the analysis model is 2D truss or 3D truss, disable shear and 
% flexural deformation checkbox. If the analysis model is grillage,
% disable axial deformation checkbox.
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 1) || (anm == 4)
    set(handles.checkbox_AxialDeformation,'Value',1,'Enable','off')
    set(handles.checkbox_FlexuralDeformation,'Value',0,'Enable','off')
    set(handles.checkbox_ShearDeformation,'Value',0,'Enable','off')
elseif anm == 3
    set(handles.checkbox_AxialDeformation,'Value',0,'Enable','off')
    set(handles.checkbox_ShearDeformation,'Value',0,'Enable','on')
%    set(handles.checkbox_FlexuralDeformation,'Value',1,'Enable','on')
else  
   set(handles.checkbox_AxialDeformation,'Value',1,'Enable','off')
    set(handles.checkbox_ShearDeformation,'Value',0,'Enable','on')
%    set(handles.checkbox_FlexuralDeformation,'Value',1,'Enable','on')
end

% If the analysis model is a 3D truss or a 3D frame, enable Vz setting option
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if (anm == 4) || (anm == 5)
    set(handles.radiobutton_SetZ,'Enable','on')
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
function varargout = GUI_Elements_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Executes on button press in "Add" pushbutton.
% Adds an Elemt object with input properties to the list of elements
function pushbutton_Add_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% Check if end nodes are not the same of a previously created element
include_constants;
equal = 0;
nel = getappdata(0,'nel');
n1 = get(handles.popupmenu_Node1,'Value');
n2 = get(handles.popupmenu_Node2,'Value');
if nel ~= 0
    elems = getappdata(0,'elems');
    for i = 1:nel
        n1j = elems(i).nodes(1).id;
        n2j = elems(i).nodes(2).id;
        if ((n1 == n1j) && (n2 == n2j)) || ((n1 == n2j) && (n2 == n1j))
            equal = 1;
        end
        break
    end
end

 % Get elem mass_consideration flag
    mass_consideration = get(handles.checkbox_MassConsideration,'Value');
   
if equal == 1
    msgbox('There is already an element with this nodal incidence.', 'Error','error');
elseif n1 == n2
    msgbox('Initial and final nodes are the same.', 'Error','error');
else
    mdata = guidata(findobj('Tag','GUI_Main'));
    model = getappdata(0,'model');
    nodes = getappdata(0,'nodes');
    
    % Get material and cross-section IDs
    matid = get(handles.popupmenu_Material,'Value');
    mat = model.materials(matid);
    secid = get(handles.popupmenu_Section,'Value');
    sec = model.sections(secid);
    
    % Get hinge information
    hingei = get(handles.popupmenu_Hinge1,'Value');
    if hingei == 1
        hi = 1; % Continuous
        kri = [];
    elseif hingei == 2
        hi = 0; % Hinged
        kri = [];
    elseif hingei == 3
        hi = 2; % Semi-rigid joint
        
        k1x = str2double(get(handles.k1x_edit,'String'));
        if isnan(k1x)
            k1x = 0;
        end
        k1y = str2double(get(handles.k1y_edit,'String'));
        if isnan(k1y)
            k1y = 0;
        end
        k1z = str2double(get(handles.k1z_edit,'String'));
        if isnan(k1z)
            k1z = 0;
        end
        kri = [abs(k1x), abs(k1y), abs(k1z)];
    end
    
    hingef = get(handles.popupmenu_Hinge2,'Value');
    if hingef == 1
        hf = 1; % Continuous
        krf = [];
    elseif hingef == 2
        hf = 0; % Hinged
        krf = [];
    elseif hingef == 3
        hf = 2; % Semi-rigid joint
        
        k2x = str2double(get(handles.k2x_edit,'String'));
        if isnan(k2x)
            k2x = 0;
        end
        k2y = str2double(get(handles.k2y_edit,'String'));
        if isnan(k2y)
            k2y = 0;
        end
        k2z = str2double(get(handles.k2z_edit,'String'));
        if isnan(k2z)
            k2z = 0;
        end
        krf = [abs(k2x), abs(k2y), abs(k2z)];
    end
    
    % Set initial and final nodes
    ni = nodes(n1);
    nf = nodes(n2);
    
    % Get Vz coordinates
    vzx = str2double(get(handles.edit_Vzx,'String'));
    vzy = str2double(get(handles.edit_Vzy,'String'));
    vzz = str2double(get(handles.edit_Vzz,'String'));
    vz = [vzx vzy vzz];
    
    % Verify if Vz is in the same direction of local axis X:
    % Get nodal coordinates
    xi = ni.coord(1);
    yi = ni.coord(2);
    zi = ni.coord(3);
    xf = nf.coord(1);
    yf = nf.coord(2);
    zf = nf.coord(3);
    % Calculate element local axis X orientation versor
    x = [xf-xi, yf-yi, zf-zi];
	z = [vzx, vzy, vzz];
    % Compare vectors 'x' and 'z'
    w = cross(x,z);
    if (abs(w(1)) < 1e-10) && (abs(w(2)) < 1e-10) && (abs(w(3)) < 1e-10)
        msgbox('Local axes X and Z are parallels. Please, change vector Z coordinates.', 'Error','error');
        return
    end
    
    % Check if model is 2D
    anm = get(mdata.popupmenu_Anm,'Value');
    elemConnect = [];
    newNodes = [];
    if anm <= 3
        coords = [xi yi;
                  xf yf];
        
        crossNodesOutput = auxModelFctn('getCrossNodePoints',coords);
        crossNodePoints = crossNodesOutput{1};
        collinearElems = crossNodesOutput{2};
        elemConnect = crossNodesOutput{3};
        
        mouse = getappdata(0,'mouse');
        cnvs = mouse.getMouseProperty('Canvas');
        % Get canvas borders
        dfltUnits = get(cnvs,'units');
        set(cnvs,'units','normalized');
        limits = get(cnvs,'Position');
        set(cnvs,'units',dfltUnits);
        axisWidth = limits(3);
        tol = axisWidth/50;
        newNodes = auxModelFctn('getNewNodes',{crossNodePoints,collinearElems,elemConnect,tol});
    else
        tol = 10^(-10);
    end
    
    % Get elem type
    if get(handles.checkbox_ShearDeformation,'Value') == 1 % TIMOSHENKO
        type = 1;
    else % NAVIER
        type = 0;
    end
    
    % Create Elem objects and alocate them in elems vector
    if ~isempty(elemConnect)
        h = ones(1,size(elemConnect,1)*2);
        h(1)   = hi;
        h(end) = hf;
        kr = cell(1,size(elemConnect,1)*2);
        kr{1}   = kri;
        kr{end} = krf;
        for e = 1:size(elemConnect,1)
            elem = Elem(type,model.anm,mat,sec,[nodes(elemConnect(e,1)) nodes(elemConnect(e,2))],h(2*e-1),h(2*e),vz,kr{2*e-1},kr{2*e},mass_consideration);
            nel = nel + 1;
            elems(nel) = elem;
        end
        newElements = size(elemConnect,1);
    else
        elem = Elem(type,model.anm,mat,sec,[ni nf],hi,hf,vz,kri,krf,mass_consideration);
        nel = nel + 1;
        elems(nel) = elem;
        newElements = 1;
    end
    
    % Compute number of new semi-rigid joints
    n_srj = 0;
    if hi == 2
        n_srj = n_srj + 1;
    end
    if hf == 2
        n_srj = n_srj + 1;
    end
    
    % Update list of elements to be deleted
    if nel ~= 0
        elm = 1:nel;
        elm = num2str(elm,'%d\n');
        set(handles.popupmenu_Delete,'Value',nel,'string',elm,'Max',nel)
    else
        set(handles.popupmenu_Delete,'Value',1,'string','No itens','Max',1)
    end
    
    % Update elements uitable
    newElems = cell(newElements,18);
    for i = 1:newElements
        e = nel - (i-1);
        id = e;
        etype = elems(e).type;
        if anm == 1 || anm == 4
            type = 'Truss';
        else
            if etype == 1
                type = 'Timoshenko';
            elseif etype == 0
                type = 'Navier';
            end
        end
        
        if (mass_consideration)
            mc = 'Yes';
        else
            mc='No';
        end
        
        mat = elems(e).material.id;
        sec = elems(e).section.id;
        ni = elems(e).nodes(1).id;
        nf = elems(e).nodes(2).id;
        
        hi = elems(e).hingei;
        if hi == 1
            hi = 'Continuous';
        elseif hi == 0
            hi = 'Hinged';
        elseif hi == 2
            hi = 'Semi-Rigid';
        end
        
        hf = elems(e).hingef;
        if hf == 1
            hf = 'Continuous';
        elseif hf == 0
            hf = 'Hinged';
        elseif hf == 2
            hf = 'Semi-Rigid';
        end
        
        v = [elems(e).vz(1), elems(e).vz(2), elems(e).vz(3)];
        v = v / norm(v);
        
        if isempty(elems(e).kri) 
            srj_i = zeros(1,3);
        else
            srj_i = [elems(e).kri(1), elems(e).kri(2), elems(e).kri(3)];
        end
        
        if isempty(elems(e).krf)
            srj_f = zeros(1,3);
        else
            srj_f = [elems(e).krf(1), elems(e).krf(2), elems(e).krf(3)];
        end
        
        aux_str = '         -';
        if anm == 1  % TRUSS 2D
            k1 = {aux_str,aux_str,aux_str};
            k2 = {aux_str,aux_str,aux_str};
        elseif anm == 2  % FRAME 2D
            k1 = {aux_str,aux_str,srj_i(3)};
            k2 = {aux_str,aux_str,srj_f(3)};
        elseif anm == 3  % GRILLAGE
            k1 = {srj_i(1),srj_i(2),aux_str};
            k2 = {srj_f(1),srj_f(2),aux_str};
        elseif anm == 4  % TRUSS 3D
            k1 = {aux_str,aux_str,aux_str};
            k2 = {aux_str,aux_str,aux_str};
        else  % FRAME 3D
            k1 = {srj_i(1),srj_i(2),srj_i(3)};
            k2 = {srj_f(1),srj_f(2),srj_f(3)};
        end
              
        newElems(end-(i-1),:) = {id,type,mc,mat,sec,ni,nf,hi,hf,k1{1},k1{2},k1{3},k2{1},k2{2},k2{3},v(1),v(2),v(3)};
    end

    if anm == 1 % TRUSS 2D
        cEdit = [false(1,2) true(1,3) false(1,2) true(1,2) false(1,9)];
    elseif anm == 2  % FRAME 2D
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,2) false(1,2) true(1,1) false(1,2) true(1,1) false(1,3)];
    elseif anm == 3  % GRILLAGE
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,4) false(1,1) true(1,2) false(1,4)];
    elseif anm == 4 % TRUSS 3D
        cEdit = [false(1,2) true(1,3) false(1,10) true(1,3)];
    else  % FRAME 3D
        cEdit = [false(1,1) true(1,4) false(1,2) true(1,11)];
    end
    tableData = get(handles.uitable_Elements,'Data');
    if isempty(tableData)
        % Reset table format
        nmat = getappdata(0,'nmat');
        mats = cell(1,nmat);
        for nm = 1:nmat
            mats(nm) = {num2str(nm)};
        end
        
        nsec = getappdata(0,'nsec');
        secs = cell(1,nsec);
        for ns = 1:nsec
            secs(ns) = {num2str(ns)};
        end
        if anm == 1 || anm == 4
            cFormat = {[] [] {'Yes' 'No'} mats secs [] [] [] [] [] [] [] [] [] [] [] [] []};
        else
            cFormat = {[] {'Navier' 'Timoshenko'} {'Yes' 'No'} mats secs [] [] {'Hinged' 'Continuous' 'Semi-Rigid'} {'Hinged' 'Continuous' 'Semi-Rigid'} [] [] [] [] [] [] [] [] []};
        end
        set(handles.uitable_Elements,'Data',newElems,'ColumnEditable',cEdit,'ColumnFormat',cFormat)
    elseif strcmp(tableData{1},'')
        % Reset table format
        nmat = getappdata(0,'nmat');
        mats = cell(1,nmat);
        for nm = 1:nmat
            mats(nm) = {num2str(nm)};
        end
        
        nsec = getappdata(0,'nsec');
        secs = cell(1,nsec);
        for ns = 1:nsec
            secs(ns) = {num2str(ns)};
        end
        if anm == 1 || anm == 4
            cFormat = {[] [] {'Yes' 'No'} mats secs [] [] [] [] [] [] [] [] [] [] [] [] []};
        else
            cFormat = {[] {'Navier' 'Timoshenko'} {'Yes' 'No'} mats secs [] [] {'Hinged' 'Continuous' 'Semi-Rigid'} {'Hinged' 'Continuous' 'Semi-Rigid'} [] [] [] [] [] [] [] [] []};
        end
        set(handles.uitable_Elements,'Data',newElems,'ColumnEditable',cEdit,'ColumnFormat',cFormat)
    else
        set(handles.uitable_Elements,'Data',vertcat(tableData,newElems),'ColumnEditable',cEdit)
    end
    
    % Set model object properties
    model.elems = elems;
    model.nel = nel;
    
    % Enable "Process Data" button in main GUI
    set(mdata.pushbutton_ProcessData,'Enable','on');
    
    % Disable result buttons
    if get(mdata.popupmenu_Results,'value') ~= 1
        loadsNeedToBeRedrawn = true;
    else
        loadsNeedToBeRedrawn = false;
    end
    set(mdata.popupmenu_Results,'Enable','off','value',1);
    set(mdata.pushbutton_Textual,'Enable','off');
    set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
    set(mdata.text_Element,'string','Elements');
    set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(mdata.pushbutton_DynamicResults,'enable','off');
    set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    % Update information panel in GUI_Main
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(4,:) = {'Elements',nel};
    infoPanelData(5,:) = {'DOFs',infoPanelData{5,2}+n_srj*model.anm.nrdof};
    infoPanelData(6,:) = {'Free DOFs',infoPanelData{6,2}+n_srj*model.anm.nrdof};
    infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',infoPanelData{9,2}+n_srj*model.anm.nrdof};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Update mouse property
    mouse = getappdata(0,'mouse');
    if ~isempty(mouse.originalData)
        mouse.originalData = infoPanelData;
    end
    setappdata(0,'mouse',mouse)
    
    % Return variables to root
    setappdata(0,'resultType',0);
    setappdata(0,'elems',elems);
    setappdata(0,'nel',nel);
    setappdata(0,'model',model);
    setappdata(0,'nodes',nodes);
    
    % Check if there are crossing points
    if ~isempty(newNodes)
        for nn = 1:size(newNodes,1)
            % Get crossing point coordinates
            x = newNodes(nn,1);
            y = newNodes(nn,2);
            z = 0;

            % Get id of intersected elements
            whichElems = auxModelFctn('isPointInElem',{[x y z],tol});

            % Update vector of handles to intersections (does not
            % create new nodes)
            if ~isempty(whichElems)
                intersections = getappdata(0,'intersections');
                existingIntSect = 0;
                for nis = 1:size(intersections,2)
                   if norm(intersections(nis).coord - [x y z]) <= tol
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
                end
                setappdata(0,'intersections',intSects);
            end
        end
    end
    
    % Enable/disable solve intersections pushbutton (toolbar)
    if size(getappdata(0,'intersections'),2) >= 1
        set(mdata.pushbutton_SolveIntSects,'enable','on')
    else
        set(mdata.pushbutton_SolveIntSects,'enable','off')
    end
    
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    % Draw updated model
    redraw(mdata,'Elements')
    if loadsNeedToBeRedrawn == true
        redraw(mdata,'Loads')
    end
    
    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
end

%--------------------------------------------------------------------------
% Executes on button press in "Delete" pushbutton.
% Deletes an Elem object from the list of elements.
function pushbutton_Delete_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
draw = getappdata(0,'draw');
elems = getappdata(0,'elems');
nel = getappdata(0,'nel');

if ~strcmp(get(handles.popupmenu_Delete,'String'),'No itens')
    % Get ID of deleted element
    del_elem = get(handles.popupmenu_Delete,'Value');
    
    % Check if deleted element had distributed or thermal loads
    if ~isempty(elems(del_elem).load.uniformGbl) || ~isempty(elems(del_elem).load.uniformLcl) ||...
       ~isempty(elems(del_elem).load.linearGbl) || ~isempty(elems(del_elem).load.linearLcl) ||...
       elems(del_elem).load.tempVar_X ~= 0 || elems(del_elem).load.tempVar_Y ~= 0 ||...
       elems(del_elem).load.tempVar_Z ~= 0
        loadsNeedToBeRedrawn = true;
    else
        loadsNeedToBeRedrawn = false;
    end
    
    % Remove deleted element from vector of elements
    elems(del_elem) = [];
    
    % Update number of elements
    nel = nel - 1;
    
    % Update list of elements to be deleted
    if nel ~= 0
        elm = 1:nel;
        elm = num2str(elm,'%d\n');
        set(handles.popupmenu_Delete,'value',nel,'string',elm,'Max',nel)
    else
        set(handles.popupmenu_Delete,'Value',1,'string','No itens','Max',1)
    end
    
    % Set updated elements uitable data
    tableData = get(handles.uitable_Elements,'Data');
    if nel ~= 0
        tableData(del_elem,:) = [];
        for ne = 1:nel
            tableData(ne,1) = {ne};
        end
        set(handles.uitable_Elements,'Data',tableData)
    else
        set(handles.uitable_Elements,'Data',{'','','','','','','','','','','',...
            '','','','','',''})
    end
    
    % Set model object properties
    model.elems = elems;
    model.nel = nel;
    
    % Check if element had any unresolved intersections
    intersections = getappdata(0,'intersections');
    if ~isempty(intersections)
        delIntSec = zeros(1,size(intersections,2));
        for nis = 1:size(intersections,2)
            if ~all(intersections(nis).elems ~= del_elem)
                intersections(nis).elems = nonzeros(intersections(nis).elems - del_elem)' + del_elem;
                if size(intersections(nis).elems,2) <= 1
                    delIntSec(nis) = nis;
                end
            end
            for nis_e = 1:size(intersections(nis).elems,2)
                if intersections(nis).elems(nis_e) > del_elem
                    intersections(nis).elems(nis_e) = intersections(nis).elems(nis_e) - 1;
                end
            end
        end
        if ~all(delIntSec == 0)
            intersections(nonzeros(delIntSec)') = [];
        end
        setappdata(0,'intersections',intersections)
    end
    
    % Enable/disable solve intersections pushbutton (toolbar)
    if size(getappdata(0,'intersections'),2) >= 1
        set(mdata.pushbutton_SolveIntSects,'enable','on')
    else
        set(mdata.pushbutton_SolveIntSects,'enable','off')
    end
    
    % Enable "Process Data" button in main GUI
    set(mdata.pushbutton_ProcessData,'Enable','on');
    
    % Disable result buttons
    if get(mdata.popupmenu_Results,'value') ~= 1
        allLoadsNeedToBeRedrawn = true;
    else
        allLoadsNeedToBeRedrawn = false;
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
    
    % Update information panel in GUI_Main
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(4,:) = {'Elements',nel};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Update mouse property
    mouse = getappdata(0,'mouse');
    if ~isempty(mouse.originalData)
        mouse.originalData = infoPanelData;
    end
    setappdata(0,'mouse',mouse)
    
    % Return variables to root
    setappdata(0,'resultType',0);
    setappdata(0,'elems',elems);
    setappdata(0,'nel',nel);
    setappdata(0,'model',model);
    draw.mdl = model;
    setappdata(0,'draw',draw);
    
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    % Draw updated model
    redraw(mdata,'Elements')
    if loadsNeedToBeRedrawn == true && allLoadsNeedToBeRedrawn == false
        redraw(mdata,'Element Loads')
    elseif allLoadsNeedToBeRedrawn == true
        redraw(mdata,'Loads')
    end

    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
end

%--------------------------------------------------------------------------
% Executes on button press in "Cancel" pushbutton.
% Return to main GUI without creating any Element object.
function pushbutton_Cancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
delete(gcf)

%--------------------------------------------------------------------------
% Executes during popupmenu_Node1 creation, after setting all properties.
function popupmenu_Node1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Node1.
function popupmenu_Node1_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Node2 creation, after setting all properties.
function popupmenu_Node2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Node2.
function popupmenu_Node2_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Hinge1 creation, after setting all properties.
function popupmenu_Hinge1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Hinge1.
function popupmenu_Hinge1_Callback(hObject, eventdata, handles)
hinge = get(hObject,'Value');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if hinge == 3
    if anm == 2
        set(handles.k1x_edit,'Enable','Off','String','0.0');
        set(handles.k1y_edit,'Enable','Off','String','0.0');
        set(handles.k1z_edit,'Enable','On');
    elseif anm == 3
        set(handles.k1x_edit,'Enable','On');
        set(handles.k1y_edit,'Enable','On');
        set(handles.k1z_edit,'Enable','Off','String','0.0');
    elseif anm == 5
        set(handles.k1x_edit,'Enable','On');
        set(handles.k1y_edit,'Enable','On');
        set(handles.k1z_edit,'Enable','On');
    end
else
    set(handles.k1x_edit,'Enable','Off','String','0.0');
    set(handles.k1y_edit,'Enable','Off','String','0.0');
    set(handles.k1z_edit,'Enable','Off','String','0.0');
end
%--------------------------------------------------------------------------
% Executes during popupmenu_Hinge2 creation, after setting all properties.
function popupmenu_Hinge2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Hinge2.
function popupmenu_Hinge2_Callback(hObject, eventdata, handles)
hinge = get(hObject,'Value');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value');
if hinge == 3
    if anm == 2
        set(handles.k2x_edit,'Enable','Off','String','0.0');
        set(handles.k2y_edit,'Enable','Off','String','0.0');
        set(handles.k2z_edit,'Enable','On');
    elseif anm == 3
        set(handles.k2x_edit,'Enable','On');
        set(handles.k2y_edit,'Enable','On');
        set(handles.k2z_edit,'Enable','Off','String','0.0');
    elseif anm == 5
        set(handles.k2x_edit,'Enable','On');
        set(handles.k2y_edit,'Enable','On');
        set(handles.k2z_edit,'Enable','On');
    end
else
    set(handles.k2x_edit,'Enable','Off','String','0.0');
    set(handles.k2y_edit,'Enable','Off','String','0.0');
    set(handles.k2z_edit,'Enable','Off','String','0.0');
end
%--------------------------------------------------------------------------
% Executes during popupmenu_Material creation, after setting all properties.
function popupmenu_Material_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Material.
function popupmenu_Material_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Section creation, after setting all properties.
function popupmenu_Section_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Section.
function popupmenu_Section_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during popupmenu_Delete creation, after setting all properties.
function popupmenu_Delete_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Executes on selection change in popupmenu_Delete.
function popupmenu_Delete_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Vzx creation, after setting all properties.
function edit_Vzx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vzx_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Vzy creation, after setting all properties.
function edit_Vzy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vzy_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes during edit_Vzz creation, after setting all properties.
function edit_Vzz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vzz_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on button press in radiobutton_DefaultZ.
function radiobutton_DefaultZ_Callback(hObject, eventdata, handles)
set(handles.edit_Vzx,'Enable','off','String','0')
set(handles.edit_Vzy,'Enable','off','String','0')
set(handles.edit_Vzz,'Enable','off','String','1')

%--------------------------------------------------------------------------
% Executes on button press in radiobutton_SetZ.
function radiobutton_SetZ_Callback(hObject, eventdata, handles)
set(handles.edit_Vzx,'Enable','on')
set(handles.edit_Vzy,'Enable','on')
set(handles.edit_Vzz,'Enable','on')

%--------------------------------------------------------------------------
% Executes on button press in checkbox_ShearDeformation.
function checkbox_ShearDeformation_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on button press in checkbox_FlexuralDeformation.
function checkbox_FlexuralDeformation_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes on button press in checkbox_AxialDeformation.
function checkbox_AxialDeformation_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
% Executes when cell is selected in uitable_Elements
function uitable_Elements_CellSelectionCallback(hObject, eventdata, handles)
if isempty(eventdata.Indices)
    return
elseif size(eventdata.Indices,1) > 1
    id = [eventdata.Indices(end,1),1];
else
    id = eventdata.Indices;
end
include_constants;
mdata = guidata(findobj('tag','GUI_Elements'));
model = getappdata(0,'model');
anm = model.anm.analysis_type;

if getappdata(0,'nel') == 0
    set(mdata.pushbutton_ApplyMultiElem,'Enable','off');
else
    if anm == TRUSS2D_ANALYSIS || anm == TRUSS3D_ANALYSIS
        if id(2) == 1 || id(2) == 3 || id(2) == 4 || id(2) == 5
            set(mdata.pushbutton_ApplyMultiElem,'Enable','on');
        else
            set(mdata.pushbutton_ApplyMultiElem,'Enable','off');
        end
    else
        if id(2) <= 5 || id(2) == 8 || id(2) == 9
            set(mdata.pushbutton_ApplyMultiElem,'Enable','on');
        else
            set(mdata.pushbutton_ApplyMultiElem,'Enable','off');
        end
    end
end
set(mdata.uitable_Elements,'UserData',id)

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_Elements
function uitable_Elements_CellEditCallback(hObject, eventdata, handles)
mdata = guidata(findobj('tag','GUI_Elements'));
if getappdata(0,'nel') > 0
    set(mdata.pushbutton_Apply,'Enable','on')
else
    set(mdata.pushbutton_Apply,'Enable','off')
end


%--------------------------------------------------------------------------
% Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(hObject, eventdata, handles)
% Get handles and intialize new hinge flag and elemOrientFlag
nel = getappdata(0,'nel');
model = getappdata(0,'model');
elems = getappdata(0,'elems');
tableData = get(handles.uitable_Elements,'Data');
newHingeFlag = false;
elemOrientFlag = false;

% Get handle to GUI_Main
mdata = guidata(findobj('Tag','GUI_Main'));

% Initialize number of semi-rigid joints added / removed
n_srj = 0;

% Get new info and set it to element objects
for e = 1:nel
   % Type
   newType = char(tableData(e,2));
   switch newType
       case 'Timoshenko'
           elems(e).type = 1;
       case 'Navier'
           elems(e).type = 0;
   end
   
   %Mass Consideration
   newMc = char(tableData(e,3));
   switch newMc
       case 'Yes'
           elems(e).mass_consideration = 1;
       case 'No'
           elems(e).mass_consideration = 0;
   end
   
   
   % Material
   newMat = cell2mat(tableData(e,4));
   elems(e).material = model.materials(newMat);
   tableData(e,4) = num2cell(newMat);
   
   % Section
   newSec = cell2mat(tableData(e,5));
   elems(e).section = model.sections(newSec);
   tableData(e,5) = num2cell(newSec);
   
   % Hinge_i
   newHingei = char(tableData(e,8));
   hingei = elems(e).hingei;
   switch newHingei
       case 'Hinged'
           elems(e).hingei = 0;
       case 'Continuous'
           elems(e).hingei = 1;
       case 'Semi-Rigid'
           elems(e).hingei = 2;
   end
   
   % Hinge_f
   newHingef = char(tableData(e,9));
   hingef = elems(e).hingef;
   switch newHingef
       case 'Hinged'
           elems(e).hingef = 0;
       case 'Continuous'
           elems(e).hingef = 1;
       case 'Semi-Rigid'
           elems(e).hingef = 2;
   end
   
   % Check if hinges were added or removed
   if elems(e).hingei ~= hingei || elems(e).hingef ~= hingef
       newHingeFlag = true;
   end
   
   % Check if semi-rigid joints were added / removed
   if elems(e).hingei == 2 && hingei ~= 2
       n_srj = n_srj + 1;
   elseif elems(e).hingei ~= 2 && hingei == 2
       n_srj = n_srj - 1;
   end
   if elems(e).hingef == 2 && hingef ~= 2
       n_srj = n_srj + 1;
   elseif elems(e).hingef ~= 2 && hingef == 2
       n_srj = n_srj - 1;
   end
   
   if elems(e).hingei == 2
       % K1_x
       newk1_x = cell2mat(tableData(e,10));
       if size(newk1_x,2) > 1
           if isempty(elems(e).kri)
               newk1_x = 0;
           else
               newk1_x = elems(e).kri(1);
           end
       elseif isnan(newk1_x)
           if isempty(elems(e).kri)
               newk1_x = 0;
           else
               newk1_x = elems(e).kri(1);
           end
       end
       newk1_x = abs(newk1_x);
       
       % K1_y
       newk1_y = cell2mat(tableData(e,11));
       if size(newk1_y,2) > 1
           if isempty(elems(e).kri)
               newk1_y = 0;
           else
               newk1_y = elems(e).kri(2);
           end
       elseif isnan(newk1_y)
           if isempty(elems(e).kri)
               newk1_y = 0;
           else
               newk1_y = elems(e).kri(2);
           end
       end
       newk1_y = abs(newk1_y);
       
       % K1_z
       newk1_z = cell2mat(tableData(e,12));
       if size(newk1_z,2) > 1
           if isempty(elems(e).kri)
               newk1_z = 0;
           else
               newk1_z = elems(e).kri(3);
           end
       elseif isnan(newk1_z)
           if isempty(elems(e).kri)
               newk1_z = 0;
           else
               newk1_z = elems(e).kri(3);
           end
       end
       newk1_z = abs(newk1_z);
       
       % Assemble vector
       newK1 = [newk1_x, newk1_y, newk1_z];
   
   else
       newK1 = [];
   end
   
   if elems(e).hingef == 2
       % K2_x
       newk2_x = cell2mat(tableData(e,13));
       if size(newk2_x,2) > 1
           if isempty(elems(e).krf)
               newk2_x = 0;
           else
               newk2_x = elems(e).krf(1);
           end
       elseif isnan(newk2_x)
           if isempty(elems(e).krf)
               newk2_x = 0;
           else
               newk2_x = elems(e).krf(1);
           end
       end
       newk2_x = abs(newk2_x);

       % K2_y
       newk2_y = cell2mat(tableData(e,14));
       if size(newk2_y,2) > 1
           if isempty(elems(e).krf)
               newk2_y = 0;
           else
               newk2_y = elems(e).krf(2);
           end
       elseif isnan(newk2_y)
           if isempty(elems(e).krf)
               newk2_y = 0;
           else
               newk2_y = elems(e).krf(2);
           end
       end
       newk2_y = abs(newk2_y);

       % K2_z
       newk2_z = cell2mat(tableData(e,15));
       if size(newk2_z,2) > 1
           if isempty(elems(e).krf)
               newk2_z = 0;
           else
               newk2_z = elems(e).krf(3);
           end
       elseif isnan(newk2_z)
           if isempty(elems(e).krf)
               newk2_z = 0;
           else
               newk2_z = elems(e).krf(3);
           end
       end
       newk2_z = abs(newk2_z);
       
       % Assemble vector
       newK2 = [newk2_x, newk2_y, newk2_z];
   
   else
       newK2 = [];
   end
   
   % Allocate new semi-rigid joint info
   % Node_i
   if isempty(newK1)
       if ~isempty(elems(e).kri)
           elems(e).kri     =  [];
       end
       newK1 =[0,0,0];
       tableData(e,10)  = {0};
       tableData(e,11)  = {0};
       tableData(e,12)  = {0};
   else
       elems(e).kri     =  newK1;
       tableData(e,10)  =  num2cell(newK1(1));
       tableData(e,11)  =  num2cell(newK1(2));
       tableData(e,12)  =  num2cell(newK1(3));
       
       newHingeFlag = true;
   end
   
   % Node_f
   if isempty(newK2)
       if ~isempty(elems(e).krf)
           elems(e).krf     =  [];
       end
       newK2 =[0,0,0];
       tableData(e,13)  = {0};
       tableData(e,14)  = {0};
       tableData(e,15)  = {0};
   else
       elems(e).krf     =  newK2;
       tableData(e,13)  =  num2cell(newK2(1));
       tableData(e,14)  =  num2cell(newK2(2));
       tableData(e,15)  =  num2cell(newK2(3));
       
       newHingeFlag = true;
   end
   
   anm = get(mdata.popupmenu_Anm,'Value');
   
   aux_str = '         -';
   if anm == 1  % TRUSS 2D
       k1 = {aux_str,aux_str,aux_str};
       k2 = {aux_str,aux_str,aux_str};
   elseif anm == 2  % FRAME 2D
       k1 = {aux_str,aux_str,newK1(3)};
       k2 = {aux_str,aux_str,newK2(3)};
   elseif anm == 3  % GRILLAGE
       k1 = {newK1(1),newK1(2),aux_str};
       k2 = {newK2(1),newK2(2),aux_str};
   elseif anm == 4  % TRUSS 3D
       k1 = {aux_str,aux_str,aux_str};
       k2 = {aux_str,aux_str,aux_str};
   else  % FRAME 3D
       k1 = {newK1(1),newK1(2),newK1(3)};
       k2 = {newK2(1),newK2(2),newK2(3)};
   end
   
   tableData(e,10)   = k1(1);
   tableData(e,11)  = k1(2);
   tableData(e,12)  = k1(3);
   tableData(e,13)  = k2(1);
   tableData(e,14)  = k2(2);
   tableData(e,15)  = k2(3);
   
   % vz_x
   newVz_x = cell2mat(tableData(e,16));
   
   % vz_y
   newVz_y = cell2mat(tableData(e,17));
   
   % vz_z
   newVz_z = cell2mat(tableData(e,18));
   
   newVz = [newVz_x, newVz_y, newVz_z];
    
   % Verify if Vz is in the same direction of local axis X:
   % Get nodal coordinates
   ni = elems(e).nodes(1);
   nf = elems(e).nodes(2);
   xi = ni.coord(1);
   yi = ni.coord(2);
   zi = ni.coord(3);
   xf = nf.coord(1);
   yf = nf.coord(2);
   zf = nf.coord(3);
   % Calculate element local axis X orientation versor
   x = [xf-xi, yf-yi, zf-zi];
   % Compare vectors 'x' and 'newVz'
   w = cross(x,newVz);
   if (abs(w(1)) > 1e-10) || (abs(w(2)) > 1e-10) || (abs(w(3)) > 1e-10)
       mdata = guidata(findobj('tag','GUI_Main'));
       if newVz(1) ~= elems(e).vz(1) || newVz(2) ~= elems(e).vz(2) ||...
          newVz(3) ~= elems(e).vz(3)
          if strcmp(get(mdata.orientationButton,'Checked'),'on')
            elemOrientFlag = 1;
          end
       end
       elems(e).vz = newVz;
   end
   tableData(e,16) = num2cell(elems(e).vz(1));
   tableData(e,17) = num2cell(elems(e).vz(2));
   tableData(e,18) = num2cell(elems(e).vz(3));
end

% Set new info in model and save them in root
model.elems = elems;
setappdata(0,'elems',elems);
setappdata(0,'model',model);

% Reset table data
set(handles.uitable_Elements,'Data',tableData)

% Update information panel in GUI_Main
infoPanelData = get(mdata.uitable_infoPanel,'Data');
infoPanelData(5,:) = {'DOFs',infoPanelData{5,2}+n_srj*model.anm.nrdof};
infoPanelData(6,:) = {'Free DOFs',infoPanelData{6,2}+n_srj*model.anm.nrdof};
infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',infoPanelData{9,2}+n_srj*model.anm.nrdof};
set(mdata.uitable_infoPanel,'Data',infoPanelData)

% Disable button
set(hObject,'Enable','off')

% Disable results and enable process data pushbutton
resType = 0;
if strcmp(get(mdata.popupmenu_Results,'Enable'),'on')
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
    
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    % Redraw
    if resType ~= 1
        redraw(mdata,'Elements');
        redraw(mdata,'Loads');
    end
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
    
    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
end

% If any hinge was added or removed, or element orientaion changed,
% redraw model.
if (newHingeFlag ~= 0 || elemOrientFlag ~= 0) && (resType == 0 || resType == 1)
    
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    redraw(mdata,'Elements');
    
    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
end

%--------------------------------------------------------------------------
% Executes on button press in pushbutton_ApplyMultiElem.
function pushbutton_ApplyMultiElem_Callback(hObject, eventdata, handles)
mdata = guidata(findobj('tag','GUI_Main'));

% Get which cell was selected by user
tableIndex = get(handles.uitable_Elements,'UserData');

% Determine what should be asked of user, based on selected cell
switch tableIndex(2)
    case 1  % ELEMENT ID (TYPE, MATERIAL AND SECTION)
        inputStr = sprintf('Apply properties of element %i (type, mass consideration, material and section) to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 2  % ELEMENT TYPE
        inputStr = sprintf('Apply type of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 3  % ELEMENT MASS 
        inputStr = sprintf('Apply mass consideration of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 4  % ELEMENT MATERIAL
        inputStr = sprintf('Apply material of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 5  % ELEMENT SECTION
        inputStr = sprintf('Apply cross-section of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 8  % ELEMENT JOINT 1
        inputStr = sprintf('Apply initial joint type of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
    case 9  % ELEMENT JOINT 2
        inputStr = sprintf('Apply final joint type of element %i to multiple elements:\nExample: "1;3-5;7"  or  "all"',...
                          tableIndex(1));
        sz = [1 65];
end

% Get string input to determine the elements to have properties modified
e_str = char(inputdlg(inputStr,'Apply change to multiple elements',sz));

% Read string entered by user
if isempty(e_str)
    return
elseif strcmp(e_str,'all')
    e_ID = 1:getappdata(0,'nel');
elseif strcmp(e_str,'ALL')
    e_ID = 1:getappdata(0,'nel');
elseif strcmp(e_str,'All')
    e_ID = 1:getappdata(0,'nel');
else
    [flag,e_ID] = readStr(e_str,getappdata(0,'nel'));
    if ~flag
        msgbox('Invalid input data!', 'Error','error');
        return
    end
end

% Get handles to the model object and elem objects
model = getappdata(0,'model');
elems = getappdata(0,'elems');

% Initialize flag for changed joints
newJointFlag = false;
n_srj = 0;

% Get uitable_Elements data
tableData = get(handles.uitable_Elements,'Data');

% Assign information to elements, based on selected cell
switch tableIndex(2)
    case 1  % ELEMENT ID (TYPE, MATERIAL AND SECTION)
        switch tableData{tableIndex(1),2}
            case 'Timoshenko'
                type = 1;
            case 'Navier'
                type = 0;
            case 'Truss'
                type = 0;
        end
        
        switch tableData{tableIndex(1),3}
            case 'Yes'
                mass_consideration = 1;
            case 'No'
                mass_consideration = 0;        
        end
        
        material = model.materials(tableData{tableIndex(1),4});
        section = model.sections(tableData{tableIndex(1),5});
        
        for e = e_ID
            elems(e).type = type;
            elems(e).mass_consideration = mass_consideration;
            elems(e).material = material;
            elems(e).section = section;
            tableData(e,2) = tableData(tableIndex(1),2);
            tableData(e,3) = tableData(tableIndex(1),3);
            tableData(e,4) = tableData(tableIndex(1),4);
            tableData(e,5) = tableData(tableIndex(1),5);
        end
        
    case 2  % ELEMENT TYPE
        switch tableData{tableIndex(1),2}
            case 'Timoshenko'
                type = 1;
            case 'Navier'
                type = 0;
            case 'Truss'
                type = 0;
        end
        for e = e_ID
            elems(e).type = type;
            tableData(e,2) = tableData(tableIndex(1),2);
        end
        
     case 3  % ELEMENT MASS CONSIDERATION
        switch tableData{tableIndex(1),3}
            case 'Yes'
                mass_consideration  = 1;
            case 'No'
                mass_consideration  = 0;           
        end
        for e = e_ID
            elems(e).mass_consideration  = mass_consideration ;
            tableData(e,3) = tableData(tableIndex(1),3);
        end
        
    case 4  % ELEMENT MATERIAL
        material = model.materials(tableData{tableIndex(1),4});
        for e = e_ID
            elems(e).material = material;
            tableData(e,4) = tableData(tableIndex(1),4);
        end
    case 5  % ELEMENT SECTION
        section = model.sections(tableData{tableIndex(1),5});
        for e = e_ID
            elems(e).section = section;
            tableData(e,5) = tableData(tableIndex(1),5);
        end
    case 8 % JOINT 1
        switch tableData{tableIndex(1),8}
            case 'Hinged'
                joint = 0;
            case 'Continuous'
                joint = 1;
            case 'Semi-Rigid'
                joint = 2;
        end
        for e = e_ID
            old_joint = elems(e).hingei;
            elems(e).hingei = joint;
            if joint ~= old_joint
                newJointFlag = true;
            end
            if joint == 2
                if old_joint ~= 2
                    n_srj = n_srj + 1;
                end
                elems(e).kri = [0,0,0];
                newJointFlag = true;
            else
                if old_joint == 2
                    n_srj = n_srj - 1;
                end
                elems(e).kri = [];
            end
            anm = get(mdata.popupmenu_Anm,'Value');
            aux_str = '         -';
            tableData(e,8) = tableData(tableIndex(1),8);
            if anm == 1  % TRUSS 2D
                stiff = {aux_str,aux_str,aux_str};
            elseif anm == 2  % FRAME 2D
                stiff = {aux_str,aux_str,0};
            elseif anm == 3  % GRILLAGE
                stiff = {0,0,aux_str};
            elseif anm == 4  % TRUSS 3D
                stiff = {aux_str,aux_str,aux_str};
            else  % FRAME 3D
                stiff = {0,0,0};
            end
            tableData(e,10) = stiff(1);
            tableData(e,11) = stiff(2);
            tableData(e,12) = stiff(3);
        end
    case 9 % JOINT 2
        switch tableData{tableIndex(1),9}
            case 'Hinged'
                joint = 0;
            case 'Continuous'
                joint = 1;
            case 'Semi-Rigid'
                joint = 2;
        end
        for e = e_ID
            old_joint = elems(e).hingef;
            elems(e).hingef = joint;
            if joint ~= old_joint
                newJointFlag = true;
            end
            if joint == 2
                if old_joint ~= 2
                    n_srj = n_srj + 1;
                end
                elems(e).krf = [0,0,0];
                newJointFlag = true;
            else
                if old_joint == 2
                    n_srj = n_srj - 1;
                end
                elems(e).krf = [];
            end
            anm = get(mdata.popupmenu_Anm,'Value');
            aux_str = '         -';
            tableData(e,9) = tableData(tableIndex(1),9);
            if anm == 1  % TRUSS 2D
                stiff = {aux_str,aux_str,aux_str};
            elseif anm == 2  % FRAME 2D
                stiff = {aux_str,aux_str,0};
            elseif anm == 3  % GRILLAGE
                stiff = {0,0,aux_str};
            elseif anm == 4  % TRUSS 3D
                stiff = {aux_str,aux_str,aux_str};
            else  % FRAME 3D
                stiff = {0,0,0};
            end
            tableData(e,13) = stiff(1);
            tableData(e,14) = stiff(2);
            tableData(e,15) = stiff(3);
        end
end

% Set new info in model and save them in root
model.elems = elems;
setappdata(0,'elems',elems);
setappdata(0,'model',model);

% Reset graphic components
set(handles.uitable_Elements,'Data',tableData)
set(hObject,'Enable','off')
set(handles.pushbutton_Apply,'Enable','off')

% Update information panel in GUI_Main
if n_srj ~= 0
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(5,:) = {'DOFs',infoPanelData{5,2}+n_srj*model.anm.nrdof};
    infoPanelData(6,:) = {'Free DOFs',infoPanelData{6,2}+n_srj*model.anm.nrdof};
    infoPanelData(9,:) = {'Semi-Rigid Joint DOFs',infoPanelData{9,2}+n_srj*model.anm.nrdof};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
end

% Disable results and enable process data pushbutton
resType = 0;
if strcmp(get(mdata.popupmenu_Results,'Enable'),'on')
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
    
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    % Redraw
    if resType ~= 1
        redraw(mdata,'Elements');
        redraw(mdata,'Loads');
    end
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
    
    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
end

% If any hinge was added or removed, redraw model
if newJointFlag && (resType == 0 || resType == 1)
    % Make GUI a normal window
    gui = findobj('Tag','GUI_Elements');
    set(gui,'WindowStyle','normal');
    
    redraw(mdata,'Elements');
    
    % Make GUI a modal window
    set(gui,'WindowStyle','modal');
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


function k1x_edit_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function k1x_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k1y_edit_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function k1y_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k1z_edit_Callback(hObject, eventdata, handles)

function k1z_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function k2x_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2y_edit_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function k2y_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k2z_edit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function k2z_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox_MassConsideration_Callback(hObject,~,handles)
mass_flag = get(handles.checkbox_MassConsideration,'Value');

if mass_flag == 1
    set(handles.checkbox_MassConsideration,'Value',0);
else
    set(handles.checkbox_MassConsideration,'Value',1);
end

