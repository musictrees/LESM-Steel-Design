%% Auxiliary Mouse Function
%
% This file contains functions that are called by an Emouse object to do
% tasks, such as:
% -> Draw nodes
% -> Draw elements
% -> Snap to nodes
% -> Snap to intersections
% -> Snap to elements
% -> Snap to grid
% -> Write selected node/element info on uitables at info panel (GUI_Main)
%
% Input:
% -> whichFunction: string that indicates which function is being called
% -> this: handle to mouse object
% -> fctnArgIn: function input arguments
%
%% ------------------------------------------------------------------------
% Works as a switch, calls other functions.
function  fctnArgOut = auxMouseFctn(whichFunction,this,fctnArgIn)
    fctnArgOut = [];
    if nargin == 2
        fctnArgIn = [];
    end
    switch whichFunction
        case 'dynResIntsctn'
            fctnArgOut = dynResIntsctn(this,fctnArgIn);
        case 'orthoPosition'
            fctnArgOut = orthoPosition(this,fctnArgIn);
        case 'snapToGridPosition'
            fctnArgOut = snapToGridPosition(this,fctnArgIn);
        case 'snapToGridDraw'
            snapToGridDraw(this,fctnArgIn)
        case 'writeNodeInfoPanel'
            writeNodeInfoPanel(this)
        case 'writeElemInfoPanel'
            writeElemInfoPanel(this)
        case 'spatial2Plane'
            fctnArgOut = spatial2Plane(this,fctnArgIn);
        case 'setElemNodes3D'
            setElemNodes3D(this,fctnArgIn);
        case 'selectNodes3D'
            selectNodes3D(this,fctnArgIn)
        case 'selectElems3D'
            selectElems3D(this,fctnArgIn)
        case 'snapToNodes'
            fctnArgOut = snapToNodes(this,fctnArgIn);
        case 'snapToNodes3D'
            snapToNodes3D(this,fctnArgIn);
        case 'snapToIntSects'
            fctnArgOut = snapToIntSects(this,fctnArgIn);
        case 'snapToIntSects3D'
            snapToIntSects3D(this,fctnArgIn);
        case 'snapToElems'
            fctnArgOut = snapToElems(this,fctnArgIn);
        case 'snapToElems3D'
            snapToElems3D(this,fctnArgIn);
        case 'drawNodes'
            drawNodes(this,fctnArgIn)
        case 'drawElements'
            drawElements(this,fctnArgIn)
        case 'drawElements3D'
            drawElements3D(this,fctnArgIn)
    end
end
%% ------------------------------------------------------------------------
function pt = dynResIntsctn(this,coords)
    % Initialize varaible to bu returned
    pt = [];
    
    % Get all plots on this canvas
    plots = findobj('Tag',sprintf('Dynamic%s',this.resStr));
    
    % Check if there aren't any plots
    if isempty(plots)
        return
    end
    
    % Check if coords are within solution interval
    if coords(1) < this.XLim(1) || coords(1) > this.XLim(2)
        return
    end
    
    % Get plots x coordinates (same for all plots)
    X = plots(1).XData;
    
    if length(X) > 2
        % Compute time step
        step = X(2) - X(1);
        
        % Get remainder of dx/step
        remX = rem(coords(1),step);
        
        % Get closer grid point to current point (x)
        if abs(remX) <= step/2
            dx = coords(1) - remX;
        elseif remX < 0
            dx = coords(1) - remX - step;
        else
            dx = coords(1) - remX + step;
        end

        i_1 = find(X == dx);
        clear X
        
        if isempty(i_1)
            return
        end
    else
        return
    end
    
    plots = findobj('Tag','DynamicDispl');
    plots = vertcat(plots,findobj('Tag','DynamicVeloc'));
    plots = vertcat(plots,findobj('Tag','DynamicAccel'));
    pt = zeros(length(plots),2);
    for i = 1:length(plots)
%         Y = plots(i).YData(i_1);
%         if abs(coords(2) - Y) <= this.numTol
%             pt = [dx, Y];
%             return
%         end
        pt(i,:) = [dx, plots(i).YData(i_1)];
    end
    
    
end

%--------------------------------------------------------------------------
function coords = orthoPosition(this,fctnArgIn)
    if nargin == 1
        fctnArgIn = [];
    end
    
    % Get element initial node coords
    a = this.elemCoords(1,1:2);
    
    % Get cursor postion
    if isempty(fctnArgIn)
        b = this.getMouseProperty('CurrentPosition');
    else
        b = fctnArgIn;
    end
    
    % Get distance between points
    dx = b(1) - a(1);
    dy = b(2) - a(2);
    L = norm([dx dy]);
    
    % Get angle between node_i-to-cursor and x axis
    theta = acos(abs(dx)/L);
    
    % Check which angle is closer to theta (0, pi/4, pi/2)
    if theta <= pi/6
        coords = [a(1)+dx a(2)];
    elseif theta <= pi/3
        L_45 = L * cos(abs(theta-pi/4));
        if dx <= 0 && dy <= 0
            coords = a - (L_45 * [cos(pi/4) sin(pi/4)]);
        elseif dx > 0 && dy <= 0
            coords = a + (L_45 * [cos(pi/4) -sin(pi/4)]);
        elseif dx <= 0 && dy > 0
            coords = a + (L_45 * [-cos(pi/4) sin(pi/4)]);
        else
            coords = a + (L_45 * [cos(pi/4) sin(pi/4)]);
        end
    else % theta <= pi/2
        coords = [a(1) a(2)+dy];
    end
end

%--------------------------------------------------------------------------
% Get snap to grid position
function stgPosition = snapToGridPosition(this,fctnArgIn)
    if nargin == 1
        fctnArgIn = [];
    end
    % Get cursor postion
    if isempty(fctnArgIn)
        currentPosition = this.getMouseProperty('CurrentPosition');
    else
        currentPosition = fctnArgIn;
    end
    
    if isempty(currentPosition)
        stgPosition = [];
        return
    end
    
    xcp = currentPosition(1);
    ycp = currentPosition(2);

    % Get remainder of currentPosition/snapPrecision
    remX = rem(xcp,this.snapPrecision);
    remY = rem(ycp,this.snapPrecision);

    % Get closer grid point to current point (x)
    if abs(remX) <= this.snapPrecision/2
        x = xcp - remX;
    elseif remX < 0
        x = xcp - remX - this.snapPrecision;
    else
        x = xcp - remX + this.snapPrecision;
    end

    % Get closer grid point to current point (y)
    if abs(remY) <= this.snapPrecision/2
        y = ycp - remY;
    elseif remY < 0
        y = ycp - remY - this.snapPrecision;
    else
        y = ycp - remY + this.snapPrecision;
    end
    
    % Assemble vector to be returned by function
    stgPosition = [x y];
end

%--------------------------------------------------------------------------
function snapToGridDraw(this,fctnArgIn)
    % Get function arguments
    coords = fctnArgIn{1};
    snapProp = fctnArgIn{2};
    if size(fctnArgIn,2) == 3
        is3D = fctnArgIn{3};
    else
        is3D = false;
    end
    
    % Get coordinates
    x = coords(1);
    y = coords(2);
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Get analysis model
    anm = snapProp(1);
    
    % Get draw properties
    sz = snapProp(2);
    
    % Flag for model still without nodes
    nnp = getappdata(0,'nnp');
    if nnp ~= 0
        noNodesFlag = false;
    else
        noNodesFlag = true;
    end
    
    % Check if size factor is zero (or there are no nodes) and update it.
    if sz == 0 || noNodesFlag == true
        sz = 5;
        draw = getappdata(0,'draw');
        draw.size = sz;
        this.sizeFlag = sz;
        setappdata(0,'draw',draw);
    end
    
    % Draw dynamic snap to grid symbol
    if anm == 1 % TRUSS_2D
        circ = 0 : pi/50 : 2*pi;
        r = sz/125;
        xcirc = x + r * cos(circ);
        ycirc = y + r * sin(circ);
        lim = axis;
        plot(xcirc, ycirc, 'color', [0.9 0.2 0], 'tag','snapGrid');
        if noNodesFlag == true
            xlim([lim(1),lim(2)])
            ylim([lim(3),lim(4)])
            % Turn grid on/off
            if strcmp(get(mdata.gridButton,'Checked'),'on') == 1
                grid on
            else
                grid off
            end
        end
        hold on
    elseif anm == 2 || (anm == 3 && is3D == false) % FRAME_2D or GRILLAGE 2D
        s = sz/200;
        xsq = [x - s , x + s , x + s , x - s];
        ysq = [y - s , y - s , y + s , y + s];
        lim = axis;
        fill(xsq, ysq, [0.9 0.2 0],'tag','snapGrid');
        if noNodesFlag == true
            xlim([lim(1),lim(2)])
            ylim([lim(3),lim(4)])
            % Turn grid on/off
            if strcmp(get(mdata.gridButton,'Checked'),'on') == 1
                grid on
            else
                grid off
            end
        end
        hold on
    elseif anm == 3 && is3D == true % GRILLAGE 3D
        canvas = this.getMouseProperty('Canvas');
        origXLim = get(canvas,'XLim');
        origYLim = get(canvas,'YLim');
        origZLim = get(canvas,'ZLim');
        draw = getappdata(0,'draw');
        draw.cube(coords(1),coords(2),0,sz/230,[0.9 0.2 0],'snapGrid');
        hold on
        if nnp == 0
            axis equal
        end
        set(canvas,'Clipping','off');
        set(canvas, 'XLim', origXLim);
        set(canvas, 'YLim', origYLim);
        set(canvas, 'ZLim', origZLim);
        if nnp == 0
            xlabel('X');
            ylabel('Y');
            zlabel(' ');
        end
        if ~strcmp(get(mdata.rulerButton,'Checked'),'on')
            mdata.axes_Canvas.XAxis.Visible = 'off';
            mdata.axes_Canvas.YAxis.Visible = 'off';
            mdata.axes_Canvas.ZAxis.Visible = 'off';
        end
        if strcmp(get(mdata.gridButton,'Checked'),'on') == 1
            grid on
        else
            grid off
        end
    end
end

%--------------------------------------------------------------------------
% Write selected node info on uitables at info panel
function writeNodeInfoPanel(this)
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Analysis type
    anl = get(mdata.popupmenu_AnalysisType,'Value');
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Get node info
    id = nodes(this.selectedNode).id;                  % id
    nodeX = nodes(this.selectedNode).coord(1);         % coord(x)
    nodeY = nodes(this.selectedNode).coord(2);         % coord(y)
    nodeZ = nodes(this.selectedNode).coord(3);         % coord(z)
    if anl == 1                                        % nodal loads
        if ~isempty(nodes(this.selectedNode).load.static)
            fx = nodes(this.selectedNode).load.static(1);
            fy = nodes(this.selectedNode).load.static(2);
            fz = nodes(this.selectedNode).load.static(3);
            mx = nodes(this.selectedNode).load.static(4);
            my = nodes(this.selectedNode).load.static(5);
            mz = nodes(this.selectedNode).load.static(6);
        else
            fx = 0;
            fy = 0;
            fz = 0;
            mx = 0;
            my = 0;
            mz = 0;
        end
    elseif anl == 2 % nodal mass
        if ~isempty(nodes(this.selectedNode).displMass)
            m = nodes(this.selectedNode).displMass*1000;
        else
            m = 0;
        end
    end
    restrX = nodes(this.selectedNode).ebc(1);          % restr(x)
    switch restrX
        case 0
            rx = 'No';
        case 1
            rx = 'Yes';
        case 2
            rx = 'Spring';
    end
    restrY = nodes(this.selectedNode).ebc(2);          % restr(y)
    switch restrY
        case 0
            ry = 'No';
        case 1
            ry = 'Yes';
        case 2
            ry = 'Spring';
    end
    restrZ = nodes(this.selectedNode).ebc(3);          % restr(z)
    switch restrZ
        case 0
            rz = 'No';
        case 1
            rz = 'Yes';
        case 2
            rz = 'Spring';
    end
    restrRotX = nodes(this.selectedNode).ebc(4);       % restr(rotX)
    switch restrRotX
        case -1
            rrx = 'No';
        case 0
            rrx = 'No';
        case 1
            rrx = 'Yes';
        case 2
            rrx = 'Spring';
    end
    restrRotY = nodes(this.selectedNode).ebc(5);       % restr(rotY)
    switch restrRotY
        case -1
            rry = 'No'; %#ok<NASGU>
        case 0
            rry = 'No'; %#ok<NASGU>
        case 1
            rry = 'Yes'; %#ok<NASGU>
        case 2
            rry = 'Spring'; %#ok<NASGU>
    end
    restrRotZ = nodes(this.selectedNode).ebc(6);       % restr(rotZ)
    switch restrRotZ
        case -1
            rrz = 'No';
        case 0
            rrz = 'No';
        case 1
            rrz = 'Yes';
        case 2
            rrz = 'Spring';
    end

    % Get info currently on panel and set it to originalData
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    if isempty(this.originalData)
        this.originalData = infoPanelData;
    end
    
    % Clear previous uitable data
    infoPanelData = {};
    
    % Set node info to uitable data (obs.: editPanelData = editable panel)
    anm = get(mdata.popupmenu_Anm,'value');
    
    infoPanelData(1,:) = {'Node ID',id};
    infoPanelData(2,:) = {'Coord X [m]',nodeX};
    infoPanelData(3,:) = {'Coord Y [m]',nodeY};
    
    if anm == 4 || anm == 5
        infoPanelData(4,:) = {'Coord Z [m]',nodeZ};
    end
    if strcmp(get(mdata.popupmenu_Results,'Enable'),'on') && get(mdata.popupmenu_Results,'Value') == 2 && anl == 1
        model = getappdata(0,'model');
        n = this.selectedNode;
        if anm == 1
            dx = 1e3*model.D(model.ID(1,n));
            if abs(dx) <= 10^-9
                dx = 0;
            end
            dy = 1e3*model.D(model.ID(2,n));
            if abs(dy) <= 10^-9
                dy = 0;
            end
            infoPanelData(4,:) = {'Global Displ. X [mm]',dx};
            infoPanelData(5,:) = {'Global Displ. Y [mm]',dy};
        elseif anm == 2
            dx = 1e3*model.D(model.ID(1,n));
            if abs(dx) <= 10^-9
                dx = 0;
            end
            dy = 1e3*model.D(model.ID(2,n));
            if abs(dy) <= 10^-9
                dy = 0;
            end
            rotz = model.D(model.ID(3,n));
            if abs(rotz) <= 10^-14
                rotz = 0;
            end
            infoPanelData(4,:) = {'Global Displ. X [mm]',dx};
            infoPanelData(5,:) = {'Global Displ. Y [mm]',dy};
            infoPanelData(6,:) = {'Rotation Z [rad]',rotz};
        elseif anm == 3
            rotx = model.D(model.ID(1,n));
            if abs(rotx) <= 10^-14
                rotx = 0;
            end
            roty = model.D(model.ID(2,n));
            if abs(roty) <= 10^-14
                roty = 0;
            end
            dz = 1e3*model.D(model.ID(3,n));
            if abs(dz) <= 10^-9
                dz = 0;
            end
            infoPanelData(4,:) = {'Global Displ. Z [mm]',dz};
            infoPanelData(5,:) = {'Rotation X [rad]',rotx};
            infoPanelData(6,:) = {'Rotation Y [rad]',roty};
        elseif anm == 4
            dx = 1e3*model.D(model.ID(1,n));
            if abs(dx) <= 10^-9
                dx = 0;
            end
            dy = 1e3*model.D(model.ID(2,n));
            if abs(dy) <= 10^-9
                dy = 0;
            end
            dz = 1e3*model.D(model.ID(3,n));
            if abs(dz) <= 10^-9
                dz = 0;
            end
            infoPanelData(5,:) = {'Global Displ. X [mm]',dx};
            infoPanelData(6,:) = {'Global Displ. Y [mm]',dy};
            infoPanelData(7,:) = {'Global Displ. Z [mm]',dz};
        else
            dx = 1e3*model.D(model.ID(1,n));
            if abs(dx) <= 10^-9
                dx = 0;
            end
            dy = 1e3*model.D(model.ID(2,n));
            if abs(dy) <= 10^-9
                dy = 0;
            end
            dz = 1e3*model.D(model.ID(3,n));
            if abs(dz) <= 10^-9
                dz = 0;
            end
            rotx = model.D(model.ID(4,n));
            if abs(rotx) <= 10^-14
                rotx = 0;
            end
            roty = model.D(model.ID(5,n));
            if abs(roty) <= 10^-14
                roty = 0;
            end
            rotz = model.D(model.ID(6,n));
            if abs(rotz) <= 10^-14
                rotz = 0;
            end
            infoPanelData(5,:) = {'Global Displ. X [mm]',dx};
            infoPanelData(6,:) = {'Global Displ. Y [mm]',dy};
            infoPanelData(7,:) = {'Global Displ. Z [mm]',dz};
            infoPanelData(8,:) = {'Rotation X [rad]',rotx};
            infoPanelData(9,:) = {'Rotation Y [rad]',roty};
            infoPanelData(10,:) = {'Rotation Z [rad]',rotz};
        end
    end
    
    % Clear previous uitable data
    editPanelData = {};
    
    if anm == 1 || anm == 2 % TRUSS 2D or FRAME 2D
        editPanelData(end+1,:) = {'Restr. Dx (y/n)',rx};
        editPanelData(end+1,:) = {'Restr. Dy (y/n)',ry};
        if anm == 2 % FRAME 2D
            editPanelData(end+1,:) = {'Restr. Rz (y/n)',rrz};
        end
        if anl == 1
            editPanelData(end+1,:) = {'Fx [kN]',fx};
            editPanelData(end+1,:) = {'Fy [kN]',fy};
            if anm == 2 % FRAME 2D
                editPanelData(end+1,:) = {'Mz [kNm]',mz};
            end
        elseif anl ==2
            editPanelData(end+1,:) = {'Mass [kg]',m};
        end
    elseif anm == 3 % GRILLAGE
        editPanelData(end+1,:) = {'Restr. Dz (y/n)',rz};
        editPanelData(end+1,:) = {'Restr. Rx and Ry (y/n)',rrx};
        if anl == 1
            editPanelData(end+1,:) = {'Fz [kN]',fz};
            editPanelData(end+1,:) = {'Mx [kNm]',mx};
            editPanelData(end+1,:) = {'My [kNm]',my};
        elseif anl ==2
            editPanelData(end+1,:) = {'Mass [kg]',m};
        end
    elseif anm == 4 || anm == 5 % TRUSS 3D or FRAME 3D
        editPanelData(end+1,:) = {'Restr. Dx (y/n)',rx};
        editPanelData(end+1,:) = {'Restr. Dy (y/n)',ry};
        editPanelData(end+1,:) = {'Restr. Dz (y/n)',rz};
        if anm == 5 % FRAME 3D
            editPanelData(end+1,:) = {'Restr. Rotation (y/n)',rrx};
        end
        if anl == 1
            editPanelData(end+1,:) = {'Fx [kN]',fx};
            editPanelData(end+1,:) = {'Fy [kN]',fy};
            editPanelData(end+1,:) = {'Fz [kN]',fz};
            if anm == 5 % FRAME 3D
                editPanelData(end+1,:) = {'Mx [kNm]',mx};
                editPanelData(end+1,:) = {'My [kNm]',my};
                editPanelData(end+1,:) = {'Mz [kNm]',mz};
            end
        elseif anl ==2
            editPanelData(end+1,:) = {'Mass [kg]',m};
        end
    end
    
    % Set data to uitable at info panel
    set(mdata.uitable_infoPanel,'Data',infoPanelData)

    % Check if current load case is a combination and set editable info
    model = getappdata(0,'model');
    lc = get(mdata.popupmenu_LoadCase,'Value');
    if lc > model.nlc
        set(mdata.uitable_infoPanelEditable,'Enable','on',...
            'ColumnFormat',{[] ,'numeric'},'ColumnEditable',...
            false(1,2),'Data',editPanelData)
    else
        set(mdata.uitable_infoPanelEditable,'Enable','on',...
            'ColumnFormat',{[] ,'numeric'},'ColumnEditable',...
            [false(1,1) true(1,1)],'Data',editPanelData)
    end
end

%--------------------------------------------------------------------------
% Write selected element info on uitables at info panel
function writeElemInfoPanel(this)
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Analysis type
    anl = get(mdata.popupmenu_AnalysisType,'Value');
    
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Get element info
    id     = this.selectedElem;                        % id
    node_i = elems(this.selectedElem).nodes(1).id;     % node_i
    node_f = elems(this.selectedElem).nodes(2).id;     % node_f
    mat    = elems(this.selectedElem).material.id;     % material
    sec    = elems(this.selectedElem).section.id;      % cross-section
    type   = elems(this.selectedElem).type;            % Navier/Timoshenko
    cx     = elems(this.selectedElem).cosine_X;
    cy     = elems(this.selectedElem).cosine_Y;
    switch type
    case 0
        shear = 'No';
    case 1
        shear = 'Yes';  
    end
    hi = elems(this.selectedElem).hingei;              % hinge_i
    switch hi
    case 1
        hinge_i = 'No';
    case 0
        hinge_i = 'Yes';
    case 2
        hinge_i = 'Semi-Rigid';
    end
    hf = elems(this.selectedElem).hingef;              % hinge_f
    switch hf
    case 1
        hinge_f = 'No';
    case 0
        hinge_f = 'Yes';
    case 2
        hinge_f = 'Semi-Rigid';
    end
    if (anl == 1) % static analysis only
        if ~isempty(elems(this.selectedElem).load.uniformGbl)   % unif load
            qx = elems(this.selectedElem).load.uniformGbl(1);
            qy = elems(this.selectedElem).load.uniformGbl(2);
            qz = elems(this.selectedElem).load.uniformGbl(3);
        else
            qx = 0;
            qy = 0;
            qz = 0;
        end
        if ~isempty(elems(this.selectedElem).load.linearGbl)    % linear load
            qx1 = elems(this.selectedElem).load.linearGbl(1);
            qy1 = elems(this.selectedElem).load.linearGbl(2);
            qz1 = elems(this.selectedElem).load.linearGbl(3);
            qx2 = elems(this.selectedElem).load.linearGbl(4);
            qy2 = elems(this.selectedElem).load.linearGbl(5);
            qz2 = elems(this.selectedElem).load.linearGbl(6);
        else
            qx1 = 0;
            qy1 = 0;
            qz1 = 0;
            qx2 = 0;
            qy2 = 0;
            qz2 = 0;
        end
        dtx = elems(this.selectedElem).load.tempVar_X;         % thermal load
        dty = elems(this.selectedElem).load.tempVar_Y;
        dtz = elems(this.selectedElem).load.tempVar_Z;
    end
    
    % Get info currently on panel and set it to originalData
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    if isempty(this.originalData)
        this.originalData = infoPanelData;
    end
    
    % Initialize index flag
    i = 0;
    
    % Clear previous panel data
    infoPanelData = {};
    
    % Set elem info to uitable data (obs.: editPanelData = editable panel)
    anm = get(mdata.popupmenu_Anm,'value');
    
    infoPanelData(1,:) = {'Element ID',id};
    infoPanelData(2,:) = {'Element length [m]',elems(this.selectedElem).length};
    
    if strcmp(get(mdata.popupmenu_Results,'Enable'),'on') && get(mdata.popupmenu_Results,'Value') >= 2 && anl == 1
        infoPanelData(3,:) = {'Local Coord. X [m]',this.elemResults(1,1)};
        if get(mdata.popupmenu_Results,'Value') == 2
            if anm == 1 || anm == 2
                rot = [ cx  cy;
                       -cy  cx ];
                dg = rot' * this.elemResults(1:2,2);
                infoPanelData(4,:) = {'Global Displ. X [mm]',dg(1)*1000};
                infoPanelData(5,:) = {'Global Displ. Y [mm]',dg(2)*1000};
            elseif anm == 3
                infoPanelData(4,:) = {'Global Displ. Z [mm]',this.elemResults(2,2)*1000};
            else
                rot = elems(this.selectedElem).T;
                dg = rot' * this.elemResults(1:3,2);
                infoPanelData(4,:) = {'Global Displ. X [mm]',dg(1)*1000};
                infoPanelData(5,:) = {'Global Displ. Y [mm]',dg(2)*1000};
                infoPanelData(6,:) = {'Global Displ. Z [mm]',dg(3)*1000};
            end
        elseif anm == 2
            switch get(mdata.popupmenu_Results,'Value')
                case 3
                    infoPanelData(4,:) = {'Axial Force [kN]',this.elemResults(1,3)};
                case 4
                    infoPanelData(4,:) = {'Shear Force Y [kN]',this.elemResults(2,3)};
                case 5
                    infoPanelData(4,:) = {'Bending Moment Z [kNm]',this.elemResults(3,3)};
            end
        elseif anm == 1 || anm == 4
            infoPanelData(4,:) = {'Axial Force [kN]',this.elemResults(1,3)};
        elseif anm == 3
            switch get(mdata.popupmenu_Results,'Value')
                case 3
                    infoPanelData(4,:) = {'Torsion Moment [kNm]',- elems(this.selectedElem).torsion_moment(1)};
                case 4
                    infoPanelData(4,:) = {'Shear Force Z [kN]',this.elemResults(1,3)};
                case 5
                    infoPanelData(4,:) = {'Bending Moment Y [kNm]',this.elemResults(2,3)};
            end
        else
            switch get(mdata.popupmenu_Results,'Value')
                case 3
                    infoPanelData(4,:) = {'Axial Force [kN]',this.elemResults(1,3)};
                case 4
                    infoPanelData(4,:) = {'Torsion Moment [kNm]',- elems(this.selectedElem).torsion_moment(1)};
                case 5
                    infoPanelData(4,:) = {'Shear Force Y [kN]',this.elemResults(2,3)};
                case 6
                    infoPanelData(4,:) = {'Shear Force Z [kN]',this.elemResults(3,3)};
                case 7
                    infoPanelData(4,:) = {'Bending Moment Y [kNm]',this.elemResults(5,3)};
                case 8
                    infoPanelData(4,:) = {'Bending Moment Z [kNm]',this.elemResults(4,3)};
            end
        end
    else
        infoPanelData(3,:) = {'Initial node',node_i};
        infoPanelData(4,:) = {'Final node',node_f};
    end
    if anm == 1 || anm == 2 % TRUSS 2D or FRAME 2D
        if anm == 2 % FRAME 2D
            editPanelData(1,:) = {'Shear deformation (y/n)',shear};
            editPanelData(2,:) = {'Hinge 1 (y/n)',hinge_i};
            editPanelData(3,:) = {'Hinge 2 (y/n)',hinge_f};
            i = 3;
        end
        editPanelData(1+i,:) = {'Material',mat};
        editPanelData(2+i,:) = {'Cross-Section',sec};
        if anl == 1
            editPanelData(3+i,:) = {'qx1 [kN/m] (Global)',qx1+qx};
            editPanelData(4+i,:) = {'qy1 [kN/m] (Global)',qy1+qy};
            editPanelData(5+i,:) = {'qx2 [kN/m] (Global)',qx2+qx};
            editPanelData(6+i,:) = {'qy2 [kN/m] (Global)',qy2+qy};
            editPanelData(7+i,:) = {'Temp. X [°C]',dtx};
            editPanelData(8+i,:) = {'Temp. Y [°C]',dty};
        end
    elseif anm  == 3 % GRILLAGE
        editPanelData(1,:) = {'Shear deformation (y/n)',shear};
        editPanelData(2,:) = {'Hinge 1 (y/n)',hinge_i};
        editPanelData(3,:) = {'Hinge 2 (y/n)',hinge_f};
        editPanelData(4,:) = {'Material',mat};
        editPanelData(5,:) = {'Cross-Section',sec};
        if anl == 1
            editPanelData(6,:) = {'qz1 [kN/m] (Global)',qz1+qz};
            editPanelData(7,:) = {'qz2 [kN/m] (Global)',qz2+qz};
            editPanelData(8,:) = {'Temp. Z [°C]',dtz};
        end
    elseif anm  == 4 || anm == 5 % TRUSS 3D or FRAME 3D
        if anm == 5 % FRAME 3D
            editPanelData(1,:) = {'Shear deformation (y/n)',shear};
            editPanelData(2,:) = {'Hinge 1 (y/n)',hinge_i};
            editPanelData(3,:) = {'Hinge 2 (y/n)',hinge_f};
            i = 3;
        end
        editPanelData(1+i,:) = {'Material',mat};
        editPanelData(2+i,:) = {'Cross-Section',sec};
        if anl == 1
            editPanelData(3+i,:) = {'qx1 [kN/m] (Global)',qx1+qx};
            editPanelData(4+i,:) = {'qy1 [kN/m] (Global)',qy1+qy};
            editPanelData(5+i,:) = {'qz1 [kN/m] (Global)',qz1+qz};
            editPanelData(6+i,:) = {'qx2 [kN/m] (Global)',qx2+qx};
            editPanelData(7+i,:) = {'qy2 [kN/m] (Global)',qy2+qy};
            editPanelData(8+i,:) = {'qz2 [kN/m] (Global)',qz2+qz};
            editPanelData(9+i,:) = {'Temp. X [°C]',dtx};
            editPanelData(10+i,:) = {'Temp. Y [°C]',dty};
            editPanelData(11+i,:) = {'Temp. Z [°C]',dtz};
        end
    end
    
    % Set data to uitable at info panel
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    % Check if current load case is a combination and set editable info
    model = getappdata(0,'model');
    lc = get(mdata.popupmenu_LoadCase,'Value');
    if lc > model.nlc
    set(mdata.uitable_infoPanelEditable,'Enable','on',...
        'ColumnFormat',{[] ,'numeric'},'ColumnEditable',...
        false(1,2),'Data',editPanelData)
    else
    set(mdata.uitable_infoPanelEditable,'Enable','on',...
        'ColumnFormat',{[] ,'numeric'},'ColumnEditable',...
        [false(1,1) true(1,1)],'Data',editPanelData)
    end
end

%--------------------------------------------------------------------------
% Get intersection of spatial line on a fixed plane
function planeCoords = spatial2Plane(this,fctnArgIn)
    % Get input data
    clickLine = getMouseProperty(this,'CurrentPosition')';  % 3D click line
    fixedAxis = fctnArgIn{1};        % string for which axis is fixed
    height = fctnArgIn{2};           % coord on fixed axis
    
    switch fixedAxis
        case 'x'
            ax = 1;
        case 'y'
            ax = 2;
        case 'z'
            ax = 3;
    end
    
    % Initialize array to be returned
    planeCoords = [];
    
    % Check if 3D click line is parallel to fixed plane
    if clickLine(1,ax) == clickLine(2,ax)
        return
    end
    
    % Get parametric value of intersection between line and plane
    t = (clickLine(1,ax) - height) / (clickLine(1,ax) - clickLine(2,ax));
    
    % Check if click line doesn't intersect fixed plane
    if t < 0 || t > 1
        return
    end
    
    % Get intersection point coordinates
    planeCoords = zeros(1,3);
    index = 1:3;
    planeCoords(ax) = height;
    index(ax) = [];
    planeCoords(index) = (1-t)*clickLine(1,index) + t*clickLine(2,index);
end

%--------------------------------------------------------------------------
% Set element end coordinates on 3D models
function setElemNodes3D(this,tol)
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));

    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    if get(mdata.popupmenu_Anm,'value') == 3 % GRILLAGE
        [coordsXY,nodeWasDrawn] = getElemNodes(this,tol);

        % If there is an intersection with plane XY, set its coords to element end
        if ~isempty(coordsXY)
            this.elemCoords(this.elemNode,:) = [coordsXY, 0];
            if nodeWasDrawn == true
                this.elemNodeID(this.elemNode) = - this.selectedNode; % minus signal as a flag for node being drawn, not selected
            else
                this.elemNodeID(this.elemNode) = this.selectedNode;
            end
            draw = getappdata(0,'draw');
            coords = this.elemCoords(this.elemNode,:);
            draw.cube(coords(1), coords(2), coords(3), this.sizeFlag/75, [1 0 0],'selectedNode');

        else % if there is no intersection with plane XY, reset elemNode property
            this.elemNode = this.elemNode - 1;
        end
    
    % If a node is selected, set its coords to element end
    elseif this.selectedNode ~= 0
        this.elemNodeID(this.elemNode) = this.selectedNode;
        this.elemCoords(this.elemNode,:) = nodes(this.selectedNode).coord;
        
    elseif this.whichIntSnap ~= 0 % If an intersection is selected, set its coords to element end
        intSects = getappdata(0,'intersections');
        this.elemCoords(this.elemNode,:) = intSects(this.whichIntSnap).coord;
        
        % Draw new node on intersection
        [~] = auxMouseFctn('drawNodes',this,{this.elemCoords(this.elemNode,:),intSects(this.whichIntSnap).elems});
        this.selectedNode = getappdata(0,'nnp');
        this.elemNodeID(this.elemNode) = - this.selectedNode; % minus signal as a flag for node being drawn, not selected
        draw = getappdata(0,'draw');
        coords = this.elemCoords(this.elemNode,:);
        draw.cube(coords(1), coords(2), coords(3), this.sizeFlag/75, [1 0 0],'selectedNode');
        
    elseif this.whichElemSnap ~= 0 % If an element is selected, set its coords to element end
        this.elemCoords(this.elemNode,:) = this.elemPoint;
        
        % Draw new node on element selection point
        [~] = auxMouseFctn('drawNodes',this,{this.elemPoint,this.whichElemSnap});
        this.selectedNode = getappdata(0,'nnp');
        this.elemNodeID(this.elemNode) = - this.selectedNode; % minus signal as a flag for node being drawn, not selected
        draw = getappdata(0,'draw');
        coords = this.elemCoords(this.elemNode,:);
        draw.cube(coords(1), coords(2), coords(3), this.sizeFlag/75, [1 0 0],'selectedNode');
        
    else % if no entities were clicked on, reset elemNode property
        this.elemNode = this.elemNode - 1;
    end
end

%--------------------------------------------------------------------------
% Select nodes on 3D models
function selectNodes3D(this,fctnArgIn)
    % Get function arguments
    lineCoords = fctnArgIn{1};
    tol = fctnArgIn{2};
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));

    % Get nodes near click point
    aux_selectedNode = 0;
    nodes = getappdata(0,'nodes');
    for n = 1:getappdata(0,'nnp')
        isNodeInClick = auxModelFctn('isPointInLine3D',{nodes(n).coord,lineCoords,tol,true});
        if isNodeInClick == true
            if aux_selectedNode == 0
                this.selectedNode = n;
                aux_selectedNode = n;
            elseif norm(nodes(n).coord - lineCoords(1,:)) <...
                   norm(nodes(aux_selectedNode).coord - lineCoords(1,:))
                this.selectedNode = n;
                aux_selectedNode = n;
            end
        end
    end
    
    % Get model and draw object
    model = getappdata(0,'model');
    draw = getappdata(0,'draw');

    % Update Draw object and draw node
    draw.mdl = model;
    if isempty(this.sizeFlag)
        draw.setSize();
        this.sizeFlag = draw.size;
    elseif this.sizeFlag == 0
        draw.setSize();
        this.sizeFlag = draw.size;
    else
        draw.size = this.sizeFlag;
    end 
    sz = draw.size;
    setappdata(0,'draw',draw)
    
    if aux_selectedNode ~= 0
        coords = nodes(this.selectedNode).coord;
        draw.cube(coords(1), coords(2), coords(3), sz/75, [1 0 0],'selectedNode');

        if strcmp(this.mouseCursor,'on')
            % Write node info on uitables at info panel
            writeNodeInfoPanel(this)
        
            % Enable 'delete entities' button
            set(mdata.pushbutton_DeleteEntities,'enable','on')
        end
    else
        if strcmp(this.mouseCursor,'on')
            % Reset info panel uitable with model data
            if ~isempty(this.originalData)
                set(mdata.uitable_infoPanel,'Data',this.originalData)
                set(mdata.uitable_infoPanelEditable,'Data',{},'enable','off')
                set(mdata.pushbutton_ApplyInfoPanel,'enable','off')
                this.originalData = {};
            end
        end
        
        % Update canvas and reinitialize this.selectedNode
        if this.selectedNode ~= 0
            this.whichNodeSnap = this.selectedNode;
            this.selectedNode = 0;
            this.moveAction();
        end
        
        % Disable 'delete entities' button
        set(mdata.pushbutton_DeleteEntities,'enable','off')
    end
end

%--------------------------------------------------------------------------
% Select elements on 3D models
function selectElems3D(this,fctnArgIn)
    % Get function arguments
    lineCoords = fctnArgIn{1}; % click line
    tol = fctnArgIn{2};
    
    % Get click line length
    line_length = norm(lineCoords(2,:) - lineCoords(1,:));
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));

    % Get elems near click point
    aux_selectedElem = [];
    elems = getappdata(0,'elems');
    for e = 1:getappdata(0,'nel')
        elemCoords = [elems(e).nodes(1).coord;
                      elems(e).nodes(2).coord];

        fctnArgOut = auxModelFctn('areLinesCrossed3D',...
                                  {elemCoords,lineCoords,tol});
                              
        isElemInClick = fctnArgOut{1};
        clickPoint = fctnArgOut{2};
                              
        if isElemInClick == true
            norm_clickPoint = norm((clickPoint(2,:) - lineCoords(1,:))/line_length);
            if isempty(aux_selectedElem)
                this.selectedElem = e;
                aux_selectedElem = clickPoint(1,:);
                prev_norm = norm_clickPoint;
            elseif norm_clickPoint < prev_norm
                this.selectedElem = e;
                aux_selectedElem = clickPoint(1,:);
            end
        end
    end
    
    elemPoint = aux_selectedElem;
    e = this.selectedElem;
    if ~isempty(aux_selectedElem)
        elemCoords = [elems(e).nodes(1).coord;
                      elems(e).nodes(2).coord];
        plot3(elemCoords(:,1), elemCoords(:,2), elemCoords(:,3), 'color', [1 0 0],'tag','selectedElem','LineWidth',1.3);
        
        if strcmp(this.mouseCursor,'on') && strcmp(get(mdata.popupmenu_Results,'Enable'),'on') && get(mdata.popupmenu_AnalysisType,'Value') == 1
            scatter3(elemPoint(1),elemPoint(2),elemPoint(3),15,[1 0 0],'filled', 'tag','selectedElem');
            this.elemResults = auxModelFctn('getElemPointDisplAndStress',{elemPoint,e});
        end

        % Write elem info on uitables at info panel
        writeElemInfoPanel(this)
        
        % Enable 'delete entities' button
        set(mdata.pushbutton_DeleteEntities,'enable','on')
    else
        % Reset info panel uitable with model data
        if ~isempty(this.originalData)
            set(mdata.uitable_infoPanel,'Data',this.originalData)
            set(mdata.uitable_infoPanelEditable,'Data',{},'enable','off')
            set(mdata.pushbutton_ApplyInfoPanel,'enable','off')
            this.originalData = {};
        end
        
        % Update canvas and reinitialize this.selectedNode
        if this.selectedElem ~= 0
            this.whichNodeSnap = this.selectedNode;
            this.selectedNode = 0;
            this.moveAction();
        end
        
        % Disable 'delete entities' button
        set(mdata.pushbutton_DeleteEntities,'enable','off')
    end
end

%--------------------------------------------------------------------------
% Check if there is a node close to cursor and draw dynamic snap symbol.
% Create an invisible square around each node, if mouse is inside one of
% the squares, the related node will be drawn in red.
function nodeFlag = snapToNodes(this,fctnArgIn)
    % Get function arguments
    coords = fctnArgIn{1};
    snapProp = fctnArgIn{2};
    
    % Get coordinates
    x = coords(1);
    y = coords(2);
    
    % Initialize flag
    nodeFlag = 0;  % flag for snapped node (0 = no nodes near, 1 = a node has been snapped)
    
    % Get number of nodes
    nnp = getappdata(0,'nnp');
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Get analysis model
    anm = snapProp(1);
    
    % Get draw properties
    sz = snapProp(2);
    axisWidth = snapProp(3);
    
    
    for n = 1:nnp
        % Check if a node has not already been snapped
        if nodeFlag ~= 1
            
            % Get node coordinates
            xn = nodes(n).coord(1);
            yn = nodes(n).coord(2);

            % Get scale factor
            if anm == 1
                scf = axisWidth/12;
            else
                scf = axisWidth/18;
            end

            % Define ivisible square around node
            snapNodeX = [(xn - scf) (xn - scf)...
                         (xn + scf) (xn + scf)...
                         (xn - scf)];
            snapNodeY = [(yn - scf) (yn + scf)...
                         (yn + scf) (yn - scf)...
                         (yn - scf)];

            % Check if cursor is inside said square
            if inpolygon(x,y,snapNodeX,snapNodeY) == 1

                % Set snap to node properties
                this.nodeSnap = [xn yn];
                previousNode = this.whichNodeSnap;
                this.whichNodeSnap = n;
                nodeFlag = 1;

                % Draw snap to node symbol
                if anm == 1 % TRUSS_2D
                    circ = 0 : pi/50 : 2*pi;
                    r = sz/125;
                    x = xn;
                    y = yn;
                    xcirc = x + r * cos(circ);
                    ycirc = y + r * sin(circ);
                    if this.selectedNode ~= 0
%                         scatter(x, y, 4*sz, [0 0.62 0], 'filled','tag', 'snapNode2')
                        plot(xcirc, ycirc, 'color', [0.9 0.2 0], 'tag','snapNode2');
                    else
%                         scatter(x, y, 4*sz, [0 0.62 0], 'filled','tag', 'snapNode')
                        plot(xcirc, ycirc, 'color', [0.9 0.2 0], 'tag','snapNode');
                    end
                    hold on
                elseif anm == 2 || anm == 3 % FRAME_2D or GRILLAGE
                    s = sz/200;
                    x = xn;
                    y = yn;
                    xsq = [x - s , x + s , x + s , x - s];
                    ysq = [y - s , y - s , y + s , y + s];
                    if this.selectedNode ~= 0
%                         scatter(x, y, 4.5*sz, [0 0.62 0], 'filled','tag', 'snapNode2')
                        fill(xsq, ysq, [0.9 0.2 0],'tag','snapNode2');
                    else
%                         scatter(x, y, 4.5*sz, [0 0.62 0], 'filled','tag', 'snapNode')
                        fill(xsq, ysq, [0.9 0.2 0],'tag','snapNode');
                    end
                    hold on
                end

            elseif n ~= this.selectedNode && n == this.whichNodeSnap  % check if node away from cursor is
                if ~isempty(findobj('tag', 'snapNode2'))              % currently snapped, if so, delete dynamic plot.
                    delete(findobj('tag', 'snapNode2'));
                end
                if this.selectedNode == 0
                    delete(findobj('tag', 'snapNode'));
                end
            end
        elseif n ~= this.selectedNode && n == previousNode   % check if node away from cursor is
            if ~isempty(findobj('tag', 'snapNode2'))         % currently snapped, if so, delete dynamic plot.
                delete(findobj('tag', 'snapNode2'));
            end
            if this.selectedNode == 0
                delete(findobj('tag', 'snapNode'));
            end
        end
    end
end

%--------------------------------------------------------------------------
% Snap to nodes on 3D models
function snapToNodes3D(this,fctnArgIn)
    % Get function arguments
    lineCoords = fctnArgIn{1};
    tol = fctnArgIn{2};

    % Get nodes near cursor current point
    aux_snapNode = 0;
    nodes = getappdata(0,'nodes');
    for n = 1:getappdata(0,'nnp')
        isNodeInCursor = auxModelFctn('isPointInLine3D',...
                                     {nodes(n).coord,lineCoords,tol,true});
        if isNodeInCursor == true
            if aux_snapNode == 0
                this.whichNodeSnap = n;
                aux_snapNode = n;
            elseif norm(nodes(n).coord - lineCoords(1,:)) <...
                   norm(nodes(aux_snapNode).coord - lineCoords(1,:))
                this.whichNodeSnap = n;
                aux_snapNode = n;
            end
        end
    end
    
    % Get draw object
    draw = getappdata(0,'draw');

    % Get size
    if isempty(this.sizeFlag)
        draw.setSize();
        this.sizeFlag = draw.size;
    elseif this.sizeFlag == 0
        draw.setSize();
        this.sizeFlag = draw.size;
    else
        draw.size = this.sizeFlag;
    end
    sz = draw.size;
    
    % Return draw object to root
    setappdata(0,'draw',draw)
    
    if aux_snapNode ~= 0
        coords = nodes(this.whichNodeSnap).coord;
        draw.cube(coords(1), coords(2), coords(3), sz/75, [1 0 0],'snapNode');
    else      
        % Update canvas and reinitialize this.whichNodeSnap
        this.whichNodeSnap = 0;
    end
end

%--------------------------------------------------------------------------
function intSectFlag = snapToIntSects(this,fctnArgIn)
    % Get function arguments
    coords = fctnArgIn{1};
    snapProp = fctnArgIn{2};
    
    % Get coordinates
    x = coords(1);
    y = coords(2);
    
    % Initialize flag
    intSectFlag = 0;  % flag for snapped intersection (0 = no intsect near, 1 = an intsect has been snapped)
    
    % Get vector of handles to elemIntersection objects
    intersections = getappdata(0,'intersections');
    
    % Get number of intsects
    nis = size(intersections,2);
    
    % Get analysis model
    anm = snapProp(1);
    
    % Get draw properties
    sz = snapProp(2);
    axisWidth = snapProp(3);
    
    for n = 1:nis
        % Check if an intsect has not already been snapped
        if intSectFlag ~= 1
            
            % Get intsect coordinates
            xi = intersections(n).coord(1);
            yi = intersections(n).coord(2);

            % Get scale factor
            if anm == 1
                scf = axisWidth/12;
            else
                scf = axisWidth/15;
            end

            % Define ivisible square around intsect
            snapNodeX = [(xi - scf) (xi - scf)...
                         (xi + scf) (xi + scf)...
                         (xi - scf)];
            snapNodeY = [(yi - scf) (yi + scf)...
                         (yi + scf) (yi - scf)...
                         (yi - scf)];

            % Check if cursor is inside said square
            if inpolygon(x,y,snapNodeX,snapNodeY) == 1

                % Set snap to intsect properties
                this.intSnap = [xi yi];
                this.whichIntSnap = n;
                intSectFlag = 1;

                % Draw snap to node symbol
                if anm == 1 % TRUSS_2D
                    circ = 0 : pi/50 : 2*pi;
                    r = sz/125;
                    x = xi;
                    y = yi;
                    xcirc = x + r * cos(circ);
                    ycirc = y + r * sin(circ);
                    plot(xcirc, ycirc, 'color', [0.9 0.2 0], 'tag','snapIntSect');
%                     scatter(x, y, 8*sz, [0 0.62 0], 'filled','tag', 'snapIntSect')
                    hold on
                elseif anm == 2 || anm == 3 % FRAME_2D or GRILLAGE
                    s = sz/200;
                    x = xi;
                    y = yi;
                    xsq = [x - s , x + s , x + s , x - s];
                    ysq = [y - s , y - s , y + s , y + s];
                    fill(xsq, ysq, [0.9 0.2 0],'tag','snapIntSect');
%                     scatter(x, y, 8*sz, [0 0.62 0], 'filled','tag', 'snapIntSect')
                    hold on
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% Snap to nodes on 3D models
function snapToIntSects3D(this,fctnArgIn)
    % Get function arguments
    lineCoords = fctnArgIn{1};
    tol = fctnArgIn{2};

    % Get nodes near cursor current point
    aux_snapIntSect = 0;
    ints = getappdata(0,'intersections');
    for n = 1:size(ints,2)
        isIntSectInCursor = auxModelFctn('isPointInLine3D',...
                                         {ints(n).coord,lineCoords,tol,true});
        if isIntSectInCursor == true
            if aux_snapIntSect == 0
                this.whichIntSnap = n;
                aux_snapIntSect = n;
            elseif norm(ints(n).coord - lineCoords(1,:)) <...
                   norm(ints(aux_snapIntSect).coord - lineCoords(1,:))
                this.whichIntSnap = n;
                aux_snapIntSect = n;
            end
        end
    end
    
    % Get model and draw object
    draw = getappdata(0,'draw');

    % Get size
    sz = this.sizeFlag;
    
    if aux_snapIntSect ~= 0
        coords = ints(this.whichIntSnap).coord;
        draw.cube(coords(1), coords(2), coords(3), sz/75, [1 0 0],'snapIntSect');
    else      
        % Update canvas and reinitialize this.whichNodeSnap
        this.whichIntSnap = 0;
    end
end

%--------------------------------------------------------------------------
% Check if there is an element close to cursor and draw dynamic snap symbol.
% Create an invisible rectangle around each element, if
% mouse is inside a rectangle, the related element will be
% drawn in red.
function inElemFlag = snapToElems(this,fctnArgIn)
    % Get function arguments
    coords = fctnArgIn{1};
    snapProp = fctnArgIn{2};
    
    % Get coordinates
    x = coords(1);
    y = coords(2);
    
    % Initialize flag
    inElemFlag = 0;  % flag for snapped element (0 = no elements near, 1 = an element has been snapped)
    
    % Get number of elements
    nel = getappdata(0,'nel');
    
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Get draw properties
    axisWidth = snapProp(3);
    
    for e = 1:nel
        % Check if an element has not already been snapped
        if inElemFlag ~= 1
            
            % Get element end coordinates
            xe1 = elems(e).nodes(1).coord(1);
            ye1 = elems(e).nodes(1).coord(2);
            xe2 = elems(e).nodes(2).coord(1);
            ye2 = elems(e).nodes(2).coord(2);
            
            % Get element cosines
            cx = elems(e).cosine_X;
            cy = elems(e).cosine_Y;
            
            % Get scale factor (draw)
            scf = axisWidth/22;
            
            % Define ivisible rectangle around element
            snapElemX = [(xe1 - scf*(cx-cy)) (xe1 - scf*(cx+cy))...
                         (xe2 + scf*(cx-cy)) (xe2 + scf*(cx+cy))...
                         (xe1 - scf*(cx-cy))];

            snapElemY = [(ye1 - scf*(cx+cy)) (ye1 + scf*(cx-cy))...
                         (ye2 + scf*(cx+cy)) (ye2 - scf*(cx-cy))...
                         (ye1 - scf*(cx+cy))];
                     
            % Check if cursor is inside said rectangle
            if inpolygon(x,y,snapElemX,snapElemY) == 1
                
                % Set snap to element properties
                this.elemSnap = [xe1 ye1;
                                 xe2 ye2];
                previousElem = this.whichElemSnap;
                this.whichElemSnap = e;
                inElemFlag = 1;
                
                % Draw snap to element symbol
                if this.selectedElem ~= 0
                    plot(this.elemSnap(:,1),this.elemSnap(:,2),'color',[0.9 0.2 0],'linewidth',1.2,'tag','snapElem2');
                else
                    plot(this.elemSnap(:,1),this.elemSnap(:,2),'color',[0.9 0.2 0],'linewidth',1.2,'tag','snapElem');
                end
                hold on
                
                % Get element angle with X axis (plane XY)
                alpha = atan((ye2-ye1)/(xe2-xe1));
                if isnan(alpha)
                    if (ye2 - ye1) >= 0
                        alpha = pi * 0.5;
                    else
                        alpha = - pi * 0.5;
                    end
                end

                % Get node_i-to-point angle with X axis (plane XY)
                alphaP = atan((y-ye1)/(x-xe1));
                if isnan(alphaP)
                    if (y - ye1) >= 0
                        alphaP = pi * 0.5;
                    else
                        alphaP = - pi * 0.5;
                    end
                end

                % Get angle between node_i-to-point and element
                theta = alphaP - alpha;

                % Get node_i-to-point length
                Lp = sqrt((x-xe1)^2 + (y-ye1)^2);

                % Projection of node_i-to-point on element
                Lpe = Lp * abs(cos(theta));

                % Coordinates of specific point inside element
                xpe = xe1 + Lpe * cx;
                ype = ye1 + Lpe * cy;
                this.elemPoint = [xpe ype];
                
                % Check if modelling options are on
                if strcmp(this.drawNode,'on') || strcmp(this.drawElem,'on')
                    scatter(xpe, ype,25, [0.9 0.2 0], 'filled','tag', 'snapElemPoint')
                    hold on
                end
                
            elseif e ~= this.selectedElem && e == this.whichElemSnap   % check if elem away from mouse is
                if ~isempty(findobj('tag', 'snapElem2'))               % currently snapped, if so, delete dynamic plot.
                    delete(findobj('tag', 'snapElem2'));
                end
                if this.selectedElem == 0
                    delete(findobj('tag', 'snapElem'));
                end
            end
        elseif e ~= this.selectedElem && e == previousElem      % check if elem away from mouse is
            if ~isempty(findobj('tag', 'snapElem2'))            % currently snapped, if so, delete dynamic plot.
                delete(findobj('tag', 'snapElem2'));
            end
            if this.selectedElem == 0
                delete(findobj('tag', 'snapElem'));
            end
        end
    end
end

%--------------------------------------------------------------------------
% Snap to elements on 3D models
function snapToElems3D(this,fctnArgIn)
    % Get function arguments
    lineCoords = fctnArgIn{1}; % click line
    tol = fctnArgIn{2};
    
    % Get click line length
    line_length = norm(lineCoords(2,:) - lineCoords(1,:));

    % Get elems near cursor current point
    aux_snapElem = [];
    elems = getappdata(0,'elems');
    for e = 1:getappdata(0,'nel')
        elemCoords = [elems(e).nodes(1).coord;
                      elems(e).nodes(2).coord];
        fctnArgOut = auxModelFctn('areLinesCrossed3D',...
                                  {elemCoords,lineCoords,tol});
                              
        isElemInClick = fctnArgOut{1};
        clickPoint = fctnArgOut{2};
                              
        if isElemInClick == true
            norm_clickPoint = norm((clickPoint(2,:) - lineCoords(1,:))/line_length);
            if isempty(aux_snapElem)
                this.whichElemSnap = e;
                this.elemPoint = clickPoint(1,:);
                aux_snapElem = clickPoint(1,:);
                prev_norm = norm_clickPoint;
            elseif norm_clickPoint < prev_norm
                this.whichElemSnap = e;
                this.elemPoint = clickPoint(1,:);
                aux_snapElem = clickPoint(1,:);
            end
        end
    end
    
    if ~isempty(aux_snapElem)
        scatter3(aux_snapElem(1),aux_snapElem(2),aux_snapElem(3),15,[1 0 0],'filled', 'tag','snapElem');
    else
        this.whichElemSnap = 0;
    end
end

%--------------------------------------------------------------------------
function drawNodes(this,fctnArgIn)
   % Get function arguments
   coords = fctnArgIn{1};
   inWhichElems = fctnArgIn{2};
   
   % Get coordinates
   x = coords(1);
   y = coords(2);
   z = coords(3);
   
   % Check if input coordinates are not the same of a previously created node
    equal = 0;
    nnp = getappdata(0,'nnp');
    if nnp ~= 0
        nodes = getappdata(0,'nodes');
        for i = 1:nnp
            xi = x;
            yi = y;
            zi = z;
            xj = nodes(i).coord(1);
            yj = nodes(i).coord(2);
            zj = nodes(i).coord(3);
            if (xi == xj) && (yi == yj) && (zi == zj)
                equal = 1;
            end
        end
    end
    if equal == 0
        mdata = guidata(findobj('Tag','GUI_Main'));
        % Enable process data pushbutton
        if strcmp(get(mdata.pushbutton_ProcessData,'Enable'),'off') && nnp ~= 0
            set(mdata.pushbutton_ProcessData,'Enable','on')
        end

        % Initialize support condition
        ebc = [0,0,0,0,0,0];

        % Initialize load case
        nodalLoadCase = zeros(12,1);

        % Increment number of nodes and create a Node object
        nnp = nnp + 1;
        n = Node(nnp,[x y z],ebc,nodalLoadCase,[],[],[]);

        % Insert created Node object in a vector of nodes
        nodes(nnp) = n;

        % Set model object properties
        model = getappdata(0,'model');
        model.nodes = nodes;
        model.nnp = nnp;
        
        % Return variables to root
        setappdata(0,'nodes',nodes);
        setappdata(0,'nnp',nnp);

        % Disable model type option
        set(mdata.popupmenu_Anm,'Enable','off');

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
        
        % Update originalData property
        if ~isempty(this.originalData)
            this.originalData(3,:) = {'Nodes',nnp};
            this.originalData(5,:) = {'DOFs',ndof};
            this.originalData(6,:) = {'Free DOFs',nfreedof};
        end

        % Get draw object
        draw = getappdata(0,'draw');

        % Update Draw object and draw node
        draw.mdl = model;
        if isempty(this.sizeFlag)
            draw.setSize();
            this.sizeFlag = draw.size;
        elseif this.sizeFlag == 0
            draw.setSize();
            this.sizeFlag = draw.size;
        else
            draw.size = this.sizeFlag;
        end 
        sz = draw.size;
        canvas = this.getMouseProperty('Canvas');
        origXLim = get(canvas,'XLim');
        origYLim = get(canvas,'YLim');
        origZLim = get(canvas,'ZLim');
        if anm == 1
            draw.circle(nodes(nnp).coord(1), nodes(nnp).coord(2), sz/125, [0 0 0],'drawNodes');
            hold on
            set(canvas, 'XLim', origXLim);
            set(canvas, 'YLim', origYLim);
        elseif anm == 2
            draw.square(nodes(nnp).coord(1), nodes(nnp).coord(2), sz/200, [0 0 0],'drawNodes');
            hold on
            set(canvas, 'XLim', origXLim);
            set(canvas, 'YLim', origYLim);
        elseif anm == 3
            if (strcmp(this.drawNode,'on') || strcmp(this.drawElem,'on')) && strcmp(get(mdata.togglebutton_2DView,'state'),'on')
                draw.square(nodes(nnp).coord(1), nodes(nnp).coord(2), sz/200, [0 0 0],'drawNodes');
                hold on
                set(canvas, 'XLim', origXLim);
                set(canvas, 'YLim', origYLim);
            else
                draw.cube(nodes(nnp).coord(1),nodes(nnp).coord(2),nodes(nnp).coord(3),sz/230,[0 0 0],'drawNodes');
                hold on
                if nnp == 1
                    axis equal
                end
                set(canvas, 'XLim', origXLim);
                set(canvas, 'YLim', origYLim);
                set(canvas, 'ZLim', origZLim);
                if nnp == 1
                    xlabel('X');
                    ylabel('Y');
                    zlabel(' ');
                end
                if ~strcmp(get(mdata.rulerButton,'Checked'),'on')
                    mdata.axes_Canvas.XAxis.Visible = 'off';
                    mdata.axes_Canvas.YAxis.Visible = 'off';
                    mdata.axes_Canvas.ZAxis.Visible = 'off';
                end
            end
        elseif anm == 4
            draw.sphere(nodes(nnp).coord(1),nodes(nnp).coord(2),nodes(nnp).coord(3),sz/125,'drawNodes');
            hold on
            set(canvas, 'XLim', origXLim);
            set(canvas, 'YLim', origYLim);
            set(canvas, 'ZLim', origZLim);
        elseif anm == 5
            draw.cube(nodes(nnp).coord(1),nodes(nnp).coord(2),nodes(nnp).coord(3),sz/200,[0 0 0],'drawNodes');
            hold on
            set(canvas, 'XLim', origXLim);
            set(canvas, 'YLim', origYLim);
            set(canvas, 'ZLim', origZLim);
        end

        if strcmp(get(mdata.nodeIDButton,'Checked'),'on')
            draw.nodeID();
        end
        
        if nnp == 1
            if (strcmp(get(mdata.rulerButton,'Checked'),'on') == 1)
                mdata.axes_Canvas.XAxis.Visible = 'on';
                mdata.axes_Canvas.YAxis.Visible = 'on';
                mdata.axes_Canvas.ZAxis.Visible = 'on';
            else
                mdata.axes_Canvas.XAxis.Visible = 'off';
                mdata.axes_Canvas.YAxis.Visible = 'off';
                mdata.axes_Canvas.ZAxis.Visible = 'off';
            end
        end
        
        if (strcmp(get(mdata.gridButton,'Checked'),'on') == 1)
            grid on
        else
            grid off
        end
        
        if anm == 3
            canvas.Clipping = 'off';
        end
        
        % Check if node is in an intersection
        intersections = getappdata(0,'intersections');
        if ~isempty(intersections)
            for nis = 1:size(intersections,2)
                if norm(intersections(nis).coord - [x y z]) <= 10^-12
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
            for i = 1:size(inWhichElems,2)
                e = inWhichElems(i);
                % Divide element
                [~] = auxModelFctn('divideElement',{e,nnp});
                % Update originalData property
                if ~isempty(this.originalData)
                    this.originalData(4,:) = {'Elements',getappdata(0,'nel')};
                end
            end
        end

        % Redraw model and disable results
        redrawFlag = 0;
        set(mdata.popupmenu_Results,'value',1,'enable','off')
        set(mdata.text_Element,'string','Elements');
        set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(mdata.pushbutton_Textual,'enable','off')
        reactionValue = get(mdata.checkbox_Reactions,'value');
        set(mdata.checkbox_Reactions,'value',0,'enable','off')
        set(mdata.pushbutton_DynamicResults,'enable','off');
        set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
        
        % If user is not viewing model, redraw
        if get(mdata.popupmenu_Results,'value') ~= 1
            redraw(mdata)
            this.sizeFlag = draw.size;
            redrawFlag = 1;
        end
        
        % If model was not already redrawn, and reactions were enabled, redraw.
        if reactionValue ~= 0 && redrawFlag == 0
            delete(findobj('tag','drawReactions'));
            delete(findobj('tag','textForceReactions'));
            delete(findobj('tag','textMomentReactions'));
        end

        % Return variables to root
        setappdata(0,'model',model);
        setappdata(0,'draw',draw);
    else
        msgbox('There is already a node with these coordinates.', 'Error','error');
        return
    end
end

%--------------------------------------------------------------------------
function drawElements(this,tol)
    % Get node coordinates and flag if a new node was drawn
    [nodeCoords,nodeWasDrawn] = getElemNodes(this,tol);
    x = nodeCoords(1);
    y = nodeCoords(2);
    
    % Get updated nodes info
    nodes = getappdata(0,'nodes');

    % Set elemNodeID (id of the selected/created node)
    % Obs.: elemNode = node_i or node_f (1 or 2).
    this.elemNodeID(this.elemNode) = this.selectedNode;

    % Set element end coordinates
    this.elemCoords(this.elemNode,1) = x;
    this.elemCoords(this.elemNode,2) = y;

    % Check if initial and final nodes are the same (if two nodes have been
    % selected/created).
    if size(this.elemCoords,1) == 2
        if this.elemCoords(1,1) == this.elemCoords(2,1) && this.elemCoords(1,2) == this.elemCoords(2,2)
            this.elemNode = 0;
            this.elemNodeID = [];
            this.elemCoords = [];
            this.selectedNode = 0;
            this.moveAction();
            %msgbox('Initial and final nodes are the same', 'Error','error');
            return
        end
    end
    
    % Check if there is an element with the same end coords as the one being drawn
    if this.elemNode == 2 && nodeWasDrawn == false
        ni = nodes(this.elemNodeID(1));
        nf = nodes(this.elemNodeID(2));
        coincidentElem = auxModelFctn('isElemEqual',{ni,nf});
        if coincidentElem ~= false
            this.elemNode = 0;
            this.elemNodeID = [];
            this.elemCoords = [];
            this.selectedNode = 0;
            this.moveAction();
            %msgbox('There is already an element with this nodal incidence', 'Error','error');
            return
        end
    end

    % Get draw object
    draw = getappdata(0,'draw');
    
    % Size scale factor
    sz = draw.size;

    % Delete dynamic snap symbol (if there is one)
    if ~isempty(findobj('tag','snapNode'))
        delete(findobj('tag','snapNode'))
    end
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Get analysis model
    anm = get(mdata.popupmenu_Anm,'value');
    
    % Draw snap node symbol (on the selected/created node)
    if this.elemNode ~= 2
        if anm == 1 % TRUSS_2D
            circ = 0 : pi/50 : 2*pi;
            r = sz/125;
            xc = this.elemCoords(this.elemNode,1);
            yc = this.elemCoords(this.elemNode,2);
            xcirc = xc + r * cos(circ);
            ycirc = yc + r * sin(circ);
            plot(xcirc, ycirc, 'color', [0.9 0.2 0], 'tag','snapNode');
            hold on
        elseif anm == 2 || anm == 3 % FRAME_2D or GRILLAGE
            s = sz/200;
            xq = this.elemCoords(this.elemNode,1);
            yq = this.elemCoords(this.elemNode,2);
            xsq = [xq - s , xq + s , xq + s , xq - s];
            ysq = [yq - s , yq - s , yq + s , yq + s];
            fill(xsq, ysq, [0.9 0.2 0],'tag','snapNode');
            hold on
        end
        
    else % if this.elemNode == 2, create element objects, draw elements 
         % and reinitialize flags.
         
         % Enable process data pushbutton
         if strcmp(get(mdata.pushbutton_ProcessData,'Enable'),'off')
             set(mdata.pushbutton_ProcessData,'Enable','on')
         end
         
         % Update flags
         setappdata(0,'resultType',0);
         setappdata(0,'vis',1);
         
         % Get info to create Elem object
         draw = getappdata(0,'draw');
         model = getappdata(0,'model');
         mat = model.materials(1);
         sec = model.sections(1);
         nodes = getappdata(0,'nodes');
         ni = nodes(this.elemNodeID(1));
         nf = nodes(this.elemNodeID(2));
         coords = [ni.coord(1) ni.coord(2);
                   nf.coord(1) nf.coord(2)];
         
         % Check if new element crosses existing intersections (without
         % nodes on them). If so, draw nodes on their coordinates (if
         % intersectElem is on).
         if strcmp(this.intersectElem,'on')
             newIntSect = auxModelFctn('getCrossIntSectPoints',{coords,tol});
             if ~isempty(newIntSect)
                 for is = 1:size(newIntSect,1)
                     intersections = getappdata(0,'intersections');
                     drawNodes(this,{intersections(newIntSect(is,1)).coord,intersections(newIntSect(is,1)).elems})
                     aux = find((newIntSect(:,1)>newIntSect(is,1)));
                     newIntSect(aux,1) = newIntSect(aux,1) - 1;
                 end
                 intersections(newIntSect(:,1)) = [];
                 setappdata(0,'intersections',intersections)
                 nodes = getappdata(0,'nodes');
             end
         end
         
         % Get number of elements
         nel = getappdata(0,'nel');
         
         % Check if new element crosses existing nodes
         crossNodesOutput = auxModelFctn('getCrossNodePoints',{coords,tol});
         crossNodePoints = crossNodesOutput{1};
         collinearElems = crossNodesOutput{2};
         elemConnect = crossNodesOutput{3};
         
         % If elemConnect is empty, no elements need to be created
         if isempty(elemConnect)
             if strcmp(this.polyline,'on')
                 this.elemNode = 1;
                 this.elemNodeID = this.elemNodeID(2);
                 this.elemCoords = this.elemCoords(2,:);
                 this.selectedNode = this.elemNodeID;
             else
                 this.elemNode = 0;
                 this.elemNodeID = [];
                 this.elemCoords = [];
                 this.selectedNode = 0;
             end
             this.moveAction();
             return
         end
         
         % Check if new element crosses any existing elements
         newNodes = auxModelFctn('getNewNodes',{crossNodePoints,collinearElems,elemConnect,tol});
         
         % Set default hinges
         anm = get(mdata.popupmenu_Anm,'Value');
         if anm == 1 || anm == 4 % Truss
             hi = 0;
             hf = 0;
         else
             hi = 1;
             hf = 1;
         end
         
         % Get vector of handles to elem objects
         if nel ~= 0
             elems = getappdata(0,'elems');
         end
         
         % Create Elem objects and alocate them in elems vector
         for e = 1:size(elemConnect,1)
             % Default element type is Navier (type = 0)
             elem = Elem(0,model.anm,mat,sec,[nodes(elemConnect(e,1)) nodes(elemConnect(e,2))],hi,hf,[0 0 1]);
             nel = nel + 1;
             elems(nel) = elem;
             
             % Plot line between nodes
             X = [nodes(elemConnect(e,1)).coord(1), nodes(elemConnect(e,2)).coord(1)];
             Y = [nodes(elemConnect(e,1)).coord(2), nodes(elemConnect(e,2)).coord(2)];
             line(X,Y,'color',[0 0 0],'tag','drawElements');
             hold on
         end
         
         % Update model object
         model.elems = elems;
         model.nel = nel;
         draw.mdl = model;
         
         % Update variables in root
         setappdata(0,'model',model);
         setappdata(0,'draw',draw);
         setappdata(0,'nel',nel);
         setappdata(0,'elems',elems);
         
         if strcmp(get(mdata.elemIDButton,'Checked'),'on')
             draw.elementID();
         end
         if strcmp(get(mdata.orientationButton,'Checked'),'on')
             draw.elementOrientation();
         end
         
         % Check if there are crossing points and divide elements
         if ~isempty(newNodes)
             for nn = 1:size(newNodes,1)
                 % Get crossing point coordinates
                 x = newNodes(nn,1);
                 y = newNodes(nn,2);
                 z = newNodes(nn,3);
                 
                 % Get id of elements to be divided
                 whichElems = auxModelFctn('isPointInElem',{[x y z],tol});
                 
                 if strcmp(this.intersectElem,'on')
                     % Draw and create new nodes on the intersections
                     drawNodes(this,{[x y z],whichElems})
                 else
                     % Update vector of handles to intersections (does not
                     % create new node)
                     if ~isempty(whichElems)
                         intersections = getappdata(0,'intersections');
                         existingIntSect = 0;
                         for nis = 1:size(intersections,2)
                             if norm(intersections(nis).coord - [x y z]) <= 2*tol
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
             % Get updated number of elements
             nel = getappdata(0,'nel');
         end
         
         % Enable/disable solve intersections pushbutton (toolbar)
         if size(getappdata(0,'intersections'),2) >= 1
             set(mdata.pushbutton_SolveIntSects,'enable','on')
         else
             set(mdata.pushbutton_SolveIntSects,'enable','off')
         end
         
         % Update information panel in GUI_Main
         infoPanelData = get(mdata.uitable_infoPanel,'Data');
         infoPanelData(4,:) = {'Elements',nel};
         set(mdata.uitable_infoPanel,'Data',infoPanelData);
         
         % Update originalData property
         if ~isempty(this.originalData)
             this.originalData(4,:) = {'Elements',nel};
         end
         
         % Draw model and disable results
         redrawFlag = 0;
         set(mdata.popupmenu_Results,'value',1,'Enable','off')
         set(mdata.text_Element,'string','Elements');
         set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
         set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
         set(mdata.pushbutton_Textual,'enable','off')
         reactionValue = get(mdata.checkbox_Reactions,'value');
         set(mdata.checkbox_Reactions,'value',0,'enable','off')
         set(mdata.pushbutton_DynamicResults,'enable','off');
         set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
         
         % If user is not viewing model, redraw.
         if get(mdata.popupmenu_Results,'value') ~= 1
             redraw(mdata);
             this.sizeFlag = draw.size;
             redrawFlag = 1;
         end
         
         % If model was not already redrawn, and reactions were enabled,
         % redraw.
         if reactionValue ~= 0 && redrawFlag == 0
%              redraw(mdata)
%              this.sizeFlag = draw.size;
            delete(findobj('tag','drawReactions'));
            delete(findobj('tag','textForceReactions'));
            delete(findobj('tag','textMomentReactions'));
         end
         
         % Delete dynamic snap node symbols
         delete(findobj('tag', 'selectedNode'));
         delete(findobj('tag', 'snapNode'))
         delete(findobj('tag', 'snapNode2'))
         delete(findobj('tag', 'snapElem'))
         delete(findobj('tag', 'snapElem2'))
         
         % Reinitialize flags
         if strcmp(this.polyline,'on')
             this.elemNode = 1;
             this.elemNodeID = this.elemNodeID(2);
             this.elemCoords = this.elemCoords(2,:);
             this.selectedNode = this.elemNodeID;
         else
             this.elemNode = 0;
             this.elemNodeID = [];
             this.elemCoords = [];
             this.selectedNode = 0;
         end
         this.orthoCoords = [];
    end
end

%--------------------------------------------------------------------------
function drawElements3D(this,tol)
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Check if new nodes were drawn
    nodeWasDrawn = false;
    for nID = 1:size(this.elemNodeID,2)
        if this.elemNodeID(nID) < 0
            nodeWasDrawn = true;
            break
        end
    end
    this.elemNodeID = abs(this.elemNodeID);

    % Check if initial and final nodes are the same (if two nodes have been
    % selected/created).
    if this.elemCoords(1,1) == this.elemCoords(2,1) &&...
       this.elemCoords(1,2) == this.elemCoords(2,2) &&...
       this.elemCoords(1,3) == this.elemCoords(2,3)
        if strcmp(this.polyline,'on')
            this.elemNode = 1;
            this.elemNodeID = this.elemNodeID(2);
            this.elemCoords = this.elemCoords(2,:);
            this.selectedNode = this.elemNodeID;
        else
            this.elemNode = 0;
            this.elemNodeID = [];
            this.elemCoords = [];
            this.selectedNode = 0;
        end
        this.moveAction();
        %msgbox('Initial and final nodes are the same', 'Error','error');
        return
    end
    
    % Check if there is an element with the same end coords as the one being drawn
    if nodeWasDrawn == false
        ni = nodes(this.elemNodeID(1));
        nf = nodes(this.elemNodeID(2));
        coincidentElem = auxModelFctn('isElemEqual',{ni,nf});
        if coincidentElem ~= false
            if strcmp(this.polyline,'on')
                this.elemNode = 1;
                this.elemNodeID = this.elemNodeID(2);
                this.elemCoords = this.elemCoords(2,:);
                this.selectedNode = this.elemNodeID;
            else
                this.elemNode = 0;
                this.elemNodeID = [];
                this.elemCoords = [];
                this.selectedNode = 0;
            end
            this.moveAction();
            %msgbox('There is already an element with this nodal incidence', 'Error','error');
            return
        end
    end
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));

    % Enable process data pushbutton
    if strcmp(get(mdata.pushbutton_ProcessData,'Enable'),'off')
        set(mdata.pushbutton_ProcessData,'Enable','on')
    end
    
    % Update flags
    setappdata(0,'resultType',0);
    setappdata(0,'vis',1);
    
    % Get canvas borders
    dfltUnits = get(mdata.axes_Canvas,'units');
    set(mdata.axes_Canvas,'units','normalized');
    limits = get(mdata.axes_Canvas,'Position');
    set(mdata.axes_Canvas,'units',dfltUnits);
    axisWidth = limits(3);
    
    % Get info to create Elem object
    draw = getappdata(0,'draw');
    model = getappdata(0,'model');
    mat = model.materials(1);
    sec = model.sections(1);
    nodes = getappdata(0,'nodes');
    ni = nodes(this.elemNodeID(1));
    nf = nodes(this.elemNodeID(2));
    coords = [ni.coord(1) ni.coord(2) ni.coord(3);
              nf.coord(1) nf.coord(2) nf.coord(3)];
        
    % Check if new element crosses existing intersections (without
    % nodes on them). If so, draw nodes on their coordinates (if
    % intersectElem is on).
    if strcmp(this.intersectElem,'on')
        newIntSect = auxModelFctn('getCrossIntSectPoints',...
                                 {coords,(axisWidth/20)/this.currentZoom});
        if ~isempty(newIntSect)
            for is = 1:size(newIntSect,1)
                intersections = getappdata(0,'intersections');
                drawNodes(this,{intersections(newIntSect(is,1)).coord,...
                                    intersections(newIntSect(is,1)).elems})
                aux = find((newIntSect(:,1)>newIntSect(is,1)));
                newIntSect(aux,1) = newIntSect(aux,1) - 1;
            end
            intersections(newIntSect(:,1)) = [];
            setappdata(0,'intersections',intersections)
            nodes = getappdata(0,'nodes');
        end
    end
    
    % Get number of elements
    nel = getappdata(0,'nel');
    
    % Check if new element crosses existing nodes
    crossNodesOutput = auxModelFctn('getCrossNodePoints',...
                                 {coords,(axisWidth/20)/this.currentZoom});
    crossNodePoints = crossNodesOutput{1};
    collinearElems = crossNodesOutput{2};
    elemConnect = crossNodesOutput{3};
    
    % If elemConnect is empty, no elements need to be created
    if isempty(elemConnect)
        if strcmp(this.polyline,'on')
            this.elemNode = 1;
            this.elemNodeID = this.elemNodeID(2);
            this.elemCoords = this.elemCoords(2,:);
            this.selectedNode = this.elemNodeID;
        else
            this.elemNode = 0;
            this.elemNodeID = [];
            this.elemCoords = [];
            this.selectedNode = 0;
        end
        this.moveAction();
        return
    end
    
    % Check if new element crosses any existing elements
    newNodes = auxModelFctn('getNewNodes',{crossNodePoints,...
                                          collinearElems,elemConnect,tol});
    
    % Set default hinges
    anm = get(mdata.popupmenu_Anm,'Value');
    if anm == 1 || anm == 4 % Truss
        hi = 0;
        hf = 0;
    else
        hi = 1;
        hf = 1;
    end
    
    % Get vector of handles to elem objects
    if nel ~= 0
        elems = getappdata(0,'elems');
    end
    
    % Create Elem objects and alocate them in elems vector
    for e = 1:size(elemConnect,1)
        % Verify if Vz is in the same direction of local axis X:
        % Get nodal coordinates
        xi = nodes(elemConnect(e,1)).coord(1);
        yi = nodes(elemConnect(e,1)).coord(2);
        zi = nodes(elemConnect(e,1)).coord(3);
        xf = nodes(elemConnect(e,2)).coord(1);
        yf = nodes(elemConnect(e,2)).coord(2);
        zf = nodes(elemConnect(e,2)).coord(3);
        % Calculate element local axis X orientation versor
        x = [xf-xi, yf-yi, zf-zi];
        vz = [0, 0, 1];
        % Compare vectors 'x' and 'vz'
        w = cross(x,vz);
        if (abs(w(1)) < 1e-10) && (abs(w(2)) < 1e-10) && (abs(w(3)) < 1e-10)
            vz = [0 -1 0];
        end
        % Default element type is Navier (type = 0)
        elem = Elem(0,model.anm,mat,sec,[nodes(elemConnect(e,1)),...
                nodes(elemConnect(e,2))],hi,hf,vz);
        nel = nel + 1;
        elems(nel) = elem;
        
        % Plot line between nodes
        X = [nodes(elemConnect(e,1)).coord(1), nodes(elemConnect(e,2)).coord(1)];
        Y = [nodes(elemConnect(e,1)).coord(2), nodes(elemConnect(e,2)).coord(2)];
        Z = [nodes(elemConnect(e,1)).coord(3), nodes(elemConnect(e,2)).coord(3)];
        plot3(X,Y,Z,'color',[0 0 0],'tag','drawElements');
        hold on
    end
    
    % Update model object
    model.elems = elems;
    model.nel = nel;
    draw.mdl = model;
    
    % Update variables in root
    setappdata(0,'model',model);
    setappdata(0,'draw',draw);
    setappdata(0,'nel',nel);
    setappdata(0,'elems',elems);
    
    if strcmp(get(mdata.elemIDButton,'Checked'),'on')
        draw.elementID();
    end
    if strcmp(get(mdata.orientationButton,'Checked'),'on')
        draw.elementOrientation();
    end
    
    % Check if there are crossing points and divide elements
    if ~isempty(newNodes)
        for nn = 1:size(newNodes,1)
            % Get crossing point coordinates
            x = newNodes(nn,1);
            y = newNodes(nn,2);
            z = newNodes(nn,3);
            
            % Get id of elements to be divided
            whichElems = auxModelFctn('isPointInElem',{[x y z],...
                                         (axisWidth/20)/this.currentZoom});
            
            if strcmp(this.intersectElem,'on')
                % Draw and create new nodes on the intersections
                drawNodes(this,{[x y z],whichElems})
            else
                % Update vector of handles to intersections (does not
                % create new nodes)
                if ~isempty(whichElems)
                    intersections = getappdata(0,'intersections');
                    existingIntSect = 0;
                    for nis = 1:size(intersections,2)
                        if norm(intersections(nis).coord - [x y z]) <=...
                                          2*(axisWidth/20)/this.currentZoom
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
        % Get updated number of elements
        nel = getappdata(0,'nel');
    end
    
    % Enable/disable solve intersections pushbutton (toolbar)
    if size(getappdata(0,'intersections'),2) >= 1
        set(mdata.pushbutton_SolveIntSects,'enable','on')
    else
        set(mdata.pushbutton_SolveIntSects,'enable','off')
    end
    
    % Update information panel in GUI_Main
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(4,:) = {'Elements',nel};
    set(mdata.uitable_infoPanel,'Data',infoPanelData);
    
    % Update originalData property
    if ~isempty(this.originalData)
        this.originalData(4,:) = {'Elements',nel};
    end
    
    % Draw model and disable results
    set(mdata.popupmenu_Results,'value',1,'Enable','off')
    set(mdata.text_Element,'string','Elements');
    set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(mdata.pushbutton_Textual,'enable','off')
    reactionValue = get(mdata.checkbox_Reactions,'value');
    set(mdata.checkbox_Reactions,'value',0,'enable','off')
    set(mdata.pushbutton_DynamicResults,'enable','off');
    set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    % If user is not viewing model, redraw.
    if get(mdata.popupmenu_Results,'value') ~= 1
        redraw(mdata,'Loads');
        this.sizeFlag = draw.size;
    end
    
    % If model was not already redrawn, and reactions were enabled,
    % redraw.
    if reactionValue ~= 0
        delete(findobj('tag','drawReactions'));
        delete(findobj('tag','textForceReactions'));
        delete(findobj('tag','textMomentReactions'));
    end
    
    % Reinitialize flags
    if strcmp(this.polyline,'on')
        this.elemNode = 1;
        this.elemNodeID = this.elemNodeID(2);
        this.elemCoords = this.elemCoords(2,:);
        this.selectedNode = this.elemNodeID;
    else
        this.elemNode = 0;
        this.elemNodeID = [];
        this.elemCoords = [];
        this.selectedNode = 0;
    end
    this.orthoCoords = [];
end

%% ------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO drawElements.
% NOT CALLED BY auxMouseFctn.
%
% Get nodes for new elements.
% Input: this = Emouse object ; tol = graphic tolerance
function [nodeCoords,nodeWasDrawnFlag] = getElemNodes(this,tol)
    % Get nodes info
    nnp = getappdata(0,'nnp');
    nodes = getappdata(0,'nodes');
    
    % Initialize flag
    nodeWasDrawnFlag = false;
    
    % Get graphic tolerance for check if new node is inside element
    draw = getappdata(0,'draw');
    if isempty(this.sizeFlag)
        draw.setSize();
        this.sizeFlag = draw.size;
    elseif this.sizeFlag == 0
        draw.setSize();
        this.sizeFlag = draw.size;
    else
        draw.size = this.sizeFlag;
    end
    sz = draw.size;
    if draw.mdl.anm.analysis_type == 0 || draw.mdl.anm.analysis_type == 3 % TRUSS
        inElemTol = sz/125;
    else
        inElemTol = sz/400;
    end

    % Check if ortho is on
    orthoIsOn = false;
    if ~isempty(this.orthoCoords) && strcmp(this.ortho,'on')
        x = this.orthoCoords(1);
        y = this.orthoCoords(2);
        orthoIsOn = true;
    end
    
    % Check if user selected an existing node
    if this.whichNodeSnap ~= 0
        nodeNeedsToBeDrawn = false;
        n = this.whichNodeSnap;
        xns = nodes(n).coord(1);
        yns = nodes(n).coord(2);
        
        if orthoIsOn == true && this.elemNode == 2
            if abs(xns-x) <= tol && abs(yns-y) > tol
                x = xns;
                whichNode = auxModelFctn('isPointInNode',[x y 0]);
                if ~isempty(whichNode)
                    this.selectedNode = whichNode;
                else
                    nodeNeedsToBeDrawn = true;
                    this.selectedNode = nnp + 1;
                    whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    if ~isempty(whichElems)
                        this.whichElemSnap = whichElems(end);
                        this.elemPoint = [x y];
                    end
                end
            elseif abs(xns-x) > tol && abs(yns-y) <= tol
                y = yns;
                whichNode = auxModelFctn('isPointInNode',[x y 0]);
                if ~isempty(whichNode)
                    this.selectedNode = whichNode;
                else
                    nodeNeedsToBeDrawn = true;
                    this.selectedNode = nnp + 1;
                    whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    if ~isempty(whichElems)
                        this.whichElemSnap = whichElems(end);
                        this.elemPoint = [x y];
                    end
                end
            else
                this.selectedNode = n;
                x = xns;
                y = yns;
            end
        else
            this.selectedNode = n;
            x = xns;
            y = yns;
        end

        if nodeNeedsToBeDrawn == true
            % Draw nodes
            drawNodes(this,{[x y 0],whichElems});
            nodeWasDrawnFlag = true;
        end
        
    else % if nodeSnap is empty, user did not click on a node (a node will be created with the element)
        whichElems = []; % flag for snapped elements
        nodeNeedsToBeDrawn = true;

        % Draw new node on click point
        if this.whichIntSnap ~= 0 % selected an intersection
            intersections = getappdata(0,'intersections');
            xi = intersections(this.whichIntSnap).coord(1);
            yi = intersections(this.whichIntSnap).coord(2);
            if orthoIsOn == true && this.elemNode == 2
                if abs(xi-x) <= tol && abs(yi-y) > tol
                    x = xi;
                    whichNode = auxModelFctn('isPointInNode',[x y 0]);
                    if ~isempty(whichNode)
                        this.selectedNode = whichNode;
                        nodeNeedsToBeDrawn = false;
                    else
                        whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    end
                elseif abs(xi-x) > tol && abs(yi-y) <= tol
                    y = yi;
                    whichNode = auxModelFctn('isPointInNode',[x y 0]);
                    if ~isempty(whichNode)
                        this.selectedNode = whichNode;
                        nodeNeedsToBeDrawn = false;
                    else
                        whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    end
                else
                    whichElems = intersections(this.whichIntSnap).elems;
                    x = xi;
                    y = yi;
                end
            else
                whichElems = intersections(this.whichIntSnap).elems;
                x = xi;
                y = yi;
            end
        elseif this.whichElemSnap ~= 0 % point is inside an element
            xe = this.elemPoint(1);
            ye = this.elemPoint(2);
            if orthoIsOn == true && this.elemNode == 2
                if abs(xe-x) <= tol && abs(ye-y) > tol
                    x = xe;
                    whichNode = auxModelFctn('isPointInNode',[x y 0]);
                    if ~isempty(whichNode)
                        this.selectedNode = whichNode;
                        nodeNeedsToBeDrawn = false;
                    else
                        whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    end
                elseif abs(xe-x) > tol && abs(ye-y) <= tol
                    y = ye;
                    whichNode = auxModelFctn('isPointInNode',[x y 0]);
                    if ~isempty(whichNode)
                        this.selectedNode = whichNode;
                        nodeNeedsToBeDrawn = false;
                    else
                        whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                    end
                elseif abs(xe-x) <= tol && abs(ye-y) <= tol
                    e = this.whichElemSnap;
                    xe_i = this.elemCoords(1,1);
                    dx = (x - xe_i);
                    ye_i = this.elemCoords(1,2);
                    dy = (y - ye_i);
                    if dx > 0 && dy == 0
                        theta = 0;
                    elseif dx > 0 && dy > 0
                        theta = pi/4;
                    elseif dx == 0 && dy > 0
                        theta = pi/2;
                    elseif dx < 0 && dy > 0
                        theta = 3*pi/4;
                    elseif dx < 0 && dy == 0
                        theta = pi;
                    elseif dx < 0 && dy < 0
                        theta = -3*pi/4;
                    elseif dx == 0 && dy < 0
                        theta = -pi/2;
                    else % dx > 0 && dy < 0
                        theta = -pi/4;
                    end

                    L = norm([dx dy]);
                    if L <= norm([(xe - xe_i), (ye - ye_i)])
                        xe_f = xe_i + (L + tol)*cos(theta);
                        ye_f = ye_i + (L + tol)*sin(theta);
                    else
                        xe_f = xe_i + L*cos(theta);
                        ye_f = ye_i + L*sin(theta);
                    end

                    crPt = auxModelFctn('areElemsCrossed',{[xe_i, ye_i; xe_f, ye_f],e});
                    if isempty (crPt)
                        x = xe;
                        y = ye;
                    else
                        x = crPt(1);
                        y = crPt(2);
                    end
                    whichElems = e;
                else
                    x = xe;
                    y = ye;
                    whichElems = this.whichElemSnap;
                end
            else
                x = xe;
                y = ye;
                whichElems = this.whichElemSnap;
            end

        elseif strcmp(this.snapToGrid,'on') % get grid coordinates
            pos = snapToGridPosition(this);
            if isempty(pos)
                x = [];
                y = [];
                nodeNeedsToBeDrawn = false;
            else
                x = pos(1);
                y = pos(2);

                % Check if there is a node on selected grid point
                inNodeFlag = false;
                for n = 1:nnp
                    if nodes(n).coord(1) == x && nodes(n).coord(2) == y
                        this.selectedNode = n;
                        inNodeFlag = true;
                        nodeNeedsToBeDrawn = false;
                        break
                    end
                end

                % If there is no node on the selected grid point, check if it 
                % is inside an element.
                if inNodeFlag == false
                    % Check if grid point is inside an element
                    whichElems = auxModelFctn('isPointInElem',{[x y 0],inElemTol});
                end
            end

        elseif orthoIsOn == false % point is not inside any intersection or element and ortho is off
            currentPosition = this.getMouseProperty('CurrentPosition');
            if isempty(currentPosition)
                x = [];
                y = [];
                nodeNeedsToBeDrawn = false;
            else
                x = currentPosition(1);
                y = currentPosition(2);
                whichNode = auxModelFctn('isPointInNode',[x y 0]);
                if ~isempty(whichNode)
                    this.selectedNode = whichNode;
                    nodeNeedsToBeDrawn = false;
                end
            end
        end

        if nodeNeedsToBeDrawn == true
            this.selectedNode = nnp + 1;
            % Draw nodes
            drawNodes(this,{[x y 0],whichElems});
            nodeWasDrawnFlag = true;
        end
    end
    
    nodeCoords = [x y];
end
