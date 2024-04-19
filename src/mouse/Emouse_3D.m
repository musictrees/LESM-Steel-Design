%% Emouse_3D class
%
%% Description
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <emouse.html *Emouse*> in the <main.html LESM (Linear Elements
% Structure Model)> program. This sub-class implements abstract methods,
% defined in super-class *Emouse*, that deal with 3D models.
%
classdef Emouse_3D < Emouse
    %%
    % <emouse.html See documentation on *Emouse* super-class>.
  
    %% Public attributes
    properties (Access = public)
        
        % Zoom variables
        zoomIndice = [];         % Zoom indice.
        axesPos = [];            % Actual axes position on window coordinates.
        defaultAxesPos = [];     % Initial axes position on window coordinates
        zoomGrid = [];           % Vector with grid coordinates.
        defaultZoomGrid = [];    % Vetor with initial grid coordinates.
        zoomSteps = 0;           % Counter for initial steps.
        defaultZoomSteps = 0;    % Counter for steps.
        currentZoom = 1;
        
        % Rotate variables
        rotIniX = 0;             % x window coordinate when RB button is clicked to rotate.
        rotIniY = 0;             % y window coordinate when RB button is clicked to rotate.
        rotIniAz = -37.5;        % Azimute when RB button is clicked to rotate.
        rotIniEl = 30;           % Elevation when RB button is clicked to rotate.
        
        % Pan variables
        panIniX = 0;             % x window coordinate when RB button is clicked to drag.
        panIniY = 0;             % y window coordinate when RB button is clicked to drag.
        
        % Double click variables
        originalAxesPos = [];    % Original axes position on window coordinates.
        originalXLim = [];       % Original limits in axis x.
        originalYLim = [];       % Original limits in axis y.
        originalZLim = [];       % Original limits in axis z.
        originalAz = -37.5;      % Original azimute of the view.
        originalEl = 30;         % Original elevation of the view.
        
        % Mouse cursor variable
        mouseCursor = 'on';     % 'on','off', flag for pressed or unpressed cursor toggle button.
        
        % Draw variables
        sizeFlag = 0;           % flag for size change
        
        % Draw node variables
        drawNode = 'off';        % 'on','off', flag for pressed or unpressed draw nodes toggle button.
        
        % Draw element variables
        drawElem = 'off';        % 'on','off', flag for pressed or unpressed draw elements toggle button.
        polyline = 'off';        % 'on','off', flag for pressed or unpressed draw multiple elements toggle button.
        ortho = 'off';           % 'on','off', flag for pressed or unpressed ortho toggle button.
        elemNode = 0;            % flag for node selection on draw elements mode (0 = empty, 1 = node_i, 2 = node_f)
        elemNodeID = [];         % element initial and final nodes
        elemCoords = [];         % element nodal coordinates
        orthoCoords = [];
        
        % Node snap variables
        nodeSnap = [];           % node coordinates [x,y]
        whichNodeSnap = 0;       % node id
        
        % Element snap variables
        elemSnap = [];           % element end coordinates [node_i(x) node_i(y)]
                                 %                         [node_f(x) node_f(y)]                        
        whichElemSnap = 0;       % element id
        elemPoint = [];          % coords of specific point inside an element (x,y)
        
        % Intersect elements variables
        intersectElem = 'off';   % 'on','off', flag for pressed or unpressed intersect elements toggle button.
        whichIntSnap = 0;
        
        % Snap to grid precision
        snapToGrid = 'off';      % 'on','off', flag for pressed or unpressed snap to grid toggle button.
        snapPrecision = 1;
        
        % Selected canvas entity variables
        selectedNode = 0;        % node id
        selectedElem = 0;        % element id
        elemResults = [];
        originalData = {};       % info panel uitable model data
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Emouse_3D(fig,axes)
            this = this@Emouse(fig,axes);
            this.zoomGrid = unique(round(logspace(0,5,51)));
            this.zoomGrid(this.zoomGrid < 10) = [];
            this.zoomSteps = length(this.zoomGrid);
            this.zoomIndice = find(this.zoomGrid == 100);
            
            % Properties for double click on mouse
            originalUnits = get(axes, 'Units');
            set(axes, 'Units', 'pixels');
            this.originalAxesPos = get(axes, 'Position');
            set(axes, 'Units', originalUnits);
            this.originalXLim = get(axes, 'XLim');
            this.originalYLim = get(axes, 'YLim');
            this.originalZLim = get(axes, 'ZLim');
            [this.originalAz,this.originalEl] = view(axes);
        end
    end
    
    %% Abstract methods
    % Implementation of the abstract methods declared in super-class <Emouse.html *Emouse*>.
    methods        
        %------------------------------------------------------------------
        % Action executed when an mouse button is pressed.
        function downAction(this)
            % Get handle to GUI_Main
            mdata = guidata(findobj('Tag','GUI_Main'));
            anm = get(mdata.popupmenu_Anm,'value');
            
            % Get toggle buttons states ('on' or 'off')
            this.mouseCursor   = get(mdata.togglebutton_Cursor,'state');
            this.drawNode      = get(mdata.togglebutton_Node,'state');
            this.drawElem      = get(mdata.togglebutton_Element,'state');
            this.snapToGrid    = get(mdata.togglebutton_SnapToGrid,'state');
            this.ortho         = get(mdata.togglebutton_Ortho,'state');
            this.intersectElem = get(mdata.togglebutton_CrossElements,'state');
            this.polyline      = get(mdata.togglebutton_Polyline,'state');
            panState           = get(mdata.toolbar_pan,'state');
            
            if strcmp(this.whichMouseButton,'left')
                setappdata(0,'isDrawingDynamic',false)
            end
            
            %--------------------------------------------------------------
            % Draw nodes
            if strcmp(this.whichMouseButton,'left') && strcmp(this.drawNode,'on')
                coords = [];
                if this.whichNodeSnap ~= 0
                    nodes = getappdata(0,'nodes');
                    coords = nodes(this.whichNodeSnap).coord;
                    inWhichElems = [];
                elseif this.whichIntSnap ~= 0
                    intSects = getappdata(0,'intersections');
                    coords = intSects(this.whichIntSnap).coord;
                    inWhichElems = intSects(this.whichIntSnap).elems;
                elseif this.whichElemSnap ~= 0
                    coords = this.elemPoint;
                    inWhichElems = this.whichElemSnap;
                elseif anm == 3
                    coords = auxMouseFctn('spatial2Plane',this,{'z',0});
                    inWhichElems = [];
                    if ~isempty(coords)
                        cnvsXLim = get(mdata.axes_Canvas,'XLim');
                        cnvsYLim = get(mdata.axes_Canvas,'YLim');
                        if coords(1) < cnvsXLim(1) || coords(1) > cnvsXLim(2) ||...
                           coords(2) < cnvsYLim(1) || coords(2) > cnvsYLim(2)
                            coords = [];
                        end
                    end
                    if strcmp(this.snapToGrid,'on') && ~isempty(coords)
                        % Get closest grid point to cursor
                        coordsXY = auxMouseFctn('snapToGridPosition',this,coords);
                        coords = [coordsXY, 0];
                        
                        % Get graphic tolerance
                         draw = getappdata(0,'draw');
                         if isempty(this.sizeFlag) || this.sizeFlag == 0
                             draw.setSize();
                             this.sizeFlag = draw.size;
                         else
                             draw.size = this.sizeFlag;
                         end
                         sz = draw.size;
                         if draw.mdl.anm.analysis_type == 1 || draw.mdl.anm.analysis_type == 3 % TRUSS
                             tol = sz/125;
                         else
                             tol = sz/400;
                         end
                        inWhichElems = auxModelFctn('isPointInElem',{coords,tol});
                    end
                end
                
                if ~isempty(coords)
                    [~] = auxMouseFctn('drawNodes',this,{coords,inWhichElems});
                end
                
            % Draw elements
            elseif strcmp(this.whichMouseButton,'left') && strcmp(this.drawElem,'on')
                % Update element node flag (initial or final node)
                this.elemNode = this.elemNode + 1;
                
                % Check if property elemNode is 1 (initial node) or 2
                % (final node). If not, reinitialize properties and return.
                if this.elemNode ~= 1 && this.elemNode ~= 2
                    this.elemNode = 0;
                    this.elemNodeID = [];
                    this.elemCoords = [];
                    this.selectedNode = 0;
                    delete(findobj('tag', 'selectedNode'))
                    this.moveAction();
                    return
                end
                
                % Get 3D click line
                lineCoords = (this.currentPosition)';
                
                % Get axis limits
                axLim = [this.canvas.XLim;this.canvas.YLim;this.canvas.ZLim];
                
                % Get canvas borders
                dfltUnits = get(gca,'units');
                set(gca,'units','normalized');
                limits = get(gca,'Position');
                set(gca,'units',dfltUnits);
                axisWidth = limits(3);
                
                % Set graphic tolerance
                scl = max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))])/2;
                tol = (axisWidth/20) / this.currentZoom * 10^(floor(log10(scl)));
                
                % Check if 3D click line is less than 1 (avoids numeric
                % problems with the graphic tolerance)
                if norm(lineCoords(1,:)-lineCoords(2,:)) < 1 &&...
                   max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))]) >= 1
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(1,:) = lineCoords(2,:) - [cx cy cz];
                end
                
                % If model is grillage, extend 3D click line beyond z=0 to avoid numeric problems.
                if anm == 3
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(2,:) = lineCoords(2,:) + [cx cy cz];
                end
                
                [~] = auxMouseFctn('selectNodes3D',this,{lineCoords,tol});
                
                % Set elemNodeID (id of the selected/created node)
                % Obs.: elemNode = node_i or node_f (1 or 2).
                this.elemNodeID(this.elemNode) = this.selectedNode;
                
                % If model is grillage, set current position in plane XY
                if anm == 3
                    coordsXY = auxMouseFctn('spatial2Plane',this,{'z',0});
                    if ~isempty(coordsXY)
                        cnvsXLim = get(mdata.axes_Canvas,'XLim');
                        cnvsYLim = get(mdata.axes_Canvas,'YLim');
                        if coordsXY(1) < cnvsXLim(1) || coordsXY(1) > cnvsXLim(2)...
                           || coordsXY(2) < cnvsYLim(1) || coordsXY(2) > cnvsYLim(2)
                            coordsXY = [];
                        end
                    end
                    this.currentPosition = coordsXY;
                end
                
                % Set element end coordinates
                [~] = auxMouseFctn('setElemNodes3D',this,tol);
                
                % Reset currentPosition property
                if anm == 3
                    pt = get(this.canvas, 'CurrentPoint');
                    xP = pt(:, 1);
                    yP = pt(:, 2);
                    zP = pt(:, 3);
                    this.currentPosition = [xP yP zP]';
                end
                
                % Draw elements and new nodes (if there are any)
                if this.elemNode == 2
                    [~] = auxMouseFctn('drawElements3D',this,tol);
                    delete(findobj('tag', 'selectedNode'));
                    delete(findobj('tag', 'snapNode'))
                    delete(findobj('tag', 'snapElem'))
                    if this.selectedNode ~= 0
                        nodes = getappdata(0,'nodes');
                        draw = getappdata(0,'draw');
                        coords = nodes(this.selectedNode).coord;
                        draw.cube(coords(1), coords(2), coords(3), this.sizeFlag/75, [1 0 0],'selectedNode');
                    end
                end
                
            % Select nodes/elems
            elseif strcmp(this.whichMouseButton,'left') && strcmp(this.mouseCursor,'on')
                % Reinitialize selection properties
                this.selectedNode = 0;
                this.selectedElem = 0;
                
                % Clear previous selected node mark
                delete(findobj('tag', 'selectedNode'))
                delete(findobj('tag', 'selectedElem'))
                
                % Get 3D click line
                lineCoords = (this.currentPosition)';
                
                % Get axis limits
                axLim = [this.canvas.XLim;this.canvas.YLim;this.canvas.ZLim];
                
                % Get canvas borders
                dfltUnits = get(this.canvas,'units');
                set(this.canvas,'units','normalized');
                limits = get(this.canvas,'Position');
                set(this.canvas,'units',dfltUnits);
                axisWidth = limits(3);
                
                % Set graphic tolerance
                scl = max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))])/2;
                tol = (axisWidth/20) / this.currentZoom * 10^(floor(log10(scl)));
                
                % Check if 3D click line is less than 1 (avoids numeric problems with the graphic tolerance)
                if norm(lineCoords(1,:)-lineCoords(2,:)) < 1 &&...
                   max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))]) >= 1
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(1,:) = lineCoords(2,:) - [cx cy cz];
                end
                
                % If model is grillage, extend 3D click line beyond z=0 to avoid numeric problems
                if anm == 3
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(2,:) = lineCoords(2,:) + [cx cy cz];
                end
                
                [~] = auxMouseFctn('selectNodes3D',this,{lineCoords,tol});
                
                % If no nodes were selected, check if user clicked on an element
                if this.selectedNode == 0
                    [~] = auxMouseFctn('selectElems3D',this,{lineCoords,tol});
                    this.elemResults = [];
                end
            end
            
            % Rotate
            if strcmp(this.whichMouseButton,'right')
                % Catches the inicial coordinates on window.
                crd = get(this.dialog, 'CurrentPoint');
                cx = crd(1); 
                cy = crd(2);
                [az,el] = view(this.canvas);

                % Incorporate the values on the variables.
                this.rotIniAz = az;
                this.rotIniEl = el;
                this.rotIniX = cx;
                this.rotIniY = cy;
                
                % Set plane visualization button to 3D
                set(mdata.plane3D,'Checked','on')
                set(mdata.planeXY,'Checked','off')
                set(mdata.planeXZ,'Checked','off')
                set(mdata.planeYZ,'Checked','off')
            end
            
            % Pan
            if strcmp(this.whichMouseButton,'center') ||...
              (strcmp(panState,'on') && strcmp(this.whichMouseButton,'left'))
                % Catches the inicial coordinates on window.
                wpt = get(this.dialog, 'CurrentPoint');
                wx = wpt(1); 
                wy = wpt(2);
                
                % Incorporate the values on the variables.
                this.panIniX = wx;
                this.panIniY = wy;
            end
            
            % Double click
            if strcmp(this.whichMouseButton,'double click')
                this.doubleClick()
            end
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse pointer moves on the canvas.
        function moveAction(this)
            % Get handle to GUI_Main
            mdata = guidata(findobj('Tag','GUI_Main'));
            
            % Get toolbar_pan's state ('on' or 'off')
            panState = get(mdata.toolbar_pan,'state');
            
            % Rotate
            if strcmp(this.whichMouseButton,'right')
                % Catches the actual coordinates on window.
                wpt = get(this.dialog, 'CurrentPoint');
                x = wpt(1); 
                y = wpt(2);
                
                % Calculate the displacements of azimute and elevation.
                dAz = this.rotIniX - x;
                dEl = this.rotIniY - y;
                if dAz > -10 && dAz < 10
                    dAz = 0;
                end
                
                % Incorporate the displacements calculated before in the
                % inicial azimute and elevation.
                az = this.rotIniAz + dAz;
                el = this.rotIniEl + dEl;
                
                % Checking the elevation values and changing them when its 
                % necessary.
                if (el > 90)
                    el = 90;
                elseif (el < -90)
                    el = -90;
                end
                
                % Sets the new view.
                view(this.canvas,[az,el]);
            end
            
            % Pan
            if strcmp(this.whichMouseButton,'center') ||...
              (strcmp(panState,'on') && strcmp(this.whichMouseButton,'left'))
                % Catches the actual window coordinates.
                wpt = get(this.dialog, 'CurrentPoint');
                wx = wpt(1); 
                wy = wpt(2);
                
                % Calculate the displacements.
                dX = this.panIniX - wx;
                dY = this.panIniY - wy;
                
                % Gets the object position.
                originalUnits = get(this.canvas, 'Units');
                set(this.canvas, 'Units', 'pixels');
                this.axesPos = get(this.canvas, 'Position');
                set(this.canvas, 'Units', originalUnits);
                
                % Catches the new limits position.
                this.axesPos(1) = this.axesPos(1) - dX;
                this.axesPos(2) = this.axesPos(2) - dY;
                
                % Sets the new limits position.
                dfltUnits = get(this.canvas, 'Units');
                set(this.canvas, 'Units', 'pixels');
                set(this.canvas, 'Position', this.axesPos);
                set(this.canvas, 'Units', dfltUnits);
                
                % Changes the old window coordinates.
                this.panIniX = wx;
                this.panIniY = wy;
            end
            
            %------------------------------------------------------ SNAP
            % Snap to nodes/elements while modelling
            if strcmp(this.drawNode,'on') || strcmp(this.drawElem,'on')
                % Reinitialize snap properties
                this.nodeSnap = [];           % node coordinates [x,y]
                this.whichNodeSnap = 0;       % node id
                
                % Element snap variables
                this.elemSnap = [];    % element end coordinates [node_i(x) node_i(y)]
                                       %                         [node_f(x) node_f(y)]
                this.whichElemSnap = 0;   % element id
                this.elemPoint = [];      % coords of specific point inside an element (x,y)
                
                % Snap to intersection property
                this.whichIntSnap = 0;

                % Clear previous selected node mark
                delete(findobj('tag', 'snapNode'))
                delete(findobj('tag', 'snapElem'))
                delete(findobj('tag', 'snapIntSect'))
                delete(findobj('tag', 'snapGrid'))
                delete(findobj('tag', 'dynamicLine'))
                
                % Get vector of handles to node/intersection objects
                nodes = getappdata(0,'nodes');
                intSects = getappdata(0,'intersections');
                
                % Get draw object
                draw = getappdata(0,'draw');
                if ~isempty(draw)
                    sz = draw.size;
                else
                    sz = 5;
                end
                
                % Get 3D click line
                lineCoords = (this.currentPosition)';
                
                % Get axis limits
                axLim = [this.canvas.XLim;this.canvas.YLim;this.canvas.ZLim];
                
                % Get canvas borders
                dfltUnits = get(this.canvas,'units');
                set(this.canvas,'units','normalized');
                limits = get(this.canvas,'Position');
                set(this.canvas,'units',dfltUnits);
                axisWidth = limits(3);
                
                % Set graphic tolerance
                scl = max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))])/2;
                tol = (axisWidth/20) / this.currentZoom * 10^(floor(log10(scl)));
                
                % Check if 3D click line is less than 1 (avoids numeric problems with the graphic tolerance)
                if norm(lineCoords(1,:)-lineCoords(2,:)) < 1 &&...
                   max([diff(axLim(1,:)),diff(axLim(2,:)),diff(axLim(3,:))]) >= 1
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(1,:) = lineCoords(2,:) - [cx cy cz];
                end
                
                % If model is grillage, extend 3D click line beyond z=0 to avoid numeric problems.
                if get(mdata.popupmenu_Anm,'value') == 3
                    dx = lineCoords(2,1) - lineCoords(1,1);
                    dy = lineCoords(2,2) - lineCoords(1,2);
                    dz = lineCoords(2,3) - lineCoords(1,3);
                    L = norm(lineCoords(1,:) - lineCoords(2,:));
                    cx = dx/L;
                    cy = dy/L;
                    cz = dz/L;
                    lineCoords(2,:) = lineCoords(2,:) + [cx cy cz];
                end
                
                if this.elemNode == 1
                    [~] = auxMouseFctn('snapToNodes3D',this,{lineCoords,tol});
                    coordsXY = [];
                    coords = [];
                    if this.whichNodeSnap ~= 0
                        coords = [this.elemCoords; nodes(this.whichNodeSnap).coord];
                        
                    else % If no nodes were snapped, check if cursor is over an intersection or element.
                        [~] = auxMouseFctn('snapToIntSects3D',this,{lineCoords,tol});
                        if this.whichIntSnap ~= 0
                            coords = [this.elemCoords; intSects(this.whichIntSnap).coord];
                        else
                            [~] = auxMouseFctn('snapToElems3D',this,{lineCoords,tol});
                            if this.whichElemSnap ~= 0
                                coords = [this.elemCoords; this.elemPoint];
                                      
                            elseif get(mdata.popupmenu_Anm,'value') == 3
                                coordsXY = auxMouseFctn('spatial2Plane',this,{'z',0});
                                if ~isempty(coordsXY)
                                    cnvsXLim = get(mdata.axes_Canvas,'XLim');
                                    cnvsYLim = get(mdata.axes_Canvas,'YLim');
                                    if coordsXY(1) < cnvsXLim(1) || coordsXY(1) > cnvsXLim(2) ||...
                                       coordsXY(2) < cnvsYLim(1) || coordsXY(2) > cnvsYLim(2)
                                        coordsXY = [];
                                    end
                                    if strcmp(this.snapToGrid,'on') && ~isempty(coordsXY)
                                        % Get closest grid point to cursor
                                        auxCoords = auxMouseFctn('snapToGridPosition',this,coordsXY);

                                        % Draw snap to grid dynamic symbol
                                        [~] = auxMouseFctn('snapToGridDraw',this,{auxCoords,[3,sz],true});
                                        coordsXY = [auxCoords, 0];
                                    end
                                end
                                coords = [this.elemCoords; coordsXY];
                            end
                        end
                    end
                    if get(mdata.popupmenu_Anm,'value') == 3 && strcmp(this.ortho,'on')
                        coords_2 = [];
                        if isempty(coordsXY)
                            coordsXY = auxMouseFctn('spatial2Plane',this,{'z',0});
                        end
                        if ~isempty(coordsXY)
                            cnvsXLim = get(mdata.axes_Canvas,'XLim');
                            cnvsYLim = get(mdata.axes_Canvas,'YLim');
                            if coordsXY(1) < cnvsXLim(1) || coordsXY(1) > cnvsXLim(2) ||...
                               coordsXY(2) < cnvsYLim(1) || coordsXY(2) > cnvsYLim(2)
                                coordsXY = [];
                            end
                        end
                        if ~isempty(coordsXY)
                            this.orthoCoords = auxMouseFctn('orthoPosition',this,coordsXY);
                            x = this.orthoCoords(1);
                            y = this.orthoCoords(2);
                            coords_2 = [x y 0];
                        else
                            this.orthoCoords = [];
                        end
                        if ~isempty(coords_2)
                        coords = [this.elemCoords; coords_2];
                        end
                    end
                    if ~isempty(coords)
                        plot3(coords(:,1),coords(:,2),coords(:,3),'color',[1 0 0],'tag','dynamicLine');
                    end
                else
                    [~] = auxMouseFctn('snapToNodes3D',this,{lineCoords,tol});
                    if this.whichNodeSnap == 0 % If no nodes were snapped, check if cursor is over an intersection or element.
                        [~] = auxMouseFctn('snapToIntSects3D',this,{lineCoords,tol});
                        if this.whichIntSnap == 0
                            [~] = auxMouseFctn('snapToElems3D',this,{lineCoords,tol});
                            
                            if this.whichElemSnap == 0
                                if get(mdata.popupmenu_Anm,'value') == 3 && strcmp(this.snapToGrid,'on')
                                    coordsXY = auxMouseFctn('spatial2Plane',this,{'z',0});
                                    if ~isempty(coordsXY)
                                        cnvsXLim = get(mdata.axes_Canvas,'XLim');
                                        cnvsYLim = get(mdata.axes_Canvas,'YLim');
                                        if coordsXY(1) < cnvsXLim(1) || coordsXY(1) > cnvsXLim(2) ||...
                                           coordsXY(2) < cnvsYLim(1) || coordsXY(2) > cnvsYLim(2)
                                            coordsXY = [];
                                        end
                                        if ~isempty(coordsXY)
                                            % Get closest grid point to cursor
                                            coords = auxMouseFctn('snapToGridPosition',this,coordsXY);

                                            % Draw snap to grid dynamic symbol
                                            [~] = auxMouseFctn('snapToGridDraw',this,{coords,[3,sz],true});
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Action executed when an mouse button is unpressed.
        function upAction(~)
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse scroll is used.
        % Executes zoom (camera) using the mouse scroll.
        function scrollAction(this,direction,wx,wy)           
            % Gets the object position.
            originalUnits = get(this.canvas, 'Units');
            set(this.canvas, 'Units', 'pixels');
            this.defaultAxesPos = get(this.canvas, 'Position');
            set(this.canvas, 'Units', originalUnits);
            this.axesPos = this.defaultAxesPos;
            
            % Check if canvas was captured or not.
            if isempty(this.canvas)
                return
            end

            % Check the direction and increment or decrement the value of
            % zoomIndice based on this.
            aux = this.zoomIndice;
            switch direction
                case 'plus'
                    if this.zoomIndice < this.zoomSteps
                        this.zoomIndice = this.zoomIndice + 1;
                    end
                case 'minus'
                    if this.zoomIndice > 1
                        this.zoomIndice = this.zoomIndice - 1;
                    end
            end

            % Picks the value that corresponds the "zoom index" of the
            % vector "zoom grid".
            zoomPct = 100 * this.zoomGrid(this.zoomIndice)/this.zoomGrid(aux);
            this.currentZoom = this.currentZoom * zoomPct/100;

            % Calculate the new height of axis x.
            dd_x = this.defaultAxesPos(3) * zoomPct/100 - this.axesPos(3);

            % Calculate the range in x.
            range_x = abs(diff([this.axesPos(1) wx]));

            % Normalize the range in x.
            cf_x = range_x/this.axesPos(3);

            % Calculate new axis x position.
            this.axesPos(3) = this.axesPos(3) + dd_x;
            this.axesPos(1) = this.axesPos(1) - dd_x * cf_x;

            % Calculate the new height of axis y.
            dd_y = this.defaultAxesPos(4) * zoomPct/100 - this.axesPos(4);

            % Calculate the range in y.
            range_y = abs(diff([this.axesPos(2) wy]));

            % Normalize the range in y.
            cf_y = range_y/this.axesPos(4);

            % Calculate new axis y position.
            this.axesPos(4) = this.axesPos(4) + dd_y;
            this.axesPos(2) = this.axesPos(2) - dd_y * cf_y;

            % Sets the new axes position, executing the zoom.
            if (this.axesPos(3) > 30 && this.axesPos(4) > 30 &&...
                this.axesPos(3) < 250000 && this.axesPos(4) < 250000)
                dfltUnits = get(this.canvas, 'Units');
                set(this.canvas, 'Units', 'pixels');
                set(this.canvas, 'Position', this.axesPos);
                set(this.canvas, 'Units', dfltUnits);
            end
        end
        
        %------------------------------------------------------------------
        function doubleClick(this)
            mdata = guidata(findobj('Tag','GUI_Main'));
            axes(this.canvas)
            
            % Reset zoom properties
            this.currentZoom = 1;
            this.zoomIndice = find(this.zoomGrid == 100);
            
            % Set axes limits
            set(this.canvas, 'XLim', this.originalXLim);
            set(this.canvas, 'YLim', this.originalYLim);
            set(this.canvas, 'ZLim', this.originalZLim);
            
            % Set axes labels
            xlabel('X');
            ylabel('Y');
            if get(mdata.popupmenu_Anm,'value') == 3 % grillage
                zlabel(' ');
                zticks(0);
                zticklabels({' '});
            else
                zlabel('Z');
            end
            
            % Set the object position
            dfltUnits = get(this.canvas, 'Units');
            set(this.canvas, 'Units', 'pixels');
            set(this.canvas, 'Position', this.originalAxesPos);
            set(this.canvas, 'Units', dfltUnits);
            
            % Sets the original parameters of visualization
            view(gca,[this.originalAz,this.originalEl]);
            
            set(mdata.plane3D,'Checked','on')
            set(mdata.planeXY,'Checked','off','enable','on')
            set(mdata.planeXZ,'Checked','off','enable','on')
            set(mdata.planeYZ,'Checked','off','enable','on')
        end
        
        %------------------------------------------------------------------
        % Returns solicited protected property
        function output = getMouseProperty(this,whichProperty)
            output = [];
            switch whichProperty
                case 'Dialog'
                    output = this.dialog;
                case 'Canvas'
                    output = this.canvas;
                case 'MainDialogName'
                    output = this.mainDialogName;
                case 'MouseButtonMode'
                    output = this.mouseButtonMode;
                case 'VerticalScrollCount'
                    output = this.verticalScrollCount;
                case 'ScrollAllowed'
                    output = this.scrollAllowed;
                case 'CurrentPosition'
                    output = this.currentPosition;
            end
        end
    end
end
