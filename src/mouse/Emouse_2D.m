%% Emouse_2D class
%
%% Description
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <emouse.html *Emouse*> in the <main.html LESM (Linear Elements
% Structure Model)> program. This sub-class implements abstract methods,
% defined in super-class *Emouse*, that deal with 2D models.
%
classdef Emouse_2D < Emouse
    %%
    % <emouse.html See documentation on *Emouse* super-class>.
  
    %% Public attributes
    properties (Access = public)
        
        % Zoom variable
        zoomIndiceX = [];        % Zoom indice for axis x.
        zoomIndiceY = [];        % Zoom indice for axis y.
        zoomGrid = [];           % Vector with grid coordinates.
        zoomSteps = 0;           % Counter for steps.
        currentZoom = 1;
        
        % Pan variables
        axesPos = [];
        panIniX = 0;            % x window coordinate when RB button is clicked to drag.
        panIniY = 0;            % y window coordinate when RB button is clicked to drag.
        
        % Double click variables
        originalXLim = [];       % Original limits in axis x.
        originalYLim = [];       % Original limits in axis y.
        
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
        orthoCoords = [];        % ortho coordinates
        
        % Node snap variables
        nodeSnap = [];           % node coordinates [x,y]
        whichNodeSnap = 0;       % node id
        
        % Element snap variables
        elemSnap = [];           % element end coordinates [node_i(x) node_i(y)]
                                 %                         [node_f(x) node_f(y)]                        
        whichElemSnap = 0;       % element id
        elemPoint = [];          % coords of specific point inside an element (x,y)
        
        % Intersect elements variables
        intersectElem = 'off';    % 'on','off', flag for pressed or unpressed intersect elements toggle button.
        intSnap = [];
        whichIntSnap = 0;
        
        % Snap to grid precision
        snapToGrid = 'off';      % 'on','off', flag for pressed or unpressed snap to grid toggle button.
        snapPrecision = 1;       % snap precision (user input - Default = 1)
        
        % Selected canvas entity variables
        selectedNode = 0;        % node id
        selectedElem = 0;        % element id
        elemResults = [];
        originalData = {};       % info panel uitable model data
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Emouse_2D(fig,axes)
            this = this@Emouse(fig,axes);
            this.zoomGrid = unique(round(logspace(0,5,51)));
            this.zoomGrid(this.zoomGrid < 10) = [];
            this.zoomSteps = length(this.zoomGrid);
            this.zoomIndiceX = find(this.zoomGrid == 100);
            this.zoomIndiceY = this.zoomIndiceX;
            
            % Properties for double click on mouse
            mdata = guidata(findobj('Tag','GUI_Main'));
            this.originalXLim = get(mdata.axes_Canvas, 'XLim');
            this.originalYLim = get(mdata.axes_Canvas, 'YLim');
        end
    end
    
    %% Abstract methods
    % Implementation of the abstract methods declared in super-class <Emouse.html *Emouse*>.
    methods
        %------------------------------------------------------------------
        % Action executed when an mouse button is pressed.
        % Executes pan (limits) initial actions when right mouse button is 
        % pressed.
        function downAction(this)
            %--------------------------------------------------------------
            % Get handle to GUI_Main
            mdata = guidata(findobj('Tag','GUI_Main'));
            
            % Get toggle buttons states ('on' or 'off')
            this.mouseCursor = get(mdata.togglebutton_Cursor,'state');
            this.drawNode = get(mdata.togglebutton_Node,'state');
            this.drawElem = get(mdata.togglebutton_Element,'state');
            this.snapToGrid = get(mdata.togglebutton_SnapToGrid,'state');
            this.ortho = get(mdata.togglebutton_Ortho,'state');
            this.intersectElem = get(mdata.togglebutton_CrossElements,'state');
            this.polyline = get(mdata.togglebutton_Polyline,'state');
            panState = get(mdata.toolbar_pan,'state');
            
            if strcmp(this.whichMouseButton,'left')
                setappdata(0,'isDrawingDynamic',false)
            end
            
            %--------------------------------------------------------------
            % Pan
            if strcmp(this.whichMouseButton,'center') ||...
               (strcmp(panState,'on') && strcmp(this.whichMouseButton,'left'))
                % Catches the inicial coordinates on window.
                crd = get(this.dialog, 'CurrentPoint');
                wx = crd(1); 
                wy = crd(2);
                
                % Incorporate the values on the variables.
                this.panIniX = wx;
                this.panIniY = wy;
            end
            %--------------------------------------------------------------
            % Double click
            if strcmp(this.mouseCursor,'on') && strcmp(this.whichMouseButton,'double click')
                this.doubleClick();
            end
            %--------------------------------------------------------------
            % Draw nodes
            if strcmp(this.drawNode,'on') && strcmp(this.whichMouseButton,'left')
                % Initialize flag for node inside element
                inWhichElems = [];
                
                % Get graphic tolerance
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
                if draw.mdl.anm.analysis_type == 1 || draw.mdl.anm.analysis_type == 3 % TRUSS
                    tol = sz/125;
                else
                    tol = sz/400;
                end
                
                % Check if snap to grid is on
                if ~isempty(this.nodeSnap) % selected an existing node
                    x = this.nodeSnap(1);
                    y = this.nodeSnap(2);
                elseif ~isempty(this.intSnap) % selected an intersection
                    x = this.intSnap(1);
                    y = this.intSnap(2);
                    intersections = getappdata(0,'intersections');
                    inWhichElems = intersections(this.whichIntSnap).elems;
                elseif ~isempty(this.elemPoint) % node is inside an elem
                    x = this.elemPoint(1);
                    y = this.elemPoint(2);
                    inWhichElems = auxModelFctn('isPointInElem',{[x y 0],tol});
                    if ~isempty(inWhichElems)
                        this.whichElemSnap = inWhichElems(end);
                        this.elemPoint = [x y];
                    end
                elseif strcmp(this.snapToGrid,'on') % get grid coordinates
                    pos = auxMouseFctn('snapToGridPosition',this);
                    x = pos(1);
                    y = pos(2);
                    inWhichElems = auxModelFctn('isPointInElem',{[x y 0],tol});
                else % clicked on blank space
                    x = this.currentPosition(1);
                    y = this.currentPosition(2);
                end
                
                % Draw nodes
                [~] = auxMouseFctn('drawNodes',this,{[x y 0],inWhichElems});
            end
            %--------------------------------------------------------------
            % Select nodes
            if strcmp(this.mouseCursor,'on') && strcmp(this.whichMouseButton,'left')
                % Clear previous selected node mark
                delete(findobj('tag', 'selectedNode'))
                delete(findobj('tag', 'selectedElem'))
                
                % Initialize flag
                deleteButtonEnabledFlag = 0;
                
                % Check if click was on a node
                if ~isempty(this.nodeSnap)
                    % Delete previous dynamic snap nodes plots
                    if ~isempty(findobj('tag', 'snapNode2'))
                        delete(findobj('tag', 'snapNode2'));
                    end
                    delete(findobj('tag', 'snapNode'));

                    % Reinitialize selectedNode property
                    this.selectedNode = 0;
                    this.moveAction();
                    
                    % Set selectedNode as the node that was clicked on
                    this.selectedNode = this.whichNodeSnap;
                    
                    % Write node info on uitables at info panel
                    [~] = auxMouseFctn('writeNodeInfoPanel',this);
                    
                    % Enable 'delete entities' button
                    set(mdata.pushbutton_DeleteEntities,'enable','on')
                    deleteButtonEnabledFlag = 1;
                    
                else % if nodeSnap is empty, user clicked on blank canvas space
                    
                    % Reset info panel uitable with model data
                    if ~isempty(this.originalData)
                        set(mdata.uitable_infoPanel,'Data',this.originalData)
                        set(mdata.uitable_infoPanelEditable,'Data',{},'enable','off')
                        set(mdata.pushbutton_ApplyInfoPanel,'enable','off')
                        this.originalData = {};
                    end
                    
                    % Update canvas and reinitialize this.selectedNode
                    if this.selectedNode ~= 0
                        this.whichNodeSnap = this.selectedNode;
                        this.selectedNode = 0;
                        this.moveAction();
                    end
                    
                    % Update canvas and reinitialize this.selectedElem
                    if this.selectedElem ~= 0
                        this.whichElemSnap = this.selectedElem;
                        elems = getappdata(0,'elems');
                        this.elemSnap = [elems(this.selectedElem).nodes(1).coord(1), elems(this.selectedElem).nodes(1).coord(2);
                                         elems(this.selectedElem).nodes(2).coord(1), elems(this.selectedElem).nodes(2).coord(2)];
                        this.selectedElem = 0;
                        this.moveAction();
                    end
                    
                    % Disable 'delete entities' button
                    set(mdata.pushbutton_DeleteEntities,'enable','off')
                end
            end
            %--------------------------------------------------------------
            % Draw elements
            if strcmp(this.drawElem,'on') && strcmp(this.whichMouseButton,'left')
                % Update element node flag (initial or final node)
                this.elemNode = this.elemNode + 1;
                
                % Check if property elemNode is 1 (initial node) or 2
                % (final node). If not, reinitialize properties and return.
                if this.elemNode ~= 1 && this.elemNode ~= 2
                    this.elemNode = 0;
                    this.elemNodeID = [];
                    this.elemCoords = [];
                    this.selectedNode = 0;
                    this.moveAction();
                    return
                end
                
                % Get canvas borders
                dfltUnits = get(gca,'units');
                set(gca,'units','normalized');
                limits = get(gca,'Position');
                set(gca,'units',dfltUnits);
                axisWidth = limits(3);
                
                % Adjust axisWidth parameter based on axis limits
                aux = [this.canvas.XLim;this.canvas.YLim];
                scl = max([diff(aux(1,:)),diff(aux(2,:))])/2;
                axisWidth = 0.75*axisWidth * 10^(floor(log10(scl)));
                
                % Draw elements and new nodes (if there are any)
                auxMouseFctn('drawElements',this,axisWidth/18);
            end
            
            %--------------------------------------------------------------
            % Select Elements
            if strcmp(this.mouseCursor,'on') && strcmp(this.whichMouseButton,'left')
                % Clear previous selected node mark
                delete(findobj('tag', 'selectedNode'))
                delete(findobj('tag', 'selectedElem'))
                
                if ~isempty(this.elemSnap)
                    % Delete previous dynamic snap elements plots
                    if ~isempty(findobj('tag', 'snapElem2'))
                        delete(findobj('tag', 'snapElem2'));
                    end
                    delete(findobj('tag', 'snapElem'));

                    % Reinitialize selectedElem property
                    this.selectedElem = 0;
                    this.moveAction();
                    
                    % Set selectedElem as the element that was clicked on
                    this.selectedElem = this.whichElemSnap;
                    
                    if strcmp(get(mdata.popupmenu_Results,'Enable'),'on') && get(mdata.popupmenu_AnalysisType,'Value') == 1
                        % Get element result values on click point
                        this.elemResults = auxModelFctn('getElemPointDisplAndStress',{this.elemPoint,this.selectedElem});
                        % Plot mark for click point
                        scatter(this.elemPoint(1),this.elemPoint(2),15,[1 0 0],'filled', 'tag','selectedElem');
                    end
                    
                    % Write selected element info on uitables at info panel
                    [~] = auxMouseFctn('writeElemInfoPanel',this);
                    
                    % Reset element results property
                    this.elemResults =[];
                    
                    % Enable 'delete entities' button
                    set(mdata.pushbutton_DeleteEntities,'enable','on')
                    
                else % if elemSnap is empty, user clicked on blank canvas space, or a node
                    
                    % Reset info panel uitable with model data
                    if ~isempty(this.originalData) && this.selectedNode == 0
                        mdata = guidata(findobj('Tag','GUI_Main'));
                        set(mdata.uitable_infoPanel,'Data',this.originalData)
                        set(mdata.uitable_infoPanelEditable,'Data',{},'enable','off')
                        set(mdata.pushbutton_ApplyInfoPanel,'enable','off','BackgroundColor',[0.94 0.94 0.94])
                        this.originalData = {};
                    end
                    
                    % Update canvas and reinitialize this.selectedElem
                    if this.selectedElem ~= 0
                        elems = getappdata(0,'elems');
                        this.whichElemSnap = this.selectedElem;
                        this.elemSnap = [elems(this.selectedElem).nodes(1).coord(1), elems(this.selectedElem).nodes(1).coord(2);
                                         elems(this.selectedElem).nodes(2).coord(1), elems(this.selectedElem).nodes(2).coord(2)];
                        this.selectedElem = 0;
                        this.moveAction();
                    end
                    
                    % Disable 'delete entities' button
                    if deleteButtonEnabledFlag == 0 % Check if delete button was enabled by selecting a node.
                        set(mdata.pushbutton_DeleteEntities,'enable','off')
                    end
                end
            end
            
            %--------------------------------------------------------------
            % Save mouse handle in root
            setappdata(0,'mouse',this)
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse pointer moves on the canvas.
        % Executes pan (limits) intermediary actions when right mouse 
        % button is still pressed.
        function moveAction(this)
            % Get handle to GUI_Main
            mdata = guidata(findobj('Tag','GUI_Main'));
            
            % Get handles to coordinate text indicators
            hToolbar = findall(groot,'tag','toolbar');
            jToolbar = hToolbar.JavaContainer.getComponentPeer;
            jaxisx   = jToolbar.getComponent(jToolbar.getComponentCount-3);
            jaxisy   = jToolbar.getComponent(jToolbar.getComponentCount-1);
            
            % Check if cursor is inside canvas limits
            x = this.currentPosition(1);
            y = this.currentPosition(2);
            xLim = get(gca,'XLim');
            yLim = get(gca,'YLim');
            inside = x >= xLim(1) && x <= xLim(2) && y >= yLim(1) && y <= yLim(2);
            
            if inside
                if strcmp(this.snapToGrid,'on')
                    % Set pointer
                    set(gcf,'Pointer','custom','PointerShapeCData',NaN(16,16));
                else
                    % Update coordinates indicator
                    jaxisx.setText(sprintf('%.3f',x));
                    jaxisy.setText(sprintf('%.3f',y));
                end
            else
                % Set pointer
                set(gcf,'Pointer','arrow');
                
                % Update coordinates indicator
                jaxisx.setText(' ');
                jaxisy.setText(' ');
            end
            
            %--------------------------------------------------------------
            % Get toolbar_pan's state ('on' or 'off')
            panState = get(mdata.toolbar_pan,'state');
            % Pan
            if strcmp(this.whichMouseButton,'center') || (strcmp(panState,'on') && strcmp(this.whichMouseButton,'left'))
                % Catches the actual window coordinates.
                wpt = get(this.dialog, 'CurrentPoint');
                wx = wpt(1); 
                wy = wpt(2);
                
                % Calculate the displacements.
                dX = this.panIniX - wx;
                dY = this.panIniY - wy;
                
                % Changes the old window coordinates.
                this.panIniX = wx;
                this.panIniY = wy;
                
                % Picks the necessary values for calculate the factors.
                originalUnits = get(this.canvas, 'Units');
                set(this.canvas, 'Units', 'pixels');
                this.axesPos = get(this.canvas, 'Position');
                set(this.canvas, 'Units', originalUnits);
                pbar = get(gca, 'PlotBoxAspectRatio');
                
                % Calculate the factor for X.
                imAspectRatioX = pbar(2) / pbar(1);
                if (imAspectRatioX ~= 1)
                    posAspectRatioX = this.axesPos(3) / this.axesPos(4);
                    arFactorX = imAspectRatioX * posAspectRatioX;
                    if (arFactorX < 1)
                        arFactorX = 1;
                    end
                else
                    arFactorX = 1;
                end
                
                % Calculate the factor for Y.
                imAspectRatioY = pbar(1) / pbar(2);
                if (imAspectRatioY ~= 1)
                    posAspectRatioY = this.axesPos(4) / this.axesPos(3);
                    arFactorY = imAspectRatioY * posAspectRatioY;
                    if (arFactorY < 1)
                        arFactorY = 1;
                    end
                else
                    arFactorY = 1;
                end

                % For log plots, transform to linear scale.
                if strcmp(get(this.canvas, 'xscale'), 'log')
                    xLim = log10(xLim);
                    xLim = FixInfLogLimits('x', xLim);
                    isXLog = true;
                else
                    isXLog = false;
                end
                dx = dX * abs(xLim(2) - xLim(1)) / (this.axesPos(3) / arFactorX);
                xLim = xLim + dx;

                % For log plots, untransform limits.
                if isXLog
                    xLim = 10.^(xLim);
                end
                
                % For log plots, transform to linear scale.
                if strcmp(get(this.canvas, 'yscale'), 'log')
                    yLim = log10(yLim);
                    yLim = FixInfLogLimits('y', yLim);
                    isYLog = true;
                else
                    isYLog = false;
                end
                dy = dY * abs(yLim(2) - yLim(1)) / (this.axesPos(4) / arFactorY);
                yLim = yLim + dy;
                
                % For log plots, untransform limits.
                if isYLog
                    yLim = 10.^(yLim);
                end
                
                % Sets the new limits.
                set(this.canvas, 'XLim', xLim);
                set(this.canvas, 'YLim', yLim);
            end
            
            %--------------------------------------------------------------
            % Snap to entities and draw dynamic lines
            
            % Check if 'rotate' or 'pan' are not being used
            if ~strcmp(this.whichMouseButton,'right') && ~strcmp(this.whichMouseButton,'center') && ~strcmp(panState,'on')
                
                % Check if selecting/modelling options are off
                if strcmp(this.mouseCursor,'off') && strcmp(this.drawNode,'off') && strcmp(this.drawElem,'off')
                    % Delete any remaining snap to nodes symbol
                    if ~isempty(findobj('tag', 'snapNode2'))
                        delete(findobj('tag', 'snapNode2'));
                    end
                    if ~isempty(findobj('tag', 'snapNode'))
                        delete(findobj('tag', 'snapNode'));
                    end
                    % Delete any remaining snap to elements symbol
                    if ~isempty(findobj('tag', 'snapElem2'))
                       delete(findobj('tag', 'snapElem2'));
                    end
                    if ~isempty(findobj('tag', 'snapElem'))
                       delete(findobj('tag', 'snapElem'));
                    end
                    % Delete any remaining snap to grid symbol
                    if ~isempty(findobj('tag','snapGrid'))
                        delete(findobj('tag','snapGrid'));
                    end
                    % Delete any remaining dynamic line
                    if ~isempty(findobj('tag','dynamicLine'))
                        delete(findobj('tag','dynamicLine'))
                    end
                    
                else % if any selecting/modelling options are on
                    
                    % Get analysis model
                    anm = get(mdata.popupmenu_Anm,'value');

                    % Get draw object
                    draw = getappdata(0,'draw');
                    if ~isempty(draw)
                        sz = draw.size;
                    else
                        sz = 5;
                    end

                    % Get canvas borders
                    dfltUnits = get(gca,'units');
                    set(gca,'units','normalized');
                    limits = get(gca,'Position');
                    set(gca,'units',dfltUnits);
                    axisWidth = limits(3);
                    
                    % Adjust axisWidth parameter based on axis limits
                    aux = [this.canvas.XLim;this.canvas.YLim];
                    scl = max([diff(aux(1,:)),diff(aux(2,:))])/2;
                    axisWidth = 0.75*axisWidth * 10^(floor(log10(scl)));
                    
                    % Set vector of snap properties
                    snapProp = [anm, sz, axisWidth];
                    
                    %---------------------------
                    % Get number of nodes
                    nnp = getappdata(0,'nnp');
                    
                    % Check if there are nodes
                    if nnp == 0
                        % ----- SNAP TO GRID ----------------------------
                        % Check which coordinates are closer to the current point,
                        % according to decimal precision, and draw temporary node, 
                        % in red.

                        % Delete previous snap to grid symbol
                        if ~isempty(findobj('tag','snapGrid'))
                            delete(findobj('tag','snapGrid'));
                        end
                        
                        % Check if snap to grid option is on
                        if strcmp(this.snapToGrid,'on')
                            % Get closest grid point to cursor
                            coords = auxMouseFctn('snapToGridPosition',this);

                            % Draw snap to grid dynamic symbol
                            [~] = auxMouseFctn('snapToGridDraw',this,{coords,snapProp});

                            % Redefine coordinates
                            x = coords(1);
                            y = coords(2);
                            
                            % Update coordinates indicator
                            if inside
                                jaxisx.setText(sprintf('%.3f',x));
                                jaxisy.setText(sprintf('%.3f',y));
                            end
                            
                            % ----- DYNAMIC LINES ----------------------------
                            % Check if an initial element node is selected, in draw
                            % elements mode.
                            if ~isempty(findobj('tag','dynamicLine'))
                                delete(findobj('tag','dynamicLine'))
                            end
                            if strcmp(this.drawElem,'on') && this.elemNode == 1
                                X = [this.elemCoords(1,1) x];
                                Y = [this.elemCoords(1,2) y];
                                plot(X,Y,'color',[0.9 0.2 0],'tag','dynamicLine')
                                hold on
                            end
                        end
                    
                    else % if there are nodes
                        
                        % ----- SNAP TO NODES ---------------------------
                        % Create an invisible square around each node, if mouse is
                        % inside a square, the related node will be drawn in red.
                        nodeFlag = auxMouseFctn('snapToNodes',this,{[x y 0],snapProp});

                        if nodeFlag == 0  % no nodes near mouse
                            this.nodeSnap = [];
                            this.whichNodeSnap = 0;
                        end
                        
                        % ----- SNAP TO INTERSECTIONS -------------------
                        intSectFlag = 0;
                        if ~isempty(findobj('tag','snapIntSect'))
                            delete(findobj('tag','snapIntSect'))
                        end
                        
                        % Get vector of handles to elemIntersection objects
                        intersections = getappdata(0,'intersections');

                        % Get number of intsects
                        nis = size(intersections,2);
                        
                        if nodeFlag == 0 && nis ~= 0 &&...
                           (strcmp(this.drawNode,'on') || strcmp(this.drawElem,'on'))
                       
                            intSectFlag = auxMouseFctn('snapToIntSects',this,{[x y 0],snapProp});
                            if intSectFlag == 0
                                this.intSnap = [];
                                this.whichIntSnap = 0;
                            end
                        else
                            this.intSnap = [];
                            this.whichIntSnap = 0;
                        end
                        
                        % ----- SNAP TO ELEMENTS ------------------------
                        % Create an invisible rectangle around each element, if
                        % mouse is inside a rectangle, the related element will be
                        % drawn in red.
                        elemFlag = 0;

                         % Delete previous snap to element point symbol
                        if ~isempty(findobj('tag','snapElemPoint'))
                            delete(findobj('tag','snapElemPoint'));
                        end

                        nel = getappdata(0,'nel');
                        if nodeFlag == 0 && intSectFlag == 0 && nel ~= 0
                           
                            elemFlag = auxMouseFctn('snapToElems',this,{[x y 0],snapProp});

                            if elemFlag == 0  % no elements near mouse
                                this.elemSnap = [];
                                this.whichElemSnap = 0;
                                this.elemPoint = [];
                            end
                        elseif (nodeFlag == 1 || intSectFlag == 1) && nel ~= 0 && this.whichElemSnap ~= 0 % a node or intersection has been snapped
                           if ~isempty(findobj('tag', 'snapElem2'))
                               delete(findobj('tag', 'snapElem2'));
                           end
                           if this.selectedElem == 0
                               delete(findobj('tag', 'snapElem'));
                           end
                           this.elemSnap = [];
                           this.whichElemSnap = 0;
                           this.elemPoint = [];
                        elseif nel == 0
                           this.elemSnap = [];
                           this.whichElemSnap = 0;
                           this.elemPoint = [];
                        end
                        
                        % ----- SNAP TO GRID ----------------------------
                        % Check which coordinates are closer to the current point,
                        % according to decimal precision, and draw temporary node, 
                        % in red.

                        % Delete previous snap to grid symbol
                        if ~isempty(findobj('tag','snapGrid'))
                            delete(findobj('tag','snapGrid'));
                        end

                        % Check if snap to grid option is on
                        if strcmp(this.snapToGrid,'on') && nodeFlag == 0 && intSectFlag == 0 && elemFlag == 0
                       
                            % Get closest grid point to cursor
                            coords = auxMouseFctn('snapToGridPosition',this);

                            % Draw snap to grid dynamic symbol
                            [~] = auxMouseFctn('snapToGridDraw',this,{coords,snapProp});

                            % Redefine coordinates
                            x = coords(1);
                            y = coords(2);
                            
                            % Update coordinates indicator
                            if inside
                                jaxisx.setText(sprintf('%.3f',x));
                                jaxisy.setText(sprintf('%.3f',y));
                            end
                        end

                        % ----- DYNAMIC LINES ----------------------------
                        % Check if an initial element node is selected, in draw
                        % elements mode.
                        if ~isempty(findobj('tag','dynamicLine'))
                            delete(findobj('tag','dynamicLine'))
                        end
                        if strcmp(this.drawElem,'on') && this.elemNode == 1
                            if strcmp(this.ortho,'on')
                                this.orthoCoords = auxMouseFctn('orthoPosition',this);
                                x = this.orthoCoords(1);
                                y = this.orthoCoords(2);
                            else
                                this.orthoCoords = [];
                            end
                            X = [this.elemCoords(1,1) x];
                            Y = [this.elemCoords(1,2) y];
                            plot(X,Y,'color',[0.9 0.2 0],'tag','dynamicLine')
                            hold on
                        end
                    end
                end
            end
            
            %--------------------------------------------------------------
            % Save mouse handle in root
            setappdata(0,'mouse',this)
        end
        
        %------------------------------------------------------------------
        % Action executed when an mouse button is unpressed.
        function upAction(~)
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse scroll is used.
        % Executes zoom (limits) using the mouse scroll.
        function scrollAction(this,direction,cx,cy)

            % Get the axes limits.
            xLim = get(this.canvas, 'XLim');
            yLim = get(this.canvas, 'YLim');
            
            % Calculate the current zoom.
            currentZoomX = (abs(diff([this.originalXLim(1) this.originalXLim(2)])))...
                * 100 / (abs(diff([xLim(1) xLim(2)])));
            currentZoomY = (abs(diff([this.originalYLim(1) this.originalYLim(2)])))...
                * 100 / (abs(diff([yLim(1) yLim(2)])));
            
            % Obtain the zoom percent in axis x and y.
            zoomPercent_x = this.zoomGrid(this.zoomIndiceX);
            zoomPercent_y = this.zoomGrid(this.zoomIndiceY);
            
            % Obtain the zoom indice for x and y.
            if (currentZoomX ~= zoomPercent_x)
                [~,this.zoomIndiceX] = min(abs(this.zoomGrid - currentZoomX));
            end
            
            if (currentZoomY ~= zoomPercent_y)
                [~,this.zoomIndiceY] = min(abs(this.zoomGrid - currentZoomY));
            end
            
            % Increment or decrement the zoom indice for x.
            switch direction
                case 'plus'
                    if this.zoomIndiceX < this.zoomSteps
                        this.zoomIndiceX = this.zoomIndiceX + 1;
                    end
                case 'minus'
                    if this.zoomIndiceX > 1
                        this.zoomIndiceX = this.zoomIndiceX - 1;
                    end
            end
            
            % Update the zoom indice for x.
            zoomPercent_x = this.zoomGrid(this.zoomIndiceX);

            % Calculate the new limits for axis x.
            if(cx < xLim(1))
                cx_aux = xLim(1);
            elseif(cx > xLim(2))
                cx_aux = xLim(2);
            else
                cx_aux = cx;
            end
            rf_x = abs(diff([xLim(1) xLim(2)]));
            ra_x = abs(diff([xLim(1) cx_aux]));
            rb_x = abs(diff([cx_aux xLim(2)]));
            cfa_x = ra_x / rf_x;
            cfb_x = rb_x / rf_x;
            newRange_x = abs(diff([this.originalXLim(1) this.originalXLim(2)])) * 100 / zoomPercent_x;
            dRange_x = newRange_x - rf_x;
            xLim(1) = xLim(1) - dRange_x * cfa_x;
            xLim(2) = xLim(2) + dRange_x * cfb_x;
            
            % Increment or decrement the zoom indice for y.
            switch direction
                case 'plus'
                    if this.zoomIndiceY < this.zoomSteps
                        this.zoomIndiceY = this.zoomIndiceY + 1;
                    end
                case 'minus'
                    if this.zoomIndiceY > 1
                        this.zoomIndiceY = this.zoomIndiceY - 1;
                    end
            end

            % Update the zoom indice for y.
            zoomPercent_y = this.zoomGrid(this.zoomIndiceY);

            % Calculate the new limits for axis y.
            if (cy < yLim(1))
                cy_aux = yLim(1);
            elseif (cy > yLim(2))
                cy_aux = yLim(2);
            else
                cy_aux = cy;
            end
            rf_y = abs(diff([yLim(1) yLim(2)]));
            ra_y = abs(diff([yLim(1) cy_aux]));
            rb_y = abs(diff([cy_aux yLim(2)]));
            cfa_y = ra_y / rf_y;
            cfb_y = rb_y / rf_y;
            newRange_y = abs(diff([this.originalYLim(1) this.originalYLim(2)])) * 100 / zoomPercent_y;
            dRange_y = newRange_y - rf_y;
            yLim(1) = yLim(1) - dRange_y * cfa_y;
            yLim(2) = yLim(2) + dRange_y * cfb_y;
            
            % Sets the new limits
            set(this.canvas, 'XLim', xLim);
            set(this.canvas, 'YLim', yLim);
            
            % Set axes position
            dfltUnits = get(this.canvas, 'Units');
            set(this.canvas, 'Units', 'normalized');
            set(this.canvas, 'Position', [0.207,0.049,0.776,0.926]);
            set(this.canvas, 'Units', dfltUnits);
            
            % Save mouse handle in root
            setappdata(0,'mouse',this)
        end
        
        %------------------------------------------------------------------
        % Action executed when double click is used.
        % Refresh the model (fit to view and redraw).
        function doubleClick(this)
            draw = getappdata(0,'draw');
            axis equal;
            draw.setLimits(diff(this.canvas.XLim)/diff(this.canvas.YLim));
            axis equal;
            setappdata(0,'draw',draw);
            this.originalXLim = get(this.canvas,'XLim');
            this.originalYLim = get(this.canvas,'YLim');
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
