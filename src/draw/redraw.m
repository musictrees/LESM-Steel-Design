%% Redraw function
%
% This file contains a function that works as a switch to call other
% functions. Redraws requested objects of the model.
%
%% ------------------------------------------------------------------------
% Draw model based on the selected option in popupmenu_Results
% INPUT: 
% * mdata -> handle to GUI_Main
% * toBeRedrawn -> string identifier for what needs to be redrawn
function redraw(mdata,toBeRedrawn,fctnArgIn)
if nargin == 1
    resVal = get(mdata.popupmenu_Results,'Value');
    resStr = get(mdata.popupmenu_Results,'String');
    result = char(resStr(resVal,:));
    toBeRedrawn = result;
    fctnArgIn = [];
elseif nargin == 2
    fctnArgIn = [];
end

% Make sure that modifications will be made on canvas
axes(mdata.axes_Canvas);

switch toBeRedrawn
    case 'Model'
        drawModel(mdata);
    case 'Deformation'
        drawDeformation(mdata);
    case 'Axial Force'
        drawAxialForce(mdata);
    case 'Torsion Moment'
        drawTorsionMoment(mdata);
    case 'Shear Force Y'
        drawShearForceY(mdata);
    case 'Shear Force Z'
        drawShearForceZ(mdata);
    case 'Bending Moment Y'
        drawBendingMomentY(mdata);
    case 'Bending Moment Z'
        drawBendingMomentZ(mdata);
    case 'Vibration Modes'
        drawVibration(mdata);
    case 'Motion'
        drawDynamicDeform(mdata);
    case 'Nodes'
        if isempty(fctnArgIn)
            drawNodes(mdata);
        else
            drawNodes(mdata,fctnArgIn);
        end
    case 'Elements'
        drawElems(mdata);
    case 'Loads'
        drawLoads(mdata);
    case 'Nodal Loads'
        drawNodalLoads(mdata);
    case 'Element Loads'
        drawElemLoads(mdata);
    case 'Units'
        textUnits(mdata);
    case 'Decimal Precision'
        textDecPrec(mdata);
end

% Update mouse draw size flag
draw = getappdata(0,'draw');
mouse = getappdata(0,'mouse');
mouse.sizeFlag = draw.size;
setappdata(0,'mouse',mouse)
end

%% ------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws model on canvas
function drawModel(mdata)
% Get objects from root
model = getappdata(0,'model');
mouse = getappdata(0,'mouse');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get analysis model
anm = get(mdata.popupmenu_Anm,'value');
anl = get(mdata.popupmenu_AnalysisType,'Value');

% Clean canvas
axes(mdata.axes_Canvas);
cla reset

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw model
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
if getappdata(0,'vis') == 1
    if strcmp(get(mdata.nodeIDButton,'Checked'),'on')
        draw.nodeID();
    end
    if strcmp(get(mdata.elemIDButton,'Checked'),'on')
        draw.elementID();
    end
    if strcmp(get(mdata.orientationButton,'Checked'),'on')
        draw.elementOrientation();
    end
end
if get(mdata.checkbox_Reactions,'Value') == 1
    draw.reactions();
else
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

% Adjust axes
ax = mdata.axes_Canvas;
axis equal;
if anm == 1 || anm == 2 || (anm == 3 && strcmp(get(mdata.togglebutton_2DView,'state'),'on'))
    draw.setLimits(diff(mdata.axes_Canvas.XLim)/diff(mdata.axes_Canvas.YLim));
    ax.Clipping = 'on';
    axis equal;
else
    draw.setLimits();
    ax.Clipping = 'off';
end

% Set plane of visualization
if anm == 1 || anm == 2
    view(2);
    set(mdata.planeButton,'Enable','off');
elseif anm == 3 && strcmp(get(mdata.togglebutton_2DView,'state'),'on')
    view(2);
    set(mdata.planeButton,'Enable','on');
    set(mdata.plane3D,'Checked','off');
    set(mdata.planeXY,'Checked','on');
    set(mdata.planeXZ,'Checked','off','enable','off');
    set(mdata.planeYZ,'Checked','off','enable','off');
else
    view(3);
    set(mdata.planeButton,'Enable','on');
    set(mdata.plane3D,'Checked','on');
    set(mdata.planeXY,'Checked','off');
    set(mdata.planeXZ,'Checked','off');
    set(mdata.planeYZ,'Checked','off');
    
    % Set axes labels positions
    if anm == 4 || anm == 5
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
    end
end

% Turn grid on/off
if strcmp(get(mdata.gridButton,'Checked'),'on')
    grid on;
else
    grid off;
end

% Turn ruler on/off
if strcmp(get(mdata.rulerButton,'Checked'),'on')
    mdata.axes_Canvas.XAxis.Visible = 'on';
    mdata.axes_Canvas.YAxis.Visible = 'on';
    mdata.axes_Canvas.ZAxis.Visible = 'on';
else
    mdata.axes_Canvas.XAxis.Visible = 'off';
    mdata.axes_Canvas.YAxis.Visible = 'off';
    mdata.axes_Canvas.ZAxis.Visible = 'off';
end

% Reinitialize object for mouse events and save it in root
mouse.originalXLim = get(mdata.axes_Canvas,'XLim');
mouse.originalYLim = get(mdata.axes_Canvas,'YLim');
if anm == 4 || anm == 5 || (anm == 3 && strcmp(get(mdata.togglebutton_2DView,'state'),'off'))
    mouse.originalZLim = get(mdata.axes_Canvas,'ZLim');
end
mouse.sizeFlag = draw.size;
setappdata(0,'mouse',mouse);
setappdata(0,'model',model);
setappdata(0,'draw',draw);

% Ensure fit to view
mouse.doubleClick();
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws deformation on canvas
function drawDeformation(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Clean model
clearLoadAndResults()

% Enable scale indicator
set(mdata.edit_Scale,'Visible','on','Enable','on');

% Calculate scale value
dsf = getappdata(0,'deform_sf');
if isnan(str2double(get(mdata.edit_Scale,'String')))
    scale = dsf * get(mdata.slider_Scale,'Value');
    set(mdata.edit_Scale,'String',num2str(scale,3))
else    
    scale = str2double(get(mdata.edit_Scale,'String'));    
end

% Draw deformed configuration
draw.deformConfig(scale);
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws Axial force diagram on canvas
function drawAxialForce(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Calculate scale value
asf = getappdata(0,'axial_sf');
scale = asf * get(mdata.slider_Scale,'Value');

% Draw axial force diagram
switch analysis
    case 'Static'
        draw.axialForce(scale);
    case 'Dynamic'
        draw.axialForceEnvelop(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws torsion moment diagram on canvas
function drawTorsionMoment(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw torsion moment diagram
switch analysis
    case 'Static'
        draw.torsionMoment();
    case 'Dynamic'
        % Calculate scale value
        tsf = getappdata(0,'torsion_sf');
        scale = tsf * get(mdata.slider_Scale,'Value');

        draw.torsionMomentEnvelop(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws shear force Y diagram on canvas
function drawShearForceY(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Calculate scale value
ssf = getappdata(0,'shearXY_sf');
scale = ssf * get(mdata.slider_Scale,'Value');

% Draw shear force diagram
switch analysis
    case 'Static'
        draw.shearForce_XY(scale);
    case 'Dynamic'
        draw.shearForceEnvelop_XY(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws shear force Z diagram on canvas
function drawShearForceZ(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Calculate scale value
ssf = getappdata(0,'shearXZ_sf');
scale = ssf * get(mdata.slider_Scale,'Value');

% Draw shear force diagram
switch analysis
    case 'Static'
        draw.shearForce_XZ(scale);
    case 'Dynamic'
        draw.shearForceEnvelop_XZ(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws bending moment Y diagram on canvas
function drawBendingMomentY(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Calculate scale value
bsf = getappdata(0,'bendingXZ_sf');
scale = bsf * get(mdata.slider_Scale,'Value');

% Draw bending moment diagram
switch analysis
    case 'Static'
        draw.bendingMoment_XZ(scale);
    case 'Dynamic'
       draw.bendingMomentEnvelop_XZ(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws bending moment Z diagram on canvas
function drawBendingMomentZ(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get type of analysis
anlVal   = get(mdata.popupmenu_AnalysisType,'Value');
anlStr   = get(mdata.popupmenu_AnalysisType,'String');
analysis = char(anlStr(anlVal,:));

% Clean model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Calculate scale value
bsf = getappdata(0,'bendingXY_sf');
scale = bsf * get(mdata.slider_Scale,'Value');

% Draw bending moment diagram
switch analysis
    case 'Static'
        draw.bendingMoment_XY(scale);
    case 'Dynamic'
       draw.bendingMomentEnvelop_XY(scale);
end
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws requested natural vibration mode on canvas
function drawVibration(mdata)
    clearResults();
    
    % Clean model
    clearLoadAndResults();
    
    % Get draw object
    draw = getappdata(0,'draw');
    model = getappdata(0,'model');
    draw.mdl = model;
    
    % Get mode identifier
    nMode = get(mdata.popupmenu_ElementResults,'Value');
        
    % Check if 'play' button is on
    if ~get(mdata.pushbutton_PlayDynamicResults,'UserData')
        % Delete any existing previous frame of dynamic oscillation
        delete(findobj('tag', 'drawVibrationMode'))

        % Update flag in root
        setappdata(0,'isDrawingDynamic',false)
        
        % Draw single frame of dynamic oscillation
        scl = 0.5 * draw.size * get(mdata.slider_Scale,'value') / get(mdata.slider_Scale,'max');
        draw.vibrationMode(nMode,scl);
    else
        % Get scale from slider
        scl = 0.5 * draw.size * get(mdata.slider_Scale,'value') / get(mdata.slider_Scale,'max');
        speed = mdata.dynResSpeed;
        if (speed<=0), speed=4; end % defaults to 4 when realtime (by choice)
        scale = linspace(scl,-scl,60/speed);

        % Initialize counters for scale vector
        counter = 1;
        incr = 1;

        % Initialize stop flag for while operation
        stop = false;

        % Set flag to root
        setappdata(0,'isDrawingDynamic',true)

        % Get initial time
        c = clock;
        % Oscillations are hardcoded to last a max of 20s
        if c(6) >= 40  % 60 - 20 = 40
            addTime = 60;
        else
            addTime = 0;
        end
        ti = c(6);

        % Draw vibration
        while ~stop
            % Delete previous frame of dynamic oscillation
            delete(findobj('tag', 'drawVibrationMode'))

            % Draw next frame of dynamic oscillation
            draw.vibrationMode(nMode,scale(counter));

            % Increment counter according to possible scale values
            if counter == length(scale) && length(scale) ~= 1
                incr = -1;
            elseif counter == length(scale) && length(scale) == 1
                incr = 0;
            elseif counter == 1
                incr = 1;
            end
            counter = counter + incr;

            % Get time after this step
            c = clock;
            if c(6) < ti
                tf = c(6) + addTime;
            else
                tf = c(6);
            end
            
            % Check if 20s has past, or flag in root signalized to stop
            if tf-ti > 20 || ~getappdata(0,'isDrawingDynamic')
                stop = true;
            end
        end
        % Update flag in root
        setappdata(0,'isDrawingDynamic',false)
        
        % Update play button properties
        set(mdata.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])
    end
    % Return objects to root
    setappdata(0,'draw',draw)
    setappdata(0,'model',model)
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws dynamic deformation on canvas
function drawDynamicDeform(mdata)
    clearResults();
    
    % Clean model
    clearLoadAndResults();
    
    % Adjust axes
    anm = get(mdata.popupmenu_Anm,'value');
    if anm == 1 || anm == 2
        ax = mdata.axes_Canvas;
        ax.Clipping = 'on';
    else
        ax = mdata.axes_Canvas;
        ax.Clipping = 'off';
    end

    % Get draw object
    draw = getappdata(0,'draw');
    model  = getappdata(0,'model');
    draw.mdl = model;
    
    % Calculate scale value
    dsf = getappdata(0,'deform_sf');
    scale = dsf * round(get(mdata.slider_Scale,'value')/2,0);
    if scale == 0
        scale = 1;
    end
    
    % Check if 'play' button is on
    if ~get(mdata.pushbutton_PlayDynamicResults,'UserData')
        % Delete any existing previous frame of dynamic oscillation
        delete(findobj('tag', 'drawDynamicDeform'))

        % Update flag in root
        setappdata(0,'isDrawingDynamic',false)
        
        % Draw single frame of dynamic oscillation
        draw.dynamicDeform(1,scale);
    else
        % Initialize stop flag for while operation
        stop = false;
        
        % Initialize counter
        counter = 1;

        % Set flag to root
        setappdata(0,'isDrawingDynamic',true)

        if mdata.dynResSpeed > 0 % speed multipliers
            % Draw vibration
            while ~stop
                % Delete previous frame of dynamic oscillation
                delete(findobj('tag', 'drawDynamicDeform'))

                % Get next scale factor and draw next frame of dynamic oscillation
                draw.dynamicDeform(counter,scale);

                % Increment counter
                counter = counter + mdata.dynResSpeed;

                % Check if has reached final step or flag in root signalized to stop
                if counter > (model.n_steps+1) || ~getappdata(0,'isDrawingDynamic')
                    stop = true;
                end
            end
        else % real time
            
            step = model.t/model.n_steps;
            c  =  clock;
            t0 = c(4:6);
            dt0 = t0(1)*3600 + t0(2)*60 + t0(3);
            
            % Draw vibration
            while ~stop
                % Delete previous frame of dynamic oscillation
                delete(findobj('tag', 'drawDynamicDeform'))

                % Get next scale factor and draw next frame of dynamic oscillation
                draw.dynamicDeform(counter,scale);

                % Increment counter
                c  =  clock;
                tt = c(4:6);
                if tt(1)<t0(1)
                    tt(1) = tt(1) + 24;
                end
                dt = tt(1)*3600 + tt(2)*60 + tt(3);
                counter = (dt-dt0)/step + counter;
                t0  = tt;
                dt0 = dt;

                % Check if has reached final step or flag in root signalized to stop
                if counter > (model.n_steps+1) || ~getappdata(0,'isDrawingDynamic')
                    stop = true;
                end
            end
        end
        
        % Update flag in root
        setappdata(0,'isDrawingDynamic',false)
        
        % Update play button properties
        set(mdata.pushbutton_PlayDynamicResults,'UserData',false,'ForegroundColor',[0,0.7,0.2])
    end
    % Return objects to root
    setappdata(0,'draw',draw)
    setappdata(0,'model',model)
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws nodes and supports on canvas
function drawNodes(mdata,setLimitsFlag)
% Flag for need to reset axis limits
if nargin == 1
    setLimitsFlag = true;
end

% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Get 2D axis propotions
ax = mdata.axes_Canvas;
axProp_XY = diff(ax.XLim)/diff(ax.YLim);

% Clean canvas
clearNodesSpprts()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw nodes
draw.setSize();
draw.nodes();
if getappdata(0,'vis') == 1
    if (strcmp(get(mdata.nodeIDButton,'Checked'),'on') == 1) % Check if nodes ID is on
        draw.nodeID();
    end
    if (strcmp(get(mdata.gridButton,'Checked'),'on') == 1) % Check if grid is on
        grid on
    else
        grid off
    end
    if model.nnp == 0 || model.nnp == 1
        switch get(mdata.rulerButton,'Checked')
            case 'off'
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                ax.ZAxis.Visible = 'off';
            case 'on'
                ax.XAxis.Visible = 'on';
                ax.YAxis.Visible = 'on';
                ax.ZAxis.Visible = 'on';
        end
    end
end
if model.anm.analysis_type <= 1 % TRUSS_2D / FRAME_2D
    if setLimitsFlag == true
        draw.setLimits(axProp_XY);
    end
elseif model.anm.analysis_type == 2 && strcmp(get(mdata.togglebutton_2DView,'state'),'on')  % GRILLAGE_2D
    if setLimitsFlag == true
        draw.setLimits(axProp_XY);
    end
elseif model.anm.analysis_type == 2 % GRILLAGE_3D
    if setLimitsFlag == true
        axis equal
        draw.setLimits();
    end
else % TRUSS_3D / FRAME_3D
    if setLimitsFlag == true
        axis equal
        draw.setLimits();
    end
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

% Reinitialize object for mouse events and save it in root
mouse = getappdata(0,'mouse');
if setLimitsFlag == true
    mouse.originalXLim = get(mdata.axes_Canvas,'XLim');
    mouse.originalYLim = get(mdata.axes_Canvas,'YLim');
    anm = get(mdata.popupmenu_Anm,'value');
    if anm == 3 && strcmp(get(mdata.togglebutton_2DView,'state'),'off')
        mouse.originalZLim = get(mdata.axes_Canvas,'ZLim');
    elseif (anm == 4) || (anm == 5)
        mouse.originalZLim = get(mdata.axes_Canvas,'ZLim');
    end
end
mouse.sizeFlag = draw.size;
setappdata(0,'mouse',mouse);
setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws elements on canvas
function drawElems(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Clean canvas
clearElems()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw model
draw.elements();
if getappdata(0,'vis') == 1
    if (strcmp(get(mdata.elemIDButton,'Checked'),'on') == 1) % Check if elements ID is on
        draw.elementID();
    end
    if (strcmp(get(mdata.orientationButton,'Checked'),'on') == 1) % Check if elements orientation is on
        draw.elementOrientation();
    end
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws loads on canvas
function drawLoads(mdata)
include_constants;

% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Clear model
clearLoadAndResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw model
if model.whichSolver == STATIC_LINEAR
    draw.elemLoadsScaleFactor();
    draw.nodalLoads();
    draw.elemLoads();
    draw.thermalLoads();
    draw.nodalPrescDispl();
else
    draw.dynamicNodalLoads();
    draw.nodalMass();
    draw.nodalInitialConditions();
end

anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws loads on canvas
function drawNodalLoads(mdata)
include_constants;

% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Clear nodal loads
delete(findobj('tag','drawNodalLoads'))
delete(findobj('tag','drawNodalMass'))
delete(findobj('tag','textNodalLoads'))
delete(findobj('tag','textNodalMoments'))
delete(findobj('tag','textNodalMass'))
delete(findobj('tag','drawPrescDispl'))
delete(findobj('tag','textPrescDispl'))
delete(findobj('tag','textPrescRot'))
delete(findobj('tag','textInitialConditions'))
clearResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw model
if model.whichSolver == STATIC_LINEAR
    draw.nodalLoads();
    draw.nodalPrescDispl();
else
    draw.dynamicNodalLoads();
    draw.nodalMass();
    draw.nodalInitialConditions();
end

anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Draws loads on canvas
function drawElemLoads(mdata)
% Get objects from root
model = getappdata(0,'model');
draw = getappdata(0,'draw');
draw.mdl = model;

% Clear elem loads
delete(findobj('tag','drawElemLoads'))
delete(findobj('tag','textElemLoads'))
delete(findobj('tag','drawThermalLoads'))
delete(findobj('tag','textThermalLoads'))
clearResults()

% Disable scale indicator
set(mdata.edit_Scale,'Visible','off','Enable','off');

% Draw model
draw.elemLoadsScaleFactor();
draw.elemLoads();
draw.thermalLoads();
anm = get(mdata.popupmenu_Anm,'value');
if anm == 1 || anm == 2
    ax = mdata.axes_Canvas;
    ax.Clipping = 'on';
else
    ax = mdata.axes_Canvas;
    ax.Clipping = 'off';
end

% Make sure not to display reactions if checkbox is unmarked
if get(mdata.checkbox_Reactions,'Value') == 0
    delete(findobj('tag','drawReactions'));
    delete(findobj('tag','textForceReactions'));
    delete(findobj('tag','textMomentReactions'));
end

setappdata(0,'model',model);
setappdata(0,'draw',draw);
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Writes/Deletes units from texts on canvas
function textUnits(mdata)
% Get units button state - on/off
check = get(mdata.unitsButton,'Checked');

% Get texts
springStiff     = findobj('tag','textSprings');
rotSpringStiff  = findobj('tag','textRotSprings');
srjoints        = findobj('tag','textSemiRigid');
nodalLoads      = findobj('tag','textNodalLoads');
nodalMoments    = findobj('tag','textNodalMoments');
nodalMass       = findobj('tag','textNodalMass');
prescDispl      = findobj('tag','textPrescDispl');
prescRot        = findobj('tag','textPrescRot');
elemLoads       = findobj('tag','textElemLoads');
thermalLoads    = findobj('tag','textThermalLoads');
axialForce      = findobj('tag','textAxialForceDiagram');
torsionMoment   = findobj('tag','textTorsionDiagram');
shearForceXY    = findobj('tag','textShearForceXYDiagram');
shearForceXZ    = findobj('tag','textShearForceXZDiagram');
bendMomentXY    = findobj('tag','textBendMomentXYDiagram');
bendMomentXZ    = findobj('tag','textBendMomentXZDiagram');
reactions       = findobj('tag','textForceReactions');
momentReactions = findobj('tag','textMomentReactions');

% Update texts accordingly to current unitsButton 'Checked' state - on/off
switch check
    case 'on' % Units button has been turned on
        for i = 1:size(springStiff,1)
            springStiff(i).String = strcat(springStiff(i).String, ' kN/m');
        end
        for i = 1:size(rotSpringStiff,1)
            rotSpringStiff(i).String = strcat(rotSpringStiff(i).String, ' kNm/rad');
        end
        for i = 1:size(srjoints,1)
            srjoints(i).String = strcat(srjoints(i).String, ' kNm/rad');
        end
        for i = 1:size(nodalLoads,1)
            nodalLoads(i).String = strcat(nodalLoads(i).String, ' kN');
        end
        for i = 1:size(nodalMass,1)
            nodalMass(i).String = strcat(nodalMass(i).String, ' kg');
        end
        for i = 1:size(nodalMoments,1)
            nodalMoments(i).String = strcat(nodalMoments(i).String, ' kNm');
        end
        for i = 1:size(prescDispl,1)
            prescDispl(i).String = strcat(prescDispl(i).String, ' mm');
        end
        for i = 1:size(prescRot,1)
            prescRot(i).String = strcat(prescRot(i).String, ' rad');
        end
        for i = 1:size(elemLoads,1)
            elemLoads(i).String = strcat(elemLoads(i).String, ' kN/m');
        end
        for i = 1:size(thermalLoads,1)
            thermalLoads(i).String = strcat(thermalLoads(i).String, ' oC');
        end
        for i = 1:size(axialForce,1)
            axialForce(i).String = strcat(axialForce(i).String, ' kN');
        end
        for i = 1:size(torsionMoment,1)
            torsionMoment(i).String = strcat(torsionMoment(i).String, ' kNm');
        end
        for i = 1:size(shearForceXY,1)
            shearForceXY(i).String = strcat(shearForceXY(i).String, ' kN');
        end
        for i = 1:size(shearForceXZ,1)
            shearForceXZ(i).String = strcat(shearForceXZ(i).String, ' kN');
        end
        for i = 1:size(bendMomentXY,1)
            bendMomentXY(i).String = strcat(bendMomentXY(i).String, ' kNm');
        end
        for i = 1:size(bendMomentXZ,1)
            bendMomentXZ(i).String = strcat(bendMomentXZ(i).String, ' kNm');
        end
        for i = 1:size(reactions,1)
            reactions(i).String = strcat(reactions(i).String, ' kN');
        end
        for i = 1:size(momentReactions,1)
            momentReactions(i).String = strcat(momentReactions(i).String, ' kNm');
        end
        
    case 'off' % Units button has been turned off
        for i = 1:size(springStiff,1)
            springStiff(i).String = springStiff(i).String(1:end-5);
        end
        for i = 1:size(rotSpringStiff,1)
            rotSpringStiff(i).String = rotSpringStiff(i).String(1:end-8);
        end
        for i = 1:size(srjoints,1)
            srjoints(i).String = srjoints(i).String(1:end-8);
        end
        for i = 1:size(nodalLoads,1)
            nodalLoads(i).String = nodalLoads(i).String(1:end-3);
        end
        for i = 1:size(nodalMass,1)
            nodalMass(i).String = nodalMass(i).String(1:end-3);
        end
        for i = 1:size(nodalMoments,1)
            nodalMoments(i).String = nodalMoments(i).String(1:end-4);
        end
        for i = 1:size(prescDispl,1)
            prescDispl(i).String = prescDispl(i).String(1:end-3);
        end
        for i = 1:size(prescRot,1)
            prescRot(i).String = prescRot(i).String(1:end-4);
        end
        for i = 1:size(elemLoads,1)
            elemLoads(i).String = elemLoads(i).String(1:end-5);
        end
        for i = 1:size(thermalLoads,1)
            thermalLoads(i).String = thermalLoads(i).String(1:end-3);
        end
        for i = 1:size(axialForce,1)
            axialForce(i).String = axialForce(i).String(1:end-3);
        end
        for i = 1:size(torsionMoment,1)
            torsionMoment(i).String = torsionMoment(i).String(1:end-4);
        end
        for i = 1:size(shearForceXY,1)
            shearForceXY(i).String = shearForceXY(i).String(1:end-3);
        end
        for i = 1:size(shearForceXZ,1)
            shearForceXZ(i).String = shearForceXZ(i).String(1:end-3);
        end
        for i = 1:size(bendMomentXY,1)
            bendMomentXY(i).String = bendMomentXY(i).String(1:end-4);
        end
        for i = 1:size(bendMomentXZ,1)
            bendMomentXZ(i).String = bendMomentXZ(i).String(1:end-4);
        end
        for i = 1:size(reactions,1)
            reactions(i).String = reactions(i).String(1:end-3);
        end
        for i = 1:size(momentReactions,1)
            momentReactions(i).String = momentReactions(i).String(1:end-4);
        end
end
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Resets all decimal precision of texts on canvas
function textDecPrec(mdata)
% Get decimal precision
dp = getappdata(0,'decPrec');

% Get texts
springStiff     = findobj('tag','textSprings');
rotSpringStiff  = findobj('tag','textRotSprings');
srjoints        = findobj('tag','textSemiRigid');
nodalLoads      = findobj('tag','textNodalLoads');
nodalMoments    = findobj('tag','textNodalMoments');
nodalMass       = findobj('tag','textNodalMass');
prescDispl      = findobj('tag','textPrescDispl');
prescRot        = findobj('tag','textPrescRot');
elemLoads       = findobj('tag','textElemLoads');
thermalLoads    = findobj('tag','textThermalLoads');
axialForce      = findobj('tag','textAxialForceDiagram');
torsionMoment   = findobj('tag','textTorsionDiagram');
shearForceXY    = findobj('tag','textShearForceXYDiagram');
shearForceXZ    = findobj('tag','textShearForceXZDiagram');
bendMomentXY    = findobj('tag','textBendMomentXYDiagram');
bendMomentXZ    = findobj('tag','textBendMomentXZDiagram');
reactions       = findobj('tag','textForceReactions');
momentReactions = findobj('tag','textMomentReactions');

% Allocate all texts in vector of handles to text objects
springTexts = vertcat(springStiff,rotSpringStiff,srjoints);
allTexts = vertcat(nodalLoads,nodalMoments,nodalMass,prescDispl,prescRot,elemLoads,...
                   axialForce,torsionMoment,shearForceXY,shearForceXZ,...
                   bendMomentXY,bendMomentXZ,reactions,momentReactions);

% Reset spring texts
for i = 1:size(springTexts,1)
    if springTexts(i).UserData >= 1000
        %springTexts(i).String = sprintf('%+.*e',dp,springTexts(i).UserData);
        springTexts(i).String = sprintf('%.*e',dp,springTexts(i).UserData);
    else
        %springTexts(i).String = sprintf('%+.*f',dp,springTexts(i).UserData);
        springTexts(i).String = sprintf('%.*f',dp,springTexts(i).UserData);
    end
end

% Reset thermal load texts
for i = 1:size(thermalLoads,1)
    str = thermalLoads(i).UserData{1};
    dt = thermalLoads(i).UserData{2};
    thermalLoads(i).String = sprintf('%s%+.*f',str,dp,dt);
end

% Reset all other texts
for i = 1:size(allTexts,1)
    %allTexts(i).String = sprintf('%+.*f',dp,allTexts(i).UserData);
    allTexts(i).String = sprintf('%.*f',dp,allTexts(i).UserData);
end

% Set units to texts, if necessary
if strcmp(get(mdata.unitsButton,'Checked'),'on') == true
    textUnits(mdata);
end
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Clears nodes and supports
function clearNodesSpprts()
    delete(findobj('tag','drawNodes'))
    delete(findobj('tag','drawSupports'))
    delete(findobj('tag','textSprings'))
    delete(findobj('tag','textRotSprings'))
    delete(findobj('tag','textNodeID'))
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Clears elements
function clearElems()
    delete(findobj('tag','drawElements'))
    delete(findobj('tag','drawSemiRigid'))
    delete(findobj('tag','drawSemiRigidTemp'))
    delete(findobj('tag','textSemiRigid'))
    delete(findobj('tag','textElemID'))
    delete(findobj('tag','drawElemOrient'))
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Clears loads/deformed configuration/diagrams
function clearLoadAndResults()
    delete(findobj('tag','drawNodalLoads'))
    delete(findobj('tag','drawNodalMass'))
    delete(findobj('tag','textNodalLoads'))
    delete(findobj('tag','textNodalMoments'))
    delete(findobj('tag','textNodalMass'))
    delete(findobj('tag','drawPrescDispl'))
    delete(findobj('tag','textPrescDispl'))
    delete(findobj('tag','textPrescRot'))
    delete(findobj('tag','textInitialConditions'))
    delete(findobj('tag','drawElemLoads'))
    delete(findobj('tag','textElemLoads'))
    delete(findobj('tag','drawThermalLoads'))
    delete(findobj('tag','textThermalLoads'))
    clearResults()
end

%--------------------------------------------------------------------------
% Auxiliary function to the redraw function
% Clears deformed configuration/diagrams
function clearResults()
    delete(findobj('tag','drawVibrationMode'))
    delete(findobj('tag','drawDynamicDeform'))
    delete(findobj('tag','drawDeformConfig'))
    delete(findobj('tag','drawAxialForceDiagram'))
    delete(findobj('tag','textAxialForceDiagram'))
    delete(findobj('tag','drawTorsionDiagram'))
    delete(findobj('tag','textTorsionDiagram'))
    delete(findobj('tag','drawShearForceXYDiagram'))
    delete(findobj('tag','textShearForceXYDiagram'))
    delete(findobj('tag','drawShearForceXZDiagram'))
    delete(findobj('tag','textShearForceXZDiagram'))
    delete(findobj('tag','drawBendMomentXYDiagram'))
    delete(findobj('tag','textBendMomentXYDiagram'))
    delete(findobj('tag','drawBendMomentXZDiagram'))
    delete(findobj('tag','textBendMomentXZDiagram'))
end
