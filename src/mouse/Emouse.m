%% Emouse class
%
%% Description
%
% This is an abstract class to facilitate the development of applications
% that handle mouse events on canvas (axes: the drawing area of a GUI
% application in MATLAB).
%
% The abstract *Emouse* class presents, in addition to the constructor
% method, 4 private concrete methods (implemented) and 4 abstract methods
% that must be implemented by the client user. Its use is achieved by
% creating a client subclass that inherits its properties and implements
% the 4 abstract methods:
%%%
% * *downAction*: This method must be implemented with the procedures to
%                 be performed when the user presses a mouse button.
%%%
% * *moveAction*: This method must be implemented with the procedures to
%                 be performed when the user moves the mouse.
%%%
% * *upAction*: This method must be implemented with the procedures to
%               be performed when the user releases the mouse button that
%               was pressed.
%%%
% * *scrollAction*: This method must be implemented with the procedures to
%                   be performed when the user uses the mouse scroll.
%
% The constructor of the abstract *Emouse* class has 2 input arguments: 
%%%
% * The handle to the target *figure* object (dialog).
%%%
% * The handle to the initial current *axes* object (canvas) in the
%   target figure. 
%
% These arguments must be provided by the client user.
% It is possible to have more than one *axes* in the *figure*.
% The current axes is updated to the axes found at the position of the
% mouse in the button down event.
% It is assumed that the Units property of the *figure* and all their
% *axes* are consistent.
%
classdef Emouse < handle
    %% Protected attributes
    properties (Access = protected)
        dialog = [];                % dialog (figure) associated to mouse events.
        canvas = [];                % canvas (axes) associated to mouse events.
        mainDialogName = [];        % name of the dialog where Emouse object was created.
        mouseButtonMode = 'up';     % Button mouse states, 'up' or 'down'.
        whichMouseButton = 'none';  % 'none', 'left', 'right', 'center', or 'double click' at button mouse down.
		verticalScrollCount = 0;    % Counter for scroll.
		scrollAllowed = false;      % Flag to scroll events.
        currentPosition = [];       % x and y coordinates of the current pointer position.
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        % Constructor method, intended to initialize an object of this
        % class.
        % This method associates the mouse button down, mouse move,
        % and mouse button up events on the target figure with
        % the private eButtonDown, eMouseMove, and eButtonUp methods,
        % respectively.
        % Input arguments:
        %  dlg: handle to the target figure object (dialog).
        %  cnvs: handle to the initial axes (canvas) of the target figure.
        function this = Emouse(dlg,cnvs)
            this.dialog = dlg;
            this.canvas = cnvs;
            set(this.dialog, 'WindowButtonDownFcn', @this.eButtonDown);
            set(this.dialog, 'WindowButtonMotionFcn', @this.eMouseMove);
            set(this.dialog, 'WindowButtonUpFcn', @this.eButtonUp);
			set(this.dialog, 'WindowScrollWheelFcn', @this.eUseScroll);
            
            if isempty(this.mainDialogName) % Prevents unwanted changes in 'mainDialogName' property
                dlgName = get(this.dialog,'Name');
                this.mainDialogName = dlgName;
            end
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % This method must be implemented by a client subclass with the
        % procedures to be performed when the user presses a mouse button.
        downAction(this)

        %------------------------------------------------------------------
        % This method must be implemented by a client subclass with the
        % procedures to be performed when when the user moves the mouse.
        moveAction(this)

        %------------------------------------------------------------------
        % This method must be implemented by a client subclass with the
        % procedures to be performed when the when the user releases the
        % mouse button that was pressed.
        upAction(this)
		
		%------------------------------------------------------------------
        % This method must be implemented by a client subclass with the
        % procedures to be performed when the when the user utilizes the
        % mouse scroll.
        scrollAction(this)
        
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % This method is a callback function associated with mouse button
        % down event on the target canvas.
        % The method finds, in the list of axes (canvases) of the 
        % target figure (dialog), the axes (canvas) in which the button
        % down position is located.
        % The method also determines which button was pressed, updates the
        % whichMouseButton property with this information, sets the
        % mouseButtonMode property to down, sets the current position to
        % the mouse button down position, and calls the abstract
        % downAction method.
        function eButtonDown(this,~,~)
            % Check if gcf is the dialog is the one where Emouse object was 
            % created. If not, do nothing.
            if ~strcmp(get(gcf,'Name'),this.mainDialogName)
                return
            end
            
            % Change units to pixels
            set(this.dialog,'Units','pixels');
            set(gca,'Units','pixels');
            
            % Get click position coordinates inside the main dialog
            figPt = get(this.dialog, 'CurrentPoint');
            
            % Get canvas borders (in pixels)
            dfltUnits = get(this.canvas,'Units');
            set(this.canvas,'Units','pixels');
            limits = get(this.canvas,'Position');
            left = limits(1);
            right = limits(1) + limits(3);
            bottom = limits(2);
            top = limits(2) + limits(4);
            set(this.canvas,'Units',dfltUnits);
            
            % Check if the click was inside the canvas
            if (figPt(1) >= left) && (figPt(1) <= right) && (figPt(2) >= bottom) && (figPt(2) <= top)
                this.canvas = gca;
            else
                return
            end
            
            % Get which button was pressed.
            this.whichMouseButton = get(this.dialog, 'SelectionType');

            if strcmp(this.whichMouseButton,'alt')
                this.whichMouseButton = 'right';
            end
            if strcmp(this.whichMouseButton,'normal')
                this.whichMouseButton = 'left';
            end
            if strcmp(this.whichMouseButton,'extend')
                this.whichMouseButton = 'center';
            end           
            if strcmp(this.whichMouseButton,'open')
                this.whichMouseButton = 'double click';
            end
            
            % Set button mode as down, get button down location, and
            % call client (subclass) button down action method.
            this.mouseButtonMode = 'down';
            pt = get(this.canvas, 'CurrentPoint');
            xP = pt(:, 1);
            yP = pt(:, 2);
			zP = pt(:, 3);
            this.currentPosition = [xP yP zP]';
            this.downAction();
        end
        
        %------------------------------------------------------------------
        % This method is a callback function associated with mouse move
        % events on the target figure (dialog).
        % It sets the current position to the current mouse position on
        % the target axes (canvas) and calls the abstract moveAction
        % method.
        function eMouseMove(this,~,~)
            % Check if gcf is the dialog is the one where Emouse object was 
            % created. If not, do nothing.
            if ~strcmp(get(gcf,'Name'),this.mainDialogName)
                return
            end
            
            % Get current mouse location, and call client (subclass)
            % mouse move action method.
            pt = get(this.canvas, 'CurrentPoint');
            xP = pt(:, 1);
            yP = pt(:, 2);
			zP = pt(:, 3);
            this.currentPosition = [xP yP zP]';
            this.moveAction();
        end
        
        %------------------------------------------------------------------
        % This method is a callback function associated with mouse button
        % up events on the target figure (dialog).
        % It sets the mouseButtonMode property to up, sets the current
        % position to the mouse button up position on the target axes
        % (canvas), and calls the abstract upAction method.
        function eButtonUp(this,~,~)
            % Check if gcf is the dialog is the one where Emouse object was 
            % created. If not, do nothing.    
            if ~strcmp(get(gcf,'Name'),this.mainDialogName)
                return
            end
            
            % Do nothing if button down event was not on a canvas.
            if strcmp(this.whichMouseButton,'none')
                return
            end
            
            % Set button mode as up, get button up location, and
            % call client (subclass) button up action method.
            this.mouseButtonMode = 'up';
            pt = get(this.canvas, 'CurrentPoint');
            xP = pt(:, 1);
            yP = pt(:, 2);
			zP = pt(:, 3);
            this.currentPosition = [xP yP zP]';
            this.upAction();

            % Reset mouse button type for next sequence of 
            % button down - mouse move - button up events.
            this.whichMouseButton = 'none';
        end
        
		%------------------------------------------------------------------
		% This method is a callback function associated with mouse scroll
        % events on the target figure (dialog).
        % It sets the current position to the current mouse position on
        % the window coordinate system and calls the abstract scrollAction
        % method.
		function eUseScroll(this,~,event)
            % Check if gcf is the dialog is the one where Emouse object was 
            % created. If not, do nothing.    
            if ~strcmp(get(gcf,'Name'),this.mainDialogName)
                return
            end
            
            % Get the scroll intensity.
			this.verticalScrollCount = event.VerticalScrollCount;
            
            % Change units to pixels
            set(this.dialog,'Units','pixels');
            set(gca,'Units','pixels');
            
            % Get click position coordinates inside the main dialog
            figPt = get(this.dialog, 'CurrentPoint');
            
            % Get canvas borders
            dfltUnits = get(this.canvas,'Units');
            set(this.canvas,'Units','pixels');
            limits = get(this.canvas,'Position');
            left = limits(1);
            right = limits(1) + limits(3);
            bottom = limits(2);
            top = limits(2) + limits(4);
            set(this.canvas,'Units',dfltUnits);
            
            % Check if the click was inside the canvas
            if (figPt(1) >= left) && (figPt(1) <= right) && (figPt(2) >= bottom) && (figPt(2) <= top)
                this.canvas = gca;
            else
                return
            end
            
            % Normal zoom obtains negative values for this counter and vice
            % versa
            if (this.verticalScrollCount > 0)
                direction = 'minus';
                this.scrollAllowed = true;
            elseif (this.verticalScrollCount < 0)
                direction = 'plus';
                this.scrollAllowed = true;
            end
            
            % RLR: Need access to anm to check for real 2D canvas
            mdata = guidata(findobj('Tag','GUI_Main'));
            anm = get(mdata.popupmenu_Anm,'value');
            
            % Executes if scroll event happens on a canvas. Get scroll 
			% window location and call client (subclass) scroll action
			% method.
            if this.scrollAllowed
                if is2D(this.canvas) && (anm < 4) % Avoid entering here in 3D models when visualization in on XY plane
                    pt = get(this.canvas, 'CurrentPoint');
                    cx = pt(1, 1);
                    cy = pt(1, 2);
                    this.scrollAction(direction,cx,cy);
                else
                    pt = get(this.dialog, 'CurrentPoint');
                    wx = pt(1);
                    wy = pt(2);
                    this.scrollAction(direction,wx,wy);
                end
            end
            
            % Reset scroll and target canvas for next events.
            this.scrollAllowed = false;
        end
		
        %------------------------------------------------------------------
        % Initializes property values of an Emouse object.
        function this = clean(this)
            this.dialog = [];
            this.canvas = [];
            this.mainDialogName = [];
            this.mouseButtonMode = 'up';
            this.whichMouseButton = 'none';
			this.verticalScrollCount = 0;
			this.scrollAllowed = false;
            this.currentPosition = [];
        end
    end
end
