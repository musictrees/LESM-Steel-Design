%% Emouse_3D_InclSupp class
%
%% Description
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <emouse.html *Emouse*> in the <main.html LESM (Linear Elements
% Structure Model)> program. This sub-class implements abstract methods,
% defined in super-class *Emouse*, that deal with the visualization and 
% handling of 3D inclined supports, on the GUI_Supports figure.
%
classdef Emouse_3D_InclSupp < Emouse
    %%
    % <emouse.html See documentation on *Emouse* super-class>.
  
    %% Public attributes
    properties (Access = public)
        % Rotate variables
        rotIniX = 0;             % x window coordinate when RB button is clicked to rotate.
        rotIniY = 0;             % y window coordinate when RB button is clicked to rotate.
        rotIniAz = 0;            % Azimuth when RB button is clicked to rotate.
        rotIniEl = 0;            % Elevation when RB button is clicked to rotate.
        
        % Double click variables
        originalAxesPos = [];    % Original axes position on window coordinates.
        originalAz = 0;          % Original azimuth of the view.
        originalEl = 0;          % Original elevation of the view.
        
        % Pick new direction variables
        numTol = 0.1;                  % numeric tolerance for click point
        dir_x_beingCollected = false;  % flag for vx direction being collected
        dir_x_1stPt = [];              % vx initial point
        dir_x_2ndPt = [];              % vx final point
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Emouse_3D_InclSupp(fig,axes)
            this = this@Emouse(fig,axes);
            
            % Properties for double click on mouse
            originalUnits = get(axes, 'Units');
            set(axes, 'Units', 'pixels');
            this.originalAxesPos = get(axes, 'Position');
            set(axes, 'Units', originalUnits);
            [this.originalAz,this.originalEl] = view(axes);
        end
    end
    
    %% Abstract methods
    % Implementation of the abstract methods declared in super-class <Emouse.html *Emouse*>.
    methods        
        %------------------------------------------------------------------
        % Action executed when an mouse button is pressed.
        function downAction(this)
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
                return
            end
            
            % Double click
            if strcmp(this.whichMouseButton,'double click')
                this.doubleClick()
                return
            end
            
            % Pick new direction
            if strcmp(this.whichMouseButton,'left')
                % Get 3D click line
                line_coords = (this.currentPosition)';
                
                switch this.dir_x_beingCollected
                    case true
                        % Get handle to GUI_Supports
                        mdata_supp = guidata(findobj('tag','GUI_Supports'));
                        vx = [this.dir_x_2ndPt(1) - this.dir_x_1stPt(1),...
                              this.dir_x_2ndPt(2) - this.dir_x_1stPt(2),...
                              this.dir_x_2ndPt(3) - this.dir_x_1stPt(3)];
                        vx = vx/norm(vx);
                        if get(mdata_supp.radiobutton_DirectionVector,'value')
                            set(mdata_supp.edit_InclSupp_1,'string',num2str(vx(1),3))
                            set(mdata_supp.edit_InclSupp_2,'string',num2str(vx(2),3))
                            set(mdata_supp.edit_InclSupp_3,'string',num2str(vx(3),3))
                        else
                            dir_1 = vx(1:2);
                            aux_norm = norm(dir_1);
                            dir_1 = dir_1 / aux_norm;
                            thetaZ = acos(dir_1(1));
                            if dir_1(2) < 0
                                thetaZ = - thetaZ;
                            end
                            
                            % Convert from rad to degrees
                            thetaZ = thetaZ * 180 / pi;
                            if abs(thetaZ) < 10^-10 || isnan(thetaZ)
                                thetaZ = 0;
                            end
                            
                            % Compute direction vector as angle
                            dir_2 = [aux_norm, vx(3)];
                            dir_2 = dir_2 / norm(dir_2);
                            thetaY = asin(dir_2(2));
                            
                            % Convert from rad to degrees
                            thetaY = thetaY * 180 / pi;
                            if abs(thetaY) < 10^-10
                                thetaY = 0;
                            end
                            
                            if strcmp(get(mdata_supp.edit_InclSupp_1,'string'),'on')
                                thetaX = 0;
                                set(mdata_supp.edit_InclSupp_1,'string',num2str(thetaX,4))
                            end
                            set(mdata_supp.edit_InclSupp_2,'string',num2str(thetaY,4))
                            set(mdata_supp.edit_InclSupp_3,'string',num2str(thetaZ,4))
                        end
                        % Draw updated model
                        GUI_Supports('updateInclSupp3DCnvs',mdata_supp,true);
                        
                        this.dir_x_beingCollected = false;
                        this.dir_x_1stPt = [];
                        this.dir_x_2ndPt = [];
                        delete(findobj('tag','dynamicLine_InclSupp_dir_x'))
                        delete(findobj('tag','text_InclSupp_dir_x'))
                    case false
                        pt = zeros(1,3);
                        if auxModelFctn('isPointInLine3D',{pt,line_coords,this.numTol,true})
                            this.dir_x_beingCollected = true;
                            this.dir_x_1stPt = pt;
                        end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse pointer moves on the canvas.
        function moveAction(this)
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
            
            % Pick new direction
            if this.dir_x_beingCollected
                % Get 3D click line
                line_coords = (this.currentPosition)';
                aux = find(abs(line_coords(2,:)) == 1);
                if isempty(aux)
                    return
                elseif aux == 1
                    coords = auxMouseFctn('spatial2Plane',this,{'x',line_coords(2,aux)});
                    id = [2,3];
                elseif aux == 2
                    coords = auxMouseFctn('spatial2Plane',this,{'y',line_coords(2,aux)});
                    id = [1,3];
                elseif aux == 3
                    coords = auxMouseFctn('spatial2Plane',this,{'z',line_coords(2,aux)});
                    id = [1,2];
                end
                
                delete(findobj('tag','dynamicLine_InclSupp_dir_x'))
                delete(findobj('tag','text_InclSupp_dir_x'))

                X = [this.dir_x_1stPt(1), coords(1)];
                Y = [this.dir_x_1stPt(2), coords(2)];
                Z = [this.dir_x_1stPt(3), coords(3)];
                plot3(X,Y,Z,...
                     'color',[1 0 0],'LineStyle','--',...
                     'tag','dynamicLine_InclSupp_dir_x');
                scatter3(coords(1),coords(2),coords(3),25,[1,0,0],'filled',...
                         'tag','dynamicLine_InclSupp_dir_x')
                     
                this.dir_x_2ndPt = coords;
                     
                r = [rem(coords(id(1)),0.05), rem(coords(id(2)),0.05)];
                for i = 1:2
                    if abs(r(i)) <= 0.025
                        this.dir_x_2ndPt(id(i)) = coords(id(i)) - r(i);
                    elseif r(i) < 0
                        this.dir_x_2ndPt(id(i)) = coords(id(i)) - r(i) - 0.05;
                    else
                        this.dir_x_2ndPt(id(i)) = coords(id(i)) - r(i) + 0.05;
                    end
                end
                
                txt = sprintf('(%.2f,%.2f,%.2f)',this.dir_x_2ndPt(1),this.dir_x_2ndPt(2),this.dir_x_2ndPt(3));
                text(coords(1)*1.05,coords(2)*1.05,coords(3)*1.05,txt,'Color',[1,0,0],'tag','text_InclSupp_dir_x');
            end
        end
        
        %------------------------------------------------------------------
        % Action executed when an mouse button is unpressed.
        function upAction(~)
        end
        
        %------------------------------------------------------------------
        % Action executed when mouse scroll is used.
        % Executes zoom (camera) using the mouse scroll.
        function scrollAction(~,direction,~,~)
            % Get handle to GUI_Supports
            mdata_supp = guidata(findobj('tag','GUI_Supports'));
%             if get(mdata_supp.checkbox_Dy,'value') || get(mdata_supp.checkbox_Dz,'value') ||...
%                get(mdata_supp.checkbox_Ky,'value') || get(mdata_supp.checkbox_Kz,'value')
           
                if get(mdata_supp.radiobutton_Angles,'value')
                    thetaX = str2double(get(mdata_supp.edit_InclSupp_1,'string'));
                    if isnan(thetaX)
                        thetaX = 0;
                    else
                        r = rem(thetaX,5);
                        if abs(r) <= 2.5
                            thetaX = thetaX - r;
                        elseif r < 0
                            thetaX = thetaX - r - 5;
                        else
                            thetaX = thetaX - r + 5;
                        end
                    end
                    switch direction
                        case 'plus'
                            thetaX = thetaX + 5;
                        case 'minus'
                            thetaX = thetaX - 5;
                    end
                    set(mdata_supp.edit_InclSupp_1,'string',num2str(thetaX,4))
                    
                else % direction vector
                    
                    % Compute direction vector as angle
                    dx = str2double(get(mdata_supp.edit_InclSupp_1,'string'));
                    if isnan(dx)
                        dx = 0;
                    end
                    dy = str2double(get(mdata_supp.edit_InclSupp_2,'string'));
                    if isnan(dy)
                        dy = 0;
                    end
                    dz = str2double(get(mdata_supp.edit_InclSupp_3,'string'));
                    if isnan(dz)
                        dz = 0;
                    end
                    
                    if dx == 0 && dy == 0 && dz == 0
                        dir_1 = [1, 0];
                        dir_2 = [1, 0];
                        dir = [1, 0, 0];
                    else
                        dir_1 = [dx, dy];
                        dir_2 = [norm(dir_1), dz];
                        dir_1 = dir_1 / norm(dir_1);
                        dir_2 = dir_2 / norm(dir_2);
                        dir = [dx, dy, dz] / norm([dx, dy, dz]);
                    end
                    
                    % Compute thetaZ
                    thetaZ = acos(dir_1(1));
                    if dir_1(2) < 0
                        thetaZ = - thetaZ;
                    end
                    
                    % Compute thetaY
                    thetaY = asin(dir_2(2));
                    
                    x = dir;
                    if abs(abs(thetaY)-(pi/2)) < 10^-6
                        thetaZ = 0;
                        if thetaY > 0
                            z1 = [-1, 0, 0];
                        else
                            z1 = [1, 0, 0];
                        end
                    else
                        z1 = [0, 0, 1];
                    end
                    
                    % Compute local axis y and z directions
                    y1 = cross(z1,x);
                    y1 = y1 / norm(y1);
                    
                    dyx = str2double(get(mdata_supp.edit_InclSupp_4,'string'));
                    if isnan(dyx)
                        dyx = 0;
                    end
                    dyy = str2double(get(mdata_supp.edit_InclSupp_5,'string'));
                    if isnan(dyy)
                        dyy = 0;
                    end
                    dyz = str2double(get(mdata_supp.edit_InclSupp_6,'string'));
                    if isnan(dyz)
                        dyz = 0;
                    end
                    if dyx == 0 && dyy == 0 && dyz == 0
                        y2 = y1;
                    else
                        vy = [dyx, dyy, dyz] / norm([dyx, dyy, dyz]);
                        z2 = cross(x,vy);
                        z2 = z2 / norm(z2);
                        y2 = cross(z2,x);
                    end
                    
                    w = cross(y1,y2);
                    aux = x * w';
                    if abs(aux) < 10^-10
                        if (y1/y2) >= 0
                            thetaX = 0;
                        else
                            thetaX = 180;
                        end
                    else
                        thetaX = acos((y1*y2')/(norm(y1)*norm(y2))) * 180 / pi;
                        if abs(thetaX) < 10^-10
                            thetaX = 0;
                        elseif aux < 0
                            thetaX = - thetaX;
                        end
                        if abs(thetaX) >= 10
                            thetaX = round(thetaX,1);
                        end
                    end
                    
                    r = rem(thetaX,5);
                    if abs(r) <= 2.5
                        thetaX = thetaX - r;
                    elseif r < 0
                        thetaX = thetaX - r - 5;
                    else
                        thetaX = thetaX - r + 5;
                    end
                    
                    switch direction
                        case 'plus'
                            thetaX = thetaX + 5;
                        case 'minus'
                            thetaX = thetaX - 5;
                    end
                    
                    thetaX = thetaX * pi / 180;
                    vz = [sin(thetaX)*sin(thetaZ) - cos(thetaX)*cos(thetaZ)*sin(thetaY), - sin(thetaX)*cos(thetaZ) - cos(thetaX)*sin(thetaY)*sin(thetaZ), cos(thetaX)*cos(thetaY)];
                    % Compute direction vector vy
                    vy = cross(vz,x);
                    vy = vy / norm(vy);
                    set(mdata_supp.edit_InclSupp_4,'string',num2str(vy(1),3));
                    set(mdata_supp.edit_InclSupp_5,'string',num2str(vy(2),3));
                    set(mdata_supp.edit_InclSupp_6,'string',num2str(vy(3),3));
                end
                % Draw updated model
                GUI_Supports('updateInclSupp3DCnvs',mdata_supp,true);
            %end
        end
        
        %------------------------------------------------------------------
        function doubleClick(this)
            % Make sure operation is performed on this canvas
            axes(this.canvas)
            
            % Set the object position.
            dfltUnits = get(this.canvas, 'Units');
            set(this.canvas, 'Units', 'pixels');
            set(this.canvas, 'Position', this.originalAxesPos);
            set(this.canvas, 'Units', dfltUnits);
            
            % Sets the original parameters of visualization.
            view(gca,[this.originalAz,this.originalEl]);
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