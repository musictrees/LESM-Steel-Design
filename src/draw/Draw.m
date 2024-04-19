%% Draw class
%
%% Description
%
% This is a handle super-class for the definition of a drawing object.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <draw_truss2d.html 2D truss model draw>.
% * <draw_frame2d.html 2D frame model draw>.
% * <draw_grillage.html Grillage model draw>.
% * <draw_truss3d.html 3D truss model draw>.
% * <draw_frame3d.html 3D frame model draw>.
%
classdef Draw < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        mdl  = [];   % handle to an object of the Model class
        size = 0;    % drawing size parameter
        az   = 0;    % camera default viewpoint azimuth
        elev = 0;    % camera default viewpoint elevation
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw(az,elev)
            draw.az   = az;
            draw.elev = elev;
        end
    end
    
    %% Class (static) auxiliary functions
    methods (Static)
        %------------------------------------------------------------------
        % Get list of selected elements IDs from the interface.
        function ids = getSelectedElements(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            n_str = get(mdata.edit_ElementResults,'string');
            if strcmp(n_str,'all') || strcmp(n_str,'ALL') || strcmp(n_str,'All')
                ids = 1:draw.mdl.nel;
            else
                [~,ids] = draw.readStr(n_str,draw.mdl.nel);
            end
        end
        function [flag,output] = readStr(str,max)
            output    = zeros(1,max);
            count     = 0;       % counter for output index
            numFlag   = false;   % flag for number being read
            vctrFlag  = false;   % flag for vector being read
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
        end
        
        %------------------------------------------------------------------
        % Plots 2D axis symbol with defined origin coordinates and color.
        % Input arguments:
        %  x: origin coordinate on the X axis
        %  y: origin coordinate on the Y axis
        %  s: size parameter (length of the arrows)
        %  c: color (RGB array)
        %  ang: rotation [rad]
        %  tag: name to identify graphic objects
        %  upperCaseFlag: flag for text with upper case
        function axis2D(draw,x,y,s,c,ang,tag,upperCaseFlag)
            if nargin < 8
                upperCaseFlag = true;
            end

            % Get colors
            if size(c,1) == 2
                mark_clr = [0,0,0];
                clr_1 = c(1,:);
                clr_2 = c(2,:);
            elseif size(c,1) == 1
                mark_clr = c;
                clr_1 = c;
                clr_2 = c;
            end
            
            % Plot origin mark
            scatter(x,y,s/10,mark_clr,'filled','tag',tag);
            
            % Plot x axis
            draw.arrow2D(x + s*cos(ang), y + s*sin(ang), s, s/5, s/5, pi + ang, clr_1,tag);
            
            % Plot y axis
            draw.arrow2D(x - s*sin(ang), y + s*cos(ang), s, s/5, s/5, -pi/2 + ang, clr_2,tag);
            
            % Plot texts
            txt_pos = [x, y; x, y] + [s/2, -s/5; -s/5,  s/2] * [cos(ang), sin(ang); -sin(ang), cos(ang)];
            if upperCaseFlag
                txt = ['X'; 'Y'];
            else
                txt = ['x'; 'y'];
            end
            text(txt_pos(1,1), txt_pos(1,2), txt(1,:), 'Color', clr_1, 'tag', tag);
            text(txt_pos(2,1), txt_pos(2,2), txt(2,:), 'Color', clr_2, 'tag', tag);
        end
        
        %------------------------------------------------------------------
        % Plots 3D axis symbol with defined origin coordinates, direction 
        % and color.
        % Input arguments:
        %  x: origin coordinate on the X axis
        %  y: origin coordinate on the Y axis
        %  z: origin coordinate on the Z axis
        %  s: size parameter (length of the arrows)
        %  c: color (RGB array)
        %  dir: matrix of axis directions
        %  tag: name to identify graphic objects
        %  upperCaseFlag: flag for text with upper case
        function axis3D(draw,x,y,z,s,c,dir,tag,upperCaseFlag)
            if nargin < 9
                upperCaseFlag = true;
            end

            % Check if directions were provided. If not, use global dir
            if isempty(dir)
                dir = eye(3);
            end
            
            % Get axis directions
            dx = dir(1,:)/norm(dir(1,:));
            dy = dir(2,:)/norm(dir(2,:));
            dz = dir(3,:)/norm(dir(3,:));
            
            % Get colors
            if size(c,1) == 3
                mark_clr = [0,0,0];
                clr_1 = c(1,:);
                clr_2 = c(2,:);
                clr_3 = c(3,:);
            else
                mark_clr = c(1,:);
                clr_1 = c(1,:);
                clr_2 = c(1,:);
                clr_3 = c(1,:);
            end
            
            % Plot origin mark
            scatter3(x,y,z,s*50,mark_clr,'filled','tag',tag);
            
            % Plot x axis
            draw.arrow3D(draw,x+s*dx(1),y+s*dx(2),z+s*dx(3),s,s/4.5,s/5,[dx;dy;dz],clr_1,tag,true);
            
            % Plot y axis
            draw.arrow3D(draw,x+s*dy(1),y+s*dy(2),z+s*dy(3),s,s/4.5,s/5,[dy;dz;dx],clr_2,tag,true);
            
            % Plot z axis
            draw.arrow3D(draw,x+s*dz(1),y+s*dz(2),z+s*dz(3),s,s/4.5,s/5,[dz;dx;dy],clr_3,tag,true);
            
            % Plot texts
            txt_pos = [x, y, z; x, y, z; x, y, z] +...
                      [s/2.75, -s/4, 0; -s/4,  s/2.75, 0; 0, -s/5, s/2] *...
                      [dx; dy; dz];
            if upperCaseFlag
                txt = ['X'; 'Y'; 'Z'];
            else
                txt = ['x'; 'y'; 'z'];
            end
            text(txt_pos(1,1), txt_pos(1,2), txt_pos(1,3), txt(1,:), 'Color', clr_1, 'tag', tag);
            text(txt_pos(2,1), txt_pos(2,2), txt_pos(2,3), txt(2,:), 'Color', clr_2, 'tag', tag);
            text(txt_pos(3,1), txt_pos(3,2), txt_pos(3,3), txt(3,:), 'Color', clr_3, 'tag', tag);
        end
        
        %------------------------------------------------------------------
        % Plots 3D plane with defined coordinates of a central point,
        % normal direction and color.
        % Input arguments:
        %  x: central point coordinate on the X axis
        %  y: central point coordinate on the Y axis
        %  z: central point coordinate on the Z axis
        %  s: size parameter (plane sides)
        %  c: color (RGB array)
        %  dir: matrix of axis directions
        %  tag: name to identify graphic objects
        function plane3D(x,y,z,s,c,dir,tag)
            % Check if directions were provided. If not, use global dir
            if isempty(dir)
                dir = eye(3);
            end

            % Get axis directions
            dy = dir(2,:)/norm(dir(2,:));
            dz = dir(3,:)/norm(dir(3,:)); 
            
            % Compute four corner points coordinates
            center = [x, y, z];
            pt_1 = center - (s/2)*dy - (s/2)*dz;
            pt_2 = pt_1 + s*dy;
            pt_3 = pt_2 + s*dz;
            pt_4 = pt_3 - s*dy;
            
            % Assemble coordinate vectors
            X = [pt_1(1), pt_2(1), pt_3(1), pt_4(1)];
            Y = [pt_1(2), pt_2(2), pt_3(2), pt_4(2)];
            Z = [pt_1(3), pt_2(3), pt_3(3), pt_4(3)];
            
            % Plot plane
            fill3(X,Y,Z,c,'EdgeColor',[0,0,0],'FaceAlpha',0.4,'tag',tag);
        end
 
        %------------------------------------------------------------------
        % Plots a circle with defined center coordinates, radius and color.
        % This method is used to draw hinges on 2D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: circle radius
        %  c: color (RGB vector)
        function circle(x,y,r,c,tag)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            if nargin == 4
                plot(xcirc, ycirc, 'color', c);
            elseif nargin == 5
                plot(xcirc, ycirc, 'color', c, 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a sphere with defined center coordinates and radius.
        % This method is used to draw hinges on 3D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  z: center coordinate on the Z axis
        %  r: sphere radius
        function s = sphere(x,y,z,r,tag)
            [a, b, c] = sphere;
            s = surf(a * r + x, b * r + y, c * r + z);
            if nargin == 4
                set(s, 'Edgecolor', [0,0,0],'FaceColor', [1,1,1]);
            elseif nargin == 5
                set(s, 'Edgecolor', [0,0,0],'FaceColor', [1,1,1], 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a square with defined center coordinates, side length and
        % color.
        % This method is used to draw nodal points and rotation constraints
        % on 2D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  S: side length
        %  c: color (RGB vector)
        function square(x,y,S,c,tag)
            s = S/2;
            
            X = [x - s , x + s , x + s , x - s];
            Y = [y - s , y - s , y + s , y + s];
            if nargin == 4
                fill(X, Y, c);
            elseif nargin == 5
                fill(X, Y, c, 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a cube with defined center coordinates, side length and
        % color.
        % This method is used to draw nodal points and rotation constraints
        % on 3D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  z: center coordinate on the Z axis
        %  S: side length
        %  c: color (RGB vector)
        function cube(x,y,z,S,c,tag)
            s = S/2;
            
            % Plot face 1
            X = [x + s, x - s, x - s, x + s];
            Y = [y + s, y + s, y - s, y - s];
            Z = [z + s, z + s, z + s, z + s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
            
            % Plot face 2
            X = [x + s, x - s, x - s, x + s];
            Y = [y + s, y + s, y - s, y - s];
            Z = [z - s, z - s, z - s, z - s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
            
            % Plot face 3
            X = [x + s, x - s, x - s, x + s];
            Y = [y + s, y + s, y + s, y + s];
            Z = [z + s, z + s, z - s, z - s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
            
            % Plot face 4
            X = [x + s, x - s, x - s, x + s];
            Y = [y - s, y - s, y - s, y - s];
            Z = [z + s, z + s, z - s, z - s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
            
            % Plot face 5
            X = [x + s, x + s, x + s, x + s];
            Y = [y + s, y + s, y - s, y - s];
            Z = [z + s, z - s, z - s, z + s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
            
            % Plot face 6
            X = [x - s, x - s, x - s, x - s];
            Y = [y + s, y + s, y - s, y - s];
            Z = [z + s, z - s, z - s, z + s];
            if nargin == 6
                fill3(X, Y, Z, c,'tag',tag);
            else
                fill3(X, Y, Z, c);
            end
            hold on
        end
        
        %------------------------------------------------------------------
        % Plots a triangle with defined top coordinates, height, base,
        % orientation, and color.
        % This method is used to draw translation constraints on 2D models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: triangle height
        %  B: triangle base
        %  ang: angle (in radian) between the axis of symmetry and the
        %       horizontal direction (counterclockwise) - 0 rad when
        %       triangle is pointing left
        %  c: color (RGB vector)
        function triangle(x,y,h,B,ang,c,tag)
            b = B/2;
            
            cx = cos(ang);
            cy = sin(ang);
            
            X = [x, x + h * cx + b * cy, x + h * cx - b * cy];
            Y = [y, y + h * cy - b * cx, y + h * cy + b * cx];
            
            if nargin == 6
                fill(X, Y, c);
            elseif nargin == 7
                fill(X, Y, c, 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 2D displacement spring with defined top coordinates,
        % height, base, and color, in the direction of axis X.
        % This method is used to draw translation spring constraints on 2D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: spring height
        %  c: color (RGB vector)
        %  tag: graphic object identification tag
        %  ang: rotation about z axis
        function springX_2D(x,y,h,c,tag,ang)
            if nargin < 6
                ang = 0;
            end
            
            % Number of coils
            nc = 4;
            % Height per coil
            hc = h/(1.375*nc);
            % Height of the straight parts
            hr = 0.1875*nc*hc;
            % Base
            b = 1.5*h/nc;
            % Coil properties
            alpha = (-pi/2):(nc*pi/15):(nc*2*pi-pi/2);
            radius = h/13.75;
            % Spring coordinates
            xs = (x-hr):(-hc*nc/30):(x-hr-hc*nc);
            ys = y - radius*cos(alpha);
            % Base coordinates
            xb = [(x-h) (x-h)];
            yb = [(y-b/2) (y+b/2)];
            
            % Assemble coordinate vectors
            X = [x xs (x-h) xb];
            Y = [y ys y yb];
            
            % Rotate drawing coordinates, if needed
            if ang ~= 0
                coords = [X; Y; ones(1,length(X))];
                
                % Rotation transformation matrix
                R = [cos(ang) -sin(ang)  0
                     sin(ang)  cos(ang)  0
                     0          0        1];
                T1 = [1 0 x
                      0 1 y
                      0 0 1];
                
                T2 = [1 0 -x
                      0 1 -y
                      0 0  1];
                
                coords = T1*R*T2*coords;
                X = coords(1,:);
                Y = coords(2,:);
            end
            
            if nargin == 4
                plot(X,Y,'color',c,'Linewidth',1.2)
            elseif nargin == 5 || nargin == 6
                plot(X,Y,'color',c,'Linewidth',1.2,'tag',tag)
            end
               
        end    
        
        %------------------------------------------------------------------
        % Plots a 2D displacement spring with defined top coordinates,
        % height, base, and color, in the direction of axis Y.
        % This method is used to draw translation spring constraints on 2D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: spring height
        %  c: color (RGB vector)
        %  tag: graphic object identification tag
        %  ang: rotation about z axis
        function springY_2D(x,y,h,c,tag,ang)
            if nargin < 6
                ang = 0;
            end
            
            % Number of coils
            nc = 4;
            % Height per coil
            hc = h/(1.375*nc);
            % Height of the straight parts
            hr = 0.1875*nc*hc;
            % Base
            b = 1.5*h/nc;
            % Coil properties
            alpha = (-pi/2):(nc*pi/15):(nc*2*pi-pi/2);
            radius = h/13.75;
            % Spring coordinates
            xs = x - radius*cos(alpha);
            ys = (y-hr):(-hc*nc/30):(y-hr-hc*nc);
            % Base coordinates
            xb = [(x-b/2) (x+b/2)];
            yb = [(y-h) (y-h)];
            
            % Assemble coordinate vectors
            X = [x xs x xb];
            Y = [y ys (y-h) yb];
            
            % Rotate drawing coordinates, if needed
            if ang ~= 0
                coords = [X; Y; ones(1,length(X))];
                
                % Rotation transformation matrix
                R = [cos(ang) -sin(ang)  0
                     sin(ang)  cos(ang)  0
                     0          0        1];
                T1 = [1 0 x
                      0 1 y
                      0 0 1];
                
                T2 = [1 0 -x
                      0 1 -y
                      0 0  1];
                
                coords = T1*R*T2*coords;
                X = coords(1,:);
                Y = coords(2,:);
            end
            
            if nargin == 4
                plot(X,Y,'color',c,'Linewidth',1.2)
            elseif nargin == 5 || nargin == 6
                plot(X,Y,'color',c,'Linewidth',1.2,'tag',tag)
            end
            
        end    
        
        %------------------------------------------------------------------
        % Plots a 2D rotational spring with defined center coordinates,
        % height, base, and color, about the direction of axis Z.
        % This method is used to draw rotational spring constraints on 2D
        % models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  h: spring height
        %  c: color (RGB vector)
        function rotSpringZ_2D(x,y,h,c,tag)
            % Spring center height
            sch = h/2;
            % Number of coils
            nc = 3.25;
            % Coil properties
            alpha = 0:(nc*pi/25):(nc*2*pi);
            rmax = 0.8*sch;
            rmin = 0;
            radius = rmin:(rmax/50):rmax;
            % Spring coordinates
            xs = zeros(1,size(alpha,2));
            ys = zeros(1,size(alpha,2));
            for i = 1:size(alpha,2)
            xs(i) = x + radius(i)*sin(alpha(i));
            ys(i) = y + radius(i)*cos(alpha(i));
            end
            % Base coordinates
            xa = xs(size(alpha,2));
            ya = ys(size(alpha,2)) - sch;
            x0 = [-0.25*sch+xs(size(alpha,2)) 0.25*sch+xs(size(alpha,2))];
            y0 = [(ys(size(alpha,2))-sch) (ys(size(alpha,2))-sch)];

            X = [xs xa x0];
            Y = [ys ya y0];
            
            if nargin == 4
                plot(X,Y,'color',c,'Linewidth',1.1)
            elseif nargin == 5
                plot(X,Y,'color',c,'Linewidth',1.1,'tag',tag)
            end
            
        end    
        
        %------------------------------------------------------------------
        % Plots a 3D displacement spring with defined top coordinates,
        % height, base, and color, in the direction of axis X.
        % This method is used to draw displacement spring constraints on 3D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  z: top coordinate on the Z axis
        %  h: spring height
        %  c: color (RGB vector)
        function SpringX_3D(x,y,z,h,c,tag)
            % Number of coils
            nc = 4.5;
            % Spring height (per coil)
            hc = h/(1.5*nc);
            % Spring height (straight parts)
            hr = 0.25*nc*hc;
            % Base
            b = 1.2*h/nc;
            % Coil properties
            alpha = -pi/2:(nc*pi/25):(nc*2*pi-pi/2);
            radius = h/12;
            % Spring coordinates
            xs = (x-hr):(-nc*hc/50):(x-hr-nc*hc);
            ys = y - radius*cos(alpha);
            zs = z - radius*sin(alpha);
            % Base coordinates
            xb = [(x-h) (x-h) (x-h) (x-h)];
            yb = [(y-b/2) (y-b/2) (y+b/2) (y+b/2)];
            zb = [(z-b/2) (z+b/2) (z+b/2) (z-b/2)];

            X = [x (x-hr) xs (x-hr-nc*hc) (x-h)];
            Y = [y y ys y y];
            Z = [z z zs z z];
            
            if nargin == 5
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0])
            elseif nargin == 6
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2,'tag',tag)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0],'tag',tag)
            end
            
        end    
        
        %------------------------------------------------------------------
        % Plots a 3D displacement spring with defined top coordinates,
        % height, base, and color, in the direction of axis Y.
        % This method is used to draw displacement spring constraints on 3D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  z: top coordinate on the Z axis
        %  h: spring height
        %  c: color (RGB vector)
        function SpringY_3D(x,y,z,h,c,tag)
            % Number of coils
            nc = 4.5;
            % Spring height (per coil)
            hc = h/(1.5*nc);
            % Spring height (straight parts)
            hr = 0.25*nc*hc;
            % Base
            b = 1.2*h/nc;
            % Coil properties
            alpha = -pi/2:(nc*pi/25):(nc*2*pi-pi/2);
            radius = h/12;
            % Spring coordinates
            xs = x - radius*cos(alpha);
            ys = (y-hr):(-nc*hc/50):(y-hr-nc*hc);
            zs = z - radius*sin(alpha);
            % Base coordinates
            xb = [(x-b/2) (x-b/2) (x+b/2) (x+b/2)];
            yb = [(y-h) (y-h) (y-h) (y-h)];
            zb = [(z-b/2) (z+b/2) (z+b/2) (z-b/2)];

            X = [x x xs x x];
            Y = [y (y-hr) ys (y-hr-nc*hc) (y-h)];
            Z = [z z zs z z];
            
            if nargin == 5
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0])
            elseif nargin == 6
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2,'tag',tag)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0],'tag',tag)
            end
            
        end    
        
        %------------------------------------------------------------------
        % Plots a 3D displacement spring with defined top coordinates,
        % height, base, and color, in the direction of axis Z.
        % This method is used to draw displacement spring constraints on 3D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  z: top coordinate on the Z axis
        %  h: spring height
        %  c: color (RGB vector)
        function SpringZ_3D(x,y,z,h,c,tag)
            % Number of coils
            nc = 4.5;
            % Spring height (per coil)
            hc = h/(1.5*nc);
            % Spring height (straight parts)
            hr = 0.25*nc*hc;
            % Base
            b = 1.2*h/nc;
            % Coil properties
            alpha = -pi/2:(nc*pi/25):(nc*2*pi-pi/2);
            radius = h/12;
            % Spring coordinates
            xs = x - radius*cos(alpha);
            ys = y - radius*sin(alpha);
            zs = (z-hr):(-nc*hc/50):(z-hr-nc*hc);
            % Base coordinates
            xb = [(x-b/2) (x-b/2) (x+b/2) (x+b/2)];
            yb = [(y-b/2) (y+b/2) (y+b/2) (y-b/2)];
            zb = [(z-h) (z-h) (z-h) (z-h)];

            X = [x x xs x x];
            Y = [y y ys y y];
            Z = [z (z-hr) zs (z-hr-nc*hc) (z-h)];
            
            if nargin == 5
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0])
            elseif nargin == 6
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2,'tag',tag)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0],'tag',tag)
            end
            
        end
        
        %------------------------------------------------------------------
        % Plots a 3D displacement spring with defined top coordinates,
        % height, base, and color, on the specified direction.
        % This method is used to draw displacement spring constraints on 3D
        % models.
        % Input arguments:
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  z: top coordinate on the Z axis
        %  h: spring height
        %  dir: mtx of local axis direction vectors [dir_x; dir_y; dir_z]
        %  c: color (RGB vector)
        %  tag: graphic object identification tag
        function displSpring_3D(x,y,z,h,dir,c,tag)
            % Compute auxiliary parameters for drawing helix-shaped springs
            % Number of coils
            nc = 4.5;
            % Spring height (per coil)
            hc = h/(1.5*nc);
            % Spring height (straight parts - outside from helix)
            hr = 0.25*nc*hc;
            % Base
            b = 1.2*h/nc;
            % Coil properties
            alpha = (-pi/2:(nc*pi/25):(nc*2*pi-pi/2))';
            radius = h/12;
            
            % Get direction vectors
            dx = dir(1,:);
            dy = dir(2,:);
            dz = dir(3,:);
            
            % Initial point (top of the spring)
            pt_1 = [x, y, z];
            
            % Final point (bottom of the spring)
            pt_2 = pt_1 - h * dx;
            
            % Helix coordinates
            Ls = linspace(0,nc*hc,size(alpha,1))';
            aux_pt = pt_1 - hr * dx;
            
            pts = zeros(size(alpha,1),3);
            pts(:,1) = ones(size(pts,1),1) * aux_pt(1);
            pts(:,2) = ones(size(pts,1),1) * aux_pt(2);
            pts(:,3) = ones(size(pts,1),1) * aux_pt(3);
            
            pts = pts - Ls * dx - radius * cos(alpha) * dy - radius * sin(alpha) * dz;
            
            xs = pts(:,1)';
            ys = pts(:,2)';
            zs = pts(:,3)';
            
            % Assemble spirng coordinate vectors
            X = [x xs pt_2(1)];
            Y = [y ys pt_2(2)];
            Z = [z zs pt_2(3)];
            
            % Base coordinates
            pt_1_b = pt_1 - h * dx - b/2 * dy - b/2 * dz;
            pt_2_b = pt_1_b + b * dz;
            pt_3_b = pt_2_b + b * dy;
            pt_4_b = pt_3_b - b * dz;
            
            xb = [pt_1_b(1) pt_2_b(1) pt_3_b(1) pt_4_b(1)];
            yb = [pt_1_b(2) pt_2_b(2) pt_3_b(2) pt_4_b(2)];
            zb = [pt_1_b(3) pt_2_b(3) pt_3_b(3) pt_4_b(3)];
            
            % Check if a tag was provided
            if nargin == 6
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0])
            elseif nargin == 7
                plot3(X,Y,Z,'Color',c,'Linewidth',1.2,'tag',tag)
                hold on
                fill3(xb,yb,zb,[0.6,0.6,0.6],'EdgeColor',[0,0,0],'tag',tag)
            end
        end
        
        %------------------------------------------------------------------
        % Plots two spheres with defined center coordinates, radius and
        % color.
        % This method is used to draw rotational springs on 3D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  z: center coordinate on the Z axis
        %  r: sphere radius
        %  clr: color (RGB vector)
        function rotSpring_3D(x,y,z,r,clr,tag)
            [a, b, c] = sphere;
            s = surf(a * r + x, b * r + y, c * r + z);
            if nargin == 5
                set(s, 'Edgecolor', clr,'FaceColor', [1,1,1]);
            elseif nargin == 6
                set(s, 'Edgecolor', clr,'FaceColor', [1,1,1], 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a pyramid with defined apex coordinates, height, base,
        % orientaion, and color.
        % This method is used to draw translation constraints on 3D models.
        % Input arguments:
        %  X: apex coordinate on the X axis
        %  Y: apex coordinate on the Y axis
        %  Z: apex coordinate on the z axis
        %  h: pyramid height
        %  B: pyramid base
        %  dir: pyramid pointing direction (x+, x-, y+, y-, z+, z-)
        %  c: color (RGB vector)
        function pyramid(X,Y,Z,h,B,dir,c,tag)
            b = B/2;
            
            % Check if direction is char or array of direction vectors
            if ischar(dir)
                % Coodinate convertion
                if strcmp(dir,'x+')
                    x = Y;
                    y = Z;
                    z = X;
                elseif strcmp(dir,'x-')
                    x = Y;
                    y = Z;
                    z = X;
                    h = -h;
                elseif strcmp(dir,'y+')
                    x = Z;
                    y = X;
                    z = Y;
                elseif strcmp(dir,'y-')
                    x = Z;
                    y = X;
                    z = Y;
                    h = -h;
                elseif strcmp(dir,'z+')
                    x = X;
                    y = Y;
                    z = Z;
                elseif strcmp(dir,'z-')
                    x = X;
                    y = Y;
                    z = Z;
                    h = -h;
                end

                % Base coordinates
                Xb = [x + b, x + b, x - b, x - b];
                Yb = [y - b, y + b, y + b, y - b];
                Zb = [z - h, z - h, z - h, z - h];

                % Face 1 coordinates
                Xf1 = [x, x + b, x + b];
                Yf1 = [y, y - b, y + b];
                Zf1 = [z, z - h, z - h];

                % Face 2 coordinates
                Xf2 = [x, x + b, x - b];
                Yf2 = [y, y + b, y + b];
                Zf2 = [z, z - h, z - h];

                % Face 3 coordinates
                Xf3 = [x, x - b, x - b];
                Yf3 = [y, y + b, y - b];
                Zf3 = [z, z - h, z - h];

                % Face 4 coordinates
                Xf4 = [x, x - b, x + b];
                Yf4 = [y, y - b, y - b];
                Zf4 = [z, z - h, z - h];
                
            else % dir is array - compute 5 points coordinates
                
                % Get direction vectors
                dx = dir(1,:);
                dy = dir(2,:);
                dz = dir(3,:);
                
                % Pt1 - apex of the pyramid
                pt_1 = [X, Y, Z];
                
                % Pt2
                pt_2 = pt_1 - h * dx - b * dy - b * dz;
                
                % Pt3
                pt_3 = pt_2 + 2 * b * dz;
                
                % Pt4
                pt_4 = pt_3 + 2 * b * dy;
                
                % Pt5
                pt_5 = pt_4 - 2 * b * dz;
                
                % Base coordinates
                Xb = [pt_2(1), pt_3(1), pt_4(1), pt_5(1)];
                Yb = [pt_2(2), pt_3(2), pt_4(2), pt_5(2)];
                Zb = [pt_2(3), pt_3(3), pt_4(3), pt_5(3)];

                % Face 1 coordinates
                Xf1 = [pt_2(1), pt_3(1), pt_1(1)];
                Yf1 = [pt_2(2), pt_3(2), pt_1(2)];
                Zf1 = [pt_2(3), pt_3(3), pt_1(3)];

                % Face 2 coordinates
                Xf2 = [pt_3(1), pt_4(1), pt_1(1)];
                Yf2 = [pt_3(2), pt_4(2), pt_1(2)];
                Zf2 = [pt_3(3), pt_4(3), pt_1(3)];

                % Face 3 coordinates
                Xf3 = [pt_4(1), pt_5(1), pt_1(1)];
                Yf3 = [pt_4(2), pt_5(2), pt_1(2)];
                Zf3 = [pt_4(3), pt_5(3), pt_1(3)];

                % Face 4 coordinates
                Xf4 = [pt_5(1), pt_2(1), pt_1(1)];
                Yf4 = [pt_5(2), pt_2(2), pt_1(2)];
                Zf4 = [pt_5(3), pt_2(3), pt_1(3)];
                
            end
            
            % Draw pyramid
            if nargin == 7
                if ischar(dir)
                    if strcmp(dir,'z+') || strcmp(dir,'z-')
                        fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0]);

                    elseif strcmp(dir,'x+') || strcmp(dir,'x-')
                        fill3(Zb, Xb, Yb, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Zf1, Xf1, Yf1, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Zf2, Xf2, Yf2, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Zf3, Xf3, Yf3, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Zf4, Xf4, Yf4, c, 'EdgeColor', [0,0,0]);

                    elseif strcmp(dir,'y+') || strcmp(dir,'y-')
                        fill3(Yb, Zb, Xb, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Yf1, Zf1, Xf1, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Yf2, Zf2, Xf2, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Yf3, Zf3, Xf3, c, 'EdgeColor', [0,0,0]);
                        hold on
                        fill3(Yf4, Zf4, Xf4, c, 'EdgeColor', [0,0,0]);
                    end
                else
                    fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0]);
                    hold on
                    fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0]);
                    hold on
                    fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0]);
                    hold on
                    fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0]);
                    hold on
                    fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0]);
                    hold on
                end
            elseif nargin == 8
                if ischar(dir)
                    if strcmp(dir,'z+') || strcmp(dir,'z-')
                        fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0], 'tag', tag);

                    elseif strcmp(dir,'x+') || strcmp(dir,'x-')
                        fill3(Zb, Xb, Yb, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Zf1, Xf1, Yf1, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Zf2, Xf2, Yf2, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Zf3, Xf3, Yf3, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Zf4, Xf4, Yf4, c, 'EdgeColor', [0,0,0], 'tag', tag);

                    elseif strcmp(dir,'y+') || strcmp(dir,'y-')
                        fill3(Yb, Zb, Xb, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Yf1, Zf1, Xf1, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Yf2, Zf2, Xf2, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Yf3, Zf3, Xf3, c, 'EdgeColor', [0,0,0], 'tag', tag);
                        hold on
                        fill3(Yf4, Zf4, Xf4, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    end
                else
                    fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    hold on
                    fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    hold on
                    fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    hold on
                    fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    hold on
                    fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0], 'tag', tag);
                    hold on
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 2D arrow with defined beggining coordinates, length,
        % arrowhead height, arrowhead base, orientation, and color.
        % This method is used to draw load symbols on 2D models.
        % Input arguments:
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  l: arrow length
        %  h: arrowhead height
        %  B: arrowhead base
        %  ang: pointing direction (angle in radian with the horizontal
        %       direction - counterclockwise) - 0 rad when pointing left
        %  c: color (RGB vector)
        function arrow2D(x,y,l,h,B,ang,c,tag)
            b = B/2;
            
            cx = cos(ang);
            cy = sin(ang);
            
            % Draw body line
            X = [x, x + l * cx];
            Y = [y, y + l * cy];
            if nargin == 7
                line(X, Y, 'Color', c);
            elseif nargin == 8
                line(X, Y, 'Color', c, 'tag', tag);
            end
            
            % Draw arrowhead
            X = [x, x + h * cx + b * cy, x + h * cx - b * cy];
            Y = [y, y + h * cy - b * cx, y + h * cy + b * cx];
            if nargin == 7
                fill(X, Y, c);
            elseif nargin == 8
                fill(X, Y, c, 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 3D arrow with defined beggining coordinates, length,
        % arrowhead height, arrowhead base, orientation, and color.
        % This method is used to draw load symbols on grillage models.
        % Input arguments:
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  z: beggining coordinate on the Z axis
        %  l: arrow length
        %  h: arrowhead height
        %  B: arrowhead base
        %  dir: pointing direction (x+, x-, y+, y-, z+, z-)
        %  c: color (RGB vector)
        function draw = arrow3D(draw,x,y,z,l,h,B,dir,c,tag,lineWidthFlag)
            % Compute arrow coordinates
            if ~ischar(dir)
                pt = [x,y,z] - l * dir(1,:);
                X = [x, pt(1)];
                Y = [y, pt(2)];
                Z = [z, pt(3)];
            elseif strcmp(dir,'x+')
                X = [x, x - l];
                Y = [y, y];
                Z = [z, z];      
            elseif strcmp(dir,'x-')
                X = [x, x + l];
                Y = [y, y];
                Z = [z, z];
            elseif strcmp(dir,'y+')
                X = [x, x];
                Y = [y, y - l];
                Z = [z, z];
            elseif strcmp(dir,'y-')
                X = [x, x];
                Y = [y, y + l];
                Z = [z, z];
            elseif strcmp(dir,'z+')
                X = [x, x];
                Y = [y, y];
                Z = [z, z - l];
            elseif strcmp(dir,'z-')
                X = [x, x];
                Y = [y, y];
                Z = [z, z + l];
            end
            
            % Draw arrowhead
            if nargin < 10
                draw.pyramid(x,y,z,h,B,dir,c);
            else
                draw.pyramid(x,y,z,h,B,dir,c,tag);
            end
            
            % Draw body line
            if nargin == 9
                line(X, Y, Z, 'Color', c);
            elseif nargin == 10
                line(X, Y, Z, 'Color', c, 'tag', tag);
            elseif nargin == 11
                if lineWidthFlag
                    lineWidth = 1.2;
                else
                    lineWidth = 1;
                end
                line(X, Y, Z, 'Color', c, 'tag', tag,'linewidth',lineWidth);
            end
        end
        
        %------------------------------------------------------------------
        % Plots a symbol of a 3D force applied to a plane (plane
        % projection of an arrow3D).
        % This method is used to draw load symbols on grillage2D models.
        % Input arguments:
        %  x: coordinate on the X axis
        %  y: coordinate on the Y axis
        %  B: arrowhead base
        %  dir: pointing direction (z+, z-)
        %  c: color (RGB vector)
        function arrow3DPlaneProj(x,y,B,dir,c,tag)
            b = B/2;
            
            % Draw arrowhead base
            xs = [x-b, x-b, x+b, x+b, x-b];
            ys = [y-b, y+b, y+b, y-b, y-b];
            if nargin == 5
                plot(xs,ys,'color',c);
            elseif nargin == 6
                plot(xs,ys,'color',c,'tag',tag);
            end
            
            % Draw symbol to indicate if force is towards or from plane
            switch dir
                case 'z+'
                    r = b/4;
                    circ = 0 : pi/50 : 2*pi;
                    xcirc = x + r * cos(circ);
                    ycirc = y + r * sin(circ);
                    if nargin == 5
                        fill(xcirc, ycirc, c);
                    elseif nargin == 6
                        fill(xcirc, ycirc, c, 'tag', tag);
                    end
                case 'z-'
                    if nargin == 5
                        plot(xs([1 3]),ys([1 3]),'color',c)
                        plot(xs([2 4]),ys([2 4]),'color',c)
                    elseif nargin == 6
                        plot(xs([1 3]),ys([1 3]),'color',c, 'tag', tag)
                        plot(xs([2 4]),ys([2 4]),'color',c, 'tag', tag)
                    end
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 3D arrow in the direction of any element local axis.
        % This method is used to draw load symbols on 3D models.
        % Input arguments:
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  z: beggining coordinate on the Z axis
        %  l: arrow length
        %  h: arrowhead height
        %  B: arrowhead base
        %  dir: pointing direction in local system (x+, x-, y+, y-, z+, z-)
        %  c: color (RGB vector)
        %  e: element ID
        function draw = spear3D(x,y,z,l,h,B,dir,c,e,draw,tag)
            if strcmp(dir,'x-')
                [x2,y2,z2] = draw.coordTransf3D(h, B, B, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(h, -B, B, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(h, -B, -B, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(h, B, -B, x, y, z, e);
            elseif strcmp(dir,'x+')
                [x2,y2,z2] = draw.coordTransf3D(-h, B, B, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(-h, -B, B, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(-h, -B, -B, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(-h, B, -B, x, y, z, e);
            elseif strcmp(dir,'y-')
                [x2,y2,z2] = draw.coordTransf3D(B, h, B, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(B, h, -B, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(-B, h, -B, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(-B, h, B, x, y, z, e);
            elseif strcmp(dir,'y+')
                [x2,y2,z2] = draw.coordTransf3D(B, -h, B, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(B, -h, -B, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(-B, -h, -B, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(-B, -h, B, x, y, z, e);
            elseif strcmp(dir,'z-')
                [x2,y2,z2] = draw.coordTransf3D(B, B, h, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(B, -B, h, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(-B, -B, h, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(-B, B, h, x, y, z, e);
            elseif strcmp(dir,'z+')
                [x2,y2,z2] = draw.coordTransf3D(B, B, -h, x, y, z, e);
                [x3,y3,z3] = draw.coordTransf3D(B, -B, -h, x, y, z, e);
                [x4,y4,z4] = draw.coordTransf3D(-B, -B, -h, x, y, z, e);
                [x5,y5,z5] = draw.coordTransf3D(-B, B, -h, x, y, z, e);
            end
            
            % Base coordinates
            Xb = [x2, x3, x4, x5];
            Yb = [y2, y3, y4, y5];
            Zb = [z2, z3, z4, z5];
            
            % Face 1 coordinates
            Xf1 = [x, x2, x3];
            Yf1 = [y, y2, y3];
            Zf1 = [z, z2, z3];
            
            % Face 2 coordinates
            Xf2 = [x, x3, x4];
            Yf2 = [y, y3, y4];
            Zf2 = [z, z3, z4];
            
            % Face 3 coordinates
            Xf3 = [x, x4, x5];
            Yf3 = [y, y4, y5];
            Zf3 = [z, z4, z5];
            
            % Face 4 coordinates
            Xf4 = [x, x5, x2];
            Yf4 = [y, y5, y2];
            Zf4 = [z, z5, z2];
            
            % Draw pyramid
            if nargin == 10
                fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0]);
                hold on
                fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0]);
                hold on
                fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0]);
                hold on
                fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0]);
                hold on
                fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0]);
            elseif nargin == 11
                fill3(Xb, Yb, Zb, c, 'EdgeColor', [0,0,0], 'tag', tag);
                hold on
                fill3(Xf1, Yf1, Zf1, c, 'EdgeColor', [0,0,0], 'tag', tag);
                hold on
                fill3(Xf2, Yf2, Zf2, c, 'EdgeColor', [0,0,0], 'tag', tag);
                hold on
                fill3(Xf3, Yf3, Zf3, c, 'EdgeColor', [0,0,0], 'tag', tag);
                hold on
                fill3(Xf4, Yf4, Zf4, c, 'EdgeColor', [0,0,0], 'tag', tag);
            end
            
            if strcmp(dir,'x-')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(-l, 0, 0, x, y, z, e);
            elseif strcmp(dir,'x+')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(l, 0, 0, x, y, z, e);
            elseif strcmp(dir,'y-')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(0, l, 0, x, y, z, e);
            elseif strcmp(dir,'y+')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(0, -l, 0, x, y, z, e);
            elseif strcmp(dir,'z-')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(0, 0, l, x, y, z, e);
            elseif strcmp(dir,'z+')
                [xpp,ypp,zpp] = draw.coordTransf3D(0, 0, 0, x, y, z, e);
                [xsp,ysp,zsp] = draw.coordTransf3D(0, 0, -l, x, y, z, e);
            end
            
            % Shaft Coordinates
            X = [xpp,xsp];
            Y = [ypp,ysp];
            Z = [zpp,zsp];
            
            % Draw body line
            if nargin == 10
                line(X, Y, Z, 'Color', [0,0,0]);
            elseif nargin == 11
                line(X, Y, Z, 'Color', [0,0,0], 'tag', tag);
            end
        end
        
        %------------------------------------------------------------------
        % Plots an out-of-plane moment symbol with defined center
        % coordinates, radius, orientation, and color.
        % This method is used to draw applied moment symbols on 2D models.
        % Input arguments:
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: symbol radius
        %  dir: moment orientation (z+, z-)
        %  c: color (RGB vector)
        function draw = moment2D(draw,x,y,r,dir,c,tag)
            h = r/2;
            b = h;
            
            % Draw a semi-circle
            th = -pi/2 : pi/30 : pi/2;
            xp = x + r * cos(th);
            yp = y + r * sin(th);
            if nargin == 6
                plot(xp, yp, 'Color', c);
            elseif nargin == 7
                plot(xp, yp, 'Color', c, 'tag', tag);
            end
            
            % Draw moment orientation
            if nargin == 6
                if strcmp(dir,'z+')
                    draw.triangle(x, y + r, h, b, 2*pi-0.3, c)
                elseif strcmp(dir,'z-')
                    draw.triangle(x, y - r, h, b, 0.3, c)
                end
            elseif nargin == 7
                if strcmp(dir,'z+')
                    draw.triangle(x, y + r, h, b, 2*pi-0.3, c, tag)
                elseif strcmp(dir,'z-')
                    draw.triangle(x, y - r, h, b, 0.3, c, tag)
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 3D moment symbol with defined beggining coordinates,
        % length, arrowhead height, arrowhead base, orientation, and color.
        % This method is used to draw applied moment symbols on grillage
        % and 3D models.
        % Input arguments:
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  z: beggining coordinate on the Z axis
        %  l: total symbol length
        %  h: arrowhead height
        %  B: arrowhead base
        %  dir: moment orientation (x+, x-, y+, y-, z+, z-)
        %  c: color (RGB vector)
        function draw = moment3D(draw,x,y,z,l,h,B,dir,c,tag,srj_flag)
            
            if nargin < 10
                draw.arrow3D(draw,x,y,z,l,h,B,dir,c);
                hold on
            elseif nargin < 11
                draw.arrow3D(draw,x,y,z,l,h,B,dir,c,tag);
                hold on
            else
                draw.arrow3D(draw,x,y,z,l,h,B,dir,c,tag,srj_flag);
                hold on
            end
            
            if ~ischar(dir)
                pt = [x,y,z] - h * dir(1,:);
                x = pt(1);
                y = pt(2);
                z = pt(3);
            elseif strcmp(dir,'x+')
                x = x - h;
            elseif strcmp(dir,'x-')
                x = x + h;
            elseif strcmp(dir,'y+')
                y = y - h;
            elseif strcmp(dir,'y-')
                y = y + h;
            elseif strcmp(dir,'z+')
                z = z - h;
            elseif strcmp(dir,'z-')
                z = z + h;
            end
            
            if nargin < 10
                draw.pyramid(x,y,z,h,B,dir,c);
                hold on
            else
                draw.pyramid(x,y,z,h,B,dir,c,tag);
                hold on
            end
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % computes drawing size parameter according to model dimensions.
        draw = setSize(draw)
        
        %------------------------------------------------------------------
        % Sets axes limits according to model dimensions.
        draw = setLimits(draw)
        
        %------------------------------------------------------------------
        % Draws structure model with nodes, elements, supports and  hinges.
        draw = model(draw)
        
        %------------------------------------------------------------------
        % Draws nodal marks with support conditions.
        draw = nodes(draw)
        
        %------------------------------------------------------------------
        % Draws elements with hinged or continuous ends.
        draw = elements(draw)
        
        %------------------------------------------------------------------
        % Draws semi-rigid joints accordingly to analysis model.
        % Input arguments:
        %  coords: semi-rigid joint coordinates ([x, y, z])
        %  sz: size parameter (2D - Spring Height, 3D - Sphere Radius)
        %  dir: direction to wich semi-rigid joint is applied (local axis)
        %  clr: semi-rigid joint symbol color
        draw = srjoint(draw,coords,sz,dir,clr)
        
        %------------------------------------------------------------------
        % Computes element loads scale factor.
        draw = elemLoadsScaleFactor(draw)
        
        %------------------------------------------------------------------
        % Draws element distributed loads (uniform and linear).
        draw = elemLoads(draw)
        
        %------------------------------------------------------------------
        % Draws applied nodal loads and moments.
        draw = nodalLoads(draw)
        
        %------------------------------------------------------------------
        % Draws applied dynamic nodal loads and moments.
        draw = dynamicNodalLoads(draw) 
        
        %------------------------------------------------------------------
        % Draws applied concentrated nodal mass.
        draw = nodalMass(draw)
        
        %------------------------------------------------------------------
        % Draws nodal prescribed displacement representation.
        draw = nodalPrescDispl(draw)
        
        %------------------------------------------------------------------
        % Draws nodal initial condition values.
        draw = nodalInitialConditions(draw)
        
        %------------------------------------------------------------------
        % Draws thermal load representation on elements.
        draw = thermalLoads(draw)
        
        %------------------------------------------------------------------
        % Plots ID number of nodes.
        draw = nodeID(draw)
        
        %------------------------------------------------------------------
        % Plots ID number of elements.
        draw = elementID(draw)
        
        %------------------------------------------------------------------
        % Draws element orientation indication from inital to final node.
        draw = elementOrientation(draw)
        
        %------------------------------------------------------------------
        % Computes deformed configuration scale factor.
        draw = deformScaleFactor(draw)
        
        %------------------------------------------------------------------
        % Computes dynamic deformed configuration scale factor.
        draw = dynamicDeformScaleFactor(draw)
        
        %------------------------------------------------------------------
        % Draws structure deformed configuration on a given scale.
        % Input arguments:
        %  scale: deformed configuration scale factor
        draw = deformConfig(draw,scale)
        
        %------------------------------------------------------------------
        % Computes axial force diagram scale factor value.
        draw = axialScaleFactor(draw)
        
        %------------------------------------------------------------------
        % Draws resulting axial force diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        draw = axialForce(draw,scale)
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        axialForceEnvelop(draw,scale)
        
        %------------------------------------------------------------------
        % Computes torsion force diagram scale factor value (for envelop).
        draw = torsionScaleFactor(draw)
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment diagram.
        draw = torsionMoment(draw)
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment envelop diagram.
        torsionMomentEnvelop(draw,scale)
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XY plane.
        draw = shearScaleFactor_XY(draw)
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        draw = shearForce_XY(draw,scale)
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XY plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        shearForceEnvelop_XY(draw,scale)
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XZ plane.
        draw = shearScaleFactor_XZ(draw)
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XZ plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        draw = shearForce_XZ(draw,scale)
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XZ plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        shearForceEnvelop_XZ(draw,scale)
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XY plane.
        draw = bendingMomentScaleFactor_XY(draw)
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        draw = bendingMoment_XY(draw,scale)
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XY plane
        % on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        bendingMomentEnvelop_XY(draw,scale)
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XZ plane.
        draw = bendingMomentScaleFactor_XZ(draw)
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XZ plane on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        draw = bendingMoment_XZ(draw,scale)
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XZ plane
        % on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        bendingMomentEnvelop_XZ(draw,scale)
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        draw = reactions(draw)
        
        %------------------------------------------------------------------
        % Draws natural undamped free vibration of specified mode
        % Input arguments:
        %  nMode: vibration mode identifier
        %  scale: deformation scale factor
        vibrationMode(draw,nMode,scale)
        
        %------------------------------------------------------------------
        % Draws deformed configuration after dynamic analysis, on given
        % time.
        % Input arguments:
        %  step : time step identifier
        %  scale: deformation scale factor
        dynamicDeform(draw,step,scale)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Transforms coordinate components from an original axis system
        % (Xi, Yi) to a new axis system (Xf, Yf) subjected to a
        % translation and rotation.
        % Output:
        %  xf: new coordinate on the Xf axis
        %  yf: new coordinate on the Yf axis
        % Input arguments:
        %  xi: original coordinate on the Xi axis
        %  yi: original coordinate on the Yi axis
        %  dx: translation of the X axis
        %  dy: translation of the Y axis
        %  ang: rotation angle (counterclockwise)
        function [xf,yf] = coordTransf2D(~,xi,yi,dx,dy,ang)
            cx = cos(ang);
            cy = sin(ang);
            
            xf = cx * xi - cy * yi + dx;
            yf = cy * xi + cx * yi + dy;
        end
        
        %------------------------------------------------------------------
        % Transforms coordinate components from an original axis system
        % (Xi, Yi, Zi) to a new axis system (Xf, Yf, Zf) subjected to a
        % translation and rotation.
        % Output:
        %  xf: new coordinate on the Xf axis
        %  yf: new coordinate on the Yf axis
        %  zf: new coordinate on the Zf axis
        % Input arguments:
        %  xi: original coordinate on the Xi axis
        %  yi: original coordinate on the Yi axis
        %  zi: original coordinate on the Zi axis
        %  dx: translation of the X axis
        %  dy: translation of the Y axis
        %  dz: translation of the Z axis
        %  e:  element ID number
        function [xf,yf,zf] = coordTransf3D(draw,xi,yi,zi,dx,dy,dz,e)
            rot = draw.mdl.elems(e).T;
            
            old_basis = [xi,yi,zi]';
            translation = [dx,dy,dz]';
            new_basis = rot' * old_basis + translation;
            
            xf = new_basis(1);
            yf = new_basis(2);
            zf = new_basis(3);
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Draw object.
        function draw = clean(draw)
            draw.mdl = [];
            draw.size = 0;
            draw.az = 0;
            draw.elev = 0;
        end
    end
end