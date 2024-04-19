%% Draw_Truss2D class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *2D Truss* draw object.
%
classdef Draw_Truss2D < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Truss2D(mdl)
            draw = draw@Draw(0,90);
            
            if (nargin > 0)
                draw.mdl = mdl;
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <draw.html *Draw*>.
    methods
        %------------------------------------------------------------------
        % Sets viewpoint position and computes drawing size parameter
        % according to model dimensions.
        function draw = setSize(draw)
            x = zeros(draw.mdl.nnp,1);
            y = zeros(draw.mdl.nnp,1);
            for n = 1:draw.mdl.nnp
                x(n) = draw.mdl.nodes(n).coord(1);
                y(n) = draw.mdl.nodes(n).coord(2);
            end
            dx = max(x) - min(x);
            dy = max(y) - min(y);
            sz = max([dx,dy]);
            if isempty(sz) || sz == 0
                draw.size = 5;
            else
                draw.size = sz;
            end
        end
        
        %------------------------------------------------------------------
        % Sets axes limits according to model dimensions.
        function draw = setLimits(draw,axProp)
            if nargin == 2
                setBothAxis = true;
            else
                setBothAxis = false;
            end
            if draw.mdl.nnp > 0
                x = zeros(draw.mdl.nnp,1);
                y = zeros(draw.mdl.nnp,1);
                
                for n = 1:draw.mdl.nnp
                    x(n) = draw.mdl.nodes(n).coord(1);
                    y(n) = draw.mdl.nodes(n).coord(2);
                end
                
                xmin = min(x);
                xmax = max(x);
                ymin = min(y);
                ymax = max(y);
                
                dx = xmax - xmin;
                dy = ymax - ymin;
                
                if (dx == 0) && (dy == 0)
                    xlim([x - 5, x + 5])
                    if setBothAxis == true
                        ylim([y - 5/axProp, y + 5/axProp])
                    end
                elseif dx > dy
                    xlim([xmin - dx/1.7, xmax + dx/1.7])
                    if setBothAxis == true
                        aux_xlim = [xmin - dx/1.7, xmax + dx/1.7];
                        xlength = abs(diff(aux_xlim));
                        ylength = xlength/axProp;
                        ymean = mean([ymin ymax]);
                        ylim([ymean - ylength/2, ymean + ylength/2])
                    end
                else
                    if setBothAxis == true
                        aux_ylim = [ymin - dy/4, ymax + dy/4];
                        ylength = abs(diff(aux_ylim));
                        xlength = ylength * axProp;
                        xmean = mean([xmin xmax]);
                        xlim([xmean - xlength/2, xmean + xlength/2])
                    end
                    ylim([ymin - dy/4, ymax + dy/4])
                end
            else
                xlim([-5,5])
                if setBothAxis == true
                    ylim([-5/axProp, 5/axProp])
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws structure model with nodes, elements, supports and  hinges.
        function draw = model(draw)
            draw.nodes();
            hold on
            draw.elements();
        end
        
        %------------------------------------------------------------------
        % Draws nodal marks with support conditions.
        function draw = nodes(draw)
            % Get flag for support visualization option
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            
            % Parameters
            r = draw.size/125;     % hinge symbol radius
            th = draw.size/35;     % translation constraint symbol (triangle height)
            tb = draw.size/50;     % translation constraint symbol (triangle base)
            sh = draw.size/20;     % displacement spring symbol (spring height)
            nclr = [0,0,0];        % node (hinge) color
            sclr = [0.6,0.6,0.6];  % support color
            sprclr = [0.6,0,0.4];  % spring color
            dc = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                % Get nodal coordinates
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Draw hinge as a nodal point
                draw.circle(x, y, r, nclr,'drawNodes');
                hold on
                
                % Draw support conditions
                if ~strcmp(drawSupports,'on')
                    continue;
                end
                
                % Get rotation
                if draw.mdl.nodes(n).isInclinedSupp
                    dir = draw.mdl.nodes(n).inclSuppDir(1:2);
                    if dir(2) >= 0
                        thetaZ = acos(dir(1));
                    else
                        thetaZ = -acos(dir(1));
                    end
                end
                
                % Draw fixed support in X direction
                if draw.mdl.nodes(n).ebc(1) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.triangle(x-r,y,th,tb,pi,sclr,'drawSupports');
                        hold on
                    else
                        draw.triangle(x-r*cos(thetaZ),y-r*sin(thetaZ),th,tb,pi+thetaZ,sclr,'drawSupports');
                        hold on
                    end
                
                % Draw spring support in X direction
                elseif draw.mdl.nodes(n).ebc(1) == SPRING_DOF
                    kx = draw.mdl.nodes(n).springStiff(1);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if kx >= 1000
                            value = sprintf('%.*e kN/m',dc,kx);
                        else
                            value = sprintf('%.*f kN/m',dc,kx);
                        end
                    else
                        if kx >= 1000
                            value = sprintf('%.*e',dc,kx);
                        else
                            value = sprintf('%.*f',dc,kx);
                        end
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.springX_2D(x-r,y,sh,sprclr,'drawSupports');
                        hold on;
                        text(x-r-sh/2,y+sh/10,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',kx);
                    else
                        draw.springX_2D(x-r*cos(thetaZ),y-r*sin(thetaZ),sh,sprclr,'drawSupports',thetaZ);
                        hold on;
                        text(x-(r+sh)*cos(thetaZ),y-(r+sh)*sin(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',kx);
                    end
                end
                
                % Draw fixed support in Y direction
                if draw.mdl.nodes(n).ebc(2) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.triangle(x,y-r,th,tb,3*pi/2,sclr,'drawSupports');
                        hold on
                    else
                        draw.triangle(x+r*sin(thetaZ),y-r*cos(thetaZ),th,tb,3*pi/2+thetaZ,sclr,'drawSupports');
                        hold on
                    end
                
                % Draw spring support in Y direction
                elseif draw.mdl.nodes(n).ebc(2) == SPRING_DOF
                    ky = draw.mdl.nodes(n).springStiff(2);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if ky >= 1000
                            value = sprintf('%.*e kN/m',dc,ky);
                        else
                            value = sprintf('%.*f kN/m',dc,ky);
                        end
                    else
                        if ky >= 1000
                            value = sprintf('%.*e',dc,ky);
                        else
                            value = sprintf('%.*f',dc,ky);
                        end
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.springY_2D(x,y-r,sh,sprclr,'drawSupports');
                        hold on;
                        text(x+sh/10,y-r-sh/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',ky);
                    else
                        draw.springY_2D(x+r*sin(thetaZ),y-r*cos(thetaZ),sh,sprclr,'drawSupports',thetaZ)
                        hold on;
                        text(x+(r+sh)*sin(thetaZ),y-(r+sh)*cos(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',ky);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws elements with hinged ends.
        function draw = elements(draw)
            % Parameters
            r = draw.size/125;  % hinge symbol radius
            clr = [0,0,0];      % element color
            
            for e = 1:draw.mdl.nel
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with X and Y axes
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Set element end coordinates
                xi = x1 + r * cx;
                yi = y1 + r * cy;
                xf = x2 - r * cx;
                yf = y2 - r * cy;
                
                % Connect element end coordinates
                X = [xi, xf];
                Y = [yi, yf];
                line(X, Y, 'Color', clr, 'tag', 'drawElements');
            end
        end
        
        %------------------------------------------------------------------
        function draw = srjoint(draw,~,~,~,~)
        end
        
        %------------------------------------------------------------------
        % Computes element loads scale factor.
        function scl = elemLoadsScaleFactor(draw)
            max_load = zeros(1,draw.mdl.nel);
            
            for e = 1:draw.mdl.nel
                % Initialize load values on element ends
                qi = 0;
                qf = 0;
                
                % Add uniform load contribution
                if isempty(draw.mdl.elems(e).load.uniformLcl) == 0
                    qy = draw.mdl.elems(e).load.uniformLcl(2);
                    qi = qi + qy;
                    qf = qf + qy;
                end
                
                % Add linear load contribution
                if isempty(draw.mdl.elems(e).load.linearLcl) == 0
                    qyi = draw.mdl.elems(e).load.linearLcl(2);
                    qyf = draw.mdl.elems(e).load.linearLcl(5);
                    qi = qi + qyi;
                    qf = qf + qyf;
                end
                
                % Get maximum and minimum load value on current element
                max_load(e) = max(abs(qi),abs(qf));
            end
            
            % Get maximum and minimum load value on model
            max_val = max(max_load);
            
            % Calculate scale factor
            if isempty(max_val) || max_val == 0
                scl = 0;
            else
                scl = draw.size/(12*max_val);
            end
            
            setappdata(0,'load_sf',scl);
            setappdata(0,'max_val',max_val);
        end
        
        %------------------------------------------------------------------
        % Draws element distributed loads (uniform and linear).
        function draw = elemLoads(draw)
            % Check if elem loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewDistribLoadsButton,'Checked'),'off')
                return
            end
            include_constants;
            
            % Parameters
            r = draw.size/125;   % hinge symbol radius
            ah = draw.size/150;  % load symbol size (arrowhead height)
            ab = draw.size/150;  % load symbol size (arrowhead base)
            clr = [1,0,0];       % load symbol color
            scl = getappdata(0,'load_sf');                 % Load symbol scale
            rmappdata(0,'max_val');
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for e = 1:draw.mdl.nel
                if ((~isempty(draw.mdl.elems(e).load.uniformGbl))  &&...
                   (~all(draw.mdl.elems(e).load.uniformGbl == 0))) ||...
                   ((~isempty(draw.mdl.elems(e).load.linearGbl))   &&...
                   (~all(draw.mdl.elems(e).load.linearGbl == 0)))
               
                    % Get element length
                    L = draw.mdl.elems(e).length;
                    
                    % Get element orientation angle
                    ang = acos(abs(draw.mdl.elems(e).cosine_X));
                    
                    % Get element nodes IDs
                    n1 = draw.mdl.elems(e).nodes(1).id;
                    n2 = draw.mdl.elems(e).nodes(2).id;
                    
                    % Get nodal coordinates
                    x1 = draw.mdl.nodes(n1).coord(1);
                    y1 = draw.mdl.nodes(n1).coord(2);
                    x2 = draw.mdl.nodes(n2).coord(1);
                    y2 = draw.mdl.nodes(n2).coord(2);
                    
                    % Initialize load values on element ends
                    qxi = 0;
                    qyi = 0;
                    qxf = 0;
                    qyf = 0;
                    
                    % Add uniform load contribtuion
                    if isempty(draw.mdl.elems(e).load.uniformLcl) == 0
                        qxi = qxi + draw.mdl.elems(e).load.uniformLcl(1);
                        qyi = qyi + draw.mdl.elems(e).load.uniformLcl(2);
                        qxf = qxi;
                        qyf = qyi;
                    end
                    
                    % Add linear load contribtuion
                    if isempty(draw.mdl.elems(e).load.linearLcl) == 0
                        qxi = qxi + draw.mdl.elems(e).load.linearLcl(1);
                        qyi = qyi + draw.mdl.elems(e).load.linearLcl(2);
                        qxf = qxf + draw.mdl.elems(e).load.linearLcl(4);
                        qyf = qyf + draw.mdl.elems(e).load.linearLcl(5);
                    end
                    
                    % Axial load equation coefficients:
                    % p(x) = Ax + B
                    A = (qxf - qxi)/L;
                    B = qxi;
                    
                    % Transversal load equation coefficients:
                    % q(x) = Cx + D
                    C = (qyf - qyi)/L;
                    D = qyi;
                    
                    % Module of load parameters
                    Qyi = abs(scl * (C * r + D));
                    Qyf = abs(scl * (C * (L-r) + D));
                    
                    % Calculate new element orientation angle and text
                    % rotation angle acording to element orientation
                    if (x1 <= x2) && (y1 <= y2)
                        alpha = ang;
                        rot = ang * 180/pi;
                    elseif (x1 >= x2) && (y1 >= y2)
                        alpha = ang + pi;
                        rot = ang * 180/pi;
                    elseif (x1 <= x2) && (y1 >= y2)
                        alpha = -ang;
                        rot = -ang * 180/pi;
                    elseif (x1 >= x2) && (y1 <= y2)
                        alpha = pi - ang;
                        rot = -ang * 180/pi;
                    end
                    
                    if qyi > 0
                        [xb,yb] = draw.coordTransf2D(r, -Qyi - r, x1, y1, alpha);
                        [xh,yh] = draw.coordTransf2D(0, -(Qyi+r), x1, y1, alpha);
                    elseif qyi < 0
                        [xb,yb] = draw.coordTransf2D(r, Qyi + r, x1, y1, alpha);
                        [xh,yh] = draw.coordTransf2D(0, Qyi+r, x1, y1, alpha);
                    elseif qyi == 0
                        if qyf > 0
                            [xb,yb] = draw.coordTransf2D(r, -Qyi - r, x1, y1, alpha);
                            [xh,yh] = draw.coordTransf2D(0, -(Qyi+r), x1, y1, alpha);
                        elseif qyf < 0
                            [xb,yb] = draw.coordTransf2D(r, Qyi + r, x1, y1, alpha);
                            [xh,yh] = draw.coordTransf2D(0, Qyi + r, x1, y1, alpha);
                        end
                    end
                    
                    if qyf > 0
                        [xc,yc] = draw.coordTransf2D(L-r, -Qyf - r, x1, y1, alpha);
                        [xj,yj] = draw.coordTransf2D(L, -(Qyf+r), x1, y1, alpha);
                    elseif qyf < 0
                        [xc,yc] = draw.coordTransf2D(L-r, Qyf + r, x1, y1, alpha);
                        [xj,yj] = draw.coordTransf2D(L, Qyf+r, x1, y1, alpha);
                    elseif qyf == 0
                        if qyi > 0
                            [xc,yc] = draw.coordTransf2D(L-r, -Qyf - r, x1, y1, alpha);
                            [xj,yj] = draw.coordTransf2D(L, -(Qyf+r), x1, y1, alpha);
                        elseif qyi < 0
                            [xc,yc] = draw.coordTransf2D(L-r, Qyf + r, x1, y1, alpha);
                            [xj,yj] = draw.coordTransf2D(L, Qyf+r, x1, y1, alpha);
                        end
                    end
                    
                    % Draw load symbol on a number cross-sections along element local axis X
                    step = (L-2*r) / round(20 * (L-2*r)/draw.size);
                    for x = r:step:(L-r) 
                        % Calculate load values on current cross-section
                        qx = A * x + B;
                        qy = C * x + D;
                        Qy = abs(scl * qy);
                        
                        % Calculate current cross-section local coordinates
                        [xs,ys] = draw.coordTransf2D(x, 0, x1, y1, alpha);
                        
                        % Coordinates of transversal load symbol Qy
                        if qy > 0
                            [xa,ya] = draw.coordTransf2D(x, -r, x1, y1, alpha);
                            [xd,yd] = draw.coordTransf2D(L/2, -(Qy+r), x1, y1, alpha);
                            if (alpha == ang)
                                dir_y = ang + 3 * pi/2;
                            elseif (alpha == ang + pi)
                                dir_y = ang + pi/2;
                            elseif (alpha == -ang)
                                dir_y = 3 * pi/2 - ang;
                            elseif (alpha == pi - ang)
                                dir_y = pi/2 - ang;
                            end
                        elseif qy < 0
                            [xa,ya] = draw.coordTransf2D(x, r, x1, y1, alpha);
                            [xd,yd] = draw.coordTransf2D(L/2, Qy+r, x1, y1, alpha);
                            if (alpha == ang)
                                dir_y = ang + pi/2;
                            elseif (alpha == ang + pi)
                                dir_y = ang + 3 * pi/2;
                            elseif (alpha == -ang)
                                dir_y = pi/2 - ang;
                            elseif (alpha == pi - ang)
                                dir_y = 3 * pi/2 - ang;
                            end
                        end
                        
                        % Coordinates of axial load symbol Qx
                        if qx > 0
                            if (alpha == ang)
                                dir_x = ang + pi;
                            elseif (alpha == ang + pi)
                                dir_x = ang;
                            elseif (alpha == -ang)
                                dir_x = pi - ang;
                            elseif (alpha == pi - ang)
                                dir_x = 2*pi - ang;
                            end
                        elseif qx < 0
                            if (alpha == ang)
                                dir_x = ang;
                            elseif (alpha == ang + pi)
                                dir_x = ang + pi;
                            elseif (alpha == -ang)
                                dir_x = 2*pi - ang;
                            elseif (alpha == pi - ang)
                                dir_x = pi - ang;
                            end
                        end
                        
                        % Draw axial load symbols
                        if ((x ~= r) && ((x+step) < L)) && ((qxi ~= 0) || (qxf ~= 0))
                            draw.triangle(xs, ys, 2*ah, 2*ab, dir_x, clr,'drawElemLoads');
                        end
                                                
                        % Draw transversal load symbol
                        if Qy >= ah
                            draw.arrow2D(xa, ya, Qy, ah, ab, dir_y, clr,'drawElemLoads');
                        elseif (abs(x-r) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [x_ini,y_ini] = draw.coordTransf2D(r, 0, x1, y1, alpha); 
                            X = [xb, x_ini];
                            Y = [yb, y_ini];
                            line(X, Y, 'Color', clr,'tag','drawElemLoads');
                        elseif (abs(x-L+r) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [x_end,y_end] = draw.coordTransf2D(L-r, 0, x1, y1, alpha);
                            X = [x_end, xc];
                            Y = [y_end, yc];
                            line(X, Y, 'Color', clr,'tag','drawElemLoads');
                        end
                    end
                    
                    % Connect diagram extremities
                    if  (qyi < 0) && (qyf > 0)
                        x0 = abs((Qyi*(L-2*r))/(Qyf+Qyi));
                        [xu,yu] = draw.coordTransf2D(x0, r, x1, y1, alpha);
                        [xl,yl] = draw.coordTransf2D(x0, -r, x1, y1, alpha);
                        X = [xb, xu, xl, xc];
                        Y = [yb, yu, yl, yc];
                        line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                    elseif (qyi > 0) && (qyf < 0)
                        x0 = abs((Qyi*(L-2*r))/(Qyf+Qyi));
                        [xu,yu] = draw.coordTransf2D(x0, r, x1, y1, alpha);
                        [xl,yl] = draw.coordTransf2D(x0, -r, x1, y1, alpha);
                        X = [xb, xl, xu, xc];
                        Y = [yb, yl, yu, yc];
                        line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                    elseif (qyi ~= 0) || (qyf ~= 0)
                        X = [xb, xc];
                        Y = [yb, yc];
                        line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                    end
                    
                    % Write load values:
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value_x1 = sprintf('%.*f kN/m',dc,abs(qxi));
                        value_x2 = sprintf('%.*f kN/m',dc,abs(qxf));
                        value_y1 = sprintf('%.*f kN/m',dc,abs(qyi));
                        value_y2 = sprintf('%.*f kN/m',dc,abs(qyf));
                    else
                        value_x1 = sprintf('%.*f',dc,abs(qxi));
                        value_x2 = sprintf('%.*f',dc,abs(qxf));
                        value_y1 = sprintf('%.*f',dc,abs(qyi));
                        value_y2 = sprintf('%.*f',dc,abs(qyf));
                    end
                    
                    % Write axial load values
                    if (qxi ~= 0) || (qxf ~= 0)
                        if qxi == qxf
                            [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,alpha);
                            text(xm,ym,value_x1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qxi));
                        else
                            [xm,ym] = draw.coordTransf2D(0,0,x1,y1,alpha);
                            [xn,yn] = draw.coordTransf2D(0,0,x2,y2,alpha);
                            if abs(alpha) <= pi/2
                                text(xm,ym,value_x1,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qxi));
                                text(xn,yn,value_x2,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qxf));
                            else
                                text(xm,ym,value_x1,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qxi));
                                text(xn,yn,value_x2,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qxf));
                            end
                        end
                    end
                    
                    % Write transversal load values
                    if (qyi ~= 0) || (qyf ~= 0)
                        if qyi == qyf
                            if abs(alpha) <= pi/2
                                if qyi > 0
                                    text(xd,yd,value_y1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                else
                                    text(xd,yd,value_y1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                end
                            else
                                if qyi > 0
                                    text(xd,yd,value_y1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                else
                                    text(xd,yd,value_y1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                end
                            end
                        else
                            if qyi ~= 0
                                if abs(alpha) <= pi/2
                                    if qyi > 0
                                        text(xh,yh,value_y1,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                    else
                                        text(xh,yh,value_y1,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                    end
                                else
                                    if qyi > 0
                                        text(xh,yh,value_y1,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                    else
                                        text(xh,yh,value_y1,'HorizontalAlignment','right','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyi));
                                    end
                                end
                            end
                            if qyf ~= 0
                                if abs(alpha) <= pi/2
                                    if qyf > 0
                                        text(xj,yj,value_y2,'HorizontalAlignment','right','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyf));
                                    else
                                        text(xj,yj,value_y2,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyf));
                                    end
                                else
                                    if qyf > 0
                                        text(xj,yj,value_y2,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyf));
                                    else
                                        text(xj,yj,value_y2,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',abs(qyf));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws applied nodal loads.
        function draw = nodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            r = draw.size/125;   % hinge symbol radius
            al = draw.size/12;   % load symbol size (arrow length)
            ah = draw.size/60;   % load symbol size (arrowhead height)
            ab = draw.size/60;   % load symbol size (arrowhead base)
            clr = [1,0,0];       % load symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.static(1);
                    fy = draw.mdl.nodes(n).load.static(2);
                    
                    % Draw horizontal load component
                    if fx > 0
                        draw.arrow2D(x - r, y, al, ah, ab, pi, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-r-(ah+al)/3,y,value,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow2D(x + r, y, al, ah, ab, 0, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+r+(ah+al)/3,y,value,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw vertical load component
                    if fy > 0
                        draw.arrow2D(x, y - r, al, ah, ab, 3 * pi/2, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y-r-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow2D(x, y + r, al, ah, ab, pi/2, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y+r+(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws applied dynamic nodal loads.
        function draw = dynamicNodalLoads(draw)
            % Check if nodal loads visualization is on and analysis is not modal
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off') ||...
               draw.mdl.drv.whichResponse ~= 1
                return
            end
            
            % Parameters
            r = draw.size/125;   % hinge symbol radius
            m = draw.size/60;    % concentrated mass symbol radius
            al = draw.size/12;   % load symbol size (arrow length)
            ah = draw.size/60;   % load symbol size (arrowhead height)
            ab = draw.size/60;   % load symbol size (arrowhead base)
            clr = [0,0.7,0];       % load symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.dynamic) == 0
                    if isempty(draw.mdl.nodes(n).displMass) == 0 && draw.mdl.nodes(n).displMass > 0
                        R = m;
                    else
                        R = r;
                    end
                    
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.dynamic(1);
                    fy = draw.mdl.nodes(n).load.dynamic(2);
                    
                    % Draw horizontal load component
                    if fx > 0
                        draw.arrow2D(x - R, y, al, ah, ab, pi, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-R-(ah+al)/3,y,value,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow2D(x + R, y, al, ah, ab, 0, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+R+(ah+al)/3,y,value,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw vertical load component
                    if fy > 0
                        draw.arrow2D(x, y - R, al, ah, ab, 3 * pi/2, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y-R-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow2D(x, y + R, al, ah, ab, pi/2, clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y+R+(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws concentrated nodal mass.
        function draw = nodalMass(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalMassButton,'Checked'),'off')
                return
            end
            
            % Parameters
            r = draw.size/60;               % concentrated mass symbol radius
            clr = [0,0,1];                  % mass color
            dc = getappdata(0,'decPrec');   % decimal precision
            
             for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).displMass)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal concenctrated mass value
                    mass = draw.mdl.nodes(n).displMass;
                    
                    % Draw mass
                    if mass > 0
                        s = draw.sphere(x, y, 0, r, 'drawNodalMass');
                        set(s,'Edgecolor',clr,'FaceColor',clr);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kg',dc,abs(mass)*1000);
                        else
                            value = sprintf('%.*f',dc,abs(mass)*1000);
                        end
                        text(x-r,y,value,'HorizontalAlignment','right','VerticalAlignment','top','Color',clr,'tag','textNodalMass','UserData',abs(mass)*1000);
                    end
                end
             end
        end
        
        %------------------------------------------------------------------
        % Draws nodal prescribed displacement representation.
        function draw = nodalPrescDispl(draw)
            % Check if presc displ visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewPrescDisplButton,'Checked'),'off')
                return
            end
            
            % Parameters
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            shift = 0;              % distance between presc. displ. symbol and support
            r = draw.size/125;      % hinge symbol radius
            al = draw.size/12;      % presc. displ. symbol size (arrow length)
            ah = draw.size/60;      % presc. displ. symbol size (arrowhead height)
            ab = draw.size/60;      % presc. displ. symbol size (arrowhead base)
            clr = [1,0,1];          % prec. displ. symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (triangle height)
            if strcmp(drawSupports,'on')
                th = draw.size/35;
            else
                th = 0;
            end
            
            for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).prescDispl)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get prescribed displacement component values and convert it to milimeter
                    dx = 1000 * draw.mdl.nodes(n).prescDispl(1);
                    dy = 1000 * draw.mdl.nodes(n).prescDispl(2);
                    
                    % Get rotation
                    if draw.mdl.nodes(n).isInclinedSupp
                        s_rot = (th + shift + r);
                        dir = draw.mdl.nodes(n).inclSuppDir(1:2);
                        if dir(2) >= 0
                            thetaZ = acos(dir(1));
                        else
                            thetaZ = - acos(dir(1));
                        end
                    end
                    
                    % Check if horizontal translation is really fixed and draw prescribed displacement indication
                    if (draw.mdl.nodes(n).ebc(1) == 1) && (dx ~= 0)
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dx));
                        else
                            value = sprintf('%.*f',dc,abs(dx));
                        end
                        if ~draw.mdl.nodes(n).isInclinedSupp
                            if dx > 0
                                draw.arrow2D(x-r-th,y,al,ah,ab,pi,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x-r-th-al,y,al,ah,ab,0,clr,'drawPrescDispl');
                            end
                            text(x-r-th-(al+ah)/2,y+al/20,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                        else
                            if dx > 0
                                draw.arrow2D(x-(s_rot)*cos(thetaZ),y-s_rot*sin(thetaZ),al,ah,ab,pi+thetaZ,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x-(s_rot+al)*cos(thetaZ),y-(s_rot+al)*sin(thetaZ),al,ah,ab,thetaZ,clr,'drawPrescDispl');
                            end
                            text(x-(s_rot+al)*cos(thetaZ),y-(s_rot+al)*sin(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                        end
                    end
                    
                    % Check if vertical translation is really fixed and draw prescribed displacement indication
                    if (draw.mdl.nodes(n).ebc(2) == 1) && (dy ~= 0)
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dy));
                        else
                            value = sprintf('%.*f',dc,abs(dy));
                        end
                        if ~draw.mdl.nodes(n).isInclinedSupp
                            if dy > 0
                                draw.arrow2D(x,y-r-th,al,ah,ab,-pi/2,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x,y-r-th-al,al,ah,ab,pi/2,clr,'drawPrescDispl');
                            end
                            text(x+al/25,y-r-th-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color', clr,'tag','textPrescDispl','UserData',abs(dy));
                        else
                            if dy > 0
                                draw.arrow2D(x+s_rot*sin(thetaZ),y-(s_rot)*cos(thetaZ),al,ah,ab,-pi/2+thetaZ,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x+(s_rot+al)*sin(thetaZ),y-(s_rot+al)*cos(thetaZ),al,ah,ab,pi/2+thetaZ,clr,'drawPrescDispl');
                            end
                            text(x+(s_rot+al)*sin(thetaZ),y-(s_rot+al)*cos(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws nodal initial condition values.
        function draw = nodalInitialConditions(draw)
            % Check if initial condition visualization is on and analysis is transient
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewInitialConditionsButton,'Checked'),'off') ||...
               draw.mdl.drv.whichResponse ~= 1
                return
            end
            
            % Parameters
            r   = draw.size/125;  % hinge symbol radius
            m   = draw.size/60;   % concentrated mass symbol radius
            clr = [0,0.7,0];      % color
            
            for n = 1:draw.mdl.nnp
                flag = 0;
                str = '';
                if isempty(draw.mdl.nodes(n).initCond)
                    continue;
                end
                if draw.mdl.nodes(n).initCond(1,1) ~= 0 || draw.mdl.nodes(n).initCond(2,1) ~= 0
                    flag = flag + 1;
                end
                if draw.mdl.nodes(n).initCond(1,2) ~= 0 || draw.mdl.nodes(n).initCond(2,2) ~= 0
                    flag = flag + 2;
                end
                if flag == 0
                    continue
                elseif flag == 1
                    str = 'd0';
                elseif flag == 2
                    str = 'v0';
                elseif flag == 3
                    str = 'd0,v0';
                end
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                if ~isempty(draw.mdl.nodes(n).displMass) &&...
                    draw.mdl.nodes(n).displMass > 0      &&...
                    strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x+m,y-0.3*r,str,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'tag','textInitialConditions');
                else
                    text(x+1.1*r,y-0.3*r,str,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'tag','textInitialConditions');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws thermal load representation on elements.
        function draw = thermalLoads(draw)
            % Check if thermal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewThermalLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            d = draw.size/150;     % distance between temperature grad symbol and element
            heatClr = [1,0,0];     % heat color
            coldClr = [0,0,1];     % cold color
            dc = getappdata(0,'decPrec');  % decimal precision
            
            for e = 1:draw.mdl.nel
                if (draw.mdl.elems(e).load.tempVar_X ~= 0) || ...
                   (draw.mdl.elems(e).load.tempVar_Y ~= 0) || ...
                   (draw.mdl.elems(e).load.tempVar_Z ~= 0)
                    % Get nodal coordinates
                    x1 = draw.mdl.elems(e).nodes(1).coord(1);
                    y1 = draw.mdl.elems(e).nodes(1).coord(2);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    
                    % Get element inclination
                    cx = draw.mdl.elems(e).cosine_X;
                    cy = draw.mdl.elems(e).cosine_Y;
                    
                    % Get temperature gradient values
                    dtx = draw.mdl.elems(e).load.tempVar_X;
                    dty = draw.mdl.elems(e).load.tempVar_Y;
                    
                    % Calculate element orientation angle and text rotation angle acording to element orientation
                    ang = acos(abs(draw.mdl.elems(e).cosine_X));
                    if (x1 <= x2) && (y1 <= y2)
                        rot = ang * 180/pi;
                    elseif (x1 >= x2) && (y1 >= y2)
                        rot = ang * 180/pi;
                    elseif (x1 <= x2) && (y1 >= y2)
                        rot = -ang * 180/pi;
                    elseif (x1 >= x2) && (y1 <= y2)
                        rot = -ang * 180/pi;
                    end
                    
                    % Check text position
                    dtx_pos = 'top';
                    if (dtx ~= 0) && (dty ~= 0)
                        dty_pos = 'bottom';
                    else
                        dty_pos = 'top';
                    end
                    
                    % Check if units are enabled
                    unitsAreOn = strcmp(get(mdata.unitsButton,'Checked'),'on');
                    
                    % Draw temperature variation symbols
                    if dtx > 0
                        line([x1,x2],[y1,y2],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTx = %.*f oC',dc,dtx),'HorizontalAlignment','center','VerticalAlignment',dtx_pos,'Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTx = %.*f',dc,dtx),'HorizontalAlignment','center','VerticalAlignment',dtx_pos,'Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        end
                    elseif dtx < 0
                        line([x1,x2],[y1,y2],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTx = %.*f oC',dc,dtx),'HorizontalAlignment','center','VerticalAlignment',dtx_pos,'Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTx = %.*f',dc,dtx),'HorizontalAlignment','center','VerticalAlignment',dtx_pos,'Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        end
                    end
                    
                    if dty > 0
                        line([x1-d*cy,x2-d*cy],[y1+d*cx,y2+d*cx],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        line([x1+d*cy,x2+d*cy],[y1-d*cx,y2-d*cx],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTy = %.*f oC',dc,dty),'HorizontalAlignment','center','VerticalAlignment',dty_pos,'Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTy = %.*f',dc,dty),'HorizontalAlignment','center','VerticalAlignment',dty_pos,'Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        end
                    elseif dty < 0
                        line([x1-d*cy,x2-d*cy],[y1+d*cx,y2+d*cx],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        line([x1+d*cy,x2+d*cy],[y1-d*cx,y2-d*cx],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTy = %.*f oC',dc,dty),'HorizontalAlignment','center','VerticalAlignment',dty_pos,'Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTy = %.*f',dc,dty),'HorizontalAlignment','center','VerticalAlignment',dty_pos,'Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots ID number of nodes.
        function draw = nodeID(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            anl = get(mdata.popupmenu_AnalysisType,'Value');
            r = draw.size/125;   % hinge symbol radius
            m = draw.size/60;    % concentrated mass symbol radius
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                id = sprintf('%d',n);
                
                if anl == 2                             &&...
                  ~isempty(draw.mdl.nodes(n).displMass) &&...
                   draw.mdl.nodes(n).displMass > 0      &&...
                   strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x+0.75*m,y+0.75*m,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                else
                    text(x+r,y+r,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots ID number of elements.
        function draw = elementID(draw)
            for e = 1:draw.mdl.nel
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                id = sprintf('%d',e);
                text((x1+x2)/2,(y1+y2)/2,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textElemID');
            end
        end
        
        %------------------------------------------------------------------
        % Draws element orientation indication from inital to final node.
        function draw = elementOrientation(draw)
            clr = [0,0.7,0]; % orientation symbol color
            for e = 1:draw.mdl.nel
                % Calculate spear length
                l = draw.size/20;
                
                % Get nodal coordinates
                xi = draw.mdl.elems(e).nodes(1).coord(1);
                yi = draw.mdl.elems(e).nodes(1).coord(2);
                xf = draw.mdl.elems(e).nodes(2).coord(1);
                yf = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Calculate element local axis X orientation vector
                x = [xf-xi, yf-yi, 0];
                x = l * x / norm(x);
                
                % Calculate element local axis Y orientation vector
                y = cross([0,0,1],x);
                y = l * y / norm(y);
                
                % Draw orientation symbol
                xm = (xi + xf)/2;
                ym = (yi + yf)/2;
                
                X = [xm, xm + x(1)];
                Y = [ym, ym + x(2)];
                line(X,Y,'Color',clr,'Linewidth',1.2,'tag','drawElemOrient');
                text(xm+x(1),ym+x(2),'X','HorizontalAlignment','left','VerticalAlignment','baseline','Color',clr,'FontWeight','bold','tag','drawElemOrient');
                
                X = [xm, xm + y(1)];
                Y = [ym, ym + y(2)];
                line(X,Y,'Color',clr,'Linewidth',1.2,'tag','drawElemOrient');
                text(xm+y(1),ym+y(2),'Y','HorizontalAlignment','left','VerticalAlignment','baseline','Color',clr,'FontWeight','bold','tag','drawElemOrient');
            end
        end
        
        %------------------------------------------------------------------
        % Computes deformed configuration scale factor.
        function draw = deformScaleFactor(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
            % Calculate maximum nodal displacement
            m = zeros(1,draw.mdl.nnp);
            for n = 1:draw.mdl.nnp
                dx = draw.mdl.D(draw.mdl.ID(1,n));
                dy = draw.mdl.D(draw.mdl.ID(2,n));
                d = [dx,dy];
                m(n) = max(abs(d));
            end
            max_ndisp = max(m);
            
            % Set adapted scale value
            if isempty(max_ndisp) || max_ndisp == 0
                dsf = 0;
            else
                dsf = draw.size/(4*slider*max_ndisp);
            end
            setappdata(0,'deform_sf',dsf);
        end
        
        %------------------------------------------------------------------
        % Computes dynamic deformed configuration scale factor.
        function draw = dynamicDeformScaleFactor(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Estimate maximum displacement of each element
            % (nodal and internal)
            m = zeros(1,draw.mdl.nnp);
            for n = 1:draw.mdl.nnp
                % Get maximum nodal displacements
                ids = [ draw.mdl.ID(1,n) ;
                        draw.mdl.ID(2,n) ];
                       
                nodeDispl = sum(draw.mdl.results.dynamicDispl(ids,:,:),3);
                m(n)      = max(max(abs(nodeDispl)));
            end
            
            % Get maximum estimated model displacement
            max_disp = max(m);
            
            % Set adapted scale value
            if isempty(max_disp) || max_disp == 0
                dsf = 0;
            else
                dsf = draw.size/(4*sliderm*max_disp);
            end
            setappdata(0,'deform_sf',dsf);
        end
        
        %------------------------------------------------------------------
        % Draws structure deformed configuration on a given scale.
        % Input arguments:
        %  scale: deformed configuration scale factor
        function draw = deformConfig(draw,scale)
            % Parameters
            clr = [1,0,0];  % deformed configuration line color
            
            for e = 1:draw.mdl.nel
                % Get element end nodes IDs
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                
                % Get nodal coordinates
                x1 = draw.mdl.nodes(n1).coord(1);
                y1 = draw.mdl.nodes(n1).coord(2);
                x2 = draw.mdl.nodes(n2).coord(1);
                y2 = draw.mdl.nodes(n2).coord(2);
                
                % Get nodal displacements
                dx1 = draw.mdl.D(draw.mdl.ID(1,n1));
                dy1 = draw.mdl.D(draw.mdl.ID(2,n1));
                dx2 = draw.mdl.D(draw.mdl.ID(1,n2));
                dy2 = draw.mdl.D(draw.mdl.ID(2,n2));
                
                % Calculate displaced nodal coordinates
                xd1 = x1 + scale * dx1;
                yd1 = y1 + scale * dy1;
                xd2 = x2 + scale * dx2;
                yd2 = y2 + scale * dy2;
                
                % Connect displaced nodal coordinates
                X = [xd1, xd2];
                Y = [yd1, yd2];
                line(X, Y, 'Color', clr, 'tag', 'drawDeformConfig');
            end
        end
        
        %------------------------------------------------------------------
        % Computes axial force diagram scale factor value.
        function draw = axialScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function draw = axialForce(draw,~)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element internal axial force value (always uniform
                % for 2D truss models) and convert it to string
                N = -draw.mdl.elems(e).axial_force(1);
                if N ~= 0
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,N);
                    else
                        value = sprintf('%+.*f',dc,N);
                    end
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%.*f kN',dc,abs(N));
                    else
                        value = sprintf('%.*f',dc,abs(N));
                    end
                    N = abs(N);
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element inclination and element length
                alpha = acos(abs(draw.mdl.elems(e).cosine_X));
                ang = alpha*180/pi;
                L = draw.mdl.elems(e).length;
                
                % Write axial force value
                if ((x1 <= x2) && (y1 <= y2))
                    [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,alpha);
                    text(xm,ym,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',N);
                elseif ((x1 >= x2) && (y1 >= y2))
                    [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,alpha+pi);
                    text(xm,ym,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',N);
                elseif ((x1 <= x2) && (y1 >= y2))
                    [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,-alpha);
                    text(xm,ym,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',-ang,'tag','textAxialForceDiagram','UserData',N);
                elseif ((x1 >= x2) && (y1 <= y2))
                    [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,pi-alpha);
                    text(xm,ym,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',-ang,'tag','textAxialForceDiagram','UserData',N);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        function axialForceEnvelop(draw,~)
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                
                % Calculate text coordinates according to elem orientation
                if (x1 <= x2) && (y1 <= y2)
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    ang = -180*(acos(cx))/pi;
                else
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Get internal forces envelop values
                Nmax = draw.mdl.elems(e).intForcesEnvelop(1,:,1);
                Nmin = draw.mdl.elems(e).intForcesEnvelop(2,:,1);
                
                % Plot text containing values
                if strcmp(get(mdata.unitsButton,'Checked'),'on')
                    value_max = sprintf('%+.*f kN',dc,Nmax(1));
                else
                    value_max = sprintf('%+.*f',dc,Nmax(1));
                end
                if strcmp(get(mdata.unitsButton,'Checked'),'on')
                    value_min = sprintf('%+.*f kN',dc,Nmin(1));
                else
                    value_min = sprintf('%+.*f',dc,Nmin(1));
                end
                text((x1+x2)/2,(y1+y2)/2,value_max,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmax(1));
                text((x1+x2)/2,(y1+y2)/2,value_min,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmin(1));
            end
        end
        
        %------------------------------------------------------------------
        % Computes torsion force diagram scale factor value (used for
        % envelop).
        function draw = torsionScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment diagram.
        function draw = torsionMoment(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment envelop diagram.
        function torsionMomentEnvelop(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XY plane.
        function draw = shearScaleFactor_XY(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XY plane on a given scale.
        function draw = shearForce_XY(draw,~)
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XY plane on a
        % given scale.
        function shearForceEnvelop_XY(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XZ plane.
        function draw = shearScaleFactor_XZ(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XZ plane on a given scale.
        function draw = shearForce_XZ(draw,~)
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XY plane on a
        % given scale.
        function shearForceEnvelop_XZ(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XY plane.
        function draw = bendingMomentScaleFactor_XY(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XY plane on a given scale.
        function draw = bendingMoment_XY(draw,~)
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XY plane on a
        % given scale.
        function bendingMomentEnvelop_XY(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XZ plane.
        function draw = bendingMomentScaleFactor_XZ(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XZ plane on a given scale.
        function draw = bendingMoment_XZ(draw,~)
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XZ plane on a
        % given scale.
        function bendingMomentEnvelop_XZ(~,~)
        end
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        function draw = reactions(draw)
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            r = draw.size/125;      % hinge symbol radius
            al = draw.size/12;      % reaction symbol size (arrow length)
            ah = draw.size/60;      % reaction symbol size (arrowhead height)
            ab = draw.size/60;      % reaction symbol size (arrowhead base)
            clr = [0,0,1];          % reaction symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (triangle/spring height)
            if strcmp(drawSupports,'on')
                th = draw.size/35;
                sh = draw.size/20;
            else
                th = 0;
                sh = 0;
            end
            
            for n = 1:draw.mdl.nnp
                % Get nodal coordinates
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Get reactions values
                rx = draw.mdl.F(draw.mdl.ID(1,n));
                ry = draw.mdl.F(draw.mdl.ID(2,n));
                
                % Support (or spring) height in the direction of axis X and Y
                if draw.mdl.nodes(n).ebc(1)== FIXED_DOF
                    hx = r+th;
                elseif draw.mdl.nodes(n).ebc(1)== SPRING_DOF
                    hx = r+sh;
                end
                if draw.mdl.nodes(n).ebc(2)== FIXED_DOF
                    hy = r+th;
                elseif draw.mdl.nodes(n).ebc(2)== SPRING_DOF
                    hy = r+sh;
                end 
                
                % Get rotation
                if draw.mdl.nodes(n).isInclinedSupp
                    dir = draw.mdl.nodes(n).inclSuppDir(1:2);
                    if dir(2) >= 0
                        thetaZ = acos(dir(1));
                    else
                        thetaZ = -acos(dir(1));
                    end
                end
                
                % Check if horizontal translation is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(1)== FIXED_DOF) || (draw.mdl.nodes(n).ebc(1)== SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(rx));
                    else
                        value = sprintf('%.*f',dc,abs(rx));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if rx >= 0
                            draw.arrow2D(x-hx,y,al,ah,ab,pi,clr,'drawReactions');
                        else
                            draw.arrow2D(x-hx-al,y,al,ah,ab,0,clr,'drawReactions');
                        end
                        text(x-hx-(al+ah)/2,y-al/20,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                    else
                        if rx >= 0
                            draw.arrow2D(x-hx*cos(thetaZ),y-hx*sin(thetaZ),al,ah,ab,pi+thetaZ,clr,'drawReactions');
                        else
                            draw.arrow2D(x-(hx+al)*cos(thetaZ),y-(hx+al)*sin(thetaZ),al,ah,ab,thetaZ,clr,'drawReactions');
                        end
                        text(x-(hx+al)*cos(thetaZ),y-(hx+al)*sin(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                    end
                end
                
                % Check if vertical translation is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(2) == FIXED_DOF) || (draw.mdl.nodes(n).ebc(2) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(ry));
                    else
                        value = sprintf('%.*f',dc,abs(ry));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if ry >= 0
                            draw.arrow2D(x,y-hy,al,ah,ab,-pi/2,clr,'drawReactions');
                        else
                            draw.arrow2D(x,y-(hy+al),al,ah,ab,pi/2,clr,'drawReactions');
                        end
                        text(x-al/25,y-hy-(ah+al)/2,value,'HorizontalAlignment','right','VerticalAlignment','middle','Color', clr,'tag','textForceReactions','UserData',abs(ry));
                    else
                        if ry >= 0
                            draw.arrow2D(x+hy*sin(thetaZ),y-(hy)*cos(thetaZ),al,ah,ab,-pi/2+thetaZ,clr,'drawReactions');
                        else
                            draw.arrow2D(x+(hy+al)*sin(thetaZ),y-(hy+al)*cos(thetaZ),al,ah,ab,pi/2+thetaZ,clr,'drawReactions');
                        end
                        text(x+(hy+al)*sin(thetaZ),y-(hy+al)*cos(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws natural undamped free vibration of specified mode
        % Input arguments:
        %  nMode: vibration mode identifier
        %  scale: deformation scale factor
        function vibrationMode(draw,nMode,scale)
            %  Parameters
            clr = [1,0,0];  % deformed configuration line color
            
            % Initialize plotting matrix
            d = zeros(2,draw.mdl.nel*3);
            
            % Loop over elems to concat nodal displacements
            for e = 1:draw.mdl.nel
                % Get end point coordinates
                coords = draw.mdl.elems(e).intCoords(:,[1,end]);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = draw.mdl.elems(e).natVibration(:,[1,end],nMode);
                
                % Get element orientation cosines
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Assemble rotation transformation matrix
                rot = [ cx  cy;
                       -cy  cx ];
                    
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords(1:2,:) + scale * dg;
                
                % Concatenate to displ mtx
                d(:,(e-1)*3+1:e*3) = [dfg nan(2,1)];
            end
            % Plot deformed configuration
            plot(d(1,:), d(2,:), 'Color', clr, 'tag', 'drawVibrationMode');
            drawnow
        end
        
        %------------------------------------------------------------------
        % Draws deformed configuration after dynamic analysis, on given
        % time.
        % Input arguments:
        %  step : time step identifier
        %  scale: deformation scale factor
        function dynamicDeform(draw,step,scale)
            %  Parameters
            clr = [1,0,0];  % deformed configuration line color
            
            % Check if step is not an integer
            dt = rem(step,1);
            if step >= draw.mdl.n_steps + 1
                step = draw.mdl.n_steps;
                dt = 1;
            end
            
            % Initialize plotting matrix
            d = zeros(2,draw.mdl.nel*3);
            
            % Loop over elems to concat nodal displacements
            for e = 1:draw.mdl.nel
                % Get end point coordinates
                coords = draw.mdl.elems(e).intCoords(:,[1,end]);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = (1-dt) * draw.mdl.elems(e).dynamicIntDispl(:,[1,end],floor(step)) +...
                        dt  * draw.mdl.elems(e).dynamicIntDispl(:,[1,end],floor(step)+1);
                
                % Get element orientation cosines
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Assemble rotation transformation matrix
                rot = [ cx  cy;
                       -cy  cx ];
                
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords(1:2,:) + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*3+1:e*3) = [dfg nan(2,1)];
            end
            % Plot deformed configuration
            plot(d(1,:), d(2,:), 'Color', clr, 'tag', 'drawDynamicDeform');
            drawnow
        end
    end
end