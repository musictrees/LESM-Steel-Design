%% Draw_Frame2D class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *2D Frame* draw object.
%
classdef Draw_Frame2D < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Frame2D(mdl)
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
            hold on;
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
            nm = draw.size/200;    % node mark symbol (square side)
            r = draw.size/125;     % hinge symbol radius
            th = draw.size/35;     % translation constraint symbol (triangle height)
            tb = draw.size/50;     % translation constraint symbol (triangle base)
            ss = draw.size/35;     % rotation constraint symbol (square side)
            sh = draw.size/20;     % displacement spring symbol (spring height)
            nclr = [0,0,0];        % node and hinge color
            sclr = [0.6,0.6,0.6];  % support color
            sprclr = [0.6,0,0.4];  % spring color
            dc = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Distance between translation constraint support symbol and nodal point
                shift = 0;
                
                if draw.mdl.nodes(n).ebc(6) == FREE_DOF
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng == tot && tot > 0
                        shift = r;
                    else
                        draw.square(x,y,nm,nclr,'drawNodes');
                        hold on;
                    end
                elseif draw.mdl.nodes(n).ebc(6) == FIXED_DOF
                    if strcmp(drawSupports,'on')
                        shift = ss/2;
                    else
                        draw.square(x,y,nm,nclr,'drawNodes');
                        hold on;
                    end
                elseif draw.mdl.nodes(n).ebc(6) == SPRING_DOF
                    if ~strcmp(drawSupports,'on')
                        draw.square(x,y,nm,nclr,'drawNodes');
                        hold on;
                    end
                end
                
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
                        draw.triangle(x - shift, y, th, tb, pi, sclr, 'drawSupports');
                        hold on;
                    else
                        draw.triangle(x-shift*cos(thetaZ),y-shift*sin(thetaZ),th,tb,pi+thetaZ,sclr,'drawSupports');
                        hold on;
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
                        draw.springX_2D(x-shift,y,sh,sprclr,'drawSupports');
                        hold on;
                        text(x-shift-sh/2,y+sh/10,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',kx);
                    else
                        draw.springX_2D(x-shift*cos(thetaZ),y-shift*sin(thetaZ),sh,sprclr,'drawSupports',thetaZ);
                        hold on;
                        text(x-(shift+sh)*cos(thetaZ),y-(shift+sh)*sin(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',kx);
                    end
                end
                
                % Draw fixed support in Y direction
                if draw.mdl.nodes(n).ebc(2) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.triangle(x, y - shift, th, tb, 3*pi/2, sclr, 'drawSupports');
                        hold on;
                    else
                        draw.triangle(x+shift*sin(thetaZ),y-shift*cos(thetaZ),th,tb,3*pi/2+thetaZ,sclr,'drawSupports');
                        hold on;
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
                        draw.springY_2D(x,y-shift,sh,sprclr,'drawSupports');
                        hold on;
                        text(x+sh/10,y-shift-sh/1.5,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',ky);
                    else
                        draw.springY_2D(x+shift*sin(thetaZ),y-shift*cos(thetaZ),sh,sprclr,'drawSupports',thetaZ);
                        hold on;
                        text(x+(shift+sh)*sin(thetaZ),y-(shift+sh)*cos(thetaZ),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',ky);
                    end
                end
                
                % Draw fixed support in Z direction
                if draw.mdl.nodes(n).ebc(6) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.square(x,y,ss,sclr,'drawSupports');
                        hold on;
                    else
                        % Rotated square
                        if dir(2) >= 0
                            cz = cos(thetaZ + pi/4);
                            sz = sin(thetaZ + pi/4);
                        else
                            cz = cos(thetaZ - pi/4);
                            sz = sin(thetaZ - pi/4);
                        end
                        S = sqrt(2)*ss/2;
                        X = [x-S*cz,x+S*sz,x+S*cz,x-S*sz];
                        Y = [y-S*sz,y-S*cz,y+S*sz,y+S*cz];
                        fill(X,Y,sclr,'tag','drawSupports');
                        hold on;
                    end
                
                % Draw spring support in Z direction
                elseif draw.mdl.nodes(n).ebc(6) == SPRING_DOF
                    krz = draw.mdl.nodes(n).springStiff(6);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if krz >= 1000
                            value = sprintf('%.*e kNm/rad',dc,krz);
                        else
                            value = sprintf('%.*f kNm/rad',dc,krz);
                        end
                    else
                        if krz >= 1000
                            value = sprintf('%.*e',dc,krz);
                        else
                            value = sprintf('%.*f',dc,krz);
                        end
                    end
                    draw.rotSpringZ_2D(x,y,sh,sprclr,'drawSupports')
                    hold on;
                    text(x+0.45*sh,y,value,'HorizontalAlignment','left','VerticalAlignment','top','Color',sprclr,'Fontsize',8.5,'tag','textRotSprings','UserData',krz);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws elements with hinged or continuous ends.
        function draw = elements(draw)
            % Get flag for semi-rigid joint visualization option
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSrj = get(mdata.viewSemiRigidButton,'Checked');
            
            % Parameters
            r      = draw.size/125;            % hinge symbol radius
            sh     = draw.size/35;             % spring height
            clr    = [0,0,0];                  % element color
            nclr   = [0,0,0];                  % node and hinge color
            sprclr = [0.6,0,0.4];              % spring color
            dc     = getappdata(0,'decPrec');  % decimal precision
            
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with X and Y axes
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get element incidence information on nodes
                [toti,hei] = draw.mdl.nodes(n1).elemsIncidence(draw.mdl);
                [totf,hef] = draw.mdl.nodes(n2).elemsIncidence(draw.mdl);
                
                % Set element end coordinates
                xi = x1;
                yi = y1;
                xf = x2;
                yf = y2;
                
                % Set element initial coordinates by checking if there is a hinge on nodal point position or on element end
                if hei == toti && draw.mdl.nodes(n1).ebc(6) == 0 % Hinge on node
                    draw.circle(x1,y1,r,nclr,'drawElements');
                    hold on;
                    xi = x1 + r * cx;
                    yi = y1 + r * cy;
                elseif draw.mdl.elems(e).hingei == 0 % Hinge on element end
                    draw.circle(x1 + r * cx, y1 + r * cy, r, clr,'drawElements');
                    xi = x1 + 2 * r * cx;
                    yi = y1 + 2 * r * cy;
                elseif draw.mdl.elems(e).hingei == 2
                    xi = x1 + sh * cx * 0.45;
                    yi = y1 + sh * cy * 0.45;
                    if strcmp(drawSrj,'on')
                        draw.srjoint([xi,yi],sh,[0 0 1],sprclr);
                        hold on
                        krzi = draw.mdl.elems(e).kri(3);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            if krzi >= 1000
                                value = sprintf('%.*e kNm/rad',dc,krzi);
                            else
                                value = sprintf('%.*f kNm/rad',dc,krzi);
                            end
                        else
                            if krzi >= 1000
                                value = sprintf('%.*e',dc,krzi);
                            else
                                value = sprintf('%.*f',dc,krzi);
                            end
                        end
                        text(xi+0.35*sh,yi,value,'HorizontalAlignment','left','VerticalAlignment','top','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krzi);
                    else
                        % Connect element end coordinates
                        X_srj = [xi,x1];
                        Y_srj = [yi,y1];
                        line(X_srj,Y_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on
                    end
                end
                
                % Set element final coordinates by checking if there is a hinge on nodal point position or on element end
                if hef == totf && draw.mdl.nodes(n2).ebc(6) == 0 % Hinge on node
                    draw.circle(x2,y2,r,nclr,'drawElements');
                    hold on;
                    xf = x2 - r * cx;
                    yf = y2 - r * cy;
                elseif draw.mdl.elems(e).hingef == 0 % Hinge on element end
                    draw.circle(x2 - r * cx, y2 - r * cy, r, clr,'drawElements');
                    xf = x2 - 2 * r * cx;
                    yf = y2 - 2 * r * cy;
                elseif draw.mdl.elems(e).hingef == 2
                    xf = x2 - sh * cx * 0.45;
                    yf = y2 - sh * cy * 0.45;
                    if strcmp(drawSrj,'on')
                        draw.srjoint([xf,yf],sh,[0 0 1],sprclr);
                        hold on
                        krzf = draw.mdl.elems(e).krf(3);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            if krzf >= 1000
                                value = sprintf('%.*e kNm/rad',dc,krzf);
                            else
                                value = sprintf('%.*f kNm/rad',dc,krzf);
                            end
                        else
                            if krzf >= 1000
                                value = sprintf('%.*e',dc,krzf);
                            else
                                value = sprintf('%.*f',dc,krzf);
                            end
                        end
                        text(xf+0.35*sh,yf,value,'HorizontalAlignment','left','VerticalAlignment','top','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krzf);
                    else
                        % Connect element end coordinates
                        X_srj = [xf,x2];
                        Y_srj = [yf,y2];
                        line(X_srj,Y_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on
                    end
                end
                
                % Connect element end coordinates
                X = [xi,xf];
                Y = [yi,yf];
                line(X,Y,'Color',clr,'tag','drawElements');
                hold on
            end
        end
        
        %------------------------------------------------------------------
        % Plots a 2D semi-rigid joint using a spiral
        % This method is used to draw semi-rigid joints on 2D models.
        % Input arguments:
        %  coords: semi-rigid joint coordinates ([x, y, z])
        %  sz: size parameter (2D - Spring Height, 3D - Sphere Radius)
        %  dir: direction to wich semi-rigid joint is applied (local axis)
        %  clr: semi-rigid joint symbol color (RGB vector)
        function draw = srjoint(draw,coords,sz,~,clr)
            % Get coordinates
            x = coords(1);
            y = coords(2);
            % Spring center height
            sch = sz/2;
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
            plot(xs,ys,'color',clr,'Linewidth',1.1,'tag','drawSemiRigid')
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
            scl = getappdata(0,'load_sf');   % Load symbol scale
            rmappdata(0,'max_val')
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
                    Qyi = abs(scl * qyi);
                    Qyf = abs(scl * qyf);
                    
                    % Calculate parameters that depends on the existence of hinges
                    if ((draw.mdl.elems(e).hingei == CONTINUOUS_END) && (draw.mdl.elems(e).hingef == CONTINUOUS_END))
                        step = L / round(20 * L/draw.size);
                        p1 = 0;
                        p2 = L;
                        r1 = 0;
                        r2 = 0;
                    elseif((draw.mdl.elems(e).hingei == HINGED_END || draw.mdl.elems(e).hingei == SEMIRIGID_END) && (draw.mdl.elems(e).hingef == CONTINUOUS_END))
                        step = (L-r) / round(20 * (L-r)/draw.size);
                        p1 = r;
                        p2 = L;
                        r1 = r;
                        r2 = 0;
                    elseif((draw.mdl.elems(e).hingei == CONTINUOUS_END) && (draw.mdl.elems(e).hingef == HINGED_END || draw.mdl.elems(e).hingef == SEMIRIGID_END))
                        step = (L-r) / round(20 * (L-r)/draw.size);
                        p1 = 0;
                        p2 = L - r;
                        r1 = 0;
                        r2 = r;
                    elseif((draw.mdl.elems(e).hingei == HINGED_END || draw.mdl.elems(e).hingei == SEMIRIGID_END) &&...
                           (draw.mdl.elems(e).hingef == HINGED_END || draw.mdl.elems(e).hingef == SEMIRIGID_END))
                        step = (L-2*r) / round(20 * (L-2*r)/draw.size);
                        p1 = r;
                        p2 = L - r;
                        r1 = r;
                        r2 = r;
                    end
                    
                    % Draw load symbol on a number cross-sections along element local axis X
                    for x = p1:step:p2
                        % Calculate load values on current cross-section
                        qx = A * x + B;
                        qy = C * x + D;
                        Qy = abs(scl * qy);
                        
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
                            [xd,yd] = draw.coordTransf2D(L/2, (Qy+r), x1, y1, alpha);
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
                        
                        if qyi > 0
                            [xb,yb] = draw.coordTransf2D(r1, -Qyi - r, x1, y1, alpha);
                            [xh,yh] = draw.coordTransf2D(0, -(Qyi+r), x1, y1, alpha);
                        elseif qyi < 0
                            [xb,yb] = draw.coordTransf2D(r1, Qyi + r, x1, y1, alpha);
                            [xh,yh] = draw.coordTransf2D(0, Qyi+r, x1, y1, alpha);
                        elseif qyi == 0
                            if qyf > 0
                                [xb,yb] = draw.coordTransf2D(r1, -r, x1, y1, alpha);
                                [xh,yh] = draw.coordTransf2D(0, -(Qyi+r), x1, y1, alpha);
                            elseif qyf < 0
                                [xb,yb] = draw.coordTransf2D(r1, r, x1, y1, alpha);
                                [xh,yh] = draw.coordTransf2D(0, Qyi+r, x1, y1, alpha);
                            end
                        end
                        
                        if qyf > 0
                            [xc,yc] = draw.coordTransf2D(L-r2, -Qyf - r, x1, y1, alpha);
                            [xj,yj] = draw.coordTransf2D(L, -(Qyf+r), x1, y1, alpha);
                        elseif qyf < 0
                            [xc,yc] = draw.coordTransf2D(L-r2, Qyf + r, x1, y1, alpha);
                            [xj,yj] = draw.coordTransf2D(L, Qyf+r, x1, y1, alpha);
                        elseif qyf == 0
                            if qyi > 0
                                [xc,yc] = draw.coordTransf2D(L-r2, -r, x1, y1, alpha);
                                [xj,yj] = draw.coordTransf2D(L, -(Qyf+r), x1, y1, alpha);
                            elseif qyi < 0
                                [xc,yc] = draw.coordTransf2D(L-r2, r, x1, y1, alpha);
                                [xj,yj] = draw.coordTransf2D(L, Qyf+r, x1, y1, alpha);
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
                        if (x ~= 0) && (x ~= L) &&((qxi ~= 0) || (qxf ~= 0))
                            draw.triangle(xs, ys, 2*ah, 2*ab, dir_x, clr,'drawElemLoads');
                        end
                        
                        % Draw transversal load symbol
                        if Qy >= ah
                            draw.arrow2D(xa, ya, Qy, ah, ab, dir_y, clr,'drawElemLoads');
                        elseif (abs(x-p1) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [x_init,y_init] = draw.coordTransf2D(p1, 0, x1, y1, alpha);
                            X = [x_init, xb];
                            Y = [y_init, yb];
                            line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                        elseif (abs(x-p2) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [x_end,y_end] = draw.coordTransf2D(p2, 0, x1, y1, alpha);
                            X = [x_end, xc];
                            Y = [y_end, yc];
                            line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                        end
                    end
                    
                    % Connect diagram extremities
                    if  (qyi < 0) && (qyf > 0)
                        x0 = abs((Qyi*(L-r1-r2))/(Qyf+Qyi));
                        [xu,yu] = draw.coordTransf2D(x0, r, x1, y1, alpha);
                        [xl,yl] = draw.coordTransf2D(x0, -r, x1, y1, alpha);
                        
                        X = [xb, xu, xl, xc];
                        Y = [yb, yu, yl, yc];
                        line(X, Y, 'Color', clr, 'tag','drawElemLoads');
                    elseif (qyi > 0) && (qyf < 0)
                        x0 = abs((Qyi*(L-r1-r2))/(Qyf+Qyi));
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
        % Draws applied nodal loads and moments.
        function draw = nodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            lshift = draw.size/100;  % distance between load symbol and nodal point
            mshift = draw.size/75;   % distance between applied moment symbol and nodal point
            al = draw.size/12;       % load symbol size (arrow length)
            ah = draw.size/60;       % load symbol size (arrowhead height)
            ab = draw.size/60;       % load symbol size (arrowhead base)
            mr = draw.size/30;       % moment load symbol radius
            clr = [1,0,0];           % load color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                % Check if current node has a nodal load
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.static(1);
                    fy = draw.mdl.nodes(n).load.static(2);
                    mz = draw.mdl.nodes(n).load.static(6);
                    
                    % Draw horizontal load component
                    if fx > 0
                        draw.arrow2D(x - lshift, y, al, ah, ab, pi, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-lshift-(ah+al)/3,y,value,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow2D(x + lshift, y, al, ah, ab, 0, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+lshift+(ah+al)/3,y,value,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw vertical load component
                    if fy > 0
                        draw.arrow2D(x, y - lshift, al, ah, ab, 3 * pi/2, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y-lshift-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow2D(x, y + lshift, al, ah, ab, pi/2, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y+lshift+(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                    
                    % Draw applied moment component
                    if mz > 0
                        draw.moment2D(draw, x + mshift, y, mr, 'z+', clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x+2.8*mshift,y+mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mz));
                        
                    elseif mz < 0
                        draw.moment2D(draw, x + mshift, y, mr, 'z-', clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x+2.8*mshift,y-mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mz));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws applied dynamic nodal loads and moments.
        function draw = dynamicNodalLoads(draw)
            % Get handle to GUI_Main
            mdata = guidata(findobj('Tag','GUI_Main'));
            
            % Check if nodal loads visualization is on and analysis is not modal
            include_constants;
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off') ||...
               draw.mdl.drv.whichResponse == MODAL_ANALYSIS
                return
            end
            
            % Parameters
            lshift = draw.size/100;  % distance between load symbol and nodal point
            mshift = draw.size/75;   % distance between applied moment symbol and nodal point
            m = draw.size/60;        % concentrated mass symbol radius
            al = draw.size/12;       % load symbol size (arrow length)
            ah = draw.size/60;       % load symbol size (arrowhead height)
            ab = draw.size/60;       % load symbol size (arrowhead base)
            mr = draw.size/30;       % moment load symbol radius
            clr = [0,0.7,0];         % load color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                % Check if current node has a nodal load
                if isempty(draw.mdl.nodes(n).load.dynamic) == 0
                    if isempty(draw.mdl.nodes(n).displMass) == 0 && draw.mdl.nodes(n).displMass > 0
                        lshift = m;
                    end
                    
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.dynamic(1);
                    fy = draw.mdl.nodes(n).load.dynamic(2);
                    mz = draw.mdl.nodes(n).load.dynamic(6);
                    
                    % Draw horizontal load component
                    if fx > 0
                        draw.arrow2D(x - lshift, y, al, ah, ab, pi, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-lshift-(ah+al)/3,y,value,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow2D(x + lshift, y, al, ah, ab, 0, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+lshift+(ah+al)/3,y,value,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw vertical load component
                    if fy > 0
                        draw.arrow2D(x, y - lshift, al, ah, ab, 3 * pi/2, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y-lshift-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow2D(x, y + lshift, al, ah, ab, pi/2, clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x+al/25,y+lshift+(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                    
                    % Draw applied moment component
                    if mz > 0
                        draw.moment2D(draw, x + mshift, y, mr, 'z+', clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x+2.8*mshift,y+mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mz));
                        
                    elseif mz < 0
                        draw.moment2D(draw, x + mshift, y, mr, 'z-', clr, 'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x+2.8*mshift,y-mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mz));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws concentrated nodal mass.
        function draw = nodalMass(draw)
            % Check if nodal mass visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalMassButton,'Checked'),'off')
                return
            end
            
            % Parameters
            r   = draw.size/60;             % concentrated mass symbol radius
            clr = [0,0,1];                  % mass color
            dc  = getappdata(0,'decPrec');  % decimal precision
            
             for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).displMass)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    mass = draw.mdl.nodes(n).displMass;
                    
                    % Draw mass
                    if mass > 0
                        s = draw.sphere(x,y,0,r,'drawNodalMass');
                        set(s,'Edgecolor',clr,'FaceColor',clr);
                        
                        % Write mass value
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
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
            shift  = 0;                        % distance between presc. displ. symbol and support
            mshift = draw.size/75;             % distance between moment symbol and nodal point
            al     = draw.size/12;             % presc. displ. symbol size (arrow length)
            ah     = draw.size/60;             % presc. displ. symbol size (arrowhead height)
            ab     = draw.size/60;             % presc. displ. symbol size (arrowhead base)
            pr     = draw.size/30;             % presc. rotation symbol size (radius)
            clr    = [1,0,1];                  % presc. displ. symbol color
            dc     = getappdata(0,'decPrec');  % decimal precision
            
            % Translation constraint symbol (triangle height)
            if strcmp(drawSupports,'on')
                th = draw.size/35;
            else
                th = draw.size/100;
            end
            
            for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).prescDispl)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get prescribed displacement component values and convert it to mm and rad
                    dx = 1000 * draw.mdl.nodes(n).prescDispl(1);
                    dy = 1000 * draw.mdl.nodes(n).prescDispl(2);
                    rz = draw.mdl.nodes(n).prescDispl(6);
                    
                    % Check if rotation is really fixed and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(6) == 1
                        ss = draw.size/70; % rotation constraint square side
                        if rz ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f rad',dc,abs(rz));
                            else
                                value = sprintf('%.*f',dc,abs(rz));
                            end
                            if rz > 0
                                draw.moment2D(draw,x+mshift,y,pr,'z+',clr,'drawPrescDispl');
                                text(x+2.8*mshift,y+pr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textPrescRot','UserData',abs(rz));
                            else
                                draw.moment2D(draw,x+mshift,y,pr,'z-',clr,'drawPrescDispl');                                
                                text(x+2.8*mshift,y-pr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textPrescRot','UserData',abs(rz));
                            end
                        end
                    else
                        [tot,he] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                        if he == tot && tot > 0
                            ss = draw.size/125; % hinge symbol radius
                        else
                            ss = 0;
                        end
                    end
                    
                    % Get rotation
                    if draw.mdl.nodes(n).isInclinedSupp
                        s_rot = ss + th + shift;
                        dir = draw.mdl.nodes(n).inclSuppDir(1:2);
                        if dir(2) >= 0
                            thetaZ = acos(dir(1));
                        else
                            thetaZ = -acos(dir(1));
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
                                draw.arrow2D(x-(ss+th+shift),y,al,ah,ab,pi,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x-(ss+th+shift+al),y,al,ah,ab,0,clr,'drawPrescDispl');
                            end
                            text(x-(ss+th+shift)-(al+ah)/2,y+al/20,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                        else
                            if dx > 0
                                draw.arrow2D(x-s_rot*cos(thetaZ),y-s_rot*sin(thetaZ),al,ah,ab,pi+thetaZ,clr,'drawPrescDispl');
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
                                draw.arrow2D(x,y-(ss+th+shift),al,ah,ab,-pi/2,clr,'drawPrescDispl');
                            else
                                draw.arrow2D(x,y-(ss+th+shift+al),al,ah,ab,pi/2,clr,'drawPrescDispl');
                            end
                            text(x+al/25,y-(ss+th+shift)-(ah+al)/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color', clr,'tag','textPrescDispl','UserData',abs(dy));
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
                return;
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
                if draw.mdl.nodes(n).initCond(1,1) ~= 0 || draw.mdl.nodes(n).initCond(2,1) ~= 0 || draw.mdl.nodes(n).initCond(6,1) ~= 0
                    flag = flag + 1;
                end
                if draw.mdl.nodes(n).initCond(1,2) ~= 0 || draw.mdl.nodes(n).initCond(2,2) ~= 0 || draw.mdl.nodes(n).initCond(6,2) ~= 0
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
                    text(x+r,y-0.3*r,str,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'tag','textInitialConditions');
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
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            anl   = get(mdata.popupmenu_AnalysisType,'Value');
            nm    = draw.size/200;   % node mark symbol (square side)
            ss    = draw.size/35;    % rotation constraint symbol (square side)
            r     = draw.size/125;   % hinge symbol radius
            m     = draw.size/60;    % concentrated mass symbol radius
            
            for n = 1:draw.mdl.nnp
                [tot,he] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                id = sprintf('%d',n);
                
                if draw.mdl.nodes(n).ebc(6) == FIXED_DOF && strcmp(get(mdata.viewSupportsButton,'Checked'),'on')
                    text(x+0.55*ss,y+0.55*ss,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                elseif anl == 2                             &&...
                      ~isempty(draw.mdl.nodes(n).displMass) &&...
                       draw.mdl.nodes(n).displMass > 0      &&...
                       strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x+0.75*m,y+0.75*m,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                elseif he == tot && tot > 0
                    text(x+r,y+r,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                else
                    text(x+0.8*nm,y+0.8*nm,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
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
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Estimate maximum displacement of each element (nodal and internal)
            m = zeros(1,draw.mdl.nel);
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                
                % Get maximum nodal displacements
                dx1 = draw.mdl.D(draw.mdl.ID(1,n1));
                dy1 = draw.mdl.D(draw.mdl.ID(2,n1));
                dx2 = draw.mdl.D(draw.mdl.ID(1,n2));
                dy2 = draw.mdl.D(draw.mdl.ID(2,n2));
                nodeDispl = [dx1,dy1,dx2,dy2];
                maxNode = max(abs(nodeDispl));
                
                % Get maximum estimated internal displacement
                maxInt = max(max(abs(draw.mdl.elems(e).intDispl)));
                
                % Get maximum element displacement
                m(e) = max(maxNode,maxInt);
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
        % Computes dynamic deformed configuration scale factor.
        function draw = dynamicDeformScaleFactor(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Estimate maximum displacement of each element (nodal and internal)
            m = zeros(1,draw.mdl.nel);
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                
                % Get maximum nodal displacements
                ids = [ draw.mdl.ID(1,n1) ;
                        draw.mdl.ID(2,n1) ;
                        draw.mdl.ID(1,n2) ;
                        draw.mdl.ID(2,n2) ];
                
                nodeDispl = sum(draw.mdl.results.dynamicDispl(ids,:,:),3);
                maxNode   = max(max(abs(nodeDispl)));
                
                % Get maximum estimated internal displacement
                maxInt = max(max(max(abs(draw.mdl.elems(e).dynamicIntDispl))));
                
                % Get maximum element displacement
                m(e) = max(maxNode,maxInt);
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
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get element length and orientation cosines
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get 50 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords;
                
                % Get element axial and transversal internal displacements in local system
                dl = draw.mdl.elems(e).intDispl;
                
                % Assemble rotation transformation matrix
                rot = [ cx  cy;
                       -cy  cx ];
                    
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords(1:2,:) + scale * dg;
                
                % Plot deformed configuration
                line(dfg(1,:),dfg(2,:),'Color',clr,'tag','drawDeformConfig');
            end
        end
        
        %------------------------------------------------------------------
        % Computes axial force diagram scale factor value.
        function draw = axialScaleFactor(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Get maximum axial internal force value of each element
            max_elem = zeros(1,draw.mdl.nel);
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get maximum value at element ends
                N1 = draw.mdl.elems(e).axial_force(:,1);
                N2 = draw.mdl.elems(e).axial_force(:,2);
                max_end = max(vertcat(abs(N1),abs(N2)));
                
                % Get maximum internal value
                if ~isempty(draw.mdl.elems(e).maxAxialForce)
                    max_int = abs(draw.mdl.elems(e).maxAxialForce(1));
                else
                    max_int = 0;
                end
                
                max_elem(e) = max(max_end,max_int);
            end
            
            % Get maximum axial internal force value of model
            max_val = max(max_elem);
            
            % Set adapted scale value
            if isempty(max_val) || max_val == 0
                asf = 0;
            else
                asf = draw.size/(2.5*sliderm*max_val);
            end
            setappdata(0,'axial_sf',asf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function draw = axialForce(draw,scale)
            % Parameters
            include_constants;
            clr = [1,0,0];     % diagram line color
            mdata = guidata(findobj('Tag','GUI_Main'));    % handle to main GUI
            dc = getappdata(0,'decPrec'); % decimal precision
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get element internal axial force value at both ends
                N1 = -draw.mdl.elems(e).axial_force(1);
                N2 = draw.mdl.elems(e).axial_force(2);
                
                % Avoid plotting garbage
                if abs(N1) < 10e-10
                    N1 = 0;
                end
                if abs(N2) < 10e-10
                    N2 = 0;
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Calculate diagram coordinates according to element orientation
                if (x1 <= x2) && (y1 <= y2)
                    xd1 = x1 - scale * N1 * abs(cy);
                    yd1 = y1 + scale * N1 * abs(cx);
                    xd2 = x2 - scale * N2 * abs(cy);
                    yd2 = y2 + scale * N2 * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    xd1 = x1 + scale * N1 * abs(cy);
                    yd1 = y1 - scale * N1 * abs(cx);
                    xd2 = x2 + scale * N2 * abs(cy);
                    yd2 = y2 - scale * N2 * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    xd1 = x1 + scale * N1 * abs(cy);
                    yd1 = y1 + scale * N1 * abs(cx);
                    xd2 = x2 + scale * N2 * abs(cy);
                    yd2 = y2 + scale * N2 * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1 = x1 - scale * N1 * abs(cy);
                    yd1 = y1 - scale * N1 * abs(cx);
                    xd2 = x2 - scale * N2 * abs(cy);
                    yd2 = y2 - scale * N2 * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Draw diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                line(Xi, Yi, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                line(Xf, Yf, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                
                % Text position
                if (xd1 == x1 && yd1 == y1) || ((xd1 >= x1 || yd1 >= y1) && (xd1 < x1 || yd1 > y1))
                    txt_pos1 = 'bottom';
                else
                    txt_pos1 = 'top';
                end
                if (xd2 == x2 && yd2 == y2) || ((xd2 >= x2 || yd2 >= y2) && (xd2 < x2 || yd2 > y2))
                    txt_pos2 = 'bottom';
                else
                    txt_pos2 = 'top';
                end
                
                % Write force values
                if abs(N1-N2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,N1);
                    else
                        value = sprintf('%+.*f',dc,N1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,value,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',N1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kN',dc,N1);
                        value2 = sprintf('%+.*f kN',dc,N2);
                    else
                        value1 = sprintf('%+.*f',dc,N1);
                        value2 = sprintf('%+.*f',dc,N2);
                    end
                    text(xd1,yd1,value1,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Rotation',ang,'tag','textAxialForceDiagram','UserData',N1);
                    text(xd2,yd2,value2,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos2,'Rotation',ang,'tag','textAxialForceDiagram','UserData',N2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate axial force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    N = draw.mdl.elems(e).intStresses(1,:);
                    
                    % Avoid plotting numeric garbage
                    if ~all(abs(N) < 10e-10)
                        % Get element 50 point division coords and internal stress values
                        coords = draw.mdl.elems(e).intCoords;
                        Nmax = draw.mdl.elems(e).maxAxialForce;
                        
                        % Get element basis transformation matrix
                        rot = draw.mdl.elems(e).T;

                        % Compute diagram coordinates
                        diagramCoords = rot' * [zeros(1,size(N,2)); scale*N; zeros(1,size(N,2))] + coords;

                        % Plot diagram
                        X = diagramCoords(1,:);
                        Y = diagramCoords(2,:);
                        line(X, Y, 'Color', clr, 'tag', 'drawAxialForceDiagram');

                        % Check if there is a maximum value within the diagram
                        if ~isempty(Nmax)
                            % Compute maximum stress value position, in global coordinates
                            NmaxGblCoords = [x1 + Nmax(2) * cx;
                                             y1 + Nmax(2) * cy;
                                             0];

                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0; scale*Nmax(1); 0] + NmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            scatter(xp, yp, 50, clr, '.', 'tag', 'drawAxialForceDiagram')
                            
                            % Avoid plotting text too close to element end
                            if abs(Nmax(2) - L) >= L/25 && Nmax(2) >= L/25
                                if (xp < (x1+x2)/2 && yp < (y1+y2)/2) || (xp >= (x1+x2)/2 && yp <= (y1+y2)/2)
                                    txt_pos = 'top';
                                else
                                    txt_pos = 'bottom';
                                end
                                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                    value = sprintf('%+.*f kN',dc,Nmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Nmax(1));
                                end
                                text(xp,yp,value,'HorizontalAlignment','center','VerticalAlignment',txt_pos,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmax(1));
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        line(X, Y, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    line(X, Y, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function axialForceEnvelop(draw,scale)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Nmax = draw.mdl.elems(e).intForcesEnvelop(1,:,1);
                Nmin = draw.mdl.elems(e).intForcesEnvelop(2,:,1);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Calculate text coordinates according to elem orientation
                if (x1 <= x2) && (y1 <= y2)
                    xd1_max = x1 - scale * Nmax(1) * abs(cy);
                    yd1_max = y1 + scale * Nmax(1) * abs(cx);
                    xd1_min = x1 - scale * Nmin(1) * abs(cy);
                    yd1_min = y1 + scale * Nmin(1) * abs(cx);
                    xd2_max = x2 - scale * Nmax(end) * abs(cy);
                    yd2_max = y2 + scale * Nmax(end) * abs(cx);
                    xd2_min = x2 - scale * Nmin(end) * abs(cy);
                    yd2_min = y2 + scale * Nmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    xd1_max = x1 + scale * Nmax(1) * abs(cy);
                    yd1_max = y1 - scale * Nmax(1) * abs(cx);
                    xd1_min = x1 + scale * Nmin(1) * abs(cy);
                    yd1_min = y1 - scale * Nmin(1) * abs(cx);
                    xd2_max = x2 + scale * Nmax(end) * abs(cy);
                    yd2_max = y2 - scale * Nmax(end) * abs(cx);
                    xd2_min = x2 + scale * Nmin(end) * abs(cy);
                    yd2_min = y2 - scale * Nmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    xd1_max = x1 + scale * Nmax(1) * abs(cy);
                    yd1_max = y1 + scale * Nmax(1) * abs(cx);
                    xd1_min = x1 + scale * Nmin(1) * abs(cy);
                    yd1_min = y1 + scale * Nmin(1) * abs(cx);
                    xd2_max = x2 + scale * Nmax(end) * abs(cy);
                    yd2_max = y2 + scale * Nmax(end) * abs(cx);
                    xd2_min = x2 + scale * Nmin(end) * abs(cy);
                    yd2_min = y2 + scale * Nmin(end) * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1_max = x1 - scale * Nmax(1) * abs(cy);
                    yd1_max = y1 - scale * Nmax(1) * abs(cx);
                    xd1_min = x1 - scale * Nmin(1) * abs(cy);
                    yd1_min = y1 - scale * Nmin(1) * abs(cx);
                    xd2_max = x2 - scale * Nmax(end) * abs(cy);
                    yd2_max = y2 - scale * Nmax(end) * abs(cx);
                    xd2_min = x2 - scale * Nmin(end) * abs(cy);
                    yd2_min = y2 - scale * Nmin(end) * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Null value
                if all(abs(Nmax) < 10e-10) && all(abs(Nmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],'Color',clr,'tag','drawAxialForceDiagram');
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',0);
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',0);
                    continue;
                end
                
                % Get coords of 27 points along element
                coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                
                % Get element basis transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Compute diagram coordinates
                maxCoords = rot' * [zeros(1,size(Nmax,2)); scale*Nmax; zeros(1,size(Nmax,2))] + coords;
                minCoords = rot' * [zeros(1,size(Nmin,2)); scale*Nmin; zeros(1,size(Nmin,2))] + coords;
                
                % Plot envelop
                Xmax = maxCoords(1,:);
                Ymax = maxCoords(2,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                line(Xmax, Ymax, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                line(Xmin, Ymin, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    line(xx, yy, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                end
                
                % Text position
                if (xd1_max == x1 && yd1_max == y1) || ((xd1_max >= x1 || yd1_max >= y1) && (xd1_max < x1 || yd1_max > y1))
                    txt_pos_max1 = 'bottom';
                else
                    txt_pos_max1 = 'top';
                end
                if (xd1_min == x1 && yd1_min == y1) || ((xd1_min >= x1 || yd1_min >= y1) && (xd1_min < x1 || yd1_min > y1))
                    txt_pos_min1 = 'bottom';
                else
                    txt_pos_min1 = 'top';
                end
                if (xd2_max == x2 && yd2_max == y2) || ((xd2_max >= x2 || yd2_max >= y2) && (xd2_max < x2 || yd2_max > y2))
                    txt_pos_max2 = 'bottom';
                else
                    txt_pos_max2 = 'top';
                end
                if (xd2_min == x2 && yd2_min == y2) || ((xd2_min >= x2 || yd2_min >= y2) && (xd2_min < x2 || yd2_min > y2))
                    txt_pos_min2 = 'bottom';
                else
                    txt_pos_min2 = 'top';
                end
                
                % Plot text containing values
                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                    value_Nmax1 = sprintf('%+.*f kN',dc,Nmax(1));
                    value_Nmin1 = sprintf('%+.*f kN',dc,Nmin(1));
                    value_Nmax2 = sprintf('%+.*f kN',dc,Nmax(end));
                    value_Nmin2 = sprintf('%+.*f kN',dc,Nmin(end));
                else
                    value_Nmax1 = sprintf('%+.*f',dc,Nmax(1));
                    value_Nmin1 = sprintf('%+.*f',dc,Nmin(1));
                    value_Nmax2 = sprintf('%+.*f',dc,Nmax(end));
                    value_Nmin2 = sprintf('%+.*f',dc,Nmin(end));
                end
                text(xd1_max,yd1_max,value_Nmax1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max1,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmax(1));
                text(xd1_min,yd1_min,value_Nmin1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min1,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmin(1));
                text(xd2_max,yd2_max,value_Nmax2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max2,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmax(end));
                text(xd2_min,yd2_min,value_Nmin2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min2,'Color',clr,'Rotation',ang,'tag','textAxialForceDiagram','UserData',Nmin(end));
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
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Get maximum internal shear force value of each element
            max_elem = zeros(1,draw.mdl.nel);
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get maximum value at element ends
                Q1_Y = draw.mdl.elems(e).shear_force_Y(:,1);
                Q2_Y = draw.mdl.elems(e).shear_force_Y(:,2);
                max_end = max(vertcat(abs(Q1_Y),abs(Q2_Y)));
                
                % Get maximum internal value
                if ~isempty(draw.mdl.elems(e).maxShearForce_XY)
                    max_int = abs(draw.mdl.elems(e).maxShearForce_XY(1));
                else
                    max_int = 0;
                end
                
                max_elem(e) = max(max_end,max_int);
            end
            
            % Get maximum shear internal force value of model
            max_val = max(max_elem);
            
            % Set adapted scale value
            if isempty(max_val) || max_val == 0
                ssf = 0;
            else
                ssf = draw.size/(2.5*sliderm*max_val);
            end
            setappdata(0,'shearXY_sf',ssf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function draw = shearForce_XY(draw,scale)
            % Parameters
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get element internal shear force value at both ends
                Q1 = draw.mdl.elems(e).shear_force_Y(1);
                Q2 = -draw.mdl.elems(e).shear_force_Y(2);
                
                % Avoid plotting garbage
                if abs(Q1) < 10e-10
                    Q1 = 0;
                end
                if abs(Q2) < 10e-10
                    Q2 = 0;
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Calculate diagram coordinates according to element orientation
                if (x1 <= x2) && (y1 <= y2)
                    xd1 = x1 - scale * Q1 * abs(cy);
                    yd1 = y1 + scale * Q1 * abs(cx);
                    xd2 = x2 - scale * Q2 * abs(cy);
                    yd2 = y2 + scale * Q2 * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    xd1 = x1 + scale * Q1 * abs(cy);
                    yd1 = y1 - scale * Q1 * abs(cx);
                    xd2 = x2 + scale * Q2 * abs(cy);
                    yd2 = y2 - scale * Q2 * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    xd1 = x1 + scale * Q1 * abs(cy);
                    yd1 = y1 + scale * Q1 * abs(cx);
                    xd2 = x2 + scale * Q2 * abs(cy);
                    yd2 = y2 + scale * Q2 * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1 = x1 - scale * Q1 * abs(cy);
                    yd1 = y1 - scale * Q1 * abs(cx);
                    xd2 = x2 - scale * Q2 * abs(cy);
                    yd2 = y2 - scale * Q2 * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Draw diagram extremities
                Xi = [x1,xd1];
                Yi = [y1,yd1];
                Xf = [x2,xd2];
                Yf = [y2,yd2];
                line(Xi,Yi,'Color',clr,'tag','drawShearForceXYDiagram');
                line(Xf,Yf,'Color',clr,'tag','drawShearForceXYDiagram');
                
                % Text position
                if (xd1 == x1 && yd1 == y1) || ((xd1 >= x1 || yd1 >= y1) && (xd1 < x1 || yd1 > y1))
                    txt_pos1 = 'bottom';
                else
                    txt_pos1 = 'top';
                end
                if (xd2 == x2 && yd2 == y2) || ((xd2 >= x2 || yd2 >= y2) && (xd2 < x2 || yd2 > y2))
                    txt_pos2 = 'bottom';
                else
                    txt_pos2 = 'top';
                end
                
                % Write force values
                if abs(Q1-Q2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,Q1);
                    else
                        value = sprintf('%+.*f',dc,Q1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,value,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Q1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kN',dc,Q1);
                        value2 = sprintf('%+.*f kN',dc,Q2);
                    else
                        value1 = sprintf('%+.*f',dc,Q1);
                        value2 = sprintf('%+.*f',dc,Q2);
                    end
                    text(xd1,yd1,value1,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Q1);
                    text(xd2,yd2,value2,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos2,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Q2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    Q = draw.mdl.elems(e).intStresses(2,:);
                    
                    % Avoid plotting numeric garbage
                    if ~all(abs(Q) < 10e-10)
                        % Get element 50 point division coords and internal stress values
                        coords = draw.mdl.elems(e).intCoords;
                        Qmax = draw.mdl.elems(e).maxShearForce_XY;
                        
                        % Get element basis transformation matrix
                        rot = draw.mdl.elems(e).T;
                        
                        % Compute diagram coordinates
                        diagramCoords = rot' * [zeros(1,size(Q,2)); scale*Q; zeros(1,size(Q,2))] + coords;
                        
                        % Plot diagram
                        X = diagramCoords(1,:);
                        Y = diagramCoords(2,:);
                        line(X,Y,'Color',clr,'tag','drawShearForceXYDiagram');
                        
                        % Check if there is a maximum value within the diagram
                        if ~isempty(Qmax)
                            % Compute maximum stress value position, in global coordinates
                            QmaxGblCoords = [x1 + Qmax(2) * cx;
                                             y1 + Qmax(2) * cy;
                                             0];

                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0; scale*Qmax(1); 0] + QmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            scatter(xp,yp,50,clr,'.','tag','drawShearForceXYDiagram')
                            
                            % Plot maximum value (avoid plotting text too close to element end)
                            if abs(Qmax(2) - L) >= L/25 && Qmax(2) >= L/25
                                if (xp < (x1+x2)/2 && yp < (y1+y2)/2) || (xp >= (x1+x2)/2 && yp <= (y1+y2)/2)
                                    txt_pos = 'top';
                                else
                                    txt_pos = 'bottom';
                                end
                                if strcmp(get(mdata.unitsButton,'Checked'),'on')
                                    value = sprintf('%+.*f kN',dc,Qmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Qmax(1));
                                end
                                text(xp,yp,value,'HorizontalAlignment','center','VerticalAlignment',txt_pos,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Qmax(1));
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1,xd2];
                        Y = [yd1,yd2];
                        line(X,Y,'Color',clr,'tag','drawShearForceXYDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1,xd2];
                    Y = [yd1,yd2];
                    line(X,Y,'Color',clr,'tag','drawShearForceXYDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XY plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function shearForceEnvelop_XY(draw,scale)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Qmax = draw.mdl.elems(e).intForcesEnvelop(1,:,2);
                Qmin = draw.mdl.elems(e).intForcesEnvelop(2,:,2);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Calculate text coordinates according to elem orientation
                if (x1 <= x2) && (y1 <= y2)
                    xd1_max = x1 - scale * Qmax(1) * abs(cy);
                    yd1_max = y1 + scale * Qmax(1) * abs(cx);
                    xd1_min = x1 - scale * Qmin(1) * abs(cy);
                    yd1_min = y1 + scale * Qmin(1) * abs(cx);
                    xd2_max = x2 - scale * Qmax(end) * abs(cy);
                    yd2_max = y2 + scale * Qmax(end) * abs(cx);
                    xd2_min = x2 - scale * Qmin(end) * abs(cy);
                    yd2_min = y2 + scale * Qmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    xd1_max = x1 + scale * Qmax(1) * abs(cy);
                    yd1_max = y1 - scale * Qmax(1) * abs(cx);
                    xd1_min = x1 + scale * Qmin(1) * abs(cy);
                    yd1_min = y1 - scale * Qmin(1) * abs(cx);
                    xd2_max = x2 + scale * Qmax(end) * abs(cy);
                    yd2_max = y2 - scale * Qmax(end) * abs(cx);
                    xd2_min = x2 + scale * Qmin(end) * abs(cy);
                    yd2_min = y2 - scale * Qmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    xd1_max = x1 + scale * Qmax(1) * abs(cy);
                    yd1_max = y1 + scale * Qmax(1) * abs(cx);
                    xd1_min = x1 + scale * Qmin(1) * abs(cy);
                    yd1_min = y1 + scale * Qmin(1) * abs(cx);
                    xd2_max = x2 + scale * Qmax(end) * abs(cy);
                    yd2_max = y2 + scale * Qmax(end) * abs(cx);
                    xd2_min = x2 + scale * Qmin(end) * abs(cy);
                    yd2_min = y2 + scale * Qmin(end) * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1_max = x1 - scale * Qmax(1) * abs(cy);
                    yd1_max = y1 - scale * Qmax(1) * abs(cx);
                    xd1_min = x1 - scale * Qmin(1) * abs(cy);
                    yd1_min = y1 - scale * Qmin(1) * abs(cx);
                    xd2_max = x2 - scale * Qmax(end) * abs(cy);
                    yd2_max = y2 - scale * Qmax(end) * abs(cx);
                    xd2_min = x2 - scale * Qmin(end) * abs(cy);
                    yd2_min = y2 - scale * Qmin(end) * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Null value
                if all(abs(Qmax) < 10e-10) && all(abs(Qmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],'Color',clr,'tag','drawShearForceXYDiagram');
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',0);
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',0);
                    continue;
                end
                
                 % Get coords of 27 points along element
                 coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                 
                 % Get element basis transformation matrix
                 rot = draw.mdl.elems(e).T;
                 
                 % Compute diagram coordinates
                 maxCoords = rot' * [zeros(1,size(Qmax,2)); scale*Qmax; zeros(1,size(Qmax,2))] + coords;
                 minCoords = rot' * [zeros(1,size(Qmin,2)); scale*Qmin; zeros(1,size(Qmin,2))] + coords;
                 
                 % Plot envelop
                 Xmax = maxCoords(1,:);
                 Ymax = maxCoords(2,:);
                 Xmin = minCoords(1,:);
                 Ymin = minCoords(2,:);
                 line(Xmax, Ymax, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                 line(Xmin, Ymin, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                 
                 % Plot lines to mark 6 internal points along envelop diagram
                 for i = plotValId
                     xx = [ minCoords(1,i) maxCoords(1,i) ];
                     yy = [ minCoords(2,i) maxCoords(2,i) ];
                     line(xx, yy, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                 end
                 
                 % Text position
                 if (xd1_max == x1 && yd1_max == y1) || ((xd1_max >= x1 || yd1_max >= y1) && (xd1_max < x1 || yd1_max > y1))
                    txt_pos_max1 = 'bottom';
                else
                    txt_pos_max1 = 'top';
                end
                if (xd1_min == x1 && yd1_min == y1) || ((xd1_min >= x1 || yd1_min >= y1) && (xd1_min < x1 || yd1_min > y1))
                    txt_pos_min1 = 'bottom';
                else
                    txt_pos_min1 = 'top';
                end
                if (xd2_max == x2 && yd2_max == y2) || ((xd2_max >= x2 || yd2_max >= y2) && (xd2_max < x2 || yd2_max > y2))
                    txt_pos_max2 = 'bottom';
                else
                    txt_pos_max2 = 'top';
                end
                if (xd2_min == x2 && yd2_min == y2) || ((xd2_min >= x2 || yd2_min >= y2) && (xd2_min < x2 || yd2_min > y2))
                    txt_pos_min2 = 'bottom';
                else
                    txt_pos_min2 = 'top';
                end
                 
                 % Plot text containing values
                 if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                    value_Qmax1 = sprintf('%+.*f kN',dc,Qmax(1));
                    value_Qmin1 = sprintf('%+.*f kN',dc,Qmin(1));
                    value_Qmax2 = sprintf('%+.*f kN',dc,Qmax(end));
                    value_Qmin2 = sprintf('%+.*f kN',dc,Qmin(end));
                else
                    value_Qmax1 = sprintf('%+.*f',dc,Qmax(1));
                    value_Qmin1 = sprintf('%+.*f',dc,Qmin(1));
                    value_Qmax2 = sprintf('%+.*f',dc,Qmax(end));
                    value_Qmin2 = sprintf('%+.*f',dc,Qmin(end));
                 end
                 text(xd1_max,yd1_max,value_Qmax1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max1,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Qmax(1));
                 text(xd1_min,yd1_min,value_Qmin1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min1,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Qmin(1));
                 text(xd2_max,yd2_max,value_Qmax2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max2,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Qmax(end));
                 text(xd2_min,yd2_min,value_Qmin2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min2,'Color',clr,'Rotation',ang,'tag','textShearForceXYDiagram','UserData',Qmin(end));
            end
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
        % Draws resulting shear force diagram in XZ plane on a given scale.
        function shearForceEnvelop_XZ(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XY plane.
        function draw = bendingMomentScaleFactor_XY(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Get maximum internal bending moment value of each element
            max_elem = zeros(1,draw.mdl.nel);
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get maximum value at element ends
                M1_Z = draw.mdl.elems(e).bending_moment_Z(:,1);
                M2_Z = draw.mdl.elems(e).bending_moment_Z(:,2);
                max_end = max(vertcat(abs(M1_Z),abs(M2_Z)));
                
                % Get maximum internal value
                if ~isempty(draw.mdl.elems(e).maxBendMoment_XY)
                    max_int = max(abs(draw.mdl.elems(e).maxBendMoment_XY(:,1)));
                else
                    max_int = 0;
                end
                
                max_elem(e) = max(max_end,max_int);
            end
            
            % Get maximum axial internal force value of model
            max_val = max(max_elem);
            
            % Set adapted scale value
            if isempty(max_val) || max_val == 0
                bsf = 0;
            else
                bsf = draw.size/(2.5*sliderm*max_val);
            end
            setappdata(0,'bendingXY_sf',bsf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function draw = bendingMoment_XY(draw,scale)
            % Parameters
            include_constants;
            clr = [1,0,0];     % diagram line color
            mdata = guidata(findobj('Tag','GUI_Main'));   % handle to main GUI
            dc = getappdata(0,'decPrec'); % decimal precision
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle with X axis
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get element internal force values ate both ends
                M1 = -draw.mdl.elems(e).bending_moment_Z(1);
                M2 = draw.mdl.elems(e).bending_moment_Z(2);
                
                % Avoid plotting garbage
                if abs(M1) < 10e-10
                    M1 = 0;
                end
                if abs(M2) < 10e-10
                    M2 = 0;
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Calculate diagram coordinates according to element orientation
                if ((x1 <= x2) && (y1 <= y2))
                    xd1 = x1 + scale * M1 * abs(cy);
                    yd1 = y1 - scale * M1 * abs(cx);
                    xd2 = x2 + scale * M2 * abs(cy);
                    yd2 = y2 - scale * M2 * abs(cx);
                    ang = 180*(acos(abs(cx)))/pi;
                elseif ((x1 >= x2) && (y1 >= y2))
                    xd1 = x1 - scale * M1 * abs(cy);
                    yd1 = y1 + scale * M1 * abs(cx);
                    xd2 = x2 - scale * M2 * abs(cy);
                    yd2 = y2 + scale * M2 * abs(cx);
                    ang = 180*(acos(abs(cx)))/pi;
                elseif ((x1 <= x2) && (y1 >= y2))
                    xd1 = x1 - scale * M1 * abs(cy);
                    yd1 = y1 - scale * M1 * abs(cx);
                    xd2 = x2 - scale * M2 * abs(cy);
                    yd2 = y2 - scale * M2 * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1 = x1 + scale * M1 * abs(cy);
                    yd1 = y1 + scale * M1 * abs(cx);
                    xd2 = x2 + scale * M2 * abs(cy);
                    yd2 = y2 + scale * M2 * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Draw diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                line(Xi, Yi, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                line(Xf, Yf, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                
                % Text position
                if (xd1 == x1 && yd1 == y1) || ((xd1 >= x1 || yd1 >= y1) && (xd1 < x1 || yd1 > y1))
                    txt_pos1 = 'bottom';
                else
                    txt_pos1 = 'top';
                end
                if (xd2 == x2 && yd2 == y2) || ((xd2 >= x2 || yd2 >= y2) && (xd2 < x2 || yd2 > y2))
                    txt_pos2 = 'bottom';
                else
                    txt_pos2 = 'top';
                end
                
                % Write bending moment values
                if abs(M1-M2) < 10e-10 && isempty(draw.mdl.elems(e).load.uniformGbl) && isempty(draw.mdl.elems(e).load.linearGbl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,M1);
                    else
                        value = sprintf('%+.*f',dc,M1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,value,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',M1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kNm',dc,M1);
                        value2 = sprintf('%+.*f kNm',dc,M2);
                    else
                        value1 = sprintf('%+.*f',dc,M1);
                        value2 = sprintf('%+.*f',dc,M2);
                    end
                    text(xd1,yd1,value1,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos1,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',M1);
                    text(xd2,yd2,value2,'Color',clr,'HorizontalAlignment','center','VerticalAlignment',txt_pos2,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',M2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if (isempty(draw.mdl.elems(e).load.uniformLcl) == 0) || (isempty(draw.mdl.elems(e).load.linearLcl) == 0)
                    M = draw.mdl.elems(e).intStresses(3,:);
                    
                    % Avoid plotting numeric garbage
                    if ~all(abs(M) < 10e-10)
                        % Get element 50 point division coords and internal stress values
                        coords = draw.mdl.elems(e).intCoords;
                        Mmax = draw.mdl.elems(e).maxBendMoment_XY;
                        
                        % Get element basis transformation matrix
                        rot = draw.mdl.elems(e).T;
                        
                        % Compute diagram coordinates
                        diagramCoords = rot' * [zeros(1,size(M,2)); -scale*M; zeros(1,size(M,2))] + coords;
                        
                        % Plot diagram
                        X = diagramCoords(1,:);
                        Y = diagramCoords(2,:);
                        line(X, Y, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                        
                        % Check if there is a maximum value within the diagram
                        if ~isempty(Mmax)
                            % Get maximum stress values
                            Mmax_val = (Mmax(:,1))';
                            
                            % Get number of maximum values within diagram
                            nm = size(Mmax_val,2);
                            
                            % Compute maximum stress value position, in global coordinates
                            MmaxGblCoords = [x1 + (Mmax(:,2))' * cx;
                                             y1 + (Mmax(:,2))' * cy;
                                             zeros(1,nm)];
                            
                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [zeros(1,nm); -scale*Mmax_val; zeros(1,nm)] + MmaxGblCoords;
                            xp = maxPointCoords(1,:);
                            yp = maxPointCoords(2,:);
                            scatter(xp, yp, 50, clr, '.', 'tag', 'drawBendMomentXYDiagram')
                            
                            % Plot maximum value in text (avoid plotting text too close to element end)
                            for np = 1:nm % 2 rounds max
                                if abs(Mmax(np,2) - L) >= L/25 && Mmax(np,2) >= L/25
                                    if (xp(np) < (x1+x2)/2 && yp(np) < (y1+y2)/2) || (xp(np) >= (x1+x2)/2 && yp(np) <= (y1+y2)/2)
                                        txt_pos = 'top';
                                    else
                                        txt_pos = 'bottom';
                                    end
                                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                        value = sprintf('%+.*f kNm',dc,Mmax_val(np));
                                    else
                                        value = sprintf('%+.*f',dc,Mmax_val(np));
                                    end
                                    text(xp(np),yp(np),value,'HorizontalAlignment','center','VerticalAlignment',txt_pos,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',Mmax_val(np));
                                end
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        line(X, Y, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    line(X, Y, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XY plane on a
        % given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function bendingMomentEnvelop_XY(draw,scale)
            % Parameters
            clr = [1,0,0];     % diagram line color
            mdata = guidata(findobj('Tag','GUI_Main')); % handle to main GUI
            dc = getappdata(0,'decPrec'); % decimal precision
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Mmax = draw.mdl.elems(e).intForcesEnvelop(1,:,3);
                Mmin = draw.mdl.elems(e).intForcesEnvelop(2,:,3);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Get element orientation angle cosine with axes X and Y
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Calculate text coordinates according to elem orientation
                if (x1 <= x2) && (y1 <= y2)
                    xd1_max = x1 + scale * Mmax(1) * abs(cy);
                    yd1_max = y1 - scale * Mmax(1) * abs(cx);
                    xd1_min = x1 + scale * Mmin(1) * abs(cy);
                    yd1_min = y1 - scale * Mmin(1) * abs(cx);
                    xd2_max = x2 + scale * Mmax(end) * abs(cy);
                    yd2_max = y2 - scale * Mmax(end) * abs(cx);
                    xd2_min = x2 + scale * Mmin(end) * abs(cy);
                    yd2_min = y2 - scale * Mmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 >= x2) && (y1 >= y2)
                    xd1_max = x1 - scale * Mmax(1) * abs(cy);
                    yd1_max = y1 + scale * Mmax(1) * abs(cx);
                    xd1_min = x1 - scale * Mmin(1) * abs(cy);
                    yd1_min = y1 + scale * Mmin(1) * abs(cx);
                    xd2_max = x2 - scale * Mmax(end) * abs(cy);
                    yd2_max = y2 + scale * Mmax(end) * abs(cx);
                    xd2_min = x2 - scale * Mmin(end) * abs(cy);
                    yd2_min = y2 + scale * Mmin(end) * abs(cx);
                    ang = 180*acos(abs(cx))/pi;
                elseif (x1 <= x2) && (y1 >= y2)
                    xd1_max = x1 - scale * Mmax(1) * abs(cy);
                    yd1_max = y1 - scale * Mmax(1) * abs(cx);
                    xd1_min = x1 - scale * Mmin(1) * abs(cy);
                    yd1_min = y1 - scale * Mmin(1) * abs(cx);
                    xd2_max = x2 - scale * Mmax(end) * abs(cy);
                    yd2_max = y2 - scale * Mmax(end) * abs(cx);
                    xd2_min = x2 - scale * Mmin(end) * abs(cy);
                    yd2_min = y2 - scale * Mmin(end) * abs(cx);
                    ang = -180*(acos(cx))/pi;
                else
                    xd1_max = x1 + scale * Mmax(1) * abs(cy);
                    yd1_max = y1 + scale * Mmax(1) * abs(cx);
                    xd1_min = x1 + scale * Mmin(1) * abs(cy);
                    yd1_min = y1 + scale * Mmin(1) * abs(cx);
                    xd2_max = x2 + scale * Mmax(end) * abs(cy);
                    yd2_max = y2 + scale * Mmax(end) * abs(cx);
                    xd2_min = x2 + scale * Mmin(end) * abs(cy);
                    yd2_min = y2 + scale * Mmin(end) * abs(cx);
                    ang = 180*(acos((cx)))/pi + 180;
                end
                
                % Null value
                if all(abs(Mmax) < 10e-10) && all(abs(Mmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],'Color',clr,'tag','drawBendMomentXYDiagram');
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',0);
                    text((x1+x2)/2,(y1+y2)/2,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',0);
                    continue;
                end
                
                % Get coords of 27 points along element
                coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                
                % Get element basis transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Compute diagram coordinates
                maxCoords = rot' * [zeros(1,size(Mmax,2)); -scale*Mmax; zeros(1,size(Mmax,2))] + coords;
                minCoords = rot' * [zeros(1,size(Mmin,2)); -scale*Mmin; zeros(1,size(Mmin,2))] + coords;
                
                % Plot envelop
                Xmax = maxCoords(1,:);
                Ymax = maxCoords(2,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                line(Xmax, Ymax, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                line(Xmin, Ymin, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    line(xx, yy, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                end
                
                % Text position
                if (xd1_max == x1 && yd1_max == y1) || ((xd1_max >= x1 || yd1_max >= y1) && (xd1_max < x1 || yd1_max > y1))
                    txt_pos_max1 = 'bottom';
                else
                    txt_pos_max1 = 'top';
                end
                if (xd1_min == x1 && yd1_min == y1) || ((xd1_min >= x1 || yd1_min >= y1) && (xd1_min < x1 || yd1_min > y1))
                    txt_pos_min1 = 'bottom';
                else
                    txt_pos_min1 = 'top';
                end
                if (xd2_max == x2 && yd2_max == y2) || ((xd2_max >= x2 || yd2_max >= y2) && (xd2_max < x2 || yd2_max > y2))
                    txt_pos_max2 = 'bottom';
                else
                    txt_pos_max2 = 'top';
                end
                if (xd2_min == x2 && yd2_min == y2) || ((xd2_min >= x2 || yd2_min >= y2) && (xd2_min < x2 || yd2_min > y2))
                    txt_pos_min2 = 'bottom';
                else
                    txt_pos_min2 = 'top';
                end
                
                % Plot text containing values
                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                    value_Mmax1 = sprintf('%+.*f kNm',dc,Mmax(1));
                    value_Mmin1 = sprintf('%+.*f kNm',dc,Mmin(1));
                    value_Mmax2 = sprintf('%+.*f kNm',dc,Mmax(end));
                    value_Mmin2 = sprintf('%+.*f kNm',dc,Mmin(end));
                else
                    value_Mmax1 = sprintf('%+.*f',dc,Mmax(1));
                    value_Mmin1 = sprintf('%+.*f',dc,Mmin(1));
                    value_Mmax2 = sprintf('%+.*f',dc,Mmax(end));
                    value_Mmin2 = sprintf('%+.*f',dc,Mmin(end));
                end
                text(xd1_max,yd1_max,value_Mmax1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max1,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',Mmax(1));
                text(xd1_min,yd1_min,value_Mmin1,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min1,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',Mmin(1));
                text(xd2_max,yd2_max,value_Mmax2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_max2,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',Mmax(end));
                text(xd2_min,yd2_min,value_Mmin2,'HorizontalAlignment','center','VerticalAlignment',txt_pos_min2,'Color',clr,'Rotation',ang,'tag','textBendMomentXYDiagram','UserData',Mmin(end));
            end
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
        % Draws resulting bending moment diagram in XZ plane on a given scale.
        function bendingMomentEnvelop_XZ(~,~)
        end
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        function draw = reactions(draw)
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            mshift = draw.size/75;  % distance between moment symbol and nodal point
            mr = draw.size/30;      % moment load symbol radius
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
                th = draw.size/100;
                sh = draw.size/100;
            end
            
            for n = 1:draw.mdl.nnp
                % Get nodal coordinates
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Get reactions values
                rx = draw.mdl.F(draw.mdl.ID(1,n));
                ry = draw.mdl.F(draw.mdl.ID(2,n));
                rz = draw.mdl.F(draw.mdl.ID(3,n));
                
                % Check if rotation is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(6)== FIXED_DOF) || (draw.mdl.nodes(n).ebc(6)== SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kNm',dc,abs(rz));
                    else
                        value = sprintf('%.*f',dc,abs(rz));
                    end
                    if rz >= 0
                        draw.moment2D(draw,x+mshift,y,mr,'z+',clr,'drawReactions');
                        text(x+2.8*mshift,y+mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textMomentReactions','UserData',abs(rz));
                    else
                        draw.moment2D(draw,x+mshift,y,mr,'z-',clr,'drawReactions');
                        text(x+2.8*mshift,y-mr,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textMomentReactions','UserData',abs(rz));
                    end
                end
                
                if draw.mdl.nodes(n).ebc(6)== FIXED_DOF
                    ss = draw.size/70; % Rotation constraint symbol (square side)
                else
                    [tot,he] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if he == tot && tot > 0
                        ss = draw.size/125; % hinge symbol radius
                    else
                        ss = 0;
                    end
                end
                
                % Support or spring height in the direction of axis X and Y
                if draw.mdl.nodes(n).ebc(1)== FIXED_DOF
                    hx = ss+th;
                elseif draw.mdl.nodes(n).ebc(1)== SPRING_DOF
                    hx = ss+sh;
                end
                if draw.mdl.nodes(n).ebc(2)== FIXED_DOF
                    hy = ss+th;
                elseif draw.mdl.nodes(n).ebc(2)== SPRING_DOF
                    hy = ss+sh;
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
            % Parameters
            clr = [1,0,0];  % deformed configuration line color
            aux_id = [1, 2:2:48, 49, 50]; % 27 auxiliary indexes for local coords
            
            % Get number of points per element
            npe = size(draw.mdl.elems(1).natVibration,2);
            
            % Initialize plotting matrix
            d = zeros(2,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 27 cross-sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get element length and orientation cosines
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get 28 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = draw.mdl.elems(e).natVibration(:,:,nMode);
                
                % Assemble rotation transformation matrix
                rot = [ cx  cy;
                       -cy  cx ];
                    
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords(1:2,:) + scale * dg;
                
                % Concatenate to displ mtx
                d(:,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan(2,1)];
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
            % Parameters
            clr = [1,0,0];  % deformed configuration line color
            aux_id = [1, 2:2:48, 49, 50]; % 27 auxiliary indexes for local coords
            
            % Check if step is not an integer
            dt = rem(step,1);
            if step >= draw.mdl.n_steps + 1
                step = draw.mdl.n_steps;
                dt = 1;
            end
            
            % Get number of points per element
            npe = size(draw.mdl.elems(1).dynamicIntDispl,2);
            
            % Initialize plotting matrix
            d = zeros(2,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 27 cross-sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get element length and orientation cosines
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                
                % Get 28 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = (1-dt) * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)) +...
                        dt  * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)+1);
                
                % Assemble rotation transformation matrix
                rot = [ cx  cy;
                       -cy  cx ];
                    
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords(1:2,:) + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan(2,1)];
            end
            % Plot deformed configuration
            plot(d(1,:), d(2,:), 'Color', clr, 'tag', 'drawDynamicDeform');
            drawnow
        end
    end
end