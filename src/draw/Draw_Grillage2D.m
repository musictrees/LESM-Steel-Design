%% Draw_Grillage2D class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *Grillage* draw object for 2D visualization.
%
classdef Draw_Grillage2D < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Grillage2D(mdl)
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
            nm     = draw.size/200;    % node mark symbol (square side)
            tb     = draw.size/50;     % translation constraint symbol (triangle base)
            ss     = draw.size/35;     % rotation constraint symbol (square side)
            sh     = draw.size/15;     % displacement spring symbol (spring height)
            sr     = draw.size/115;    % rotational spring symbol radius
            nclr   = [0,0,0];          % node and hinge color
            sclr   = [0.6,0.6,0.6];    % support color
            sprclr = [0.6,0,0.4];      % spring color
            dc     = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Draw nodal mark
                if draw.mdl.nodes(n).ebc(4) == FREE_DOF || draw.mdl.nodes(n).ebc(5) == FREE_DOF
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng < tot || tot == 0
                        draw.square(x,y,nm,nclr,'drawNodes');
                        hold on;
                    end
                end
                
                % Draw support conditions
                if ~strcmp(drawSupports,'on')
                    continue;
                end
                
                % Draw fixed support in Z direction
                if draw.mdl.nodes(n).ebc(3) == FIXED_DOF
                    draw.square(x,y,tb,sclr,'drawSupports');
                    hold on;
                    plot([x x-tb/2],[y y-tb/2],'color',[0 0 0],'tag','drawSupports');
                    plot([x x+tb/2],[y y-tb/2],'color',[0 0 0],'tag','drawSupports');
                    plot([x x-tb/2],[y y+tb/2],'color',[0 0 0],'tag','drawSupports');
                    plot([x x+tb/2],[y y+tb/2],'color',[0 0 0],'tag','drawSupports');
                    
                % Draw spring support in Z direction
                elseif draw.mdl.nodes(n).ebc(3) == SPRING_DOF
                    kz = draw.mdl.nodes(n).springStiff(3);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if kz >= 1000
                            value = sprintf('%.*e kN/m',dc,kz);
                        else
                            value = sprintf('%.*f kN/m',dc,kz);
                        end
                    else
                        if kz >= 1000
                            value = sprintf('%.*e',dc,kz);
                        else
                            value = sprintf('%.*f',dc,kz);
                        end
                    end
                    draw.square(x,y,tb,sclr,'drawSupports');
                    hold on;
                    draw.circle(x,y,sh/10,sprclr,'drawSupports');
                    hold on;
                    text(x-sh/5,y+sh/5,value,'HorizontalAlignment','right','VerticalAlignment','middle','Color',sprclr,'tag','textSprings','UserData',kz);
                end
                
                % Draw fixed rotation support
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF
                    draw.square(x, y, ss, sclr,'drawSupports');
                    hold on;
                
                % Draw spring rotation support
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF
                    krxy = draw.mdl.nodes(n).springStiff(4);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if krxy >= 1000
                            value = sprintf('%.*e kNm/rad',dc,krxy);
                        else
                            value = sprintf('%.*f kNm/rad',dc,krxy);
                        end
                    else
                        if krxy >= 1000
                            value = sprintf('%.*e',dc,krxy);
                        else
                            value = sprintf('%.*f',dc,krxy);
                        end
                    end
                    circ = 0:pi/50:2*pi;
                    xcirc = x + sr * cos(circ);
                    ycirc = y + sr * sin(circ);
                    fill(xcirc,ycirc,sprclr,'tag','drawSupports');
                    hold on;
                    text(x-sh/5,y-sh/5,value,'HorizontalAlignment','right','VerticalAlignment','middle','Color',sprclr,'tag','textRotSprings','UserData',krxy);
                end    
            end
        end
        
        %------------------------------------------------------------------
        % Draws elements with hinged or continuous ends.
        function draw = elements(draw)
            % Get flag for semi-rigid joint visualization option
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSrj = get(mdata.viewSemiRigidButton,'Checked');
            
            % Parameters
            r      = draw.size/125; % hinge symbol radius
            sr     = draw.size/35;  % semi-rigid joint symbol radius
            clr    = [0,0,0];       % element color
            nclr   = [0,0,0];       % node and hinge color
            sprclr = [0.6,0,0.4];   % spring color
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
                if hei == toti && (draw.mdl.nodes(n1).ebc(4) == 0 || draw.mdl.nodes(n1).ebc(5) == 0) % Hinge on node
                    draw.circle(x1,y1,r,nclr,'drawElements');
                    hold on;
                    xi = x1 + r * cx;
                    yi = y1 + r * cy;
                elseif draw.mdl.elems(e).hingei == 0 % Hinge on element end
                    draw.circle(x1+r*cx,y1+r*cy,r,clr,'drawElements');
                    xi = x1 + 2 * r * cx;
                    yi = y1 + 2 * r * cy;
                elseif draw.mdl.elems(e).hingei == 2
                    xi = x1 + sr * cx * 0.45;
                    yi = y1 + sr * cy * 0.45;
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),~] = draw.mdl.elems(e).locAxis;
                        draw.srjoint([xi,yi,0],sr,dir,sprclr);
                        hold on;
                        krxi = draw.mdl.elems(e).kri(1);
                        kryi = draw.mdl.elems(e).kri(2);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            if krxi >= 1000
                                value_x = sprintf('%.*e kNm/rad',dc,krxi);
                            else
                                value_x = sprintf('%.*f kNm/rad',dc,krxi);
                            end
                            if kryi >= 1000
                                value_y = sprintf('%.*e kNm/rad',dc,kryi);
                            else
                                value_y = sprintf('%.*f kNm/rad',dc,kryi);
                            end
                        else
                            if krxi >= 1000
                                value_x = sprintf('%.*e',dc,krxi);
                            else
                                value_x = sprintf('%.*f',dc,krxi);
                            end
                            if kryi >= 1000
                                value_y = sprintf('%.*e',dc,kryi);
                            else
                                value_y = sprintf('%.*f',dc,kryi);
                            end
                        end
                        text(xi+1.1*sr*cx,yi+1.1*sr*cy,value_x,'HorizontalAlignment','center','VerticalAlignment','top','Color',[0 0.65 0],'Fontsize',8.5,'tag','textSemiRigid','UserData',krxi);
                        text(xi-1.1*sr*cy,yi+1.1*sr*cx,value_y,'HorizontalAlignment','center','VerticalAlignment','top','Color',[0 0 1],'Fontsize',8.5,'tag','textSemiRigid','UserData',kryi);
                    else
                        % Connect element end coordinates
                        X_srj = [xi,x1];
                        Y_srj = [yi,y1];
                        line(X_srj, Y_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Set element final coordinates by checking if there is a hinge on nodal point position or on element end
                if hef == totf && (draw.mdl.nodes(n2).ebc(4) == 0 || draw.mdl.nodes(n2).ebc(5) == 0) % Hinge on node
                    draw.circle(x2,y2,r,nclr,'drawElements');
                    hold on;
                    xf = x2 - r * cx;
                    yf = y2 - r * cy;
                elseif draw.mdl.elems(e).hingef == 0 % Hinge on element end
                    draw.circle(x2-r*cx,y2-r*cy,r,clr,'drawElements');
                    xf = x2 - 2 * r * cx;
                    yf = y2 - 2 * r * cy;
                elseif draw.mdl.elems(e).hingef == 2
                    xf = x2 - sr * cx * 0.45;
                    yf = y2 - sr * cy * 0.45;
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),~] = draw.mdl.elems(e).locAxis;
                        draw.srjoint([xf,yf,0],sr,dir,sprclr);
                        hold on;
                        krxf = draw.mdl.elems(e).krf(1);
                        kryf = draw.mdl.elems(e).krf(2);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            if krxf >= 1000
                                value_x = sprintf('%.*e kNm/rad',dc,krxf);
                            else
                                value_x = sprintf('%.*f kNm/rad',dc,krxf);
                            end
                            if kryf >= 1000
                                value_y = sprintf('%.*e kNm/rad',dc,kryf);
                            else
                                value_y = sprintf('%.*f kNm/rad',dc,kryf);
                            end
                        else
                            if krxf >= 1000
                                value_x = sprintf('%.*e',dc,krxf);
                            else
                                value_x = sprintf('%.*f',dc,krxf);
                            end
                            if kryf >= 1000
                                value_y = sprintf('%.*e',dc,kryf);
                            else
                                value_y = sprintf('%.*f',dc,kryf);
                            end
                        end
                        text(xf+sr*cx,yf+sr*cy,value_x,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0.65 0],'Fontsize',8.5,'tag','textSemiRigid','UserData',krxf);
                        text(xf-sr*cy,yf+sr*cx,value_y,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 1],'Fontsize',8.5,'tag','textSemiRigid','UserData',kryf);
                    else
                        % Connect element end coordinates                        
                        X_srj = [xf,x2];
                        Y_srj = [yf,y2];
                        line(X_srj,Y_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Update element drawing end coordinates
                X = [xi,xf];
                Y = [yi,yf];
                line(X,Y,'Color',clr,'tag','drawElements');
                hold on
            end
        end
        
        %------------------------------------------------------------------
        % Plots sphere and moment 3D arrows to represent semi-rigid joint.
        % This method is used to draw semi-rigid joints on 3D models.
        % Input arguments:
        %  coords: semi-rigid joint coordinates ([x, y, z])
        %  sz: size parameter (2D - Spring Height, 3D - Sphere Radius)
        %  dir: direction to wich semi-rigid joint is applied (local axis)
        %  clr: semi-rigid joint symbol color (RGB vector)
        function draw = srjoint(draw,coords,sz,dir,clr)
            % Get coordinates
            x = coords(1);
            y = coords(2);
            
            % local axis
            dx = dir(1,:);
            dy = dir(2,:);

            % arrow properties
            l = 3 * sz;
            h = 0.8 * sz;
            B = 0.5 * sz;
            
            % Determine angle of arrow rotation
            ang = abs(acos(dx(1)));
            if dx(1,2) < 0
                ang = -ang;
            end
            
            % Draw double arrows to indicate the direction of local axis
            aux = [x,y,0] + l * dx;
            draw.arrow2D(aux(1),aux(2),l,h,B,pi+ang,[0 0.7 0],'drawSemiRigid');
            hold on;
            aux = [x,y,0] + (1-(h/l)) * l * dx;
            draw.arrow2D(aux(1),aux(2),(1-(h/l))*l,h,B,pi+ang,[0 0.7 0],'drawSemiRigid');
            hold on;
            
            aux = [x,y,0] + l * dy;
            draw.arrow2D(aux(1),aux(2),l,h,B,-pi/2+ang,[0 0 1],'drawSemiRigid');
            hold on;
            aux = [x,y,0] + (1-(h/l)) * l * dy;
            draw.arrow2D(aux(1),aux(2),(1-(h/l))*l,h,B,-pi/2+ang,[0 0 1],'drawSemiRigid');
            hold on;
            
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
        % Plots ID number of nodes.
        function draw = nodeID(draw)
            for n = 1:draw.mdl.nnp
                x  = draw.mdl.nodes(n).coord(1);
                y  = draw.mdl.nodes(n).coord(2);
                id = sprintf('%d',n);
                text(x+draw.size/150,y+draw.size/150,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
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
        % Computes element loads scale factor.
        function draw = elemLoadsScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws element distributed loads (uniform and linear).
        function draw = elemLoads(draw)
            % Check if elem loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewDistribLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            pb    = draw.size/125;           % load symbol size (pyramid base)
            clr   = [1,0,0];                 % load symbol color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            for e = 1:draw.mdl.nel
                if ((~isempty(draw.mdl.elems(e).load.uniformGbl))  &&...
                   (~all(draw.mdl.elems(e).load.uniformGbl == 0))) ||...
                   ((~isempty(draw.mdl.elems(e).load.linearGbl))   &&...
                   (~all(draw.mdl.elems(e).load.linearGbl == 0)))
               
                    % Get element length
                    L = draw.mdl.elems(e).length;
                    
                    % Get element orientation angle
                    ang = acos(abs(draw.mdl.elems(e).cosine_X));
                    
                    % Get element orientation angle cosine with axes X and Y
                    cx = draw.mdl.elems(e).cosine_X;
                    cy = draw.mdl.elems(e).cosine_Y;
                    
                    % Get element end nodes IDs
                    n1 = draw.mdl.elems(e).nodes(1).id;
                    n2 = draw.mdl.elems(e).nodes(2).id;
                    
                    % Get nodal coordinates
                    x1 = draw.mdl.nodes(n1).coord(1);
                    y1 = draw.mdl.nodes(n1).coord(2);
                    x2 = draw.mdl.nodes(n2).coord(1);
                    y2 = draw.mdl.nodes(n2).coord(2);
                    
                    % Calculate new element orientation angle and text
                    % rotation angle acording to element orientation
                    if (x1 <= x2) && (y1 <= y2)
                        alpha = ang;
                        rot   = ang * 180/pi;
                    elseif (x1 >= x2) && (y1 >= y2)
                        alpha = ang + pi;
                        rot   = ang * 180/pi;
                    elseif (x1 <= x2) && (y1 >= y2)
                        alpha = -ang;
                        rot   = -ang * 180/pi;
                    elseif (x1 >= x2) && (y1 <= y2)
                        alpha = pi - ang;
                        rot   = -ang * 180/pi;
                    end
                    
                    % Initialize load values on element ends
                    qi = 0;
                    qf = 0;
                    
                    % Add uniform load contribtuion
                    if isempty(draw.mdl.elems(e).load.uniformGbl) == 0
                        qz = draw.mdl.elems(e).load.uniformGbl(3);
                        qi = qi + qz;
                        qf = qf + qz;
                    end
                    
                    % Add linear load contribtuion
                    if isempty(draw.mdl.elems(e).load.linearGbl) == 0
                        qzi = draw.mdl.elems(e).load.linearGbl(3);
                        qzf = draw.mdl.elems(e).load.linearGbl(6);
                        qi = qi + qzi;
                        qf = qf + qzf;
                    end
                    
                    % Calculate load quation coefficients:
                    % q(x) = Ax + B
                    A = (qf - qi)/L;
                    B = qi;
                    
                    % Draw load symbol on a number cross-sections along element local axis X
                    step = L / round(20 * L/draw.size);
                    for x = 0:step:L
                        % Calculate current cross-section coordinates
                        xs = x1 + x * cx;
                        ys = y1 + x * cy;
                        
                        % Calculate load value on current cross-section
                        q = A * x + B;
                        
                        % Draw load symbol on current cross-section
                        if q > 0
                            dir = 'z+';
                        elseif q < 0
                            dir = 'z-';
                        end
                        
                        if q ~= 0
                            draw.arrow3DPlaneProj(xs,ys,pb,dir,clr,'drawElemLoads')
                            hold on
                        end
                    end
                    
                    % Plot line
                    plot([x1 x2],[y1 y2],':','color',clr,'LineWidth',1.35,'tag','drawElemLoads')

                    % Write load values:
                    % If load is uniform, draw load value in the middle of the element
                    if qi == qf
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN/m',dc,qi);
                        else
                            value = sprintf('%.*f',dc,qi);
                        end
                        [xm,ym] = draw.coordTransf2D(L/2,0,x1,y1,alpha);
                        text(xm,ym,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',qi);
                        
                    % If load is linear, draw initial and final load values
                    else
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value1 = sprintf('%.*f kN/m',dc,qi);
                            value2 = sprintf('%.*f kN/m',dc,qf);
                        else
                            value1 = sprintf('%.*f',dc,qi);
                            value2 = sprintf('%.*f',dc,qf);
                        end
                        [xm,ym] = draw.coordTransf2D(0,0,x1,y1,alpha);
                        [xn,yn] = draw.coordTransf2D(0,0,x2,y2,alpha);
                        if abs(alpha) <= pi/2
                            text(xm,ym,value1,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',qi);
                            text(xn,yn,value2,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',qf);
                        else
                            text(xm,ym,value1,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',qi);
                            text(xn,yn,value2,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'Rotation',rot,'tag','textElemLoads','UserData',qf);
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
            shift = draw.size/75;  % distance between load symbol and nodal point
            al    = draw.size/12;  % load symbol size (arrow lenght)
            ah    = draw.size/60;  % load symbol size (arrowhead height)
            ab    = draw.size/60;  % load symbol size (arrowhead base)
            clr   = [1,0,0];       % load color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                % Check if current node has a nodal load
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fz = draw.mdl.nodes(n).load.static(3);
                    mx = draw.mdl.nodes(n).load.static(4);
                    my = draw.mdl.nodes(n).load.static(5);
                    
                    % Draw load component in the Z direction
                    if fz > 0
                        draw.arrow3DPlaneProj(x,y,ab,'z+',clr,'drawNodalLoads')
                        hold on       
                    elseif fz < 0
                        draw.arrow3DPlaneProj(x,y,ab,'z-',clr,'drawNodalLoads')
                        hold on
                    end
                    if fz ~= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,fz);
                        else
                            value = sprintf('%.*f',dc,fz);
                        end
                        text(x+shift,y+1.2*shift,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',fz);
                    end
                    
                    % Draw moment component in the X direction 
                    if mx > 0
                        draw.arrow2D(x-shift,y,al,ah,ab,pi,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x-shift-0.25*al,y,0.75*al,ah,ab,pi,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x-shift-al,y,value,'HorizontalAlignment','right','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                        
                    elseif mx < 0
                        draw.arrow2D(x+shift,y,al,ah,ab,0,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x+shift+0.25*al,y,0.75*al,ah,ab,0,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x+shift+al,y,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    end
                    
                    % Draw moment component in the Y direction 
                    if my > 0
                        draw.arrow2D(x,y-shift,al,ah,ab,-pi/2,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x,y-shift-0.25*al,0.75*al,ah,ab,-pi/2,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y-shift-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textNodalMoments','UserData',abs(my));
                        
                    elseif my < 0
                        draw.arrow2D(x,y+shift,al,ah,ab,pi/2,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x,y+shift+0.25*al,0.75*al,ah,ab,pi/2,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y+shift+al,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
                    end
                end
            end
        end
        
         %------------------------------------------------------------------
        % Draws applied dynamic nodal loads and moments.
        function draw = dynamicNodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            shift = draw.size/75;  % distance between load symbol and nodal point
            al    = draw.size/12;  % load symbol size (arrow lenght)
            ah    = draw.size/60;  % load symbol size (arrowhead height)
            ab    = draw.size/60;  % load symbol size (arrowhead base)
            clr   = [0,0.7,0];     % load color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.dynamic) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fz = draw.mdl.nodes(n).load.dynamic(3);
                    mx = draw.mdl.nodes(n).load.dynamic(4);
                    my = draw.mdl.nodes(n).load.dynamic(5);
                    
                    % Draw load component in the Z direction
                    if fz > 0
                        draw.arrow3DPlaneProj(x,y,ab,'z+',clr,'drawNodalLoads')
                        hold on       
                    elseif fz < 0
                        draw.arrow3DPlaneProj(x,y,ab,'z-',clr,'drawNodalLoads')
                        hold on
                    end
                    if fz ~= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,fz);
                        else
                            value = sprintf('%.*f',dc,fz);
                        end
                        text(x+shift,y+1.2*shift,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalLoads','UserData',fz);
                    end
                    
                    % Draw moment component in the X direction
                    if mx > 0
                        draw.arrow2D(x-shift,y,al,ah,ab,pi,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x-shift-0.25*al,y,0.75*al,ah,ab,pi,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x-shift-al,y,value,'HorizontalAlignment','right','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                        
                    elseif mx < 0
                        draw.arrow2D(x+shift,y,al,ah,ab,0,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x+shift+0.25*al,y,0.75*al,ah,ab,0,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x+shift+al,y,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    end
                    
                    % Draw moment component in the Y direction
                    if my > 0
                        draw.arrow2D(x,y-shift,al,ah,ab,-pi/2,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x,y-shift-0.25*al,0.75*al,ah,ab,-pi/2,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y-shift-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textNodalMoments','UserData',abs(my));
                        
                    elseif my < 0
                        draw.arrow2D(x,y+shift,al,ah,ab,pi/2,clr,'drawNodalLoads');
                        hold on
                        draw.arrow2D(x,y+shift+0.25*al,0.75*al,ah,ab,pi/2,clr,'drawNodalLoads');
                        hold on
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y+shift+al,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
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
            r   = draw.size/60;            % concentrated mass symbol radius
            clr = [0,0,1];                 % mass color
            dc  = getappdata(0,'decPrec'); % decimal precision
            
             for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).displMass) == 0
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    mass = draw.mdl.nodes(n).displMass;
                                      
                    % Draw mass symbol
                    if mass > 0
                        s = draw.sphere(x,y,0,r,'drawNodalMass');
                        set(s,'Edgecolor',clr,'FaceColor',clr);
                        
                        % Write concentrated mass value
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kg',dc,abs(mass)*1000);
                        else
                            value = sprintf('%.*f',dc,abs(mass)*1000);
                        end
                        text(x+0.8*r,y-0.8*r,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',clr,'tag','textNodalMass','UserData',abs(mass)*1000);
                     end
                end
             end
        end
        
        %------------------------------------------------------------------
        % Draws nodal prescribed displacement representation.
        function draw = nodalPrescDispl(draw)
        end
        
        %------------------------------------------------------------------
        % Draws nodal initial condition values.
        function draw = nodalInitialConditions(draw)
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
            d       = draw.size/150;  % distance between temperature grad symbol and element
            heatClr = [1,0,0];        % heat color
            coldClr = [0,0,1];        % cold color
            dc      = getappdata(0,'decPrec'); % decimal precision
            
            for e = 1:draw.mdl.nel
                if draw.mdl.elems(e).load.tempVar_Z ~= 0
                    % Get nodal coordinates
                    x1 = draw.mdl.elems(e).nodes(1).coord(1);
                    y1 = draw.mdl.elems(e).nodes(1).coord(2);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    
                    % Get element orientation angle cosine with axes X and Y
                    cx = draw.mdl.elems(e).cosine_X;
                    cy = draw.mdl.elems(e).cosine_Y;
                    
                    % Get temperature variation values
                    dtz = draw.mdl.elems(e).load.tempVar_Z;
                    
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
                    
                    % Check if units are enabled
                    unitsAreOn = strcmp(get(mdata.unitsButton,'Checked'),'on');
                    
                    % Draw temperature variation symbols
                    if dtz > 0
                        line([x1-d*cy,x2-d*cy],[y1+d*cx,y2+d*cx],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        line([x1+d*cy,x2+d*cy],[y1-d*cx,y2-d*cx],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTz = %.*f oC',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        end
                    elseif dtz < 0
                        line([x1-d*cy,x2-d*cy],[y1+d*cx,y2+d*cx],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        line([x1+d*cy,x2+d*cy],[y1-d*cx,y2-d*cx],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTz = %.*f oC',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'Rotation',rot,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Computes deformed configuration scale factor.
        function draw = deformScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Computes dynamic deformed configuration scale factor.
        function draw = dynamicDeformScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws structure deformed configuration on a given scale.
        % Input arguments:
        %  scale: deformed configuration scale factor
        function draw = deformConfig(draw,~)
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
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function axialForceEnvelop(~,~)
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
        % Input arguments:
        %  scale: shear force diagram scale factor
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
        % Input arguments:
        %  scale: shear force diagram scale factor
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
        % Input arguments:
        %  scale: bending moment diagram scale factor
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
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function bendingMomentEnvelop_XZ(~,~)
        end
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        function draw = reactions(draw)
        end
        
        %------------------------------------------------------------------
        % Draws natural undamped free vibration of specified mode
        % Input arguments:
        %  nMode: vibration mode identifier
        %  scale: deformation scale factor
        function vibrationMode(~,~,~)
        end
        
        %------------------------------------------------------------------
        % Draws deformed configuration after dynamic analysis, on given
        % time.
        % Input arguments:
        %  step : time step identifier
        %  scale: deformation scale factor
        function dynamicDeform(~,~,~)
        end
    end
end