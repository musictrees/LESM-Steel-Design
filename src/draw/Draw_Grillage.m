%% Draw_Grillage class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *Grillage* draw object.
%
classdef Draw_Grillage < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Grillage(mdl)
            draw = draw@Draw(37.5,30);
            
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
                draw.size = 10;
            else
                draw.size = sz;
            end
        end
        
        %------------------------------------------------------------------
        % Sets axes limits according to model dimensions.
        function draw = setLimits(draw)
            mdata = guidata(findobj('Tag','GUI_Main'));
            ax = mdata.axes_Canvas;
            ax.Clipping = 'off';
            if draw.mdl.nnp == 0
                xmin = 0;
                xmax = 0;
                ymin = 0;
                ymax = 0;
                d = 5;
            elseif draw.mdl.nnp == 1
                xmin = draw.mdl.nodes(1).coord(1);
                xmax = draw.mdl.nodes(1).coord(1);
                ymin = draw.mdl.nodes(1).coord(2);
                ymax = draw.mdl.nodes(1).coord(2);
                d = 5;
            else
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
                d = draw.size/5;
            end
            xlim([xmin - d, xmax + d]);
            ylim([ymin - d, ymax + d]);
            zlim([0,((max(xlim)-min(xlim))+(max(ylim)-min(ylim)))/1000]);
            xlabel('X');
            ylabel('Y');
            zlabel(' ');
            zticks(0);
            zticklabels({' '});
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
            nm      = draw.size/230;  % node mark symbol (cube side)
            r       = draw.size/125;  % hinge symbol radius
            ph      = draw.size/40;   % translation constraint symbol (pyramid height)
            pb      = draw.size/80;   % translation constraint symbol (pyramid base)
            cs      = draw.size/75;   % rotation constraint symbol (cube side)
            sh      = draw.size/20;   % spring symbol height
            sr      = draw.size/125;  % rotational spring symbol radius
            nclr    = [0.0,0.0,0.0];  % node and hinge color
            sclr    = [0.6,0.6,0.6];  % support color
            sprclr  = [0.6,0.0,0.4];  % spring color
            rotsclr = [0.7,0.0,0.5];  % rotational spring color
            dc      = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Distance between translation constraint support symbol and nodal point
                shift = 0;
                
                if draw.mdl.nodes(n).ebc(4) == FREE_DOF || draw.mdl.nodes(n).ebc(5) == FREE_DOF
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng == tot && tot > 0
                        shift = r;
                    else
                        draw.cube(x,y,0,nm,nclr,'drawNodes');
                        hold on;
                    end
                elseif strcmp(drawSupports,'on')
                    if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF
                        shift = cs/2;
                    elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF
                        shift = sr;
                    end
                else
                    draw.cube(x,y,0,nm,nclr,'drawNodes');
                    hold on;
                end
                
                % Draw support conditions
                if ~strcmp(drawSupports,'on')
                    continue;
                end
                
                % Draw fixed support in Z direction
                if draw.mdl.nodes(n).ebc(3) == FIXED_DOF
                    draw.pyramid(x,y,-shift,ph,pb,'z+',sclr,'drawSupports');
                    hold on;
                    
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
                    draw.SpringZ_3D(x,y,-shift,sh,sprclr,'drawSupports');
                    hold on;
                    text(x+sh/10,y+sh/10,-shift-sh/1.5,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSprings','UserData',kz);
                end
                
                % Draw fixed rotation support
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF
                    draw.cube(x,y,0,cs,sclr,'drawSupports');
                    hold on;
                
                % Draw spring rotation support
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF
                    kr = draw.mdl.nodes(n).springStiff(4);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        if kr >= 1000
                            value = sprintf('%.*e kNm/rad',dc,kr);
                        else
                            value = sprintf('%.*f kNm/rad',dc,kr);
                        end
                    else
                        if kr >= 1000
                            value = sprintf('%.*e',dc,kr);
                        else
                            value = sprintf('%.*f',dc,kr);
                        end
                    end
                    draw.rotSpring_3D(x,y,0,sr,rotsclr,'drawSupports');
                    hold on;
                    text(x,y,1.1*shift,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textRotSprings','UserData',kr);
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
            nm     = draw.size/230;            % node mark symbol (cube side)
            r      = draw.size/125;            % hinge symbol radius
            sr     = draw.size/125;            % semi-rigid joint symbol radius
            sprclr = [0.6,0,0.4];              % spring color
            clr    = [0,0,0];                  % element color
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
                    draw.sphere(x1,y1,0,r,'drawElements');
                    hold on;
                    xi = x1 + r * cx;
                    yi = y1 + r * cy;
                elseif draw.mdl.elems(e).hingei == 0 % Hinge on element end
                    draw.sphere(x1+(r+nm/2)*cx,y1+(r+nm/2)*cy,0,r,'drawElements');
                    xi = x1 + 2 * r * cx;
                    yi = y1 + 2 * r * cy;
                elseif draw.mdl.elems(e).hingei == 2
                    xi = x1 + (sr + nm/2) * cx;
                    yi = y1 + (sr + nm/2) * cy;
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),dir(3,:)] = draw.mdl.elems(e).locAxis;
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
                        if krxi == kryi
                            text(xi,yi,sr,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krxi);
                        else
                            text(xi+1.6*sr*cx,yi+1.6*sr*cy,sr,value_x,'HorizontalAlignment','center','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krxi);
                            text(xi-1.6*sr*cy,yi+1.6*sr*cx,sr,value_y,'HorizontalAlignment','center','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',kryi);
                        end
                    else
                        % Connect element end coordinates
                        X_srj = [xi,x1];
                        Y_srj = [yi,y1];
                        Z_srj = [0,0];
                        plot3(X_srj,Y_srj,Z_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Set element final coordinates by checking if there is a hinge on nodal point position or on element end
                if hef == totf && (draw.mdl.nodes(n2).ebc(4) == 0 || draw.mdl.nodes(n2).ebc(5) == 0) % Hinge on node
                    draw.sphere(x2,y2,0,r,'drawElements');
                    hold on;
                    xf = x2 - r * cx;
                    yf = y2 - r * cy;
                elseif draw.mdl.elems(e).hingef == 0 % Hinge on element end
                    draw.sphere(x2-(r+nm/2)*cx,y2-(r+nm/2)*cy,0,r,'drawElements');
                    xf = x2 - 2 * r * cx;
                    yf = y2 - 2 * r * cy;
                elseif draw.mdl.elems(e).hingef == 2
                    xf = x2 - (sr + nm/2) * cx;
                    yf = y2 - (sr + nm/2) * cy;
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),dir(3,:)] = draw.mdl.elems(e).locAxis;
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
                        if krxf == kryf
                            text(xf,yf,sr,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krxf);
                        else
                            text(xf+1.6*sr*cx,yf+1.6*sr*cy,sr,value_x,'HorizontalAlignment','center','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',krxf);
                            text(xf-1.6*sr*cy,yf+1.6*sr*cx,sr,value_y,'HorizontalAlignment','center','VerticalAlignment','middle','Color',sprclr,'Fontsize',8.5,'tag','textSemiRigid','UserData',kryf);
                        end
                    else
                        % Connect element end coordinates
                        X_srj = [xf,x2];
                        Y_srj = [yf,y2];
                        Z_srj = [0,0];
                        plot3(X_srj,Y_srj,Z_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Connect element end coordinates
                X = [xi,xf];
                Y = [yi,yf];
                Z = [0,0];
                plot3(X,Y,Z,'Color',clr,'tag','drawElements');
                hold on;
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
            z = coords(3);
            
            % Semi-rigid joint symbol radius
            r = sz;
            
            % local axis
            dx = dir;
            dy = [dir(2,:);dir(3,:);dir(1,:)];

            % arrow properties
            l = 4 * r;
            h = 0.8 * r;
            B = 0.5 * r;
            
            % Draw double arrows to indicate the direction of local axis
            srjclr = [0.6,0,0.4];
            aux = [x,y,z] + l * dir(1,:);
            draw.moment3D(draw,aux(1),aux(2),aux(3),l,h,B,dx,srjclr,'drawSemiRigid',true);
            hold on
            aux = [x,y,z] + l * dir(2,:);
            draw.moment3D(draw,aux(1),aux(2),aux(3),l,h,B,dy,srjclr,'drawSemiRigid',true);
            hold on
            
            % Draw sphere to represent semi-rigd joint
            [a, b, c] = sphere;
            s = surf(a * r + x, b * r + y, c * r + z);
            set(s, 'Edgecolor', clr,'FaceColor', [1,1,1], 'tag', 'drawSemiRigid');
        end
        
        %------------------------------------------------------------------
        % Computes element loads scale factor.
        function scl = elemLoadsScaleFactor(draw)
            max_elem = zeros(1,draw.mdl.nel);
            
            for e = 1:draw.mdl.nel
                % Initialize load values on element ends
                qi = 0;
                qf = 0;
                
                % Add uniform load contribtuion
                if ~isempty(draw.mdl.elems(e).load.uniformGbl)
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
                
                % Get maximum load value on current element
                max_elem(e) = max(abs(qi),abs(qf));
            end
            
            % Get maximum load value on model
            max_val = max(max_elem);
            
            % Calculate scale factor
            if isempty(max_val) || max_val == 0
                scl = 0;
            else
                scl = draw.size/(12*max_val);
            end
            
            setappdata(0,'load_sf',scl);
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
            shift = draw.size/130;            % distance between load symbol and element
            ph    = draw.size/200;            % load symbol size (pyramid height)
            pb    = draw.size/300;            % load symbol size (pyramid base)
            clr   = [1,0,0];                  % load symbol color
            scl   = getappdata(0,'load_sf');  % Load symbol scale
            dc    = getappdata(0,'decPrec');  % decimal precision
            
            for e = 1:draw.mdl.nel
                if ((~isempty(draw.mdl.elems(e).load.uniformGbl))  &&...
                   (~all(draw.mdl.elems(e).load.uniformGbl == 0))) ||...
                   ((~isempty(draw.mdl.elems(e).load.linearGbl))   &&...
                   (~all(draw.mdl.elems(e).load.linearGbl == 0)))
               
                    % Get element length
                    L = draw.mdl.elems(e).length;
                    
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
                    
                    % Initialize load values on element ends
                    qi = 0;
                    qf = 0;
                    
                    % Add uniform load contribtuion
                    if isempty(draw.mdl.elems(e).load.uniformGbl) == 0
                        qz = draw.mdl.elems(e).load.uniformGbl(3);
                        qi = qz;
                        qf = qz;
                    end
                    
                    % Add linear load contribtuion
                    if isempty(draw.mdl.elems(e).load.linearGbl) == 0
                        qi = qi + draw.mdl.elems(e).load.linearGbl(3);
                        qf = qf + draw.mdl.elems(e).load.linearGbl(6);
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
                        if abs(scl * q) >= ph
                            if q > 0
                                draw.arrow3D(draw,xs,ys,-shift,scl*q,ph,pb,'z+',clr,'drawElemLoads');
                            elseif q < 0
                                draw.arrow3D(draw,xs,ys,shift,-scl*q,ph,pb,'z-',clr,'drawElemLoads');
                            end
                        end
                    end
                    
                    % Connect load symbols
                    if  (qi < 0) && (qf > 0)
                        x0 = (abs(qi)*(L))/(abs(qf)+abs(qi)); 
                        [xu,yu,zu] = draw.coordTransf3D(x0,0,shift,x1,y1,0,e);
                        [xd,yd,zd] = draw.coordTransf3D(x0,0,-shift,x1,y1,0,e);
                        X = [x1,xu,xd,x2];
                        Y = [y1,yu,yd,y2];
                        Z = [shift-scl*qi,zu,zd,-shift-scl*qf];
                        line(X,Y,Z,'Color',clr,'tag','drawElemLoads');
                    elseif (qi > 0) && (qf < 0)
                        x0 = (abs(qi)*(L))/(abs(qf)+abs(qi));
                        [xu,yu,zu] = draw.coordTransf3D(x0,0,shift,x1,y1,0,e);
                        [xd,yd,zd] = draw.coordTransf3D(x0,0,-shift,x1,y1,0,e);
                        X = [x1,xd,xu,x2];
                        Y = [y1,yd,yu,y2];
                        Z = [-shift-scl*qi,zd,zu,shift-scl*qf];
                        line(X,Y,Z,'Color',clr,'tag','drawElemLoads');
                    elseif (qi~=0) || (qf~=0)
                        X = [x1,x2];
                        Y = [y1,y2];
                        if qi >= 0 && qf >= 0
                            Z = [-shift-scl*qi,-shift-scl*qf];
                        elseif qi <= 0 && qf <= 0
                            Z = [shift-scl*qi,shift-scl*qf];
                        end
                        line(X,Y,Z,'Color',clr,'tag','drawElemLoads');
                        if qi == 0 && qf > 0
                            line([x1,x1],[y1,y1],[0,-shift],'Color',clr,'tag','drawElemLoads');
                        elseif qi == 0 && qf < 0
                            line([x1,x1],[y1,y1],[0,shift],'Color',clr,'tag','drawElemLoads');
                        end
                        if qf == 0 && qi > 0
                            line([x2,x2],[y2,y2],[0,-shift],'Color',clr,'tag','drawElemLoads');
                        elseif qf == 0 && qi < 0
                            line([x2,x2],[y2,y2],[0,shift],'Color',clr,'tag','drawElemLoads');
                        end
                    end
                    
                    % Write load values
                    if qi == qf % If load is uniform, draw load value in the middle of the element
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kN/m',dc,abs(qi));
                        else
                            value = sprintf('%.*f',dc,abs(qi));
                        end
                        if qi > 0
                            text((x1+x2)/2,(y1+y2)/2,-shift-scl*qi,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textElemLoads','UserData',abs(qi));
                        elseif qi < 0
                            text((x1+x2)/2,(y1+y2)/2,shift-scl*qi,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textElemLoads','UserData',abs(qi));
                        end
                        
                    else % If load is linear, draw initial and final load values
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value1 = sprintf('%.*f kN/m',dc,abs(qi));
                            value2 = sprintf('%.*f kN/m',dc,abs(qf));
                        else
                            value1 = sprintf('%.*f',dc,abs(qi));
                            value2 = sprintf('%.*f',dc,abs(qf));
                        end
                        if qi > 0
                            text(x1,y1,-shift-scl*qi,value1,'HorizontalAlignment','left','VerticalAlignment','top','Color',clr,'tag','textElemLoads','UserData',abs(qi));
                        elseif qi < 0
                            text(x1,y1,shift-scl*qi,value1,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',clr,'tag','textElemLoads','UserData',abs(qi));
                        end
                        if qf > 0
                            text(x2,y2,-shift-scl*qf,value2,'HorizontalAlignment','right','VerticalAlignment','top','Color',clr,'tag','textElemLoads','UserData',abs(qf));
                        elseif qf < 0
                            text(x2,y2,shift-scl*qf,value2,'HorizontalAlignment','right','VerticalAlignment','bottom','Color',clr,'tag','textElemLoads','UserData',abs(qf));
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
            shift = draw.size/150;            % distance between load symbol and nodal point
            al    = draw.size/17;             % load symbol size (arrow lenght)
            ah    = draw.size/60;             % load symbol size (arrowhead height)
            ab    = draw.size/130;            % load symbol size (arrowhead base)
            clr   = [1,0,0];                  % load color
            dc    = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fz = draw.mdl.nodes(n).load.static(3);
                    mx = draw.mdl.nodes(n).load.static(4);
                    my = draw.mdl.nodes(n).load.static(5);
                    
                    % Z load
                    if fz > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        draw.arrow3D(draw,x,y,-shift,al,ah,ab,'z+',clr,'drawNodalLoads');
                        text(x,y,-shift-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                    elseif fz < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        draw.arrow3D(draw,x,y,shift,al,ah,ab,'z-',clr,'drawNodalLoads');
                        text(x,y,shift+al,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                    end
                    
                    % X moment
                    if mx > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        draw.moment3D(draw,x-shift,y,0,al,ah,ab,'x+',clr,'drawNodalLoads');
                        text(x-shift-al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    elseif mx < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        draw.moment3D(draw,x+shift,y,0,al,ah,ab,'x-',clr,'drawNodalLoads');
                        text(x+shift+al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    end
                    
                    % Y moment
                    if my > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        draw.moment3D(draw,x,y-shift,0,al,ah,ab,'y+',clr,'drawNodalLoads');
                        text(x,y-shift-al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
                    elseif my < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        draw.moment3D(draw,x,y+shift,0,al,ah,ab,'y-',clr,'drawNodalLoads');
                        text(x,y+shift+al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
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
            shift = draw.size/150;            % distance between load symbol and nodal point
            al    = draw.size/16;             % load symbol size (arrow lenght)
            ah    = draw.size/60;             % load symbol size (arrowhead height)
            ab    = draw.size/110;            % load symbol size (arrowhead base)
            clr   = [0,0.7,0];                % load color
            dc    = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.dynamic) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get nodal load components
                    fz = draw.mdl.nodes(n).load.dynamic(3);
                    mx = draw.mdl.nodes(n).load.dynamic(4);
                    my = draw.mdl.nodes(n).load.dynamic(5);
                    
                    % Z load
                    if fz > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        draw.arrow3D(draw,x,y,-shift,al,ah,ab,'z+',clr,'drawNodalLoads');
                        text(x,y,-shift-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                    elseif fz < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        draw.arrow3D(draw,x,y,shift,al,ah,ab,'z-',clr,'drawNodalLoads');
                        text(x,y,shift+al,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                    end
                    
                    % X moment
                    if mx > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        draw.moment3D(draw,x-shift,y,0,al,ah,ab,'x+',clr,'drawNodalLoads');
                        text(x-shift-al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    elseif mx < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        draw.moment3D(draw,x+shift,y,0,al,ah,ab,'x-',clr,'drawNodalLoads');
                        text(x+shift+al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(mx));
                    end
                    
                    % Y moment
                    if my > 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        draw.moment3D(draw,x,y-shift,0,al,ah,ab,'y+',clr,'drawNodalLoads');
                        text(x,y-shift-al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
                    elseif my < 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        draw.moment3D(draw,x,y+shift,0,al,ah,ab,'y-',clr,'drawNodalLoads');
                        text(x,y+shift+al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMoments','UserData',abs(my));
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
            r   = draw.size/100;            % concentrated mass symbol radius
            clr = [0,0,1];                  % mass color
            dc  = getappdata(0,'decPrec');  % decimal precision
            
             for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).displMass)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    mass = draw.mdl.nodes(n).displMass;
                    
                    % Draw horizontal load component
                    if mass > 0
                        s = draw.sphere(x,y,0,r,'drawNodalMass');
                        set(s,'Edgecolor',clr,'FaceColor',clr);
                        
                        % Write mass value
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value = sprintf('%.*f kg',dc,abs(mass)*1000);
                        else
                            value = sprintf('%.*f',dc,abs(mass)*1000);
                        end
                        text(x,y,1.2*r,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalMass','UserData',abs(mass)*1000);
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
            al  = draw.size/17;             % presc. displ. symbol size (arrow length)
            ah  = draw.size/60;             % presc. displ. symbol size (pyramid height)
            ab  = draw.size/130;            % presc. displ. symbol size (pyramid base)
            clr = [1,0,1];                  % presc. displ. symbol color
            dc  = getappdata(0,'decPrec');  % decimal precision
            
            % translation constraint symbol (pyramid height)
            if strcmp(drawSupports,'on')
                ph = draw.size/40;
            else
                ph = draw.size/150;
            end
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).prescDispl) == 0
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    
                    % Get prescribed displacement component values (convert to mm and rad)
                    rx = draw.mdl.nodes(n).prescDispl(4);
                    ry = draw.mdl.nodes(n).prescDispl(5);
                    dz = 1000 * draw.mdl.nodes(n).prescDispl(3);
                    
                    % Check if rotation is really fixed and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(4) == 1 && draw.mdl.nodes(n).ebc(5) == 1
                        cs = draw.size/150;
                        if rx ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on')
                                value = sprintf('%.*f rad',dc,abs(rx));
                            else
                                value = sprintf('%.*f',dc,abs(rx));
                            end
                            if rx > 0
                                draw.moment3D(draw,x-cs,y,0,al,ah,ab,'x+',clr,'drawPrescDispl');
                                text(x-cs-al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescRot','UserData',abs(rx));
                            else
                                draw.moment3D(draw,x+cs,y,0,al,ah,ab,'x-',clr,'drawPrescDispl');
                                text(x+cs+al,y,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescRot','UserData',abs(rx));
                            end
                        end
                        if ry ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on')
                                value = sprintf('%.*f rad',dc,abs(ry));
                            else
                                value = sprintf('%.*f',dc,abs(ry));
                            end
                            if ry > 0
                                draw.moment3D(draw,x,y-cs,0,al,ah,ab,'y+',clr,'drawPrescDispl');
                                text(x,y-cs-al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescRot','UserData',abs(ry));
                            else
                                draw.moment3D(draw,x,y+cs,0,al,ah,ab,'y-',clr,'drawPrescDispl');
                                text(x,y+cs+al,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescRot','UserData',abs(ry));
                            end
                        end
                    elseif draw.mdl.nodes(n).ebc(4) == 2 && draw.mdl.nodes(n).ebc(5) == 2
                        cs = draw.size/125; % rotational spring symbol radius
                    else
                        [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                        if hng == tot && tot > 0
                            cs = draw.size/125; % hinge symbol radius
                        else
                            cs = 0;
                        end
                    end
                    
                    % Check if translation in the Z axis is really fixed and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(3) == 1 && dz ~= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dz));
                        else
                            value = sprintf('%.*f',dc,abs(dz));
                        end
                        if dz > 0
                            draw.arrow3D(draw,x,y,-cs-ph,al,ah,ab,'z+',clr,'drawPrescDispl');
                        elseif dz < 0
                            draw.arrow3D(draw,x,y,-cs-ph-al,al,ah,ab,'z-',clr,'drawPrescDispl');
                        end
                        text(x,y,-cs-ph-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
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
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            r   = draw.size/125;   % hinge symbol radius
            m   = draw.size/100;   % concentrated mass symbol radius
            clr = [0,0.7,0];       % color
            
            for n = 1:draw.mdl.nnp
                flag = 0;
                str = '';
                if isempty(draw.mdl.nodes(n).initCond)
                    continue;
                end
                if draw.mdl.nodes(n).initCond(3,1) ~= 0 || draw.mdl.nodes(n).initCond(4,1) ~= 0 || draw.mdl.nodes(n).initCond(5,1) ~= 0
                    flag = flag + 1;
                end
                if draw.mdl.nodes(n).initCond(3,2) ~= 0 || draw.mdl.nodes(n).initCond(4,2) ~= 0 || draw.mdl.nodes(n).initCond(5,2) ~= 0
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
                    text(x,y,2*m,str,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textInitialConditions');
                elseif draw.mdl.nodes(n).ebc(4) == 2 && draw.mdl.nodes(n).ebc(5) == 2 && strcmp(drawSupports,'on')
                    text(x,y,2*r,str,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textInitialConditions');
                else
                    text(x,y,1.2*r,str,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textInitialConditions');
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
            d       = draw.size/150;            % distance between temperature grad symbol and element
            heatClr = [1,0,0];                  % heat color
            coldClr = [0,0,1];                  % cold color
            dc      = getappdata(0,'decPrec');  % decimal precision
            
            for e = 1:draw.mdl.nel
                if draw.mdl.elems(e).load.tempVar_Z ~= 0
                    x1 = draw.mdl.elems(e).nodes(1).coord(1);
                    y1 = draw.mdl.elems(e).nodes(1).coord(2);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    
                    % Get temperature variation values
                    dtz = draw.mdl.elems(e).load.tempVar_Z;
                    
                    % Check if units are enabled
                    unitsAreOn = strcmp(get(mdata.unitsButton,'Checked'),'on');
                    
                    % Draw temperature variation symbols
                    if dtz > 0
                        line([x1,x2],[y1,y2],[d,d],'Color',coldClr,'LineStyle','-.', 'tag','drawThermalLoads');
                        line([x1,x2],[y1,y2],[-d,-d],'Color',heatClr,'LineStyle','-.', 'tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,-d,sprintf('dTz = %.*f ºC',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,-d,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        end
                    else
                        line([x1,x2],[y1,y2],[d,d],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        line([x1,x2],[y1,y2],[-d,-d],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,-d,sprintf('dTz = %.*f ºC',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,-d,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots ID number of nodes.
        function draw = nodeID(draw)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            anl   = get(mdata.popupmenu_AnalysisType,'Value');
            tz    = draw.size/120;  % distance between text and nodal point (axis Z)
            
            for n = 1:draw.mdl.nnp
                x  = draw.mdl.nodes(n).coord(1);
                y  = draw.mdl.nodes(n).coord(2);
                id = sprintf('%d',n);
                
                if anl == 2                             &&...
                  ~isempty(draw.mdl.nodes(n).displMass) &&...
                   draw.mdl.nodes(n).displMass > 0      &&...
                   strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                   text(x,y,2*tz,id,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','tag','textNodeID');
                elseif draw.mdl.nodes(n).ebc(4) ~= 0 || draw.mdl.nodes(n).ebc(5) ~= 0 && get(mdata.viewSupportsButton,'Checked')
                    text(x,y,1.5*tz,id,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','tag','textNodeID');
                else
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng == tot && tot > 0
                        text(x,y,1.5*tz,id,'FontWeight','bold','tag','textNodeID');
                    else
                        text(x,y,tz,id,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','tag','textNodeID');
                    end
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
                text((x1+x2)/2,(y1+y2)/2,0,id,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold','tag','textElemID');
            end
        end
        
        %------------------------------------------------------------------
        % Draws element orientation indication from inital to final node.
        function draw = elementOrientation(draw)
            clr = [0,.7,0]; % orientation symbol color
            for e = 1:draw.mdl.nel
                % Calculate spear length
                l = draw.size/20;
                
                % Get nodal coordinates
                xi = draw.mdl.elems(e).nodes(1).coord(1);
                yi = draw.mdl.elems(e).nodes(1).coord(2);
                xf = draw.mdl.elems(e).nodes(2).coord(1);
                yf = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Calculate element local axis X orientation vector
                x = [xf-xi,yf-yi,0];
                x = l*x/norm(x);
                
                % Get orientation vector on the local XZ plane
                z = [draw.mdl.elems(e).vz(1),draw.mdl.elems(e).vz(2),draw.mdl.elems(e).vz(3)];
                
                % Calculate element local axis Y orientation vector
                y = cross(z,x);
                y = l*y/norm(y);
                
                % Calculate element local axis Z orientation vector
                z = cross(x,y);
                z = l*z/norm(z);
                
                % Draw orientation symbol
                xm = (xi+xf)/2;
                ym = (yi+yf)/2;
                
                X = [xm,xm+x(1)];
                Y = [ym,ym+x(2)];
				Z = [0,x(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.5,'tag','drawElemOrient');
                text(xm+x(1),ym+x(2),x(3),'X','HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'FontWeight','bold','tag','drawElemOrient');
                
                X = [xm,xm+y(1)];
                Y = [ym,ym+y(2)];
				Z = [0,y(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.5,'tag','drawElemOrient');
                text(xm+y(1),ym+y(2),y(3),'Y','HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'FontWeight','bold','tag','drawElemOrient');
                
                X = [xm,xm+z(1)];
                Y = [ym,ym+z(2)];
				Z = [0,z(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.5,'tag','drawElemOrient');
                text(xm+z(1),ym+z(2),z(3),'Z','HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'FontWeight','bold','tag','drawElemOrient');
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
                dz1 = draw.mdl.D(draw.mdl.ID(3,n1));
                dz2 = draw.mdl.D(draw.mdl.ID(3,n2));
                nodeDispl = [dz1,dz2];
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
                dsf = draw.size/(5*sliderm*max_disp);
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
                ids = [ draw.mdl.ID(3,n1) ;
                        draw.mdl.ID(3,n2) ];
                       
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
                % Get 50 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords;
                
                % Get element transversal internal displacements in local system
                dl = draw.mdl.elems(e).intDispl;
                
                % Deformed configuration global coordinates
                dfg = coords(3,:) + scale * dl(2,:);
                
                % Plot deformed configuration
                line(coords(1,:),coords(2,:),dfg,'Color',clr,'tag','drawDeformConfig');
            end
        end
        
        %------------------------------------------------------------------
        % Computes axial force diagram scale factor value.
        function draw = axialScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force diagram on a given scale.
        function draw = axialForce(draw,~)
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        function axialForceEnvelop(~,~)
        end
        
        %------------------------------------------------------------------
        % Computes torsion force diagram scale factor value (for envelop).
        function draw = torsionScaleFactor(draw)
        end
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment diagram.
        function draw = torsionMoment(draw)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element internal tortion moment value (always uniform
                % for grillage models) and convert it to string
                T = -draw.mdl.elems(e).torsion_moment(1);
                if T ~= 0
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kNm',dc,T);
                    else
                        value = sprintf('%+.*f',dc,T);
                    end
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%.*f kNm',dc,abs(T));
                    else
                        value = sprintf('%.*f',dc,abs(T));
                    end
                    T = abs(T);
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                
                % Write tortion moment value above the element
                text((x1+x2)/2,(y1+y2)/2,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textTorsionDiagram','UserData',T);
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment envelop diagram.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function torsionMomentEnvelop(draw,~)
            % Parameters
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
                
                % Get internal forces envelop values
                Tmax = draw.mdl.elems(e).intForcesEnvelop(1,:,3);
                Tmin = draw.mdl.elems(e).intForcesEnvelop(2,:,3);
                
                % Plot text containing values
                if strcmp(get(mdata.unitsButton,'Checked'),'on')
                    value_max = sprintf('%+.*f kNm',dc,Tmax(1));
                else
                    value_max = sprintf('%+.*f',dc,Tmax(1));
                end
                if strcmp(get(mdata.unitsButton,'Checked'),'on')
                    value_min = sprintf('%+.*f kNm',dc,Tmin(1));
                else
                    value_min = sprintf('%+.*f',dc,Tmin(1));
                end
                text((x1+x2)/2,(y1+y2)/2,0,value_max,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textTorsionDiagram','UserData',Tmax(1));
                text((x1+x2)/2,(y1+y2)/2,0,value_min,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textTorsionDiagram','UserData',Tmin(1));
            end
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
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Get maximum internal shear force value of each element
            max_elem = zeros(1,draw.mdl.nel);
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get maximum value at element ends
                Q1_Z = draw.mdl.elems(e).shear_force_Z(:,1);
                Q2_Z = draw.mdl.elems(e).shear_force_Z(:,2);
                max_end = max(vertcat(abs(Q1_Z),abs(Q2_Z)));
                
                % Get maximum internal value
                if ~isempty(draw.mdl.elems(e).maxShearForce_XZ)
                    max_int = abs(draw.mdl.elems(e).maxShearForce_XZ(1));
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
            setappdata(0,'shearXZ_sf',ssf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XZ plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function draw = shearForce_XZ(draw,scale)
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
                Q1 = draw.mdl.elems(e).shear_force_Z(1);
                Q2 = -draw.mdl.elems(e).shear_force_Z(2);
                
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
                
                % Calculate diagram beggining coordinates
                zd1 = scale * Q1;
                zd2 = scale * Q2;
                
                % Draw diagram extremities
                Xi = [x1,x1];
                Yi = [y1,y1];
                Zi = [0,zd1];
                Xf = [x2,x2];
                Yf = [y2,y2];
                Zf = [0,zd2];
                line(Xi,Yi,Zi,'Color',clr,'tag','drawShearForceXZDiagram');
                line(Xf,Yf,Zf,'Color',clr,'tag','drawShearForceXZDiagram');
                
                % Write force values
                if abs(Q1-Q2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,Q1);
                    else
                        value = sprintf('%+.*f',dc,Q1);
                    end
                    if Q1 >= 0
                        text((x1+x2)/2,(y1+y2)/2,zd1,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                    else
                        text((x1+x2)/2,(y1+y2)/2,zd1,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                    end
                else
                    if Q1 >= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value1 = sprintf('%+.*f kN',dc,Q1);
                        else
                            value1 = sprintf('%+.*f',dc,Q1);
                        end
                        text(x1,y1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                    else
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value1 = sprintf('%+.*f kN',dc,Q1);
                        else
                            value1 = sprintf('%+.*f',dc,Q1);
                        end
                        text(x1,y1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                    end
                    if Q2 >= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value2 = sprintf('%+.*f kN',dc,Q2);
                        else
                            value2 = sprintf('%+.*f',dc,Q2);
                        end
                        text(x2,y2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Q2);
                    else
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value2 = sprintf('%+.*f kN',dc,Q2);
                        else
                            value2 = sprintf('%+.*f',dc,Q2);
                        end
                        text(x2,y2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Q2);
                    end
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    Q = draw.mdl.elems(e).intStresses(1,:);
                    
                    % Avoid plotting numeric garbage
                    if ~all(abs(Q) < 10e-10)
                        % Get element 50 point division coords and internal stress values
                        coords = draw.mdl.elems(e).intCoords;
                        Qmax = draw.mdl.elems(e).maxShearForce_XZ;
                        
                        % Get element basis transformation matrix
                        rot = draw.mdl.elems(e).T;
                        
                        % Compute diagram coordinates
                        diagramCoords = rot' * [zeros(1,size(Q,2)); zeros(1,size(Q,2)); scale*Q] + coords;
                        
                        % Plot diagram
                        X = diagramCoords(1,:);
                        Y = diagramCoords(2,:);
                        Z = diagramCoords(3,:);
                        line(X,Y,Z,'Color',clr,'tag','drawShearForceXZDiagram');
                        
                        % Check if there is a maximum value within the diagram
                        if ~isempty(Qmax)
                            % Compute maximum stress value position, in global coordinates
                            QmaxGblCoords = [x1 + Qmax(2) * cx;
                                             y1 + Qmax(2) * cy;
                                             0];
                            
                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0;0;scale*Qmax(1)] + QmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            zp = maxPointCoords(3);
                            scatter3(xp,yp,zp,50,clr,'.','tag','drawShearForceXZDiagram')
                            
                            % Plot maximum value (avoid plotting text too close to element end)
                            if abs(Qmax(2) - L) >= L/25 && Qmax(2) >= L/25
                                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                    value = sprintf('%+.*f kN',dc,Qmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Qmax(1));
                                end
                                if zp >= 0
                                    text(xp,yp,zp,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                                else
                                    text(xp,yp,zp,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                                end
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1,xd2];
                        Y = [yd1,yd2];
                        Z = [zd1,zd2];
                        line(X,Y,Z,'Color',clr,'tag','drawShearForceXZDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [x1,x2];
                    Y = [y1,y2];
                    Z = [zd1,zd2];
                    line(X,Y,Z,'Color',clr,'tag','drawShearForceXZDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XZ plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function shearForceEnvelop_XZ(draw,scale)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Qmax = draw.mdl.elems(e).intForcesEnvelop(1,:,1);
                Qmin = draw.mdl.elems(e).intForcesEnvelop(2,:,1);
                
                if all(abs(Qmax) < 10e-10) && all(abs(Qmin) < 10e-10)
                    x1 = draw.mdl.elems(e).nodes(1).coord(1);
                    y1 = draw.mdl.elems(e).nodes(1).coord(2);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[0,0],'Color',clr,'tag','drawShearForceXZDiagram');
                    text((x1+x2)/2,(y1+y2)/2,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',0);
                    text((x1+x2)/2,(y1+y2)/2,0,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',0);
                else
                    % Get coords of 27 points along element
                    coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                    
                    % Compute diagram coordinates
                    maxCoords = [zeros(1,size(Qmax,2)); zeros(1,size(Qmax,2)); scale*Qmax] + coords;
                    minCoords = [zeros(1,size(Qmin,2)); zeros(1,size(Qmin,2)); scale*Qmin] + coords;
                    
                    % Plot envelop
                    Xmax = maxCoords(1,:);
                    Ymax = maxCoords(2,:);
                    Zmax = maxCoords(3,:);
                    Xmin = minCoords(1,:);
                    Ymin = minCoords(2,:);
                    Zmin = minCoords(3,:);
                    line(Xmax,Ymax,Zmax,'Color',clr,'tag','drawShearForceXZDiagram');
                    line(Xmin,Ymin,Zmin,'Color',clr,'tag','drawShearForceXZDiagram');
                    
                    % Plot lines to mark 6 internal points along envelop diagram
                    for i = plotValId
                        xx = [ minCoords(1,i) maxCoords(1,i) ];
                        yy = [ minCoords(2,i) maxCoords(2,i) ];
                        zz = [ minCoords(3,i) maxCoords(3,i) ];
                        line(xx, yy, zz, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                    end
                    
                    % Plot text containing values
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
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
                    if Qmax(1) >= 0
                        text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Qmax1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                    else
                        text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Qmax1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                    end
                    if Qmin(1) >= 0
                        text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Qmin1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(1));
                    else
                        text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Qmin1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(1));
                    end
                    if Qmax(end) >= 0
                        text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Qmax2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(end));
                    else
                        text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Qmax2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(end));
                    end
                    if Qmin(end) >= 0
                        text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Qmin2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(end));
                    else
                        text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Qmin2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(end));
                    end
                end
            end
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
            mdata = guidata(findobj('Tag','GUI_Main'));
            sliderm = get(mdata.slider_Scale,'Max');
            
            % Get maximum internal bending moment value of each element
            max_elem = zeros(1,draw.mdl.nel);
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get maximum value at element ends
                M1_Y = draw.mdl.elems(e).bending_moment_Y(:,1);
                M2_Y = draw.mdl.elems(e).bending_moment_Y(:,2);
                max_end = max(vertcat(abs(M1_Y),abs(M2_Y)));
                
                % Get maximum internal value
                if ~isempty(draw.mdl.elems(e).maxBendMoment_XZ)
                    max_int = max(abs(draw.mdl.elems(e).maxBendMoment_XZ(:,1)));
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
            setappdata(0,'bendingXZ_sf',bsf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function draw = bendingMoment_XZ(draw,scale)
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
                
                % Get element internal force values ate both ends
                M1 = draw.mdl.elems(e).bending_moment_Y(1);
                M2 = -draw.mdl.elems(e).bending_moment_Y(2);
                
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
                
                % Calculate diagram beginning coordinates
                zd1 = -scale * M1;
                zd2 = -scale * M2;
                
                % Draw diagram extremities
                Xi = [x1,x1];
                Yi = [y1,y1];
                Zi = [0,zd1];
                Xf = [x2,x2];
                Yf = [y2,y2];
                Zf = [0,zd2];
                line(Xi,Yi,Zi,'Color',clr,'tag','drawBendMomentXZDiagram');
                line(Xf,Yf,Zf,'Color',clr,'tag','drawBendMomentXZDiagram');
                
                % Write force values
                if abs(M1-M2) < 10e-10 && isempty(draw.mdl.elems(e).load.uniformGbl) && isempty(draw.mdl.elems(e).load.linearGbl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kNm',dc,M1);
                    else
                        value = sprintf('%+.*f',dc,M1);
                    end
                    if M1 >= 0
                        text((x1+x2)/2,(y1+y2)/2,zd1,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                    else
                        text((x1+x2)/2,(y1+y2)/2,zd1,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                    end
                else
                    if M1 >= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value1 = sprintf('%+.*f kNm',dc,M1);
                        else
                            value1 = sprintf('%+.*f',dc,M1);
                        end
                        text(x1,y1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                    else
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value1 = sprintf('%+.*f kNm',dc,M1);
                        else
                            value1 = sprintf('%+.*f',dc,M1);
                        end
                        text(x1,y1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                    end
                    if M2 >= 0
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value2 = sprintf('%+.*f kNm',dc,M2);
                        else
                            value2 = sprintf('%+.*f',dc,M2);
                        end
                        text(x2,y2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',M2);
                    else
                        if strcmp(get(mdata.unitsButton,'Checked'),'on')
                            value2 = sprintf('%+.*f kNm',dc,M2);
                        else
                            value2 = sprintf('%+.*f',dc,M2);
                        end
                        text(x2,y2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',M2);
                    end
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    M = draw.mdl.elems(e).intStresses(2,:);
                    
                    % Avoid plotting numeric garbage
                    if ~all(abs(M) < 10e-10)
                        % Get element 50 point division coords and internal stress values
                        coords = draw.mdl.elems(e).intCoords;
                        Mmax = draw.mdl.elems(e).maxBendMoment_XZ;
                        
                        % Get element basis transformation matrix
                        rot = draw.mdl.elems(e).T;
                        
                        % Compute diagram coordinates
                        diagramCoords = rot' * [zeros(1,size(M,2)); zeros(1,size(M,2)); -scale*M] + coords;
                        
                        % Plot diagram
                        X = diagramCoords(1,:);
                        Y = diagramCoords(2,:);
                        Z = diagramCoords(3,:);
                        line(X,Y,Z,'Color',clr,'tag','drawBendMomentXZDiagram');
                        
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
                            maxPointCoords = rot' * [zeros(1,nm);zeros(1,nm);-scale*Mmax_val] + MmaxGblCoords;
                            xp = maxPointCoords(1,:);
                            yp = maxPointCoords(2,:);
                            zp = maxPointCoords(3,:);
                            scatter3(xp,yp,zp,50,clr,'.','tag','drawBendMomentXZDiagram')
                            
                            % Plot maximum value (avoid plotting text too close to element end)
                            for np = 1:nm % 2 max
                                if abs(Mmax(np,2) - L) >= L/25 && Mmax(np,2) >= L/25
                                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                                        value = sprintf('%+.*f kNm',dc,Mmax_val(np));
                                    else
                                        value = sprintf('%+.*f',dc,Mmax_val(np));
                                    end
                                    if zp(np) >= 0
                                        text(xp(np),yp(np),zp(np),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax_val(np));
                                    else
                                        text(xp(np),yp(np),zp(np),value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax_val(np));
                                    end
                                end
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1,xd2];
                        Y = [yd1,yd2];
                        Z = [zd1,zd2];
                        line(X,Y,Z,'Color',clr,'tag','drawBendMomentXZDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [x1,x2];
                    Y = [y1,y2];
                    Z = [zd1,zd2];
                    line(X,Y,Z,'Color',clr,'tag','drawBendMomentXZDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XZ plane on a
        % given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function bendingMomentEnvelop_XZ(draw,scale)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Mmax = draw.mdl.elems(e).intForcesEnvelop(1,:,2);
                Mmin = draw.mdl.elems(e).intForcesEnvelop(2,:,2);
                
                if all(abs(Mmax) < 10e-10) && all(abs(Mmin) < 10e-10)
                    x1 = draw.mdl.elems(e).nodes(1).coord(1);
                    y1 = draw.mdl.elems(e).nodes(1).coord(2);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%+.*f kNm',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[0,0],'Color',clr,'tag','drawBendMomentXZDiagram');
                    text((x1+x2)/2,(y1+y2)/2,0,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',0);
                    text((x1+x2)/2,(y1+y2)/2,0,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',0);
                else
                    % Get coords of 27 points along element
                    coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                    
                    % Compute diagram coordinates
                    maxCoords = [zeros(1,size(Mmax,2)); zeros(1,size(Mmax,2)); -scale*Mmax] + coords;
                    minCoords = [zeros(1,size(Mmin,2)); zeros(1,size(Mmin,2)); -scale*Mmin] + coords;
                    
                    % Plot envelop
                    Xmax = maxCoords(1,:);
                    Ymax = maxCoords(2,:);
                    Zmax = maxCoords(3,:);
                    Xmin = minCoords(1,:);
                    Ymin = minCoords(2,:);
                    Zmin = minCoords(3,:);
                    line(Xmax,Ymax,Zmax,'Color',clr,'tag','drawBendMomentXZDiagram');
                    line(Xmin,Ymin,Zmin,'Color',clr,'tag','drawBendMomentXZDiagram');
                    
                    % Plot lines to mark 6 internal points along envelop diagram
                    for i = plotValId
                        xx = [ minCoords(1,i) maxCoords(1,i) ];
                        yy = [ minCoords(2,i) maxCoords(2,i) ];
                        zz = [ minCoords(3,i) maxCoords(3,i) ];
                        line(xx, yy, zz, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                    end
                    
                    % Plot text containing values
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
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
                    if Mmax(1) > 0
                        text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Mmax1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(1));
                    else
                        text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Mmax1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(1));
                    end
                    if Mmin(1) > 0
                        text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Mmin1,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(1));
                    else
                        text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Mmin1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(1));
                    end
                    if Mmax(end) > 0
                        text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Mmax2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(end));
                    else
                        text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Mmax2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(end));
                    end
                    if Mmin(end) > 0
                        text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Mmin2,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(end));
                    else
                        text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Mmin2,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(end));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        function draw = reactions(draw)
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            al = draw.size/17;      % reaction symbol size (arrow length)
            ah = draw.size/60;      % reaction symbol size (pyramid height)
            ab = draw.size/130;     % reaction symbol size (pyramid base)
            clr = [0,0,1];          % reaction symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (pyramid/spring height)
            if strcmp(drawSupports,'on')
                ph = draw.size/40;
                sh = draw.size/20;
            else
                ph = draw.size/150;
                sh = draw.size/150;
            end
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                
                % Get reactions values
                rx = draw.mdl.F(draw.mdl.ID(1,n));
                ry = draw.mdl.F(draw.mdl.ID(2,n));
                rz = draw.mdl.F(draw.mdl.ID(3,n));
                
                % Constraint symbol sizes
                cs = 0; hz = 0;
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5)== FIXED_DOF
                    cs = draw.size/150;
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF
                    cs = draw.size/125;
                else
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng == tot && tot > 0
                        cs = draw.size/125;
                    end
                end    
                if draw.mdl.nodes(n).ebc(3)== FIXED_DOF
                    hz = ph;
                elseif draw.mdl.nodes(n).ebc(3)== SPRING_DOF
                    hz = sh;
                end
                
                % Check if rotation is fixed and draw reaction indication
                if ((draw.mdl.nodes(n).ebc(4) == FIXED_DOF)  && (draw.mdl.nodes(n).ebc(5)== FIXED_DOF)) ||...
                   ((draw.mdl.nodes(n).ebc(4) == SPRING_DOF) && (draw.mdl.nodes(n).ebc(5)== SPRING_DOF))     
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value_x = sprintf('%.*f kNm',dc,abs(rx));
                        value_y = sprintf('%.*f kNm',dc,abs(ry));
                    else
                        value_x = sprintf('%.*f',dc,abs(rx));
                        value_y = sprintf('%.*f',dc,abs(ry));
                    end
                    if rx >= 0
                        draw.moment3D(draw,x-cs,y,0,al,ah,ab,'x+',clr,'drawReactions');
                        text(x-cs-al,y,0,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(rx));
                    else
                        draw.moment3D(draw,x+cs,y,0,al,ah,ab,'x-',clr,'drawReactions');
                        text(x+cs+al,y,0,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(rx));
                    end
                    if ry >= 0
                        draw.moment3D(draw,x,y-cs,0,al,ah,ab,'y+',clr,'drawReactions');
                        text(x,y-cs-al,0,value_y,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(ry));
                    else
                        draw.moment3D(draw,x,y+cs,0,al,ah,ab,'y-',clr,'drawReactions');
                        text(x,y+cs+al,0,value_y,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(ry));
                    end
                end
                
                % Check if translation in the Z axis is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(3) == FIXED_DOF) || (draw.mdl.nodes(n).ebc(3) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on')
                        value = sprintf('%.*f kN',dc,abs(rz));
                    else
                        value = sprintf('%.*f',dc,abs(rz));
                    end
                    if rz >= 0
                        draw.arrow3D(draw,x,y,-cs-hz,al,ah,ab,'z+',clr,'drawReactions');
                    else
                        draw.arrow3D(draw,x,y,-cs-hz-al,al,ah,ab,'z-',clr,'drawReactions');
                    end
                    text(x,y,-cs-hz-al,value,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textForceReactions','UserData',abs(rz));
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
            aux_id = [1, 2:2:48, 49, 50]; % 27 auxiliary indexes for local coords
            
            % Get number of points per element
            npe = size(draw.mdl.elems(1).natVibration,2);
            
            % Initialize plotting arrays
            d = zeros(1,draw.mdl.nel*(npe+1));
            coords = zeros(2,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get 28 cross-sections coordinates
                crds = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = draw.mdl.elems(e).natVibration(:,:,nMode);
                
                % Deformed configuration global coordinates
                dfg = crds(3,:) + scale * dl(2,:);
                
                % Concatenate to dislp mtx
                d(1,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan];
                coords(:,(e-1)*(npe+1)+1:e*(npe+1)) = [crds(1:2,:) nan(2,1)];
            end
            % Plot deformed configuration
            plot3(coords(1,:), coords(2,:), d,'Color', clr, 'tag', 'drawVibrationMode');
            drawnow;
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
            aux_id = [1, 2:2:48, 49, 50]; % 27 auxiliary indexes for local coords
            
            % Check if step is not an integer
            dt = rem(step,1);
            if step >= draw.mdl.n_steps + 1
                step = draw.mdl.n_steps;
                dt = 1;
            end
            
            % Get number of points per element
            npe = size(draw.mdl.elems(1).dynamicIntDispl,2);
            
            % Initialize plotting arrays
            d = zeros(1,draw.mdl.nel*(npe+1));
            coords = zeros(2,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get 28 cross-sections coordinates
                crds = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = (1-dt) * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)) +...
                        dt  * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)+1);
                
                % Deformed configuration global coordinates
                dfg = crds(3,:) + scale * dl(2,:);
                
                % Concatenate to dislp mtx
                d(1,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan];
                coords(:,(e-1)*(npe+1)+1:e*(npe+1)) = [crds(1:2,:) nan(2,1)];
            end
            % Plot deformed configuration
            plot3(coords(1,:), coords(2,:), d,'Color', clr, 'tag', 'drawDynamicDeform');
            drawnow
        end
    end
end