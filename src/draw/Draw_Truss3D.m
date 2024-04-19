%% Draw_Truss3D class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *3D Truss* draw object.
%
classdef Draw_Truss3D < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Truss3D(mdl)
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
            z = zeros(draw.mdl.nnp,1);
            for n = 1:draw.mdl.nnp
                x(n) = draw.mdl.nodes(n).coord(1);
                y(n) = draw.mdl.nodes(n).coord(2);
                z(n) = draw.mdl.nodes(n).coord(3);
            end
            dx = max(x) - min(x);
            dy = max(y) - min(y);
            dz = max(z) - min(z);
            sz = max([dx,dy,dz]);
            if isempty(sz) || sz == 0
                draw.size = 5;
            else
                draw.size = sz;
            end
        end
        
        %------------------------------------------------------------------
        % Sets axes limits according to model dimensions.
        function draw = setLimits(draw)
            if draw.mdl.nnp > 0
                x = zeros(draw.mdl.nnp,1);
                y = zeros(draw.mdl.nnp,1);
                z = zeros(draw.mdl.nnp,1);
                
                for n = 1:draw.mdl.nnp
                    x(n) = draw.mdl.nodes(n).coord(1);
                    y(n) = draw.mdl.nodes(n).coord(2);
                    z(n) = draw.mdl.nodes(n).coord(3);
                end
                
                xmin = min(x);
                xmax = max(x);
                ymin = min(y);
                ymax = max(y);
                zmin = min(z);
                zmax = max(z);
                
                dx = xmax - xmin;
                dy = ymax - ymin;
                dz = zmax - zmin;
                
                if (dx == 0) && (dy == 0) && (dz == 0)
                    xlim([x - 2, x + 2])
                    ylim([y - 2, y + 2])
                    zlim([z - 1, z + 1])
                    
                elseif (dx >= dy) && (dx >= dz) ||...
                       (dy >= dx) && (dy >= dz) ||...
                       (dz >= dx) && (dz >= dy)
                    d = draw.size/5.2;
                    xlim([xmin - d, xmax + d])
                    ylim([ymin - d, ymax + d])
                    zlim([zmin - d, zmax + d])
                else
                    xlim([xmin - dx/5, xmax + dx/5])
                    ylim([ymin - dy/5, ymax + dy/5])
                    zlim([zmin - dz/5, zmax + dz/5])
                end
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
            else
                xlim([0,10])
                ylim([0,10])
                zlim([-1,1])
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
            % Get flag for support visualization option status
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            
            % Parameters
            r = draw.size/125;     % hinge symbol radius
            ph = draw.size/35;     % translation constraint symbol (pyramid height)
            pb = draw.size/50;     % translation constraint symbol (pyramid base)
            sh = draw.size/20;     % spring symbol height
            sclr = [0.6,0.6,0.6];  % support color
            sprclr = [0.6,0,0.4];  % spring color
            dc = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
				z = draw.mdl.nodes(n).coord(3);
                
                % Draw hinge as a nodal point
                draw.sphere(x, y, z, r,'drawNodes');
                hold on
                
                % Draw support conditions
                if ~strcmp(drawSupports,'on')
                    continue;
                end
                
                % Get direction
                if draw.mdl.nodes(n).isInclinedSupp
                    [dir_x,dir_y,dir_z] = draw.mdl.nodes(n).getInclinedSuppLocAxis;
                end
                
                % Draw fixed support in X direction
                if draw.mdl.nodes(n).ebc(1) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.pyramid(x-r,y,z,ph,pb,'x+',sclr,'drawSupports');
                        hold on
                    else
                        pt = [x,y,z] - r * dir_x;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_x;dir_y;dir_z],sclr,'drawSupports');
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
                        draw.SpringX_3D(x-r,y,z,sh,sprclr,'drawSupports');
                        hold on;
                        text(x-r-sh/2,y,z+0.1*sh,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kx);
                    else
                        pt = [x,y,z] - r * dir_x;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_x;dir_y;dir_z],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_x;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kx);
                    end
                end
                
                % Draw fixed support in Y direction
                if draw.mdl.nodes(n).ebc(2) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.pyramid(x,y-r,z,ph,pb,'y+',sclr,'drawSupports');
                        hold on
                    else
                        pt = [x,y,z] - r * dir_y;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_y;dir_z;dir_x],sclr,'drawSupports');
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
                        draw.SpringY_3D(x,y-r,z,sh,sprclr,'drawSupports');
                        hold on;
                        text(x,y-r-sh/2,z+0.1*sh,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',ky);
                    else
                        pt = [x,y,z] - r * dir_y;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_y;dir_z;dir_x],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_y;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',ky);
                    end
                end
                
                % Draw fixed support in Z direction
                if draw.mdl.nodes(n).ebc(3) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.pyramid(x,y,z-r,ph,pb,'z+',sclr,'drawSupports');
                        hold on
                    else
                        pt = [x,y,z] - r * dir_z;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_z;dir_x;dir_y],sclr,'drawSupports');
                        hold on
                    end
                
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
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.SpringZ_3D(x,y,z-r,sh,sprclr,'drawSupports');
                        hold on;
                        text(x+0.05*sh,y+0.05*sh,z-r-sh/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',sprclr,'tag','textSprings','UserData',kz);
                    else
                        pt = [x,y,z] - r * dir_z;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_z;dir_x;dir_y],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_z;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kz);
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
				z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
				z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Get element orientation angle cosine with X, Y and Z axes
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
				cz = draw.mdl.elems(e).cosine_Z;
                
                % Set element end coordinates
                xi = x1 + r * cx;
                yi = y1 + r * cy;
				zi = z1 + r * cz;
                xf = x2 - r * cx;
                yf = y2 - r * cy;
				zf = z2 - r * cz;
                
                % Connect element end coordinates
                X = [xi, xf];
                Y = [yi, yf];
				Z = [zi, zf];
                plot3(X, Y, Z, 'Color', clr,'tag','drawElements');
            end
        end
        
        %------------------------------------------------------------------
        function draw = srjoint(draw,~,~,~,~)
        end
        
        %------------------------------------------------------------------
        % Computes element loads scale factor.
        function scl = elemLoadsScaleFactor(draw)
            max_load_y = zeros(1,draw.mdl.nel);
            max_load_z = zeros(1,draw.mdl.nel);
            
            for e = 1:draw.mdl.nel
                % Initialize load values on element ends
                qyi = 0;
                qyf = 0;
				qzi = 0;
				qzf = 0;
                
                % Add uniform load contribtuion
                if isempty(draw.mdl.elems(e).load.uniformLcl) == 0
                    qy = draw.mdl.elems(e).load.uniformLcl(2);
					qz = draw.mdl.elems(e).load.uniformLcl(3);
                    qyi = qyi + qy;
                    qyf = qyf + qy;
					qzi = qzi + qz;
					qzf = qzf + qz;
                end
                
                % Add linear load contribtuion
                if isempty(draw.mdl.elems(e).load.linearLcl) == 0
                    qy1 = draw.mdl.elems(e).load.linearLcl(2);
                    qy2 = draw.mdl.elems(e).load.linearLcl(5);
					qz1 = draw.mdl.elems(e).load.linearLcl(3);
                    qz2 = draw.mdl.elems(e).load.linearLcl(6);
                    qyi = qyi + qy1;
                    qyf = qyf + qy2;
                    qzi = qzi + qz1;
                    qzf = qzf + qz2;
                end
                
                % Get maximum load value on current element
                max_load_y(e) = max(abs(qyi),abs(qyf));
                max_load_z(e) = max(abs(qzi),abs(qzf));
            end
            
            % Get maximum load value on model
            max_val_y = max(max_load_y);
            max_val_z = max(max_load_z);
            max_val = max(max_val_y, max_val_z);
            
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
            include_constants;
            
            % Parameters
            r = draw.size/125;   % hinge symbol radius
            ah = draw.size/200;  % load symbol size (arrowhead height)
            ab = draw.size/300;  % load symbol size (arrowhead base)
            clr = [1,0,0];       % load symbol color
            scl = getappdata(0,'load_sf');                 % Load symbol scale
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for e = 1:draw.mdl.nel
                if ((~isempty(draw.mdl.elems(e).load.uniformGbl))  &&...
                   (~all(draw.mdl.elems(e).load.uniformGbl == 0))) ||...
                   ((~isempty(draw.mdl.elems(e).load.linearGbl))   &&...
                   (~all(draw.mdl.elems(e).load.linearGbl == 0)))
               
                    % Get element length
                    L = draw.mdl.elems(e).length;
                    
                    % Get node IDs
                    n1 = draw.mdl.elems(e).nodes(1).id;
                    
                    % Get nodal coordinates
                    x1 = draw.mdl.nodes(n1).coord(1);
                    y1 = draw.mdl.nodes(n1).coord(2);
					z1 = draw.mdl.nodes(n1).coord(3);
                    
                    % Initialize load values on element ends
                    qxi = 0;
                    qyi = 0;
					qzi = 0;
                    qxf = 0;
                    qyf = 0;
					qzf = 0;
                    
                    % Add uniform load contribtuion
                    if isempty(draw.mdl.elems(e).load.uniformLcl) == 0
                        qxi = qxi + draw.mdl.elems(e).load.uniformLcl(1);
                        qyi = qyi + draw.mdl.elems(e).load.uniformLcl(2);
						qzi = qzi + draw.mdl.elems(e).load.uniformLcl(3);
                        qxf = qxi;
                        qyf = qyi;
						qzf = qzi;
                    end
                    
                    % Add linear load contribtuion
                    if isempty(draw.mdl.elems(e).load.linearLcl) == 0
                        qxi = qxi + draw.mdl.elems(e).load.linearLcl(1);
                        qyi = qyi + draw.mdl.elems(e).load.linearLcl(2);
						qzi = qzi + draw.mdl.elems(e).load.linearLcl(3);
                        qxf = qxf + draw.mdl.elems(e).load.linearLcl(4);
                        qyf = qyf + draw.mdl.elems(e).load.linearLcl(5);
						qzf = qzf + draw.mdl.elems(e).load.linearLcl(6);
                    end
                    
                    % Axial load equation coefficients:
                    % p(x) = Ax + B
                    A = (qxf - qxi)/L;
                    B = qxi;
                    
                    % Transversal load equation coefficients in local axis Y:
                    % q(x) = Cx + D
                    C = (qyf - qyi)/L;
                    D = qyi;
					
                    % Transversal load equation coefficients in local axis Z:
                    % q(x) = Ex + F
                    E = (qzf - qzi)/L;
                    F = qzi;
                    
                    % Module of load parameters
                    Qyi = abs(scl * (C * r + D));
                    Qyf = abs(scl * (C * (L-r) + D));
                    Qzi = abs(scl * (E * r + F));
					Qzf = abs(scl * (E * (L-r) + F));
                    
                    % Get coordinates for text positions
                    if qyi > 0
                        [xb_y,yb_y,zb_y] = draw.coordTransf3D(r, -(Qyi+r), 0, x1, y1, z1, e);
                        [xh_y,yh_y,zh_y] = draw.coordTransf3D(0, -(Qyi+r), 0, x1, y1, z1, e);
                    elseif qyi < 0
                        [xb_y,yb_y,zb_y] = draw.coordTransf3D(r, Qyi+r, 0, x1, y1, z1, e);
                        [xh_y,yh_y,zh_y] = draw.coordTransf3D(0, Qyi+r, 0, x1, y1, z1, e);
                    elseif qyi == 0
                        if qyf > 0
                            [xb_y,yb_y,zb_y] = draw.coordTransf3D(r, -(Qyi+r), 0, x1, y1, z1, e);
                            [xh_y,yh_y,zh_y] = draw.coordTransf3D(0, -(Qyi+r), 0, x1, y1, z1, e);
                        elseif qyf < 0
                            [xb_y,yb_y,zb_y] = draw.coordTransf3D(r, Qyi+r, 0, x1, y1, z1, e);
                            [xh_y,yh_y,zh_y] = draw.coordTransf3D(0, Qyi+r, 0, x1, y1, z1, e);
                        end
                    end
                    
                    if qyf > 0
                        [xc_y,yc_y,zc_y] = draw.coordTransf3D(L-r, -(Qyf+r), 0, x1, y1, z1, e);
                        [xj_y,yj_y,zj_y] = draw.coordTransf3D(L, -(Qyf+r), 0, x1, y1, z1, e);
                    elseif qyf < 0
                        [xc_y,yc_y,zc_y] = draw.coordTransf3D(L-r, Qyf+r, 0, x1, y1, z1, e);
                        [xj_y,yj_y,zj_y] = draw.coordTransf3D(L, Qyf+r, 0, x1, y1, z1, e);
                    elseif qyf == 0
                        if qyi > 0
                            [xc_y,yc_y,zc_y] = draw.coordTransf3D(L-r, -(Qyf+r), 0, x1, y1, z1, e);
                            [xj_y,yj_y,zj_y] = draw.coordTransf3D(L, -(Qyf+r), 0, x1, y1, z1, e);
                        elseif qyi < 0
                            [xc_y,yc_y,zc_y] = draw.coordTransf3D(L-r, Qyf+r, 0, x1, y1, z1, e);
                            [xj_y,yj_y,zj_y] = draw.coordTransf3D(L, Qyf+r, 0, x1, y1, z1, e);
                        end
                    end
                    
                    if qzi > 0
                        [xb_z,yb_z,zb_z] = draw.coordTransf3D(r, 0, -(Qzi+r), x1, y1, z1, e);
                        [xh_z,yh_z,zh_z] = draw.coordTransf3D(0, 0, -(Qzi+r), x1, y1, z1, e);
                    elseif qzi < 0
                        [xb_z,yb_z,zb_z] = draw.coordTransf3D(r, 0, Qzi+r, x1, y1, z1, e);
                        [xh_z,yh_z,zh_z] = draw.coordTransf3D(0, 0, Qzi+r, x1, y1, z1, e);
                    elseif qzi == 0
                        if qzf > 0
                            [xb_z,yb_z,zb_z] = draw.coordTransf3D(r, 0, -(Qzi+r), x1, y1, z1, e);
                            [xh_z,yh_z,zh_z] = draw.coordTransf3D(0, 0, -(Qzi+r), x1, y1, z1, e);
                        elseif qzf < 0
                            [xb_z,yb_z,zb_z] = draw.coordTransf3D(r, 0, Qzi+r, x1, y1, z1, e);
                            [xh_z,yh_z,zh_z] = draw.coordTransf3D(0, 0, Qzi+r, x1, y1, z1, e);
                        end
                    end
                    
                    if qzf > 0
                        [xc_z,yc_z,zc_z] = draw.coordTransf3D(L-r, 0, -(Qzf+r), x1, y1, z1, e);
                        [xj_z,yj_z,zj_z] = draw.coordTransf3D(L, 0, -(Qzf+r), x1, y1, z1, e);
                    elseif qzf < 0
                        [xc_z,yc_z,zc_z] = draw.coordTransf3D(L-r, 0, Qzf+r, x1, y1, z1, e);
                        [xj_z,yj_z,zj_z] = draw.coordTransf3D(L, 0, Qzf+r, x1, y1, z1, e);
                    elseif qzf == 0
                        if qzi > 0
                            [xc_z,yc_z,zc_z] = draw.coordTransf3D(L-r, 0, -(Qzf+r), x1, y1, z1, e);
                            [xj_z,yj_z,zj_z] = draw.coordTransf3D(L, 0, -(Qzf+r), x1, y1, z1, e);
                        elseif qzi < 0
                            [xc_z,yc_z,zc_z] = draw.coordTransf3D(L-r, 0, Qzf+r, x1, y1, z1, e);
                            [xj_z,yj_z,zj_z] = draw.coordTransf3D(L, 0, Qzf+r, x1, y1, z1, e);
                        end
                    end
                    
                    % Draw load symbol on a number cross-sections along element local axis X
                    step = (L-2*r) / round(20 * (L-2*r)/draw.size);
                    for x = r:step:(L-r)
                        % Calculate load values on current cross-section
                        qx = A * x + B;
                        qy = C * x + D;
						qz = E * x + F;
                        Qy = abs(scl * qy);
						Qz = abs(scl * qz);
                        
                        % Calculate current cross-section local coordinates
                        [xs,ys,zs] = draw.coordTransf3D(x, 0, 0, x1, y1, z1, e);
                        
                        % Coordinates of transversal load symbol Qy
                        if qy > 0
                            [xa_y,ya_y,za_y] = draw.coordTransf3D(x, -r, 0, x1, y1, z1, e);
                            [xd_y,yd_y,zd_y] = draw.coordTransf3D(L/2, -(Qy+r), 0, x1, y1, z1, e);
                        elseif qy < 0
                            [xa_y,ya_y,za_y] = draw.coordTransf3D(x, r, 0, x1, y1, z1, e);
                            [xd_y,yd_y,zd_y] = draw.coordTransf3D(L/2, Qy+r, 0, x1, y1, z1, e);
                        end
						
                        % Coordinates of transversal load symbol Qz
                        if qz > 0
                            [xa_z,ya_z,za_z] = draw.coordTransf3D(x, 0, -r, x1, y1, z1, e);
                            [xd_z,yd_z,zd_z] = draw.coordTransf3D(L/2, 0, -(Qz+r), x1, y1, z1, e);
                        elseif qz < 0
                            [xa_z,ya_z,za_z] = draw.coordTransf3D(x, 0, r, x1, y1, z1, e);
                            [xd_z,yd_z,zd_z] = draw.coordTransf3D(L/2, 0, Qz+r, x1, y1, z1, e);
                        end
                        
                        % Draw axial load symbols
                        if ((x ~= r) && ((x+step) < L)) && ((qxi ~= 0) || (qxf ~= 0))
                            if qx > 0
                                draw.spear3D(xs, ys, zs, scl * qx, 2*ah, ab, 'x+', clr, e, draw,'drawElemLoads');
                            elseif qx < 0
                                draw.spear3D(xs, ys, zs, scl * qx, 2*ah, ab, 'x-', clr, e, draw,'drawElemLoads');
                            end
                        end
                                                
                        % Draw transversal load symbol in local axis Y
                        if Qy >= ah
                            if qy > 0
                                draw.spear3D(xa_y, ya_y, za_y, Qy, ah, 0.8*ab, 'y+', clr, e, draw,'drawElemLoads');
                            else
                                draw.spear3D(xa_y, ya_y, za_y, Qy, ah, 0.8*ab, 'y-', clr, e, draw,'drawElemLoads');
                            end
                        elseif (abs(x-r) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [xi,yi,zi] = draw.coordTransf3D(r, 0, 0, x1, y1, z1, e);
                            X = [xi, xb_y];
                            Y = [yi, yb_y];
                            Z = [zi, zb_y];
                            line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                        elseif (abs(x-L+r) <= 10e-10) && ((Qyi ~=0) || (Qyf ~=0))
                            [xf,yf,zf] = draw.coordTransf3D(L-r, 0, 0, x1, y1, z1, e);
                            X = [xf, xc_y];
                            Y = [yf, yc_y];
                            Z = [zf, zc_y];
                            line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                        end
                        
                        % Draw transversal load symbol in local axis Z
                        if Qz >= ah
                            if qz > 0
                                draw.spear3D(xa_z, ya_z, za_z, Qz, ah, 0.8*ab, 'z+', clr, e, draw,'drawElemLoads');
                            else
                                draw.spear3D(xa_z, ya_z, za_z, Qz, ah, 0.8*ab, 'z-', clr, e, draw,'drawElemLoads');
                            end
                        elseif (abs(x-r) <= 10e-10) && ((Qzi ~=0) || (Qzf ~=0))
                            [xi,yi,zi] = draw.coordTransf3D(r, 0, 0, x1, y1, z1, e);
                            X = [xi, xb_z];
                            Y = [yi, yb_z];
                            Z = [zi, zb_z];
                            line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                        elseif (abs(x-L+r) <= 10e-10) && ((Qzi ~=0) || (Qzf ~=0))
                            [xf,yf,zf] = draw.coordTransf3D(L-r, 0, 0, x1, y1, z1, e);
                            X = [xf, xc_z];
                            Y = [yf, yc_z];
                            Z = [zf, zc_z];
                            line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                        end
                    end
                    
                    % Connect diagram extremities in local axis Y
                    if  (qyi < 0) && (qyf > 0)
                        x0 = abs((Qyi*(L-2*r))/(Qyf+Qyi));
                        [xu,yu,zu] = draw.coordTransf3D(x0, r, 0, x1, y1, z1, e);
                        [xd,yd,zd] = draw.coordTransf3D(x0, -r, 0, x1, y1, z1, e);
                        X = [xb_y, xu, xd, xc_y];
                        Y = [yb_y, yu, yd, yc_y];
                        Z = [zb_y, zu, zd, zc_y];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    elseif (qyi > 0) && (qyf < 0)
                        x0 = abs((Qyi*(L-2*r))/(Qyf+Qyi));
                        [xu,yu,zu] = draw.coordTransf3D(x0, r, 0, x1, y1, z1, e);
                        [xd,yd,zd] = draw.coordTransf3D(x0, -r, 0, x1, y1, z1, e);
                        X = [xb_y, xd, xu, xc_y];
                        Y = [yb_y, yd, yu, yc_y];
                        Z = [zb_y, zd, zu, zc_y];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    elseif (qyi ~= 0) || (qyf ~= 0)
                        X = [xb_y, xc_y];
                        Y = [yb_y, yc_y];
						Z = [zb_y, zc_y];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    end
					
                    % Connect diagram extremities in local axis Z
                    if  (qzi < 0) && (qzf > 0)
                        x0 = abs((Qzi*(L-2*r))/(Qzf+Qzi));
                        [xu,yu,zu] = draw.coordTransf3D(x0, 0, r, x1, y1, z1, e);
                        [xd,yd,zd] = draw.coordTransf3D(x0, 0, -r, x1, y1, z1, e);
                        X = [xb_z, xu, xd, xc_z];
                        Y = [yb_z, yu, yd, yc_z];
                        Z = [zb_z, zu, zd, zc_z];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    elseif (qzi > 0) && (qzf < 0)
                        x0 = abs((Qzi*(L-2*r))/(Qzf+Qzi));
                        [xu,yu,zu] = draw.coordTransf3D(x0, 0, r, x1, y1, z1, e);
                        [xd,yd,zd] = draw.coordTransf3D(x0, 0, -r, x1, y1, z1, e);
                        X = [xb_z, xd, xu, xc_z];
                        Y = [yb_z, yd, yu, yc_z];
                        Z = [zb_z, zd, zu, zc_z];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    elseif (qzi ~= 0) || (qzf ~= 0)
                        X = [xb_z, xc_z];
                        Y = [yb_z, yc_z];
						Z = [zb_z, zc_z];
                        line(X, Y, Z, 'Color', clr,'tag','drawElemLoads');
                    end
                    
                    % Write load values:
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value_x1 = sprintf('%.*f kN/m',dc,abs(qxi));
                        value_x2 = sprintf('%.*f kN/m',dc,abs(qxf));
                        value_y1 = sprintf('%.*f kN/m',dc,abs(qyi));
                        value_y2 = sprintf('%.*f kN/m',dc,abs(qyf));
						value_z1 = sprintf('%.*f kN/m',dc,abs(qzi));
                        value_z2 = sprintf('%.*f kN/m',dc,abs(qzf));
                    else
                        value_x1 = sprintf('%.*f',dc,abs(qxi));
                        value_x2 = sprintf('%.*f',dc,abs(qxf));
                        value_y1 = sprintf('%.*f',dc,abs(qyi));
                        value_y2 = sprintf('%.*f',dc,abs(qyf));
						value_z1 = sprintf('%.*f',dc,abs(qzi));
                        value_z2 = sprintf('%.*f',dc,abs(qzf));
                    end

                    % Write axial load values
                    if (qxi ~= 0) || (qxf ~= 0)
                        if qxi == qxf
                            [xm,ym,zm] = draw.coordTransf3D(L/2,0,0,x1,y1,z1,e);
                            text(xm,ym,zm,value_x1,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textElemLoads','UserData',abs(qxi));
                        else
                            [xm,ym,zm] = draw.coordTransf3D(L/50,0,0,x1,y1,z1,e);
                            [xn,yn,zn] = draw.coordTransf3D(L-L/50,0,0,x1,y1,z1,e);
                            text(xm,ym,zm,value_x1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qxi));
                            text(xn,yn,zn,value_x2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qxf));
                        end
                    end
                    
                    % Write transversal load values in local axis Y
                    if (qyi ~= 0) || (qyf ~= 0)
                        if qyi == qyf
                            text(xd_y,yd_y,zd_y,value_y1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qyi));
                        else
                            if qyi ~= 0
                                text(xh_y,yh_y,zh_y,value_y1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qyi));
                            end
                            if qyf ~= 0
                                text(xj_y,yj_y,zj_y,value_y2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qyf));
                            end
                        end
                    end
                    
                    % Write transversal load values in local axis Z
                    if (qzi ~= 0) || (qzf ~= 0)
                        if qzi == qzf
                            text(xd_z,yd_z,zd_z,value_z1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qzi));
                        else
                            if qzi ~= 0
                                text(xh_z,yh_z,zh_z,value_z1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qzi));
                            end
                            if qzf ~= 0
                                text(xj_z,yj_z,zj_z,value_z2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textElemLoads','UserData',abs(qzf));
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
            ah = draw.size/40;   % load symbol size (arrowhead height)
            ab = draw.size/80;   % load symbol size (arrowhead base)
            clr = [1,0,0];       % load symbol color
            dc = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
					z = draw.mdl.nodes(n).coord(3);
                
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.static(1);
                    fy = draw.mdl.nodes(n).load.static(2);
					fz = draw.mdl.nodes(n).load.static(3);
                    
                    % Draw load component in X axis
                    if fx > 0
                        draw.arrow3D(draw,x-r,y,z,al,ah,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-r-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow3D(draw,x+r,y,z,al,ah,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+r+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw load component in Y axis
                    if fy > 0
                        draw.arrow3D(draw,x,y-r,z,al,ah,ab,'y+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y-r-(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow3D(draw,x,y+r,z,al,ah,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y+r+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                    
                    % Draw load component in Z axis
                    if fz > 0
                        draw.arrow3D(draw,x,y,z-r,al,ah,ab,'z+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z-r-(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                        
                    elseif fz < 0
                        draw.arrow3D(draw,x,y,z+r,al,ah,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z+r+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws applied dynamic nodal loads.
        function draw = dynamicNodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off') || draw.mdl.drv.whichResponse ~= 1
                return
            end
            
            % Parameters
            r   = draw.size/125;  % hinge symbol radius
            m   = draw.size/60;   % concentrated mass symbol radius
            al  = draw.size/12;   % load symbol size (arrow length)
            ah  = draw.size/40;   % load symbol size (arrowhead height)
            ab  = draw.size/80;   % load symbol size (arrowhead base)
            clr = [0,0.7,0];      % load symbol color
            dc  = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).load.dynamic)
                    if ~isempty(draw.mdl.nodes(n).displMass) && draw.mdl.nodes(n).displMass > 0
                        R = m;
                    else
                        R = r;
                    end
                    
                    % Get nodal coordinates
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
					z = draw.mdl.nodes(n).coord(3);
                
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.dynamic(1);
                    fy = draw.mdl.nodes(n).load.dynamic(2);
					fz = draw.mdl.nodes(n).load.dynamic(3);
                    
                    % Draw load component in X axis
                    if fx > 0
                        draw.arrow3D(draw,x-R,y,z,al,ah,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-R-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        
                    elseif fx < 0
                        draw.arrow3D(draw,x+R,y,z,al,ah,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+R+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                    end
                    
                    % Draw load component in Y axis
                    if fy > 0
                        draw.arrow3D(draw,x,y-R,z,al,ah,ab,'y+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y-R-(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        
                    elseif fy < 0
                        draw.arrow3D(draw,x,y+R,z,al,ah,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y+R+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                    end
                    
                    % Draw load component in Z axis
                    if fz > 0
                        draw.arrow3D(draw,x,y,z-R,al,ah,ab,'z+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z-R-(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                        
                    elseif fz < 0
                        draw.arrow3D(draw,x,y,z+R,al,ah,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z+R+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
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
                    z = draw.mdl.nodes(n).coord(3);
                    
                    % Get nodal concenctrated mass value
                    mass = draw.mdl.nodes(n).displMass;
                    
                    % Draw mass
                    if mass > 0
                        s = draw.sphere(x, y, z, r, 'drawNodalMass');
                        set(s, 'Edgecolor', clr,'FaceColor', clr);
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kg',dc,abs(mass)*1000);
                        else
                            value = sprintf('%.*f',dc,abs(mass)*1000);
                        end
                        text(x-r,y,z,value,'HorizontalAlignment','right','VerticalAlignment','top','Color',clr,'tag','textNodalMass','UserData',abs(mass)*1000);
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
            r     = draw.size/125;  % hinge symbol radius
            al    = draw.size/12;   % presc. displ. symbol size (arrow length)
            ah    = draw.size/40;   % presc. displ. symbol size (arrowhead height)
            ab    = draw.size/80;   % presc. displ. symbol size (arrowhead base)
            clr   = [1,0,1];        % presc. displ. symbol color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (pyramid height)
            if strcmp(drawSupports,'on')
                ph = draw.size/35;
            else
                ph = 0;
            end
            
            for n = 1:draw.mdl.nnp
                if ~isempty(draw.mdl.nodes(n).prescDispl)
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    z = draw.mdl.nodes(n).coord(3);
                    
                    % Get prescribed displacement component values and convert it to milimeter
                    dx = 1000 * draw.mdl.nodes(n).prescDispl(1);
                    dy = 1000 * draw.mdl.nodes(n).prescDispl(2);
                    dz = 1000 * draw.mdl.nodes(n).prescDispl(3);
                    
                    % Get direction
                    if draw.mdl.nodes(n).isInclinedSupp
                        [dir_x,dir_y,dir_z] = draw.mdl.nodes(n).getInclinedSuppLocAxis;
                    end
                    
                    % Check if translation is really fixed in X axis direction and draw prescribed displacement indication
                    if (draw.mdl.nodes(n).ebc(1) == 1) && (dx ~= 0)
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dx));
                        else
                            value = sprintf('%.*f',dc,abs(dx));
                        end
                        if ~draw.mdl.nodes(n).isInclinedSupp
                            if dx > 0
                                draw.arrow3D(draw,x-r-shift-ph,y,z,al,ah,ab,'x+',clr,'drawPrescDispl');
                                hold on
                                text(x-r-ph-shift-(ah+al)/2,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                            else
                                draw.arrow3D(draw,x-r-shift-ph-al,y,z,al,ah,ab,'x-',clr,'drawPrescDispl');
                                hold on
                                text(x-r-ph-shift-(ah+al)/3,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                            end
                        else
                            if dx > 0
                                pt  = [x,y,z] - (r+ph+shift) * dir_x;
                                txt = [x,y,z] - (r+ph+shift+0.7*al) * dir_x;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                            else
                                pt  = [x,y,z] - (r+ph+shift+al) * dir_x;
                                txt = [x,y,z] - (r+ph+shift+0.5*al) * dir_x;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                            end
                            hold on
                            text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                        end
                    end
                    
                    % Check if translation is really fixed in Y axis direction and draw prescribed displacement indication
                    if (draw.mdl.nodes(n).ebc(2) == 1) && (dy ~= 0)
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dy));
                        else
                            value = sprintf('%.*f',dc,abs(dy));
                        end
                        if ~draw.mdl.nodes(n).isInclinedSupp
                            if dy > 0
                                draw.arrow3D(draw,x,y-r-shift-ph,z,al,ah,ab,'y+',clr,'drawPrescDispl');
                                hold on
                                text(x,y-r-ph-shift-(ah+al)/2,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                            else
                                draw.arrow3D(draw,x,y-r-shift-ph-al,z,al,ah,ab,'y-',clr,'drawPrescDispl');
                                hold on
                                text(x,y-r-ph-shift-(ah+al)/3,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                            end
                        else
                            if dy > 0
                                pt  = [x,y,z] - (r+ph+shift) * dir_y;
                                txt = [x,y,z] - (r+ph+shift+0.7*al) * dir_y;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                            else
                                pt  = [x,y,z] - (r+ph+shift+al) * dir_y;
                                txt = [x,y,z] - (r+ph+shift+0.5*al) * dir_y;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                            end
                            hold on
                            text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                        end
                    end
                    
                    % Check if translation is really fixed in Z axis direction and draw prescribed displacement indication
                    if (draw.mdl.nodes(n).ebc(3) == 1) && (dz ~= 0)
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f mm',dc,abs(dz));
                        else
                            value = sprintf('%.*f',dc,abs(dz));
                        end
                        if ~draw.mdl.nodes(n).isInclinedSupp
                            if dz > 0
                                draw.arrow3D(draw,x,y,z-r-shift-ph,al,ah,ab,'z+',clr,'drawPrescDispl');
                                hold on
                                text(x,y,z-r-ph-shift-(ah+al)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
                            else
                                draw.arrow3D(draw,x,y,z-r-shift-ph-al,al,ah,ab,'z-',clr,'drawPrescDispl');
                                hold on
                                text(x,y,z-r-ph-shift-(ah+al)/3,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
                            end
                        else
                            if dz > 0
                                pt  = [x,y,z] - (r+ph+shift) * dir_z;
                                txt = [x,y,z] - (r+ph+shift+0.7*al) * dir_z;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                            else
                                pt  = [x,y,z] - (r+ph+shift+al) * dir_z;
                                txt = [x,y,z] - (r+ph+shift+0.5*al) * dir_z;
                                draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                            end
                            hold on
                            text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
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
            if strcmp(get(mdata.viewInitialConditionsButton,'Checked'),'off') || draw.mdl.drv.whichResponse ~= 1
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
                if draw.mdl.nodes(n).initCond(1,1) ~= 0 || draw.mdl.nodes(n).initCond(2,1) ~= 0 || draw.mdl.nodes(n).initCond(3,1) ~= 0
                    flag = flag + 1;
                end
                if draw.mdl.nodes(n).initCond(1,2) ~= 0 || draw.mdl.nodes(n).initCond(2,2) ~= 0 || draw.mdl.nodes(n).initCond(3,2) ~= 0
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
                z = draw.mdl.nodes(n).coord(3);
                if ~isempty(draw.mdl.nodes(n).displMass) &&...
                    draw.mdl.nodes(n).displMass > 0      &&...
                    strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x-m,y-m,z+r,str,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'FontWeight','bold','tag','textInitialConditions');
                else
                    text(x-r,y-r,z+r,str,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'FontWeight','bold','tag','textInitialConditions');
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
                    z1 = draw.mdl.elems(e).nodes(1).coord(3);
                    x2 = draw.mdl.elems(e).nodes(2).coord(1);
                    y2 = draw.mdl.elems(e).nodes(2).coord(2);
                    z2 = draw.mdl.elems(e).nodes(2).coord(3);
                    
                    % Get temperature variation values
                    dtx = draw.mdl.elems(e).load.tempVar_X;
                    dty = draw.mdl.elems(e).load.tempVar_Y;
                    dtz = draw.mdl.elems(e).load.tempVar_Z;
                    
                    % New coordinates of the temperature variation
                    [xi_ty1,yi_ty1,zi_ty1] = draw.coordTransf3D(0, d, 0, x1, y1, z1, e);
                    [xf_ty1,yf_ty1,zf_ty1] = draw.coordTransf3D(0, d, 0, x2, y2, z2, e);
                    [xi_ty2,yi_ty2,zi_ty2] = draw.coordTransf3D(0, -d, 0, x1, y1, z1, e);
                    [xf_ty2,yf_ty2,zf_ty2] = draw.coordTransf3D(0, -d, 0, x2, y2, z2, e);
                    
                    [xi_tz1,yi_tz1,zi_tz1] = draw.coordTransf3D(0, 0, d, x1, y1, z1, e);
                    [xf_tz1,yf_tz1,zf_tz1] = draw.coordTransf3D(0, 0, d, x2, y2, z2, e);
                    [xi_tz2,yi_tz2,zi_tz2] = draw.coordTransf3D(0, 0, -d, x1, y1, z1, e);
                    [xf_tz2,yf_tz2,zf_tz2] = draw.coordTransf3D(0, 0, -d, x2, y2, z2, e);
                    
                    % Check if units are enabled
                    unitsAreOn = strcmp(get(mdata.unitsButton,'Checked'),'on');
                    
                    % Draw temperature variation symbols
                    if dtx > 0
                        line([x1,x2],[y1,y2],[z1,z2],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTx = %.*f oC',dc,dtx),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTx = %.*f',dc,dtx),'HorizontalAlignment','center','VerticalAlignment','top','Color',heatClr,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        end
                    elseif dtx < 0
                        line([x1,x2],[y1,y2],[z1,z2],'Color',coldClr,'LineStyle','-.', 'tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTx = %.*f oC',dc,dtx),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTx = %.*f',dc,dtx),'HorizontalAlignment','center','VerticalAlignment','top','Color',coldClr,'tag','textThermalLoads','UserData',{'dTx = ',dtx});
                        end
                    end
                    
                    if dty > 0
                        plot3([xi_ty1,xf_ty1],[yi_ty1,yf_ty1],[zi_ty1,zf_ty1],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        plot3([xi_ty2,xf_ty2],[yi_ty2,yf_ty2],[zi_ty2,zf_ty2],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTy = %.*f oC',dc,dty),'HorizontalAlignment','left','VerticalAlignment','bottom','Color',heatClr,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTy = %.*f',dc,dty),'HorizontalAlignment','left','VerticalAlignment','bottom','Color',heatClr,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        end
                    elseif dty < 0
                        plot3([xi_ty1,xf_ty1],[yi_ty1,yf_ty1],[zi_ty1,zf_ty1],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        plot3([xi_ty2,xf_ty2],[yi_ty2,yf_ty2],[zi_ty2,zf_ty2],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTy = %.*f oC',dc,dty),'HorizontalAlignment','left','VerticalAlignment','bottom','Color',coldClr,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTy = %.*f',dc,dty),'HorizontalAlignment','left','VerticalAlignment','bottom','Color',coldClr,'tag','textThermalLoads','UserData',{'dTy = ',dty});
                        end
                    end
                    
                    if dtz > 0
                        plot3([xi_tz1,xf_tz1],[yi_tz1,yf_tz1],[zi_tz1,zf_tz1],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        plot3([xi_tz2,xf_tz2],[yi_tz2,yf_tz2],[zi_tz2,zf_tz2],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTz = %.*f oC',dc,dtz),'HorizontalAlignment','right','VerticalAlignment','bottom','Color',heatClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','right','VerticalAlignment','bottom','Color',heatClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        end
                    elseif dtz < 0
                        plot3([xi_tz1,xf_tz1],[yi_tz1,yf_tz1],[zi_tz1,zf_tz1],'Color',heatClr,'LineStyle','-.','tag','drawThermalLoads');
                        plot3([xi_tz2,xf_tz2],[yi_tz2,yf_tz2],[zi_tz2,zf_tz2],'Color',coldClr,'LineStyle','-.','tag','drawThermalLoads');
                        if unitsAreOn
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTz = %.*f oC',dc,dtz),'HorizontalAlignment','right','VerticalAlignment','bottom','Color',coldClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
                        else
                            text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,sprintf('dTz = %.*f',dc,dtz),'HorizontalAlignment','right','VerticalAlignment','bottom','Color',coldClr,'tag','textThermalLoads','UserData',{'dTz = ',dtz});
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
                z = draw.mdl.nodes(n).coord(3);
                id = sprintf('%d',n);
                
                if anl == 2                             &&...
                  ~isempty(draw.mdl.nodes(n).displMass) &&...
                   draw.mdl.nodes(n).displMass > 0      &&...
                   strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x+0.8*m,y,z+0.8*m,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                else
                    text(x+r,y,z+r,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Plots ID number of elements.
        function draw = elementID(draw)
            for e = 1:draw.mdl.nel
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                id = sprintf('%d',e);
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,id,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold','tag','textElemID');
            end
        end
        
        %------------------------------------------------------------------
        % Draws element orientation indication from inital to final node.
        function draw = elementOrientation(draw)
            clr = [0,.7,0]; % orientation symbol color
            for e = 1:draw.mdl.nel
                % Calculate spear length
                l = draw.size/25;
                
                % Get nodal coordinates
                xi = draw.mdl.elems(e).nodes(1).coord(1);
                yi = draw.mdl.elems(e).nodes(1).coord(2);
                zi = draw.mdl.elems(e).nodes(1).coord(3);
                xf = draw.mdl.elems(e).nodes(2).coord(1);
                yf = draw.mdl.elems(e).nodes(2).coord(2);
                zf = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate element local axis X orientation vector
                x = [xf-xi, yf-yi, zf-zi];
                x = l * x / norm(x);
                
                % Get orientation vector on the local XZ plane
                z = [draw.mdl.elems(e).vz(1), draw.mdl.elems(e).vz(2), draw.mdl.elems(e).vz(3)];
                
                % Calculate element local axis Y orientation vector
                y = cross(z,x);
                y = l * y / norm(y);
                
                % Calculate element local axis Z orientation vector
                z = cross(x,y);
                z = l * z/norm(z);
                
                % Draw orientation symbol
                xm = (xi + xf)/2;
                ym = (yi + yf)/2;
                zm = (zi + zf)/2;
                
                X = [xm, xm + x(1)];
                Y = [ym, ym + x(2)];
				Z = [zm, zm + x(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.2,'tag','drawElemOrient');
                text(xm+x(1),ym+x(2),zm+x(3),'X','HorizontalAlignment','left','VerticalAlignment','baseline','Color',clr,'FontWeight','bold','tag','drawElemOrient');
                
                X = [xm, xm + y(1)];
                Y = [ym, ym + y(2)];
				Z = [zm, zm + y(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.2,'tag','drawElemOrient');
                text(xm+y(1),ym+y(2),zm+y(3),'Y','HorizontalAlignment','left','VerticalAlignment','baseline','Color',clr,'FontWeight','bold','tag','drawElemOrient');
                
                X = [xm, xm + z(1)];
                Y = [ym, ym + z(2)];
				Z = [zm, zm + z(3)];
                line(X,Y,Z,'Color',clr,'Linewidth',1.2,'tag','drawElemOrient');
                text(xm+z(1),ym+z(2),zm+z(3),'Z','HorizontalAlignment','left','VerticalAlignment','baseline','Color',clr,'FontWeight','bold','tag','drawElemOrient');
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
				dz = draw.mdl.D(draw.mdl.ID(3,n));
                d = [dx,dy,dz];
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
                        draw.mdl.ID(2,n) ;
                        draw.mdl.ID(3,n) ];
                       
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
				z1 = draw.mdl.nodes(n1).coord(3);
                x2 = draw.mdl.nodes(n2).coord(1);
                y2 = draw.mdl.nodes(n2).coord(2);
				z2 = draw.mdl.nodes(n2).coord(3);
                
                % Get nodal displacements
                dx1 = draw.mdl.D(draw.mdl.ID(1,n1));
                dy1 = draw.mdl.D(draw.mdl.ID(2,n1));
				dz1 = draw.mdl.D(draw.mdl.ID(3,n1));
                dx2 = draw.mdl.D(draw.mdl.ID(1,n2));
                dy2 = draw.mdl.D(draw.mdl.ID(2,n2));
				dz2 = draw.mdl.D(draw.mdl.ID(3,n2));
                
                % Calculate displaced nodal coordinates
                xd1 = x1 + scale * dx1;
                yd1 = y1 + scale * dy1;
				zd1 = z1 + scale * dz1;
                xd2 = x2 + scale * dx2;
                yd2 = y2 + scale * dy2;
				zd2 = z2 + scale * dz2;
                
                % Connect displaced nodal coordinates
                X = [xd1, xd2];
                Y = [yd1, yd2];
				Z = [zd1, zd2];
                line(X, Y, Z, 'Color', clr,'tag','drawDeformConfig');
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
                % for 3D truss models) and convert it to string
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Write axial force value
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textAxialForceDiagram','UserData',N);                
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function axialForceEnvelop(draw,~)
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr   = [1,0,0];
            dc    = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
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
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value_max,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textAxialForceDiagram','UserData',Nmax(1));
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value_min,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textAxialForceDiagram','UserData',Nmin(1));
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
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            r     = draw.size/125;  % hinge symbol radius
            al    = draw.size/12;   % reaction symbol size (arrow length)
            ah    = draw.size/40;   % reaction symbol size (arrowhead height)
            ab    = draw.size/80;   % reaction symbol size (arrowhead base)
            clr   = [0,0,1];        % reaction symbol color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (pyramid/spring height)
            if strcmp(drawSupports,'on')
                ph = draw.size/35;
                sh = draw.size/20;
            else
                ph = 0;
                sh = 0;
            end
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
				z = draw.mdl.nodes(n).coord(3);
                
                % Get reactions values
                rx = draw.mdl.F(draw.mdl.ID(1,n));
                ry = draw.mdl.F(draw.mdl.ID(2,n));
				rz = draw.mdl.F(draw.mdl.ID(3,n));
                
                % Get direction
                if draw.mdl.nodes(n).isInclinedSupp
                    [dir_x,dir_y,dir_z] = draw.mdl.nodes(n).getInclinedSuppLocAxis;
                end
                
                % Support (or spring) height in the direction of each axis
                if draw.mdl.nodes(n).ebc(1)== FIXED_DOF
                    hx = ph;
                elseif draw.mdl.nodes(n).ebc(1)== SPRING_DOF
                    hx = sh;
                end
                if draw.mdl.nodes(n).ebc(2)== FIXED_DOF
                    hy = ph;
                elseif draw.mdl.nodes(n).ebc(2)== SPRING_DOF
                    hy = sh;
                end
                if draw.mdl.nodes(n).ebc(3)== FIXED_DOF
                    hz = ph;
                elseif draw.mdl.nodes(n).ebc(3)== SPRING_DOF
                    hz = sh;
                end
                
                % Check if translation is fixed in local axis X and draw reaction indication
                if (draw.mdl.nodes(n).ebc(1) == FIXED_DOF) || (draw.mdl.nodes(n).ebc(1) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(rx));
                    else
                        value = sprintf('%.*f',dc,abs(rx));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if rx >= 0
                            draw.arrow3D(draw,x-r-hx,y,z,al,ah,ab,'x+',clr,'drawReactions');
                            hold on
                            text(x-r-hx-(ah+al)/2,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                        else
                            draw.arrow3D(draw,x-r-hx-al,y,z,al,ah,ab,'x-',clr,'drawReactions');
                            hold on
                            text(x-r-hx-(ah+al)/3,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                        end
                    else
                        if rx >= 0
                            pt  = [x,y,z] - (r+hx) * dir_x;
                            txt = [x,y,z] - (r+hx+0.7*al) * dir_x;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_x;dir_y;dir_z],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (r+hx+al) * dir_x;
                            txt = [x,y,z] - (r+hx+0.5*al) * dir_x;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_x;dir_y;dir_z],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                    end
                end
                
                % Check if translation is fixed in local axis Y and draw reaction indication
                if (draw.mdl.nodes(n).ebc(2) == FIXED_DOF) || (draw.mdl.nodes(n).ebc(2) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(ry));
                    else
                        value = sprintf('%.*f',dc,abs(ry));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if ry >= 0
                            draw.arrow3D(draw,x,y-r-hy,z,al,ah,ab,'y+',clr,'drawReactions');
                            hold on
                            text(x,y-r-hy-(ah+al)/2,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                        else
                            draw.arrow3D(draw,x,y-r-hy-al,z,al,ah,ab,'y-',clr,'drawReactions');
                            hold on
                            text(x,y-r-hy-(ah+al)/3,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                        end
                    else
                        if ry >= 0
                            pt  = [x,y,z] - (r+hy) * dir_y;
                            txt = [x,y,z] - (r+hy+0.7*al) * dir_y;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_y;dir_z;dir_x],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (r+hy+al) * dir_y;
                            txt = [x,y,z] - (r+hy+0.5*al) * dir_y;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_y;dir_z;dir_x],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                    end
                end
                
                % Check if vertical translation in z is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(3) == FIXED_DOF) || (draw.mdl.nodes(n).ebc(3) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(rz));
                    else
                        value = sprintf('%.*f',dc,abs(rz));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if rz >= 0
                            draw.arrow3D(draw,x,y,z-r-hz,al,ah,ab,'z+',clr,'drawReactions');
                            hold on
                            text(x,y,z-r-hz-(ah+al)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
                        else
                            draw.arrow3D(draw,x,y,z-r-hz-al,al,ah,ab,'z-',clr,'drawReactions');
                            hold on
                            text(x,y,z-r-hz-(ah+al)/3,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
                        end
                    else
                        if rz >= 0
                            pt  = [x,y,z] - (r+hz) * dir_z;
                            txt = [x,y,z] - (r+hz+0.7*al) * dir_z;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_z;dir_x;dir_y],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (r+hz+al) * dir_z;
                            txt = [x,y,z] - (r+hz+0.5*al) * dir_z;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_z;dir_x;dir_y],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
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
            d = zeros(3,draw.mdl.nel*3);
            
            % Loop over elems to concat nodal displacements
            for e = 1:draw.mdl.nel
                % Get end point coordinates
                coords = draw.mdl.elems(e).intCoords(:,[1,end]);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = draw.mdl.elems(e).natVibration(:,[1,end],nMode);
                
                % Get rotation transformation matrix
                rot = draw.mdl.elems(e).T;
                    
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*3+1:e*3) = [dfg nan(3,1)];
            end
            % Plot deformed configuration
            plot3(d(1,:), d(2,:), d(3,:),'Color', clr, 'tag', 'drawVibrationMode');
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
            d = zeros(3,draw.mdl.nel*3);
            
            % Loop over elems to concat nodal displacements
            for e = 1:draw.mdl.nel
                % Get end point coordinates
                coords = draw.mdl.elems(e).intCoords(:,[1,end]);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = (1-dt) * draw.mdl.elems(e).dynamicIntDispl(:,[1,end],floor(step)) +...
                        dt  * draw.mdl.elems(e).dynamicIntDispl(:,[1,end],floor(step)+1);
                
                % Get rotation transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Rotate displacements vector to global system
                dg = rot' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*3+1:e*3) = [dfg nan(3,1)];
            end
            % Plot deformed configuration
            plot3(d(1,:), d(2,:), d(3,:),'Color', clr, 'tag', 'drawDynamicDeform');
            drawnow
        end
    end
end