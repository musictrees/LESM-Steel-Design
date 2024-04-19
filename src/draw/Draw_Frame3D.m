%% Draw_Frame3D class
%
%% Description
%
% This is a sub-class of the <draw.html *Draw*> class for the
% implementation of the *3D Frame* draw object.
%
classdef Draw_Frame3D < Draw
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = Draw_Frame3D(mdl)
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
            nm = draw.size/200;       % node mark symbol (cube side)
            r  = draw.size/125;       % hinge symbol radius
            ph = draw.size/35;        % translation constraint symbol (pyramid height)
            pb = draw.size/50;        % translation constraint symbol (pyramid base)
            cs = draw.size/50;        % rotation constraint symbol (cube side)
            sh = draw.size/20;        % spring symbol height
            sr = draw.size/105;       % rotational spring symbol radius
            nclr    = [0,0,0];        % node and hinge color
            sclr    = [0.6,0.6,0.6];  % support color
            sprclr  = [0.6,0,0.4];    % spring color
            rotsclr = [0.7,0,0.5];    % rotational spring color
            dc = getappdata(0,'decPrec');  % decimal precision
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                z = draw.mdl.nodes(n).coord(3);
                
                % Distance between translation constraint support symbol and nodal point
                shift = 0;
                
                if draw.mdl.nodes(n).ebc(4) == FREE_DOF || draw.mdl.nodes(n).ebc(5) == FREE_DOF || draw.mdl.nodes(n).ebc(6) == FREE_DOF
                    [tot,hng] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                    if hng == tot && tot > 0
                        shift = r;
                    else
                        draw.cube(x,y,z,nm,nclr,'drawNodes');
                        hold on;
                    end
                elseif draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF && draw.mdl.nodes(n).ebc(6) == FIXED_DOF
                    if strcmp(drawSupports,'on')
                        shift = cs/2;
                    else
                        draw.cube(x,y,z,nm,nclr,'drawNodes');
                        hold on;
                    end
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF && draw.mdl.nodes(n).ebc(6) == SPRING_DOF
                    if ~strcmp(drawSupports,'on')
                        draw.cube(x,y,z,nm,nclr,'drawNodes');
                        hold on;
                    end
                end
                
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
                        draw.pyramid(x-shift,y,z,ph,pb,'x+',sclr,'drawSupports');
                        hold on;
                    else
                        pt = [x,y,z] - shift * dir_x;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_x;dir_y;dir_z],sclr,'drawSupports');
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
                        draw.SpringX_3D(x-shift,y,z,sh,sprclr,'drawSupports');
                        hold on;
                        text(x-shift-sh/2,y,z+0.1*sh,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kx);
                    else
                        pt = [x,y,z] - shift * dir_x;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_x;dir_y;dir_z],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_x;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kx);
                    end
                end
                
                % Draw fixed support in Y direction
                if draw.mdl.nodes(n).ebc(2) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.pyramid(x,y-shift,z,ph,pb,'y+',sclr,'drawSupports');
                        hold on;
                    else
                        pt = [x,y,z] - shift * dir_y;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_y;dir_z;dir_x],sclr,'drawSupports');
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
                        draw.SpringY_3D(x,y-shift,z,sh,sprclr,'drawSupports');
                        hold on;
                        text(x,y-shift-sh/2,z+0.1*sh,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',ky);
                    else
                        pt = [x,y,z] - shift * dir_y;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_y;dir_z;dir_x],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_y;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',ky);
                    end
                end
                
                % Draw fixed support in Z direction
                if draw.mdl.nodes(n).ebc(3) == FIXED_DOF
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        draw.pyramid(x,y,z-shift,ph,pb,'z+',sclr,'drawSupports');
                        hold on;
                    else
                        pt = [x,y,z] - shift * dir_z;
                        draw.pyramid(pt(1),pt(2),pt(3),ph,pb,[dir_z;dir_x;dir_y],sclr,'drawSupports');
                        hold on;
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
                        draw.SpringZ_3D(x,y,z-shift,sh,sprclr,'drawSupports');
                        hold on;
                        text(x+0.05*sh,y+0.05*sh,z-shift-sh/2,value,'HorizontalAlignment','left','VerticalAlignment','middle','Color',sprclr,'tag','textSprings','UserData',kz);
                    else
                        pt = [x,y,z] - shift * dir_z;
                        draw.displSpring_3D(pt(1),pt(2),pt(3),sh,[dir_z;dir_x;dir_y],sprclr,'drawSupports');
                        hold on;
                        txt = [x,y,z] - 0.8*sh * dir_z;
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',sprclr,'tag','textSprings','UserData',kz);
                    end
                end
                
                % Draw fixed rotational support
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF && draw.mdl.nodes(n).ebc(6) == FIXED_DOF
                    draw.cube(x,y,z,cs,sclr,'drawSupports');
                    hold on;
                
                % Draw spring rotational support
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF && draw.mdl.nodes(n).ebc(6) == SPRING_DOF
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
                    draw.rotSpring_3D(x,y,z,sr,rotsclr,'drawSupports');
                    hold on;
                    text(x+sr,y,z,value,'HorizontalAlignment','left','VerticalAlignment','bottom','Color',sprclr,'tag','textRotSprings','UserData',kr);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws elements with hinged or continuous ends.
        function draw = elements(draw)
            % Get flag for semi-rigid joint visualization option status
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSrj =  get(mdata.viewSemiRigidButton,'Checked');
            
            % Parameters
            nm     = draw.size/200;            % node mark symbol (cube side)
            r      = draw.size/125;            % hinge symbol radius
            sr     = draw.size/100;            % semi-rigid joint symbol radius
            clr    = [0,0,0];                  % element color
            sprclr = [0.6,0,0.4];              % spring color
            dc     = getappdata(0,'decPrec');  % decimal precision
            
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
				z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
				z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Get element orientation angle cosine with X and Y axes
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
				cz = draw.mdl.elems(e).cosine_Z;					
                
                % Get element incidence information on nodes
                [toti,hei] = draw.mdl.nodes(n1).elemsIncidence(draw.mdl);
                [totf,hef] = draw.mdl.nodes(n2).elemsIncidence(draw.mdl);
                
                % Set element end coordinates
                xi = x1;
                yi = y1;
				zi = z1;
                xf = x2;
                yf = y2;
				zf = z2;
                
                % Set element initial coordinates by checking if there is a hinge on nodal point position or on element end
                if (hei == toti) && (draw.mdl.nodes(n1).ebc(4) == 0 || draw.mdl.nodes(n1).ebc(5) == 0 || draw.mdl.nodes(n1).ebc(6) == 0) % Hinge on node
                    draw.sphere(x1,y1,z1,r,'drawElements');
                    hold on;
                    xi = x1 + r * cx;
                    yi = y1 + r * cy;
                    zi = z1 + r * cz;
                elseif draw.mdl.elems(e).hingei == 0 % Hinge on element end
                    draw.sphere(x1+r*cx,y1+r*cy,z1+r*cz,r,'drawElements');
                    xi = x1 + 2 * r * cx;
                    yi = y1 + 2 * r * cy;
                    zi = z1 + 2 * r * cz;
                elseif draw.mdl.elems(e).hingei == 2
                    xi = x1 + (sr + nm/2) * cx;
                    yi = y1 + (sr + nm/2) * cy;
                    zi = z1 + (sr + nm/2) * cz;
                    
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),dir(3,:)] = draw.mdl.elems(e).locAxis;
                        draw.srjoint([xi,yi,zi],sr,dir,sprclr);
                        hold on;
                        
                        krxi = draw.mdl.elems(e).kri(1);
                        kryi = draw.mdl.elems(e).kri(2);
                        krzi = draw.mdl.elems(e).kri(3);
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
                            if krzi >= 1000
                                value_z = sprintf('%.*e kNm/rad',dc,krzi);
                            else
                                value_z = sprintf('%.*f kNm/rad',dc,krzi);
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
                            if krzi >= 1000
                                value_z = sprintf('%.*e',dc,krzi);
                            else
                                value_z = sprintf('%.*f',dc,krzi);
                            end
                        end
                        
                        txt_pos_x = [xi,yi,zi] + 4.5*sr*dir(1,:);
                        txt_pos_y = [xi,yi,zi] + 4.5*sr*dir(2,:);
                        txt_pos_z = [xi,yi,zi] + 4.5*sr*dir(3,:);
                        
                        text(txt_pos_x(1),txt_pos_x(2),txt_pos_x(3),value_x,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[1 0.5 0.1],'tag','textSemiRigid','UserData',krxi);
                        text(txt_pos_y(1),txt_pos_y(2),txt_pos_y(3),value_y,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0.65 0],'tag','textSemiRigid','UserData',kryi);
                        text(txt_pos_z(1),txt_pos_z(2),txt_pos_z(3),value_z,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 1],'tag','textSemiRigid','UserData',krzi);
                    else
                        % Connect element end coordinates
                        X_srj = [xi,x1];
                        Y_srj = [yi,y1];
                        Z_srj = [zi,z1];
                        plot3(X_srj,Y_srj,Z_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Set element final coordinates by checking if there is a hinge on nodal point position or on element end
                if (hef == totf) && (draw.mdl.nodes(n2).ebc(4) == 0 || draw.mdl.nodes(n2).ebc(5) == 0 || draw.mdl.nodes(n2).ebc(6) == 0) % Hinge on node
                    draw.sphere(x2,y2,z2,r,'drawElements');
                    hold on;
                    xf = x2 - r * cx;
                    yf = y2 - r * cy;
                    zf = z2 - r * cz;
                elseif draw.mdl.elems(e).hingef == 0 % Hinge on element end
                    draw.sphere(x2-r*cx,y2-r*cy,z2-r*cz,r,'drawElements');
                    xf = x2 - 2 * r * cx;
                    yf = y2 - 2 * r * cy;
                    zf = z2 - 2 * r * cz;
                elseif draw.mdl.elems(e).hingef == 2
                    xf = x2 - (sr + nm/2) * cx;
                    yf = y2 - (sr + nm/2) * cy;
                    zf = z2 - (sr + nm/2) * cz;
                    
                    if strcmp(drawSrj,'on')
                        [dir(1,:),dir(2,:),dir(3,:)] = draw.mdl.elems(e).locAxis;
                        draw.srjoint([xf,yf,zf],sr,dir,sprclr);
                        hold on;
                        
                        krxf = draw.mdl.elems(e).krf(1);
                        kryf = draw.mdl.elems(e).krf(2);
                        krzf = draw.mdl.elems(e).krf(3);
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
                            if krzf >= 1000
                                value_z = sprintf('%.*e kNm/rad',dc,krzf);
                            else
                                value_z = sprintf('%.*f kNm/rad',dc,krzf);
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
                            if krzf >= 1000
                                value_z = sprintf('%.*e',dc,krzf);
                            else
                                value_z = sprintf('%.*f',dc,krzf);
                            end
                        end
                        
                        txt_pos_x = [xf,yf,zf] + 4.5*sr*dir(1,:);
                        txt_pos_y = [xf,yf,zf] + 4.5*sr*dir(2,:);
                        txt_pos_z = [xf,yf,zf] + 4.5*sr*dir(3,:);
                        
                        text(txt_pos_x(1),txt_pos_x(2),txt_pos_x(3),value_x,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[1 0.5 0.1],'tag','textSemiRigid','UserData',krxf);
                        text(txt_pos_y(1),txt_pos_y(2),txt_pos_y(3),value_y,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0.65 0],'tag','textSemiRigid','UserData',kryf);
                        text(txt_pos_z(1),txt_pos_z(2),txt_pos_z(3),value_z,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 1],'tag','textSemiRigid','UserData',krzf);
                    else
                        % Connect element end coordinates
                        X_srj = [xf,x2];
                        Y_srj = [yf,y2];
                        Z_srj = [zf,z2];
                        plot3(X_srj,Y_srj,Z_srj,'Color',sprclr,'tag','drawSemiRigidTemp','linewidth',3.5);
                        hold on;
                    end
                end
                
                % Connect element end coordinates
                X = [xi,xf];
                Y = [yi,yf];
                Z = [zi,zf];
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
            
            % Compute semi-rigid joint symbol radius
            r = sz;
            
            % local axis
            dx = dir;
            dy = [dir(2,:);dir(3,:);dir(1,:)];
            dz = [dir(3,:);dir(1,:);dir(2,:)];

            % arrow properties
            l = 4.0 * r;
            h = 0.8 * r;
            B = 0.5 * r;
            
            % Draw double arrows to indicate the direction of local axis
            aux = [x,y,z] + l * dir(1,:);
            draw.moment3D(draw,aux(1),aux(2),aux(3),l,h,B,dx,[1 0.5 0.1],'drawSemiRigid',true);
            hold on
            aux = [x,y,z] + l * dir(2,:);
            draw.moment3D(draw,aux(1),aux(2),aux(3),l,h,B,dy,[0 0.7 0],'drawSemiRigid',true);
            hold on
            aux = [x,y,z] + l * dir(3,:);
            draw.moment3D(draw,aux(1),aux(2),aux(3),l,h,B,dz,[0 0 1],'drawSemiRigid',true);
            hold on
            
            % Draw sphere to represent semi-rigd joint
            [a, b, c] = sphere;
            s = surf(a * r + x, b * r + y, c * r + z);
            set(s, 'Edgecolor', clr,'FaceColor', [1,1,1], 'tag', 'drawSemiRigid');

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
        % Draws applied nodal loads and moments.
        function draw = nodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off')
                return
            end
            
            % Parameters
            r   = draw.size/125;           % distance between load and nodal point (hinge radius)
            al  = draw.size/12;            % load symbol size (arrow length)
            ah  = draw.size/40;            % load symbol size (arrowhead height)
            ab  = draw.size/80;            % load symbol size (arrowhead base)
            clr = [1,0,0];                 % load color
            dc  = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.static) == 0
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    z = draw.mdl.nodes(n).coord(3);
                    
                    % Get nodal load components
                    fx = draw.mdl.nodes(n).load.static(1);
                    fy = draw.mdl.nodes(n).load.static(2);
                    fz = draw.mdl.nodes(n).load.static(3);
                    mx = draw.mdl.nodes(n).load.static(4);
                    my = draw.mdl.nodes(n).load.static(5);
                    mz = draw.mdl.nodes(n).load.static(6);
                    
                    % Initialize distance between moment symbol and nodal point
                    mshift_x1 = r;
                    mshift_x2 = r;
                    mshift_y1 = r;
                    mshift_y2 = r;
                    mshift_z1 = r;
                    mshift_z2 = r;
                    
                    % Draw load component in X axis
                    if fx > 0
                        draw.arrow3D(draw,x-r,y,z,al,ah,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-r-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        mshift_x1 = r+1.1*al;
                        
                    elseif fx < 0
                        draw.arrow3D(draw,x+r,y,z,al,ah,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+r+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        mshift_x2 = r+1.1*al;
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
                        mshift_y1 = r+1.1*al;
                        
                    elseif fy < 0
                        draw.arrow3D(draw,x,y+r,z,al,ah,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y+r+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        mshift_y2 = r+1.1*al;
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
                        mshift_z1 = r+1.1*al;
                        
                    elseif fz < 0
                        draw.arrow3D(draw,x,y,z+r,al,ah,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z+r+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                        mshift_z2 = r+1.1*al;
                    end
                    
                    % Draw momment component in X axis
                    if mx > 0
                        draw.moment3D(draw,x-mshift_x1,y,z,al,ah/1.4,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x-mshift_x1-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mx));
                        
                    elseif mx < 0
                        draw.moment3D(draw,x+mshift_x2,y,z,al,ah/1.4,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x+mshift_x2+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mx));
                    end
                    
                    % Draw momment component in Y axis
                    if my > 0
                        draw.moment3D(draw,x,y-mshift_y1,z,al,ah/1.4,ab,'y+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y-mshift_y1-(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(my));
                        
                    elseif my < 0
                        draw.moment3D(draw,x,y+mshift_y2,z,al,ah/1.4,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y+mshift_y2+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(my));
                    end
                    
                    % Draw momment component in Z axis
                    if mz > 0
                        draw.moment3D(draw,x,y,z-mshift_z1,al,ah/1.4,ab,'z+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x,y,z-mshift_z1-(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mz));
                        
                    elseif mz < 0
                        draw.moment3D(draw,x,y,z+mshift_z2,al,ah/1.4,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x,y,z+mshift_z2+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mz));
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws applied dynamic nodal loads and moments.
        function draw = dynamicNodalLoads(draw)
            % Check if nodal loads visualization is on
            mdata = guidata(findobj('Tag','GUI_Main'));
            if strcmp(get(mdata.viewNodalLoadsButton,'Checked'),'off') || draw.mdl.drv.whichResponse ~= 1
                return
            end
            
            % Parameters
            r   = draw.size/125;           % distance between load and nodal point (hinge radius)
            m   = draw.size/60;            % concentrated mass symbol radius
            al  = draw.size/12;            % load symbol size (arrow length)
            ah  = draw.size/40;            % load symbol size (arrowhead height)
            ab  = draw.size/80;            % load symbol size (arrowhead base)
            clr = [0,0.7,0];               % load color
            dc  = getappdata(0,'decPrec'); % decimal precision
            
            for n = 1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).load.dynamic) == 0
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
                    mx = draw.mdl.nodes(n).load.dynamic(4);
                    my = draw.mdl.nodes(n).load.dynamic(5);
                    mz = draw.mdl.nodes(n).load.dynamic(6);
                    
                    % Initialize distance between moment symbol and nodal point
                    mshift_x1 = R;
                    mshift_x2 = R;
                    mshift_y1 = R;
                    mshift_y2 = R;
                    mshift_z1 = R;
                    mshift_z2 = R;
                    
                    % Draw load component in X axis
                    if fx > 0
                        draw.arrow3D(draw,x-R,y,z,al,ah,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x-R-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        mshift_x1 = R+1.1*al;
                        
                    elseif fx < 0
                        draw.arrow3D(draw,x+R,y,z,al,ah,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fx));
                        else
                            value = sprintf('%.*f',dc,abs(fx));
                        end
                        text(x+R+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fx));
                        mshift_x2 = R+1.1*al;
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
                        mshift_y1 = R+1.1*al;
                        
                    elseif fy < 0
                        draw.arrow3D(draw,x,y+R,z,al,ah,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fy));
                        else
                            value = sprintf('%.*f',dc,abs(fy));
                        end
                        text(x,y+R+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fy));
                        mshift_y2 = R+1.1*al;
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
                        mshift_z1 = R+1.1*al;
                        
                    elseif fz < 0
                        draw.arrow3D(draw,x,y,z+R,al,ah,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kN',dc,abs(fz));
                        else
                            value = sprintf('%.*f',dc,abs(fz));
                        end
                        text(x,y,z+R+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(fz));
                        mshift_z2 = R+1.1*al;
                    end
                    
                    % Draw momment component in X axis
                    if mx > 0
                        draw.moment3D(draw,x-mshift_x1,y,z,al,ah/1.4,ab,'x+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x-mshift_x1-(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mx));
                        
                    elseif mx < 0
                        draw.moment3D(draw,x+mshift_x2,y,z,al,ah/1.4,ab,'x-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mx));
                        else
                            value = sprintf('%.*f',dc,abs(mx));
                        end
                        text(x+mshift_x2+(ah+al)/1.5,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mx));
                    end
                    
                    % Draw momment component in Y axis
                    if my > 0
                        draw.moment3D(draw,x,y-mshift_y1,z,al,ah/1.4,ab,'y+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y-mshift_y1-(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(my));
                        
                    elseif my < 0
                        draw.moment3D(draw,x,y+mshift_y2,z,al,ah/1.4,ab,'y-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(my));
                        else
                            value = sprintf('%.*f',dc,abs(my));
                        end
                        text(x,y+mshift_y2+(ah+al)/1.5,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(my));
                    end
                    
                    % Draw momment component in Z axis
                    if mz > 0
                        draw.moment3D(draw,x,y,z-mshift_z1,al,ah/1.4,ab,'z+',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x,y,z-mshift_z1-(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mz));
                        
                    elseif mz < 0
                        draw.moment3D(draw,x,y,z+mshift_z2,al,ah/1.4,ab,'z-',clr,'drawNodalLoads');
                        if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                            value = sprintf('%.*f kNm',dc,abs(mz));
                        else
                            value = sprintf('%.*f',dc,abs(mz));
                        end
                        text(x,y,z+mshift_z2+(ah+al)/1.5,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textNodalLoads','UserData',abs(mz));
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
            al    = draw.size/12;   % presc. displ. symbol size (arrow length)
            ah    = draw.size/40;   % presc. displ. symbol size (pyramid height)
            ab    = draw.size/80;   % presc. displ. symbol size (pyramid base)
            clr   = [1,0,1];        % presc. displ. symbol color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (pyramid/spring height)
            if strcmp(drawSupports,'on')
                ph = draw.size/35;
                sh = draw.size/20;
            else
                ph = draw.size/125;
                sh = draw.size/125;
            end
            
            for n =1:draw.mdl.nnp
                if isempty(draw.mdl.nodes(n).prescDispl) == 0
                    x = draw.mdl.nodes(n).coord(1);
                    y = draw.mdl.nodes(n).coord(2);
                    z = draw.mdl.nodes(n).coord(3);
                    
                    % Get prescribed displacement component values and convert it to milimeter and radian
                    dx = 1000 * draw.mdl.nodes(n).prescDispl(1);
                    dy = 1000 * draw.mdl.nodes(n).prescDispl(2);
                    dz = 1000 * draw.mdl.nodes(n).prescDispl(3);
                    rx = draw.mdl.nodes(n).prescDispl(4);
                    ry = draw.mdl.nodes(n).prescDispl(5);
                    rz = draw.mdl.nodes(n).prescDispl(6);
                    
                    % Get direction
                    if draw.mdl.nodes(n).isInclinedSupp
                        [dir_x,dir_y,dir_z] = draw.mdl.nodes(n).getInclinedSuppLocAxis;
                    end
                    
                    % Set rotation constraint symbol (cube side)
                    if draw.mdl.nodes(n).ebc(4) == 1 && draw.mdl.nodes(n).ebc(5) == 1 && draw.mdl.nodes(n).ebc(6) == 1
                        cs = draw.size/100;
                    else
                        cs = 0;
                    end
                    
                    % Initialize distance between moment symbol and nodal point
                    mshift_x = 0;
                    mshift_y = 0;
                    mshift_z = 0;
                    
                    % Check if translation is really fixed in X axis direction and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(1) == 1
                        if dx ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f mm',dc,abs(dx));
                            else
                                value = sprintf('%.*f',dc,abs(dx));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if dx > 0
                                    draw.arrow3D(draw,x-cs-ph,y,z,al,ah,ab,'x+',clr,'drawPrescDispl');
                                    hold on
                                    text(x-cs-ph-(ah+al)/2,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                                else
                                    draw.arrow3D(draw,x-cs-ph-al,y,z,al,ah,ab,'x-',clr,'drawPrescDispl');
                                    hold on
                                    text(x-cs-ph-(ah+al)/3,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                                end
                            else
                                if dx > 0
                                    pt  = [x,y,z] - (cs+ph) * dir_x;
                                    txt = [x,y,z] - (cs+ph+0.7*al) * dir_x;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+ph+al) * dir_x;
                                    txt = [x,y,z] - (cs+ph+0.5*al) * dir_x;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dx));
                            end
                            mshift_x = mshift_x + 1.1*al;
                        end
                        mshift_x = mshift_x + ph;
                        
                    elseif draw.mdl.nodes(n).ebc(1) == 2
                        if draw.mdl.nodes(n).springStiff(1) ~= 0
                            mshift_x = mshift_x + sh;
                        end
                    end
                    
                    % Check if translation is really fixed in Y axis direction and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(2) == 1
                        if dy ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f mm',dc,abs(dy));
                            else
                                value = sprintf('%.*f',dc,abs(dy));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if dy > 0
                                    draw.arrow3D(draw,x,y-cs-ph,z,al,ah,ab,'y+',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y-cs-ph-(ah+al)/2,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                                else
                                    draw.arrow3D(draw,x,y-cs-ph-al,z,al,ah,ab,'y-',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y-cs-ph-(ah+al)/3,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                                end
                            else
                                if dy > 0
                                    pt  = [x,y,z] - (cs+ph) * dir_y;
                                    txt = [x,y,z] - (cs+ph+0.7*al) * dir_y;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+ph+al) * dir_y;
                                    txt = [x,y,z] - (cs+ph+0.5*al) * dir_y;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dy));
                            end
                            mshift_y = mshift_y + 1.1*al;
                        end
                        mshift_y = mshift_y + ph;
                        
                    elseif draw.mdl.nodes(n).ebc(2) == 2
                        if draw.mdl.nodes(n).springStiff(2) ~= 0
                            mshift_y = mshift_y + sh;
                        end
                    end
                    
                    % Check if translation is really fixed in Z axis direction and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(3) == 1
                        if dz ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f mm',dc,abs(dz));
                            else
                                value = sprintf('%.*f',dc,abs(dz));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if dz > 0
                                    draw.arrow3D(draw,x,y,z-cs-ph,al,ah,ab,'z+',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y,z-cs-ph-(ah+al)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
                                else
                                    draw.arrow3D(draw,x,y,z-cs-ph-al,al,ah,ab,'z-',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y,z-cs-ph-(ah+al)/3,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
                                end
                            else
                                if dz > 0
                                    pt  = [x,y,z] - (cs+ph) * dir_z;
                                    txt = [x,y,z] - (cs+ph+0.7*al) * dir_z;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+ph+al) * dir_z;
                                    txt = [x,y,z] - (cs+ph+0.5*al) * dir_z;
                                    draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(dz));
                            end
                            mshift_z = mshift_z + 1.1*al;
                        end
                        mshift_z = mshift_z + ph;
                        
                    elseif draw.mdl.nodes(n).ebc(3) == 2
                        if draw.mdl.nodes(n).springStiff(3) ~= 0
                            mshift_z = mshift_z + sh;
                        end
                    end
                    
                    % Check if rotation is really fixed and draw prescribed displacement indication
                    if draw.mdl.nodes(n).ebc(4) == 1 && draw.mdl.nodes(n).ebc(5) == 1 && draw.mdl.nodes(n).ebc(6) == 1
                        % Draw prescribed displacement indication in X axis direction
                        if rx ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f rad',dc,abs(rx));
                            else
                                value = sprintf('%.*f',dc,abs(rx));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if rx > 0
                                    draw.moment3D(draw,x-cs-mshift_x,y,z,al,ah/1.4,ab,'x+',clr,'drawPrescDispl');
                                    hold on
                                    text(x-cs-mshift_x-(ah+al)/2,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rx));
                                else
                                    draw.moment3D(draw,x-cs-mshift_x-al,y,z,al,ah/1.4,ab,'x-',clr,'drawPrescDispl');
                                    hold on
                                    text(x-cs-mshift_x-(ah+al)/3,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rx));
                                end
                            else
                                if rx > 0
                                    pt  = [x,y,z] - (cs+mshift_x) * dir_x;
                                    txt = [x,y,z] - (cs+mshift_x+0.7*al) * dir_x;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+mshift_x+al) * dir_x;
                                    txt = [x,y,z] - (cs+mshift_x+0.5*al) * dir_x;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_x;dir_y;dir_z],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rx));
                            end
                        end
                        
                        % Draw prescribed displacement indication in Y axis direction
                        if ry ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f rad',dc,abs(ry));
                            else
                                value = sprintf('%.*f',dc,abs(ry));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if ry > 0
                                    draw.moment3D(draw,x,y-cs-mshift_y,z,al,ah/1.4,ab,'y+',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y-cs-mshift_y-(ah+al)/2,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(ry));
                                else
                                    draw.moment3D(draw,x,y-cs-mshift_y-al,z,al,ah/1.4,ab,'y-',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y-cs-mshift_y-(ah+al)/3,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(ry));
                                end
                            else
                                if ry > 0
                                    pt  = [x,y,z] - (cs+mshift_y) * dir_y;
                                    txt = [x,y,z] - (cs+mshift_y+0.7*al) * dir_y;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+mshift_y+al) * dir_y;
                                    txt = [x,y,z] - (cs+mshift_y+0.5*al) * dir_y;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_y;dir_z;dir_x],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(ry));
                            end
                        end
                        
                        % Draw prescribed displacement indication in Z axis direction
                        if rz ~= 0
                            if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                value = sprintf('%.*f rad',dc,abs(rz));
                            else
                                value = sprintf('%.*f',dc,abs(rz));
                            end
                            if ~draw.mdl.nodes(n).isInclinedSupp
                                if rz > 0
                                    draw.moment3D(draw,x,y,z-cs-mshift_z,al,ah/1.4,ab,'z+',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y,z-cs-mshift_z-(ah+al)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rz));
                                else
                                    draw.moment3D(draw,x,y,z-cs-mshift_z-al,al,ah/1.4,ab,'z-',clr,'drawPrescDispl');
                                    hold on
                                    text(x,y,z-cs-mshift_z-(ah+al)/3,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rz));
                                end
                            else
                                if rz > 0
                                    pt  = [x,y,z] - (cs+mshift_z) * dir_z;
                                    txt = [x,y,z] - (cs+mshift_z+0.7*al) * dir_z;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                                else
                                    pt  = [x,y,z] - (cs+mshift_z+al) * dir_z;
                                    txt = [x,y,z] - (cs+mshift_z+0.5*al) * dir_z;
                                    draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_z;dir_x;dir_y],clr,'drawPrescDispl');
                                end
                                hold on
                                text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textPrescDispl','UserData',abs(rz));
                            end
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
                if draw.mdl.nodes(n).initCond(1,1) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(2,1) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(3,1) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(4,1) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(5,1) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(6,1) ~= 0
                    flag = flag + 1;
                end
                if draw.mdl.nodes(n).initCond(1,2) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(2,2) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(3,2) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(4,2) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(5,2) ~= 0 ||...
                   draw.mdl.nodes(n).initCond(6,2) ~= 0
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
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            anl   = get(mdata.popupmenu_AnalysisType,'Value');
            nm    = draw.size/200;   % node mark symbol (cube side)
            cs    = draw.size/50;    % rotation constraint symbol (cube side)
            sr    = draw.size/105;   % rotational spring symbol radius
            r     = draw.size/125;   % hinge symbol radius
            m     = draw.size/60;    % concentrated mass symbol radius
            
            for n = 1:draw.mdl.nnp
                [tot,he] = draw.mdl.nodes(n).elemsIncidence(draw.mdl);
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                z = draw.mdl.nodes(n).coord(3);
                id = sprintf('%d',n);
                
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF &&...
                   draw.mdl.nodes(n).ebc(5) == FIXED_DOF &&...
                   draw.mdl.nodes(n).ebc(6) == FIXED_DOF &&...
                   strcmp(get(mdata.viewSupportsButton,'Checked'),'on')
                    text(x+0.55*cs,y,z+0.55*cs,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF &&...
                       draw.mdl.nodes(n).ebc(5) == SPRING_DOF &&...
                       draw.mdl.nodes(n).ebc(6) == SPRING_DOF &&...
                       strcmp(get(mdata.viewSupportsButton,'Checked'),'on')
                    text(x+0.75*sr,y,z+0.75*sr,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                elseif anl == 2                             &&...
                      ~isempty(draw.mdl.nodes(n).displMass) &&...
                       draw.mdl.nodes(n).displMass > 0      &&...
                       strcmp(get(mdata.viewNodalMassButton,'Checked'),'on')
                    text(x+m,y,z+m,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                elseif he == tot && tot > 0
                    text(x+r,y,z+r,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
                else
                    text(x+0.8*nm,y,z+0.8*nm,id,'HorizontalAlignment','left','VerticalAlignment','baseline','FontWeight','bold','tag','textNodeID');
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
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
            % Estimate maximum displacement of each element (nodal and internal)
            m = zeros(1,draw.mdl.nel);
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                
                % Get maximum nodal displacements
                dx1 = draw.mdl.D(draw.mdl.ID(1,n1));
                dy1 = draw.mdl.D(draw.mdl.ID(2,n1));
                dz1 = draw.mdl.D(draw.mdl.ID(3,n1));
                dx2 = draw.mdl.D(draw.mdl.ID(1,n2));
                dy2 = draw.mdl.D(draw.mdl.ID(2,n2));
                dz2 = draw.mdl.D(draw.mdl.ID(3,n2));
                nodeDispl = [dx1,dy1,dz1,dx2,dy2,dz2];
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
                dsf = draw.size/(4*slider*max_disp);
            end
            setappdata(0,'deform_sf',dsf);
        end
        
        %------------------------------------------------------------------
        % Computes dynamic deformed configuration scale factor.
        function draw = dynamicDeformScaleFactor(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
            % Estimate maximum displacement of each element (nodal and internal)
            m = zeros(1,draw.mdl.nel);
            for e = 1:draw.mdl.nel
                n1 = draw.mdl.elems(e).nodes(1).id;
                n2 = draw.mdl.elems(e).nodes(2).id;
                
                % Get maximum nodal displacements
                ids = [ draw.mdl.ID(1,n1) ;
                        draw.mdl.ID(2,n1) ;
                        draw.mdl.ID(3,n1) ;
                        draw.mdl.ID(1,n2) ;
                        draw.mdl.ID(2,n2) ;
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
                dsf = draw.size/(4*slider*max_disp);
            end
            setappdata(0,'deform_sf',dsf);
        end
        
        %------------------------------------------------------------------
        % Draws structure deformed configuration on a given scale.
        % Input arguments:
        %  scale: deformed configuration scale factor
        function draw = deformConfig(draw,scale)
            %  Parameters
            clr = [1,0,0];  % deformed configuration line color
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get 50 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords;
                
                % Get element axial and transversal internal displacements in local system
                dl = draw.mdl.elems(e).intDispl;
                
                % Rotate displacements vector to global system
                dg = draw.mdl.elems(e).T' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords + scale * dg;
                
                % Plot deformed configuration
                line(dfg(1,:), dfg(2,:), dfg(3,:), 'Color', clr, 'tag', 'drawDeformConfig');
            end
        end
        
        %------------------------------------------------------------------
        % Computes axial force diagram scale factor value.
        function draw = axialScaleFactor(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
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
                asf = draw.size/(2.5*slider*max_val);
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
            clr = [1,0,0];  % diagram line color
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
                cz = draw.mdl.elems(e).cosine_Z;
                
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate diagram coordinates according to element orientation
                [xd1,yd1,zd1] = draw.coordTransf3D(0, scale*N1, 0, x1, y1, z1, e);
                [xd2,yd2,zd2] = draw.coordTransf3D(L, scale*N2, 0, x1, y1, z1, e);
                
                % Draw diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Zi = [z1, zd1];
                line(Xi, Yi, Zi, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                Zf = [z2, zd2];
                line(Xf, Yf, Zf, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                
                % Write force values
                if abs(N1-N2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,N1);
                    else
                        value = sprintf('%+.*f',dc,N1);
                    end
                    [xt,yt,zt] = draw.coordTransf3D(L/2, scale*N1, 0, x1, y1, z1, e);
                    text(xt,yt,zt,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',N1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kN',dc,N1);
                        value2 = sprintf('%+.*f kN',dc,N2);
                    else
                        value1 = sprintf('%+.*f',dc,N1);
                        value2 = sprintf('%+.*f',dc,N2);
                    end
                    text(xd1,yd1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',N1);
                    text(xd2,yd2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',N2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
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
                        Z = diagramCoords(3,:);
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawAxialForceDiagram');

                        % Check if there is a maximum value within the diagram
                        if ~isempty(Nmax)
                            % Compute maximum stress value position, in global coordinates
                            NmaxGblCoords = [x1 + Nmax(2) * cx;
                                             y1 + Nmax(2) * cy;
                                             z1 + Nmax(2) * cz];

                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0; scale*Nmax(1); 0] + NmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            zp = maxPointCoords(3);
                            scatter3(xp, yp, zp, 50, clr, '.', 'tag', 'drawAxialForceDiagram')
                            
                            % Avoid plotting text too close to element end
                            if abs(Nmax(2) - L) >= L/25 && Nmax(2) >= L/25
                                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                    value = sprintf('%+.*f kN',dc,Nmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Nmax(1));
                                end
                                text(xp,yp,zp,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',Nmax(1));
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        Z = [zd1, zd2];
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    Z = [zd1, zd2];
                    line(X, Y, Z, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting axial force envelop diagram on a given scale.
        % Input arguments:
        %  scale: axial force diagram scale factor
        function axialForceEnvelop(draw,scale)
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main')); % handle to main GUI
            clr = [1,0,0];     % diagram line color
            dc = getappdata(0,'decPrec'); % decimal precision
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Null value
                if all(abs(Nmax) < 10e-10) && all(abs(Nmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[z1,z2],'Color',clr,'tag','drawAxialForceDiagram');
                    text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',0);
                end
                
                % Get coords of 27 points along element
                coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                
                % Get element basis transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Compute diagram coordinates
                maxCoords = rot' * [zeros(1,size(Nmax,2)); zeros(1,size(Nmax,2)); scale*Nmax] + coords;
                minCoords = rot' * [zeros(1,size(Nmin,2)); zeros(1,size(Nmin,2)); scale*Nmin] + coords;
                
                % Plot envelop
                Xmax = maxCoords(1,:);
                Ymax = maxCoords(2,:);
                Zmax = maxCoords(3,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                Zmin = minCoords(3,:);
                line(Xmax, Ymax, Zmax, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                line(Xmin, Ymin, Zmin, 'Color', clr, 'tag', 'drawAxialForceDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    zz = [ minCoords(3,i) maxCoords(3,i) ];
                    line(xx, yy, zz, 'Color', clr, 'tag', 'drawAxialForceDiagram');
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
                text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Nmax1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',Nmax(1));
                text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Nmin1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',Nmin(1));
                text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Nmax2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',Nmax(end));
                text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Nmin2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textAxialForceDiagram','UserData',Nmin(end));
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
            % Parameters
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element internal tortion moment value (always uniform) and convert it to string
                T = -draw.mdl.elems(e).torsion_moment(1);
                if T ~= 0
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,T);
                    else
                        value = sprintf('%+.*f',dc,T);
                    end
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kNm',dc,abs(T));
                    else
                        value = sprintf('%.*f',dc,abs(T));
                    end
                    T = abs(T);
                end
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Write tortion moment value in the middle of the element
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textTorsionDiagram','UserData',T);
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting torsion moment envelop diagram.
        function torsionMomentEnvelop(draw,~)
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
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value_max,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textTorsionDiagram','UserData',Tmax(1));
                text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value_min,'HorizontalAlignment','center','VerticalAlignment','top','Color',clr,'tag','textTorsionDiagram','UserData',Tmin(1));
            end
        end
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XY plane.
        function draw = shearScaleFactor_XY(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
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
                ssf = draw.size/(2.5*slider*max_val);
            end
            setappdata(0,'shearXY_sf',ssf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function draw = shearForce_XY(draw,scale)
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle cosine with axes X, Y and Z
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                cz = draw.mdl.elems(e).cosine_Z;
                
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate diagram extrimities coordinates
                [xd1,yd1,zd1] = draw.coordTransf3D(0, scale*Q1, 0, x1, y1, z1, e);
                [xd2,yd2,zd2] = draw.coordTransf3D(0, scale*Q2, 0, x2, y2, z2, e);
                
                % Draw diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Zi = [z1, zd1];
                line(Xi, Yi, Zi, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                Zf = [z2, zd2];
                line(Xf, Yf, Zf, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                
                % Write force values
                if abs(Q1-Q2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,Q1);
                    else
                        value = sprintf('%+.*f',dc,Q1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Q1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kN',dc,Q1);
                        value2 = sprintf('%+.*f kN',dc,Q2);
                    else
                        value1 = sprintf('%+.*f',dc,Q1);
                        value2 = sprintf('%+.*f',dc,Q2);
                    end
                    text(xd1,yd1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Q1);
                    text(xd2,yd2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Q2);
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
                        Z = diagramCoords(3,:);
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXYDiagram');

                        % Check if there is a maximum value within the diagram
                        if ~isempty(Qmax)
                            % Compute maximum stress value position, in global coordinates
                            QmaxGblCoords = [x1 + Qmax(2) * cx;
                                             y1 + Qmax(2) * cy;
                                             z1 + Qmax(2) * cz];

                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0; scale*Qmax(1); 0] + QmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            zp = maxPointCoords(3);
                            scatter3(xp, yp, zp, 50, clr, '.', 'tag', 'drawShearForceXYDiagram')

                            % Avoid plotting text too close to element end
                            if abs(Qmax(2) - L) >= L/25 && Qmax(2) >= L/25
                                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                    value = sprintf('%+.*f kN',dc,Qmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Qmax(1));
                                end
                                text(xp,yp,zp,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Qmax(1));
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        Z = [zd1, zd2];
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    Z = [zd1, zd2];
                    line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XY plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function shearForceEnvelop_XY(draw,scale)
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Null value
                if all(abs(Qmax) < 10e-10) && all(abs(Qmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[z1,z2],'Color',clr,'tag','drawShearForceXYDiagram');
                    text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',0);
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
                Zmax = maxCoords(3,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                Zmin = minCoords(3,:);
                line(Xmax, Ymax, Zmax, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                line(Xmin, Ymin, Zmin, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    zz = [ minCoords(3,i) maxCoords(3,i) ];
                    line(xx, yy, zz, 'Color', clr, 'tag', 'drawShearForceXYDiagram');
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
                text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Qmax1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Qmax(1));
                text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Qmin1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Qmin(1));
                text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Qmax2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Qmax(end));
                text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Qmin2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXYDiagram','UserData',Qmin(end));
            end
        end
        
        %------------------------------------------------------------------
        % Computes shear force diagram scale factor value in XZ plane.
        function draw = shearScaleFactor_XZ(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
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
                ssf = draw.size/(2.5*slider*max_val);
            end
            setappdata(0,'shearXZ_sf',ssf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force diagram in XZ plane on a given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function draw = shearForce_XZ(draw,scale)
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle cosine with axes X, Y and Z
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                cz = draw.mdl.elems(e).cosine_Z;
                
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate diagram extrimities coordinates
                [xd1,yd1,zd1] = draw.coordTransf3D(0, 0, scale*Q1, x1, y1, z1, e);
                [xd2,yd2,zd2] = draw.coordTransf3D(0, 0, scale*Q2, x2, y2, z2, e);
                
                % Draw diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Zi = [z1, zd1];
                line(Xi, Yi, Zi, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                Zf = [z2, zd2];
                line(Xf, Yf, Zf, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                
                % Write force values
                if abs(Q1-Q2) < 10e-10 && isempty(draw.mdl.elems(e).load.linearLcl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,Q1);
                    else
                        value = sprintf('%+.*f',dc,Q1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kN',dc,Q1);
                        value2 = sprintf('%+.*f kN',dc,Q2);
                    else
                        value1 = sprintf('%+.*f',dc,Q1);
                        value2 = sprintf('%+.*f',dc,Q2);
                    end
                    text(xd1,yd1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Q1);
                    text(xd2,yd2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Q2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    Q = draw.mdl.elems(e).intStresses(3,:);
                    
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
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXZDiagram');

                        % Check if there is a maximum value within the diagram
                        if ~isempty(Qmax)
                            % Compute maximum stress value position, in global coordinates
                            QmaxGblCoords = [x1 + Qmax(2) * cx;
                                             y1 + Qmax(2) * cy;
                                             z1 + Qmax(2) * cz];

                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [0; 0; scale*Qmax(1)] + QmaxGblCoords;
                            xp = maxPointCoords(1);
                            yp = maxPointCoords(2);
                            zp = maxPointCoords(3);
                            scatter3(xp, yp, zp, 50, clr, '.', 'tag', 'drawShearForceXZDiagram')

                            % Avoid plotting text too close to element end
                            if abs(Qmax(2) - L) >= L/25 && Qmax(2) >= L/25
                                if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                    value = sprintf('%+.*f kN',dc,Qmax(1));
                                else
                                    value = sprintf('%+.*f',dc,Qmax(1));
                                end
                                text(xp,yp,zp,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        Z = [zd1, zd2];
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    Z = [zd1, zd2];
                    line(X, Y, Z, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting shear force envelop diagram in XZ plane on a
        % given scale.
        % Input arguments:
        %  scale: shear force diagram scale factor
        function shearForceEnvelop_XZ(draw,scale)
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Qmax = draw.mdl.elems(e).intForcesEnvelop(1,:,3);
                Qmin = draw.mdl.elems(e).intForcesEnvelop(2,:,3);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Null value
                if all(abs(Qmax) < 10e-10) && all(abs(Qmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kN',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[z1,z2],'Color',clr,'tag','drawShearForceXZDiagram');
                    text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',0);
                end
                
                % Get coords of 27 points along element
                coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                
                % Get element basis transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Compute diagram coordinates
                maxCoords = rot' * [zeros(1,size(Qmax,2)); zeros(1,size(Qmax,2)); scale*Qmax] + coords;
                minCoords = rot' * [zeros(1,size(Qmin,2)); zeros(1,size(Qmin,2)); scale*Qmin] + coords;
                
                % Plot envelop
                Xmax = maxCoords(1,:);
                Ymax = maxCoords(2,:);
                Zmax = maxCoords(3,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                Zmin = minCoords(3,:);
                line(Xmax, Ymax, Zmax, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                line(Xmin, Ymin, Zmin, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    zz = [ minCoords(3,i) maxCoords(3,i) ];
                    line(xx, yy, zz, 'Color', clr, 'tag', 'drawShearForceXZDiagram');
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
                text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Qmax1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(1));
                text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Qmin1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(1));
                text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Qmax2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmax(end));
                text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Qmin2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textShearForceXZDiagram','UserData',Qmin(end));
            end
        end
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XY plane.
        function draw = bendingMomentScaleFactor_XY(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');

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
                bsf = draw.size/(2.5*slider*max_val);
            end
            setappdata(0,'bendingXY_sf',bsf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XY plane on a given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function draw = bendingMoment_XY(draw,scale)
            include_constants
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle with X axis
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                cz = draw.mdl.elems(e).cosine_Z;
                
                % Get element internal forces value at both ends
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate diagram extrimities coordinates
                [xd1,yd1,zd1] = draw.coordTransf3D(0, -scale*M1, 0, x1, y1, z1, e);
                [xd2,yd2,zd2] = draw.coordTransf3D(0, -scale*M2, 0, x2, y2, z2, e);
                
                % Draw Mz diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Zi = [z1, zd1];
                line(Xi, Yi, Zi, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                Zf = [z2, zd2];
                line(Xf, Yf, Zf, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                
                % Write moment values
                if abs(M1-M2) < 10e-10 && isempty(draw.mdl.elems(e).load.uniformGbl) && isempty(draw.mdl.elems(e).load.linearGbl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,M1);
                    else
                        value = sprintf('%+.*f',dc,M1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',M1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kNm',dc,M1);
                        value2 = sprintf('%+.*f kNm',dc,M2);
                    else
                        value1 = sprintf('%+.*f',dc,M1);
                        value2 = sprintf('%+.*f',dc,M2);
                    end
                    text(xd1,yd1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',M1);
                    text(xd2,yd2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',M2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    M = draw.mdl.elems(e).intStresses(4,:);
                    
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
                        Z = diagramCoords(3,:);
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                        
                        % Check if there is a maximum value within the diagram
                        if ~isempty(Mmax)
                            % Compute maximum stress value position, in global coordinates
                            MmaxGblCoords = [x1 + (Mmax(:,2))' * cx;
                                             y1 + (Mmax(:,2))' * cy;
                                             z1 + (Mmax(:,2))' * cz];
                            
                            % Get maximum stress values
                            Mmax_val = (Mmax(:,1))';
                            
                            % Get number of maximum values within diagram
                            nm = size(Mmax_val,2);
                            
                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [zeros(1,nm); -scale*Mmax_val; zeros(1,nm)] + MmaxGblCoords;
                            xp = maxPointCoords(1,:);
                            yp = maxPointCoords(2,:);
                            zp = maxPointCoords(3,:);
                            scatter3(xp, yp, zp, 50, clr, '.', 'tag', 'drawBendMomentXYDiagram')
                            
                            % Plot maximum value in text
                            for np = 1:nm % 2 rounds max
                                % Avoid plotting text too close to element end
                                if abs(Mmax(np,2) - L) >= L/25 && Mmax(np,2) >= L/25
                                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                        value = sprintf('%+.*f kNm',dc,Mmax_val(np));
                                    else
                                        value = sprintf('%+.*f',dc,Mmax_val(np));
                                    end
                                    text(xp(np),yp(np),zp(np),value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',Mmax_val(np));
                                end
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        Z = [zd1, zd2];
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    Z = [zd1, zd2];
                    line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XY plane on a
        % given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function bendingMomentEnvelop_XY(draw,scale)
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Mmax = draw.mdl.elems(e).intForcesEnvelop(1,:,4);
                Mmin = draw.mdl.elems(e).intForcesEnvelop(2,:,4);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Null value
                if all(abs(Mmax) < 10e-10) && all(abs(Mmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[z1,z2],'Color',clr,'tag','drawBendMomentXYDiagram');
                    text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',0);
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
                Zmax = maxCoords(3,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                Zmin = minCoords(3,:);
                line(Xmax, Ymax, Zmax, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                line(Xmin, Ymin, Zmin, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    zz = [ minCoords(3,i) maxCoords(3,i) ];
                    line(xx, yy, zz, 'Color', clr, 'tag', 'drawBendMomentXYDiagram');
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
                text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Mmax1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',Mmax(1));
                text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Mmin1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',Mmin(1));
                text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Mmax2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',Mmax(end));
                text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Mmin2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXYDiagram','UserData',Mmin(end));
            end
        end
        
        %------------------------------------------------------------------
        % Computes bending moment diagram scale factor value in XZ plane.
        function draw = bendingMomentScaleFactor_XZ(draw)
            mdata  = guidata(findobj('Tag','GUI_Main'));
            slider = get(mdata.slider_Scale,'Max');
            
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
                bsf = draw.size/(2.5*slider*max_val);
            end
            setappdata(0,'bendingXZ_sf',bsf);
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment diagram in XZ plane on a given scale.
        function draw = bendingMoment_XZ(draw,scale)
            include_constants
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get element length
                L = draw.mdl.elems(e).length;
                
                % Get element orientation angle with X axis
                cx = draw.mdl.elems(e).cosine_X;
                cy = draw.mdl.elems(e).cosine_Y;
                cz = draw.mdl.elems(e).cosine_Z;
                
                % Get element internal forces value at both ends
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
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Calculate diagram extrimities coordinates
                [xd1,yd1,zd1] = draw.coordTransf3D(0, 0, -scale*M1, x1, y1, z1, e);
                [xd2,yd2,zd2] = draw.coordTransf3D(0, 0, -scale*M2, x2, y2, z2, e);
                
                % Draw Mz diagram extremities
                Xi = [x1, xd1];
                Yi = [y1, yd1];
                Zi = [z1, zd1];
                line(Xi, Yi, Zi, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                
                Xf = [x2, xd2];
                Yf = [y2, yd2];
                Zf = [z2, zd2];
                line(Xf, Yf, Zf, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                
                % Write force values
                if abs(M1-M2) < 10e-10 && isempty(draw.mdl.elems(e).load.uniformGbl) && isempty(draw.mdl.elems(e).load.linearGbl)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,M1);
                    else
                        value = sprintf('%+.*f',dc,M1);
                    end
                    text((xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                else
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value1 = sprintf('%+.*f kNm',dc,M1);
                        value2 = sprintf('%+.*f kNm',dc,M2);
                    else
                        value1 = sprintf('%+.*f',dc,M1);
                        value2 = sprintf('%+.*f',dc,M2);
                    end
                    text(xd1,yd1,zd1,value1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',M1);
                    text(xd2,yd2,zd2,value2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',M2);
                end
                
                % Connect diagram extremities:
                % Check if element has distributed load.
                % -If so, calculate shear force value along element length
                % -If not, connect both ends with a straight line
                if ~isempty(draw.mdl.elems(e).load.uniformLcl) || ~isempty(draw.mdl.elems(e).load.linearLcl)
                    M = draw.mdl.elems(e).intStresses(5,:);
                    
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
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                        
                        % Check if there is a maximum value within the diagram
                        if ~isempty(Mmax)
                            % Compute maximum stress value position, in global coordinates
                            MmaxGblCoords = [x1 + (Mmax(:,2))' * cx;
                                             y1 + (Mmax(:,2))' * cy;
                                             z1 + (Mmax(:,2))' * cz];

                            % Get maximum stress values
                            Mmax_val = (Mmax(:,1))';

                            % Get number of maximum values within diagram
                            nm = size(Mmax_val,2);
                            
                            % Plot point indicating maximum value
                            maxPointCoords = rot' * [zeros(1,nm); zeros(1,nm); -scale*Mmax_val] + MmaxGblCoords;
                            xp = maxPointCoords(1,:);
                            yp = maxPointCoords(2,:);
                            zp = maxPointCoords(3,:);
                            scatter3(xp, yp, zp, 50, clr, '.', 'tag', 'drawBendMomentXZDiagram')

                            % Plot maximum value in text
                            for np = 1:nm % 2 rounds max
                                % Avoid plotting text too close to element end
                                if abs(Mmax(np,2) - L) >= L/25 && Mmax(np,2) >= L/25
                                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                                        value = sprintf('%+.*f kNm',dc,Mmax_val(np));
                                    else
                                        value = sprintf('%+.*f',dc,Mmax_val(np));
                                    end
                                    text(xp(np),yp(np),zp(np),value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax_val(np));
                                end
                            end
                        end
                    else
                        % Connect both ends with a straight line
                        X = [xd1, xd2];
                        Y = [yd1, yd2];
                        Z = [zd1, zd2];
                        line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                    end
                else
                    % Connect both ends with a straight line
                    X = [xd1, xd2];
                    Y = [yd1, yd2];
                    Z = [zd1, zd2];
                    line(X, Y, Z, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                end
            end
        end
        
        %------------------------------------------------------------------
        % Draws resulting bending moment envelop diagram in XZ plane on a
        % given scale.
        % Input arguments:
        %  scale: bending moment diagram scale factor
        function bendingMomentEnvelop_XZ(draw,scale)
            mdata = guidata(findobj('Tag','GUI_Main'));
            clr = [1,0,0];
            dc = getappdata(0,'decPrec');
            plotValId = [1,6,11,16,21,27];
            
            % Loop over selected elements
            ids = draw.getSelectedElements(draw);
            for e = ids
                % Get internal forces envelop values
                Mmax = draw.mdl.elems(e).intForcesEnvelop(1,:,5);
                Mmin = draw.mdl.elems(e).intForcesEnvelop(2,:,5);
                
                % Get nodal coordinates
                x1 = draw.mdl.elems(e).nodes(1).coord(1);
                y1 = draw.mdl.elems(e).nodes(1).coord(2);
                z1 = draw.mdl.elems(e).nodes(1).coord(3);
                x2 = draw.mdl.elems(e).nodes(2).coord(1);
                y2 = draw.mdl.elems(e).nodes(2).coord(2);
                z2 = draw.mdl.elems(e).nodes(2).coord(3);
                
                % Null value
                if all(abs(Mmax) < 10e-10) && all(abs(Mmin) < 10e-10)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%+.*f kNm',dc,0);
                    else
                        value = sprintf('%+.*f',dc,0);
                    end
                    line([x1,x2],[y1,y2],[z1,z2],'Color',clr,'tag','drawBendMomentXZDiagram');
                    text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,value,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',0);
                end
                
                % Get coords of 27 points along element
                coords = draw.mdl.elems(e).intCoords(:,[1, 2:2:48, 49, 50]);
                
                % Get element basis transformation matrix
                rot = draw.mdl.elems(e).T;
                
                % Compute diagram coordinates
                maxCoords = rot' * [zeros(1,size(Mmax,2)); zeros(1,size(Mmax,2)); -scale*Mmax] + coords;
                minCoords = rot' * [zeros(1,size(Mmin,2)); zeros(1,size(Mmin,2)); -scale*Mmin] + coords;
                
                % Plot envelop
                Xmax = maxCoords(1,:);
                Ymax = maxCoords(2,:);
                Zmax = maxCoords(3,:);
                Xmin = minCoords(1,:);
                Ymin = minCoords(2,:);
                Zmin = minCoords(3,:);
                line(Xmax, Ymax, Zmax, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                line(Xmin, Ymin, Zmin, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
                
                % Plot lines to mark 6 internal points along envelop diagram
                for i = plotValId
                    xx = [ minCoords(1,i) maxCoords(1,i) ];
                    yy = [ minCoords(2,i) maxCoords(2,i) ];
                    zz = [ minCoords(3,i) maxCoords(3,i) ];
                    line(xx, yy, zz, 'Color', clr, 'tag', 'drawBendMomentXZDiagram');
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
                text(maxCoords(1,1),maxCoords(2,1),maxCoords(3,1),value_Mmax1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(1));
                text(minCoords(1,1),minCoords(2,1),minCoords(3,1),value_Mmin1,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(1));
                text(maxCoords(1,end),maxCoords(2,end),maxCoords(3,end),value_Mmax2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmax(end));
                text(minCoords(1,end),minCoords(2,end),minCoords(3,end),value_Mmin2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr,'tag','textBendMomentXZDiagram','UserData',Mmin(end));
            end
        end
        
        %------------------------------------------------------------------
        % Draws reactions indication next to nodal supports.
        function draw = reactions(draw)
            % Parameters
            include_constants;
            mdata = guidata(findobj('Tag','GUI_Main'));
            drawSupports = get(mdata.viewSupportsButton,'Checked');
            al    = draw.size/12;      % reaction symbol size (arrow length)
            ah    = draw.size/40;      % reaction symbol size (pyramid height)
            ab    = draw.size/80;      % reaction symbol size (pyramid base)
            clr   = [0,0,1];           % reaction symbol color
            dc    = getappdata(0,'decPrec'); % decimal precision
            
            % Translation constraint symbol (pyramid/spring height)
            if strcmp(drawSupports,'on')
                ph = draw.size/35;
                sh = draw.size/20;
            else
                ph = draw.size/125;
                sh = draw.size/125;
            end
            
            for n = 1:draw.mdl.nnp
                x = draw.mdl.nodes(n).coord(1);
                y = draw.mdl.nodes(n).coord(2);
                z = draw.mdl.nodes(n).coord(3);
                
                % Get reactions values
                rx = draw.mdl.F(draw.mdl.ID(1,n));
                ry = draw.mdl.F(draw.mdl.ID(2,n));
                rz = draw.mdl.F(draw.mdl.ID(3,n));
                mx = draw.mdl.F(draw.mdl.ID(4,n));
                my = draw.mdl.F(draw.mdl.ID(5,n));
                mz = draw.mdl.F(draw.mdl.ID(6,n));
                
                % Get direction
                if draw.mdl.nodes(n).isInclinedSupp
                    [dir_x,dir_y,dir_z] = draw.mdl.nodes(n).getInclinedSuppLocAxis;
                end
                
                % Set rotation constraint symbol size (cube/spring size)
                cs = 0;
                if draw.mdl.nodes(n).ebc(4) == FIXED_DOF && draw.mdl.nodes(n).ebc(5) == FIXED_DOF && draw.mdl.nodes(n).ebc(6) == FIXED_DOF
                    cs = draw.size/100;
                elseif draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF && draw.mdl.nodes(n).ebc(6) == SPRING_DOF
                    cs = draw.size/70;
                end
                
                % Support (or spring) height in the direction of each axis
                hx = 0;
                hy = 0;
                hz = 0;
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
                
                % Initialize distance between moment symbol and nodal point
                mshift_x = cs;
                mshift_y = cs;
                mshift_z = cs;
                
                % Check if translation is fixed in local axis X and draw reaction indication
                if draw.mdl.nodes(n).ebc(1) == FIXED_DOF || draw.mdl.nodes(n).ebc(1) == SPRING_DOF
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(rx));
                    else
                        value = sprintf('%.*f',dc,abs(rx));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if rx >= 0
                            draw.arrow3D(draw,x-cs-hx,y,z,al,ah,ab,'x+',clr,'drawReactions');
                            hold on
                            text(x-cs-hx-(ah+al)/2,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                        else
                            draw.arrow3D(draw,x-cs-hx-al,y,z,al,ah,ab,'x-',clr,'drawReactions');
                            hold on
                            text(x-cs-hx-(ah+al)/3,y,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                        end
                    else
                        if rx >= 0
                            pt  = [x,y,z] - (cs+hx) * dir_x;
                            txt = [x,y,z] - (cs+hx+0.7*al) * dir_x;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_x;dir_y;dir_z],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (cs+hx+al) * dir_x;
                            txt = [x,y,z] - (cs+hx+0.5*al) * dir_x;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_x;dir_y;dir_z],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rx));
                    end
                    mshift_x = mshift_x + hx + 1.1*al;
                end
                
                % Check if translation is fixed in local axis Y and draw reaction indication
                if draw.mdl.nodes(n).ebc(2) == FIXED_DOF || draw.mdl.nodes(n).ebc(2) == SPRING_DOF
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(ry));
                    else
                        value = sprintf('%.*f',dc,abs(ry));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if ry >= 0
                            draw.arrow3D(draw,x,y-cs-hy,z,al,ah,ab,'y+',clr,'drawReactions');
                            hold on
                            text(x,y-cs-hy-(ah+al)/2,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                        else
                            draw.arrow3D(draw,x,y-cs-hy-al,z,al,ah,ab,'y-',clr,'drawReactions');
                            hold on
                            text(x,y-cs-hy-(ah+al)/3,z,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                        end
                    else
                        if ry >= 0
                            pt  = [x,y,z] - (cs+hy) * dir_y;
                            txt = [x,y,z] - (cs+hy+0.7*al) * dir_y;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_y;dir_z;dir_x],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (cs+hy+al) * dir_y;
                            txt = [x,y,z] - (cs+hy+0.5*al) * dir_y;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_y;dir_z;dir_x],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(ry));
                    end
                    mshift_y = mshift_y + hy + 1.1*al;
                end
                
                % Check if translation is fixed in local axis Z and draw reaction indication
                if draw.mdl.nodes(n).ebc(3) == FIXED_DOF || draw.mdl.nodes(n).ebc(3) == SPRING_DOF
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value = sprintf('%.*f kN',dc,abs(rz));
                    else
                        value = sprintf('%.*f',dc,abs(rz));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if rz >= 0
                            draw.arrow3D(draw,x,y,z-cs-hz,al,ah,ab,'z+',clr,'drawReactions');
                            hold on
                            text(x,y,z-cs-hz-(ah+al)/2,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
                        else
                            draw.arrow3D(draw,x,y,z-cs-hz-al,al,ah,ab,'z-',clr,'drawReactions');
                            hold on
                            text(x,y,z-cs-hz-(ah+al)/3,value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
                        end
                    else
                        if rz >= 0
                            pt  = [x,y,z] - (cs+hz) * dir_z;
                            txt = [x,y,z] - (cs+hz+0.7*al) * dir_z;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,[dir_z;dir_x;dir_y],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (cs+hz+al) * dir_z;
                            txt = [x,y,z] - (cs+hz+0.5*al) * dir_z;
                            draw.arrow3D(draw,pt(1),pt(2),pt(3),al,ah,ab,-[dir_z;dir_x;dir_y],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textForceReactions','UserData',abs(rz));
                    end
                    mshift_z = mshift_z + hz + 1.1*al;
                end
                
                % Check if rotation is fixed and draw reaction indication
                if (draw.mdl.nodes(n).ebc(4) == FIXED_DOF  && draw.mdl.nodes(n).ebc(5) == FIXED_DOF  && draw.mdl.nodes(n).ebc(6) == FIXED_DOF) ||...
                   (draw.mdl.nodes(n).ebc(4) == SPRING_DOF && draw.mdl.nodes(n).ebc(5) == SPRING_DOF && draw.mdl.nodes(n).ebc(6) == SPRING_DOF)
                    if strcmp(get(mdata.unitsButton,'Checked'),'on') == 1
                        value_x = sprintf('%.*f kNm',dc,abs(mx));
                        value_y = sprintf('%.*f kNm',dc,abs(my));
                        value_z = sprintf('%.*f kNm',dc,abs(mz));
                    else
                        value_x = sprintf('%.*f',dc,abs(mx));
                        value_y = sprintf('%.*f',dc,abs(my));
                        value_z = sprintf('%.*f',dc,abs(mz));
                    end
                    if ~draw.mdl.nodes(n).isInclinedSupp
                        if mx >= 0
                            draw.moment3D(draw,x-mshift_x,y,z,al,ah/1.4,ab,'x+',clr,'drawReactions');
                            hold on
                            text(x-mshift_x-(ah+al)/2,y,z,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mx));
                        else
                            draw.moment3D(draw,x-mshift_x-al,y,z,al,ah/1.4,ab,'x-',clr,'drawReactions');
                            hold on
                            text(x-mshift_x-(ah+al)/3,y,z,value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mx));
                        end
                        if my >= 0
                            draw.moment3D(draw,x,y-mshift_y,z,al,ah/1.4,ab,'y+',clr,'drawReactions');
                            hold on
                            text(x,y-mshift_y-(ah+al)/2,z,value_y,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(my));
                        else
                            draw.moment3D(draw,x,y-mshift_y-al,z,al,ah/1.4,ab,'y-',clr,'drawReactions');
                            hold on
                            text(x,y-mshift_y-(ah+al)/3,z,value_y,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(my));
                        end
                        if mz >= 0
                            draw.moment3D(draw,x,y,z-mshift_z,al,ah/1.4,ab,'z+',clr,'drawReactions');
                            hold on
                            text(x,y,z-mshift_z-(ah+al)/2,value_z,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mz));
                        else
                            draw.moment3D(draw,x,y,z-mshift_z-al,al,ah/1.4,ab,'z-',clr,'drawReactions');
                            hold on
                            text(x,y,z-mshift_z-(ah+al)/3,value_z,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mz));
                        end
                    else
                        if mx >= 0
                            pt  = [x,y,z] - (mshift_x) * dir_x;
                            txt = [x,y,z] - (mshift_x+0.7*al) * dir_x;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_x;dir_y;dir_z],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (mshift_x+al) * dir_x;
                            txt = [x,y,z] - (mshift_x+0.5*al) * dir_x;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_x;dir_y;dir_z],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value_x,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mx));
                        if my >= 0
                            pt  = [x,y,z] - (mshift_y) * dir_y;
                            txt = [x,y,z] - (mshift_y+0.7*al) * dir_y;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_y;dir_z;dir_x],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (mshift_y+al) * dir_y;
                            txt = [x,y,z] - (mshift_y+0.5*al) * dir_y;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_y;dir_z;dir_x],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value_y,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(my));
                        if mz >= 0
                            pt  = [x,y,z] - (mshift_z) * dir_z;
                            txt = [x,y,z] - (mshift_z+0.7*al) * dir_z;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,[dir_z;dir_x;dir_y],clr,'drawReactions');
                        else
                            pt  = [x,y,z] - (mshift_z+al) * dir_z;
                            txt = [x,y,z] - (mshift_z+0.5*al) * dir_z;
                            draw.moment3D(draw,pt(1),pt(2),pt(3),al,ah/1.4,ab,-[dir_z;dir_x;dir_y],clr,'drawReactions');
                        end
                        hold on
                        text(txt(1),txt(2),txt(3),value_z,'HorizontalAlignment','center','VerticalAlignment','bottom','Color',clr,'tag','textMomentReactions','UserData',abs(mz));
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
            aux_id = [1, 2:2:48, 49, 50]; % 27 auxiliary indexes for local coords
            
            % Get number of points per element
            npe = size(draw.mdl.elems(1).natVibration,2);
            
            % Initialize plotting matrix
            d = zeros(3,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get 28 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = draw.mdl.elems(e).natVibration(:,:,nMode);
                    
                % Rotate displacements vector to global system
                dg = draw.mdl.elems(e).T' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan(3,1)];
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
            d = zeros(3,draw.mdl.nel*(npe+1));
            
            % Calculate deformed configuration coordinates of 50 cross-
            % sections along element local axis X and connect them
            for e = 1:draw.mdl.nel
                % Get 28 cross-sections coordinates
                coords = draw.mdl.elems(e).intCoords(:,aux_id);
                
                % Get element axial and transversal internal normalized
                % displacements due vibration mode in local system
                dl = (1-dt) * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)) +...
                        dt  * draw.mdl.elems(e).dynamicIntDispl(:,:,floor(step)+1);
                    
                % Rotate displacements vector to global system
                dg = draw.mdl.elems(e).T' * dl;
                
                % Deformed configuration global coordinates
                dfg = coords + scale * dg;
                
                % Concatenate to dislp mtx
                d(:,(e-1)*(npe+1)+1:e*(npe+1)) = [dfg nan(3,1)];
            end
            % Plot deformed configuration
            plot3(d(1,:), d(2,:), d(3,:),'Color', clr, 'tag', 'drawDynamicDeform');
            drawnow
        end
    end
end