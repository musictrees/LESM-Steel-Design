%% Print_Truss3D Class
%
%% Description
%
% This is a sub-class of the <print.html *Print*> class for the
% implementation of the *3D Truss* print object.
%
classdef Print_Truss3D < Print
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print_Truss3D(model)
            print = print@Print(1);
            
            if (nargin > 0)
                print.model = model;
            end
        end
    end

    %% Public methods
    % Implementation of the abstract methods declared in super-class <print.html *Print*>.
    methods
        %------------------------------------------------------------------
        % Prints analyis results.
        function results(print,lc,currentLc,GUI_Mode)
            if nargin < 4
                GUI_Mode = false;
            end

            if GUI_Mode
                % Create a wait bar to show the progress of the printing process
                h = waitbar(0,'Working on text file...','WindowStyle','modal');
            end

            print.header();
            if GUI_Mode
                waitbar(1/16);
            end

            fprintf(print.txt, '\n\n\n____________ M O D E L  I N F O R M A T I O N ____________\n');
            print.modelLabel();
            [n_srj,n_ns,n_nl,n_pd,n_ul,n_ll,n_tv] = print.modelDescrip(lc,currentLc);
            if GUI_Mode
                waitbar(2/16);
            end

            print.material();
            if GUI_Mode
                waitbar(3/16);
            end

            print.section();
            if GUI_Mode
                waitbar(4/16);
            end

            print.nodalCoords();
            if GUI_Mode
                waitbar(5/16);
            end

            print.nodalSupport();
            if GUI_Mode
                waitbar(6/16);
            end

            print.spring(n_ns);
            if GUI_Mode
                waitbar(7/16);
            end

            print.nodalLoads(n_nl);
            if GUI_Mode
                waitbar(8/16);
            end

            print.nodalPrescDisp(n_pd);
            if GUI_Mode
                waitbar(9/16);
            end

            print.elements();
            if GUI_Mode
                waitbar(10/16);
            end

            print.srjoints(n_srj);
            if GUI_Mode
                waitbar(11/16);
            end

            print.unifElementLoads(n_ul);
            print.linearElementLoads(n_ll);
            print.temperatureVariation(n_tv);
            if GUI_Mode
                waitbar(12/16);
            end

            fprintf(print.txt, '\n\n\n\n_____________ A N A L Y S I S  R E S U L T S _____________\n');
            print.nodalDisplRot();
            if GUI_Mode
                waitbar(13/16);
            end

            print.reactions();
            if GUI_Mode
                waitbar(14/16);
            end

            print.intForces();
            if GUI_Mode
                waitbar(15/16);
            end

            print.elemDispl();
            if GUI_Mode
                waitbar(1);
            end

            if GUI_Mode
                % Close waitbar
                close(h)
            end
        end

        %------------------------------------------------------------------
        % Prints model type.
        function modelLabel(print)
            fprintf(print.txt, '\n\n----------------------------\n');
            fprintf(print.txt, 'MODEL TYPE: ');
            fprintf(print.txt, '3D TRUSS\n');
            fprintf(print.txt, '----------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints model description.
        % Output:
        %  n_srj: number of semi-rigid joints
        %  n_ns: number of nodes with spring supports
        %  n_nl: number of nodes with applied loads
        %  n_pd: number of nodes with prescribed displacement
        %  n_ul: number of elements with uniformly distributed loads
        %  n_ll: number of elements with linearly distributed loads
        %  n_tv: number of elements with temperature variation
        function [n_srj,n_ns,n_nl,n_pd,n_ul,n_ll,n_tv] = modelDescrip(print,lc,currentLc)
            include_constants;
            % Initialize variables
            neq = print.model.nnp * print.model.anm.ndof;
            neqfixed = 0;
            neqspring = 0;
            n_srj = 0;
            n_ns = 0;
            n_nl = 0;
            n_pd = 0;
            n_ul  = 0;
            n_ll = 0;
            n_tv = 0;
            
            % Loop over all nodes
            for n = 1:print.model.nnp
                % Increment number of fixed d.o.f.
                for i = 1:3
                    if print.model.nodes(n).ebc(i) == FIXED_DOF
                        neqfixed = neqfixed + 1;
                    elseif print.model.nodes(n).ebc(i) == SPRING_DOF
                        neqspring = neqspring + 1;
                    end
                end
                
                % Increment number of nodes with spring supports
                if ~isempty(print.model.nodes(n).springStiff)
                    n_ns = n_ns + 1;
                end
                
                % Increment number of nodes with applied loads
                if ~isempty(print.model.nodes(n).load.static) && ~all(print.model.nodes(n).load.static == 0)
                    n_nl = n_nl + 1;
                end
                
                % Increment number of nodes with prescribed displacement
                if ~isempty(print.model.nodes(n).prescDispl) && ~all(print.model.nodes(n).prescDispl == 0)
                    n_pd = n_pd + 1;
                end
            end
            neqfree = neq - neqfixed - neqspring; % number of free d.o.f.
            
            % Loop over all elements
            for e = 1:print.model.nel
                % Increment number of elements with uniformly distributed loads
                if ~isempty(print.model.elems(e).load.uniformGbl) && ~all(print.model.elems(e).load.uniformGbl == 0)
                    n_ul = n_ul + 1;
                end
                
                % Increment number of elements with linearly distributed loads
                if ~isempty(print.model.elems(e).load.linearGbl) && ~all(print.model.elems(e).load.linearGbl == 0)
                    n_ll = n_ll + 1;
                end
                
                % Increment number of elements with temperature variation
                if (print.model.elems(e).load.tempVar_X ~= 0) || ...
                   (print.model.elems(e).load.tempVar_Y ~= 0) || ...
                   (print.model.elems(e).load.tempVar_Z ~= 0)
                    n_tv = n_tv + 1;
                end
            end

            % Check if there are load case combinations
            if isempty(print.model.loadComb) == 1
                ncomb = 0;
            else
                ncomb = size(print.model.loadComb,2);
            end    
            
            % Check if current load cas/comb is case or comb
            if lc > print.model.nlc
                str = 'COMB';
                value = lc - print.model.nlc;
            else
                str = 'CASE';
                value = lc;
            end    

            fprintf(print.txt, '\n\n----------------------------------------\n' );
            fprintf(print.txt, 'M O D E L  D E S C R I P T I O N:\n' );
            fprintf(print.txt, '----------------------------------------\n');
            fprintf(print.txt, 'NUMBER OF NODES.......................:%4d\n', print.model.nnp);
            fprintf(print.txt, 'NUMBER OF ELEMENTS ...................:%4d\n', print.model.nel);
            fprintf(print.txt, 'NUMBER OF DEGREES OF FREEDOM..........:%4d\n', neq);
            fprintf(print.txt, 'NUMBER OF FREE DEGREES OF FREEDOM.....:%4d\n', neqfree);
            fprintf(print.txt, 'NUMBER OF FIXED DEGREES OF FREEDOM....:%4d\n', neqfixed);
            fprintf(print.txt, 'NUMBER OF SPRINGS.....................:%4d\n', neqspring);
            fprintf(print.txt, 'NUMBER OF NODES W/ SPRINGS............:%4d\n', n_ns);
            fprintf(print.txt, 'NUMBER OF MATERIALS...................:%4d\n', print.model.nmat);
            fprintf(print.txt, 'NUMBER OF CROSS-SECTIONS..............:%4d\n', print.model.nsec);
            fprintf(print.txt, 'NUMBER OF LOAD CASES..................:%4d\n', print.model.nlc);
            fprintf(print.txt, 'NUMBER OF LOAD CASE COMBINATIONS......:%4d\n', ncomb);
            fprintf(print.txt, 'CURRENT LOAD CASE / COMBINATION.......:%10s (%4s %2d)\n', char(currentLc), char(str), value);
            fprintf(print.txt, 'NUMBER OF NODES W/ APPLIED LOADS......:%4d\n', n_nl);
            fprintf(print.txt, 'NUMBER OF NODES W/ PRESC. DISPL.......:%4d\n', n_pd);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ UNIFORM LOAD....:%4d\n', n_ul);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ LINEAR LOAD.....:%4d\n', n_ll);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ TEMP. VAR.......:%4d\n', n_tv);
            fprintf(print.txt, '----------------------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints cross-section properties.
        function section(print)
            fprintf(print.txt, '\n\n----------------------------------------------\n');
            fprintf(print.txt, 'C R O S S  S E C T I O N  P R O P E R T I E S\n');
            fprintf(print.txt, '----------------------------------------------\n');

            if print.model.nsec > 0
                fprintf(print.txt, ' SECTION     AREA [cm²]      HEIGHT Y [cm]      HEIGHT Z [cm]\n');

                for s = 1:print.model.nsec
                    fprintf(print.txt, '%4d   %14.2f   %13.2f     %13.2f\n', s,...
                            1e4*print.model.sections(s).area_x,...
                            1e2*print.model.sections(s).height_y,...
							1e2*print.model.sections(s).height_z);
                end
            else
                fprintf(print.txt, ' NO CROSS-SECTION\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal support conditions.
        function nodalSupport(print)
            include_constants;
            
            fprintf(print.txt, '\n\n------------------------------\n');
            fprintf(print.txt, 'N O D A L  R E S T R A I N T S \n');
            fprintf(print.txt, '------------------------------\n');
            fprintf(print.txt, ' NODE   DISPL X    DISPL Y    DISPL Z   |   INCL.SUPP    ROT X [°]    ROT Y [°]    ROT Z [°]    Vx_x     Vx_y     Vx_z    Vy_x     Vy_y     Vy_z\n');

            for n = 1:print.model.nnp
                if(print.model.nodes(n).ebc(1) == FIXED_DOF)
                    node_restr1 = 'FIXED';
                elseif (print.model.nodes(n).ebc(1) == SPRING_DOF)
                    node_restr1 = 'SPRING';
                else
                    node_restr1 = 'FREE';
                end

                if(print.model.nodes(n).ebc(2) == FIXED_DOF)
                    node_restr2 = 'FIXED';
                elseif (print.model.nodes(n).ebc(2) == SPRING_DOF)
                    node_restr2 = 'SPRING';
                else
                    node_restr2 = 'FREE';
                end
                
                if(print.model.nodes(n).ebc(3) == FIXED_DOF)
                    node_restr3 = 'FIXED';
                elseif (print.model.nodes(n).ebc(3) == SPRING_DOF)
                    node_restr3 = 'SPRING';
                else
                    node_restr3 = 'FREE';
                end

                if ~print.model.nodes(n).isInclinedSupp
                    incl = ' NO';
                    rot_x = '   --- ';
                    rot_y = '   --- ';
                    rot_z = '   --- ';
                    dir_x1 = '  ---';
                    dir_x2 = '  ---';
                    dir_x3 = '  ---';
                    dir_y1 = '  ---';
                    dir_y2 = '  ---';
                    dir_y3 = '  ---';
                else
                    % Compute direction vector as angle
                    dir_1 = print.model.nodes(n).inclSuppDir(1:2);
                    aux_norm = norm(dir_1);
                    dir_1 = dir_1 / aux_norm;
                    thetaZ = acos(dir_1(1));
                    if dir_1(2) < 0
                        thetaZ = - thetaZ;
                    end
                    
                    % Convert from rad to degrees
                    thetaZ = thetaZ * 180 / pi;
                    if abs(thetaZ) < 10^-10
                        thetaZ = 0;
                    end
                    
                    % Compute direction vector as angle
                    dir_2 = [aux_norm, print.model.nodes(n).inclSuppDir(3)];
                    dir_2 = dir_2 / norm(dir_2);
                    thetaY = asin(dir_2(2));
                    
                    % Convert from rad to degrees
                    thetaY = thetaY * 180 / pi;
                    if abs(thetaY) < 10^-10
                        thetaY = 0;
                    end
                    
                    % Get rot X, from vy direction vector
                    if ~isempty(print.model.nodes(n).inclSupp_vy)
                        if abs(abs(thetaY)-90) < 10^-10 && abs(thetaZ) > 10^-10
                            thetaX = 0;
                        else
                            [x,dy1,~] = print.model.nodes(n).getInclinedSuppLocAxis(true);
                            [~,dy2,~] = print.model.nodes(n).getInclinedSuppLocAxis;
                            w = cross(dy1,dy2);
                            aux = x * w';
                            thetaX = acos((dy1*dy2')/(norm(dy1)*norm(dy2))) * 180 / pi;
                            if abs(thetaX) < 10^-10
                                thetaX = 0;
                            elseif aux < 0
                                thetaX = - thetaX;
                            end
                        end
                    else
                        thetaX = [];
                    end
                    
                    % Assemble strings to be printed
                    incl = 'YES';
                    
                    if ~isempty(thetaX)
                        if abs(thetaX) >= 100
                            rot_x = sprintf('%.1f',thetaX);
                        elseif abs(thetaX) >= 10
                            rot_x = sprintf('%.2f',thetaX);
                        else
                            rot_x = sprintf('%.3f',thetaX);
                        end
                    else
                        rot_x = '   --- ';
                    end
                    
                    if abs(thetaY) >= 100
                        rot_y = sprintf('%.1f',thetaY);
                    elseif abs(thetaY) >= 10
                        rot_y = sprintf('%.2f',thetaY);
                    else
                        rot_y = sprintf('%.3f',thetaY);
                    end
                    
                    if abs(thetaZ) >= 100
                        rot_z = sprintf('%.1f',thetaZ);
                    elseif abs(thetaZ) >= 10
                        rot_z = sprintf('%.2f',thetaZ);
                    else
                        rot_z = sprintf('%.3f',thetaZ);
                    end
                    
                    dir_x1 = sprintf('%.2f',print.model.nodes(n).inclSuppDir(1));
                    dir_x2 = sprintf('%.2f',print.model.nodes(n).inclSuppDir(2));
                    dir_x3 = sprintf('%.2f',print.model.nodes(n).inclSuppDir(3));
                    
                    if ~isempty(print.model.nodes(n).inclSupp_vy)
                        dir_y1 = sprintf('%.2f',print.model.nodes(n).inclSupp_vy(1));
                        dir_y2 = sprintf('%.2f',print.model.nodes(n).inclSupp_vy(2));
                        dir_y3 = sprintf('%.2f',print.model.nodes(n).inclSupp_vy(3));
                    else
                        dir_y1 = '  ---';
                        dir_y2 = '  ---';
                        dir_y3 = '  ---';
                    end
                end
                
                fprintf(print.txt, '%4d     %6s     %6s     %6s   |   %9s    %7s      %7s      %7s     %5s    %5s    %5s   %5s    %5s    %5s\n', n,...
                        node_restr1,...
                        node_restr2,...
                        node_restr3,...
                        incl,...
                        rot_x,...
                        rot_y,...
                        rot_z,...
                        dir_x1,...
                        dir_x2,...
                        dir_x3,...
                        dir_y1,...
                        dir_y2,...
                        dir_y3);
            end
        end
        
        %------------------------------------------------------------------
        % Prints spring supports.
        % Input arguments:
        % n_ns: number of nodes with spring supports
        function spring(print,n_ns)
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'S P R I N G  S T I F F N E S S\n');
            fprintf(print.txt, '--------------------------------\n');
            
            if n_ns ~= 0
                fprintf(print.txt, ' NODE    KDX [kN/m]      KDY [kN/m]      KDZ [kN/m]\n');

                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).springStiff)
                        fprintf(print.txt, '%4d   %12.3e   %13.3e   %13.3e\n', n,...
                                print.model.nodes(n).springStiff(1),...
                                print.model.nodes(n).springStiff(2),...
                                print.model.nodes(n).springStiff(3) );
                    end
                end
            else
                fprintf(print.txt, ' NO SPRINGS\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal loads.
        % Input arguments:
        %  n_nl: number of nodes with applied loads
        function nodalLoads(print,n_nl)
            fprintf(print.txt, '\n\n--------------------\n');
            fprintf(print.txt, 'N O D A L  L O A D S \n');
            fprintf(print.txt, '--------------------\n');
            
            if n_nl ~= 0
                aux = 0;
                fprintf(print.txt, ' NODE     FX [kN]          FY [kN]          FZ [kN]\n');
                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).load.static) && ~all(print.model.nodes(n).load.static == 0)
                        fprintf(print.txt, '%4d   %10.3f   %14.3f   %14.3f\n', n,...
                                print.model.nodes(n).load.static(1),...
                                print.model.nodes(n).load.static(2),...
                                print.model.nodes(n).load.static(3));
                        aux = aux + 1;
                        if aux == n_nl
                            break
                        end
                    end
                end
            else
                fprintf(print.txt, ' NO NODAL LOAD\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal prescribed displacements.
        % Input arguments:
        %  n_pd: number of nodes with prescribed displacement
        function nodalPrescDisp(print,n_pd)
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'N O D A L  P R E S C.  D I S P L. \n');
            fprintf(print.txt, '--------------------------------\n');
            
            if n_pd ~= 0
                aux = 0;
                fprintf(print.txt, ' NODE     DX [mm]          DY [mm]          DZ [mm]\n');
                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).prescDispl) && ~all(print.model.nodes(n).prescDispl == 0)
                        fprintf(print.txt, '%4d   %8.1f   %14.1f   %14.1f\n', n,...
                                1e3*print.model.nodes(n).prescDispl(1),...
                                1e3*print.model.nodes(n).prescDispl(2),...
                                1e3*print.model.nodes(n).prescDispl(3));
                        aux = aux + 1;
                        if aux == n_pd
                            break
                        end
                    end
                end
            else
                fprintf(print.txt, ' NO PRESCRIBED DISPLACEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints elements information.
        function elements(print)
            fprintf(print.txt, '\n\n---------------\n');
            fprintf(print.txt, 'E L E M E N T S\n');
            fprintf(print.txt, '---------------\n');

            if print.model.nel > 0
                fprintf(print.txt, ' ELEMENT      TYPE        MAT  SEC  JOINTi    JOINTf    NODEi NODEf  LENGTH [m]   vz_X     vz_Y    vz_Z\n');

                for e = 1:print.model.nel
                    type = 'Truss element ';
                    mat = print.model.elems(e).material.id;
                    sec = print.model.elems(e).section.id;
                    hingei = 'hinged';
                    hingef = 'hinged';
                    nodei = print.model.elems(e).nodes(1).id;
                    nodef = print.model.elems(e).nodes(2).id;
                    L = print.model.elems(e).length;
                    v = [print.model.elems(e).vz(1), print.model.elems(e).vz(2), print.model.elems(e).vz(3)];
                    v = v / norm(v);
                    
                    fprintf(print.txt, '%5d   %15s %4d %4d   %6s    %6s    %3d  %4d  %9.3f     %6.3f   %6.3f  %6.3f\n', ...
                            e, type, mat, sec, hingei, hingef, nodei, nodef, L, v(1), v(2), v(3));
                end
            else
                fprintf(print.txt, ' NO ELEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints semi-rigid joints information.
        function srjoints(~,~)
        end

        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        % Input arguments:
        %  n_ul: number of elements with uniformly distributed loads
        function unifElementLoads(print,n_ul)
            fprintf(print.txt, '\n\n----------------------------------------\n');
            fprintf(print.txt, 'U N I F O R M  E L E M E N T  L O A D S\n');
            fprintf(print.txt, '----------------------------------------\n');
            
            if n_ul ~= 0
                aux = 0;
                fprintf(print.txt, ' ELEMENT  DIRECTION    QX [kN/m]      QY [kN/m]     QZ [kN/m]\n');
                for e = 1:print.model.nel
                    if ~isempty(print.model.elems(e).load.uniformGbl) && ~all(print.model.elems(e).load.uniformGbl == 0)
                        if print.model.elems(e).load.uniformDir == 0
                            dir = 'GLOBAL';
                            qx = print.model.elems(e).load.uniformGbl(1);
                            qy = print.model.elems(e).load.uniformGbl(2);
                            qz = print.model.elems(e).load.uniformGbl(3);
                        else
                            dir = 'LOCAL';
                            qx = print.model.elems(e).load.uniformLcl(1);
                            qy = print.model.elems(e).load.uniformLcl(2);
                            qz = print.model.elems(e).load.uniformLcl(3);
                        end
                        fprintf(print.txt, '%5d  %9s %13.3f %13.3f  %13.3f\n',e,dir,qx,qy,qz);
                        aux = aux + 1;
                        if aux == n_ul
                            break
                        end
                    end
                end
            else
                fprintf(print.txt, ' NO UNIFORM LOAD\n');
            end
        end

        %------------------------------------------------------------------
        % Prints linearly distributed loads information.
        % Input arguments:
        %  n_ll: number of elements with linearly distributed loads
        function linearElementLoads(print,n_ll)
            fprintf(print.txt, '\n\n--------------------------------------\n');
            fprintf(print.txt, 'L I N E A R  E L E M E N T  L O A D S\n');
            fprintf(print.txt, '--------------------------------------\n');
            
            if n_ll ~= 0
                aux = 0;
                fprintf(print.txt, ' ELEMENT  DIRECTION     QXi [kN/m]     QXf [kN/m]     QYi [kN/m]     QYf [kN/m]     QZi [kN/m]     QZf [kN/m]\n');
                for e = 1:print.model.nel
                    if ~isempty(print.model.elems(e).load.linearGbl) && ~all(print.model.elems(e).load.linearGbl == 0)
                        if print.model.elems(e).load.linearDir == 0
                            dir = 'GLOBAL';
                            qxi = print.model.elems(e).load.linearGbl(1);
                            qyi = print.model.elems(e).load.linearGbl(2);
                            qzi = print.model.elems(e).load.linearGbl(3);
                            qxf = print.model.elems(e).load.linearGbl(4);
                            qyf = print.model.elems(e).load.linearGbl(5);
                            qzf = print.model.elems(e).load.linearGbl(6);
                        else
                            dir = 'LOCAL ';
                            qxi = print.model.elems(e).load.linearLcl(1);
                            qyi = print.model.elems(e).load.linearLcl(2);
                            qzi = print.model.elems(e).load.linearLcl(3);
                            qxf = print.model.elems(e).load.linearLcl(4);
                            qyf = print.model.elems(e).load.linearLcl(5);
                            qzf = print.model.elems(e).load.linearLcl(6);
                        end
                        fprintf(print.txt, '%5d   %9s %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f\n',...
                                e, dir, qxi, qxf, qyi, qyf, qzi, qzf);
                        aux = aux + 1;
                        if aux == n_ll
                            break
                        end
                    end
                end
            else
                fprintf(print.txt, ' NO LINEAR LOAD\n');
            end
        end

        %------------------------------------------------------------------
        % Prints thermal loads information.
        % Input arguments:
        %  n_tv: number of elements with temperature variation
        function temperatureVariation(print,n_tv)
            fprintf(print.txt, '\n\n-----------------------------------------\n');
            fprintf(print.txt, 'T E M P E R A T U R E  V A R I A T I O N\n');
            fprintf(print.txt, '-----------------------------------------\n');
            
            if n_tv ~= 0
                fprintf(print.txt, ' ELEMENT     dTx [°C]     dTy [°C]     dTz [°C]\n');

                for e = 1:print.model.nel
                    if (print.model.elems(e).load.tempVar_X ~= 0) || ...
                       (print.model.elems(e).load.tempVar_Y ~= 0) || ...
                       (print.model.elems(e).load.tempVar_Z ~= 0)
                        dtx = print.model.elems(e).load.tempVar_X;
                        dty = print.model.elems(e).load.tempVar_Y;
                        dtz = print.model.elems(e).load.tempVar_Z;

                        fprintf(print.txt, '%5d   %12.3f %12.3f %12.3f\n', e, dtx, dty, dtz);
                    end
                end
            else
                fprintf(print.txt, ' NO TEMPERATURE VARIATION\n');
            end
        end

        %------------------------------------------------------------------
        % Prints results of nodal displacement/rotation.
        function nodalDisplRot(print)
            fprintf(print.txt, '\n\n---------------------------------------------------------\n');
            fprintf(print.txt, 'N O D A L  D I S P L A C E M E N T S / R O T A T I O N S \n');
            fprintf(print.txt, '---------------------------------------------------------\n');
            fprintf(print.txt, ' NODE      DISPL X [mm]       DISPL Y [mm]       DISPL Z [mm]\n');

            for n = 1:print.model.nnp
                dx = 1e3*print.model.D(print.model.ID(1,n));
                if abs(dx) < 10^-15
                    dx = 0;
                end
                dy = 1e3*print.model.D(print.model.ID(2,n));
                if abs(dy) < 10^-15
                    dy = 0;
                end
                dz = 1e3*print.model.D(print.model.ID(3,n));
                if abs(dz) < 10^-15
                    dz = 0;
                end
                fprintf(print.txt, '%4d     %11.2e     %14.2e     %14.2e\n', n, dx, dy, dz);
            end
        end

        %------------------------------------------------------------------
        % Prints results of support reactions.
        function reactions(print)
            fprintf(print.txt, '\n\n---------------------------------\n');
            fprintf(print.txt, 'S U P P O R T  R E A C T I O N S\n');
            fprintf(print.txt, '---------------------------------\n');
            fprintf(print.txt, ' NODE       FORCE X [kN]     FORCE Y [kN]     FORCE Z [kN]\n');

            for n = 1:print.model.nnp
                if( (print.model.ID(1,n) > print.model.neqfree) || ...
                    (print.model.ID(2,n) > print.model.neqfree) || ...
                    (print.model.ID(3,n) > print.model.neqfree) )

                    if print.model.ID(1,n) > print.model.neqfree
                        node_reaction1 = print.model.F(print.model.ID(1,n));
                    else
                        node_reaction1 = 0.0;
                    end

                    if print.model.ID(2,n) > print.model.neqfree
                        node_reaction2 = print.model.F(print.model.ID(2,n));
                    else
                        node_reaction2 = 0.0;
                    end
                    
                    if print.model.ID(3,n) > print.model.neqfree
                        node_reaction3 = print.model.F(print.model.ID(3,n));
                    else
                        node_reaction3 = 0.0;
                    end

                    fprintf(print.txt, '%4d     %13.3f   %14.3f   %14.3f\n', n, ...
                            node_reaction1, node_reaction2, node_reaction3);
                end
            end
        end

        %------------------------------------------------------------------
        % Prints results of internal forces at element nodes.
        function intForces(print)
            fprintf(print.txt, '\n\n------------------------------------------------------------\n');
            fprintf(print.txt, 'I N T E R N A L  F O R C E S  A T  E L E M E N T  N O D E S\n' );
            fprintf(print.txt, '------------------------------------------------------------\n');
            fprintf(print.txt, ' ELEM      AXIAL FORCE [kN]\n');
            fprintf(print.txt, '           Nodei     Nodef\n');

            for e = 1:print.model.nel
                fprintf(print.txt, '%4d   %10.3f %9.3f\n', e,...
                        print.model.elems(e).axial_force(1),...
                        print.model.elems(e).axial_force(2));
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of elements internal displacements.
        function elemDispl(print)
            fprintf(print.txt, '\n\n-----------------------------------------------------------------------------------------\n');
            fprintf(print.txt, 'E L E M E N T S  I N T E R N A L  D I S P L A C E M E N T S  I N  L O C A L  S Y S T E M\n' );
            fprintf(print.txt, '-----------------------------------------------------------------------------------------\n');
            fprintf(print.txt, 'Axial and transversal displacements in 10 cross-sections from x = 0 to x = L\n\n');

            for e = 1:print.model.nel
                L = print.model.elems(e).length;
                l = zeros(10,1);
                l(1) = 0;
                for i = 2:10
                    l(i) = l(i-1) + L/9;
                end
                
                du = zeros(10,1);
                dv_XY = zeros(10,1);
                dv_XZ = zeros(10,1);
                j = 1;
                for x = 0:L/9:L
                    del_gblAnl = 1000*print.model.elems(e).gblAnlIntDispl(print.model,x);
                    del_lclAnl = 1000*print.model.anm.lclAnlIntDispl(print.model.elems(e),x);
                    d = del_gblAnl + del_lclAnl;
                    du(j) = d(1);
                    dv_XY(j) = d(2);
                    dv_XZ(j) = d(3);
                    j = j + 1;
                end
                
                fprintf(print.txt, ' ELEM %d\n', e);
                
                fprintf(print.txt, '   X [m]     %9.3f     %6.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f   %8.3f\n',...
                                                 l(1),     l(2),     l(3),    l(4),    l(5),    l(6),    l(7),    l(8),    l(9),   l(10));
                
                fprintf(print.txt, '  du [mm]    %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n',...
                                                  du(1),  du(2),  du(3),  du(4),  du(5),  du(6),  du(7),  du(8),  du(9),  du(10));
                                           
                fprintf(print.txt, ' dv_XY [mm]  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n',...
                                                dv_XY(1),dv_XY(2),dv_XY(3),dv_XY(4),dv_XY(5),dv_XY(6),dv_XY(7),dv_XY(8),dv_XY(9),dv_XY(10));
                                           
                fprintf(print.txt, ' dv_XZ [mm]  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n\n',...
                                                 dv_XZ(1),dv_XZ(2),dv_XZ(3),dv_XZ(4),dv_XZ(5),dv_XZ(6),dv_XZ(7),dv_XZ(8),dv_XZ(9),dv_XZ(10));
            end
        end
    end
end