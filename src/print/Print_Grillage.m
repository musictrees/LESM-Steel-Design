%% Print_Grillage Class
%
%% Description
%
% This is a sub-class of the <print.html *Print*> class for the
% implementation of the *Grillage* print object.
%
classdef Print_Grillage < Print
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print_Grillage(model)
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
            fprintf(print.txt, 'GRILLAGE\n');
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
                for i = 3:5
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
                % Increment number of semi-rigid joints
                if print.model.elems(e).hingei == SEMIRIGID_END
                    n_srj = n_srj + 1;
                end
                if print.model.elems(e).hingef == SEMIRIGID_END
                    n_srj = n_srj + 1;
                end
                
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
            fprintf(print.txt, 'NUMBER OF SEMI-RIGID JOINTS...........:%4d\n', n_srj);
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
                fprintf(print.txt, ' SECTION     FULL AREA [cm²]     SHEAR AREA [cm²]     INERTIA X [cm4]      INERTIA Y [cm4]      HEIGHT [cm]\n');

                for s = 1:print.model.nsec
                    fprintf(print.txt, '%4d   %15.2f   %17.2f    %18.2f   %18.2f   %15.2f\n', s,...
                            1e4*print.model.sections(s).area_x,...
                            1e4*print.model.sections(s).area_z,...
                            1e8*print.model.sections(s).inertia_x,...
                            1e8*print.model.sections(s).inertia_y,...
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
            fprintf(print.txt, ' NODE   ROTAT X    ROTAT Y    DISPL Z\n');

            for n = 1:print.model.nnp
                if(print.model.nodes(n).ebc(4) == FIXED_DOF)
                    node_restr1 = 'FIXED';
                elseif (print.model.nodes(n).ebc(4) == SPRING_DOF)
                    node_restr1 = 'SPRING';
                else
                    node_restr1 = 'FREE';
                end

                if(print.model.nodes(n).ebc(5) == FIXED_DOF)
                    node_restr2 = 'FIXED';
                elseif (print.model.nodes(n).ebc(5) == SPRING_DOF)
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

                fprintf(print.txt, '%4d     %6s     %6s     %6s\n', n,...
                        node_restr1,...
                        node_restr2,...
                        node_restr3);
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
                fprintf(print.txt, ' NODE     KRX [kNm/rad]      KRY [kNm/rad]      KDZ [kN/m]\n');

                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).springStiff)
                        fprintf(print.txt, '%4d   %16.3e   %16.3e   %13.3e\n', n,...
                                print.model.nodes(n).springStiff(4),...
                                print.model.nodes(n).springStiff(5),...
                                print.model.nodes(n).springStiff(3));
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
                aux  = 0;
                fprintf(print.txt, ' NODE     MX [kNm]         MY [kNm]         FZ [kN]\n');
                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).load.static) && ~all(print.model.nodes(n).load.static == 0)
                        fprintf(print.txt, '%4d   %10.3f   %14.3f   %14.3f\n', n,...
                                print.model.nodes(n).load.static(4),...
                                print.model.nodes(n).load.static(5),...
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
                fprintf(print.txt, ' NODE     RX [rad]          RY [rad]          DZ [mm]\n');
                for n = 1:print.model.nnp
                    if ~isempty(print.model.nodes(n).prescDispl) && ~all(print.model.nodes(n).prescDispl == 0)
                        fprintf(print.txt, '%4d   %10.3f   %15.3f   %13.1f\n', n,...
                                print.model.nodes(n).prescDispl(1),...
                                print.model.nodes(n).prescDispl(2),...
                                1e3*print.model.nodes(n).prescDispl(6));
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
            include_constants;

            fprintf(print.txt, '\n\n---------------\n');
            fprintf(print.txt, 'E L E M E N T S\n');
            fprintf(print.txt, '---------------\n');

            if print.model.nel > 0
                etype = zeros(1,print.model.nel);
                for e = 1:print.model.nel
                    etype(e) = print.model.elems(e).type;
                end    
                    
                if all(etype(:) == MEMBER_NAVIER)
                    fprintf(print.txt, ' ELEMENT      TYPE        MAT  SEC    JOINTi        JOINTf      NODEi NODEf   LENGTH [m]   vz_X     vz_Y    vz_Z\n');
                else
                    fprintf(print.txt, ' ELEMENT      TYPE            MAT  SEC    JOINTi        JOINTf      NODEi NODEf   LENGTH [m]   vz_X     vz_Y    vz_Z\n');
                end
            
                for e = 1:print.model.nel
                    element_type = print.model.elems(e).type;

                    switch element_type
                        case MEMBER_NAVIER
                            type = 'Navier element';
                        case MEMBER_TIMOSHENKO
                            type = 'Timoshenko element';
                    end

                    if print.model.elems(e).hingei == HINGED_END
                        hingei = 'hinged    ';
                    elseif print.model.elems(e).hingei == SEMIRIGID_END
                        hingei = 'semi-rigid';
                    else
                        hingei = 'continuous';
                    end

                    if print.model.elems(e).hingef == HINGED_END
                        hingef = 'hinged    ';
                    elseif print.model.elems(e).hingef == SEMIRIGID_END
                        hingef = 'semi-rigid';
                    else
                        hingef = 'continuous';
                    end

                    mat = print.model.elems(e).material.id;
                    sec = print.model.elems(e).section.id;
                    nodei = print.model.elems(e).nodes(1).id;
                    nodef = print.model.elems(e).nodes(2).id;
                    L = print.model.elems(e).length;
                    v = [print.model.elems(e).vz(1), print.model.elems(e).vz(2), print.model.elems(e).vz(3)];
                    v = v / norm(v);
                    
                    if all(etype(:) == MEMBER_NAVIER) == 1
                        fprintf(print.txt, '%5d   %15s %4d %4d   %10s    %10s    %3d  %4d  %9.3f     %6.3f   %6.3f  %6.3f\n', ...
                            e, type, mat, sec, hingei, hingef, nodei, nodef, L, v(1), v(2), v(3));
                    else
                        fprintf(print.txt, '%5d   %19s %4d %4d   %10s    %10s    %3d  %4d  %9.3f     %6.3f   %6.3f  %6.3f\n', ...
                            e, type, mat, sec, hingei, hingef, nodei, nodef, L, v(1), v(2), v(3));
                    end    
                    
                end
            else
                fprintf(print.txt, ' NO ELEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints semi-rigid joints information.
        % Input arguments:
        %  n_srj: number of semi-rigid joints
        function srjoints(print,n_srj)
            include_constants;
            
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'S E M I - R I G I D  J O I N T S \n');
            fprintf(print.txt, '--------------------------------\n');
            
            if n_srj > 0
                aux = 0;
                fprintf(print.txt, ' ELEMENT   JOINTi_KRX [kNm/rad]  JOINTi_KRY [kNm/rad]  JOINTi_KRZ [kNm/rad]  JOINTf_KRX [kNm/rad]  JOINTf_KRY [kNm/rad]  JOINTf_KRZ [kNm/rad]\n');
                for e = 1:print.model.nel
                    
                    elemHasSrj = false;
                    
                    if print.model.elems(e).hingei == SEMIRIGID_END
                        
                        krxi = sprintf('%9.3e',print.model.elems(e).kri(1));
                        kryi = sprintf('%9.3e',print.model.elems(e).kri(2));
                        krzi = '   --    ';
                        
                        aux = aux + 1;
                        elemHasSrj = true;
                    else
                        krxi = '   --    ';
                        kryi = '   --    ';
                        krzi = '   --    ';
                    end
                    
                    if print.model.elems(e).hingef == SEMIRIGID_END
                        
                        krxf = sprintf('%9.3e',print.model.elems(e).krf(1));
                        kryf = sprintf('%9.3e',print.model.elems(e).krf(2));
                        krzf = '   --    ';
                        
                        aux = aux + 1;
                        elemHasSrj = true;
                    else
                        krxf = '   --    ';
                        kryf = '   --    ';
                        krzf = '   --    ';
                    end
                        
                        
                    if elemHasSrj
                        fprintf(print.txt, ' %4d            %s             %s             %s             %s             %s             %s       \n', e,...
                                krxi, kryi, krzi,...
                                krxf, kryf, krzf);
            
                        if aux >= n_srj
                            break
                        end
                    end
                end
            else
                fprintf(print.txt, ' NO SEMI-RIGID JOINTS\n');
            end
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
                aux  = 0;
                fprintf(print.txt, ' ELEMENT  DIRECTION     QZ [kN/m]\n');
                for e = 1:print.model.nel
                    if ~isempty(print.model.elems(e).load.uniformGbl) && ~all(print.model.elems(e).load.uniformGbl == 0)
                        dir = 'GLOBAL';
                        qz = print.model.elems(e).load.uniformGbl(3);
                        fprintf(print.txt, '%5d   %9s %13.3f\n',e,dir,qz);
                        aux = aux + 1;
                        if aux ==  n_ul
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
                fprintf(print.txt, ' ELEMENT  DIRECTION     QZi [kN/m]     QZf [kN/m]\n');
                for e = 1:print.model.nel
                    if ~isempty(print.model.elems(e).load.linearGbl) && ~all(print.model.elems(e).load.linearGbl == 0)
                        dir = 'GLOBAL';
                        qzi = print.model.elems(e).load.linearGbl(3);
                        qzf = print.model.elems(e).load.linearGbl(6);
                        fprintf(print.txt, '%5d   %9s %14.3f %14.3f\n',...
                                e, dir, qzi, qzf);
                        aux  = aux + 1;
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
                fprintf(print.txt, ' ELEMENT     dTx [°C]     dTz [°C]\n');

                for e = 1:print.model.nel
                    if (print.model.elems(e).load.tempVar_X ~= 0) || ...
                       (print.model.elems(e).load.tempVar_Z ~= 0)
                        dtx = print.model.elems(e).load.tempVar_X;
                        dtz = print.model.elems(e).load.tempVar_Z;

                        fprintf(print.txt, '%5d   %12.3f %12.3f\n', e, dtx, dtz);
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
            fprintf(print.txt, ' NODE      ROTAT X [rad]      ROTAT Y [rad]      DISPL Z [mm]\n');

            for n = 1:print.model.nnp
                rx = print.model.D(print.model.ID(1,n));
                if abs(rx) < 10^-15
                    rx = 0;
                end
                ry = print.model.D(print.model.ID(2,n));
                if abs(ry) < 10^-15
                    ry = 0;
                end
                dz = 1e3*print.model.D(print.model.ID(3,n));
                if abs(dz) < 10^-15
                    dz = 0;
                end
                fprintf(print.txt, '%4d     %14.5e     %14.5e     %13.2e\n', n, rx, ry, dz);
            end
        end

        %------------------------------------------------------------------
        % Prints results of support reactions.
        function reactions(print)
            fprintf(print.txt, '\n\n---------------------------------\n');
            fprintf(print.txt, 'S U P P O R T  R E A C T I O N S\n');
            fprintf(print.txt, '---------------------------------\n');
            fprintf(print.txt, ' NODE       MOMENT X [kNm]     MOMENT Y [kNm]    FORCE Z [kN]\n');

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

                    fprintf(print.txt, '%4d     %13.3f   %15.3f   %15.3f\n', n, ...
                            node_reaction1, node_reaction2, node_reaction3);
                end
            end
        end

        %------------------------------------------------------------------
        % Prints results of internal forces at element nodes.
        function intForces(print)
            fprintf(print.txt, '\n\n------------------------------------------------------------\n');
            fprintf(print.txt, 'I N T E R N A L  F O R C E S  A T  E L E M E N T  N O D E S\n');
            fprintf(print.txt, '------------------------------------------------------------\n');
            fprintf(print.txt, ' ELEM      SHEAR FORCE             BENDING MOMENT          TORSION MOMENT \n');
            fprintf(print.txt, '           Nodei      Nodef        Nodei      Nodef        Nodei      Nodef\n');

            for e = 1:print.model.nel
                fprintf(print.txt, '%4d   %10.3f %10.3f   %10.3f %10.3f   %10.3f %10.3f\n', e, ...
                        print.model.elems(e).shear_force_Z(1),...
                        print.model.elems(e).shear_force_Z(2), ...
                        print.model.elems(e).bending_moment_Y(1),...
                        print.model.elems(e).bending_moment_Y(2), ...
                        print.model.elems(e).torsion_moment(1),...
                        print.model.elems(e).torsion_moment(2));
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
                
                dv = zeros(10,1);
                j = 1;
                for x = 0:L/9:L
                    del_gblAnl = 1000*print.model.elems(e).gblAnlIntDispl(print.model,x);
                    del_lclAnl = 1000*print.model.anm.lclAnlIntDispl(print.model.elems(e),x);
                    d = del_gblAnl + del_lclAnl;
                    dv(j) = d(2);
                    j = j + 1;
                end
                
                fprintf(print.txt, ' ELEM %d\n', e);
                
                fprintf(print.txt, '   X [m]     %6.3f     %6.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f   %8.3f\n',...
                                                 l(1),     l(2),     l(3),    l(4),    l(5),    l(6),    l(7),    l(8),    l(9),   l(10));
                                             
                fprintf(print.txt, ' dv [mm]  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n\n',...
                                               dv(1), dv(2),  dv(3),  dv(4),  dv(5),  dv(6),  dv(7),  dv(8),  dv(9),  dv(10));
            end
        end
    end
end