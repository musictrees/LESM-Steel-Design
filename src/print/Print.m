%% Print Class
%
%% Description
%
% This is a handle super-class for the definition of a printing object.
%
% Essentially, this super-class declares abstract methods and define
% public methods that print model information and analysis results in a
% text file or in the default output (MATLAB command window).
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <print_truss2d.html Print 2D truss model and results>.
% * <print_frame2d.html Print 2D frame model and results>.
% * <print_grillage.html Print grillage model and results>.
% * <print_truss3d.html Print 3D truss model and results>.
% * <print_frame3d.html Print 3D frame model and results>.
%
classdef Print < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        txt = 0;    % output identifier (1 = MATLAB command window)
        model = [];   % handle to an object of the model class
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print(txt)
            print.txt = txt;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Prints static analyis results.
        results(print)
        
        %------------------------------------------------------------------
        % Prints model type.
        modelLabel(print)
        
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
        [n_srj,n_ns,n_nl,n_pd,n_ul,n_ll,n_tv] = modelDescrip(print,lc,currentLc)
        
        %------------------------------------------------------------------
        % Prints cross-section properties.
        section(print)
        
        %------------------------------------------------------------------
        % Prints nodal support conditions.
        nodalSupport(print)
        
         %------------------------------------------------------------------
        % Prints spring information.
        % Input arguments:
        %  n_ns: number of nodes with applied loads
        spring(print,n_ns)
        
        %------------------------------------------------------------------
        % Prints nodal loads.
        % Input arguments:
        %  n_nl: number of nodes with applied loads
        nodalLoads(print,n_nl)
        
        %------------------------------------------------------------------
        % Prints nodal prescribed displacements.
        % Input arguments:
        %  n_pd: number of nodes with prescribed displacement
        nodalPrescDisp(print,n_pd)
        
        %------------------------------------------------------------------
        % Prints elements information.
        elements(print)
        
        %------------------------------------------------------------------
        % Prints semi-rigid joints information.
        % Input arguments:
        %  n_srj: number of semi-rigid joints
        srjoints(print,n_srj)
        
        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        % Input arguments:
        %  n_ul: number of elements with uniformly distributed loads
        unifElementLoads(print,n_ul)
        
        %------------------------------------------------------------------
        % Prints linearly distributed loads information.
        % Input arguments:
        %  n_ll: number of elements with linearly distributed loads
        linearElementLoads(print,n_ll)
        
        %------------------------------------------------------------------
        % Prints thermal loads information.
        % Input arguments:
        %  n_tv: number of elements with temperature variation
        temperatureVariation(print,n_tv)
        
        %------------------------------------------------------------------
        % Prints results of nodal displacement/rotation.
        nodalDisplRot(print)
        
        %------------------------------------------------------------------
        % Prints results of support reactions.
        reactions(print)
        
        %------------------------------------------------------------------
        % Prints results of internal forces at element nodes.
        intForces(print)
        
        %------------------------------------------------------------------
        % Prints results of elements internal displacements.
        elemDispl(print)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Prints header of analysis results.
        function header(print)
            fprintf(print.txt, '\n=========================================================\n');
            fprintf(print.txt, ' LESM - Linear Elements Structure Model analysis program\n');
            fprintf(print.txt, '    PONTIFICAL CATHOLIC UNIVERSITY OF RIO DE JANEIRO\n');
            fprintf(print.txt, '    DEPARTMENT OF CIVIL AND ENVIRONMENTAL ENGINEERING\n');
            fprintf(print.txt, '                          AND\n');
            fprintf(print.txt, '               TECGRAF/PUC-RIO INSTITUTE\n');
            fprintf(print.txt, '=========================================================\n');
        end
        
        %------------------------------------------------------------------
        % Prints material properties.
        function material(print)
            fprintf(print.txt, '\n\n-------------------------------------\n');
            fprintf(print.txt, 'M A T E R I A L  P R O P E R T I E S\n');
            fprintf(print.txt, '-------------------------------------\n');
            if print.model.nmat > 0
                fprintf(print.txt, ' MATERIAL     E [MPa]      G [MPa]     POISSON      THERMAL EXP. COEFF. [/°C]\n');
                
                for m = 1:print.model.nmat
                    fprintf(print.txt, '%4d   %13.0f   %10.0f   %8.2f   %20.3d\n', m,...
                            1e-3*print.model.materials(m).elasticity,...
                            1e-3*print.model.materials(m).shear,...
                            print.model.materials(m).poisson,...
                            print.model.materials(m).thermExp);
                end
            else
                fprintf(print.txt, ' NO MATERIAL\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal coordinates.
        function nodalCoords(print)
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'N O D A L  C O O R D I N A T E S \n');
            fprintf(print.txt, '--------------------------------\n');
            fprintf(print.txt, ' NODE     COORD X [m]      COORD Y [m]      COORD Z [m]    \n');
            
            for n = 1:print.model.nnp
                fprintf(print.txt, '%4d   %11.3f   %14.3f   %14.3f\n', n,...
                        print.model.nodes(n).coord(1),...
                        print.model.nodes(n).coord(2),...
                        print.model.nodes(n).coord(3));
            end
        end
        
        %------------------------------------------------------------------
        % Prints modal dynamic analyis results.
        function results_modalDynamic(print)
            print.header();
            fprintf(print.txt, '\n\n\n____________ M O D E L  I N F O R M A T I O N ____________\n');
            print.modelLabel();
            
            fprintf(print.txt, '\n\n----------------------------------------\n' );
            fprintf(print.txt, 'M O D E L  D E S C R I P T I O N:\n' );
            fprintf(print.txt, '----------------------------------------\n');
            fprintf(print.txt, 'NUMBER OF NODES.......................:%4d\n', print.model.nnp);
            fprintf(print.txt, 'NUMBER OF ELEMENTS ...................:%4d\n', print.model.nel);
            fprintf(print.txt, 'NUMBER OF DEGREES OF FREEDOM..........:%4d\n', print.model.neq);
            fprintf(print.txt, 'NUMBER OF FREE DEGREES OF FREEDOM.....:%4d\n', print.model.neqfree);
            fprintf(print.txt, 'NUMBER OF FIXED DEGREES OF FREEDOM....:%4d\n', print.model.neqfixed);
            fprintf(print.txt, 'NUMBER OF SPRINGS.....................:%4d\n', print.model.neqspring);
            fprintf(print.txt, 'NUMBER OF SEMI-RIGID JOINTS...........:%4d\n', print.model.njoints);
            fprintf(print.txt, 'NUMBER OF MATERIALS...................:%4d\n', print.model.nmat);
            fprintf(print.txt, 'NUMBER OF CROSS-SECTIONS..............:%4d\n', print.model.nsec);
            
            print.material();
            print.section();
            print.nodalCoords();
            print.nodalSupport();
            print.spring(print.model.neqspring);
            print.elements();
            print.srjoints(print.model.njoints);
            
            fprintf(print.txt, '\n\n\n\n_____________ A N A L Y S I S  R E S U L T S _____________\n');
            fprintf(print.txt, '\n\n---------------------------------\n');
            fprintf(print.txt, '   M O D A L  A N A L Y S I S \n');
            fprintf(print.txt, '---------------------------------\n');
            fprintf(print.txt, ' MODE       EIGENVALUE      w [rad/s]       f [Hz]          T [s]       VIBRATION MODE\n');
            for i = 1:length(print.model.W)
                fprintf(print.txt, '%4d       %10.5e     %9.5e    %8.5e    %7.5e       ',i,...
                        print.model.W(i)^2, print.model.W(i), print.model.W(i)/(2*pi), (2*pi)/print.model.W(i));
                for j = 1:size(print.model.V,2)
                    fprintf(print.txt,'  %5.3e  ',print.model.V(j,i));
                end
                fprintf(print.txt,'\n');
            end
            fprintf(print.txt, '\n\n---------------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Print object.
        function clean(print)
            print.txt = [];
            print.model = [];
        end
    end
end