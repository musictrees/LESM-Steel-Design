%% Linear Elastic Dynamic Driver Class
%
%% Description
%
% This is a sub-class of the <drv.html *Drv*> class for the
% implementation of the *Linear-Elastic Dynamic* analysis type.
%
classdef Drv_LED < Drv
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        whichResponse = 1;  % flag for response to be obtained by dynamic analysis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function drv = Drv_LED(solver_flag,graphical_flag,model)
            drv = drv@Drv(solver_flag,graphical_flag,model);
        end
    end
    
    %% Private methods
    % Implementation of methods that are exclusively called by dynamic
    % analysis modules.
    methods (Access = private)
        %------------------------------------------------------------------
        function status = eigenValues(drv)
            % Assemble free-free global stiffness matrix
            Kff = drv.model.K(1:drv.model.neqfree+drv.model.neqspring,...
                              1:drv.model.neqfree+drv.model.neqspring);
            
            % Check for instablity
            status = 1;
            if (rcond(full(Kff)) < 10e-12)
                status = 0;
                return;
            end
            
            % Assemble free-free global mass matrix
            Mff = drv.model.M(1:drv.model.neqfree+drv.model.neqspring,...
                              1:drv.model.neqfree+drv.model.neqspring);
            
            % Compute eigenvectors and eigenvalues
            %--------------------------------------------------------------
            % ATTENTION
            % The following code would be more efficient, but does not work 
            % on all recent Matlab versions
            %
            %[v,w2] = eigs(full(Kff),full(Mff),drv.n_modes,'smallestabs');
            %[drv.W,i] = sort(sqrt(diag(w2)));
            %drv.V = v(:,i);
            %
            % Instead of eigs, for now, the program uses eig, gets all
            % eigenvalues, sorts them and stores the requested number of
            % modes (drv.n_modes)
            %--------------------------------------------------------------
            [v,w2] = eig(full(Kff),full(Mff));
            [w,i]  = sort(sqrt(diag(w2)));
            drv.model.W  = w(1:drv.model.n_modes);
            drv.model.V  = v(:,i(1:drv.model.n_modes));
            
            % Normalize eigenvectors
            for i = 1:size(drv.model.V,2)
                drv.model.V(:,i) = drv.model.V(:,i) / norm(drv.model.V(:,i));
            end
        end
        
        %------------------------------------------------------------------
        function dampingMtx(drv)
            % Compute damping coefficients of propotionality to mass and
            % stiffness
            drv.computeDampingCoeffs();

            % Computes damping matrix
            drv.model.C = drv.model.massDampCoeff  * drv.model.M +...
                          drv.model.stiffDampCoeff * drv.model.K;
        end
        
        %------------------------------------------------------------------
        % Computes damping coefficients of propotionality to mass and
        % stiffness.
        % There are 3 possible cases:
        % * Coefficients are already provided (method does nothing)
        % * Only 1st mode critical damping ratio is provided (solve system)
        % * 1st and 2nd modes critical damping ratio are provided (solve 
        %   system)
        function computeDampingCoeffs(drv)
            include_constants;
            
            % Check if rayleigh coefficients were provided already
            if drv.model.damping == RAYLEIGH_COEFFS
                return
            end
            
            % Check if only one vibration mode was found. If so, consider
            % only mass proportional coefficient
            if length(drv.model.W) == 1
                drv.model.massDampCoeff  = 2 * drv.model.W * drv.model.xi(1);
                drv.model.stiffDampCoeff = 0;
            else
                % Assemble critical damping vector and matrix of first two
                % angular frequencies
                Xi = drv.model.xi'.*ones(2,1);
                A  =  0.5 * [ 1/drv.model.W(1)   drv.model.W(1) ;
                              1/drv.model.W(2)   drv.model.W(2) ];
                          
                if (rcond(full(A)) < 10e-12)
                    drv.model.massDampCoeff  = 2 * drv.model.W(1) * drv.model.xi(1);
                    drv.model.stiffDampCoeff = 0;
                    return;
                end

                % Compute damping coefficients
                coeffs = A \ Xi;
                drv.model.massDampCoeff  = coeffs(1);
                drv.model.stiffDampCoeff = coeffs(2);
            end
        end
        
        %------------------------------------------------------------------
        function initCondMtx(drv)
            % Initialize initial conditions matrix
            drv.model.c0 = zeros(drv.model.neq,2);
            
            % Set nodal initial conditions to global initial conditions mtx
            drv.model.anm.initialConditions(drv.model);
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements related to normalized
        % eigenvectors of vibration modes
        function elemIntVbrtnMode(drv)
            % Initialize step vector
            aux_x = linspace(0,1,50);
            x = aux_x([1, 2:2:48, 49, 50]);
            clear aux_x
            
            % Get eigenvectors matrix, already appended with fixed dofs and
            % rotated (inclined supports)
            v = drv.model.anm.locNodeToGblVbrtnModes(drv.model);
            
            for e = 1:drv.model.nel
                % Get element length
                L = drv.model.elems(e).length;
                
                % Compute element normalizad internal displacements
                % from each found vibration mode
                drv.model.elems(e).vibrtnModes(v,x*L);
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <drv.html *Drv*>.
    methods
        %------------------------------------------------------------------
        % Dimensions and initializes global matrices and vectors needed for
        % the current type of analysis
        function dimMtxVctrs(drv)
            drv.model.K = sparse(drv.model.neq,drv.model.neq);
            drv.model.M = sparse(drv.model.neq,drv.model.neq);
            drv.model.F = sparse(drv.model.neq,drv.model.n_steps + 1);
        end
        
        %------------------------------------------------------------------
        % Assembles global matrices
        function gblMtx(drv)
            for e = 1:drv.model.nel
                % Compute element stiffness matrix in global system
                keg = drv.model.elems(e).gblStiffMtx();
                
                % Assemble element stiffness matrix to global stiffness matrix
                drv.model.assembleElemStiffMtx(keg,e);
                
                % Check if this element should have its mass considered
                if (drv.model.elems(e).mass_consideration)
                    
                % Compute element mass matrix in global system
                drv.model.elems(e).mass_type = drv.model.mass_type;
                drv.model.elems(e).mass_mi = drv.model.mass_mi;
                meg = drv.model.elems(e).gblMassMtx();                
                               
                % Assemble element mass matrix to global mass matrix
                drv.model.assembleElemMassMtx(meg,e);
                
                end
            end
            
            for j = 1:drv.model.njoints
                % Compute semi-rigid joint stiffness matrix in global system
                kjg = drv.model.srjoints(j).gblStiffMtx();
                
                % Assemble semi-rigid joint stiffness matrix to global matrix
                drv.model.assembleSrjointStiffMtx(kjg,j);
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal forces in several cross-section
        % positions.
        function elemIntForces(drv)
            % Initialize step vector
            aux_x = linspace(0,1,50);
            x = aux_x([1, 2:2:48, 49, 50]);
            clear aux_x

            for e = 1:drv.model.nel
                % Internal forces on element ends
                %----------------------------------------------------------
                % Initialize element internal force arrays with null values
                drv.model.anm.initIntForce(drv.model.elems(e),drv.model.n_steps);
                
                % Compute element internal forces from global analysis
                % (from nodal displacements and rotations)
                fel = drv.model.elems(e).gblDynamicAnlIntForce(drv.model);
                
                % Assemble element internal force arrays
                drv.model.anm.assembleIntForce(drv.model.elems(e),fel);
                %----------------------------------------------------------
                
                % Internal forces throughout element length
                %----------------------------------------------------------
                % Get element length
                L = drv.model.elems(e).length;
                
                % Compute element internal forces from global
                % analysis
                forceValues = drv.model.anm.intDynamicStress(drv.model.elems(e),x*L);
                
                % Assemble element internal force values matrix
                drv.model.elems(e).intStresses = forceValues;
                %----------------------------------------------------------
                
                % Compute elements internal forces envelop
                drv.model.elems(e).forcesEnvelop();
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements in several cross-section
        % positions.
        function elemIntDispl(drv)
            % Initialize step vector
            aux_x = linspace(0,1,50);
            x = aux_x([1, 2:2:48, 49, 50]);
            clear aux_x
            
            % Sum of dynamic displacement over all modes
            d = sum(drv.model.results.dynamicDispl,3);
            
            for e = 1:drv.model.nel
                % Get element length
                L = drv.model.elems(e).length;
                
                % Compute element internal dynamic displacements
                drv.model.elems(e).dynamicDispl(d,x*L);
            end
        end
        
        %------------------------------------------------------------------
        % Processes current model data according to the current analysis
        % type
        % Considers the direct stiffness method.
        % Assembles global equilibrium system of equations, calculates
        % nodal displacements and rotations, and computes support reactions
        % and elements internal forces.
        % Graphical version.
        % Displays the processing stages in a waitbar and the error
        % message on a message box.
        % Output:
        %  status: flag for stable free-free global matrix (0 = unstable, 1 = stable)
        function status = process(drv)
            include_constants;
            
            % Check if provided materials have non-null density (TEMPORARY)
            for aux = 1:length(drv.model.materials)
                if drv.model.materials(aux).density <= 0
                    if drv.graphical
                        msgbox('Invalid material property. Density must be positive.', 'Error','error');
                    else
                        fprintf(1,'Invalid material property. Density must be positive.\n');
                    end
                    status = false;
                    return
                end
            end
            
            % Check if there is a Result object. If not, create one
            if isempty(drv.model.results)
                drv.model.results = Result();
            else
                % Clean previous results
                drv.model.results.clean();
            end
            drv.model.results.type = drv.analysis;
            
            % Create a wait bar to show the progress of the analysis process
            if drv.graphical
                h = waitbar(0,'Processing data...','WindowStyle','modal');
            else
                fprintf(1,'Preparing analysis data...\n');
            end
            
            % Create ficticious rotation constraints
            drv.model.fictRotConstraint(1);
            
            % Compute number of equations
            drv.model.computeNeq();
            if drv.graphical
                waitbar(1/15);
            end
            
            % Dimension and initialize global stiffness and mass matrices
            drv.dimMtxVctrs();
            if drv.graphical
                waitbar(2/15);
            end
            
            % Generate global d.o.f. numbering matrix
            drv.model.anm.setupDOFNum(drv.model);
            if drv.graphical
                waitbar(3/15);
            end
            
            % Check if all dofs are fixed (no dynamic analysis would be possible)
            if drv.model.neqfixed == drv.model.neq
                status = false;
                if drv.graphical
                    msgbox('All DOFs are fixed. Dynamic analysis cannot be performed', 'Error','error');
                else
                    fprintf(1,'All DOFs are fixed. Dynamic analysis cannot be performed.\n');
                end
                return
            elseif drv.model.n_modes <= 0 || drv.model.n_modes > drv.model.neqfree + drv.model.neqspring
                drv.model.n_modes = drv.model.neqfree + drv.model.neqspring;
            end
            
            % Assemble global d.o.f. numbering matrix
            drv.model.assembleDOFNum();
            if drv.graphical
                waitbar(4/15);
            end
            
            % Dimesion and generate srjoint objects
            drv.model.createSrjoints();
            if drv.graphical
                waitbar(5/15);
            end
            
            % Assemble element gather vector
            drv.model.assembleGatherVectors();
            if drv.graphical
                waitbar(6/15);
            end
            
            % Store spring stiffness coefficients in global stiffness matrix
            drv.model.anm.setupSpringStiff(drv.model);
            
            % Store concentrated mass coefficients in global mass matrix
            drv.model.anm.setupConcentratedMass(drv.model);
            
            % Assemble global stiffness and mass matrices
            if drv.graphical
                waitbar(7/15,h,'Assembling stiffness and mass matrices...');
            else
                fprintf(1,'Assembling stiffness and mass matrices...\n');
            end
            drv.gblMtx();
            
            % Compute eigenvalues and eigenvectors
            if drv.graphical
                waitbar(8/15,h,'Computing eigenvalues...');
            else
                fprintf(1,'Computing eigenvalues...\n');
            end
            status = drv.eigenValues();
            
            % Check structural model stability
            if status
                if drv.graphical
                    waitbar(9/15,h,'Computing vibration modes...');
                else
                    fprintf(1,'Computing vibration modes...\n');
                end
                
                % Compute elem internal displacement due to natural frequencies
                drv.elemIntVbrtnMode();
                
                if drv.model.n_steps > 0 && drv.whichResponse ~= MODAL_ANALYSIS
                    if drv.graphical
                        waitbar(10/15,h,'Assembling dynamic loads...');
                    else
                        fprintf(1,'Assembling dynamic loads...\n');
                    end
                    
                    % Assemble global forcing matrix
                    drv.model.anm.dynamicNodalLoads(drv.model);
                    if drv.graphical
                        waitbar(11/15,h,'Assembling damping matrix...');
                    else
                        fprintf(1,'Assembling damping matrix...\n');
                    end
            
                    % Assemble damping matrix
                    drv.dampingMtx();
                    if drv.graphical
                        waitbar(12/15,h,'Assembling initial conditions matrix...');
                    else
                        fprintf(1,'Assembling initial conditions matrix...\n');
                    end
                    
                    % Assemble initial conditions matrix
                    drv.initCondMtx();
                    if drv.graphical
                        waitbar(13/15,h,'Solving system of transient equations...');
                    else
                        fprintf(1,'Solving system of transient equations...\n');
                    end
                    
                    % Switch between solvers for dynamic analysis
                    [~] = drv.solver.solve();
                    % Rotate results back to global system, in case there were inclined supports
                    drv.model.anm.locNodeToGblDynamicRes(drv.model);
                    
                    if drv.graphical
                        waitbar(14/15,h,'Computing internal forces and displacements...');
                    else
                        fprintf(1,'Computing element internal forces and displacements...\n');
                    end
                    
                    % Compute elements internal forces
                    drv.elemIntForces();
                    
                    % Compute elem internal displacement due to forced vibrations
                    drv.elemIntDispl();
                    
                    if drv.graphical
                        waitbar(1);
                    else
                        fprintf(1,'Dynamic structural analysis successfully finished.\n');
                    end
                end
                
                % Remove ficticious rotation constraints
                drv.model.fictRotConstraint(0);
            else
                if drv.graphical
                    msgbox('Unstable structure or invalid input data.', 'Error','error');
                else
                    fprintf(1,'Unstable structure or invalid input data.\n');
                end
                
                % Remove ficticious rotation constraints
                drv.model.fictRotConstraint(0);
            end
            
            % Close waitbar
            if drv.graphical
                close(h)
            end
        end
    end
end
