%% Linear Elastic Static Driver Class
%
%% Description
%
% This is a sub-class of the <drv.html *Drv*> class for the
% implementation of the *Linear-Elastic Static* analysis type.
%
classdef Drv_LES < Drv
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function drv = Drv_LES(graphical_flag,model)
            include_constants;
            drv = drv@Drv(STATIC_LINEAR,graphical_flag,model);
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
            drv.model.F = sparse(drv.model.neq,1);
            drv.model.D = sparse(drv.model.neq,1);
        end
        
        %------------------------------------------------------------------
        % Assembles global matrices
        function gblMtx(drv)
            for e = 1:drv.model.nel
                % Compute element stiffness matrix in global system
                keg = drv.model.elems(e).gblStiffMtx();
                
                % Assemble element stiffness matrix to global matrix
                drv.model.assembleElemStiffMtx(keg,e);
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
            x = linspace(0,1,50);

            for e = 1:drv.model.nel
                % Internal forces on element ends
                %----------------------------------------------------------
                % Initialize element internal force arrays with null values
                drv.model.anm.initIntForce(drv.model.elems(e));
                
                % Compute element internal forces from global analysis
                % (from nodal displacements and rotations)
                fel_gblAnl = drv.model.elems(e).gblAnlIntForce(drv.model);
                
                % Get element internal forces from local analysis
                % (fixed end forces of distributed loads and thermal loads)
                fel_lclAnl = drv.model.elems(e).fel_distribLoad +...
                             drv.model.elems(e).fel_thermalLoad;
                
                % Add contribution of element internal forces from global
                % analysis and from local analysis
                fel = fel_gblAnl + fel_lclAnl;
                
                % Assemble element internal force arrays
                drv.model.anm.assembleIntForce(drv.model.elems(e),fel);
                %----------------------------------------------------------
                
                % Internal forces throughout element length
                %----------------------------------------------------------
                % Get element length
                L = drv.model.elems(e).length;
                
                % Compute element internal forces from global
                % analysis
                forceValues = drv.model.anm.intStress(drv.model.elems(e),x*L);
                
                % Assemble element internal force values matrix
                drv.model.elems(e).intStresses = forceValues;
                %----------------------------------------------------------
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements in several cross-section
        % positions.
        function elemIntDispl(drv)
            % Initialize step vector
            x = linspace(0,1,50);
            
            for e = 1:drv.model.nel
                % Get element length
                L = drv.model.elems(e).length;
                
                % Compute element internal displacements from global
                % analysis (from nodal displacements and rotations)
                del_gblAnl = drv.model.elems(e).gblAnlIntDispl(drv.model,x*L);
                
                % Compute element internal displacements from local
                % analysis (internal displacements from distributed
                % loads and thermal loads)
                del_lclAnl = drv.model.anm.lclAnlIntDispl(drv.model.elems(e),x*L);
                
                % Assemble element internal displacements matrix
                drv.model.elems(e).intDispl = del_gblAnl + del_lclAnl;
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
                waitbar(1/14);
            end
            
            % Dimension and initialize global matrix and vectors
            drv.dimMtxVctrs();
            if drv.graphical
                waitbar(2/14);
            end
            
            % Generate global d.o.f. numbering matrix
            drv.model.anm.setupDOFNum(drv.model);
            if drv.graphical
                waitbar(3/14);
            end
            
            % Assemble global d.o.f. numbering matrix
            drv.model.assembleDOFNum();
            if drv.graphical
                waitbar(4/14);
            end
            
            % Dimesion and generate srjoint objects
            drv.model.createSrjoints();
            if drv.graphical
                waitbar(5/14);
            end
            
            % Assemble element gather vector
            drv.model.assembleGatherVectors();
            if drv.graphical
                waitbar(6/14);
            end
            
            % Store prescribed displacements in global displacement vector
            drv.model.anm.setupPrescDispl(drv.model);
            if drv.graphical
                waitbar(7/14);
            end
            
            % Store spring stiffness coefficients in global stiffness matrix
            drv.model.anm.setupSpringStiff(drv.model);
            if drv.graphical
                waitbar(8/14);
            else
                fprintf(1,'Assembling stiffness matrix...\n');
            end
            
            % Assemble global stiffness matrix
            drv.gblMtx();
            if drv.graphical
                waitbar(9/14);
            else
                fprintf(1,'Assembling forcing vector...\n');
            end
            
            % Assemble global forcing vector
            drv.model.anm.nodalLoads(drv.model); % Initialize forcing vector with nodal loads
            if drv.graphical
                waitbar(10/14);
            end
            drv.model.elemLoads();         % Add element equivalent nodal loads to forcing vector
            if drv.graphical
                waitbar(11/14);
            else
                fprintf(1,'Solving system of equations...\n');
            end
            
            % Partition and solve the system of equations
            status = drv.solver.solve();
            % Rotate displacement vector back to global system, in case there were inclined supports
            drv.model.anm.locNodeToGblDispl(drv.model);
            if drv.graphical
                waitbar(12/14);
            end
            
            % Check structural model stability
            if (status == 1)
                if ~drv.graphical
                    fprintf(1,'Computing element internal forces...\n');
                end
                
                % Compute elements internal forces
                drv.elemIntForces();
                if drv.graphical
                    waitbar(13/14);
                else
                    fprintf(1,'Computing element internal displacements...\n');
                end
                
                % Compute elements internal displacements
                drv.elemIntDispl();
                if drv.graphical
                    waitbar(1);
                else
                    fprintf(1,'Structural analysis successfully finished.\n');
                end
                
                % Remove ficticious rotation constraints
                drv.model.fictRotConstraint(0);
            else
                % Remove ficticious rotation constraints
                drv.model.fictRotConstraint(0);
    
                if drv.graphical
                    msgbox('Unstable structure or invalid input data.', 'Error','error');
                else
                    fprintf(1,'Unstable structure or invalid input data.\n');
                end
            end
            
            % Close waitbar
            if drv.graphical
                close(h)
            end
        end
    end
end
