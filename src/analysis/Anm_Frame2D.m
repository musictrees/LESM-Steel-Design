%% Anm_Frame2D Class
%
%% Description
%
% This is a sub-class of the <anm.html *Anm*> class for the
% implementation of the *2D Frame* analysis model.
%
% A 2D frame model has the following assumptions:
%
% * Frame elements are usually rigidly connected at joints.
%   However, a frame element might have a hinge (rotation liberation)
%   at an end or hinges at both ends.
% * It is assumed that a hinge in a 2D frame element releases continuity of
%   rotation in the out-of-plane direction.
% * A 2D frame model is considered to be laid in the XY-plane, with only
%   in-plane behavior, that is, there is no displacement transversal to
%   the frame plane.
% * Internal forces at any cross-section of a 2D frame element are:
%   axial force (local x direction), shear force (local y direction),
%   and bending moment (about local z direction).
% * Each node of a 2D frame model has three d.o.f.'s: a horizontal
%   displacement in X direction, a vertical displacement in Y
%   direction, and a rotation about the Z-axis.
%
classdef Anm_Frame2D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame2D()
            include_constants;
            anm = anm@Anm(FRAME2D_ANALYSIS,3,1,[1,2,6],3);
            
            anm.gln_axl    = [1,4];
            anm.gln_flx_XY = [2,3,5,6];
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <anm.html *Anm*>.
    methods
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function rot = gblToLocElemRotMtx(~,elem)
            % Get 3x3 basis rotation transformation matrix
            T = elem.T;
            
            % Assemble element d.o.f. rotation transformation matrix
            % rot = [ T 0
            %         0 T ]
            rot = blkdiag(T,T);
        end
        
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to the local inclined nodal support
        % system.
        % Output:
        %  rot: rotation transformation matrix
        %  flag: flag for there being inclined supports
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function [rot,flag] = gbltoLocNodeMtx(~,elem)
            
            % Initialize flag
            flag = true;
            
            % Get the info if the nodes have inclined support
            isInc1 = elem.nodes(1).isInclinedSupp;
            isInc2 = elem.nodes(2).isInclinedSupp;
            
            % Neither nodes have inclined supports
            if ~isInc1 && ~isInc2
                
                rot = 1;
                flag = false;
                return
                
            % The first node is inclined support and the second is not.
            elseif isInc1 && ~isInc2
                
                % Get node basis transformation matrix
                T1 = elem.nodes(1).T(1:2,1:2);
                T2 = eye(2);
                
            % The first node is not an inclined support and the second is.
            elseif ~isInc1 && isInc2
                
                % Get node basis transformation matrix
                T1 = eye(2);
                T2 = elem.nodes(2).T(1:2,1:2);
                
            % Both nodes have inclined support
            else
                
                % Get node basis transformation matrix
                T1 = elem.nodes(1).T(1:2,1:2);
                T2 = elem.nodes(2).T(1:2,1:2);
                
            end
            
            % Assemble element d.o.f. rotation transformation matrix
            % RZ = [ T1 0
            %        0 T2]
            rot = blkdiag(T1,1,T2,1);
        end
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint d.o.f. (degree of freedom) rotation
        % transformation matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        function rot = gblToLocSrjointRotMtx(~,~)
            rot = eye(2);
        end

        %------------------------------------------------------------------
        % Assembles element stiffness matrix in local system.
        % Output:
        %  kel: target element stiffness matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function kel = elemLocStiffMtx(anm,elem)
            % Initialize element stiffness matrix in local system
            kel = zeros(elem.nen * anm.ndof);
            
            % Compute axial stiffness coefficients
            kel(anm.gln_axl,anm.gln_axl) = elem.axialStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralStiffCoeff_XY();
        end
        
        %------------------------------------------------------------------
        % Assembles element mass matrix in local system.
        % Output:
        %  mel: target element mass matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function mel = elemLocMassMtx(anm,elem)
            % Initialize element mass matrix in local system
            mel = zeros(elem.nen * anm.ndof);
            
            % Compute axial mass coefficients
            mel(anm.gln_axl,anm.gln_axl) = elem.axialMassCoeff();
            
            % Compute flexural mass coefficients
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralMassCoeff_XY();
        end
        
        %------------------------------------------------------------------
        % Assembles element lumped mass matrix in local system.
        % Output:
        %  mel: target element mass matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function mel = elemLocLumpedMassMtx(anm,elem)
            % Initialize element lumped mass matrix in local system
            mel = zeros(elem.nen * anm.ndof);
            
            % Compute axial lumped mass coefficients
            mel(anm.gln_axl,anm.gln_axl) = elem.axialLumpedMassCoeff();
            
            % Compute flexural lumped mass coefficients
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralLumpedMassCoeff_XY();
        end
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint stiffness matrix in local system.
        % Output:
        %  kjl: target joint stiffness matrix in local system
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        function kjl = jointLocStiffMtx(~,srj)
            kjl = [  srj.krz  -srj.krz;
                    -srj.krz   srj.krz ];
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        function nodalLoads(~,model)
            for n = 1:model.nnp
                if ~isempty(model.nodes(n).load.static)
                    
                    % Check if node has inclined support and assemble
                    % forcing vector to be added
                    if model.nodes(n).isInclinedSupp
                        newF = model.nodes(n).T  * [model.nodes(n).load.static(1);...
                                                    model.nodes(n).load.static(2);...
                                                    0];
                    else
                        newF = [model.nodes(n).load.static(1);...
                                model.nodes(n).load.static(2)];
                    end
                    
                    % Add applied force in global X direction
                    id = model.ID(1,n);
                    model.F(id) = model.F(id) + newF(1);
                    
                    % Add applied force in global Y direction
                    id = model.ID(2,n);
                    model.F(id) = model.F(id) + newF(2);
                    
                    % Add applied moment about global Z direction
                    id = model.ID(3,n);
                    model.F(id) = model.F(id) + model.nodes(n).load.static(6);
                    
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing matrix.
        % Each column consists on a forcing vector on a time instant
        % Input arguments:
        %  model: handle to an object of the model class
        function dynamicNodalLoads(~,model)             
            % Loop through all nodes
            for n = 1:model.nnp
                % Check if node has dynamic loads
                if ~isempty(model.nodes(n).load.dynamic)
                    if ~all(model.nodes(n).load.dynamic == 0) && ~isempty(model.nodes(n).load.getFcn())
                        % Evalute dynamic loads
                        f = model.nodes(n).load.evalDynamicLoad();

                        % Check if node has inclined support and assemble
                        % forcing vector to be added
                        if model.nodes(n).isInclinedSupp
                            newF = model.nodes(n).T  * [f(1,:);...
                                                        f(2,:);...
                                                        zeros(1,model.n_steps+1)];
                        else
                            newF = [f(1,:);...
                                    f(2,:)];
                        end

                        % Add applied force in global X direction
                        id = model.ID(1,n);
                        model.F(id,:) = model.F(id,:) + newF(1,:);

                        % Add applied force in global Y direction
                        id = model.ID(2,n);
                        model.F(id,:) = model.F(id,:) + newF(2,:);

                        % Add applied moment about global Z direction
                        id = model.ID(3,n);
                        model.F(id,:) = model.F(id,:) + f(6,:);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute nodal displacement obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        function locNodeToGblDispl(~,model)
            for n = 1:model.nnp
                % Check if node has inclined support
                if model.nodes(n).isInclinedSupp

                % Get node displacement DOF ids
                id_1 = model.ID(1,n);
                id_2 = model.ID(2,n);
                id   = [id_1, id_2];
                
                % Rotate displacement to global system from local
                % inclined support system
                model.D(id) = (model.nodes(n).T(1:2,1:2))' * model.D(id);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute vibration modes in global system,
        % rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        %  v: [OPTIONAL] vibration modes matrix
        % Output:
        %  v: rotated vibration modes matrix
        function v = locNodeToGblVbrtnModes(~,model,v)
            if nargin < 3
              v = [model.V; zeros(model.neqfixed,model.n_modes)];
            end
            for n = 1:model.nnp
                % Check if node has inclined support
                if model.nodes(n).isInclinedSupp

                % Get node displacement DOF ids
                id_1 = model.ID(1,n);
                id_2 = model.ID(2,n);
                id   = [id_1, id_2];
                
                % Rotate displacement to global system from local
                % inclined support system
                v(id,:) = (model.nodes(n).T(1:2,1:2))' * v(id,:);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute nodal results obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        function locNodeToGblDynamicRes(~,model)
            for n = 1:model.nnp
                % Check if node has inclined support
                if model.nodes(n).isInclinedSupp

                % Get node displacement DOF ids
                id_1 = model.ID(1,n);
                id_2 = model.ID(2,n);
                id   = [id_1, id_2];
                
                % Rotate results to global system from local
                % inclined support system 
                for m = 1:size(model.results.dynamicDispl,3)
                    model.results.dynamicDispl(id,:,m)       = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicDispl(id,:,m);
                    model.results.dynamicVeloc(id,:,m)       = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicVeloc(id,:,m);
                    model.results.dynamicAccel(id,:,m)       = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicAccel(id,:,m);
                    model.results.dynamicDisplForced(id,:,m) = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicDisplForced(id,:,m);
                    model.results.dynamicVelocForced(id,:,m) = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicVelocForced(id,:,m);
                    model.results.dynamicAccelForced(id,:,m) = (model.nodes(n).T(1:2,1:2))' * model.results.dynamicAccelForced(id,:,m);
                end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied distributed load.
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        function fel = elemLocDistribLoadFEF(anm,load)
            % Initialize load vector
            fel = zeros(load.elem.nen*anm.ndof,1);
            
            % Compute axial fixed end force components
            fel(anm.gln_axl) = load.axialDistribLoadFEF();
            
            % Compute flexural (transversal) fixed end force components
            fel(anm.gln_flx_XY) = load.flexuralDistribLoadFEF_XY();
        end
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied thermal load (temperature variation).
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        function fel = elemLocThermalLoadFEF(anm,load)
            % Initialize load vector
            fel = zeros(load.elem.nen*anm.ndof,1);
            
            % Compute axial fixed end force components
            fel(anm.gln_axl) = load.axialThermalLoadFEF();
            
            % Compute flexural (transversal) fixed end force components
            fel(anm.gln_flx_XY) = load.flexuralThermalLoadFEF_XY();
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  axial_force(nsteps+1,2)
        %   Ni = axial_force(:,1) - init value throughout time steps
        %   Nf = axial_force(:,2) - final value throughout time steps
        %  shear_force(nsteps+1,2)
        %   Qi = shear_force(:,1) - init value throughout time steps
        %   Qf = shear_force(:,2) - final value throughout time steps
        %  bending_moment(nsteps+1,2)
        %   Mi = bending_moment(:,1) - init value throughout time steps
        %   Mf = bending_moment(:,2) - final value throughout time steps
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  nsteps: number of steps for transient analysis
        function initIntForce(~,elem,nsteps)
            if (nargin < 3)
                nsteps = 0;
            end
            elem.axial_force      = zeros(nsteps+1,elem.nen);
            elem.shear_force_Y    = zeros(nsteps+1,elem.nen);
            elem.bending_moment_Z = zeros(nsteps+1,elem.nen);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.axial_force(:,1)      = elem.axial_force(:,1)      + fel(1,:)';
            elem.axial_force(:,2)      = elem.axial_force(:,2)      + fel(4,:)';
            
            elem.shear_force_Y(:,1)    = elem.shear_force_Y(:,1)    + fel(2,:)';
            elem.shear_force_Y(:,2)    = elem.shear_force_Y(:,2)    + fel(5,:)';
            
            elem.bending_moment_Z(:,1) = elem.bending_moment_Z(:,1) + fel(3,:)';
            elem.bending_moment_Z(:,2) = elem.bending_moment_Z(:,2) + fel(6,:)';
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dv (transversal displacement in local y-axis)
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntDispl(~,elem)
            elem.intDispl = zeros(2,50);
        end
        
        %------------------------------------------------------------------
        % Assembles displacement shape function matrix evaluated at a 
        % given cross-section position.
        % Output:
        %  N: displacement shape function matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function N = displShapeFcnMtx(anm,elem,x)
            % Get number of points
            np = size(x,2);
            
            % Compute axial displacement shape functions vector
            Nu = elem.axialDisplShapeFcnVector(x);
            
            % Compute transversal displacement shape functions vector
            Nv = elem.flexuralDisplShapeFcnVector_XY(x);
            
            % Initialize displacement shape function matrix
            N = zeros(2*np,elem.nen * anm.ndof);
            
            % Assemble displacement shape function matrix
            N(2*(1:np)-1, anm.gln_axl)    = Nu;
            N(2*(1:np)  , anm.gln_flx_XY) = Nv;
        end
        
        %------------------------------------------------------------------
        % Computes internal displacements matrix from global analysis
        % Output:
        %  del: element internal displacements matrix
        % Input arguments:
        %  N: displacement shape function matrix
        %  dl: nodal displacements, on local coordinates
        function del = gblAnlIntDisplMtx(~,N,dl)
            % Get number of points
            np = size(N,1)/2;
            
            % Compute internal displacements matrix
            del_aux = N * dl;
            
            % Predimension del as a 3D matrix (uncoupled dynamic analysis)
            del = zeros(2,np,size(dl,2));
            
            % Rearrange internal displacements matrix
            for i = 1:size(dl,2)
                del(:,:,i) = [ (del_aux(2*(1:np)-1,i))';
                               (del_aux(2*(1:np),i))'  ];
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: a 2x1 vector of with element internal displacements:
        %      du -> axial displacement
        %      dv -> transversal displacement in local y-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dv(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Get number of points
            np = size(x,2);
            
            % Initialize displacements vector resulting from local analysis
            del = zeros(2,np);
            
            % Add the contribution of axial and transversal displacements
            % resulting from distributed loads
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(1,:) = elem.load.axialDistribLoadDispl(x);
                del(2,:) = elem.load.flexuralDistribLoadDispl_XY(x);
            end
            
            % Add the contribution of axial and transversal displacements
            % resulting from thermal loads
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Y ~= 0)
                del(1,:) = del(1,:) + elem.load.axialThermalLoadDispl(x);
                del(2,:) = del(2,:) + elem.load.flexuralThermalLoadDispl_XY(x);
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal stresses vector in local system,
        % in a given cross-section position.
        % Output:
        %  stressValues: element internal stresses vector in a given
        %                cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function stressValues = intStress(~,elem,x)
            % Get axial force values
            [N,elem.maxAxialForce] = elem.intAxialForce(x);
            
            % Get shear force values
            [Q,elem.maxShearForce_XY] = elem.intShearForce_XY(x);
            
            % Get bending moment values
            [M,elem.maxBendMoment_XY] = elem.intBendingMoment_XY(x);
            
            % Assemble stress values matrix
            stressValues = [N;
                            Q;
                            M];
        end
        
        %------------------------------------------------------------------
        % Computes element internal stresses vector in local system,
        % in a given cross-section position, throughout time steps.
        % Output:
        %  stressValues: element internal stresses vector in a given
        %                cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function stressValues = intDynamicStress(~,elem,x)
            % Get axial force values
            [N,elem.maxAxialForce] = elem.intDynamicAxialForce(x);
            
            % Get shear force values
            [Q,elem.maxShearForce_XY] = elem.intDynamicShearForce_XY(x);
            
            % Get bending moment values
            [M,elem.maxBendMoment_XY] = elem.intDynamicBendingMoment_XY(x);
            
            % Assemble stress values matrix
            stressValues(:,:,1) = N;
            stressValues(:,:,2) = Q;
            stressValues(:,:,3) = M;
        end
    end
end
