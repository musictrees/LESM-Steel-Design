%% Anm_Frame3D Class
%
%% Description
%
% This is a sub-class of the <anm.html *Anm*> class for the
% implementation of the *3D Frame* analysis model.
%
% A 3D frame model has the following assumptions:
%
% * Frame elements are usually rigidly connected at joints.
%   However, a frame element might have a hinge (rotation liberation)
%   at an end or hinges at both ends.
% * It is assumed that a hinge in a 3D frame element releases continuity
%   of rotation in all directions.
% * Internal forces at any cross-section of a 3D frame element are:
%   axial force, shear forces (local y direction and local z direction),
%   bending moments (about local y direction and local z direction),
%   and torsion moment (about x direction).
% * Each node of a 3D frame model has six d.o.f.'s: displacements in
%   X, Y and Z directions, and rotations about X, Y and Z directions.
%
classdef Anm_Frame3D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame3D()
            include_constants;
            anm = anm@Anm(FRAME3D_ANALYSIS,6,3,1:6,4:6);
            
            anm.gln_axl    = [1,7];
            anm.gln_flx_XY = [2,6,8,12];
            anm.gln_flx_XZ = [3,5,9,11];
            anm.gln_tor    = [4,10];
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
            % rot = [ T 0 0 0
            %         0 T 0 0
            %         0 0 T 0
            %         0 0 0 T ]
            rot = blkdiag(T,T,T,T);
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
                T1 = elem.nodes(1).T;
                T2 = eye(3);
                
            % The first node is not an inclined support and the second is.
            elseif ~isInc1 && isInc2
                
                % Get node basis transformation matrix
                T1 = eye(3);
                T2 = elem.nodes(2).T;
                
            % Both nodes have inclined support
            else
                
                % Get node basis transformation matrix
                T1 = elem.nodes(1).T;
                T2 = elem.nodes(2).T;
                
            end
            
            % Assemble element d.o.f. rotation transformation matrix
            % RZ = [ T1 0
            %        0 T2]
            rot = blkdiag(T1,eye(3),T2,eye(3));
        end
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint d.o.f. (degree of freedom) rotation
        % transformation matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        function rot = gblToLocSrjointRotMtx(~,srj)
            % Get 3x3 basis rotation transformation matrix
            T = srj.elem.T;
            
            % Assemble element d.o.f. rotation transformation matrix
            % rot = [ T 0
            %         0 T ]
            rot = blkdiag(T,T);
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
            
            % Compute torsion stiffness coefficients
            kel(anm.gln_tor,anm.gln_tor) = elem.torsionStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralStiffCoeff_XY();
            kel(anm.gln_flx_XZ,anm.gln_flx_XZ) = elem.flexuralStiffCoeff_XZ();
            
%  kel = [ kea(1,1)  0            0           0         0            0           kea(1,2)  0            0           0         0            0;
%          0         kef_XY(1,1)  0           0         0            kef_XY(1,2) 0         kef_XY(1,3)  0           0         0            kef_XY(1,4);
%          0         0            kef_XZ(1,1) 0         kef_XZ(1,2)  0           0         0            kef_XZ(1,3) 0         kef_XZ(1,4)  0;
%          0         0            0           ket(1,1)  0            0           0         0            0           ket(1,2)  0            0;
%          0         0            kef_XZ(2,1) 0         kef_XZ(2,2)  0           0         0            kef_XZ(2,3) 0         kef_XZ(2,4)  0;
%          0         kef_XY(2,1)  0           0         0            kef_XY(2,2) 0         kef_XY(2,3)  0           0         0            kef_XY(2,4);
%          kea(2,1)  0            0           0         0            0           kea(2,2)  0            0           0         0            0;
%          0         kef_XY(3,1)  0           0         0            kef_XY(3,2) 0         kef_XY(3,3)  0           0         0            kef_XY(3,4);
%          0         0            kef_XZ(3,1) 0         kef_XZ(3,2)  0           0         0            kef_XZ(3,3) 0         kef_XZ(3,4)  0;
%          0         0            0           ket(2,1)  0            0           0         0            0           ket(2,2)  0            0;
%          0         0            kef_XZ(4,1) 0         kef_XZ(4,2)  0           0         0            kef_XZ(4,3) 0         kef_XZ(4,4)  0;
%          0         kef_XY(4,1)  0           0         0            kef_XY(4,2) 0         kef_XY(4,3)  0           0         0            kef_XY(4,4) ];
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
            
            % Compute torsion mass coefficients
            mel(anm.gln_tor,anm.gln_tor) = elem.torsionMassCoeff();
            
            % Compute flexural mass coefficients
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralMassCoeff_XY();
            mel(anm.gln_flx_XZ,anm.gln_flx_XZ) = elem.flexuralMassCoeff_XZ();
            
%  mel = [ mea(1,1)  0            0           0         0            0           mea(1,2)  0            0           0         0            0;
%          0         mef_XY(1,1)  0           0         0            mef_XY(1,2) 0         mef_XY(1,3)  0           0         0            mef_XY(1,4);
%          0         0            mef_XZ(1,1) 0         mef_XZ(1,2)  0           0         0            mef_XZ(1,3) 0         mef_XZ(1,4)  0;
%          0         0            0           met(1,1)  0            0           0         0            0           met(1,2)  0            0;
%          0         0            mef_XZ(2,1) 0         mef_XZ(2,2)  0           0         0            mef_XZ(2,3) 0         mef_XZ(2,4)  0;
%          0         mef_XY(2,1)  0           0         0            mef_XY(2,2) 0         mef_XY(2,3)  0           0         0            mef_XY(2,4);
%          mea(2,1)  0            0           0         0            0           mea(2,2)  0            0           0         0            0;
%          0         mef_XY(3,1)  0           0         0            mef_XY(3,2) 0         mef_XY(3,3)  0           0         0            mef_XY(3,4);
%          0         0            mef_XZ(3,1) 0         mef_XZ(3,2)  0           0         0            mef_XZ(3,3) 0         mef_XZ(3,4)  0;
%          0         0            0           met(2,1)  0            0           0         0            0           met(2,2)  0            0;
%          0         0            mef_XZ(4,1) 0         mef_XZ(4,2)  0           0         0            mef_XZ(4,3) 0         mef_XZ(4,4)  0;
%          0         mef_XY(4,1)  0           0         0            mef_XY(4,2) 0         mef_XY(4,3)  0           0         0            mef_XY(4,4) ];
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
            
            % Compute torsion lumped mass coefficients
            mel(anm.gln_tor,anm.gln_tor) = elem.torsionLumpedMassCoeff();
            
            % Compute flexural lumped mass coefficients
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.flexuralLumpedMassCoeff_XY();
            mel(anm.gln_flx_XZ,anm.gln_flx_XZ) = elem.flexuralLumpedMassCoeff_XZ();
        end
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint stiffness matrix in local system.
        % Output:
        %  kjl: target joint stiffness matrix in local system
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        function kjl = jointLocStiffMtx(~,srj)
            kjlx = [  srj.krx  -srj.krx;
                     -srj.krx   srj.krx ];
            kjly = [  srj.kry  -srj.kry;
                     -srj.kry   srj.kry ];
            kjlz = [  srj.krz  -srj.krz;
                     -srj.krz   srj.krz ];
            kjl = zeros(6);
            kjl([1 4],[1 4]) = kjlx;
            kjl([2 5],[2 5]) = kjly;
            kjl([3 6],[3 6]) = kjlz;
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
                        newF = model.nodes(n).T  * [model.nodes(n).load.static(1);
                                                    model.nodes(n).load.static(2);
                                                    model.nodes(n).load.static(3)];
                                                
                        newM = model.nodes(n).T  * [model.nodes(n).load.static(4);
                                                    model.nodes(n).load.static(5);
                                                    model.nodes(n).load.static(6)];
                    else
                        newF = [model.nodes(n).load.static(1);
                                model.nodes(n).load.static(2);
                                model.nodes(n).load.static(3)];
                            
                        newM = [model.nodes(n).load.static(4);
                                model.nodes(n).load.static(5);
                                model.nodes(n).load.static(6)];
                    end
                    
                    
                    
                    % Add applied force in global X direction
                    id = model.ID(1,n);
                    model.F(id) = model.F(id) + newF(1);
                    
                    % Add applied force in global Y direction
                    id = model.ID(2,n);
                    model.F(id) = model.F(id) + newF(2);
                    
                    % Add applied force in global Z direction
                    id = model.ID(3,n);
                    model.F(id) = model.F(id) + newF(3);
                    
                    % Add applied moment about global X direction
                    id = model.ID(4,n);
                    model.F(id) = model.F(id) + newM(1);
                    
                    % Add applied moment about global Y direction
                    id = model.ID(5,n);
                    model.F(id) = model.F(id) + newM(2);
                    
                    % Add applied moment about global Z direction
                    id = model.ID(6,n);
                    model.F(id) = model.F(id) + newM(3);
                    
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
                                                        f(3,:)];
                        else
                            newF = [f(1,:);...
                                    f(2,:);...
                                    f(3,:)];
                        end

                        % Add applied force in global X direction
                        id = model.ID(1,n);
                        model.F(id,:) = model.F(id,:) + newF(1,:);

                        % Add applied force in global Y direction
                        id = model.ID(2,n);
                        model.F(id,:) = model.F(id,:) + newF(2,:);

                        % Add applied force in global Z direction
                        id = model.ID(3,n);
                        model.F(id,:) = model.F(id,:) + newF(3,:);

                        % Add applied moment about global X direction
                        id = model.ID(4,n);
                        model.F(id,:) = model.F(id,:) + f(4,:);

                        % Add applied moment about global Y direction
                        id = model.ID(5,n);
                        model.F(id,:) = model.F(id,:) + f(5,:);

                        % Add applied moment about global Z direction
                        id = model.ID(6,n);
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
                id_3 = model.ID(3,n);
                id   = [id_1, id_2, id_3];
                
                % Rotate displacement to global system from local
                % inclined support system
                model.D(id) = model.nodes(n).T' * model.D(id);
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
                id_3 = model.ID(3,n);
                id   = [id_1, id_2, id_3];
                
                % Rotate displacement to global system from local
                % inclined support system
                v(id,:) = model.nodes(n).T' * v(id,:);
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
                id_3 = model.ID(3,n);
                id   = [id_1, id_2, id_3];
                
                % Rotate results to global system from local
                % inclined support system 
                for m = 1:size(model.results.dynamicDispl,3)
                    model.results.dynamicDispl(id,:,m)       = model.nodes(n).T' * model.results.dynamicDispl(id,:,m);
                    model.results.dynamicVeloc(id,:,m)       = model.nodes(n).T' * model.results.dynamicVeloc(id,:,m);
                    model.results.dynamicAccel(id,:,m)       = model.nodes(n).T' * model.results.dynamicAccel(id,:,m);
                    model.results.dynamicDisplForced(id,:,m) = model.nodes(n).T' * model.results.dynamicDisplForced(id,:,m);
                    model.results.dynamicVelocForced(id,:,m) = model.nodes(n).T' * model.results.dynamicVelocForced(id,:,m);
                    model.results.dynamicAccelForced(id,:,m) = model.nodes(n).T' * model.results.dynamicAccelForced(id,:,m);
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
            fel(anm.gln_flx_XZ) = load.flexuralDistribLoadFEF_XZ();
            
%             fel = [ fea(1);
%                     fef_XY(1);
%                     fef_XZ(1);
%                     0;
%                     fef_XZ(2);
%                     fef_XY(2);
%                     fea(2);
%                     fef_XY(3);
%                     fef_XZ(3);
%                     0;
%                     fef_XZ(4);
%                     fef_XY(4)];
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
            fel(anm.gln_flx_XZ) = load.flexuralThermalLoadFEF_XZ();
            
%             fel = [ fea(1);
%                     fef_XY(1);
%                     fef_XZ(1);
%                     0;
%                     fef_XZ(2);
%                     fef_XY(2);
%                     fea(2);
%                     fef_XY(3);
%                     fef_XZ(3);
%                     0;
%                     fef_XZ(4);
%                     fef_XY(4)];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  axial_force(nsteps+1,2)
        %   Ni = axial_force(:,1) - init value throughout time steps
        %   Nf = axial_force(:,2) - final value throughout time steps
        %  shear_force(nsteps+1,2)
        %   Qi = shear_force(:,1) - init value throughout time steps
        %   Qf = shear_force(:,2) - final value throughout time steps
        %  torsion moment(nsteps+1,2)
        %   Ti = torsion_moment(:,1) - init value throughout time steps
        %   Tf = torsion_moment(:,2) - final value throughout time steps
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
            elem.axial_force      = zeros(nsteps+1,2);
            elem.shear_force_Y    = zeros(nsteps+1,2);
            elem.shear_force_Z    = zeros(nsteps+1,2);
            elem.torsion_moment   = zeros(nsteps+1,2);
            elem.bending_moment_Y = zeros(nsteps+1,2);
            elem.bending_moment_Z = zeros(nsteps+1,2);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.axial_force(:,1)      = elem.axial_force(:,1)      + fel(1,:)';
            elem.axial_force(:,2)      = elem.axial_force(:,2)      + fel(7,:)';
            
            elem.shear_force_Y(:,1)    = elem.shear_force_Y(:,1)    + fel(2,:)';
            elem.shear_force_Y(:,2)    = elem.shear_force_Y(:,2)    + fel(8,:)';
            elem.shear_force_Z(:,1)    = elem.shear_force_Z(:,1)    + fel(3,:)';
            elem.shear_force_Z(:,2)    = elem.shear_force_Z(:,2)    + fel(9,:)';
            
            elem.torsion_moment(:,1)   = elem.torsion_moment(:,1)   + fel( 4,:)';
            elem.torsion_moment(:,2)   = elem.torsion_moment(:,2)   + fel(10,:)';
            
            elem.bending_moment_Y(:,1) = elem.bending_moment_Y(:,1) + fel( 5,:)';
            elem.bending_moment_Y(:,2) = elem.bending_moment_Y(:,2) + fel(11,:)';
            elem.bending_moment_Z(:,1) = elem.bending_moment_Z(:,1) + fel( 6,:)';
            elem.bending_moment_Z(:,2) = elem.bending_moment_Z(:,2) + fel(12,:)';
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dv (transversal displacement in local y-axis)
        %  intDispl(3,:) -> dw (transversal displacement in local z-axis)
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntDispl(~,elem)
            elem.intDispl = zeros(3,50);
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
            Nw = elem.flexuralDisplShapeFcnVector_XZ(x);
            
            % Initialize displacement shape function matrix
            N = zeros(3*np,elem.nen * anm.ndof);
            
            % Assemble displacement shape function matrix
            N(3*(1:np)-2,anm.gln_axl)    = Nu;
            N(3*(1:np)-1,anm.gln_flx_XY) = Nv;
            N(3*(1:np)  ,anm.gln_flx_XZ) = Nw;
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
            np = size(N,1)/3;
            
            % Compute internal displacements matrix
            del_aux = N * dl;
            
            % Predimension del as a 3D matrix (uncoupled dynamic analysis)
            del = zeros(3,np,size(dl,2));
            
            % Rearrange internal displacements matrix
            for i = 1:size(dl,2)
                del(:,:,i) = [ (del_aux(3*(1:np)-2,i))';
                               (del_aux(3*(1:np)-1,i))';
                               (del_aux(3*(1:np),i))'  ];
            end              
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: a 3x1 vector of with element internal displacements:
        %      du -> axial displacement
        %      dv -> transversal displacement in local y-axis direction
        %      dw -> transversal displacement in local z-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dv(x);
        %          dw(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Get number of points
            np = size(x,2);
            
            % Initialize displacements vector resulting from local analysis
            del = zeros(3,np);
            
            % Add the contribution of axial and transversal displacements
            % resulting from distributed loads
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(1,:) = elem.load.axialDistribLoadDispl(x);
                del(2,:) = elem.load.flexuralDistribLoadDispl_XY(x);
                del(3,:) = elem.load.flexuralDistribLoadDispl_XZ(x);
            end
            
            % Add the contribution of axial and transversal displacements
            % resulting from thermal loads
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Y ~= 0) || (elem.load.tempVar_Z ~= 0)
                del(1,:) = del(1,:) + elem.load.axialThermalLoadDispl(x);
                del(2,:) = del(2,:) + elem.load.flexuralThermalLoadDispl_XY(x);
                del(3,:) = del(3,:) + elem.load.flexuralThermalLoadDispl_XZ(x);
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
            [Q_XY,elem.maxShearForce_XY] = elem.intShearForce_XY(x);
            [Q_XZ,elem.maxShearForce_XZ] = elem.intShearForce_XZ(x);
            
            % Get bending moment values
            [M_XY,elem.maxBendMoment_XY] = elem.intBendingMoment_XY(x);
            [M_XZ,elem.maxBendMoment_XZ] = elem.intBendingMoment_XZ(x);
            
            % Assemble stress values matrix
            stressValues = [N;
                            Q_XY;
                            Q_XZ;
                            M_XY;
                            M_XZ];
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
            [Q_XY,elem.maxShearForce_XY] = elem.intDynamicShearForce_XY(x);
            [Q_XZ,elem.maxShearForce_XZ] = elem.intDynamicShearForce_XZ(x);
            
            % Get bending moment values
            [M_XY,elem.maxBendMoment_XY] = elem.intDynamicBendingMoment_XY(x);
            [M_XZ,elem.maxBendMoment_XZ] = elem.intDynamicBendingMoment_XZ(x);
            
            % Assemble stress values matrix         
            stressValues(:,:,1) = N;
            stressValues(:,:,2) = Q_XY;
            stressValues(:,:,3) = Q_XZ;
            stressValues(:,:,4) = M_XY;
            stressValues(:,:,5) = M_XZ;
            %stressValues(:,:,6) = ones(size(N,1),size(N,2)).* elem.torsion_moment(:,1);
            T = elem.torsion_moment(:,1);
            stressValues(:,:,6) = sparse(1:size(T,1),1:size(T,1),T') * ones(size(N,1),size(N,2));
        end
    end
end
