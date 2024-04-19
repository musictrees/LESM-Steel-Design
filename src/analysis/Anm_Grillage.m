%% Anm_Grillage Class
%
%% Description
%
% This is a sub-class of the <anm.html *Anm*> class for the
% implementation of the *Grillage* analysis model.
%
% A grillage model has the following assumptions:
%
% * It is a 2D model, which in LESM is considered in the XY-plane.
% * Beam elements are laid out in a grid pattern in a single plane,
%   rigidly connected at nodes. However, a grillage element might have
%   a hinge (rotation liberation) at an end or hinges at both ends.
% * It is assumed that a hinge in a grillage element releases continuity of
%   both bending and torsion rotations.
% * By assumption, there is only out-of-plane behavior, which
%   includes displacement transversal to the grillage plane, and
%   rotations about in-plane axes.
% * Internal forces at any cross-section of a grillage element are:
%   shear force (local z direction), bending moment (about local y direction),
%   and torsion moment (about local x direction).
%   By assumption, there is no axial force in a grillage element.
% * Each node of a grillage model has three d.o.f.'s: a transversal
%   displacement in Z direction, and rotations about the X and Y directions.
%
classdef Anm_Grillage < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Grillage()
            include_constants;
            anm = anm@Anm(GRILLAGE_ANALYSIS,3,2,[4,5,3],1:2);

            anm.gln_flx_XZ = [3,2,6,5];
            anm.gln_tor    = [1,4];
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
        function [rot,flag] = gbltoLocNodeMtx(~,~)
            % There are no inclined displacement supports on grillage models
            rot = 1;
            flag = false;
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
            % rot = [ T11 T12  0   0
            %         T21 T22  0   0
            %         0   0    T11 T12
            %         0   0    T21 T22 ]
            rot = blkdiag(T([1 2],[1 2]),T([1 2],[1 2]));
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
            
            % Compute torsion stiffness coefficients
            kel(anm.gln_tor,anm.gln_tor) = elem.torsionStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kel(anm.gln_flx_XZ,anm.gln_flx_XZ) = elem.flexuralStiffCoeff_XZ();

%             kel = [ ket(1,1)  0         0         ket(1,2)  0         0;
%                     0         kef(2,2)  kef(2,1)  0         kef(2,4)  kef(2,3);
%                     0         kef(1,2)  kef(1,1)  0         kef(1,4)  kef(1,3);
%                     ket(2,1)  0         0         ket(2,2)  0         0;
%                     0         kef(4,2)  kef(4,1)  0         kef(4,4)  kef(4,3);
%                     0         kef(3,2)  kef(3,1)  0         kef(3,4)  kef(3,3) ];
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
            
            % Compute torsion mass coefficients
            mel(anm.gln_tor,anm.gln_tor) = elem.torsionMassCoeff();
            
            % Compute flexural mass coefficients
            mel(anm.gln_flx_XZ,anm.gln_flx_XZ) = elem.flexuralMassCoeff_XZ();

%             mel = [ met(1,1)  0         0         met(1,2)  0         0;
%                     0         mef(2,2)  mef(2,1)  0         mef(2,4)  mef(2,3);
%                     0         mef(1,2)  mef(1,1)  0         mef(1,4)  mef(1,3);
%                     met(2,1)  0         0         met(2,2)  0         0;
%                     0         mef(4,2)  mef(4,1)  0         mef(4,4)  mef(4,3);
%                     0         mef(3,2)  mef(3,1)  0         mef(3,4)  mef(3,3) ];
        end
        
        %------------------------------------------------------------------
        % Assembles element mass matrix in local system.
        % Output:
        %  mel: target element mass matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function mel = elemLocLumpedMassMtx(anm,elem)
            % Initialize element lumped mass matrix in local system
            mel = zeros(elem.nen * anm.ndof);
            
            % Compute torsion lumped mass coefficients
            mel(anm.gln_tor,anm.gln_tor) = elem.torsionLumpedMassCoeff();
            
            % Compute flexural lumped mass coefficients
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
            kjl = zeros(4);
            kjl([1 3],[1 3]) = kjlx;
            kjl([2 4],[2 4]) = kjly;
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        function nodalLoads(~,model)
            for n = 1:model.nnp
                if ~isempty(model.nodes(n).load.static)
                    
                    % Add applied moment about global X direction
                    id = model.ID(1,n);
                    model.F(id) = model.F(id) + model.nodes(n).load.static(4);
                    
                    % Add applied moment about global Y direction
                    id = model.ID(2,n);
                    model.F(id) = model.F(id) + model.nodes(n).load.static(5);
                    
                    % Add applied force in global Z direction
                    id = model.ID(3,n);
                    model.F(id) = model.F(id) + model.nodes(n).load.static(3);
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

                        % Add applied moment about global X direction
                        id = model.ID(1,n);
                        model.F(id,:) = model.F(id,:) + f(4,:);

                        % Add applied moment about global Y direction
                        id = model.ID(2,n);
                        model.F(id,:) = model.F(id,:) + f(5,:);

                        % Add applied force in global Z direction
                        id = model.ID(3,n);
                        model.F(id,:) = model.F(id,:) + f(3,:);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute nodal displacement obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % DOES NOTHING ON GRILLAGE MODELS
        function locNodeToGblDispl(~,~)
        end
        
        %------------------------------------------------------------------
        % Compute vibration modes in global system,
        % rotating results from nodes with inclined supports.
        % DOES NOTHING ON GRILLAGE MODELS
        function v = locNodeToGblVbrtnModes(~,model,v)
            if nargin < 3
              v = [model.V; zeros(model.neqfixed,model.n_modes)];
            end
        end
        
        %------------------------------------------------------------------
        % Compute nodal results obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % DOES NOTHING ON GRILLAGE MODELS
        function locNodeToGblDynamicRes(~,~)
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
            
            % Compute flexural (transversal) fixed end force components
            fel(anm.gln_flx_XZ) = load.flexuralDistribLoadFEF_XZ();

%             fel = [ 0;
%                     fef(2);
%                     fef(1);
%                     0;
%                     fef(4);
%                     fef(3) ];
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
            
            % Compute flexural (transversal) fixed end force components
            fel(anm.gln_flx_XZ) = load.flexuralThermalLoadFEF_XZ();

%             fel = [ 0;
%                     fef(2);
%                     fef(1);
%                     0;
%                     fef(4);
%                     fef(3) ];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  torsion moment(nsteps+1,2)
        %   Ti = torsion_moment(:,1) - init value throughout time steps
        %   Tf = torsion_moment(:,2) - final value throughout time steps
        %  bending_moment(nsteps+1,2)
        %   Mi = bending_moment(:,1) - init value throughout time steps
        %   Mf = bending_moment(:,2) - final value throughout time steps
        %  shear_force(nsteps+1,2)
        %   Qi = shear_force(:,1) - init value throughout time steps
        %   Qf = shear_force(:,2) - final value throughout time steps
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  nsteps: number of steps for transient analysis
        function initIntForce(~,elem,nsteps)
            if (nargin < 3)
                nsteps = 0;
            end
            elem.torsion_moment   = zeros(nsteps+1,2);
            elem.bending_moment_Y = zeros(nsteps+1,2);
            elem.shear_force_Z    = zeros(nsteps+1,2);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.torsion_moment(:,1)   = elem.torsion_moment(:,1)   + fel(1,:)';
            elem.torsion_moment(:,2)   = elem.torsion_moment(:,2)   + fel(4,:)';
            
            elem.bending_moment_Y(:,1) = elem.bending_moment_Y(:,1) + fel(2,:)';
            elem.bending_moment_Y(:,2) = elem.bending_moment_Y(:,2) + fel(5,:)';
            
            elem.shear_force_Z(:,1)    = elem.shear_force_Z(:,1)    + fel(3,:)';
            elem.shear_force_Z(:,2)    = elem.shear_force_Z(:,2)    + fel(6,:)';
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dw (transversal displacement in local z-axis)
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
            
            % Compute transversal displacement shape functions vector
            Nw = elem.flexuralDisplShapeFcnVector_XZ(x);
            
            % Initialize displacement shape function matrix
            N = zeros(2*np,elem.nen * anm.ndof);
            
            % Assemble displacement shape function matrix
            N(2*(1:np), anm.gln_flx_XZ) = Nw;
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
        %      du -> axial displacement (always zero)
        %      dw -> transversal displacement in local z-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dw(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Get number of points
            np = size(x,2);
            
            % Initialize displacements vector resulting from local analysis
            del = zeros(2,np);
            
            % Add the contribution of transversal displacements resulting
            % from distributed loads (there is no axial displacement in a
            % grillage model)
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(2,:) = elem.load.flexuralDistribLoadDispl_XZ(x);
            end
            
            % Add the contribution of transversal displacements resulting
            % from thermal loads (there is no axial displacement in a
            % grillage model)
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Z ~= 0)
                del(2,:) = del(2,:) + elem.load.flexuralThermalLoadDispl_XZ(x);
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
            % Get shear force values
            [Q,elem.maxShearForce_XZ] = elem.intShearForce_XZ(x);
            
            % Get bending moment values
            [M,elem.maxBendMoment_XZ] = elem.intBendingMoment_XZ(x);
            
            % Assemble stress values matrix
            stressValues = [Q;
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
            % Get shear force values
            [Q,elem.maxShearForce_XZ] = elem.intDynamicShearForce_XZ(x);
            
            % Get bending moment values
            [M,elem.maxBendMoment_XZ] = elem.intDynamicBendingMoment_XZ(x);
            
            % Assemble stress values matrix
            stressValues(:,:,1) = Q;
            stressValues(:,:,2) = M;
            %stressValues(:,:,3) = ones(size(Q,1),size(Q,2)).* elem.torsion_moment(:,1);
            T = elem.torsion_moment(:,1);
            stressValues(:,:,3) = sparse(1:size(T,1),1:size(T,1),T') * ones(size(Q,1),size(Q,2));
        end
    end
end
