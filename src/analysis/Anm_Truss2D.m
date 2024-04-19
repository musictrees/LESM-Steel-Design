%% Anm_Truss2D Class
%
%% Description
%
% This is a sub-class of the <anm.html *Anm*> class for the
% implementation of the *2D Truss* analysis model.
%
% A 2D truss model has the following assumptions:
%
% * Truss elements are bars connected at their ends only, and they are
%   connected by friction-less pins. Therefore, a truss element does not
%   present any secondary bending moment or torsion moment induced by
%   rotation continuity at joints.
% * A truss model is loaded only at nodes.
%   Any load action along an element, such as self weight, is statically
%   transferred as concentrated forces to the element end nodes.
% * Local bending of elements due to internal loads is neglected, when
%   compared to the effect of global load acting on the truss.
% * Therefore, there is only one type of internal force in a truss
%   element: axial internal force, which may be tension or compression.
% * A 2D truss model is considered to be laid in the XY-plane, with only
%   in-plane behavior, that is, there is no displacement transversal to
%   the truss plane.
% * Each node of a 2D truss model has two d.o.f.'s (degrees of
%   freedom): a horizontal displacement in X direction and a vertical
%   displacement in Y direction.
%
classdef Anm_Truss2D < Anm
    %% Private attributes
    properties (SetAccess = private, GetAccess = private)
        gln_flx2axl = [1,3]; % Gather vector of flexural to axial dofs (no rotational dofs)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Truss2D()
            include_constants;
            anm = anm@Anm(TRUSS2D_ANALYSIS,2,0,1:2,[]);
            
            anm.gln_axl    = [1,3];
            anm.gln_flx_XY = [2,4]; % used for distributed loads FEF and
        end                         % displacement functions matrix and 
                                    % local mass matrix
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
            % Get 2x2 basis rotation transformation matrix
            T = elem.T(1:2,1:2);
            
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
            rot = blkdiag(T1,T2);
        end
        
        %------------------------------------------------------------------
        function gblToLocSrjointRotMtx(~,~)
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
            
%             kel = [ kea(1,1)  0     kea(1,2)  0;
%                     0         0     0         0;
%                     kea(2,1)  0     kea(2,2)  0;
%                     0         0     0         0 ];
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
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.axialMassCoeff();
            
%             mel = [ mea(1,1)  0         mea(1,2)         0;
%                     0         mea(1,1)     0         mea(1,2);
%                     mea(2,1)  0         mea(2,2)         0;
%                     0         mea(1,2)     0         mea(2,2)];
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
            mel(anm.gln_flx_XY,anm.gln_flx_XY) = elem.axialLumpedMassCoeff();
        end
        
        %------------------------------------------------------------------
        function jointLocStiffMtx(~,~)
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
            fef = load.flexuralDistribLoadFEF_XY();
            fel(anm.gln_flx_XY) = fef(anm.gln_flx2axl);
            
%             fel = [ fea(1);
%                     fef(1);
%                     fea(2);
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
            
            % Compute axial fixed end force components
            fel(anm.gln_axl) = load.axialThermalLoadFEF();
            
            % Compute flexural (transversal) fixed end force components
            fef = load.flexuralThermalLoadFEF_XY();
            fel(anm.gln_flx_XY) = fef(anm.gln_flx2axl);

%             fel = [ fea(1);
%                     fef(1);
%                     fea(2);
%                     fef(3) ];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  axial_force(nsteps+1,2)
        %   Ni = axial_force(:,1) - init value throughout time steps
        %   Nf = axial_force(:,2) - final value throughout time steps
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  nsteps: number of steps for transient analysis
        function initIntForce(~,elem,nsteps)
            if (nargin < 3)
                nsteps = 0;
            end
            elem.axial_force = zeros(nsteps+1,2);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.axial_force(:,1) = elem.axial_force(:,1) + fel(1,:)';
            elem.axial_force(:,2) = elem.axial_force(:,2) + fel(3,:)';
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
            N(2*(1:np)  , anm.gln_flx_XY) = Nv(:,anm.gln_flx2axl);
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
        function del = lclAnlIntDispl(~,~,x)
            % Get number of points
            np = size(x,2);
            
            % Internal displacements from local analysis are neglected in
            % truss models.
            del = zeros(2,np);
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
            [N,~] = elem.intAxialForce(x);
            
            % Assemble stress values matrix
            stressValues = N;
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
            [N,~] = elem.intDynamicAxialForce(x);
            
            % Assemble stress values matrix
            stressValues = N;
        end
    end
end
