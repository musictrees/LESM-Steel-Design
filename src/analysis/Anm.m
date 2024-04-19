%% Anm (Analysis Model) Class
%
%% Description
%
% This is a handle super-class for the definition of an analysis model.
%
% An object of the *Anm* class "projects" the generic 3D behavior of
% objects of the <elem.html *Elem*> and <node.html *Node*> classes to a
% specific model behavior.
% Abstract methods of this class access specific properties related to each
% analysis model behavior.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <anm_truss2d.html 2D truss analysis model>.
% * <anm_frame2d.html 2D frame analysis model>.
% * <anm_truss3d.html 3D truss analysis model>.
% * <anm_frame3d.html 3D frame analysis model>.
% * <anm_grillage.html Grillage analysis model>.
%
classdef Anm < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        analysis_type =  0;   % flag for type of analysis model
        ndof          =  0;   % number of degrees-of-freedom per node
        nrdof         =  0;   % number of rotational degrees-of-freedom per node
        gln           = [];   % local nodal dofs gather vector
        gln_displ     = [];   % local nodal displacement dofs gather vector (inside gle)
        gln_rot       = [];   % local nodal rotational dofs gather vector (inside gle)
        gln_axl       = [];   % axial gather vector
        gln_flx_XY    = [];   % flexural (XY) gather vector
        gln_flx_XZ    = [];   % flexural (XZ) gather vector
        gln_tor       = [];   % torsional gather vector
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm(type,ndof,nrdof,gln,gln_rot)
            anm.analysis_type =     type;
            anm.ndof          =     ndof;
            anm.nrdof         =    nrdof;
            anm.gln           =      gln;
            anm.gln_rot       =  gln_rot;
            
            % Compute nodal displacement dof gather vector
            if ~isempty(gln_rot)
                aux = zeros(1,length(gln));
                aux(gln_rot) = gln(gln_rot);
                [~,anm.gln_displ,~] = find(gln-aux);
            else
                anm.gln_displ = anm.gln;
            end
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        rot = gblToLocElemRotMtx(~,elem)
        
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to the local inclined nodal support
        % system.
        % Output:
        %  rot: rotation transformation matrix
        %  flag: flag for there being inclined supports
        % Input arguments:
        %  elem: handle to an object of the Elem class
        [rot,flag] = gbltoLocNodeMtx(~,elem)
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint d.o.f. (degree of freedom) rotation
        % transformation matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        rot = gblToLocSrjointRotMtx(~,srj)
        
        %------------------------------------------------------------------
        % Assembles element stiffness matrix in local system.
        % Output:
        %  kel: target element stiffness matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        kel = elemLocStiffMtx(anm,elem)
        
        %------------------------------------------------------------------
        % Assembles element mass matrix in local system.
        % Output:
        %  mel: target element mass matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        mel = elemLocMassMtx(anm,elem)
        
        %------------------------------------------------------------------
        % Assembles element lumped mass matrix in local system.
        % Output:
        %  mel: target element mass matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        mel = elemLocLumpedMassMtx(anm,elem)
        
        %------------------------------------------------------------------
        % Assembles semi-rigid joint stiffness matrix in local system.
        % Output:
        %  kjl: target joint stiffness matrix in local system
        % Input arguments:
        %  srj: handle to an object of the Srjoint class
        kjl = jointLocStiffMtx(~,srj)
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        nodalLoads(~,model)
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing matrix.
        % Each column consists on a forcing vector on a time instant
        % Input arguments:
        %  model: handle to an object of the model class
        dynamicNodalLoads(~,model)
        
        %------------------------------------------------------------------
        % Compute nodal displacement obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        locNodeToGblDispl(~,model)
        
        %------------------------------------------------------------------
        % Compute vibration modes in global system,
        % rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        %  v: [OPTIONAL] vibration modes matrix
        % Output:
        %  v: rotated vibration modes matrix
        v = locNodeToGblVbrtnModes(~,model,v);
        
        %------------------------------------------------------------------
        % Compute nodal results obtained from solution in the global
        % system, rotating results from nodes with inclined supports.
        % Input:
        %  model: handle to an object of the model class
        locNodeToGblDynamicRes(~,model)
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied distributed load.
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        fel = elemLocDistribLoadFEF(anm,load)
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied thermal load (temperature variation).
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        fel = elemLocThermalLoadFEF(anm,load)
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        initIntForce(~,elem)
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        assembleIntForce(~,elem,fel)
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        initIntDispl(~,elem)
        
        %------------------------------------------------------------------
        % Assembles displacement shape functions matrix evaluated at a 
        % given cross-section position.
        % Output:
        %  N: displacement shape function matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position
        N = displShapeFcnMtx(anm,elem,x)
        
        %------------------------------------------------------------------
        % Computes internal displacements matrix from global analysis
        % Output:
        %  del: element internal displacements matrix
        % Input arguments:
        %  N: displacement shape function matrix
        %  dl: nodal displacements, on local coordinates
        del = gblAnlIntDisplMtx(~,N,dl)
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: element internal displacements vector in a given
        %       cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        del = lclAnlIntDispl(anm,elem,x)
        
        %------------------------------------------------------------------
        % Computes element internal stresses vector in local system,
        % in a given cross-section position.
        % Output:
        %  stressValues: element internal stresses vector in a given
        %                cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        stressValues = intStress(anm,elem,x)
        
        %------------------------------------------------------------------
        % Computes element internal stresses vector in local system,
        % in a given cross-section position, throughout time steps.
        % Output:
        %  stressValues: element internal stresses vector in a given
        %                cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        stressValues = intDynamicStress(anm,elem,x)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Initializes global d.o.f (degree of freedom) numbering ID matrix
        % with ones and zeros, and counts total number of equations of free
        % d.o.f.'s, total number of equations of fixed d.o.f.'s and total
		% number of d.o.f's constrained by spring support.
        %  ID matrix initialization:
        %  if ID(k,n) =  0, d.o.f. k of node n is free.
        %  if ID(k,n) =  1, d.o.f. k of node n is constrained by support.
		%  if ID(k,n) =  2, d.o.f. k of node n is constrained by spring.
        %  if ID(k,n) = -1, d.o.f. k of node n is constrained by ficticious support.
        % Input arguments:
        %  model: handle to an object of the model class
        function setupDOFNum(anm,model)
            include_constants;
            
            % Dimension global d.o.f. numbering ID matrix
            model.ID = zeros(anm.ndof,model.nnp);
            
            % Initialize number of fixed d.o.f. and d.o.f. constrained by spring
            model.neqfixed = 0;
            model.neqspring = 0;

            % Count number of fixed d.o.f. and setup ID matrix
            for n = 1:model.nnp
                i = 0;
                for dof = anm.gln
                    i = i + 1;
                    % Check for fixed dof
                    if (model.nodes(n).ebc(dof) == FIXED_DOF || model.nodes(n).ebc(dof) == FICTFIXED_DOF)
                        model.neqfixed = model.neqfixed + 1;
                        model.ID(i,n) = FIXED_DOF;
                    % Check for dof associated to spring
                    elseif model.nodes(n).ebc(dof) == SPRING_DOF
                        model.neqspring = model.neqspring + 1;
                        model.ID(i,n) = SPRING_DOF;
                    end
                end
            end
                
            % Compute total number of free d.o.f.
            model.neqfree = model.neq - model.neqfixed - model.neqspring;
        end
        
        %------------------------------------------------------------------
        % Adds prescribed displacements (known support settlement values)
        % to global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of the global displacement vector that corresponds to a free d.o.f.,
		% or a d.o.f. constrained by spring.
        % Input arguments:
        %  model: handle to an object of the model class
        function setupPrescDispl(anm,model)
            for n = 1:model.nnp
                if ~isempty(model.nodes(n).prescDispl)
                    i = 0;
                    for dof = anm.gln
                        i = i + 1;
                        % Add prescribed displacement
                        id = model.ID(i,n);
                        if (id > (model.neqfree + model.neqspring)) && (model.nodes(n).prescDispl(dof) ~= 0)
                            model.D(id) = model.nodes(n).prescDispl(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adds spring stiffness coefficients to global stiffness matrix.
        % Avoids storing a spring stiffness component in a position of the
        % global stiffness matrix that corresponds to a free d.o.f.
        % Input arguments:
        %  model: handle to an object of the model class
        function setupSpringStiff(anm,model)
            model.gblSprStiff = zeros(1,model.neqspring);
            for n = 1:model.nnp
                if ~isempty(model.nodes(n).springStiff)
                    i = 0;
                    for dof = anm.gln
                        i = i + 1;
                        % Add spring stiffness to global stiffness matrix
                        id = model.ID(i,n);
                        if (id > model.neqfree) && (id <= model.neqfree + model.neqspring) && (model.nodes(n).springStiff(dof) ~= 0)
                            model.K(id,id) = model.K(id,id) + model.nodes(n).springStiff(dof);
                            model.gblSprStiff(id-model.neqfree) = model.nodes(n).springStiff(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adds concentrated mass coefficients to global mass matrix.
        % Input arguments:
        %  model: handle to an object of the model class
        function setupConcentratedMass(anm,model)
            for n = 1:model.nnp
                if model.nodes(n).displMass > 0
                    for i = 1:length(anm.gln_displ)
                        % Add concentrated mass to global mass matrix
                        id = model.ID(i,n);
                        model.M(id,id) = model.M(id,id) + model.nodes(n).displMass;
                    end
                end
                if model.nodes(n).rotMass > 0
                    for i = (1:length(anm.gln_rot))+length(anm.gln_displ)
                        % Add concentrated mass to global mass matrix
                        id = model.ID(i,n);
                        model.M(id,id) = model.M(id,id) + model.nodes(n).rotMass;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element gather vector (gle) that stores element d.o.f.'s
        % equation numbers.
        function assembleGle(anm,model,elem)
            include_constants;
            
            % Initialize element gather vector
            elem.gle = zeros(elem.nen * anm.ndof, 1);
            
            for n = 1:elem.nen
                % Assemble translation d.o.f.'s of initial node to element gather vector
                for dof = 1:(anm.ndof - anm.nrdof)
                    id = anm.gln_displ(dof);
                    elem.gle(id+anm.ndof*(n-1)) = model.ID(id,elem.nodes(n).id);
                end
                
                % Assemble rotation d.o.f.'s of initial node to element gather vector
                for dof = 1:anm.nrdof
                    id = anm.gln_rot(dof);
                    if n == 1
                        if elem.hingei == SEMIRIGID_END
                            elem.gle(id+anm.ndof*(n-1)) = elem.srjointi.eqs(dof);
                        else
                            elem.gle(id+anm.ndof*(n-1)) = model.ID(id,elem.nodes(n).id);
                        end
                    else % if n == 2
                        if elem.hingef == SEMIRIGID_END
                            elem.gle(id+anm.ndof*(n-1)) = elem.srjointf.eqs(dof);
                        else
                            elem.gle(id+anm.ndof*(n-1)) = model.ID(id,elem.nodes(n).id);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles semi-rigid gather vector (gle) that stores joint
        % d.o.f.'s equation numbers.
        function assembleGlj(anm,model,srjoint)
            % Initialize joint gather vector
            srjoint.glj = zeros(2 * anm.nrdof, 1);
            
            for i = 1:anm.nrdof
                % Get local rotaional dof id
                id = anm.gln_rot(i);
                
                % First equation number corresponds to end node rotation
                srjoint.glj(i) = model.ID(id,srjoint.node.id);
                
                % Second equation number corresponds to the additional
                % rotation equation number of semi-rigid joint
                srjoint.glj(i+anm.nrdof) = srjoint.eqs(i);
            end
        end
        
        %------------------------------------------------------------------
        % Sets nodal initial conditions to global initial conditions matrix
        % Input arguments:
        %  model: handle to an object of the model class
        function initialConditions(anm,model)
            for n = 1:model.nnp
                i = 0;
                for dof = anm.gln
                    i = i + 1;
                    % Set initial conditions
                    id = model.ID(i,n);
                    model.c0(id,:) = model.nodes(n).initCond(dof,:);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of an Anm object.
        function clean(anm)
            anm.analysis_type =  0;
            anm.ndof          =  0;
            anm.nrdof         =  0;
            anm.gln           = [];
            anm.gln_displ     = [];
            anm.gln_rot       = [];
            anm.gln_axl       = [];
            anm.gln_flx_XY    = [];
            anm.gln_flx_XZ    = [];
            anm.gln_tor       = [];
        end
    end
end
