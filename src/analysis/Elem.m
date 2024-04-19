%% Elem (Element) Class
%
%% Description
%
% This is a handle class for the definition of an element.
%
% This class generically handles a three-dimensional behavior of a linear
% element. An object of the <anm.html *Anm*> class is responsible to
% "project" this generic 3D behavior to a specific model behavior, such as
% 2D frame, 2D truss, grillage, 3D truss or 3D frame model.
%
% Two different types of flexural behavior are considered:
%
% In Euler-Bernoulli flexural behavior, it is assumed that there is no
% shear deformation. As consequence, bending of a linear element is such
% that its cross-section remains plane and normal to the element
% longitudinal axis.
%
% In Timoshenko flexural behavior, shear deformation is considered in an
% approximated manner. Bending of a linear structure element is such that
% its cross-section remains plane but it is not normal to the element
% longitudinal axis.
%
% In truss models, the two types of elements may be used indistinguishably,
% since there is no bending behavior of a truss element, and
% Euler-Bernoulli elements and Timoshenko elements are equivalent for the
% axial behavior.
%
% To handle different types of flexural behavior, the formulation is based
% on the Timoshenko´s theory and it is expressed in a generic way, in which
% the shear parameter Omega assumes null value for the Navier´s theory.
%
%% Element Local Coordinate System
%
% In 2D models of LESM, the local axes of an element are defined
% uniquely in the following manner:
%
% * The local z-axis of an element is always in the direction of the
%   global Z-axis, which is perpendicular to the model plane and its
%   positive direction points out of the screen.
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The local y-axis of an element lays on the global XY-plane and is
%   perpendicular to the element x-axis in such a way that the
%   cross-product (x-axis * y-axis) results in a vector in the
%   global Z direction.
%
% In 3D models of LESM, the local y-axis and z-axis are defined by an
% auxiliary vector vz = (vzx,vzy,vzz), which is an element property and
% should be specified as an input data of each element:
%
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The auxiliary vector lays in the local xz-plane of an element, and the
%   cross-product (vz * x-axis) defines the the local y-axis vector.
% * The direction of the local z-axis is then calculated with the
%   cross-product (x-axis * y-axis).
%
% In 2D models, the auxiliary vector vz is automatically set to (0,0,1).
% In 3D models, it is important that the auxiliary vector is not parallel
% to the local x-axis; otherwise, the cross-product (vz * x-axis) is zero
% and local y-axis and z-axis will not be defined.
%
classdef Elem < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nen                = 2;    % number of nodes (always equal to 2)
        type               = 0;    % flag for type of element (Navier = 0, Timoshenko = 1)
        anm                = [];   % handle to an object of the Anm class
        load               = [];   % handle to an object of the Lelem class
        material           = [];   % handle to an object of the Material class
        section            = [];   % handle to an object of the Section class
        nodes              = [];   % vector of handles to objects of the Node class [initial_node final_node]
        hingei             = 0;    % flag for rotation liberation at initial node (hinged = 0, continuous = 1, semirigid = 2)
        hingef             = 0;    % flag for rotation liberation at final node   (hinged = 0, continuous = 1, semirigid = 2)
        srjointi           = [];   % handle to Srjoint object at initial node
        srjointf           = [];   % handle to Srjoint object at final node
        kri                = [];   % initial rotational stiffness vector [krx kry krz]
        krf                = [];   % final rotational stiffness vector [krx kry krz]
        length             = 0;    % element length
        cosine_X           = 0;    % cosine of orientation angle with global X-axis
        cosine_Y           = 0;    % cosine of orientation angle with global Y-axis
        cosine_Z           = 0;    % cosine of orientation angle with global Z-axis
        vz                 = [];   % auxiliary vector in local xz-plane [vz_X vz_Y vz_Z]
        gle                = [];   % gather vector (stores element d.o.f. eqn. numbers)
        T                  = [];   % basis rotation transformation matrix between two coordinate systems
        rot                = [];   % rotation transformation matrix from global system to local system
        rot_inclSupp       = [];   % rotation transformation matrix from global system to inclined support local system
        kel                = [];   % stiffness matrix in local system
        mel                = [];   % mass matrix in local system
        mass_type          = 1;    % flag for mass matrx type (0 = lumped, 1 = consistent, 2 = mixed)
        mass_consideration = 1;    % flag for mass consideration (0 = element mass is not considered, 1 = mass is considered)
        mass_mi            = 1;    % mixed mass matrix proportion coefficient
        fel_distribLoad    = [];   % fixed end forces in local system for a distributed load
        fel_thermalLoad    = [];   % fixed end forces in local system for a thermal load
        axial_force        = [];   % array of axial internal forces at element ends
        shear_force_Y      = [];   % array of shear internal forces at element ends in local y-axis
        shear_force_Z      = [];   % array of shear internal forces at element ends in local z-axis
        bending_moment_Y   = [];   % array of bending internal moments at element ends in local y-axis
        bending_moment_Z   = [];   % array of bending internal moments at element ends in local z-axis
        torsion_moment     = [];   % array of torsional internal moments at element ends
        intDispl           = [];   % array of axial and transversal internal displacements
        natVibration       = [];   % 3D matrix of normalized vibration modes behavior [dof, intCoords, modes]
        dynamicIntDispl    = [];   % 3D matrix of internal dynamic displacement [dof, intCoords, time steps]
        intCoords          = [];   % array of coordinates for 50 points inside element
        intStresses        = [];   % array of stress values on 50 cross-sections
        intForcesEnvelop   = [];   % 3D matrix of internal forces envelop [{max,min}, intCoords, internal force]
        maxAxialForce      = [];   % array of maximum axial force values and its local coordinates
        maxShearForce_XY   = [];   % array of maximum shear force values on the xy-plane and its local coordinates
        maxShearForce_XZ   = [];   % array of maximum shear force values on the xz-plane and its local coordinates
        maxBendMoment_XY   = [];   % array of maximum bending moment values about the xy-plane and its local coordinates
        maxBendMoment_XZ   = [];   % array of maximum bending moment values about the xz-plane and its local coordinates
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem(type,anm,mat,sec,nodes,hi,hf,vz,kri,krf,mc,gle,af,sfy,sfz,bmy,bmz,tm,intdispl)
            if (nargin > 0)
                % Get nodal coordinates
                xi = nodes(1).coord(1);
                yi = nodes(1).coord(2);
                zi = nodes(1).coord(3);
                xf = nodes(2).coord(1);
                yf = nodes(2).coord(2);
                zf = nodes(2).coord(3);
                
                % Calculate element length
                dx = xf - xi;
                dy = yf - yi;
                dz = zf - zi;
                L  = sqrt(dx^2 + dy^2 + dz^2);
                
                % Calculate element cosines with global axes
                cx = dx/L;
                cy = dy/L;
                cz = dz/L;
                
                % Set properties
                elem.type = type;
                elem.anm = anm;
                elem.material = mat;
                elem.section = sec;
                elem.nodes = nodes;
                elem.hingei = hi;
                elem.hingef = hf;
                elem.length = L;
                elem.cosine_X = cx;
                elem.cosine_Y = cy;
                elem.cosine_Z = cz;
                elem.vz = vz;
                
                if (nargin > 8)
                elem.kri = kri;
                elem.krf = krf;
                end      
                                
                if (nargin > 10)
                elem.mass_consideration = mc;
                end
                
                if (nargin > 11)
                elem.gle = gle;
                end
                
                if (nargin > 12)
                elem.axial_force = af;
                elem.shear_force_Y = sfy;
                elem.shear_force_Z = sfz;
                elem.bending_moment_Y = bmy;
                elem.bending_moment_Z = bmz;
                elem.torsion_moment = tm;
                elem.intDispl = intdispl;
                end
                
                elem.intStresses = [];
                elem.maxAxialForce = [];
                elem.maxShearForce_XY = [];
                elem.maxShearForce_XZ = [];
                elem.maxBendMoment_XY = [];
                elem.maxBendMoment_XZ = [];
                
                % Calculate deformed configuration coordinates of 50 cross-
                % sections along element local axis X.
                intCoords = ones(3,50);
                i = linspace(0,L,size(intCoords,2));
                intCoords(1,:) = intCoords(1,:) * xi + i * cx;
                intCoords(2,:) = intCoords(2,:) * yi + i * cy;
                intCoords(3,:) = intCoords(3,:) * zi + i * cz;
                elem.intCoords = intCoords;
                
                % Compute and set basis rotation transformation matrix and
                % d.o.f. rotation transformation matrix
                elem.T = elem.rotTransMtx();
                elem.rot = anm.gblToLocElemRotMtx(elem);
                
                % Create and set load object
                elem.load = Lelem(elem);
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes elements local axis direction vectors.
        % Output:
        %  x: local axis x direction vector
        %  y: local axis y direction vector
        %  z: local axis z direction vector
        function [x,y,z] = locAxis(elem)
            % Get nodal coordinates
            xi = elem.nodes(1).coord(1);
            yi = elem.nodes(1).coord(2);
            zi = elem.nodes(1).coord(3);
            xf = elem.nodes(2).coord(1);
            yf = elem.nodes(2).coord(2);
            zf = elem.nodes(2).coord(3);
            
            % Calculate element local x-axis
            x = [xf-xi, yf-yi, zf-zi];
            x = x / norm(x);
            
            % Calculate element local y-axis
            y = cross(elem.vz,x);
            y = y / norm(y);
            
            % Calculate element local z-axis
            z = cross(x,y);
        end
        
        %------------------------------------------------------------------
        % Computes basis rotation transformation matrix between local and
        % global coordinate systems of an element.
        % Output:
        %  T: basis rotation transformation matrix
        function T = rotTransMtx(elem)
            % Get element local axis
            [x,y,z] = locAxis(elem);
            
            % Global axes
            X = [1,0,0];
            Y = [0,1,0];
            Z = [0,0,1];
            
            % Calculate angle cosines between local x-axis and global axes
            cxX = dot(x,X);
            cxY = dot(x,Y);
            cxZ = dot(x,Z);
            
            % Calculate angle cosines between local y-axis and global axes
            cyX = dot(y,X);
            cyY = dot(y,Y);
            cyZ = dot(y,Z);
            
            % Calculate angle cosines between local z-axis and global axes
            czX = dot(z,X);
            czY = dot(z,Y);
            czZ = dot(z,Z);
            
            % Assemble basis transformation matrix
            T = [ cxX cxY cxZ;
                  cyX cyY cyZ;
                  czX czY czZ ];
        end
        
        %------------------------------------------------------------------
        % Computes element stiffness matrix in global system.
        % Output:
        %  keg: element stiffness matrix in global system
        function keg = gblStiffMtx(elem)
            % Compute and store element stiffness matrix in local system
            elem.kel = elem.anm.elemLocStiffMtx(elem);
            
            % Transform element stiffness matrix from local to global system
            keg = elem.rot' * elem.kel * elem.rot;
            
            % Get inclined supports rotation matrix
            [elem.rot_inclSupp,thereAreInclSupps] = elem.anm.gbltoLocNodeMtx(elem);
            
            % Transform element stiffness matrix from global to inclined support local system
            % obs.: This affects only displacement DOFs.
            if thereAreInclSupps
                keg  = elem.rot_inclSupp * keg * elem.rot_inclSupp';
            end
        end
        
        %------------------------------------------------------------------
        % Computes element mass matrix in global system.
        % Output:
        %  meg: element mass matrix in global system
        function meg = gblMassMtx(elem)
            include_constants;
            
            % Switch between lumped and consistent mass matrices
            switch elem.mass_type
                case LUMPED_MASS
                    % Compute and store element lumped mass matrix in local system
                    elem.mel = elem.anm.elemLocLumpedMassMtx(elem);
                case CONSISTENT_MASS
                    % Compute and store element consistent mass matrix in local system
                    elem.mel = elem.anm.elemLocMassMtx(elem);
                case MIXED_MASS
                    % Compute element lumped mass matrix in local system
                    mel_lumped = elem.anm.elemLocLumpedMassMtx(elem);
                    
                    % Compute element consistent mass matrix in local system
                    mel_consis = elem.anm.elemLocMassMtx(elem);
                    
                    % Compute and store element mixed mass matrix in local system
                    elem.mel = (1 - elem.mass_mi) * mel_lumped +...
                                    elem.mass_mi  * mel_consis;
            end
            
            % Transform element mass matrix from local to global system
            meg = elem.rot' * elem.mel * elem.rot;
            
            % Get inclined supports rotation matrix
            [elem.rot_inclSupp,thereAreInclSupps] = elem.anm.gbltoLocNodeMtx(elem);
            
            % Transform element mass matrix from global to inclined support local system
            % obs.: This affects only displacement DOFs.
            if thereAreInclSupps
                meg  = elem.rot_inclSupp * meg * elem.rot_inclSupp';
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal force vector in local system at end nodes,
        % for the global analysis (from nodal displacements and rotations).
        % Output:
        %  fel: element internal force vector in local system at end nodes
        % Input arguments:
        %  model: handle to an object of the model class
        function fel = gblAnlIntForce(elem,model)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = model.D(elem.gle);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Compute internal force vector in local system
            fel = elem.kel * dl;
        end
        
        %------------------------------------------------------------------
        % Computes element internal force array in local system at end nodes,
        % for the global analysis (from nodal displacements and rotations).
        % Output:
        %  fel: element internal force vector in local system at end nodes
        %       throughout time steps of analysis
        % Input arguments:
        %  model: handle to an object of the model class
        function fel = gblDynamicAnlIntForce(elem,model)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = sum(model.results.dynamicDispl(elem.gle,:,:),3);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Compute internal force vector in local system
            fel = elem.kel * dl;
        end
        
        %------------------------------------------------------------------
        % Computes element axial and transversal internal displacements
        % vector in local system at a given cross-section position,
        % for the global analysis (from nodal displacements and rotations).
        % Output:
        %  del: element internal displacements vector at given
        %       cross-section position
        % Input arguments:
        %  model: handle to an object of the model class
        %  x: cross-section position on element local x-axis
        function del = gblAnlIntDispl(elem,model,x)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = model.D(elem.gle);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Evaluate element shape function matrix at given cross-section
            N = elem.anm.displShapeFcnMtx(elem,x);
            
            % Compute internal displacements matrix from global analysis
            del = elem.anm.gblAnlIntDisplMtx(N,dl);
        end
        
        %------------------------------------------------------------------
        % Input arguments:
        %  v: eigenvetors matrix
        %  x: cross-section position on element local x-axis
        function vibrtnModes(elem,v,x)
            % Get nodal displacements and rotations at element end nodes
            % in global system (from normalized eigenvectors)
            dg = v(elem.gle,:);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Evaluate element shape function matrix at given cross-section
            N = elem.anm.displShapeFcnMtx(elem,x);
            
            % Compute internal displacements matrix from global analysis
            elem.natVibration = elem.anm.gblAnlIntDisplMtx(N,dl);
        end
        
        %------------------------------------------------------------------
        % Input arguments:
        %  v: (displacement x time) matrix
        %  x: cross-section position on element local x-axis
        function dynamicDispl(elem,d,x)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = d(elem.gle,:);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Evaluate element shape function matrix at given cross-section
            N = elem.anm.displShapeFcnMtx(elem,x);
            
            % Compute internal displacements matrix from global analysis
            elem.dynamicIntDispl = elem.anm.gblAnlIntDisplMtx(N,dl);
        end
        
        %------------------------------------------------------------------
        % Generates element axial stiffness coefficient matrix.
        % Output:
        % kea: a 2x2 matrix with axial stiffness coefficients
        function kea = axialStiffCoeff(elem)
            E = elem.material.elasticity;
            A = elem.section.area_x;
            L = elem.length;
            
            kea = [ E*A/L  -E*A/L;
                   -E*A/L   E*A/L ];
        end
        
        %------------------------------------------------------------------
        % Generates element torsion stiffness coefficient matrix.
        % Output:
        %  ket: a 2x2 matrix with torsion stiffness coefficients
        function ket = torsionStiffCoeff(elem)
            include_constants;
            
            G  = elem.material.shear;
            Jt = elem.section.inertia_x;
            L  = elem.length;
            
            if ((elem.hingei == CONTINUOUS_END || elem.hingei == SEMIRIGID_END)...
               && (elem.hingef == CONTINUOUS_END || elem.hingef == SEMIRIGID_END))
                
                ket = [ G*Jt/L  -G*Jt/L;
                       -G*Jt/L   G*Jt/L ];
                
            else
                
                ket = [ 0  0;
                        0  0 ];
                
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural stiffness coefficient matrix in local
        % xy-plane.
        % Output:
        %  kef: a 4x4 matrix with flexural stiffness coefficients:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   kef = [ k_Dini_Dini  k_Dini_Rini  k_Dini_Dend  k_Dini_Rend;
        %           k_Rini_Dini  k_Rini_Rini  k_Rini_Dend  k_Dini_Rend;
        %           k_Dend_Dini  k_Dend_Rini  k_Dend_Dend  k_Dend_Rend;
        %           k_Rend_Dini  k_Rend_Rini  k_Rend_Dend  k_Rend_Rend ]
        function kef = flexuralStiffCoeff_XY(elem)
            include_constants;
            
            % Basic element properties
            E  = elem.material.elasticity;
            I  = elem.section.inertia_z;
            L  = elem.length;
            EI  = E * I;
            
            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                G  = elem.material.shear;
                As = elem.section.area_y;
                
                Omega  = EI / (G * As * L * L);
            end
            
            % Auxiliary parameters
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            laL  = lambda * L;
            laL2 = laL    * L;
            laL3 = laL2   * L;
            
            muL  = mu   * L;
            muL2 = muL  * L;
            muL3 = muL2 * L;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                
                kef = [ 12*EI/muL3   6*EI/muL2        -12*EI/muL3   6*EI/muL2;
                         6*EI/muL2   4*lambda*EI/muL   -6*EI/muL2   2*gamma*EI/muL;
                       -12*EI/muL3  -6*EI/muL2         12*EI/muL3  -6*EI/muL2;
                         6*EI/muL2   2*gamma*EI/muL    -6*EI/muL2   4*lambda*EI/muL ];
                
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                
                kef = [  3*EI/laL3   0                 -3*EI/laL3   3*EI/laL2;
                         0           0                  0           0;
                        -3*EI/laL3   0                  3*EI/laL3  -3*EI/laL2;
                         3*EI/laL2   0                 -3*EI/laL2   3*EI/laL ];
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                
                kef = [  3*EI/laL3   3*EI/laL2         -3*EI/laL3   0;
                         3*EI/laL2   3*EI/laL          -3*EI/laL2   0;
                        -3*EI/laL3  -3*EI/laL2          3*EI/laL3   0;
                         0           0                  0           0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                
                % kef = [  0           0                  0           0;
                %          0           0                  0           0;
                %          0           0                  0           0;
                %          0           0                  0           0 ];
                kef  = sparse(4,4);
                
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural stiffness coefficient matrix in local
        % xz-plane.
        % Output:
        %  kef: a 4x4 matrix with flexural stiffness coefficients:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   kef = [ k_Dini_Dini  k_Dini_Rini  k_Dini_Dend  k_Dini_Rend;
        %           k_Rini_Dini  k_Rini_Rini  k_Rini_Dend  k_Dini_Rend;
        %           k_Dend_Dini  k_Dend_Rini  k_Dend_Dend  k_Dend_Rend;
        %           k_Rend_Dini  k_Rend_Rini  k_Rend_Dend  k_Rend_Rend ]
        function kef = flexuralStiffCoeff_XZ(elem)
            include_constants;
            
            % Basic element properties
            E  = elem.material.elasticity;
            I  = elem.section.inertia_y;
            L  = elem.length;
            EI  = E * I;

            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                G  = elem.material.shear;
                As = elem.section.area_z;
                
                Omega  = EI / (G * As * L * L);
            end
            
            % Auxiliary parameters
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            laL  = lambda * L;
            laL2 = laL    * L;
            laL3 = laL2   * L;
            
            muL  = mu   * L;
            muL2 = muL  * L;
            muL3 = muL2 * L;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                
                kef = [ 12*EI/muL3  -6*EI/muL2        -12*EI/muL3  -6*EI/muL2;
                        -6*EI/muL2   4*lambda*EI/muL    6*EI/muL2   2*gamma*EI/muL;
                       -12*EI/muL3   6*EI/muL2         12*EI/muL3   6*EI/muL2;
                        -6*EI/muL2   2*gamma*EI/muL     6*EI/muL2   4*lambda*EI/muL ];
                
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                
                kef = [  3*EI/laL3   0                 -3*EI/laL3  -3*EI/laL2;
                         0           0                  0           0;
                        -3*EI/laL3   0                  3*EI/laL3   3*EI/laL2;
                        -3*EI/laL2   0                  3*EI/laL2   3*EI/laL ];
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                
                kef = [  3*EI/laL3  -3*EI/laL2         -3*EI/laL3   0;
                        -3*EI/laL2   3*EI/laL           3*EI/laL2   0;
                        -3*EI/laL3   3*EI/laL2          3*EI/laL3   0;
                         0           0                  0           0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                
                % kef = [  0           0                  0           0;
                %          0           0                  0           0;
                %          0           0                  0           0;
                %          0           0                  0           0 ];
                kef  = sparse(4,4);
                
            end
        end
        
        %------------------------------------------------------------------
        % Generates element axial mass coefficient matrix.
        % Output:
        % mea: a 2x2 matrix with axial mass coefficients
        function mea = axialMassCoeff(elem)
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            mea = [ M/3   M/6;
                    M/6   M/3 ];
        end
        
        %------------------------------------------------------------------
        % Generates element axial lumped mass coefficient matrix.
        % Output:
        % mea: a 2x2 matrix with axial mass coefficients
        function mea = axialLumpedMassCoeff(elem)
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end

            mea = (M/2) * [ 1 0 ;
                            0 1 ];
        end
        
        %------------------------------------------------------------------
        % Generates element torsion mass coefficient matrix.
        % Output:
        %  met: a 2x2 matrix with torsion mass coefficients
        function met = torsionMassCoeff(elem)
            include_constants;
            
            rho  = elem.material.density;
            Jt   = elem.section.inertia_x;
            A    = elem.section.area_x;
            L    = elem.length;
            M    = rho * A * L;
            At   = Jt / A;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            if ((elem.hingei == CONTINUOUS_END || elem.hingei == SEMIRIGID_END)...
               && (elem.hingef == CONTINUOUS_END || elem.hingef == SEMIRIGID_END))
                
                met = (M * At / 6) *[ 2  1;
                                      1  2 ];
                
            else
                
                met = [ 0  0;
                        0  0 ];
                
            end
        end
        
        %------------------------------------------------------------------
        % Generates element torsion lumped mass coefficient matrix.
        % Output:
        %  met: a 2x2 matrix with torsion mass coefficients
        function met = torsionLumpedMassCoeff(~)
            met = zeros(2);
        end
        
        %------------------------------------------------------------------
        % Generates element flexural mass coefficient matrix in local
        % xy-plane.
        % Output:
        %  mef: a 4x4 matrix with flexural mass coefficients:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   mef = [ m_Dini_Dini  m_Dini_Rini  m_Dini_Dend  m_Dini_Rend;
        %           m_Rini_Dini  m_Rini_Rini  m_Rini_Dend  m_Dini_Rend;
        %           m_Dend_Dini  m_Dend_Rini  m_Dend_Dend  m_Dend_Rend;
        %           m_Rend_Dini  m_Rend_Rini  m_Rend_Dend  m_Rend_Rend ]
        function mef = flexuralMassCoeff_XY(elem)
            include_constants;
            
            % Basic element properties
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            L2   = L * L;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                E  = elem.material.elasticity;
                I  = elem.section.inertia_z;
                G  = elem.material.shear;
                As = elem.section.area_y;
                
                Omega  = E * I / (G * As * L * L);
            end
            
            % Auxiliary parameters
            mu     = 1 + 12 * Omega;
            theta  = 1 + 4  * Omega;
            lambda = 1 + 3  * Omega;
            
            om2 = Omega * Omega;
            mu2 = mu * mu;
            la2 = lambda * lambda;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
               
                m_1 =  264 + 48  * theta / mu - 48 * Omega / mu2;
                m_2 = (44  - 108 * Omega / mu - 24 * Omega / mu2) * L;
                m_3 =  108 + 360 * Omega / mu + 72 * Omega * theta / mu2;
                m_4 = (26  + 108 * Omega / mu + 24 * Omega / mu2) * L;
                m_5 = (7 + 1 / mu2) * L2;
                m_6 = (7 - 1 / mu2) * L2;
               
                mef = (M/840) * [ m_1   m_2   m_3  -m_4 ;
                                  m_2   m_5   m_4  -m_6 ;
                                  m_3   m_4   m_1  -m_2 ;
                                 -m_4  -m_6  -m_2   m_5 ];
           
%                mef = (M/210) * [ 78       11*L     27      -13*L/2;
%                                  11*L     2*L2     13*L/2  -3*L2/2;
%                                  27       13*L/2   78      -11*L;
%                                 -13*L/2  -3*L2/2  -11*L     2*L2 ];
                 
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
               
                m_1 =  198 + 198 * Omega / lambda + 48  * om2 / la2;
                m_2 =  117 + 117 * Omega / lambda - 144 * om2 / la2;
                m_3 =  408 - 432 * Omega / lambda + 144 * om2 / la2;
                m_4 = (33 / lambda + 48 * Omega / la2) * L;
                m_5 = (72 / lambda - 48 * Omega / la2) * L;
                m_6 = (16 / la2) * L2;
                
                mef = (M/840) * [ m_1   0   m_2  -m_4 ;
                                    0   0     0     0 ;
                                  m_2   0   m_3  -m_5 ;
                                 -m_4   0  -m_5   m_6 ];
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                
                m_1 =  198 + 198 * Omega / lambda + 48  * om2 / la2;
                m_2 =  117 + 117 * Omega / lambda - 144 * om2 / la2;
                m_3 =  408 - 432 * Omega / lambda + 144 * om2 / la2;
                m_4 = (33 / lambda + 48 * Omega / la2) * L;
                m_5 = (72 / lambda - 48 * Omega / la2) * L;
                m_6 = (16 / la2) * L2;
                
                mef = (M/840) * [ m_3   m_5   m_2   0 ;
                                  m_5   m_6   m_4   0 ;
                                  m_2   m_4   m_1   0 ;
                                    0     0     0   0 ];
                                 
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                mef = [ M/3  0  M/6  0;
                          0  0    0  0;
                        M/6  0  M/3  0;
                          0  0    0  0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural lumped mass coefficient mtx in local
        % xy-plane.
        % Output:
        %  mef: a 4x4 matrix with flexural mass coefficients:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   mef = [ m_Dini_Dini  m_Dini_Rini  m_Dini_Dend  m_Dini_Rend;
        %           m_Rini_Dini  m_Rini_Rini  m_Rini_Dend  m_Dini_Rend;
        %           m_Dend_Dini  m_Dend_Rini  m_Dend_Dend  m_Dend_Rend;
        %           m_Rend_Dini  m_Rend_Rini  m_Rend_Dend  m_Rend_Rend ]
        function mef = flexuralLumpedMassCoeff_XY(elem)
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            mef = (M/2) * [ 1 0 0 0 ;
                            0 0 0 0 ;
                            0 0 1 0 ;
                            0 0 0 0 ];
        end
        
        %------------------------------------------------------------------
        % Generates element flexural mass coefficient matrix in local
        % xz-plane.
        % Output:
        %  mef: a 4x4 matrix with flexural mass coefficients:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   mef = [ m_Dini_Dini  m_Dini_Rini  m_Dini_Dend  m_Dini_Rend;
        %           m_Rini_Dini  m_Rini_Rini  m_Rini_Dend  m_Dini_Rend;
        %           m_Dend_Dini  m_Dend_Rini  m_Dend_Dend  m_Dend_Rend;
        %           m_Rend_Dini  m_Rend_Rini  m_Rend_Dend  m_Rend_Rend ]
        function mef = flexuralMassCoeff_XZ(elem)
            include_constants;
            
            % Basic element properties
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            L2   = L * L;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                E  = elem.material.elasticity;
                I  = elem.section.inertia_y;
                G  = elem.material.shear;
                As = elem.section.area_z;
                
                Omega  = E * I / (G * As * L * L);
            end
            
            % Auxiliary parameters
            mu     = 1 + 12 * Omega;
            theta  = 1 + 4  * Omega;
            lambda = 1 + 3  * Omega;
            
            om2 = Omega * Omega;
            mu2 = mu * mu;
            la2 = lambda * lambda;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                
                m_1 =  264 + 48  * theta / mu - 48 * Omega / mu2;
                m_2 = (44  - 108 * Omega / mu - 24 * Omega / mu2) * L;
                m_3 =  108 + 360 * Omega / mu + 72 * Omega * theta / mu2;
                m_4 = (26  + 108 * Omega / mu + 24 * Omega / mu2) * L;
                m_5 = (7 + 1 / mu2) * L2;
                m_6 = (7 - 1 / mu2) * L2;
               
                mef = (M/840) * [ m_1  -m_2   m_3   m_4 ;
                                 -m_2   m_5  -m_4  -m_6 ;
                                  m_3  -m_4   m_1   m_2 ;
                                  m_4  -m_6   m_2   m_5 ];

%                mef = (M/210) * [ 78      -11*L     27       13*L/2;
%                                 -11*L     2*L2    -13*L/2  -3*L2/2;
%                                  27      -13*L/2   78       11*L;
%                                  13*L/2  -3*L2/2   11*L     2*L2 ];

              
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
               
                m_1 =  198 + 198 * Omega / lambda + 48  * om2 / la2;
                m_2 =  117 + 117 * Omega / lambda - 144 * om2 / la2;
                m_3 =  408 - 432 * Omega / lambda + 144 * om2 / la2;
                m_4 = (33 / lambda + 48 * Omega / la2) * L;
                m_5 = (72 / lambda - 48 * Omega / la2) * L;
                m_6 = (16 / la2) * L2;
                
                mef = (M/840) * [ m_1   0   m_2   m_4 ;
                                    0   0     0     0 ;
                                  m_2   0   m_3   m_5 ;
                                  m_4   0   m_5   m_6 ];
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                
                m_1 =  198 + 198 * Omega / lambda + 48  * om2 / la2;
                m_2 =  117 + 117 * Omega / lambda - 144 * om2 / la2;
                m_3 =  408 - 432 * Omega / lambda + 144 * om2 / la2;
                m_4 = (33 / lambda + 48 * Omega / la2) * L;
                m_5 = (72 / lambda - 48 * Omega / la2) * L;
                m_6 = (16 / la2) * L2;
                
                mef = (M/840) * [ m_3  -m_5   m_2   0 ;
                                 -m_5   m_6  -m_4   0 ;
                                  m_2  -m_4   m_1   0 ;
                                    0     0     0   0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                mef = [ M/3  0  M/6  0;
                          0  0    0  0;
                        M/6  0  M/3  0;
                          0  0    0  0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural lumped mass coefficient mtx in local
        % xz-plane.
        % Output:
        %  mef: a 4x4 matrix with flexural mass coefficients:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   mef = [ m_Dini_Dini  m_Dini_Rini  m_Dini_Dend  m_Dini_Rend;
        %           m_Rini_Dini  m_Rini_Rini  m_Rini_Dend  m_Dini_Rend;
        %           m_Dend_Dini  m_Dend_Rini  m_Dend_Dend  m_Dend_Rend;
        %           m_Rend_Dini  m_Rend_Rini  m_Rend_Dend  m_Rend_Rend ]
        function mef = flexuralLumpedMassCoeff_XZ(elem)
            rho  = elem.material.density;
            A    = elem.section.area_x;
            L    = elem.length;
            M    = rho * A * L;
            
            if elem.mass_consideration == 0
                M =0; 
            end
            
            mef = (M/2) * [ 1 0 0 0 ;
                            0 0 0 0 ;
                            0 0 1 0 ;
                            0 0 0 0 ];
        end
        
        %------------------------------------------------------------------
        % Evaluates element axial displacement shape functions at a given
        % cross-section position of given element.
        % Axial shape functions give the longitudinal displacement of an
        % element subjected to unit nodal displacements in axial directions.
        % Output:
        %  Nu: a 1x2 vector with axial displacement shape functions:
        %      D -> axial displacement in local x-axis
        %      ini -> initial node
        %      end -> end node
        %   Nu = [ N_Dini  N_Dend ]
        % Input arguments:
        %  x: cross-section position on element local x-axis
        function Nu = axialDisplShapeFcnVector(elem,x)
            L = elem.length;
            
            Nu1 = 1 - x./L;
            Nu2 = x./L;
            
            Nu = [ Nu1'  Nu2' ];
        end
        
        %------------------------------------------------------------------
        % Evaluates element flexural displacement shape functions in local
        % xy-plane at a given cross-section position.
        % Flexural shape functions give the transversal displacement of an
        % element subjected to unit nodal displacements and rotations.
        % Output:
        %  Nv: a 1x4 vector with flexural displacement shape functions:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   Nv = [ N_Dini  N_Rini  N_Dend  N_Rend ]
        % Input arguments:
        %  x: element cross-section position in local x-axis
        function Nv = flexuralDisplShapeFcnVector_XY(elem,x)
            include_constants;
            
            % Get number of points
            np = size(x,2);
            
            % Basic element properties
            L  = elem.length;
            L2  = L^2;
            L3  = L2 * L;         
            
            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                E  = elem.material.elasticity;
                G  = elem.material.shear;
                As = elem.section.area_y;
                I  = elem.section.inertia_z;
                
                Omega  = E * I / (G * As * L2);
            end
            
            % Auxiliary parameters
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                Nv1 = 1 - 12*Omega*x./(L*mu) - 3*x.^2/(L2*mu) + 2*x.^3/(L3*mu);
                Nv2 = x - 6*Omega*x./mu - 2*lambda*x.^2/(L*mu) + x.^3/(L2*mu);
                Nv3 = 12*Omega*x./(L*mu) + 3*x.^2/(L2*mu) - 2*x.^3/(L3*mu);
                Nv4 = -6*Omega*x./mu - gamma*x.^2/(L*mu) + x.^3/(L2*mu);
                
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                Nv1 = 1 - 3*x./(2*L*lambda) - 3*Omega*x./(L*lambda) + x.^3/(2*L3*lambda);
                Nv2 = zeros(1,np);
                Nv3 = 3*x./(2*L*lambda) + 3*Omega*x./(L*lambda) - x.^3/(2*L3*lambda);
                Nv4 = -gamma*x./(2*lambda) - 3*Omega*x./lambda + x.^3/(2*L2*lambda);
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                Nv1 = 1 - 3*Omega*x./(L*lambda) - 3*x.^2/(2*L2*lambda) + x.^3/(2*L3*lambda);
                Nv2 = x - 3*Omega*x./lambda - 3*x.^2/(2*L*lambda) + x.^3/(2*L2*lambda);
                Nv3 = 3*Omega*x./(L*lambda) + 3*x.^2/(2*L2*lambda) - x.^3/(2*L3*lambda);
                Nv4 = zeros(1,np);
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                Nv1 = 1 - x./L;
                Nv2 = zeros(1,np);
                Nv3 = x./L;
                Nv4 = zeros(1,np);
            end
            
            Nv = [ Nv1'  Nv2'  Nv3'  Nv4' ];
        end
        
        %------------------------------------------------------------------
        % Evaluates element flexural displacement shape functions in local
        % xz-plane at a given cross-section position.
        % Flexural shape functions give the transversal displacement of an
        % element subjected to unit nodal displacements and rotations.
        % Output:
        %  Nw: a 1x4 vector with flexural displacement shape functions:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   Nw = [ N_Dini  N_Rini  N_Dend  N_Rend ]
        % Input arguments:
        %  x: element cross-section position in local x-axis
        function Nw = flexuralDisplShapeFcnVector_XZ(elem,x)
            include_constants;
            
            % Get number of points
            np = size(x,2);
            
            % Basic element properties
            L  = elem.length;
            L2  = L^2;
            L3  = L2 * L;
            
            % Timoshenko parameter
            if elem.type == 0     % Navier element
                Omega  = 0;
            elseif elem.type == 1 % Timoshenko element
                E  = elem.material.elasticity;
                G  = elem.material.shear;
                As = elem.section.area_z;
                I  = elem.section.inertia_y;
                
                Omega  = E * I / (G * As * L2);
            end
            
            % Auxiliary parameters
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            if ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ...
               ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                Nw1 = 1 - 12*Omega*x./(L*mu) - 3*x.^2/(L2*mu) + 2*x.^3/(L3*mu);
                Nw2 = -x + 6*Omega*x./mu + 2*lambda*x.^2/(L*mu) - x.^3/(L2*mu);
                Nw3 = 12*Omega*x./(L*mu) + 3*x.^2/(L2*mu) - 2*x.^3/(L3*mu);
                Nw4 = 6*Omega*x./mu + gamma*x.^2/(L*mu) - x.^3/(L2*mu);
                
            elseif (elem.hingei == HINGED_END) && ...
                   ((elem.hingef == CONTINUOUS_END) || (elem.hingef == SEMIRIGID_END))
                Nw1 = 1 - 3*x./(2*L*lambda) - 3*Omega*x./(L*lambda) + x.^3/(2*L3*lambda);
                Nw2 = zeros(1,np);
                Nw3 = 3*x./(2*L*lambda) + 3*Omega*x./(L*lambda) - x.^3/(2*L3*lambda);
                Nw4 = gamma*x./(2*lambda) + 3*Omega*x./lambda - x.^3/(2*L2*lambda);
                
            elseif ((elem.hingei == CONTINUOUS_END) || (elem.hingei == SEMIRIGID_END)) && ... 
                   (elem.hingef == HINGED_END)
                Nw1 = 1 - 3*Omega*x./(L*lambda) - 3*x.^2/(2*L2*lambda) + x.^3/(2*L3*lambda);
                Nw2 = -x + 3*Omega*x./lambda + 3*x.^2/(2*L*lambda) - x.^3/(2*L2*lambda);
                Nw3 = 3*Omega*x./(L*lambda) + 3*x.^2/(2*L2*lambda) - x.^3/(2*L3*lambda);
                Nw4 = zeros(1,np);
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                Nw1 = 1 - x./L;
                Nw2 = zeros(1,np);
                Nw3 = x./L;
                Nw4 = zeros(1,np);
            end
            
            Nw = [ Nw1'  Nw2'  Nw3'  Nw4' ];
        end
        
        %------------------------------------------------------------------
        function [N,Nmax] = intAxialForce(elem,x)
            include_constants;
            
            % Get element length
            L = elem.length;
            
            % Get element internal axial force value at both ends
            N1 = -elem.axial_force(1);
            N2 = elem.axial_force(2);
            
            % Avoid numeric garbage
            if abs(N1) < 1e-12
                N1 = 0;
            end
            if abs(N2) < 1e-12
                N2 = 0;
            end

            if ((isempty(elem.load.uniformLcl) == 0) || ...
               (isempty(elem.load.linearLcl) == 0)) &&...
               elem.anm.analysis_type ~= TRUSS2D_ANALYSIS &&...
               elem.anm.analysis_type ~= TRUSS3D_ANALYSIS
                qui = 0;
                quf = 0;
                qli = 0;
                qlf = 0;
                
                if isempty(elem.load.uniformLcl) == 0
                    qui = elem.load.uniformLcl(1);
                    quf = qui;
                end
                
                if isempty(elem.load.linearLcl) == 0
                    qli = elem.load.linearLcl(1);
                    qlf = elem.load.linearLcl(4);
                end
                
                % Calculate total load value at both ends
                qi = qui + qli;
                qf = quf + qlf;
                
                % Load equation coefficients:
                % p(x) = Ax + B
                A = (qf - qi)/L;
                B = qi;
                
                if A ~= 0
                    % Axial force equation coefficients:
                    % N(x) = Ax^2 + Bx + C
                    A = -A/2;
                    B = -B;
                    C = N1;
                    
                    % Calculate force intensity on current
                    % cross-section using shear force equation
                    N = A * x.^2 + B * x + C;
                    
                    % Get maximum coord
                    xNmax = roots([2*A B]);
                    
                    % Get maximum axial force
                    if xNmax > 0 && xNmax < L
                        Nmx = A * xNmax^2 + B * xNmax + C;
                        Nmax = [Nmx xNmax];
                    else
                        Nmax = [];
                    end
                else
                    N = linspace(N1,N2,size(x,2));
                    Nmax = [];
                end
            elseif elem.anm.analysis_type == TRUSS2D_ANALYSIS ||...
                   elem.anm.analysis_type == TRUSS3D_ANALYSIS
                % On truss analysis, axial force is constant
                N = N1;
                Nmax = [];
            else
                % Axial force is constant
                N = ones(1,size(x,2)) * N1;
                Nmax = [];
            end
        end
        
        %------------------------------------------------------------------
        function [N,Nmax] = intDynamicAxialForce(elem,x)
            % Get element internal axial force value at initial end
            N1 = -elem.axial_force(:,1);
            
            % Obs.: For now, only dynamic nodal loads are being considered,
            % thus axial force is constant along elements.
            
            % Avoid numeric garbage
            if all(abs(N1) < 1e-12)
                N1 = zeros(size(N1,1),1);
            end

            % Axial force is constant for nodal loads
            %N = ones(size(N1,1),size(x,2)).* N1;
            N = sparse(1:size(N1,1),1:size(N1,1),N1') * ones(size(N1,1),size(x,2));
            Nmax = [];
        end
        
        %------------------------------------------------------------------
        function [Q,Qmax] = intShearForce_XY(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal shear force value at both ends
            Q1 = elem.shear_force_Y(1);
            Q2 = -elem.shear_force_Y(2);
            
            % Avoid numeric garbage
            if abs(Q1) < 1e-10
                Q1 = 0;
            end
            if abs(Q2) < 1e-10
                Q2 = 0;
            end
            
            if (isempty(elem.load.uniformLcl) == 0) || ...
               (isempty(elem.load.linearLcl) == 0)
                qui = 0;
                quf = 0;
                qli = 0;
                qlf = 0;
                
                if isempty(elem.load.uniformLcl) == 0
                    qui = elem.load.uniformLcl(2);
                    quf = qui;
                end
                
                if isempty(elem.load.linearLcl) == 0
                    qli = elem.load.linearLcl(2);
                    qlf = elem.load.linearLcl(5);
                end
                
                % Calculate total load value at both ends
                qi = qui + qli;
                qf = quf + qlf;
                
                % Load equation coefficients:
                % q(x) = Ax + B
                A = (qf - qi)/L;
                B = qi;
                
                if A ~= 0
                    % Shear force equation coefficients:
                    % Q(x) = Ax^2 + Bx + C
                    A = A/2;
                    C = Q1;
                    
                    % Calculate force intensity on current
                    % cross-section using shear force equation
                    Q = A * x.^2 + B * x + C;
                    
                    % Get maximum coord
                    xQmax = roots([2*A B]);
                    
                    % Get maximum shear force
                    if xQmax > 0 && xQmax < L
                        Qmx = A * xQmax^2 + B * xQmax + C;
                        Qmax = [Qmx xQmax];
                    else
                        Qmax = [];
                    end
                else
                    Q = linspace(Q1,Q2,size(x,2));
                    Qmax = [];
                end
                
            else
                Q = ones(1,size(x,2)) * Q1;
                Qmax = [];
            end
        end
        
        %------------------------------------------------------------------
        function [Q,Qmax] = intDynamicShearForce_XY(elem,x)
            % Get element internal shear force value at both ends
            Q1 = elem.shear_force_Y(:,1);
            
            % Obs.: For now, only dynamic nodal loads are being considered,
            % thus shear force is constant along elements.
            
            % Avoid numeric garbage
            if all(abs(Q1) < 1e-12)
                Q1 = zeros(size(Q1,1),1);
            end
            
            %Q = ones(size(Q1,1),size(x,2)).* Q1;
            Q = sparse(1:size(Q1,1),1:size(Q1,1),Q1') * ones(size(Q1,1),size(x,2));
            Qmax = [];
        end
        
        %------------------------------------------------------------------
        function [Q,Qmax] = intShearForce_XZ(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal shear force value at both ends
            Q1 = elem.shear_force_Z(1);
            Q2 = -elem.shear_force_Z(2);
            
            % Avoid numeric garbage
            if abs(Q1) < 1e-10
                Q1 = 0;
            end
            if abs(Q2) < 1e-10
                Q2 = 0;
            end
            
            if (isempty(elem.load.uniformLcl) == 0) || ...
               (isempty(elem.load.linearLcl) == 0)
                qui = 0;
                quf = 0;
                qli = 0;
                qlf = 0;
                
                if isempty(elem.load.uniformLcl) == 0
                    qui = elem.load.uniformLcl(3);
                    quf = qui;
                end
                
                if isempty(elem.load.linearLcl) == 0
                    qli = elem.load.linearLcl(3);
                    qlf = elem.load.linearLcl(6);
                end
                
                % Calculate total load value at both ends
                qi = qui + qli;
                qf = quf + qlf;
                
                % Load equation coefficients:
                % q(x) = Ax + B
                A = (qf - qi)/L;
                B = qi;
                
                if A ~= 0
                    % Shear force equation coefficients:
                    % Q(x) = Ax^2 + Bx + C
                    A = A/2;
                    C = Q1;
                    
                    % Calculate force intensity on current
                    % cross-section using shear force equation
                    Q = A * x.^2 + B * x + C;
                    
                    % Get maximum coord
                    xQmax = roots([2*A B]);
                    
                    % Get maximum shear force
                    if xQmax > 0 && xQmax < L
                        Qmx = A * xQmax^2 + B * xQmax + C;
                        Qmax = [Qmx xQmax];
                    else
                        Qmax = [];
                    end
                else
                    Q = linspace(Q1,Q2,size(x,2));
                    Qmax = [];
                end
                
            else
                Q = ones(1,size(x,2)) * Q1;
                Qmax = [];
            end
        end
        
        %------------------------------------------------------------------
        function [Q,Qmax] = intDynamicShearForce_XZ(elem,x)
            % Get element internal shear force value at both ends
            Q1 = elem.shear_force_Z(:,1);
            
            % Obs.: For now, only dynamic nodal loads are being considered,
            % thus shear force is constant along elements.
            
            % Avoid numeric garbage
            if all(abs(Q1) < 1e-12)
                Q1 = zeros(size(Q1,1),1);
            end
            
            %Q = ones(size(Q1,1),size(x,2)).* Q1;
            Q = sparse(1:size(Q1,1),1:size(Q1,1),Q1') * ones(size(Q1,1),size(x,2));
            Qmax = [];
        end
        
        %------------------------------------------------------------------
        function [M,Mmax] = intBendingMoment_XY(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal forces value at both ends
            M1 = -elem.bending_moment_Z(1);
            M2 = elem.bending_moment_Z(2);
            Q1 = elem.shear_force_Y(1);
            
            % Avoid numeric garbage
            if abs(M1) < 1e-10
                M1 = 0;
            end
            if abs(M2) < 1e-10
                M2 = 0;
            end
            
            if (isempty(elem.load.uniformLcl) == 0) || ...
               (isempty(elem.load.linearLcl) == 0)
                qui = 0;
                quf = 0;
                qli = 0;
                qlf = 0;
                
                if isempty(elem.load.uniformLcl) == 0
                    qui = elem.load.uniformLcl(2);
                    quf = qui;
                end
                
                if isempty(elem.load.linearLcl) == 0
                    qli = elem.load.linearLcl(2);
                    qlf = elem.load.linearLcl(5);
                end
                
                % Calculate total load value at both ends
                qi = qui + qli;
                qf = quf + qlf;
                
                % Load equation coefficients:
                % q(x) = Ax + B
                A = (qf - qi)/L;
                B = qi;
                
                % Bending moment equation coefficients:
                % M(x) = Ax^3 + Bx^2 + Cx + D
                A = A/6;
                B = B/2;
                C = Q1;
                D = M1;
                
                % Calculate internal force intensity on current
                % cross-section using bending moment equation
                M = A * x.^3 + B * x.^2 + C * x + D;
                
                % Get maximum coord
                xMmax = roots([3*A 2*B C]);
                
                % Initialize maximum bending moment matrix
                Mmax = zeros(size(xMmax,1),2);
                
                % Assemble maximum bending moment matrix
                for i = 1:size(xMmax,1)  % this loop is of 2 rounds max
                    if isreal(xMmax(i)) && xMmax(i) > 0 && xMmax(i) < L
                        xx = xMmax(i);
                        Mmax(i,:) = [(A*xx^3 + B*xx^2 + C*xx + D), xx];
                    end
                end
                
                % Avoid getting the same maximum bending moment value twice
                if size(Mmax,1) == 2
                    if all(Mmax(1,:) == 0) && all(Mmax(2,:) == 0)
                        Mmax = [];
                    elseif all(Mmax(1,:) == 0) && ~all(Mmax(2,:) == 0)
                        Mmax(1,:) = [];
                    elseif ~all(Mmax(1,:) == 0) && all(Mmax(2,:) == 0)
                        Mmax(2,:) = [];
                    elseif Mmax(1,2) == Mmax(2,2)
                        Mmax(2,:) = [];
                    end
                end
                
            else
                M = linspace(M1,M2,size(x,2));
                Mmax = [];
            end
        end
        
        %------------------------------------------------------------------
        function [M,Mmax] = intDynamicBendingMoment_XY(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal forces value at both ends
            M1 = -elem.bending_moment_Z(:,1);
            M2 =  elem.bending_moment_Z(:,2);
            
            % Obs.: For now, only dynamic nodal loads are being considered,
            % thus bending moment is varies linearly along elements.
            
            % Avoid numeric garbage
            if all(abs(M1) < 1e-12)
                M1 = zeros(size(M1,1),1);
            end
            if all(abs(M2) < 1e-12)
                M2 = zeros(size(M2,1),1);
            end
            
            M = zeros(size(M1,1),size(x,2));
            M(:,  1) = M1;
            M(:,end) = M2;
            M(:,2:end-1) = sparse(1:size(M1,1),1:size(M1,1),M1') * ones(size(M1,1),size(x,2)-2) +...
                           (M2 - M1) / L * x(2:end-1);
            Mmax = [];
        end
        
        %------------------------------------------------------------------
        function [M,Mmax] = intBendingMoment_XZ(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal forces value at both ends
            M1 = elem.bending_moment_Y(1);
            M2 = -elem.bending_moment_Y(2);
            Q1 = elem.shear_force_Z(1);
            
            % Avoid numeric garbage
            if abs(M1) < 1e-10
                M1 = 0;
            end
            if abs(M2) < 1e-10
                M2 = 0;
            end
            
            if (isempty(elem.load.uniformLcl) == 0) || ...
               (isempty(elem.load.linearLcl) == 0)
                qui = 0;
                quf = 0;
                qli = 0;
                qlf = 0;
                
                if isempty(elem.load.uniformLcl) == 0
                    qui = elem.load.uniformLcl(3);
                    quf = qui;
                end
                
                if isempty(elem.load.linearLcl) == 0
                    qli = elem.load.linearLcl(3);
                    qlf = elem.load.linearLcl(6);
                end
                
                % Calculate total load value at both ends
                qi = qui + qli;
                qf = quf + qlf;
                
                % Load equation coefficients:
                % q(x) = Ax + B
                A = (qf - qi)/L;
                B = qi;
                
                % Bending moment equation coefficients:
                % M(x) = Ax^3 + Bx^2 + Cx + D
                A = A/6;
                B = B/2;
                C = Q1;
                D = M1;
                
                % Calculate internal force intensity on current
                % cross-section using bending moment equation
                M = A * x.^3 + B * x.^2 + C * x + D;
                
                % Get maximum coord
                xMmax = roots([3*A 2*B C]);
                
                % Initialize maximum bending moment matrix
                Mmax = zeros(size(xMmax,1),2);
                
                % Assemble maximum bending moment matrix
                for i = 1:size(xMmax,1)  % this loop is of 2 rounds max
                    if isreal(xMmax(i)) && xMmax(i) > 0 && xMmax(i) < L
                        xx = xMmax(i);
                        Mmax(i,:) = [(A*xx^3 + B*xx^2 + C*xx + D), xx];
                    end
                end
                
                % Avoid getting the same maximum bending moment value twice
                if size(Mmax,1) == 2
                    if all(Mmax(1,:) == 0) && all(Mmax(2,:) == 0)
                        Mmax = [];
                    elseif all(Mmax(1,:) == 0) && ~all(Mmax(2,:) == 0)
                        Mmax(1,:) = [];
                    elseif ~all(Mmax(1,:) == 0) && all(Mmax(2,:) == 0)
                        Mmax(2,:) = [];
                    elseif Mmax(1,2) == Mmax(2,2)
                        Mmax(2,:) = [];
                    end
                end
                
            else
                M = linspace(M1,M2,size(x,2));
                Mmax = [];
            end
        end
        
        %------------------------------------------------------------------
        function [M,Mmax] = intDynamicBendingMoment_XZ(elem,x)
            % Get element length
            L = elem.length;
            
            % Get element internal forces value at both ends
            M1 =  elem.bending_moment_Y(:,1);
            M2 = -elem.bending_moment_Y(:,2);
            
            % Obs.: For now, only dynamic nodal loads are being considered,
            % thus bending moment is varies linearly along elements.
            
            % Avoid numeric garbage
            if all(abs(M1) < 1e-12)
                M1 = zeros(size(M1,1),1);
            end
            if all(abs(M2) < 1e-12)
                M2 = zeros(size(M2,1),1);
            end
            
            M = zeros(size(M1,1),size(x,2));
            M(:,  1) = M1;
            M(:,end) = M2;
            %M(:,2:end-1) = M1 + (M2 - M1) / L * x(2:end-1);
            M(:,2:end-1) = sparse(1:size(M1,1),1:size(M1,1),M1') * ones(size(M1,1),size(x,2)-2) +...
                           (M2 - M1) / L * x(2:end-1);
            Mmax = [];
        end

        %------------------------------------------------------------------
        % Computes maximum and minimum values from transient analysis to
        % obtain internal forces envelop.
        function forcesEnvelop(elem)
            
            elem.intForcesEnvelop = zeros(2,size(elem.intStresses,2),size(elem.intStresses,3));
            
            for i = 1:size(elem.intStresses,3)
                
                maxF = max(elem.intStresses(:,:,i));
                aux = find(maxF < 0);
                if ~isempty(aux)
                    maxF(aux) = 0;
                end
                
                minF = min(elem.intStresses(:,:,i));
                aux = find(minF > 0);
                if ~isempty(aux)
                    minF(aux) = 0;
                end
                
                elem.intForcesEnvelop(:,:,i) = [ maxF ;
                                                 minF ];
            end
            
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of an Elem object.
        function clean(elem)
            elem.nen                = 2;
            elem.type               = 0;
            elem.anm                = [];
            elem.material           = [];
            elem.section            = [];
            elem.nodes              = [];
            elem.hingei             = 0;
            elem.hingef             = 0;
            elem.kri                = [];
            elem.krf                = [];
            elem.length             = 0;
            elem.cosine_X           = 0;
            elem.cosine_Y           = 0;
            elem.cosine_Z           = 0;
            elem.vz                 = [];
            elem.gle                = [];
            elem.rot                = [];
            elem.T                  = [];
            elem.kel                = [];
            elem.mel                = [];
            elem.mass_type          = 1;
            elem.mass_consideration = 1;
            elem.mass_mi            = 1;
            elem.fel_distribLoad    = [];
            elem.fel_thermalLoad    = [];
            elem.load               = [];
            elem.axial_force        = [];
            elem.shear_force_Y      = [];
            elem.shear_force_Z      = [];
            elem.bending_moment_Y   = [];
            elem.bending_moment_Z   = [];
            elem.torsion_moment     = [];
            elem.intDispl           = [];
            elem.natVibration       = [];
            elem.dynamicIntDispl    = [];
            elem.intCoords          = [];
            elem.intStresses        = [];
            elem.intForcesEnvelop   = [];
            elem.maxAxialForce      = [];
            elem.maxShearForce_XY   = [];
            elem.maxShearForce_XZ   = [];
            elem.maxBendMoment_XY   = [];
            elem.maxBendMoment_XZ   = [];
        end
    end
end