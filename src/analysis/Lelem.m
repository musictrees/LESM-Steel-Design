%% Lelem (Load Element) Class
%
%% Description
%
% This is a handle class for the definition of an element load.
%
% An element load object is responsible for defining the response of a
% linear element to loading effects, such as the computation of element
% fixed end forces (FEF), element equivalent nodal loads (ENL), and element
% internal displacements. There is a mutual relation between an object of
% the *Lelem* class and an object of the <elem.html *Elem*> class, since
% each element has its own load properties and each load object is
% associated with one element.
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
classdef Lelem < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        elem         = [];   % handle to an object of the Elem class
        elemLoadCase = [];   % array of element loads for each load case
        uniformDir   = 0;    % flag for uniform load direction (global system = 0, local system = 1)
        uniformGbl   = [];   % vector of uniformly distributed load components in global system [qx qy qz]
        uniformLcl   = [];   % vector of uniformly distributed load components in local system  [qx qy qz]
        linearDir    = 0;    % flag for linear load direction (global system = 0, local system = 1)
        linearGbl    = [];   % vector of linearly distributed load components in global system [qxi qyi qzi qxf qyf qzf]
        linearLcl    = [];   % vector of linearly distributed load components in local system  [qxi qyi qzi qxf qyf qzf]
        tempVar_X    = 0;    % temperature variation on element center of gravity
        tempVar_Y    = 0;    % temperature gradient relative to local y-axis
        tempVar_Z    = 0;    % temperature gradient relative to local z-axis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function load = Lelem(elem,elemLoadCase,unifdir,unifgbl,uniflcl,lindir,lingbl,linlcl,tvx,tvy,tvz)
            load.elem = elem;
            if (nargin > 1)
                if size(unifgbl,1) == 1 && size(unifgbl,2) > 1
                    unifgbl = unifgbl';
                end
                if size(uniflcl,1) == 1 && size(uniflcl,2) > 1
                    uniflcl = uniflcl';
                end
                if size(lingbl,1) == 1 && size(lingbl,2) > 1
                    lingbl = lingbl';
                end
                if size(linlcl,1) == 1 && size(linlcl,2) > 1
                    linlcl = linlcl';
                end
                load.elemLoadCase = elemLoadCase;
                load.uniformDir = unifdir;
                load.uniformGbl = unifgbl;
                load.uniformLcl = uniflcl;
                load.linearDir  = lindir;
                load.linearGbl  = lingbl;
                load.linearLcl  = linlcl;
                load.tempVar_X  = tvx;
                load.tempVar_Y  = tvy;
                load.tempVar_Z  = tvz;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Sets uniformly distributed load in global and local system.
        % Input arguments:
        %  unifload: vector of uniformly distributed load components
        %  dir:      specified direction of the distributed load components
        function setUnifLoad(load,unifload,dir)
            include_constants
            rot = load.elem.T;
            if size(unifload,1) == 1 && size(unifload,2) > 1
                unifload = unifload';
            end
            if isempty(load.uniformGbl)
                if dir == GLOBAL_LOAD
                    load.uniformGbl = unifload;
                    load.uniformLcl = rot * unifload;

                elseif dir == LOCAL_LOAD
                    load.uniformLcl = unifload;
                    load.uniformGbl = rot' * unifload;
                end
            else
                if dir == GLOBAL_LOAD
                    load.uniformGbl = load.uniformGbl + unifload;
                    load.uniformLcl = load.uniformLcl + rot * unifload;

                elseif dir == LOCAL_LOAD
                    load.uniformLcl = load.uniformLcl + unifload;
                    load.uniformGbl = load.uniformGbl + rot' * unifload;
                end
            end 
            if all(load.uniformGbl == 0) && all(load.uniformLcl == 0)
                load.uniformGbl = [];
                load.uniformLcl = [];
            end
        end
        
        %------------------------------------------------------------------
        % Sets linearly distributed load in global and local system.
        % Input arguments:
        %  linearload: vector of linearly distributed load components
        %  dir:        specified direction of the distributed load components
        function setLinearLoad(load,linearload,dir)
            include_constants
            rot = blkdiag(load.elem.T,load.elem.T);
            if size(linearload,1) == 1 && size(linearload,2) > 1
                linearload = linearload';
            end
            if isempty(load.linearGbl)
                if dir == GLOBAL_LOAD
                    load.linearGbl = linearload;
                    load.linearLcl = rot * linearload;

                elseif dir == LOCAL_LOAD
                    load.linearLcl = linearload;
                    load.linearGbl = rot' * linearload;
                end
            else
                if dir == GLOBAL_LOAD
                    load.linearGbl = load.linearGbl + linearload;
                    load.linearLcl = load.linearLcl + rot * linearload;

                elseif dir == LOCAL_LOAD
                    load.linearLcl = load.linearLcl + linearload;
                    load.linearGbl = load.linearGbl + rot' * linearload;
                end
            end 
            if all(load.linearGbl == 0) && all(load.linearLcl == 0)
                load.linearGbl = [];
                load.linearLcl = [];
            end
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load (ENL) vector in global
        % system for an applied distributed load.
        % Output:
        %  feg: equivalent nodal load vector in global system
        function feg = gblDistribLoadENL(load)
            % Compute and store element fixed end forces (FEF) in local
            % system for a distributed load
            load.elem.fel_distribLoad = load.elem.anm.elemLocDistribLoadFEF(load);
            
            % Transform element fixed end forces in local system to
            % equivalent nodal loads in global system
            feg = -load.elem.rot' * load.elem.fel_distribLoad;
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load (ENL) vector in global
        % system for an applied thermal load.
        % Output:
        %  feg: equivalent nodal load vector in global system
        function feg = gblThermalLoadENL(load)
            % Compute and store element fixed end forces (FEF) in local
            % system for a thermal load
            load.elem.fel_thermalLoad = load.elem.anm.elemLocThermalLoadFEF(load);
            
            % Transform element fixed end forces in local system to
            % equivalent nodal loads in global system
            feg = -load.elem.rot' * load.elem.fel_thermalLoad;
        end
        
        %------------------------------------------------------------------
        % Generates element axial fixed end force (FEF) vector for an
        % applied linearly distributed load.
        % Output:
        %  fea: 2 position vector with axial FEF's:
        %       fea(1) -> axial force at initial node
        %       fea(2) -> axial force at final node
        function fea = axialDistribLoadFEF(load)
            % Initialize axial load values at end nodes
            qxi = 0;
            qxf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qxi = qxi + load.uniformLcl(1);
                qxf = qxi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qxi = qxi + load.linearLcl(1);
                qxf = qxf + load.linearLcl(4);
            end
            
            % Check if axial load is not null over element
            if (qxi ~= 0) || (qxf ~= 0)
                L = load.elem.length;
                
                % Separate uniform portion from linear portion of axial load
                qx0 = qxi;
                qx1 = qxf - qxi;
                
                % Calculate fixed end force vector
                fea = [ -(qx0*L/2 + qx1*L/6);
                        -(qx0*L/2 + qx1*L/3) ];
               
            else
                
                fea = [ 0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralDistribLoadFEF_XY(load)
            include_constants;
            
            % Initialize transversal load values at end nodes
            qyi = 0;
            qyf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qyi = qyi + load.uniformLcl(2);
                qyf = qyi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qyi = qyi + load.linearLcl(2);
                qyf = qyf + load.linearLcl(5);
            end
            
            % Check if transversal load is not null over element
            if (qyi ~= 0) || (qyf ~= 0)
                % Basic element properties
                L  = load.elem.length;
                L2 = L^2;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    E  = load.elem.material.elasticity;
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_y;
                    I  = load.elem.section.inertia_z;
                    
                    Omega = E * I / (G * As * L2);
                end
                
                % Auxiliary parameters
                mu       = 1 + 12 * Omega;
                lambda   = 1 + 3  * Omega;
                zeta     = 1 + 40 * Omega/3;
                xi       = 1 + 5  * Omega;
                eta      = 1 + 15 * Omega;
                vartheta = 1 + 4  * Omega;
                psi      = 1 + 12 * Omega/5;
                varpi    = 1 + 20 * Omega/9;
                epsilon  = 1 + 80 * Omega/7;
                varrho   = 1 + 10 * Omega;
                upsilon  = 1 + 5  * Omega/2;
                varsigma = 1 + 40 * Omega/11;
               
                % Separate uniform portion from linear portion of tranversal load
                qy0 = qyi;
                qy1 = qyf - qyi;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ -(qy0*L/2   + qy1*3*L*zeta/(20*mu));
                            -(qy0*L2/12 + qy1*L2*eta/(30*mu));
                            -(qy0*L/2   + qy1*7*L*epsilon/(20*mu));
                              qy0*L2/12 + qy1*L2*varrho/(20*mu) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ -(qy0*3*L*vartheta/(8*lambda) + qy1*L*xi/(10*lambda));
                              0;
                            -(qy0*5*L*psi/(8*lambda)      + qy1*2*L*upsilon/(5*lambda));
                              qy0*L2/(8*lambda)           + qy1*L2/(15*lambda) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qy0*5*L*psi/(8*lambda)      + qy1*9*L*varpi/(40*lambda));
                            -(qy0*L2/(8*lambda)           + qy1*7*L2/(120*lambda));
                            -(qy0*3*L*vartheta/(8*lambda) + qy1*11*L*varsigma/(40*lambda));
                              0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qy0*L/2 + qy1*L/6);
                              0;
                            -(qy0*L/2 + qy1*L/3);
                              0 ];
                end
                
            else

                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralDistribLoadFEF_XZ(load)
            include_constants;
            
            % Initialize transversal load values at end nodes
            qzi = 0;
            qzf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qzi = qzi + load.uniformLcl(3);
                qzf = qzi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qzi = qzi + load.linearLcl(3);
                qzf = qzf + load.linearLcl(6);
            end
            
            % Check if transversal load is not null over element
            if (qzi ~= 0) || (qzf ~= 0)
                % Basic element properties
                L  = load.elem.length;
                L2 = L^2;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    E  = load.elem.material.elasticity;
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_z;
                    I  = load.elem.section.inertia_y;
                    
                    Omega = E * I / (G * As * L2);
                end
                
                % Auxiliary parameters
                mu       = 1 + 12 * Omega;
                lambda   = 1 + 3  * Omega;
                zeta     = 1 + 40 * Omega/3;
                xi       = 1 + 5  * Omega;
                eta      = 1 + 15 * Omega;
                vartheta = 1 + 4  * Omega;
                psi      = 1 + 12 * Omega/5;
                varpi    = 1 + 20 * Omega/9;
                epsilon  = 1 + 80 * Omega/7;
                varrho   = 1 + 10 * Omega;
                upsilon  = 1 + 5  * Omega/2;
                varsigma = 1 + 40 * Omega/11;
                
                % Separate uniform portion from linear portion of tranversal load
                qz0 = qzi;
                qz1 = qzf - qzi;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ -(qz0*L/2   + qz1*3*L*zeta/(20*mu));
                              qz0*L2/12 + qz1*L2*eta/(30*mu);
                            -(qz0*L/2   + qz1*7*L*epsilon/(20*mu));
                            -(qz0*L2/12 + qz1*L2*varrho/(20*mu)) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ -(qz0*3*L*vartheta/(8*lambda) + qz1*L*xi/(10*lambda));
                              0;
                            -(qz0*5*L*psi/(8*lambda)      + qz1*2*L*upsilon/(5*lambda));
                            -(qz0*L2/(8*lambda)           + qz1*L2/(15*lambda)) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qz0*5*L*psi/(8*lambda)      + qz1*9*L*varpi/(40*lambda));
                              qz0*L2/(8*lambda)           + qz1*7*L2/(120*lambda);
                            -(qz0*3*L*vartheta/(8*lambda) + qz1*11*L*varsigma/(40*lambda));
                              0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qz0*L/2 + qz1*L/6);
                              0;
                            -(qz0*L/2 + qz1*L/3);
                              0 ];
                end
                
            else

                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element axial fixed end force (FEF) vector for an
        % applied thermal load.
        % Output:
        %  fea: 2 position vector with axial FEF's:
        %       fea(1,1) -> axial force at initial node
        %       fea(2,1) -> axial force at final node
        function fea = axialThermalLoadFEF(load)
            % Get temperature variation on element center of gravity
            dtx = load.tempVar_X;
            
            % Check if temperature variation is not null
            if dtx ~= 0
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                A     = load.elem.section.area_x;
                
                % Calculate fixed end force vector
                fea = [ E * A * alpha * dtx;
                       -E * A * alpha * dtx ];
            
            else
                
                fea = [ 0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralThermalLoadFEF_XY(load)
            include_constants;
            
            % Get temperature gradient relative to element local y-axis
            dty = load.tempVar_Y;
            
            % Check if temperature gradient is not null
            if dty ~= 0
                % Basic element properties
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                I     = load.elem.section.inertia_z;
                h     = load.elem.section.height_y;
                L     = load.elem.length;
                EI    = E * I;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_y;
                    
                    Omega = E * I / (G * As * L * L);
                end
                
                % Auxiliary parameter
                lambda = 1 + 3 * Omega;
                
                % Compute unitary dimensionless temperature gradient
                tg = (alpha * dty) / h;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ 0;
                            tg * (EI);
                            0;
                            tg * (-EI) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ tg * (-3*EI/(2*L*lambda));
                            0;
                            tg * (3*EI/(2*L*lambda));
                            tg * (-3*EI/(2*lambda)) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ tg * (3*EI/(2*L*lambda));
                            tg * (3*EI/(2*lambda));
                            tg * (-3*EI/(2*L*lambda));
                            0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ 0;
                            0;
                            0;
                            0 ];
                end
                
            else
                
                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralThermalLoadFEF_XZ(load)
            include_constants;
            
            % Get temperature gradient relative to element local z-axis
            dtz = load.tempVar_Z;
            
            % Check if temperature gradient is not null
            if dtz ~= 0
                % Basic element properties
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                I     = load.elem.section.inertia_y;
                h     = load.elem.section.height_z;
                L     = load.elem.length;
                EI    = E * I;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_z;
                    
                    Omega = E * I / (G * As * L * L);
                end
                
                % Auxiliary parameter
                lambda = 1 + 3 * Omega;
                
                % Compute unitary dimensionless temperature gradient
                tg = (alpha * dtz) / h;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ 0;
                            tg * (-EI);
                            0;
                            tg * (EI) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    fef = [ tg * (-3*EI/(2*L*lambda));
                            0;
                            tg * (3*EI/(2*L*lambda));
                            tg * (3*EI/(2*lambda)) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ tg * (3*EI/(2*L*lambda));
                            tg * (-3*EI/(2*lambda));
                            tg * (-3*EI/(2*L*lambda));
                            0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ 0;
                            0;
                            0;
                            0 ];
                end
                
            else
                
                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Computes element axial displacement at given cross-section position
        % for an applied linearly distributed load assuming a fixed end
        % condition.
        % Output:
        %  du: axial displacement in longitudinal direction
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function du = axialDistribLoadDispl(load,x)
            include_constants;
            
            % Initialize axial load values at end nodes
            qxi = 0;
            qxf = 0;
            
            % Add uniform load contribution in local system
            if isempty(load.uniformLcl) == 0
                qxi = qxi + load.uniformLcl(1);
                qxf = qxi;
            end
            
            % Add linear load contribution in local system
            if isempty(load.linearLcl) == 0
                qxi = qxi + load.linearLcl(1);
                qxf = qxf + load.linearLcl(4);
            end
            
            % Check if transversal load is not null over element
            if (qxi ~= 0) || (qxf ~= 0)
                E = load.elem.material.elasticity;
                A = load.elem.section.area_x;
                L = load.elem.length;

                EA = E  * A;

                % Separate uniform portion from linear portion of axial load
                qx0 = qxi;
                qx1 = qxf - qxi;

                % Calculate axial displacements from uniform load portion (up0)
                % and from linear axial load portion (up1)
                up0 = qx0/EA * (L*x./2 - x.^2/2);
                up1 = qx1/EA * (L*x./6 - x.^3/(6*L));
                
                du = up0 + up1;
                
            else
                du = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dv = flexuralDistribLoadDispl_XY(load,x)
            include_constants;
            
            % Initialize transversal load value at end nodes
            qyi = 0;
            qyf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qyi = qyi + load.uniformLcl(2);
                qyf = qyi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qyi = qyi + load.linearLcl(2);
                qyf = qyf + load.linearLcl(5);
            end
            
            % Check if transversal load is not null over element
            if (qyi ~= 0) || (qyf ~= 0)
                % Basic element properties
                E  = load.elem.material.elasticity;
                I  = load.elem.section.inertia_z;
                EI  = E  * I;
                L  = load.elem.length;
                L2  = L  * L;
                L3  = L2 * L;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_y;
                    
                    Omega = EI / (G * As * L2);
                end
                
                % Auxiliary parameters
                mu       = 1 + 12 * Omega;
                lambda   = 1 + 3  * Omega;
                zeta     = 1 + 40 * Omega/3;
                xi       = 1 + 5  * Omega;
                eta      = 1 + 15 * Omega;
                vartheta = 1 + 4  * Omega;
                psi      = 1 + 12 * Omega/5;
                varpi    = 1 + 20 * Omega/9;

                % Separate uniform portion from linear portion of tranversal load
                qy0 = qyi;
                qy1 = qyf - qyi;

                % Calculate transversal displacements from uniform load portion (vq0)
                % and from linear load portion (vq1)
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    vq0 = qy0/EI * (L3*Omega*x./2 + (L2/24 - L2*Omega/2)*x.^2 - L*x.^3/12 + x.^4/24);
                    vq1 = qy1/EI * (3*L3*Omega*zeta*x./(20*mu) + L2*eta*x.^2/(60*mu) - (L*zeta/(40*mu) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    vq0 = qy0/EI * ((L3*mu/(48*lambda) + 3*L3*Omega*vartheta/(8*lambda))*x - L2*Omega*x.^2/2 - L*vartheta*x.^3/(16*lambda) + x.^4/24);
                    vq1 = qy1/EI * ((L3*eta/(120*lambda) + L3*Omega*xi/(10*lambda))*x - (L*xi/(60*lambda) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    vq0 = qy0/EI * (5*L3*Omega*psi*x./(8*lambda) + (L2/(16*lambda) - L2*Omega/2)*x.^2 - 5*L*psi*x.^3/(48*lambda) + x.^4/24);
                    vq1 = qy1/EI * (9*L3*Omega*varpi*x./(40*lambda) + 7*L2*x.^2/(240*lambda) - (3*L*varpi/(80*lambda) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    vq0 = qy0/EI * ((L3/24 + L3*Omega/2)*x - L2*Omega*x.^2/2 - L*x.^3/12 + x.^4/24);
                    vq1 = qy1/EI * ((7*L3/360 + L3*Omega/6)*x - (L/36 + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                end
                
                dv = vq0 + vq1;
                
            else
                dv = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dw = flexuralDistribLoadDispl_XZ(load,x)
            include_constants;
            
            % Initialize transversal load value at end nodes
            qzi = 0;
            qzf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qzi = qzi + load.uniformLcl(3);
                qzf = qzi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qzi = qzi + load.linearLcl(3);
                qzf = qzf + load.linearLcl(6);
            end
            
            % Check if transversal load is not null over element
            if (qzi ~= 0) || (qzf ~= 0)
                % Basic element parameters
                E  = load.elem.material.elasticity;
                I  = load.elem.section.inertia_y;
                EI = E  * I;
                L  = load.elem.length;
                L2 = L  * L;
                L3 = L2 * L;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_z;
                    
                    Omega = EI / (G * As * L2);
                end
                
                % Auxiliary parameters
                mu       = 1 + 12 * Omega;
                lambda   = 1 + 3  * Omega;
                zeta     = 1 + 40 * Omega/3;
                xi       = 1 + 5  * Omega;
                eta      = 1 + 15 * Omega;
                vartheta = 1 + 4  * Omega;
                psi      = 1 + 12 * Omega/5;
                varpi    = 1 + 20 * Omega/9;

                % Separate uniform portion from linear portion of tranversal load
                qz0 = qzi;
                qz1 = qzf - qzi;

                % Calculate transversal displacements from uniform load portion (wq0)
                % and from linear load portion (wq1)
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    wq0 = qz0/EI * (L3*Omega*x./2 + (L2/24 - L2*Omega/2)*x.^2 - L*x.^3/12 + x.^4/24);
                    wq1 = qz1/EI * (3*L3*Omega*zeta*x./(20*mu) + L2*eta*x.^2/(60*mu) - (L*zeta/(40*mu) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    wq0 = qz0/EI * ((L3*mu/(48*lambda) + 3*L3*Omega*vartheta/(8*lambda))*x - L2*Omega*x.^2/2 - L*vartheta*x.^3/(16*lambda) + x.^4/24);
                    wq1 = qz1/EI * ((L3*eta/(120*lambda) + L3*Omega*xi/(10*lambda))*x - (L*xi/(60*lambda) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    wq0 = qz0/EI * (5*L3*Omega*psi*x./(8*lambda) + (L2/(16*lambda) - L2*Omega/2)*x.^2 - 5*L*psi*x.^3/(48*lambda) + x.^4/24);
                    wq1 = qz1/EI * (9*L3*Omega*varpi*x./(40*lambda) + 7*L2*x.^2/(240*lambda) - (3*L*varpi/(80*lambda) + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    wq0 = qz0/EI * ((L3/24 + L3*Omega/2)*x - L2*Omega*x.^2/2 - L*x.^3/12 + x.^4/24);
                    wq1 = qz1/EI * ((7*L3/360 + L3*Omega/6)*x - (L/36 + L*Omega/6)*x.^3 + x.^5/(120*L));
                    
                end
                
                dw = wq0 + wq1;
                
            else
                dw = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element axial displacement at given cross-section position
        % for an applied thermal load assuming a fixed end condition.
        % Output:
        %  du: axial displacement in longitudinal direction
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function du = axialThermalLoadDispl(~,~)
            % Axial displacement is null when an element with fixed ends is
            % subjected to a temperature variation
            du = 0;
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dv = flexuralThermalLoadDispl_XY(load,x)
            include_constants;
            
            % Get temperature gradient relative to element local y-axis
            dty = load.tempVar_Y;
            
            % Check if temperature gradient is not null
            if dty ~= 0
                % Basic element parameters
                alpha = load.elem.material.thermExp;
                h     = load.elem.section.height_y;
                L     = load.elem.length;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    E  = load.elem.material.elasticity;
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_y;
                    I  = load.elem.section.inertia_z;
                    
                    Omega = E * I / (G * As * L * L);
                end
                
                % Auxiliary parameters
                mu     = 1 + 12 * Omega;
                lambda = 1 + 3  * Omega;
                gamma  = 1 - 6  * Omega;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dty) / h;
                
                % Calculate transversal displacement
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    dv = 0;
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    dv = tg * ((-L*mu/(4*lambda) + 3*L*Omega/(2*lambda))*x + x.^2/2 - x.^3/(4*L*lambda));
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    dv = tg *(-3*L*Omega*x./(2*lambda) - gamma*x.^2/(4*lambda) + x.^3/(4*L*lambda));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    dv = tg * (-L*x./2 + x.^2/2);
                    
                end
                
            else
                dv = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dw = flexuralThermalLoadDispl_XZ(load,x)
            include_constants;
            
            % Get temperature gradient relative to element local z-axis
            dtz = load.tempVar_Z;
            
            % Check if temperature gradient is not null
            if dtz ~= 0
                % Basic element parameters
                alpha = load.elem.material.thermExp;
                h     = load.elem.section.height_z;
                L     = load.elem.length;
                
                % Timoshenko parameter
                if load.elem.type == 0     % Navier element
                    Omega = 0;
                elseif load.elem.type == 1 % Timoshenko element
                    E  = load.elem.material.elasticity;
                    G  = load.elem.material.shear;
                    As = load.elem.section.area_z;
                    I  = load.elem.section.inertia_y;
                    
                    Omega = E * I / (G * As * L * L);
                end
                
                % Auxiliary parameters
                mu     = 1 + 12 * Omega;
                lambda = 1 + 3  * Omega;
                gamma  = 1 - 6  * Omega;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dtz) / h;
                
                % Calculate transversal displacement
                if (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    dw = 0;
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END || load.elem.hingef == SEMIRIGID_END)
                    
                    dw = tg * ((-L*mu/(4*lambda) + 3*L*Omega/(2*lambda))*x + x.^2/2 - x.^3/(4*L*lambda));
                    
                elseif (load.elem.hingei == CONTINUOUS_END || load.elem.hingei == SEMIRIGID_END) && (load.elem.hingef == HINGED_END)
                    
                    dw = tg *(-3*L*Omega*x./(2*lambda) - gamma*x.^2/(4*lambda) + x.^3/(4*L*lambda));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    dw = tg * (-L*x./2 + x.^2/2);
                    
                end
                
            else
                dw = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Clean data structure of a Lelem object.
        function clean(load)
            load.elem         = [];
            load.elemLoadCase = [];
            load.uniformDir   = 0;
            load.uniformGbl   = [];
            load.uniformLcl   = [];
            load.linearDir    = 0;
            load.linearGbl    = [];
            load.linearLcl    = [];
            load.tempVar_X    = 0;
            load.tempVar_Y    = 0;
            load.tempVar_Z    = 0;
        end
    end
end