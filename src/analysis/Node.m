%% Node Class
%
%% Description
%
% This is a handle class for the definition of a node.
%
% A node is a joint between two or more elements, or any element end, used
% to discretize the model. It is always considered as a tri-dimensional
% entity which may have applied loads or prescribed displacements.
%
classdef Node < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id             =  0;    % identification number
        coord          = [];    % vector of coordinates on global system [X Y Z]
        ebc            = [];    % vector of essential boundary condition flags [dx dy dz rx ry rz]: -1 = fictitious, 0 = free, 1 = fixed, 2 = spring
        nodalLoadCase  = [];    % array of nodal loads and prescribed displacements for each load case [fx fy fz mx my mz dx dy dz rx ry rz]'
        prescDispl     = [];    % vector of prescribed displacement values [dx dy dz rx ry rz]
        springStiff    = [];    % vector of spring stiffness coefficients [kx ky kz krx kry krz]
        isInclinedSupp = false; % flag for inclined support
        inclSuppDir    = [];    % vector of inclined support direction [dir_x dir_y dir_z]
        T              = [];    % basis rotation transformation matrix between two coordinate systems (global to inclined supp)
        inclSupp_vy    = [];    % auxiliar vy node local axis y orientation versor
        initCond       = [];    % array of initial conditions for dynamic analysis
        
        % ATTENTION
        % In the future, all nodal load info will be here
        load = [];    % vector of handles to Lnode objects
        
        % Addition: Concentrated mass. Not considering concentrated rotational inertia for now.
        displMass = 0.0;
        rotMass   = 0.0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node(id,coord,ebc,nodalLoadCase,nodalLoad,presDispl,springStiff,inclSupp,dir)
            if (nargin > 0)
                node.id = id;
                node.coord = coord;
                node.ebc = ebc;
                node.nodalLoadCase = nodalLoadCase;
                node.prescDispl = presDispl;
                node.springStiff = springStiff;
                node.load = Lnode(node);
                node.load.static = nodalLoad;
                if (nargin > 7)
                    node.isInclinedSupp = inclSupp;
                    node.inclSuppDir = dir;
                    node.computeRotMtxInclinedSupport();
                end
            end
            node.initCond = sparse(6,2);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Counts total number of elements and number of hinged elements
        % connected to a node.
        % Output:
        %  tot: total number of elements connected to a node
        %  hng: number of hinged elements connected to a node
        % Input arguments:
        %  model: handle to an object of the model class
        function [tot,hng] = elemsIncidence(node,model)
            % Initialize values
            tot = 0;
            hng = 0;
            
            for e = 1:model.nel
                % Check if initial node of current element is the target one
                if model.elems(e).nodes(1).id == node.id
                    tot = tot + 1;
                    if model.elems(e).hingei == 0
                        hng = hng + 1;
                    end
                
                % Check if final node of current element is the target one
                elseif model.elems(e).nodes(2).id == node.id
                    tot = tot + 1;
                    if model.elems(e).hingef == 0
                        hng = hng + 1;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Compute inclined support local axis
        function [x,y,z] = getInclinedSuppLocAxis(node,defaultAxisFlag)
            % Check if a flag was provided
            % The defaultAxisFlag aims to signalize that the output should
            % be as if there was no specified vy direction versor
            if nargin < 2
                defaultAxisFlag = false;
            end
            
            % Get node local axis x direction
            x = node.inclSuppDir / norm(node.inclSuppDir);
            
            % Check if there isn't a specified orientation versor 
            if isempty(node.inclSupp_vy) || defaultAxisFlag
                % Auxiliar vz node local axis z orientation versor
                if abs(x(1)) <= (10^-10) && abs(x(2)) <= (10^-10) && abs(x(3)) > (10^-10)
                    if x(3) > 0
                        vz = [-1, 0, 0];
                    else
                        vz = [1, 0, 0];
                    end
                else
                    vz = [0, 0, 1];
                end

                % Compute local axis y and z directions
                y = cross(vz,x);
                y = y / norm(y);
                z = cross(x,y);
                z = z / norm(z);
                
            else
                % Compute local axis y and z directions
                z = cross(x,node.inclSupp_vy);
                z = z / norm(z);
                y = cross(z,x);
                y = y / norm(y);
            end
        end
        
        %------------------------------------------------------------------
        % Compute inclined support rotation rotation matrix 
        function computeRotMtxInclinedSupport(node)
            % Get node local axis
            [x,y,z] = getInclinedSuppLocAxis(node);
            
            % Assemble basis transformation matrix
            node.T = [x;
                      y;
                      z];
        end
        
        %------------------------------------------------------------------
        % Adds inclined displacement restraint to node
        % Input arguments:
        %  dir: vector of inclined support direction [dir_x dir_y dir_z]
        %  vy: auxiliar vy node local axis y orientation versor
        function setInclinedSupp(node,dir,vy)
            if nargin ~= 3
                vy = [];
            end
            node.isInclinedSupp = true;
            node.inclSuppDir = dir / norm(dir);
            node.inclSupp_vy = vy / norm(vy);
            node.computeRotMtxInclinedSupport();
        end
        
        %------------------------------------------------------------------
        % Removes inclined displacement restraints from node
        function removeInclinedSupp(node)
            node.isInclinedSupp = false;
            node.inclSuppDir = [];
            node.T = [];
            node.inclSupp_vy = [];
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Node object.
        function clean(node)
            node.id             =  0;
            node.coord          = [];
            node.ebc            = [];
            node.nodalLoadCase  = [];
            node.prescDispl     = [];
            node.springStiff    = [];
            node.removeInclinedSupp;
            node.initCond       = [];
            node.load           = [];
        end
    end
end