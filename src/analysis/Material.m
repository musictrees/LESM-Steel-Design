%% Material Class
%
%% Description
%
% This is a handle class for the definition of a material.
%
% All materials are considered to have linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.
%
classdef Material < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id         = 0;   % identification number
        elasticity = 0;   % elasticity modulus
        poisson    = 0;   % poisson ratio
        shear      = 0;   % shear modulus
        thermExp   = 0;   % thermal expansion coefficient
        density    = 0;   % mass density
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function material = Material(id,e,v,te,rho)
            if (nargin > 0)
                material.id         = id;
                material.elasticity = e;
                material.poisson    = v;
                material.shear      = e / (2 * (1 + v));
                material.thermExp   = te;
                if nargin == 5
                    material.density = rho;
                end
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a Material object.
        function clean(material)
            material.id         = 0;
            material.elasticity = 0;
            material.poisson    = 0;
            material.shear      = 0;
            material.thermExp   = 0;
            material.density    = 0;
        end
    end
end