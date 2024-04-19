%% Result Class
%
%% Description
%
% This is a handle class for the definition of a result storage.
%
% Objects of this class are responsible for storing results after data is
% processed by a drv object.
%
classdef Result < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type               = [];   % flag for analysis type (static, numeric dynamic, analytical dynamic, etc.)
        dynamicDispl       = [];   % matrix of total dynamic displacement results [neq, n_steps+1, n_modes]
        dynamicVeloc       = [];   % matrix of total dynamic velocity results [neq, n_steps+1, n_modes]
        dynamicAccel       = [];   % matrix of total dynamic acceleration results [neq, n_steps+1, n_modes]
        dynamicDisplForced = [];   % matrix of forced dynamic displacement results [neq, n_steps+1, n_modes]
        dynamicVelocForced = [];   % matrix of forced dynamic velocity results [neq, n_steps+1, n_modes]
        dynamicAccelForced = [];   % matrix of forced dynamic acceleration results [neq, n_steps+1, n_modes]
    end
    
    %% Constructor method
    methods
        function results = Result(type)
            if (nargin > 0)
                results.type = type;
            end
        end
    end
    
    %% Public methods
    methods 
        function clean(results)
            results.type               = [];
            results.dynamicDispl       = [];
            results.dynamicVeloc       = [];
            results.dynamicAccel       = [];
            results.dynamicDisplForced = [];
            results.dynamicVelocForced = [];
            results.dynamicAccelForced = [];
        end
    end
end