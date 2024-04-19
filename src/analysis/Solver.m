%% Solver Class
%
%% Description
%
% This is a handle super-class for the definition of a solver algorithm.
%
% The solver is created and called by an analysis driver to compute the
% solution of an assembled stuctural analysis problem.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <solver_les.html linear-elastic static solver>.
% * <solver_led_am3.html Adams-Moulton dynamic solver>.
% * <solver_led_modalsup.html modal superposition dynamic solver>.
% * <solver_led_newmark.html Newmark dynamic solver>.
% * <solver_led_rk4.html Runge-Kutta 4th order dynamic solver>.
% * <solver_led_wilson.html Wilson-theta dynamic solver>.
%
classdef Solver < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        drv = [];   % handle to a driver object
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver(drv)
            solver.drv = drv;
        end
    end
    
    %% Abstract method
    methods (Abstract)
        %------------------------------------------------------------------
        % Calls algorithm to solve structural analysis problem
        status = solve(solver)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of the algorithm object.
        function clean(solver)
            solver.drv = [];
        end
    end
end