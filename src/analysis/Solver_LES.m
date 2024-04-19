%% Linear Elastic Static Solver Class
%
%% Description
%
% This is a sub-class of the <solver.html *Solver*> class for the
% implementation of the *Linear-Elastic Static* solver algorithm.
%
classdef Solver_LES < Solver
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver_LES(drv)
            solver = solver@Solver(drv);
        end
    end
    
    %% Public method
    % Implementation of the abstract method declared in super-class <solver.html *Solver*>.
    methods
        %------------------------------------------------------------------
        % Partitions and solves the system of equilibrium equations.
        % Partitions the coefficient matrix K, forcing vector F and unknown 
        % vector D:
        %  f -> free d.o.f.'s (numbered first) - natural B.C. (unknown)
        %  s -> spring d.o.f.'s (numbred after free d.o.f.'s) - natural B.C. (unknown)
        %  c -> constrainted by support d.o.f.'s (numbered later) - essential B.C. (known)
        % 
        % Note that spring d.o.f.'s have unknown displacement values, so 
        % they are treated as free d.o.f.'s on the next steps.
        %
        % Partitioned equilibrium system:
        % [ Kff Kfc ] * [ Df ] = [ Ff ]
        % [ Kcf Kcc ]   [ Dc ] = [ Fc ]
        %
        % Output:
        %  status: flag for stable free-free global matrix (0 = unstable, 1 = stable)
        function status = solve(solver)
            % Get handle to model
            model = solver.drv.model;
            
            % Assemble free-free global matrix
            Kff = model.K(1:model.neqfree+model.neqspring,1:model.neqfree+model.neqspring);
            
            % Check for stable Kff global matrix by verifying its
            % reciprocal condition number (a very low reciprocal condition
            % number indicates that the matrix is badly conditioned and
            % may be singular)
            status = 1;
            if (rcond(full(Kff)) < 10e-12)
                status = 0;
                return;
            end
            
            % Partition system of equations
            Kfc = model.K(1:model.neqfree+model.neqspring,model.neqfree+model.neqspring+1:model.neq);
            Kcf = model.K(model.neqfree+model.neqspring+1:model.neq,1:model.neqfree+model.neqspring);
            Kcc = model.K(model.neqfree+model.neqspring+1:model.neq,model.neqfree+model.neqspring+1:model.neq);

            Ff  = model.F(1:model.neqfree+model.neqspring);
            Fc  = model.F(model.neqfree+model.neqspring+1:model.neq); 

            Dc  = model.D(model.neqfree+model.neqspring+1:model.neq);
            
            % Solve for Df
            Df = Kff \ (Ff - Kfc * Dc);

            % Recover forcing unknown values (reactions) at essential B.C.
            % It is assumed that the Fc vector currently stores combined 
            % nodal loads applied directly to fixed d.o.f's.
            % Superimpose computed reaction values to combined nodal loads,
            % with inversed direction, that were applied directly to fixed 
            % d.o.f.'s.
            Fc = -Fc + Kcf * Df + Kcc * Dc;
            
            % Initialize spring reactions vector
            Fs = zeros(1,model.neqspring);
            
            % Evaluates reactions in spring supports
            for ns = 1:model.neqspring
                Fs(ns) = - model.gblSprStiff(ns) * Df(model.neqfree+ns);
            end    
            Ff(model.neqfree+1:model.neqfree+model.neqspring) = Fs;
            
            % Reconstruct the global force vector (F) the global unknown 
            % vector D
            model.D = [ Df
                        Dc ];
            
            model.F = [ Ff
                        Fc ];
            
            % Store full displacement and forcing vectors
            % Avoids trouble with sparse form in post processing functions
            model.D = full(model.D);
            model.F = full(model.F);
        end
    end
end
