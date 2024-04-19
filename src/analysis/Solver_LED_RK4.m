%% Linear Elastic Dynamic Solver Class (RK4 Method)
%
%% Description
%
% This is a sub-class of the <solver.html *Solver*> class for the
% implementation of the *4th order Runge-Kutta* solver algorithm.
%
classdef Solver_LED_RK4 < Solver
    %% Private properties
    properties (SetAccess = private, GetAccess = private)
        nf   =  0;
        step =  0;
        Mff  = [];
        Cff  = [];
        Kff  = [];
        Ff   = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver_LED_RK4(drv)
            solver = solver@Solver(drv);
        end
    end
    
    %% Private method
    % Evaluates derivative -> "f" from RK formulation
    methods (Access = private)
        %------------------------------------------------------------------
        function f = evalDeriv(solver,t,y,forceFlag)
            if ~forceFlag
                F = 0;
            else
                % Check if step is not an integer
                dt = rem(t,1);
                if t >= size(solver.Ff,2)
                    t = size(solver.Ff,2) - 1;
                    dt = 1;
                end
                F = (1-dt) * solver.Ff(:,floor(t)) + dt * solver.Ff(:,floor(t)+1);
            end
            
            % Get displ and veloc
            d = y(1:solver.nf);
            v = y(solver.nf+1:end);
            
            % Compute derivatives
            f = zeros(2*solver.nf,1);
            f(1:solver.nf) = v;
            f(solver.nf+1:end) = solver.Mff \ ( - solver.Kff * d  - solver.Cff * v + F);
        end
    end
    
    %% Public method
    % Implementation of the abstract method declared in super-class <solver.html *Solver*>.
    methods
        %------------------------------------------------------------------
        function status = solve(solver)
            % Get handle to model object
            model = solver.drv.model;
            
            % Compute number of free dofs
            solver.nf = model.neqfree + model.neqspring;

            % Assemble free-free global stiffness matrix
            solver.Kff = model.K(1:solver.nf , 1:solver.nf);
            
            % Assemble free-free global mass matrix
            solver.Mff = model.M(1:solver.nf , 1:solver.nf);
            
            % Assemble free-free global damping matrix
            solver.Cff = model.C(1:solver.nf , 1:solver.nf);
            
            % Assemble free global forcing matrix
            solver.Ff = model.F(1:solver.nf  ,  :  );
            
            % Initialize displacement, velocity and acceleration matrices
            d  = zeros(model.neq , model.n_steps + 1);
            dp = zeros(model.neq , model.n_steps + 1);
            v  = zeros(model.neq , model.n_steps + 1);
            vp = zeros(model.neq , model.n_steps + 1);
            a  = zeros(model.neq , model.n_steps + 1);
            ap = zeros(model.neq , model.n_steps + 1);
            
             % Set initial conditions for free vibration
            d(:,1) = model.c0(:,1);
            v(:,1) = model.c0(:,2);
            a(1:solver.nf,1) = solver.Mff \ ( - solver.Cff*v(1:solver.nf,1)...
                                              - solver.Kff*d(1:solver.nf,1));
                                                   
            % Set initial conditions for free + forced vibration
            dp(:,1) = d(:,1);
            vp(:,1) = v(:,1);
            ap(1:solver.nf,1) = solver.Mff \ (solver.Ff(:,1) - solver.Cff*v(1:solver.nf,1)...
                                                             - solver.Kff*d(1:solver.nf,1));
            
            % Compute time step
            solver.step = model.t/model.n_steps;
            
            % Initialize counter
            count = 1;
            
            % Loop through all steps
            for i = solver.step:solver.step:model.t
                % Update counter
                count = count + 1;
                
                % Compute fourth order Runge-Kutta parameters
                f1  = solver.evalDeriv(count-1,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)],false);
                f1p = solver.evalDeriv(count-1,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)],true);
                
                f2  = solver.evalDeriv(count-0.5,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)]+ f1*solver.step*0.5,false);
                f2p = solver.evalDeriv(count-0.5,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)]+f1p*solver.step*0.5,true);
                
                f3  = solver.evalDeriv(count-0.5,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)]+ f2*solver.step*0.5,false);
                f3p = solver.evalDeriv(count-0.5,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)]+f2p*solver.step*0.5,true);
                
                f4  = solver.evalDeriv(count,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)]+ f3*solver.step,false);
                f4p = solver.evalDeriv(count,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)]+f3p*solver.step,true);
                
                % Compute displacement on the next step
                u  =  d(1:solver.nf,count-1) + (1/6) * solver.step * ( f1(1:solver.nf,:) +  2*f2(1:solver.nf,:) +  2*f3(1:solver.nf,:) +  f4(1:solver.nf,:));
                up = dp(1:solver.nf,count-1) + (1/6) * solver.step * (f1p(1:solver.nf,:) + 2*f2p(1:solver.nf,:) + 2*f3p(1:solver.nf,:) + f4p(1:solver.nf,:));
                
                % Compute speed on the next step
                w  =  v(1:solver.nf,count-1) + (1/6) * solver.step * ( f1(solver.nf+1:end,:) +  2*f2(solver.nf+1:end,:) +  2*f3(solver.nf+1:end,:) +  f4(solver.nf+1:end,:));
                wp = vp(1:solver.nf,count-1) + (1/6) * solver.step * (f1p(solver.nf+1:end,:) + 2*f2p(solver.nf+1:end,:) + 2*f3p(solver.nf+1:end,:) + f4p(solver.nf+1:end,:));
                
                % Compute acceleration on the next step
                z  = (1/6) * ( f1(solver.nf+1:end,:) +  2*f2(solver.nf+1:end,:) +  2*f3(solver.nf+1:end,:) +  f4(solver.nf+1:end,:));
                zp = (1/6) * (f1p(solver.nf+1:end,:) + 2*f2p(solver.nf+1:end,:) + 2*f3p(solver.nf+1:end,:) + f4p(solver.nf+1:end,:));
                
                % Store computed components to matrices
                d(1:solver.nf,count)  = u;
                dp(1:solver.nf,count) = up;
                v(1:solver.nf,count)  = w;
                vp(1:solver.nf,count) = wp;
                a(1:solver.nf,count)  = z;
                ap(1:solver.nf,count) = zp;
            end
            
            % Store results
            model.results.dynamicDispl        = dp;
            model.results.dynamicVeloc        = vp;
            model.results.dynamicAccel        = ap;
            model.results.dynamicDisplForced  = dp - d;
            model.results.dynamicVelocForced  = vp - v;
            model.results.dynamicAccelForced  = ap - a;
            
            % Clear private properties
            solver.nf   =  0;
            solver.step =  0;
            solver.Mff  = [];
            solver.Cff  = [];
            solver.Kff  = [];
            solver.Ff   = [];
            
            % Return success flag (will always be the case, unless a bug
            % occurs previously).
            status = true;
        end
    end
end