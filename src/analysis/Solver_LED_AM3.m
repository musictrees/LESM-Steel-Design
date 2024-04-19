%% Linear Elastic Dynamic Solver Class (Adams-Moulton Method)
%
%% Description
%
% This is a sub-class of the <solver.html *Solver*> class for the
% implementation of the *3rd order Adams-Moulton* solver algorithm.
%
classdef Solver_LED_AM3 < Solver
    %% Private properties
    properties (SetAccess = private, GetAccess = private)
        nf     =  0;
        step   =  0;
        Mff    = [];
        invMff = [];
        Cff    = [];
        Kff    = [];
        Ff     = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver_LED_AM3(drv)
            solver = solver@Solver(drv);
        end
    end
    
    %% Private method
    methods (Access = private)
        % Evaluates derivative -> "f" from RK formulation
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
            f(solver.nf+1:end) = solver.invMff * ( - solver.Kff * d  - solver.Cff * v + F);
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
            if (rcond(full(solver.Mff)) < 10e-12)
                solver.invMff = pinv(solver.Mff);
            else
                solver.invMff = inv(solver.Mff);
            end
            
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
            
            % Compute auxiliar matrices
            Aux_B = (2/solver.step) * solver.Cff + (2/solver.step)^2 * solver.Mff;
            Aux_A = solver.Kff + Aux_B;
            Aux_C = solver.Cff + (4/solver.step) * solver.Mff;
            
            % Initialize counter
            count = 1;
            
            % Compute first step using Newmark
            %--------------------------------------------------------------
            % Update counter
            count = count + 1;

            % Compute displacement on the next step
            u  = Aux_A \ (Aux_B * d(1:solver.nf,count-1)  + Aux_C * v(1:solver.nf,count-1)  + solver.Mff * a(1:solver.nf,count-1));
            up = Aux_A \ (Aux_B * dp(1:solver.nf,count-1) + Aux_C * vp(1:solver.nf,count-1) + solver.Mff * ap(1:solver.nf,count-1) + solver.Ff(:,count));

            % Compute speed on the next step
            w  = (2/solver.step) * (u  - d(1:solver.nf,count-1))  - v(1:solver.nf,count-1);
            wp = (2/solver.step) * (up - dp(1:solver.nf,count-1)) - vp(1:solver.nf,count-1);

            % Compute acceleration on the next step
            z  = (2/solver.step)^2 * (u  - d(1:solver.nf,count-1))  - (4/solver.step) * v(1:solver.nf,count-1)  - a(1:solver.nf,count-1);
            zp = (2/solver.step)^2 * (up - dp(1:solver.nf,count-1)) - (4/solver.step) * vp(1:solver.nf,count-1) - ap(1:solver.nf,count-1);

            % Store computed components to matrices
            d(1:solver.nf,count)  = u;
            dp(1:solver.nf,count) = up;
            v(1:solver.nf,count)  = w;
            vp(1:solver.nf,count) = wp;
            a(1:solver.nf,count)  = z;
            ap(1:solver.nf,count) = zp;
            %--------------------------------------------------------------
            
            % Compute second step using second order Adams-Moulton
            %--------------------------------------------------------------
            if (2*solver.step) <= model.t
                % Update counter
                count = count + 1;

                % Compute two previous derivatives
                f0  = solver.evalDeriv(count-2,[ d(1:solver.nf,count-2); v(1:solver.nf,count-2)],false);
                f0p = solver.evalDeriv(count-2,[dp(1:solver.nf,count-2);vp(1:solver.nf,count-2)],true);

                f1  = solver.evalDeriv(count-1,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)],false);
                f1p = solver.evalDeriv(count-1,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)],true);

                A = [  eye(solver.nf)                                 -solver.step*(5/12)*eye(solver.nf)                               ;
                    (solver.step*(5/12))*(solver.invMff*solver.Kff)   (solver.step*(5/12))*(solver.invMff*solver.Cff) + eye(solver.nf) ];

                F = [ zeros(solver.nf,1) ; solver.step*(5/12)*solver.invMff*solver.Ff(:,count) ];

                % Compute displacement and speed on the next step
                x  = [  d(1:solver.nf,count-1); v(1:solver.nf,count-1) ] + (1/12) * solver.step * (8*f1 -  f0);
                xp = [ dp(1:solver.nf,count-1);vp(1:solver.nf,count-1) ] + (1/12) * solver.step * (8*f1p -  f0p) + F;

                y  = A \  x;
                yp = A \ xp;

                u  =  y(1:solver.nf,:);
                up = yp(1:solver.nf,:);

                w  =  y(solver.nf+1:end,:);
                wp = yp(solver.nf+1:end,:);

                % Compute acceleration on the next step
                z  = solver.invMff * ( - solver.Kff * u  - solver.Cff * w );
                zp = solver.invMff * ( - solver.Kff * up - solver.Cff * wp + solver.Ff(:,count));

                % Store computed components to matrices
                d(1:solver.nf,count)  = u;
                dp(1:solver.nf,count) = up;
                v(1:solver.nf,count)  = w;
                vp(1:solver.nf,count) = wp;
                a(1:solver.nf,count)  = z;
                ap(1:solver.nf,count) = zp;
            end
            %--------------------------------------------------------------
            
            % Compute derivative matrix prior to entering loop
            A = [  eye(solver.nf)                                   -solver.step*(9/24)*eye(solver.nf)                               ;
                  (solver.step*(9/24))*(solver.invMff*solver.Kff)   (solver.step*(9/24))*(solver.invMff*solver.Cff) + eye(solver.nf) ];
            
            % Loop through all remaining steps
            for i = (3*solver.step):solver.step:model.t
                % Update counter
                count = count + 1;
                
                % Compute two previous derivatives
                f0  = solver.evalDeriv(count-3,[ d(1:solver.nf,count-3); v(1:solver.nf,count-3)],false);
                f0p = solver.evalDeriv(count-3,[dp(1:solver.nf,count-3);vp(1:solver.nf,count-3)],true);
                
                f1  = solver.evalDeriv(count-2,[ d(1:solver.nf,count-2); v(1:solver.nf,count-2)],false);
                f1p = solver.evalDeriv(count-2,[dp(1:solver.nf,count-2);vp(1:solver.nf,count-2)],true);
                
                f2  = solver.evalDeriv(count-1,[ d(1:solver.nf,count-1); v(1:solver.nf,count-1)],false);
                f2p = solver.evalDeriv(count-1,[dp(1:solver.nf,count-1);vp(1:solver.nf,count-1)],true);
                  
                F = [ zeros(solver.nf,1) ; solver.step*(9/24)*solver.invMff*solver.Ff(:,count) ];
                
                % Compute displacement and speed on the next step
                x  = [  d(1:solver.nf,count-1); v(1:solver.nf,count-1) ] + (1/24) * solver.step * ( 19*f2 -  5*f1 +  f0 );
                xp = [ dp(1:solver.nf,count-1);vp(1:solver.nf,count-1) ] + (1/24) * solver.step * (19*f2p - 5*f1p +  f0p) + F;
                
                y  = A \  x;
                yp = A \ xp;
                
                u  =  y(1:solver.nf,:);
                up = yp(1:solver.nf,:);
                
                w  =  y(solver.nf+1:end,:);
                wp = yp(solver.nf+1:end,:);
                
                % Compute acceleration on the next step
                z  = solver.invMff * ( - solver.Kff * u  - solver.Cff * w );
                zp = solver.invMff * ( - solver.Kff * up - solver.Cff * wp + solver.Ff(:,count));
                
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
            solver.nf     =  0;
            solver.step   =  0;
            solver.Mff    = [];
            solver.invMff = [];
            solver.Cff    = [];
            solver.Kff    = [];
            solver.Ff     = [];
            
            % Return success flag (will always be the case, unless a bug
            % occurs previously).
            status = true;
        end
    end
end