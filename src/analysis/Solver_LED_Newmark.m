%% Linear Elastic Dynamic Solver Class (Newmark Method)
%
%% Description
%
% This is a sub-class of the <solver.html *Solver*> class for the
% implementation of the *Newmark* solver algorithm.
%
classdef Solver_LED_Newmark < Solver
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver_LED_Newmark(drv)
            solver = solver@Solver(drv);
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
            nf = model.neqfree + model.neqspring;

            % Assemble free-free global stiffness matrix
            Kff = model.K(1:nf , 1:nf);
            
            % Assemble free-free global mass matrix
            Mff = model.M(1:nf , 1:nf);
            
            % Assemble free-free global damping matrix
            Cff = model.C(1:nf , 1:nf);
            
            % Assemble free global forcing matrix
            Ff = model.F(1:nf  ,  :  );
            
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
            a(1:nf,1) = Mff \(- Cff*v(1:nf,1)  - Kff*d(1:nf,1));
            
            % Set initial conditions for forced + free vibration
            dp(:,1) = d(:,1);
            vp(:,1) = v(:,1);
            ap(1:nf,1) = Mff \(Ff(:,1) - Cff*v(1:nf,1)  - Kff*d(1:nf,1));
            
            % Compute time step
            step = model.t/model.n_steps;
            
            % Compute auxiliar matrices
            Aux_B = (2/step) * Cff + (2/step)^2 * Mff;
            Aux_A = Kff + Aux_B;
            Aux_C = Cff + (4/step) * Mff;
            
            % Initialize counter
            count = 1;
            
            % Loop through all steps
            for i = step:step:model.t
                % Update counter
                count = count + 1;
                
                % Compute displacement on the next step
                u  = Aux_A \ (Aux_B * d(1:nf,count-1)  + Aux_C * v(1:nf,count-1)  + Mff * a(1:nf,count-1));
                up = Aux_A \ (Aux_B * dp(1:nf,count-1) + Aux_C * vp(1:nf,count-1) + Mff * ap(1:nf,count-1) + Ff(:,count));
                
                % Compute speed on the next step
                w  = (2/step) * (u  - d(1:nf,count-1))  - v(1:nf,count-1);
                wp = (2/step) * (up - dp(1:nf,count-1)) - vp(1:nf,count-1);
                
                % Compute acceleration on the next step
                z  = (2/step)^2 * (u  - d(1:nf,count-1))  - (4/step) * v(1:nf,count-1)  - a(1:nf,count-1);
                zp = (2/step)^2 * (up - dp(1:nf,count-1)) - (4/step) * vp(1:nf,count-1) - ap(1:nf,count-1);
                
                % Store computed components to matrices
                d(1:nf,count)  = u;
                dp(1:nf,count) = up;
                v(1:nf,count)  = w;
                vp(1:nf,count) = wp;
                a(1:nf,count)  = z;
                ap(1:nf,count) = zp;
            end
            
            % Store results
            model.results.dynamicDispl       = dp;
            model.results.dynamicVeloc       = vp;
            model.results.dynamicAccel       = ap;
            model.results.dynamicDisplForced = dp - d;
            model.results.dynamicVelocForced = vp - v;
            model.results.dynamicAccelForced = ap - a;
            
            % Return success flag (will always be the case, unless a bug
            % occurs previously).
            status = true;
        end
    end
end