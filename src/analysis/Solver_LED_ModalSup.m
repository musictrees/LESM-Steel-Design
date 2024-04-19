%% Linear Elastic Dynamic Solver Class (Modal Superposition - Uncoupled System of ODEs)
%
%% Description
%
% This is a sub-class of the <solver.html *Solver*> class for the
% implementation of the *Modal Superposition* solver algorithm.
%
classdef Solver_LED_ModalSup < Solver
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function solver = Solver_LED_ModalSup(drv)
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
            
            % Get free-free mass matrix
            Mff = model.M(1:nf ,1:nf);
            
            % Get free global forcing matrix
            Ff = model.F(1:nf  ,  :  );
            
            % Compute diagonal uncoupled modal stiffness, damping and mass
            Kd = diag(model.V' * model.K(1:nf , 1:nf) * model.V);
            Cd = diag(model.V' * model.C(1:nf , 1:nf) * model.V);
            Md = diag(model.V' *        Mff           * model.V);
            
            % Find initial null mass index
            % Null modal mass coefficients occur when dealing with element
            % lumped mass matrices. It means that no mass is attributed to
            % rotational dofs, and so modal behavior of such dofs should
            % not be considered.
            nullMass = find(abs(Md) <= 10^-15);
            if ~isempty(nullMass)
                n_modes =  nullMass(1) - 1;
            else
                n_modes = model.n_modes;
            end
            
            % Initialize displacement, velocity and acceleration matrices
            d = zeros(model.n_modes , model.n_steps + 1);
            v = zeros(model.n_modes , model.n_steps + 1);
            a = zeros(model.n_modes , model.n_steps + 1);            
                      
            % Compute time step
            step = model.t/model.n_steps;
            
            % Vector of time steps
            dt = step:step:model.t;
                     
            % Loop through all free DOFs
            for n = 1:n_modes
                % Get initial conditions
                d0 = model.V(:,n)' * Mff * model.c0(1:nf,1) / Md(n);
                v0 = model.V(:,n)' * Mff * model.c0(1:nf,2) / Md(n);
                a0 = (-Cd(n)*v0 - Kd(n)*d0) / Md(n);
                
                % Check if there is are initial displacement and speed
                % conditions
                if d0 == 0 && v0 == 0
                    d(n,:) = zeros(1,model.n_steps + 1);
                    v(n,:) = zeros(1,model.n_steps + 1);
                    a(n,:) = zeros(1,model.n_steps + 1);
                else
                    % Set initial conditions in vectors
                    d(n,1) = d0;
                    v(n,1) = v0;
                    a(n,1) = a0;
                    
                    % Compute critical damping ratio for mode n
                    xi_n = 0.5 * Cd(n) / (Md(n) * model.W(n));
                    
                    % Compute auxiliar parameter
                    ox2 = abs(1 - xi_n^2);
                    
                    % Compute initial phase angle
                    phi = atan(d0*model.W(n)*ox2^0.5 / (v0 + xi_n*model.W(n)*d0));
                    if isnan(phi)
                        phi = 0;
                    end
                    
                    % Compute initial displacement amplitude
                    criticalDamping = false;
                    if phi ~= 0
                        y = d0 / sin(phi);
                    elseif ox2 ~= 0
                        y_aux = ((v0/model.W(n)) + xi_n * d0) / (ox2^0.5);
                        y = (d0^2 + y_aux^2)^0.5;
                    else % critically damped system
                        criticalDamping = true;
                    end
                    
                    % Switch between non critically damped and critically
                    % damped system. This is necessary because it changes
                    % the analytical solution of the ODE of motion.
                    % Avoids going through trigonometric functions with
                    % null angular frequency on each time step.
                    if ~criticalDamping
                        
                        % Compute auxiliar parameters
                        wt   =          model.W(n) * dt;
                        xwt  =                xi_n * wt;
                        s    =  sin(ox2^0.5 * wt + phi);
                        c    =  cos(ox2^0.5 * wt + phi);

                        % Compute displacement on the next step
                        d(n,2:end) = y * exp(-xwt).* s;

                        % Compute speed on the next step
                        v(n,2:end) = model.W(n) * (y * ox2^0.5 * exp(-xwt).* c - xi_n * d(n,2:end));

                        % Compute acceleration on the next step
                        a(n,2:end) = model.W(n) *  xi_n * v(n,2:end) +...
                                   model.W(n)^2 * (-ox2 * d(n,2:end) - y * ox2^0.5 * xi_n * exp(-xwt).* c);
                               
                    else % is critically damped

                        % Compute auxiliar parameters
                        wt = model.W(n) * dt;
                        
                        % Compute displacement on the next step
                        d(n,2:end) = (d0 + (v0 + model.W(n)*d0)*dt).* exp(-wt);
                        
                        % Compute speed on the next step
                        v(n,2:end) = (v0 + model.W(n)*d0) * exp(-wt) - model.W(n) * d(n,2:end);
                        
                        % Compute acceleration on the next step
                        a(n,2:end) = model.W(n) * ((v0 + model.W(n)*d0) * exp(-wt) - v(n,2:end));
                    end
                end
                
           end
            
            % Patricular solution on the ODE (external loads) - Newmark
            %--------------------------------------------------------------
            % Get modal load matrix
            Ffd = model.V' * Ff;

            % Initialize displacement, velocity and acceleration matrices
            dp = zeros(model.n_modes , model.n_steps + 1);
            vp = zeros(model.n_modes , model.n_steps + 1);
            ap = zeros(model.n_modes , model.n_steps + 1);
            
            % Set initial conditions for forced + free vibration
            dp(:,1) = d(:,1);
            vp(:,1) = v(:,1);
            ap(:,1) = (Ffd(:,1) - Cd.*v0 - Kd.*d0) ./ Md;

            % Compute auxiliar matrices
            Aux_B = (2/step) * Cd + (2/step)^2 * Md;
            Aux_A = Kd + Aux_B;
            Aux_C = Cd + (4/step) * Md;

            % Initialize counter
            count = 1;

            % Loop through all steps
            for i = step:step:model.t
                % Update counter
                count = count + 1;

                % Compute displacement on the next step
                u = diag(Aux_A) \ (Aux_B.* dp(:,count-1) +...
                                   Aux_C.* vp(:,count-1) +...
                                   Md.* ap(:,count-1) + Ffd(:,count));

                % Compute speed on the next step
                w = (2/step) * (u - dp(:,count-1)) - vp(:,count-1);

                % Compute acceleration on the next step
                z = (2/step)^2 * (u - dp(:,count-1)) -...
                    (4/step)   * vp(:,count-1) - ap(:,count-1);

                % Store computed components to matrices
                dp(:,count) = u;
                vp(:,count) = w;
                ap(:,count) = z;
            end
            
            %--------------------------------------------------------------
            
            % Initialize 3D matrices to store results
            dynamicDispl = zeros(model.neq,model.n_steps + 1,model.n_modes);
            dynamicVeloc = zeros(model.neq,model.n_steps + 1,model.n_modes);
            dynamicAccel = zeros(model.neq,model.n_steps + 1,model.n_modes);
            dynamicDisplForced = dynamicDispl;
            dynamicVelocForced = dynamicVeloc;
            dynamicAccelForced = dynamicAccel;
            
            % Loop through eigenvectors matrix rows and results columns
            % (time steps) to compute contribution of each vibration mode
            % on each DOF on each instant
            for i = 1:nf
                for j = 1:model.n_steps + 1
                    dynamicDispl(i,j,:) = model.V(i,:).*(dp(:,j))';
                    dynamicVeloc(i,j,:) = model.V(i,:).*(vp(:,j))';
                    dynamicAccel(i,j,:) = model.V(i,:).*(ap(:,j))';
                    dynamicDisplForced(i,j,:) = model.V(i,:).*(d(:,j))';
                    dynamicDisplForced(i,j,:) = dynamicDispl(i,j,:) - dynamicDisplForced(i,j,:);
                    dynamicVelocForced(i,j,:) = model.V(i,:).*(v(:,j))';
                    dynamicVelocForced(i,j,:) = dynamicVeloc(i,j,:) - dynamicVelocForced(i,j,:);
                    dynamicAccelForced(i,j,:) = model.V(i,:).*(a(:,j))';
                    dynamicAccelForced(i,j,:) = dynamicAccel(i,j,:) - dynamicAccelForced(i,j,:);
                end
            end
            
            % Store results
            model.results.dynamicDispl       =       dynamicDispl;
            model.results.dynamicVeloc       =       dynamicVeloc;
            model.results.dynamicAccel       =       dynamicAccel;
            model.results.dynamicDisplForced = dynamicDisplForced;
            model.results.dynamicVelocForced = dynamicVelocForced;
            model.results.dynamicAccelForced = dynamicAccelForced;
            
            % Return success flag (will always be the case, unless a bug
            % occurs previously).
            status = true;
            
            % #########################################################################
            %                 ATTENTION
            %                 The following block of code consisted on an attempt to
            %                 implement a particular analytical solution of the
            %                 uncoupled modal system. This was verified to be of a high
            %                 level of complexity, due to the fact that multiple
            %                 loading frequencies are possible.
            %                 Once acknowledged that this process would lead to
            %                 numerical solutions anyway, a decision was made to solve
            %                 this as it is done by the numerical solver, using
            %                 Newmark's mehtod, but considering the uncoupled modal 
            %                 system.
            %                 In the future, this may be changed.
            %
            %                 % Get load oscilation
            %                 swt  = model.V(:,n)' * sin(Ff(:,2) * dt + Ff(:,3));
            %                 
            %                 % Compute load over time vector
            %                 ft = Ampl * swt;
            %                 
            %                 % Get load frequency
            %                 wl  = model.F(n,2);
            %                 wlt =       w*dt;
            %                 
            %                 % Get load phase angle
            %                 phL = model.F(n,3);
            %                 
            %                 % Compute load to natural frequency ratio
            %                 beta = wl / model.W(n);
            %                 
            %                 % Compute auxiliary parameters
            %                 aux_1 = (1 - beta^2);
            %                 aux_2 = (Cd(n)/Md(n)) *  wl * model.W(n);
            %                 
            %                 % Compute dynamic amplification factor
            %                 DAF = 1/(aux_1^2 + aux_2^2)^0.5;
            %                 
            %                 % Compute displacement amplitude for external loading on
            %                 % this vibration mode
            %                 yL = ft * DAF / Kd(n);
            %                 
            %                 % Compute displacement due to external dynamic loads
            %                 dp(2:end) = yL * (aux_1 * sin(wlt + phL) - aux_2 * cos(wlt + phL));
            %                 vp(2:end) = yL * wl * (aux_1 * cos(wlt + phL) + aux_2 * sin(wlt + phL));
            %                 ap(2:end) = -(wl^2) * dp(2:end);
            % #########################################################################
        end
    end
end