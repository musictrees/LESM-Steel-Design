%% Tfcn_Static (Static time function) Class
%
%% Description
%
% This is a sub-class of the <tfcn.html *Tfcn*> class for the
% implementation of the *Static (Uniform)* time function.
%
% The idea of static time functions is that, although their value will stay
% always the same, they only occur on a determined time interval.
%
classdef Tfcn_Static < Tfcn
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function fcn = Tfcn_Static(id,fctr,prev,ti,tf,t,ns,evalFlag)
            include_constants;
            fcn = fcn@Tfcn(STATIC,id,fctr,ti,tf,t,ns,prev);
            if (nargin <= 7)
                fcn.evalFcn();
            elseif evalFlag
                fcn.evalFcn();
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <tfcn.html *Tfcn*>.
    methods
        %------------------------------------------------------------------
        % Computes normalized load values throughout analysis interval
        % Input arguments:
        %  fcn: handle to this object
        function evalFcn(fcn)
            % Initialize vector to be returned with zeros
            fcn.value = zeros(fcn.nsteps+1,1);
            
            % Compute time steps
            step = fcn.totalT/fcn.nsteps;
            
            % Compute initial and final time step where function acts
            ni =  ceil(fcn.ti / step) + 1;
            if ni > fcn.nsteps + 1
                return
            end
            nf = floor(fcn.tf / step) + 1;
            if nf > fcn.nsteps + 1
                nf = fcn.nsteps + 1;
            end
            if ni == nf
                return
            end
            
            % Set normalized values to vector
            fcn.value(ni:nf) = fcn.weightFctr * ones(nf-ni+1,1);
        end
    end
end