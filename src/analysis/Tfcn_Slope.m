%% Tfcn_Slope (Slope time function) Class
%
%% Description
%
% This is a sub-class of the <tfcn.html *Tfcn*> class for the
% implementation of the *Slope* time function:
%
% FF(t)  = (t - ti) / (tf - ti)
%
classdef Tfcn_Slope < Tfcn    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function fcn = Tfcn_Slope(id,fctr,prev,ti,tf,t,ns,evalFlag)
            include_constants;
            fcn = fcn@Tfcn(SLOPE,id,fctr,ti,tf,t,ns,prev);
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
            
            % Compute time vector
            t = linspace(0,(nf-ni)*step,nf-ni+1);
            
            % Set normalized values to vector
            fcn.value(ni:nf) = fcn.weightFctr * t / ((nf-ni)*step);
        end
    end
end