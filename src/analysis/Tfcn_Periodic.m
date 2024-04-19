%% Tfcn_Periodic (Periodic time function) Class
%
%% Description
%
% This is a sub-class of the <tfcn.html *Tfcn*> class for the
% implementation of the *Periodic* time function:
%
% F(t) = cos(wt + phi)
%
classdef Tfcn_Periodic < Tfcn
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        w   = 0;   % angular speed [rad/s]
        phi = 0;   % phase angle [rad]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function fcn = Tfcn_Periodic(id,fctr,prev,w,phi,ti,tf,t,ns,evalFlag)
            include_constants;
            fcn = fcn@Tfcn(PERIODIC,id,fctr,ti,tf,t,ns,prev);
            fcn.w   =    w;
            fcn.phi =  phi;
            if (nargin <= 9)
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
            fcn.value(ni:nf) = fcn.weightFctr * cos(fcn.w * t + fcn.phi);
        end
    end
end