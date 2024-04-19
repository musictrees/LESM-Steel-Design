%% Tfcn_Table (Table time function) Class
%
%% Description
%
% This is a sub-class of the <tfcn.html *Tfcn*> class for the
% implementation of the *Table* time function.
%
classdef Tfcn_Table < Tfcn
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        x        = [];   % function domain - time [s]
        F        = [];   % function values - Amplitude []
        src_file = '';   % table sorce file name
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function fcn = Tfcn_Table(id,fctr,prev,x,F,t,ns,evalFlag,srcf)
            include_constants;
            ti = x(1);
            tf = x(end);
            if x(1) < 0
                ti = 0;                
            end
            if x(end) > t
                tf = t;
            end
            fcn = fcn@Tfcn(TABLE,id,fctr,ti,tf,t,ns,prev);
            fcn.x =  x;
            fcn.F =  F;
            if (nargin <= 8)
                fcn.evalFcn();
            elseif evalFlag
                fcn.evalFcn();
            end
            if nargin == 9
                fcn.src_file = srcf;
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
            t = linspace(0,(nf-ni)*step,nf-ni+1)+fcn.x(1);
                        
            % Set normalized values to vector
            fcn.value(ni:nf) = fcn.weightFctr * interp1(fcn.x,fcn.F,t);
            
            % For safaty - Avoid extrapolation
            fcn.value(isnan(fcn.value)) = 0;
        end
    end
end