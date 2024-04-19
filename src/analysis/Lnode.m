%% Lnode (Nodal Load) Class
%
%% Description
%
% This is a handle class for the definition of a nodal load.
%
% Nodal loads are external forcing conditions considered to be singularly 
% applied on specified DOFs.
%
classdef Lnode < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        node    = [];  % handle to respective node object
        static  = [];  % vector of applied static load components [fx fy fz mx my mz]
        dynamic = [];  % vector of amplitudes for applied dynamic load components [fx fy fz mx my mz]
    end
    
    %% Protected attributes
    properties (SetAccess = private, GetAccess = protected)
        fcn = [];  % handles to Lfcn objects
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function load = Lnode(node,fcn)
            if (nargin > 0)
                load.node = node;
                if (nargin > 1)
                    load.fcn = fcn;
                end
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Sets provided time function to this load.
        function setFcn(load,fcn)
            load.fcn = fcn;
        end
        
        %------------------------------------------------------------------
        % Gets time function.
        function ptr = getFcn(load)
            if isempty(load.fcn)
                ptr=[];
                return
            end
            if ~isvalid(load.fcn)
                load.fcn = [];
            end
            ptr = load.fcn;
        end
        
        %------------------------------------------------------------------
        % Sets time function associated to this load to another load.
        % Any changes to the time functions will affect both loads.
        function setFcn2(load,otherLoad)
            if isempty(load.fcn)
                otherLoad.fcn=[];
                return
            end
            if ~isvalid(load.fcn)
                load.fcn = [];
            end
            
            otherLoad.fcn = load.fcn;
        end
        
        %------------------------------------------------------------------
        % Evaluates dynamic load values throught analysis duration
        % Output:
        %    f: load values
        function f = evalDynamicLoad(load)
            
            if isempty(load.fcn)
                f = zeros(length(load.dynamic),1);
                return
            end
            
            if ~isvalid(load.fcn)
                load.fcn = [];
                f = zeros(length(load.dynamic),1);
                return
            end
            
            f0 = load.fcn.evalAll();
            if ~f0
                f = zeros(length(load.dynamic),1);
                return
            end
            
            f_aux = zeros(length(load.dynamic),length(f0));
            for i = 1:length(load.dynamic)
                f_aux(i,:) = f0;
            end
            f = diag(load.dynamic) * f_aux;
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Lnode object.
        function clean(load)
            load.node    = [];
            load.fcn     = [];
            load.static  = [];
            load.dynamic = [];
        end
    end
end