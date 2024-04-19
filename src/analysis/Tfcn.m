%% Tfcn (Time function) Class
%
%% Description
%
% This is a handle super-class for the definition of a time function.
%
% Time functions are considered to describe the behaviour of load values
% throughout the analysis interval.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <tfcn_static.html static (uniform) time function>.
% * <tfcn_slope.html slope time function>.
% * <tfcn_periodic.html periodic time function>.
% * <tfcn_table.html table time function>.
%
classdef Tfcn < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id         = [];   % load function identifier
        next       = [];   % handle to next Tfcn object
        prev       = [];   % handle to prev Tfcn object
        type       =  0;   % flag for load type
        weightFctr =  1;   % weightFactor (intensity) of normalized function
        nsteps     =  0;   % number of steps for analysis
        ti         =  0;   % initial instant of load action
        tf         =  0;   % final instant of load action
        totalT     =  0;   % total analysis interval
        value      = [];   % normalized function values
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function fcn = Tfcn(type,id,fctr,ti,tf,t,ns,prev)
            if (nargin > 0)
                fcn.type       = type;
                fcn.id         =   id;
                fcn.weightFctr = fctr;
                fcn.ti         =   ti;
                fcn.tf         =   tf;
                fcn.totalT     =    t;
                fcn.nsteps     =   ns;
                fcn.next       =   [];
                fcn.prev       = prev;
                if ~isempty(prev), fcn.prev.next = fcn; end
            end
        end
    end
    
    %% Private methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Computes normalized load values throughout analysis interval
        % Input arguments:
        %  fcn: handle to this object
        evalFcn(fcn)
    end
    
    %% Public methods
    methods
        % Evaluates normalized dynamic load values throught analysis
        % Output:
        %    f: normalized load values
        function f = evalAll(fcn,mustEvalFcns)
            if nargin == 1
                mustEvalFcns = false;
            end
            
            ptr = fcn;
            
            % Initialize vector to be returned with zeros
            f = zeros(ptr.nsteps+1,1);
            
            % Add contribution of each fcn
            if mustEvalFcns
                while ~isempty(ptr)
                    ptr.evalFcn();
                    f = f + ptr.value;
                    ptr = ptr.next;
                end
            else
                while ~isempty(ptr)
                    f = f + ptr.value;
                    ptr = ptr.next;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Goes through list of fcn handles. Returns the end of the list and
        % the number of handles on the list.
        function [ptr,nFcn] = goThrough(fcn)
            ptr = fcn;
            nFcn = 1;
            
            while ~isempty(ptr.next)
                ptr = ptr.next;
                nFcn = nFcn + 1;
            end
        end
        
        %------------------------------------------------------------------
        % Goes through list of fcn handles. Returns a vector of identifiers
        % for each fcn on the list
        function id = getIds(fcn)
            ptr = fcn; 
            
            [~,nFcn] = fcn.goThrough();
            id = zeros(1,nFcn);
            count = 1;
            
            while ~isempty(ptr)
                id(count) = ptr.id;
                ptr = ptr.next;
                count = count + 1;
            end
        end
        
        %------------------------------------------------------------------
        % Goes through list of fcn handles searching for fcn with given id.
        % Returns fcn handle or empty handle (in case fcn id is not found)
        function ptr = getById(fcn,id)
            ptr = fcn;
            
            while ~isempty(ptr)
                if length(ptr.id) == length(id)
                    if (ischar(ptr.id) && ischar(id)) ||...
                       (~ischar(ptr.id) && ~ischar(id))
                        if all(ptr.id == id)
                            return
                        end
                    end
                end
                ptr = ptr.next;
            end
        end
        
        %------------------------------------------------------------------
        % Creates a copy of list of time functions
        % ATTENTION: The resulting copy will have the same properties as
        % the original one, but is a different list. That is, afterwards,
        % any changes made to the original will not affect the copy.
        function copy = createCopy(fcn)
            include_constants;
            ptr = fcn;
            
            switch ptr.type
                case STATIC
                    % Create Tfcn object
                    newFcn = Tfcn_Static(ptr.id,ptr.weightFctr,[],ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                    newFcn.value = ptr.value;
                    
                case PERIODIC
                    % Create Tfcn object
                    newFcn = Tfcn_Periodic(ptr.id,ptr.weightFctr,[],ptr.w,ptr.phi,ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                    newFcn.value = ptr.value;
                    
                case SLOPE
                    % Create Tfcn object
                    newFcn = Tfcn_Slope(ptr.id,ptr.weightFctr,[],ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                    newFcn.value = ptr.value;
                    
                case TABLE
                    % Create Tfcn object
                    newFcn = Tfcn_Table(ptr.id,ptr.weightFctr,[],ptr.x,ptr.F,ptr.totalT,ptr.nsteps,false);
                    newFcn.value = ptr.value;
            end
            
            copy = newFcn;
            prev_copy = copy;
            ptr = ptr.next;
            
            while ~isempty(ptr)
                switch ptr.type
                    case STATIC
                        % Create Tfcn object
                        newFcn = Tfcn_Static(ptr.id,ptr.weightFctr,prev_copy,ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                        newFcn.value = ptr.value;

                    case PERIODIC
                        % Create Tfcn object
                        newFcn = Tfcn_Periodic(ptr.id,ptr.weightFctr,prev_copy,ptr.w,ptr.phi,ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                        newFcn.value = ptr.value;

                    case SLOPE
                        % Create Tfcn object
                        newFcn = Tfcn_Slope(ptr.id,ptr.weightFctr,prev_copy,ptr.ti,ptr.tf,ptr.totalT,ptr.nsteps,false);
                        newFcn.value = ptr.value;
                        
                    case TABLE
                        % Create Tfcn object
                        newFcn = Tfcn_Table(ptr.id,ptr.weightFctr,prev_copy,ptr.x,ptr.F,ptr.totalT,ptr.nsteps,false);
                        newFcn.value = ptr.value;
                end
                prev_copy = newFcn;
                ptr = ptr.next;
            end
        end
        
        %------------------------------------------------------------------
        % Sets new number of steps for Tfcn objects and reevaluate
        % normalized functions.
        function update_nsteps(fcn,new_nsteps)
            ptr = fcn;
            
            while ~isempty(ptr)
                ptr.nsteps = new_nsteps;
                ptr.evalFcn();
                ptr = ptr.next;
            end
        end
        %------------------------------------------------------------------
        % Sets new total time value for Tfcn objects and reevaluate
        % normalized functions.
        function update_time(fcn,new_time)
            ptr = fcn;
            
            while ~isempty(ptr)
                ptr.totalT = new_time;
                ptr.evalFcn();
                ptr = ptr.next;
            end
        end
         %------------------------------------------------------------------
        % Sets new total time value and number of steps for Tfcn objects
        % and reevaluate normalized functions.
        function update_time_nsteps(fcn,new_time,new_nsteps)
            ptr = fcn;
            
            while ~isempty(ptr)
                ptr.nsteps = new_nsteps;
                ptr.totalT = new_time;
                ptr.evalFcn();
                ptr = ptr.next;
            end
        end  
        
        %------------------------------------------------------------------
        % Cleans data structure of a Tfcn object.
        function clean(fcn)
            fcn.id         = [];
            fcn.next       = [];
            fcn.prev       = [];
            fcn.type       = [];
            fcn.weightFctr = [];
            fcn.nsteps     = [];
            fcn.ti         = [];
            fcn.tf         = [];
            fcn.totalT     = [];
            fcn.value      = [];
        end  
    end
end