%% Model Class
%
%% Description
%
% This is a handle class for the definition of a structural model.
%
% A model object is responsible for storing the global variables of a
% structural analysis problem.
%
classdef Model < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        drv            = [];      % handle to an object of the Drv class
        anm            = [];      % handle to an object of the Anm class
        materials      = [];      % vector of handles to objects of the Material class
        sections       = [];      % vector of handles to objects of the Section class
        nodes          = [];      % vector of handles to objects of the Node class
        elems          = [];      % vector of handles to objects of the Elem class
        srjoints       = [];      % vector of handles to objects of the Srjoint class
        results        = [];      % handle to an object of the Result class
        K              = [];      % global stiffness matrix
        M              = [];      % global mass matrix
        C              = [];      % global damping matrix
        F              = [];      % global forcing vector
        D              = [];      % global displacement vector
        W              = [];      % array of vibration frequencies
        V              = [];      % eigenvectors vibration modules matrix
        damping        = 0;       % flag for type of damping to be considered
        xi             = [0,0];   % critical damping coefficients
        massDampCoeff  = 0;       % mass proportional damping coefficient
        stiffDampCoeff = 0;       % stiffness proportional damping coefficient
        c0             = [];      % initial condition matrix for dynamic analysis [d0, v0]
        n_steps        = 0;       % number of steps for numerical integration method
        t              = 0;       % time interval for numerical integration method
        n_modes        = 0;       % max number of vibration modes to be computed
        whichSolver    = 0;       % flag for solver to be used when processing data
        mass_type      = 1;       % flag for type of element local mass matrix to be considered (Consistent, Lumped, Mixed)
        mass_mi        = 1;       % mixed mass matrix proportion coefficient
        ID             = [];      % global d.o.f. numbering matrix
        nmat           = 0;       % number of materials
        nsec           = 0;       % number of cross-sections
        nnp            = 0;       % number of nodes
        nel            = 0;       % number of elements
        njoints        = 0;       % number of semi-rigid joint nodes
        neq            = 0;       % number of equations
        neqfree        = 0;       % number of equations of free d.o.f.
        neqfixed       = 0;       % number of equations of fixed d.o.f.
        neqspring      = 0;       % number of equations of spring d.o.f.
        gblSprStiff    = [];      % global spring stiffness coefficients vector
        nlc            = 0;       % number of load cases
        ncomb          = 0;       % number of load combinations
        strLc          = {};      % cell array of load case names
        strComb        = {};      % cell array of load case combination names
        loadComb       = [];      % matrix of combination factors for each load case combination
        loadCombID     = [];      % matrix of flags for load cases in each combination
        timeFcns       = {};      % time fucntions list
        strTimeFcns    = {};      % cell array of time function names
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function model = Model(drv,anm,mat,sec,nodes,elems,strlc,strcomb,loadcomb)
            if (nargin > 0)
                model.drv        =             drv;
                model.anm        =             anm;
                model.materials  =             mat;
                model.sections   =             sec;
                model.nodes      =           nodes;
                model.elems      =           elems;
                model.nmat       =     length(mat);
                model.nsec       =     length(sec);
                model.nnp        =   length(nodes);
                model.nel        =   length(elems);
                
                % Store load cases and combinations
                if nargin >= 6
                    model.nlc = length(strlc);
                    model.strLc = strlc;
                    if nargin >= 7
                        model.ncomb = length(strcomb);
                        model.strComb = strcomb;
                        model.loadComb = loadcomb;
                        model.loadCombID = loadcomb;
                        for i = 1:size(loadcomb,1)
                            for j = 1:size(loadcomb,2)
                                if model.loadCombID(i,j) ~= 0
                                    model.loadCombID(i,j) = 1;
                                end
                            end
                        end
                    end
                else
                    model.nlc = 1;
                    model.strLc = {'CASE 01'};
                end
            end
        end
    end
    
    %% Private methods 
    methods(Access=private)
        function flag = checkFcnId(model,fcnId)
            flag = true;
            if fcnId < 1 || fcnId > length(model.timeFcns)
                flag = false;
                return
            end
            if isempty(model.timeFcns{fcnId})
                flag = false;
                return
            end
            if ~isvalid(model.timeFcns{fcnId})
                model.timeFcns{fcnId} = [];
                flag = false;
                return
            end
        end
    end
    
    %% Public methods
     methods
        %------------------------------------------------------------------
        % Adds a new fcn list (a load case) to the cell arrays
        function fcnId = createFcnList(model,name)
            fcnId = length(model.timeFcns)+1;
            model.timeFcns{fcnId} = [];
            model.strTimeFcns{fcnId} = name;
        end
        
        %------------------------------------------------------------------
        % deletes a fcn list (a load case) from the cell arrays
        function deleteFcnList(model,fcnId)
            if fcnId > length(model.timeFcns)
                return
            end
            model.removeFcn(fcnId,1); %removing top erases whole list
            model.timeFcns(fcnId) = [];
            model.strTimeFcns(fcnId) = [];
        end
        
        %------------------------------------------------------------------
        % Finds a fcn listid from a handle (pointer)
        function fcnId = findFcnListByHandle(model,ptr)
            fcnId = 0;
            if isempty(ptr) || isempty(model.timeFcns), return; end
            for ii = 1:length(model.timeFcns)
                if ptr == model.timeFcns{ii}
                    fcnId = ii;
                    return
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add a new handle to a Tfcn object to fcn list
        function addFcn(model,fcnId,type,fcnInput)
            include_constants;                    
            
            if ~model.checkFcnId(fcnId)
                prev = []; 
            else
               [prev,~] = model.timeFcns{fcnId}.goThrough(); 
            end
            
            switch type
                case STATIC
                    % Get function input from cell
                    id   = fcnInput{1};
                    fctr = fcnInput{2};
                    ti   = fcnInput{3};
                    tf   = fcnInput{4};
                                        
                    % Create Tfcn object
                    newFcn = Tfcn_Static(id,fctr,prev,ti,tf,model.t,model.n_steps);
                    
                case PERIODIC
                    % Get function input from cell
                    id   = fcnInput{1};
                    fctr = fcnInput{2};
                    w    = fcnInput{3};
                    phi  = fcnInput{4};
                    ti   = fcnInput{5};
                    tf   = fcnInput{6};
                                      
                    % Create load function object
                    newFcn = Tfcn_Periodic(id,fctr,prev,w,phi,ti,tf,model.t,model.n_steps);
                    
                case SLOPE
                    % Get function input from cell
                    id   = fcnInput{1};
                    fctr = fcnInput{2};
                    ti   = fcnInput{3};
                    tf   = fcnInput{4};
                                       
                    newFcn = Tfcn_Slope(id,fctr,prev,ti,tf,model.t,model.n_steps);
                    
                case TABLE 
                    % Get input from cell
                    id   = fcnInput{1};
                    fctr = fcnInput{2};
                    x    = fcnInput{3};
                    fval = fcnInput{4};
                                        
                    newFcn = Tfcn_Table(id,fctr,prev,x,fval,model.t,model.n_steps);
            end
            
            % If this is the first handle to fcn on the list,
            % store as a property of this Lnode object
            if isempty(prev)
                model.timeFcns{fcnId} = newFcn;
            end
        end
        
        %------------------------------------------------------------------
        % Remove handle to a Tfcn object from fcn list
        function removeFcn(model,fcnId,id)
            
            if ~model.checkFcnId(fcnId)
                return
            end                
            
            fcn = model.timeFcns{fcnId};
            
            % Get handle to be removed from list
            toBeRemoved = fcn.getById(id);
            
            % If fcn with given id was not found, do nothing
            if isempty(toBeRemoved)
                return
            end
            
            % Check if fcn to be removed is top of the list
            if fcn == toBeRemoved
            % Destroy whole list of time functions to aviod errors
                ptr=toBeRemoved.next;
                while (~isempty(ptr))
                    ptr.prev.next = ptr.next;
                    if (~isempty(ptr.next))
                    ptr.next.prev = ptr.prev;
                    end
                    temp = ptr;                    
                    ptr = ptr.next;
                    delete(temp);
                end    
                delete(toBeRemoved);
                % leave empty fcn, but does not exclude load case from cell
                model.timeFcns{fcnId} = [];
                return
            end
            
            % If reached this point, simply remove fcn from list
            toBeRemoved.prev.next = toBeRemoved.next;
            if ~isempty(toBeRemoved.next)
                toBeRemoved.next.prev = toBeRemoved.prev;
            end
            delete(toBeRemoved);
        end
        
        %------------------------------------------------------------------
        % Creates or removes fictitious rotation constraints on nodes where
        % all incident elements are hinged, to force global stability
        % during the analysis process by avoiding a singular stifness matrix.
        % Input arguments:
        %  fict: flag to check if constraints must be created or removed
        function fictRotConstraint(model,fict)
            include_constants;
            for n = 1:model.nnp
                
                % Create fictitious rotation constraints (before analysis process)
                if fict == 1
                    % Calculate total number of incident elements and hinged
                    % elements on current node
                    [tot,hng] = model.nodes(n).elemsIncidence(model);
                    
                    % Check if all incident elements are hinged
                    if tot == hng
                         % Insert ficticious constraint (flag -1)
                         % on free rotation d.o.f.'s. (flag 0)
                        if model.nodes(n).ebc(4) == FREE_DOF
                            model.nodes(n).ebc(4) = FICTFIXED_DOF;
                        end
                        if model.nodes(n).ebc(5) == FREE_DOF
                            model.nodes(n).ebc(5) = FICTFIXED_DOF;
                        end
                        if model.nodes(n).ebc(6) == FREE_DOF
                            model.nodes(n).ebc(6) = FICTFIXED_DOF;
                        end
                    end
                    
                % Remove fictitious rotation constraints (after analysis process)
                elseif fict == 0
                    if model.nodes(n).ebc(4) == FICTFIXED_DOF
                        model.nodes(n).ebc(4) = FREE_DOF;
                    end
                    if model.nodes(n).ebc(5) == FICTFIXED_DOF
                        model.nodes(n).ebc(5) = FREE_DOF;
                    end
                    if model.nodes(n).ebc(6) == FICTFIXED_DOF
                        model.nodes(n).ebc(6) = FREE_DOF;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Computes number of equations to be processed and updtes neq
        % property.
        function computeNeq(model)
            % Compute total number of semi-rigid joints
            model.computeNumSrjoints();
            
            % Set "neq" property of the model object
            model.neq = (model.nnp * model.anm.ndof) + (model.njoints * model.anm.nrdof);
        end
        
        %------------------------------------------------------------------
        % Computes the number of semi-rigid joints of the entire model.
        % It traverses the list of elements and, for each one, checks
        % whether it has semi-rigid joint at the element extremities. 
        function computeNumSrjoints(model)
            include_constants;
            model.njoints = 0;
            for e = 1:model.nel
                if (model.elems(e).hingei == SEMIRIGID_END)
                    model.njoints = model.njoints + 1;
                end
                if (model.elems(e).hingef == SEMIRIGID_END)
                    model.njoints = model.njoints + 1;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles global d.o.f (degree of freedom) numbering ID matrix:
        %  ID(k,n) = equation number of d.o.f. k of node n
        %  Free d.o.f.'s have the initial numbering.
        %  countF --> counts free d.o.f.'s. (numbered first)
        %  countS --> counts spring d.o.f.'s (numbered after free d.o.f.'s)
        %  countC --> counts fixed d.o.f.'s. (numbered last)
        function assembleDOFNum(model)
            include_constants;
            
            % Initialize equation numbers for free, fixed and spring d.o.f.'s.
            countF = 0;
            countS = model.neqfree;
            countC = model.neqfree + model.neqspring;
            
            % Check if each d.o.f. is free, fixed or spring to increment eqn.
            % number and store it in the ID matrix.
            for n = 1:model.nnp
                for k = 1:model.anm.ndof
                    if model.ID(k,n) == FREE_DOF
                        countF = countF + 1;
                        model.ID(k,n) = countF;
                    elseif model.ID(k,n) == SPRING_DOF
                        countS = countS + 1;
                        model.ID(k,n) = countS;
                    else
                        countC = countC + 1;
                        model.ID(k,n) = countC;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Creates the array of semi-rigid joints of the entire model.
        % It traverses the list of elements and, for each one, checks
        % whether it has semi-rigid joint at the element extremities.
        % It also numbers the eqs. of free d.o.f.'s associated to 
        % semi-rigid joints.
        function createSrjoints(model)
            include_constants;
            if model.njoints > 0
                joints(1,model.njoints) = Srjoint();
                nj = 0;   % Index to the current semi-rigid joint
                % Holds eq. number of joint rot. d.o.f.
                srjeq = model.neqfree - (model.njoints * model.anm.nrdof);
                for e = 1:model.nel
                    if (model.elems(e).hingei == SEMIRIGID_END)
                        nj = nj + 1;
                        joints(nj).setElemToSrj(model.elems(e));
                        joints(nj).krx = joints(nj).elem.kri(1);
                        joints(nj).kry = joints(nj).elem.kri(2);
                        joints(nj).krz = joints(nj).elem.kri(3);
                        joints(nj).elem.srjointi = joints(nj);
                        joints(nj).node = joints(nj).elem.nodes(1);
                        for i = 1:joints(nj).anm.nrdof
                            srjeq = srjeq + 1;
                            joints(nj).eqs(i) = srjeq;
                        end
                    end
                    if (model.elems(e).hingef == SEMIRIGID_END)
                        nj = nj + 1;
                        joints(nj).setElemToSrj(model.elems(e));
                        joints(nj).krx = joints(nj).elem.krf(1);
                        joints(nj).kry = joints(nj).elem.krf(2);
                        joints(nj).krz = joints(nj).elem.krf(3);
                        joints(nj).elem.srjointf = joints(nj);
                        joints(nj).node = joints(nj).elem.nodes(2);
                        for i = 1:joints(nj).anm.nrdof
                            srjeq = srjeq + 1;
                            joints(nj).eqs(i) = srjeq;
                        end
                    end
                end
                % Store semi-rigid objects on the model object
                model.srjoints = joints;
            end
        end
        
        %------------------------------------------------------------------
        % Assembles gather vectors of elements (gle) and semi-rigid joints
        % (glj) that stores element and joint d.o.f.'s equation numbers.
        function assembleGatherVectors(model)
            for e = 1:model.nel
                model.anm.assembleGle(model,model.elems(e));
            end
            
            for j = 1:model.njoints
                model.anm.assembleGlj(model,model.srjoints(j));
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element stiffness matrix coefficients (in global system)
        % to the correct positions of the global stiffness matrix.
        % Input arguments:
        %  keg: element stiffness matrix in global system
        %  e: element identification number
        function assembleElemStiffMtx(model,keg,e)
            % Extract element gather vector
            gle = model.elems(e).gle;
            
            % Add contribution of element stifness matrix to
            % global stiffness matrix.
            model.K(gle,gle) = model.K(gle,gle) + keg;
        end
        
        %------------------------------------------------------------------
        % Assembles element mass matrix coefficients (in global system)
        % to the correct positions of the global mass matrix.
        % Input arguments:
        %  meg: element mass matrix in global system
        %  e: element identification number
        function assembleElemMassMtx(model,meg,e)
            % Extract element gather vector
            gle = model.elems(e).gle;
            
            % Add contribution of element mass matrix to
            % global mass matrix.
            model.M(gle,gle) = model.M(gle,gle) + meg;
        end
        
        %------------------------------------------------------------------
        % Assembles semi-ridig joint stiffness matrix coefficients
        % (in global system) to the correct positions of the global
        % stiffness matrix.
        % Input arguments:
        %  kjg: semi-ridig joint stiffness matrix in global system
        %  j: semi-ridig joint identification number
        function assembleSrjointStiffMtx(model,kjg,j)
            % Extract semi-rigid joint gather vector
            glj = model.srjoints(j).glj;
            
            % Add contribution of semi-rigid joint stifness matrix to
            % global stiffness matrix.
            model.K(glj,glj) = model.K(glj,glj) + kjg;
        end
        
        %------------------------------------------------------------------
        % Assembles semi-ridig joint mass matrix coefficients
        % (in global system) to the correct positions of the global
        % mass matrix.
        % Input arguments:
        %  j: semi-ridig joint identification number
        function assembleSrjointMassMtx(model,j) %#ok<INUSD>
%             % Extract semi-rigid joint gather vector
%             glj = model.srjoints(j).glj;
%             
%             % Add contribution of semi-rigid joint mass matrix to
%             % global mass matrix.
%             model.M(glj,glj) = model.M(glj,glj) + mjg;
        end
        
        %------------------------------------------------------------------
        % Adds element equivalent nodal loads (ENL), from distributed loads
        % and thermal loads, to global forcing vector.
        function elemLoads(model)
            for e = 1:model.nel
                % Add distributed load as ENL to global forcing vector
                feg = model.elems(e).load.gblDistribLoadENL();
                model.assembleENL(feg,e);
                
                % Add thermal load as ENL to global forcing vector
                feg = model.elems(e).load.gblThermalLoadENL();
                model.assembleENL(feg,e);
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element equivalent nodal load vector (in global system)
        % to any term of the global forcing vector, including the terms that
        % correspond to constrained d.o.f.'s.
        % Input arguments:
        %  feg: element equivalent nodal load vector in global system
        %  e: element identification number
        function assembleENL(model,feg,e)
            % Extract element gather vector
            gle = model.elems(e).gle;
            
            % Add contribution of element equivalent nodal load vector to
            % global forcing vector
            model.F(gle) = model.F(gle) + model.elems(e).rot_inclSupp * feg;
        end
        %------------------------------------------------------------------
        %ronald
        function designSolver(model,standard,fid)
            model.standard=standard;
            
         %%put parts of moment xy in group(bending)
         for i=1:length(model.allStructuralGroups)
             M = zeros(model.nlc,50*length(model.allStructuralGroups(i).elems));
             for e = 1:length(model.allStructuralGroups(i).elems)
                 for lc = 1:model.nlc
                     x = linspace(0,1,50);
                     model.assembleLoadCase(lc);
                     [M(lc,((e-1)*50+1):e*50),~] = model.allStructuralGroups(i).elems(e).intBendingMoment_XY(model.allStructuralGroups(i).elems(e).length*x,lc);
                 end
                 
             end
             model.assembleLoadCase(model.current_lc);
             
             model.allStructuralGroups(i).bending_XY=M; %REVISAR
             %alterar RONALD quando trocar a funcão para bending
             %moments
             for mPos=1:1:size(M,1)
                 model.loadsInfomartions(mPos).load=M(mPos,:);
             end
             currentGroupLoadsInfomartionsLoad=transpose(model.loadsInfomartions);
           
            
             model.allStructuralGroups(i).bendingCases_XY=combinationLoads(currentGroupLoadsInfomartionsLoad);
             model.allStructuralGroups(i).bendingCases_XY;
             model.allStructuralGroups(i).getMaxBending();
         end
         %%
            %Shear       
            for i=1:length(model.allStructuralGroups)
                M = zeros(model.nlc,50*length(model.allStructuralGroups(i).elems));
                for e = 1:length(model.allStructuralGroups(i).elems)
                    for lc = 1:model.nlc
                        x = linspace(0,1,50);
                        model.assembleLoadCase(lc);
                        [M(lc,((e-1)*50+1):e*50),~] = model.allStructuralGroups(i).elems(e).intShearForce_XY(model.allStructuralGroups(i).elems(e).length*x,lc);
                    end
                    
                end
                model.assembleLoadCase(model.current_lc);
                
                model.allStructuralGroups(i).shear_XY=M; %REVISAR
                %alterar RONALD quando trocar a funcão para bending
                %moments
                for mPos=1:1:size(M,1)
                    model.loadsInfomartions(mPos).load=M(mPos,:);
                end
                currentGroupLoadsInfomartionsLoad=transpose(model.loadsInfomartions);
                model.allStructuralGroups(i).shearCases_XY=combinationLoads(currentGroupLoadsInfomartionsLoad);
                model.allStructuralGroups(i).getMaxShear();
                
            end
            %%   
            
             %%
            %AXIAL       
            for i=1:length(model.allStructuralGroups)
                M = zeros(model.nlc,50*length(model.allStructuralGroups(i).elems));
                for e = 1:length(model.allStructuralGroups(i).elems)
                    for lc = 1:model.nlc
                        x = linspace(0,1,50);
                        model.assembleLoadCase(lc);
                        [M(lc,((e-1)*50+1):e*50),~] = model.allStructuralGroups(i).elems(e).intAxialForce(model.allStructuralGroups(i).elems(e).length*x,lc);
                    end
                    
                end
                model.assembleLoadCase(model.current_lc);
                
                model.allStructuralGroups(i).axialForce=M; %REVISAR
                %alterar RONALD quando trocar a funcão para bending
                %moments
                for mPos=1:1:size(M,1)
                    model.loadsInfomartions(mPos).load=M(mPos,:);
                end
                currentGroupLoadsInfomartionsLoad=transpose(model.loadsInfomartions);
                model.allStructuralGroups(i).axialForceCases=combinationLoads(currentGroupLoadsInfomartionsLoad);
                
                model.allStructuralGroups(i).separeTractionCases();%function to separe axial(traction) forces
                model.allStructuralGroups(i).separeCompressCases();%function to separe axial(compress) forces
                model.allStructuralGroups(i).getMaxCompress();%function to get max value(compress).
                model.allStructuralGroups(i).getMaxTraction();%function to get max value(traction)
            end
            %% 
             
            setappdata(0,'allStructuralGroups',model.allStructuralGroups);
             for i=1:1:size(model.allStructuralGroups ,2)
                  model.standard.solve(model.allStructuralGroups(i),fid,i);
             end
            
        end
        %------------------------------------------------------------------
        %ronald
        % Feeds lelem and lnode objects with load case (or combination)
        % associated to index lc
        function assembleLoadCase(model,lc)
            % Need to check if lc reffers to load case of combination
            if lc <= model.nlc % load case
                 % Allocate nodal loads and prescribed displacements in each Node object
                 % obs.: nodalLoadCase(:,i) = [ fx
                 %                              fy
                 %                              fz     (nodalLoadCase is a matrix, each
                 %                              mx      column refers to one specific
                 %                              my      load case)
                 %                              mz
                 %                              dx
                 %                              dy
                 %                              dz
                 %                              rx
                 %                              ry
                 %                              rz ]
                 for n = 1:model.nnp
                     if isempty(model.nodes(n).nodalLoadCase) == 0 && lc <= size(model.nodes(n).nodalLoadCase,2)
                         if all(model.nodes(n).nodalLoadCase(1:6,lc) == 0)
                             model.nodes(n).load.static = [];
                         else %if there are any nodal loads, set them to load.static
                             model.nodes(n).load.static = model.nodes(n).nodalLoadCase(1:6,lc);
                         end
                         if size(model.nodes(n).nodalLoadCase,1) <= 6
                             model.nodes(n).prescDispl = [];
                         elseif all(model.nodes(n).nodalLoadCase(7:12,lc) == 0)
                             model.nodes(n).prescDispl = [];
                         else %if there are any prescribed displacements, set them to prescDispl
                             model.nodes(n).prescDispl = model.nodes(n).nodalLoadCase(7:12,lc);
                         end
                     else
                         model.nodes(n).load.static = [];
                         model.nodes(n).prescDispl = [];
                     end
                 end
                 
                 % Allocate element loads (distributed loads and thermal loads) in each
                 % Lelem object
                 % obs.: elemLoadCase(:,i) = [  unifDir
                 %                                qx
                 %                                qy
                 %                                qz        (elemLoadCase is a matrix,
                 %                             linearDir     each column refers to one
                 %                                qx1        specific load case)
                 %                                qy1
                 %                                qz1
                 %                                qx2
                 %                                qy2
                 %                                qz2
                 %                                dtx
                 %                                dty
                 %                                dtz   ]
                 for e = 1:model.nel
                     if isempty(model.elems(e).load.elemLoadCase) == 0 && lc <= size(model.elems(e).load.elemLoadCase,2)
                         % Clear previous uniform loads
                         model.elems(e).load.uniformGbl = [];
                         model.elems(e).load.uniformLcl = [];
                         if all(model.elems(e).load.elemLoadCase(2:4,lc) == 0)
                             model.elems(e).load.uniformDir = 0;
                         else  % if there are uniform loads, set their local and global components
                             model.elems(e).load.uniformDir = model.elems(e).load.elemLoadCase(1,lc);
                             model.elems(e).load.setUnifLoad((model.elems(e).load.elemLoadCase(2:4,lc))',model.elems(e).load.elemLoadCase(1,lc));
                         end
                         % Clear previous linear loads
                         model.elems(e).load.linearGbl = [];
                         model.elems(e).load.linearLcl = [];
                         if size(model.elems(e).load.elemLoadCase,1) <= 5
                             model.elems(e).load.linearDir = 0;
                         elseif all(model.elems(e).load.elemLoadCase(6:11,lc) == 0)
                             model.elems(e).load.linearDir = 0;
                         else  % if there are linear loads, set their local and global components
                             model.elems(e).load.linearDir = model.elems(e).load.elemLoadCase(5,lc);
                             model.elems(e).load.setLinearLoad((model.elems(e).load.elemLoadCase(6:11,lc))',model.elems(e).load.elemLoadCase(5,lc));
                             % Avoids numeric problems with new loads
                             if model.elems(e).load.linearDir == 0
                                 for i = 1:size(model.elems(e).load.linearLcl,1)
                                     if abs(model.elems(e).load.linearLcl(i)) <= 10^-10
                                         model.elems(e).load.linearLcl(i) = 0;
                                     end
                                 end
                             elseif model.elems(e).load.linearDir == 1
                                 for i = 1:size(model.elems(e).load.linearGbl,1)
                                     if abs(model.elems(e).load.linearGbl(i)) <= 10^-10
                                         model.elems(e).load.linearGbl(i) = 0;
                                     end
                                 end
                             end
                         end
                         % Set thermal loads
                         if size(model.elems(e).load.elemLoadCase,1) <= 11
                             model.elems(e).load.tempVar_X = 0;
                             model.elems(e).load.tempVar_Y = 0;
                             model.elems(e).load.tempVar_Z = 0;
                         else
                             model.elems(e).load.tempVar_X = model.elems(e).load.elemLoadCase(12,lc);
                             model.elems(e).load.tempVar_Y = model.elems(e).load.elemLoadCase(13,lc);
                             model.elems(e).load.tempVar_Z = model.elems(e).load.elemLoadCase(14,lc);
                         end
                     else
                         model.elems(e).load.uniformDir = 0;
                         model.elems(e).load.uniformGbl = [];
                         model.elems(e).load.uniformLcl = [];
                         model.elems(e).load.linearDir = 0;
                         model.elems(e).load.linearGbl = [];
                         model.elems(e).load.linearLcl = [];
                         model.elems(e).load.tempVar_X = 0;
                         model.elems(e).load.tempVar_Y = 0;
                         model.elems(e).load.tempVar_Z = 0;
                     end
                 end
                 
            else % lc is a load combination
                % Consider load cases factors for the selected combination and allocate
                % the resulting nodal loads and prescribed displacements in each Node
                % object.
                % obs.: nodalLoadCase(:,i) = [ fx
                %                              fy
                %                              fz     (nodalLoadCase is a matrix, each
                %                              mx      column refers to a specific load
                %                              my      case)
                %                              mz
                %                              dx
                %                              dy
                %                              dz
                %                              rx
                %                              ry
                %                              rz ]
                %
                % obs.2: loadComb is a matrix that holds combination factors for load
                % case combinations, where each column refers to a specific load comb.
                
                for n = 1:model.nnp
                    if isempty(model.nodes(n).nodalLoadCase) == 0
                        allLogic = all(model.nodes(n).nodalLoadCase(1:6,:) == 0);
                        if all(allLogic == 1)
                            model.nodes(n).load.static = [];
                        else  % if there are nodal loads, ,multiply them by their
                            % proper comb factors, sum the results and set them to
                            % load.static.
                            model.nodes(n).load.static(1) = model.nodes(n).nodalLoadCase(1,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).load.static(2) = model.nodes(n).nodalLoadCase(2,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).load.static(3) = model.nodes(n).nodalLoadCase(3,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).load.static(4) = model.nodes(n).nodalLoadCase(4,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).load.static(5) = model.nodes(n).nodalLoadCase(5,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).load.static(6) = model.nodes(n).nodalLoadCase(6,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            if all(model.nodes(n).load.static == 0)
                                model.nodes(n).load.static = [];
                            end
                        end
                        if size(model.nodes(n).nodalLoadCase,1) <= 6
                            allLogic = 1;
                        else
                            allLogic = all(model.nodes(n).nodalLoadCase(7:12,:) == 0);
                        end
                        if all(allLogic == 1)
                            model.nodes(n).prescDispl = [];
                        else  % if there are prescribed displacements, ,multiply them
                            % by their proper comb factors, sum the results and set
                            % them to prescDispl.
                            model.nodes(n).prescDispl(1) = model.nodes(n).nodalLoadCase(7,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).prescDispl(2) = model.nodes(n).nodalLoadCase(8,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).prescDispl(3) = model.nodes(n).nodalLoadCase(9,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).prescDispl(4) = model.nodes(n).nodalLoadCase(10,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).prescDispl(5) = model.nodes(n).nodalLoadCase(11,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            model.nodes(n).prescDispl(6) = model.nodes(n).nodalLoadCase(12,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-model.nlc);
                            if all(model.nodes(n).prescDispl == 0)
                                model.nodes(n).prescDispl = [];
                            end
                        end
                    else
                        model.nodes(n).load.static = [];
                        model.nodes(n).prescDispl = [];
                    end
                end
                
                % Consider load cases factors for the selected combination and allocate
                % the resulting element loads and prescribed displacements in each
                % Lelem object.
                % obs.: elemLoadCase(:,i) = [  unifDir
                %                                qx
                %                                qy
                %                                qz        (elemLoadCase is a matrix,
                %                             linearDir     each column refers to a
                %                                qx1        specific load case)
                %                                qy1
                %                                qz1
                %                                qx2
                %                                qy2
                %                                qz2
                %                                dtx
                %                                dty
                %                                dtz   ]
                %
                % obs.2: loadComb is a matrix that holds combination factors for load
                % case combinations, where each column refers to a specific load comb.
                
                for e = 1:model.nel
                    if isempty(model.elems(e).load.elemLoadCase) == 0
                        % Clear previous uniform loads
                        model.elems(e).load.uniformGbl = [];
                        model.elems(e).load.uniformLcl = [];
                        allLogic = all(model.elems(e).load.elemLoadCase(2:4,:) == 0);
                        if all(allLogic == 1)
                            model.elems(e).load.uniformDir = 0;
                        else  % if there are uniform loads, set their local and global
                            % components and consider their comb factors
                            unifLoad = zeros(size(model.elems(e).load.elemLoadCase,2),3);
                            for i = 1:size(model.elems(e).load.elemLoadCase,2)
                                unifLoad(i,:) = model.loadComb(i,lc-model.nlc) * model.elems(e).load.elemLoadCase(2:4,i);
                                if ~all(unifLoad(i,:) == 0)
                                    model.elems(e).load.setUnifLoad(unifLoad(i,:),model.elems(e).load.elemLoadCase(1,i));
                                end
                            end
                        end
                        
                        % Clear previous linear loads
                        model.elems(e).load.linearGbl = [];
                        model.elems(e).load.linearLcl = [];
                        if size(model.elems(e).load.elemLoadCase,1) <= 5
                            allLogic = 1;
                        else
                            allLogic = all(model.elems(e).load.elemLoadCase(6:11,:) == 0);
                        end
                        if all(allLogic == 1)
                            model.elems(e).load.linearDir = 0;
                        else  % if there are linear loads, set their local and global
                            % components and consider their comb factors
                            linLoad = zeros(size(model.elems(e).load.elemLoadCase,2),6);
                            for i = 1:size(model.elems(e).load.elemLoadCase,2)
                                linLoad(i,:) = model.loadComb(i,lc-model.nlc) * model.elems(e).load.elemLoadCase(6:11,i);
                                if ~all(linLoad(i,:) == 0)
                                    model.elems(e).load.setLinearLoad(linLoad(i,:),model.elems(e).load.elemLoadCase(5,i));
                                end
                            end
                            % Avoids numeric problems with new loads
                            for i = 1:size(model.elems(e).load.linearLcl,1)
                                if abs(model.elems(e).load.linearLcl(i)) <= 10^-10
                                    model.elems(e).load.linearLcl(i) = 0;
                                end
                                if abs(model.elems(e).load.linearGbl(i)) <= 10^-10
                                    model.elems(e).load.linearGbl(i) = 0;
                                end
                            end
                        end
                        % Set thermal loads
                        if size(model.elems(e).load.elemLoadCase,1) <= 11
                            model.elems(e).load.tempVar_X = 0;
                            model.elems(e).load.tempVar_Y = 0;
                            model.elems(e).load.tempVar_Z = 0;
                        else
                            model.elems(e).load.tempVar_X = model.elems(e).load.elemLoadCase(12,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-model.nlc);
                            model.elems(e).load.tempVar_Y = model.elems(e).load.elemLoadCase(13,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-model.nlc);
                            model.elems(e).load.tempVar_Z = model.elems(e).load.elemLoadCase(14,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-model.nlc);
                        end
                    else
                        model.elems(e).load.uniformDir = 0;
                        model.elems(e).load.uniformGbl = [];
                        model.elems(e).load.uniformLcl = [];
                        model.elems(e).load.linearDir = 0;
                        model.elems(e).load.linearGbl = [];
                        model.elems(e).load.linearLcl = [];
                        model.elems(e).load.tempVar_X = 0;
                        model.elems(e).load.tempVar_Y = 0;
                        model.elems(e).load.tempVar_Z = 0;
                    end
                end
            end
                    
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of the entire model.
        function clean(model)
            % Clean data structure of the Anm object
            if model.anm ~= 0
                model.anm.clean();
            end
            
            % Clean data structure of all Material objects
            for m = 1:model.nmat
                model.materials(m).clean();
            end
            
            % Clean data structure of all Section objects
            for s = 1:model.nsec
                model.sections(s).clean();
            end
            
            % Clean data structure of all Node objects
            for n = 1:model.nnp
                model.nodes(n).clean();
            end
            
            % Clean data structure of all Element objects and Lelem objects
            for e = 1:model.nel
                model.elems(e).load.clean();
                model.elems(e).clean();
            end
            
            % Clean data structure of Result object
            if ~isempty(model.results)
                model.results.clean();
            end
            
            % Clean data structure of the model object
            model.drv            = [];
            model.anm            = [];
            model.materials      = [];
            model.sections       = [];
            model.nodes          = [];
            model.elems          = [];
            model.srjoints       = [];
            model.results        = [];
            model.K              = [];
            model.M              = [];
            model.C              = [];
            model.F              = [];
            model.D              = [];
            model.W              = [];
            model.V              = [];
            model.damping        = 0;
            model.xi             = [];
            model.massDampCoeff  = 0;
            model.stiffDampCoeff = 0;
            model.c0             = [];
            model.n_steps        = 0;
            model.t              = 0;
            model.n_modes        = 0;
            model.whichSolver    = 0;
            model.mass_type      = 1;
            model.mass_mi        = 1;
            model.ID             = [];
            model.nmat           = 0;
            model.nsec           = 0;
            model.nnp            = 0;
            model.njoints        = 0;
            model.nel            = 0;
            model.neq            = 0;
            model.neqfree        = 0;
            model.neqfixed       = 0;
            model.neqspring      = 0;
            model.gblSprStiff    = [];
            model.nlc            = 0;
            model.ncomb          = 0;
            model.strLc          = {};
            model.strComb        = {};
            model.loadComb       = [];
            model.loadCombID     = [];
            model.timeFcns       = {};
            model.strTimeFcns    = {};
        end
    end
end