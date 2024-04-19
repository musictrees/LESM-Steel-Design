%% Read file Function (Version 2.1)
%
% This is an auxiliary file that contains functions to read a
% neutral-format file with the _.lsm_ extension 
% This neutral-format file contains all information about a linear elements
% structural model.
%
%% Main function
% Output:
%  print:    object of the Print class
%  draw:     object of the Draw class
%  nclc:     current load case identifier
% Input arguments:
%  fid:      integer identifier of the input file
%  model:    handle to an object of the model class
%  GUI_Mode: flag for the graphical version of LESM being used
function [print,draw,nclc] = readFile_v2_1(fid,model,GUI_Mode)
    % Initialize nclc
    nclc = 1;

    % Initialize flag for while loop
    keep_going = true;
    
    while keep_going && ~feof(fid)
        % Get file line
        tline = fgetl(fid);
        % Get rid of blank spaces
        string = deblank(tline);
        
        % Look for for tag strings
        switch string
            case '%HEADER.ANALYSIS'
                [print,draw] = readHeaderAnalysis(fid,model,GUI_Mode);
            case '%NODE.COORD'
                readNodeCoord(fid,model);
            case '%NODE.SUPPORT'
                readNodeSupport(fid,model);
            case '%NODE.INCLINED.SUPPORT'
                readNodeInclinedSupps(fid,model);
            case '%NODE.SUPPORT.SPRING'
                readSpringStiff(fid,model);
            case '%MATERIAL.ISOTROPIC'
                readMaterialIsotropic(fid,model);
            case '%SECTION.PROPERTY'
                readSectionProperty(fid,model);
            case '%BEAM.END.LIBERATION'
                rotlib = readBeamEndLiberation(fid);
            case '%ELEMENT.BEAM'
                readElementBeam(fid,model,rotlib);
            case '%ELEMENT.BEAM.INTERSECTIONS'
                if GUI_Mode
                    readElementIntersections(fid);
                end
            case '%ELEMENT.BEAM.JOINTS'
                readElementBeamJoints(fid,model);
            case '%LOAD'
                readLoadCase(fid,model);
            case '%LOAD.COMBINATION'
                readLoadComb(fid,model);
            case '%LOAD.CASE'
                nclc = readCurrentLoadCase(fid);
            case '%LOAD.CASE.NODAL.CASES'
                readNodalLoadCase(fid,model);
            case '%LOAD.CASE.BEAM.CASES'
                readElemLoadCase(fid,model);
            case '%END'
                setCurrentLoadCase(model,nclc);
                keep_going = false;    
        end
    end
end

%% Auxiliary functions
%--------------------------------------------------------------------------
function [print,draw] = readHeaderAnalysis(fid,model,GUI_Mode)
    % Get file line
    tline = fgetl(fid);
    % Get rid of blank spaces
    string = deblank(tline);
    
    % Get target analysis type
    switch string
        case '''TRUSS2D'''
            model.anm = Anm_Truss2D();
            print = Print_Truss2D(model);
            if GUI_Mode == true
                draw = Draw_Truss2D(model);
            else
                draw = [];
            end
        case '''FRAME2D'''
            model.anm = Anm_Frame2D();
            print = Print_Frame2D(model);
            if GUI_Mode == true
                draw = Draw_Frame2D(model);
            else
                draw = [];
            end
        case '''GRILLAGE'''
            model.anm = Anm_Grillage();
            print = Print_Grillage(model);
            if GUI_Mode == true
                draw = Draw_Grillage(model);
            else
                draw = [];
            end
        case '''TRUSS3D'''
            model.anm = Anm_Truss3D();
            print = Print_Truss3D(model);
            if GUI_Mode == true
                draw = Draw_Truss3D(model);
            else
                draw = [];
            end
        case '''FRAME3D'''
            model.anm = Anm_Frame3D();
            print = Print_Frame3D(model);
            if GUI_Mode == true
                draw = Draw_Frame3D(model);
            else
                draw = [];
            end
    end
end

%--------------------------------------------------------------------------
function readNodeCoord(fid,model)
    % Read and store total number of nodes and equations
    nnp = fscanf(fid,'%f',1);
    model.nnp = nnp;
    model.neq = nnp * model.anm.ndof;
    
    % Initialize vector of objects of the Node class
    nodes = initNodes(nnp);
    model.nodes = nodes;
    
    for n = 1:nnp
        b = fscanf(fid,'%f',4);
        % b(1) = node id
        % b(2) = X coordinate
        % b(3) = Y coordinate
        % b(4) = Z coordinate
        model.nodes(n).id = b(1);
        model.nodes(n).coord = [b(2) b(3) b(4)];
    end
end

%--------------------------------------------------------------------------
function nodes = initNodes(nnp)
    % Initialize vector of objects of the Node class
    nodes(1,nnp)  =  Node();
    
    % Initialize vector of objects of the Lnode class
    lnodes(1,nnp) = Lnode();
    
    % Set Lnode objects as properties of each respective Node object
    for n = 1:nnp
        nodes(n).load      = lnodes(n);
        nodes(n).load.node =  nodes(n);
    end
end

%--------------------------------------------------------------------------
function readNodeSupport(fid,model)
    % Read total number of nodes with essential boundary conditions (supports)
    nnebc = fscanf(fid,'%f',1);
    
    % Initialize number of fixed dof's
    model.neqfixed = 0;
    
    for n = 1:nnebc
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx
        % b(3) = Dy
        % b(4) = Dz
        % b(5) = Rx
        % b(6) = Ry
        % b(7) = Rz
        model.nodes(b(1)).ebc = [b(2) b(3) b(4) b(5) b(6) b(7)];
        for i = 1:6
           if b(i+1) == 1
               model.neqfixed = model.neqfixed + 1;
           end    
        end
    end
end

%--------------------------------------------------------------------------
function readNodeInclinedSupps(fid,model)
    % Read total number of nodes with inclined supports
    nnis = fscanf(fid,'%f',1);
    
    for n = 1:nnis
        b = fscanf(fid,'%f',[1 8]);
        % b(1) = node id
        % b(2) = dir_x
        % b(3) = dir_y
        % b(4) = dir_z
        % b(5) = flag_vy
        % b(6) = vy(1)
        % b(7) = vy(2)
        % b(8) = vy(3)
        dir = [b(2), b(3), b(4)];
        if ~b(5)
            model.nodes(b(1)).setInclinedSupp(dir);
        else
            vy = [b(6), b(7), b(8)];
            model.nodes(b(1)).setInclinedSupp(dir,vy);
        end
    end
end

%--------------------------------------------------------------------------
function readSpringStiff(fid,model)
    % Read total number of nodes with spring supports
    nnsp = fscanf(fid,'%f',1);
    
    % Initialize number of springs
    model.neqspring = 0;
    
    for n = 1:nnsp
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = kdx
        % b(3) = kdy
        % b(4) = kdz
        % b(5) = krx
        % b(6) = kry
        % b(7) = krz
        model.nodes(b(1)).springStiff = [b(2) b(3) b(4) b(5) b(6) b(7)];
        for i = 1:6
           if b(i+1) ~= 0
               model.neqspring = model.neqspring + 1;
           end
        end    
    end
end

%--------------------------------------------------------------------------
function readMaterialIsotropic(fid,model)
    % Read and store total number of materials
    nmat = fscanf(fid,'%f',1);
    model.nmat = nmat;
    
    % Initialize vector of objects of the Material class
    materials(1,nmat) = Material();
    model.materials = materials;
    
    for m = 1:nmat
        b = fscanf(fid,'%f',4);
        % b(1) = material id
        % b(2) = elasticity modulus [MPa]
        % b(3) = possion ratio
        % b(4) = thermal expansion coefficient  [/oC]
        model.materials(m) = Material(b(1),1e3*b(2),b(3),b(4));
    end
end

%--------------------------------------------------------------------------
function readSectionProperty(fid,model)
    % Read and store total number of cross-sections
    nsec = fscanf(fid,'%f',1);
    model.nsec = nsec;
    
    % Initialize vector of objects of the Section class
    sections(1,nsec) = Section();
    model.sections = sections;
    
    for s = 1:nsec
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = Cross-section id
        % b(2) = Ax [cm2]
        % b(3) = Ay [cm2]
        % b(4) = Az [cm2]
        % b(5) = Ix [cm4]
        % b(6) = Iy [cm4]
        % b(7) = Iz [cm4]
        % b(8) = Hy [cm]
        % b(9) = Hz [cm]
        model.sections(s).id = b(1);
        model.sections(s).area_x = 1e-4*b(2);
        model.sections(s).area_y = 1e-4*b(3);
        model.sections(s).area_z = 1e-4*b(4);
        model.sections(s).inertia_x = 1e-8*b(5);
        model.sections(s).inertia_y = 1e-8*b(6);
        model.sections(s).inertia_z = 1e-8*b(7);
        model.sections(s).height_y = 1e-2*b(8);
        model.sections(s).height_z = 1e-2*b(9);
    end
end

%--------------------------------------------------------------------------
function rotlib = readBeamEndLiberation(fid)
    include_constants;
    
    % Read number of beam end rotation liberation properties
    nrotlib = fscanf(fid,'%f',1);
    rotlib = zeros(2,nrotlib);
    
    % Read beam end rotation liberation (hinge) data
    for i = 1:nrotlib
        b = fscanf(fid,'%f',[1 13]);
        % b( 1) = rotation liberation property id (not used)
        % b( 2) = idx
        % b( 3) = idy
        % b( 4) = idz
        % b( 5) = irx
        % b( 6) = iry
        % b( 7) = irz
        % b( 8) = jdx
        % b( 9) = jdy
        % b(10) = jdz
        % b(11) = jrx
        % b(12) = jry
        % b(13) = jrz
        % If rotation liberation in any of the directions is free, it is
        % considered that it is free in all directions, thus it will receive a
        % hinged flag
        if (b(5)==0) || (b(6)==0) || (b(7)==0)
            rotlib(1,i) = HINGED_END;
        elseif (b(5)==2) || (b(6)==2) || (b(7)==2)
            rotlib(1,i) = SEMIRIGID_END;
        else
            rotlib(1,i) = CONTINUOUS_END;
        end
        
        if (b(11)==0) || (b(12)==0) || (b(13)==0)
            rotlib(2,i) = HINGED_END;
        elseif (b(11)==2) || (b(12)==2) || (b(13)==2)
            rotlib(2,i) = SEMIRIGID_END;
        else
            rotlib(2,i) = CONTINUOUS_END;
        end
    end
end

%--------------------------------------------------------------------------
function readElementBeam(fid,model,rotlib)
    % Read and store total number of elements
    nel = fscanf(fid,'%f',1);
    model.nel = nel;
    
    % Initialize vector of objects of the Elem class
    elems(1,nel) = Elem();
    model.elems = elems;
    
    for e = 1:nel
        b = fscanf(fid,'%f',[1 10]);
        % b(1)  = element id
        % b(2)  = element type (0 = Navier, 1 = Timoshenko)
        % b(3)  = material id
        % b(4)  = section property id
        % b(5)  = end release id
        % b(6)  = init node id
        % b(7)  = final node id
        % b(8)  = vz_X
        % b(9)  = vz_Y
        % b(10) = vz_Z
        model.elems(e) = Elem(b(2),...
                            model.anm,...
                            model.materials(b(3)),...
                            model.sections(b(4)),...
                            [model.nodes(b(6)), model.nodes(b(7))],...
                            rotlib(1,b(5)), rotlib(2,b(5)),...
                            [b(8), b(9), b(10)]);
    end
end

%--------------------------------------------------------------------------
function readElementIntersections(fid)
    % Read number of intersections
    nints = fscanf(fid,'%f',1);

    for ni = 1:nints
        b = fscanf(fid,'%f',[1 5]);
        % b(1)  = intersection id
        % b(2)  = x
        % b(3)  = y
        % b(4)  = z
        % b(5)  = number of elements
        intersections(b(1)).coord = [b(2), b(3), b(4)];        %#ok<AGROW>
        intersections(b(1)).elems = fscanf(fid,'%f',[1 b(5)]); %#ok<AGROW>
    end
    
    % Store intersections structs in root (GUI_Mode only)
    setappdata(0,'intersections',intersections)
end

%--------------------------------------------------------------------------
function readElementBeamJoints(fid,model)
    include_constants;

    % Read and store total number of semi-rigid joints
    nsrj = fscanf(fid,'%f',1);
    model.njoints = nsrj;
    
    % Update number of equations
    model.neq = model.neq + model.njoints * model.anm.nrdof;
    
    % Initilize counter
    nj = 1;
    aux = 1;
    
    while nj <= nsrj
        b = fscanf(fid,'%f',[1 9]);
        % b(1)  = element id
        % b(2)  = init joint  (0 = cont./hinged, 1 = semi-rigid)
        % b(3)  = final joint (0 = cont./hinged, 1 = semi-rigid)
        % b(4)  = krix
        % b(5)  = kriy
        % b(6)  = kriz
        % b(7)  = krfx
        % b(8)  = krfy
        % b(9)  = krfz
        
        if b(2)
            model.elems(b(1)).hingei = SEMIRIGID_END;
            model.elems(b(1)).kri = [b(4), b(5), b(6)];
            nj = nj + 1;
        end
        
        if b(3)
            model.elems(b(1)).hingef = SEMIRIGID_END;
            model.elems(b(1)).krf = [b(7), b(8), b(9)];
            nj = nj + 1;
        end
        
        % Auxiliar counter - avoids infinte loop in while statement, in
        % case there is an error in the input file
        aux = aux + 1;
        if aux > nsrj
            nj = aux;
        end
    end
end

%--------------------------------------------------------------------------
function readLoadCase(fid,model)
    % Read total number of load cases
    model.nlc = fscanf(fid,'%f',1);
    str0 = fgetl(fid); %#ok<NASGU>
    for i = 1:model.nlc
       string(i) = {fgetl(fid)}; %#ok<AGROW>
    end
    model.strLc = string';
end

%--------------------------------------------------------------------------
function nclc = readCurrentLoadCase(fid)
    % Read current load case id
    nclc = fscanf(fid,'%f',1);
     
end

%--------------------------------------------------------------------------
function readLoadComb(fid,model)
     % Read total number of load case combinations
     model.ncomb = fscanf(fid,'%f',1);
     for i = 1:model.ncomb
        str0 = fgetl(fid); %#ok<NASGU>
        string(i) = {fgetl(fid)};  %#ok<AGROW>
        b = fscanf(fid,'%f',[1 model.nlc]); 
        % b(1) = comb factor 1
        % b(2) = comb factor 2
        % (...)
        % b(nlc) = comb factor nlc
        if isempty(model.loadComb) == 1
            model.loadComb = zeros(model.nlc,model.ncomb);
        end
        model.loadComb(:,i) = b';
        if isempty(model.loadCombID) == 1
            model.loadCombID = zeros(model.nlc,model.ncomb);
        end
        model.loadCombID(:,i) = (~b' - 1) * -1;
     end
     model.strComb = string';
end

%--------------------------------------------------------------------------
function readNodalLoadCase(fid,model)
    % Read total number of nodes with applied load and/or prescribed
    % displacements in at least one load case
    nnlc = fscanf(fid,'%f',1);

    for i = 1:nnlc
        b = fscanf(fid,'%f',[1 14]);
        % b(1) = node id
        % b(2) = Fx [kN]
        % b(3) = Fy [kN]
        % b(4) = Fz [kN]
        % b(5) = Mx [kNm]
        % b(6) = My [kNm]
        % b(7) = Mz [kNm]
        % b(8) = Dx [mm]
        % b(9) = Dy [mm]
        % b(10) = Dz [mm]
        % b(11) = Rx [rad]
        % b(12) = Ry [rad]
        % b(13) = Rz [rad]
        % b(14) = number of load cases with loads in this specific node
        model.nodes(b(1)).nodalLoadCase(:,1) = [b(2); b(3); b(4); b(5);...
                                              b(6); b(7); 1e-3*b(8);...
                                              1e-3*b(9); 1e-3*b(10);...
                                              b(11); b(12); b(13)];                                 
        if b(14) ~= 1
            for j = 2:b(14)
                c = fscanf(fid,'%f',[1 12]);
                % c(1) = Fx [kN]
                % c(2) = Fy [kN]
                % c(3) = Fz [kN]
                % c(4) = Mx [kNm]
                % c(5) = My [kNm]
                % c(6) = Mz [kNm]
                % c(7) = Dx [mm]
                % c(8) = Dy [mm]
                % c(9) = Dz [mm]
                % c(10) = Rx [rad]
                % c(11) = Ry [rad]
                % c(12) = Rz [rad]
                model.nodes(b(1)).nodalLoadCase(:,j) = [c(1); c(2); c(3);...
                                                      c(4); c(5); c(6);...
                                                      1e-3*c(7); 1e-3*c(8);... 
                                                      1e-3*c(9); c(10);...
                                                      c(11); c(12)];
            end
        end
    end
end

%--------------------------------------------------------------------------
function readElemLoadCase(fid,model)
    % Read total number of nodes with applied load and/or prescribed
    % displacements in at least one load case
    nelc = fscanf(fid,'%f',1);

    for i = 1:nelc
        b = fscanf(fid,'%f',[1 16]);
        % b(1) = elem id
        % b(2) = unifDir
        % b(3) = qX [kN/m]
        % b(4) = qY [kN/m]
        % b(5) = qZ [kN/m]
        % b(6) = linearDir
        % b(7) = qX1 [kN/m]
        % b(8) = qY1 [kN/m]
        % b(9) = qZ1 [kN/m]
        % b(10) = qX2 [kN/m]
        % b(11) = qY2 [kN/m]
        % b(12) = qZ2 [kN/m]
        % b(13) = dt_x [oC]
        % b(14) = dt_y [oC]
        % b(15) = dt_z [oC]
        % b(16) = number of load cases with loads in this specific node
        model.elems(b(1)).load.elemLoadCase(:,1) = [b(2); b(3); b(4); b(5);...
                                                 b(6); b(7); b(8); b(9);...
                                                 b(10); b(11); b(12);...
                                                 b(13); b(14); b(15)];
        if b(16) ~= 1
            for j = 2:b(16)
                c = fscanf(fid,'%f',[1 14]);
                % c(1) = unifDir
                % c(2) = qX [kN/m]
                % c(3) = qY [kN/m]
                % c(4) = qZ [kN/m]
                % c(5) = linearDir
                % c(6) = qX1 [kN/m]
                % c(7) = qY1 [kN/m]
                % c(8) = qZ1 [kN/m]
                % c(9) = qX2 [kN/m]
                % c(10) = qY2 [kN/m]
                % c(11) = qZ2 [kN/m]
                % c(12) = dt_x [oC]
                % c(13) = dt_y [oC]
                % c(14) = dt_z [oC]
                model.elems(b(1)).load.elemLoadCase(:,j) = [c(1); c(2); c(3);...
                                                         c(4); c(5); c(6);...
                                                         c(7); c(8); c(9);...
                                                         c(10); c(11); c(12);...
                                                         c(13); c(14)];
            end
        end
    end
end

%--------------------------------------------------------------------------
function setCurrentLoadCase(model,lc)
% Get number of load cases
nlc = model.nlc;

% Check if currebt load case is a case or a combination
if lc <= nlc
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
            if all(model.nodes(n).nodalLoadCase(1:6,lc) == 0) == 1
                model.nodes(n).load.static = []; 
            else %if there are any nodal loads, set them to load.static
                model.nodes(n).load.static = model.nodes(n).nodalLoadCase(1:6,lc);
            end
            if size(model.nodes(n).nodalLoadCase,1) <= 6
                model.nodes(n).prescDispl = [];
            elseif all(model.nodes(n).nodalLoadCase(7:12,lc) == 0) == 1
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
            % Initialize uniform loads
            model.elems(e).load.uniformGbl = [];    
            model.elems(e).load.uniformLcl = [];
            if all(model.elems(e).load.elemLoadCase(2:4,lc) == 0) == 1
                model.elems(e).load.uniformDir = 0;  
            else  % if there are uniform loads, set their local and global components
                model.elems(e).load.uniformDir = model.elems(e).load.elemLoadCase(1,lc);
                model.elems(e).load.setUnifLoad((model.elems(e).load.elemLoadCase(2:4,lc))',model.elems(e).load.elemLoadCase(1,lc));
            end
            % Initialize linear loads
            model.elems(e).load.linearGbl = [];
            model.elems(e).load.linearLcl = [];
            if size(model.elems(e).load.elemLoadCase,1) <= 5
                model.elems(e).load.linearDir = 0;
            elseif all(model.elems(e).load.elemLoadCase(6:11,lc) == 0) == 1
                model.elems(e).load.linearDir = 0;
            else  % if there are linear loads, set their local and global components
                model.elems(e).load.linearDir = model.elems(e).load.elemLoadCase(5,lc);
                model.elems(e).load.setLinearLoad((model.elems(e).load.elemLoadCase(6:11,lc))',model.elems(e).load.elemLoadCase(5,lc));
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
    
else  % if lc > nlc, the current load case is a combination
    
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
            if all(allLogic == 1) == 1
                model.nodes(n).load.static = []; 
            else  % if there are nodal loads, ,multiply them by their 
                  % proper comb factors, sum the results and set them to 
                  % load.static.
                model.nodes(n).load.static(1) = model.nodes(n).nodalLoadCase(1,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).load.static(2) = model.nodes(n).nodalLoadCase(2,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).load.static(3) = model.nodes(n).nodalLoadCase(3,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).load.static(4) = model.nodes(n).nodalLoadCase(4,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).load.static(5) = model.nodes(n).nodalLoadCase(5,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).load.static(6) = model.nodes(n).nodalLoadCase(6,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                if all(model.nodes(n).load.static == 0) == 1
                    model.nodes(n).load.static = [];
                end    
            end
            if size(model.nodes(n).nodalLoadCase,1) <= 6
                allLogic = 1;
            else
                allLogic = all(model.nodes(n).nodalLoadCase(7:12,:) == 0);
            end
            if all(allLogic == 1) == 1
                model.nodes(n).prescDispl = [];
            else  % if there are prescribed displacements, ,multiply them 
                  % by their proper comb factors, sum the results and set  
                  % them to prescDispl.
                model.nodes(n).prescDispl(1) = model.nodes(n).nodalLoadCase(7,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).prescDispl(2) = model.nodes(n).nodalLoadCase(8,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).prescDispl(3) = model.nodes(n).nodalLoadCase(9,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).prescDispl(4) = model.nodes(n).nodalLoadCase(10,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).prescDispl(5) = model.nodes(n).nodalLoadCase(11,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                model.nodes(n).prescDispl(6) = model.nodes(n).nodalLoadCase(12,:) * model.loadComb(1:size(model.nodes(n).nodalLoadCase,2),lc-nlc);
                if all(model.nodes(n).prescDispl == 0) == 1
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
            % Initilize uniform loads
            model.elems(e).load.uniformGbl = [];    
            model.elems(e).load.uniformLcl = [];
            allLogic = all(model.elems(e).load.elemLoadCase(2:4,:) == 0);
            if all(allLogic == 1) == 1
                model.elems(e).load.uniformDir = 0;
            else  % if there are uniform loads, set their local and global
                  % components and consider their comb factors
                unifLoad = zeros(size(model.elems(e).load.elemLoadCase,2),3);
                for i = 1:size(model.elems(e).load.elemLoadCase,2)
                    unifLoad(i,:) = model.loadComb(i,lc-nlc) * model.elems(e).load.elemLoadCase(2:4,i);
                    if all(unifLoad(i,:) == 0) == 0
                        model.elems(e).load.setUnifLoad(unifLoad(i,:),model.elems(e).load.elemLoadCase(1,i));
                    end    
                end
            end  
            
            % Initilize linear loads
            model.elems(e).load.linearGbl = [];
            model.elems(e).load.linearLcl = [];
            if size(model.elems(e).load.elemLoadCase,1) <= 5
                allLogic = 1;
            else
                allLogic = all(model.elems(e).load.elemLoadCase(6:11,:) == 0);
            end
            if all(allLogic == 1) == 1
                model.elems(e).load.linearDir = 0;
            else  % if there are linear loads, set their local and global
                  % components and consider their comb factors    
                linLoad = zeros(size(model.elems(e).load.elemLoadCase,2),6);
                for i = 1:size(model.elems(e).load.elemLoadCase,2)
                    linLoad(i,:) = model.loadComb(i,lc-nlc) * model.elems(e).load.elemLoadCase(6:11,i);
                    if all(linLoad(i,:) == 0) == 0
                        model.elems(e).load.setLinearLoad(linLoad(i,:),model.elems(e).load.elemLoadCase(5,i));
                    end    
                end
            end
            % Set thermal loads
            if size(model.elems(e).load.elemLoadCase,1) <= 11
                model.elems(e).load.tempVar_X = 0;
                model.elems(e).load.tempVar_Y = 0;  
                model.elems(e).load.tempVar_Z = 0; 
            else
                model.elems(e).load.tempVar_X = model.elems(e).load.elemLoadCase(12,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-nlc);
                model.elems(e).load.tempVar_Y = model.elems(e).load.elemLoadCase(13,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-nlc);  
                model.elems(e).load.tempVar_Z = model.elems(e).load.elemLoadCase(14,:) * model.loadComb(1:size(model.elems(e).load.elemLoadCase,2),lc-nlc);
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