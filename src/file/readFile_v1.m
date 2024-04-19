%% Read file Function (Versions 1.0 and 1.1)
%
% This is an auxiliary file that contains functions to read a
% neutral-format file with the _.lsm_ extension 
% This neutral-format file contains all information about a linear element
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
function [print,draw,nclc] = readFile_v1(fid,model,GUI_Mode)
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
            case '%NODE.SUPPORT.SPRING'
                readSpringStiff(fid,model);
            case '%MATERIAL.ISOTROPIC'
                readMaterialIsotropic(fid,model);
            case '%SECTION.PROPERTY'
                readSectionProperty(fid,model);
            case '%BEAM.END.LIBERATION'
                rotlib = readBeamEndLiberation(fid);
            case '%ELEMENT.BEAM.NAVIER'
                readElementBeamNavier(fid,model,rotlib);
            case '%ELEMENT.BEAM.TIMOSHENKO'
                readElementBeamTimoshenko(fid,model,rotlib);
            case '%LOAD.CASE.NODAL.FORCE'
                readNodeAppliedForceMoment(fid,model);
            case '%LOAD.CASE.NODAL.DISPLACEMENT'
                readNodePrescDisplRot(fid,model);
            case '%LOAD.CASE.BEAM.UNIFORM'
                readElemAppliedUnifLoad(fid,model);
            case '%LOAD.CASE.BEAM.LINEAR'
                readElemAppliedLinearLoad(fid,model);
            case '%LOAD.CASE.BEAM.TEMPERATURE'
                readElemThermalLoad(fid,model);
            case '%END'
                nclc = loadCases_v1(model);
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
        else
            rotlib(1,i) = CONTINUOUS_END;
        end
        
        if (b(11)==0) || (b(12)==0) || (b(13)==0)
            rotlib(2,i) = HINGED_END;
        else
            rotlib(2,i) = CONTINUOUS_END;
        end
    end
end

%--------------------------------------------------------------------------
function readElementBeamNavier(fid,model,rotlib)
    % Read and store total number of elements
    nel = fscanf(fid,'%f',1);
    model.nel = nel;
    
    % Initialize vector of objects of the Elem_Navier sub-class
    elems(1,nel) = Elem();
    model.elems = elems;
    
    for e = 1:nel
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = element id
        % b(2) = material id
        % b(3) = section property id
        % b(4) = end release id
        % b(5) = init node id
        % b(6) = final node id
        % b(7) = vz_X
        % b(8) = vz_Y
        % b(9) = vz_Z
        model.elems(e) = Elem(0,model.anm,...
                           model.materials(b(2)),...
                           model.sections(b(3)),...
                           [model.nodes(b(5)), model.nodes(b(6))],...
                           rotlib(1,b(4)), rotlib(2,b(4)),...
                           [b(7), b(8), b(9)]);
    end
end

%--------------------------------------------------------------------------
function readElementBeamTimoshenko(fid,model,rotlib)
    % Read and store total number of elements
    nel = fscanf(fid,'%f',1);
    model.nel = nel;
    
    % Initialize vector of objects of the Elem_Timoshenko class
    elems(1,nel) = Elem();
    model.elems = elems;
    
    for e = 1:nel
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = element id
        % b(2) = material id
        % b(3) = section property id
        % b(4) = end release id
        % b(5) = init node id
        % b(6) = final node id
        % b(7) = vz_X
        % b(8) = vz_Y
        % b(9) = vz_Z
        
        % Create element object
        model.elems(e) = Elem(1,model.anm,...
                           model.materials(b(2)),...
                           model.sections(b(3)),...
                           [model.nodes(b(5)), model.nodes(b(6))],...
                           rotlib(1,b(4)), rotlib(2,b(4)),...
                           [b(7), b(8), b(9)]);
    end
end

%--------------------------------------------------------------------------
function readNodeAppliedForceMoment(fid,model)
    % Read total number of nodal loads
    n_nodalload = fscanf(fid,'%f',1);
    
    for i = 1:n_nodalload
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Fx [kN]
        % b(3) = Fy [kN]
        % b(4) = Fz [kN]
        % b(5) = Mx [kNm]
        % b(6) = My [kNm]
        % b(7) = Mz [kNm]
        model.nodes(b(1)).load.static = [b(2) b(3) b(4) b(5) b(6) b(7)];
    end
end

%--------------------------------------------------------------------------
function readNodePrescDisplRot(fid,model)
    % Read number of nodes with prescribed displacements
    nnprescdispl = fscanf(fid,'%f',1);
    
    for i = 1:nnprescdispl
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx [mm]
        % b(3) = Dy [mm]
        % b(4) = Dz [mm]
        % b(5) = Rx [rad]
        % b(6) = Ry [rad]
        % b(7) = Rz [rad]
        model.nodes(b(1)).prescDispl(1) = 1e-3*b(2);
        model.nodes(b(1)).prescDispl(2) = 1e-3*b(3);
        model.nodes(b(1)).prescDispl(3) = 1e-3*b(4);
        model.nodes(b(1)).prescDispl(4) = b(5);
        model.nodes(b(1)).prescDispl(5) = b(6);
        model.nodes(b(1)).prescDispl(6) = b(7);
    end
end

%--------------------------------------------------------------------------
function readElemAppliedUnifLoad(fid,model)
    % Read number of elements with applied uniformely distrib. load
    n_uniformload = fscanf(fid,'%f',1);
    
    for i = 1:n_uniformload
        b = fscanf(fid,'%f',[1 2]);
        % b(1) = element id
        % b(2) = direction (global or local system)
        model.elems(b(1)).load.uniformDir = b(2);
        
        c = fscanf(fid,'%f',[1 3]);
        % c(1) = Qx [kN/m]
        % c(2) = Qy [kN/m]
        % c(3) = Qz [kN/m]
        model.elems(b(1)).load.setUnifLoad([c(1),c(2),c(3)],b(2));
    end
end

%--------------------------------------------------------------------------
function readElemAppliedLinearLoad(fid,model)
    % Read number of elements with applied linearly distrib. load
    n_linearload = fscanf(fid,'%f',1);
    
    for i = 1:n_linearload
        b = fscanf(fid,'%f',[1 2]);
        % b(1) = element id
        % b(2) = direction (local or global system)
        model.elems(b(1)).load.linearDir = b(2);
        
        c = fscanf(fid,'%f',[1 3]);
        % c(1) = Qxi [kN/m]
        % c(2) = Qyi [kN/m]
        % c(3) = Qzi [kN/m]
        d = fscanf(fid,'%f',[1 3]);
        % d(1) = Qxj [kN/m]
        % d(2) = Qyj [kN/m]
        % d(3) = Qzj [kN/m]
        model.elems(b(1)).load.setLinearLoad([c(1),c(2),c(3),d(1),d(2),d(3)],b(2));
    end
end

%--------------------------------------------------------------------------
function readElemThermalLoad(fid,model)
    % Read number of elements with thermal load
    n_tempvar = fscanf(fid,'%f',1);
    
    for i = 1:n_tempvar
        b = fscanf(fid,'%f',[1 4]);
        % b(1) = element id
        % b(2) = dtx [oC]
        % b(3) = dty [oC]
        % b(4) = dtz [oC]
        model.elems(b(1)).load.tempVar_X = b(2);
        model.elems(b(1)).load.tempVar_Y = b(3);
        model.elems(b(1)).load.tempVar_Z = b(4);
    end
end    

%--------------------------------------------------------------------------
function nclc = loadCases_v1(model)
    % Set current load case as the first (and only) one
    nclc = 1;
    
    % Initialize load case parameters in the model object
    model.nlc = 1;
    model.strLc = {'CASE 01'};
    
    % Initialize load combinations parameters in the model object
    model.ncomb = 0;
    model.strComb = {' '};
    
    % Set nodal loads and prescribed displacements as one nodal load case 
    for n = 1:model.nnp
        if ~isempty(model.nodes(n).load.static)
            model.nodes(n).nodalLoadCase = zeros(12,1);
            model.nodes(n).nodalLoadCase(1:6) = model.nodes(n).load.static;
        end
        if ~isempty(model.nodes(n).prescDispl)
            if isempty(model.nodes(n).nodalLoadCase) == 1
                model.nodes(n).nodalLoadCase = zeros(12,1);
            end
            model.nodes(n).nodalLoadCase(7:12) = model.nodes(n).prescDispl;
        end
    end
    
    % Set distributed and thermal loads as one element load case
    for e = 1:model.nel
        if ~isempty(model.elems(e).load.uniformGbl)
            model.elems(e).load.elemLoadCase = zeros(14,1);
            unifDir = model.elems(e).load.uniformDir;
            model.elems(e).load.elemLoadCase(1) = unifDir;
            if unifDir == 0
                model.elems(e).load.elemLoadCase(2:4) = model.elems(e).load.uniformGbl;
            else
                model.elems(e).load.elemLoadCase(2:4) = model.elems(e).load.uniformLcl;
            end
        end
        if ~isempty(model.elems(e).load.linearGbl)
            if isempty(model.elems(e).load.elemLoadCase)
                model.elems(e).load.elemLoadCase = zeros(14,1);
            end    
            linDir = model.elems(e).load.linearDir;
            model.elems(e).load.elemLoadCase(5) = linDir;
            if linDir == 0
                model.elems(e).load.elemLoadCase(6:11) = model.elems(e).load.linearGbl;
            else
                model.elems(e).load.elemLoadCase(6:11) = model.elems(e).load.linearLcl;
            end
        end
        if isempty(model.elems(e).load.elemLoadCase)
            model.elems(e).load.elemLoadCase = zeros(14,1);
        end
        model.elems(e).load.elemLoadCase(12) = model.elems(e).load.tempVar_X;
        model.elems(e).load.elemLoadCase(13) = model.elems(e).load.tempVar_Y;
        model.elems(e).load.elemLoadCase(14) = model.elems(e).load.tempVar_Z;
    end
end 
