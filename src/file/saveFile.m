%% Save File Function
%
% This is an auxiliary function that writes all information about a
% structural model in a neutral-format file with the _.lsm_ extension.
%
%% Function Code
function saveFile(model,lsm)
include_constants;

% fprintf(lsm, 'NEUTRAL-FORMAT FILE\n');
% fprintf(lsm, 'This file stores all information about a structure model that can be read and\n');
% fprintf(lsm, 'loaded by the LESM (Linear Elements Structure Model) program.\n');
% fprintf(lsm, 'To modify the model, edit only the data below the %%TAGS\n\n');

%--------------------------------------------------------------------------
% Version header
%fprintf(lsm, '-----------------------------------------------------------------------\n');
%fprintf(lsm, 'Specify version number: ''X.XX''\n');
fprintf(lsm, '%%HEADER.VERSION\n');
fprintf(lsm, '3.0');
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Analysis header
%fprintf(lsm, '-----------------------------------------------------------------------\n');
%fprintf(lsm, 'Specify analysis model type: ''TRUSS2D'', ''FRAME2D'', ''Grillage'', ''TRUSS3D'' OR ''FRAME3D''\n');
fprintf(lsm, '%%HEADER.ANALYSIS\n');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value') - 1;
if anm == TRUSS2D_ANALYSIS
    fprintf(lsm, '''TRUSS2D''');
elseif anm == FRAME2D_ANALYSIS
    fprintf(lsm, '''FRAME2D''');
elseif anm == GRILLAGE_ANALYSIS
    fprintf(lsm, '''GRILLAGE''');
elseif anm == TRUSS3D_ANALYSIS
    fprintf(lsm, '''TRUSS3D''');
elseif anm == FRAME3D_ANALYSIS
    fprintf(lsm, '''FRAME3D''');
end
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Analysis algorithm header
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Specify analysis algorithm type: ''STATIC_LINEAR'', ''DYNAMIC_NEWMARK_LINEAR'', ''DYNAMIC_MODALSUP_LINEAR'',\n');
% fprintf(lsm, '                                 ''DYNAMIC_RK4_LINEAR'',''DYNAMIC_AM3_LINEAR'',''DYNAMIC_WILSON_LINEAR''\n');
fprintf(lsm, '%%HEADER.ANALYSIS.ALGORITHM\n');
switch model.whichSolver
    case STATIC_LINEAR
        fprintf(lsm, '''STATIC_LINEAR''');
    case DYNAMIC_NEWMARK_LINEAR
        fprintf(lsm, '''DYNAMIC_NEWMARK_LINEAR''');
    case DYNAMIC_MODALSUP_LINEAR
        fprintf(lsm, '''DYNAMIC_MODALSUP_LINEAR''');
    case DYNAMIC_RK4_LINEAR
        fprintf(lsm, '''DYNAMIC_RK4_LINEAR''');
    case DYNAMIC_AM3_LINEAR
        fprintf(lsm, '''DYNAMIC_AM3_LINEAR''');
    case DYNAMIC_WILSON_LINEAR
        fprintf(lsm, '''DYNAMIC_WILSON_LINEAR''');
end
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Analysis time steps header
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Provide number of time steps for dynamic analysis\n');
fprintf(lsm, '%%HEADER.ANALYSIS.TIME.STEPS\n');
fprintf(lsm, sprintf('%i',model.n_steps));
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Analysis time interval header
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Provide time interval [s] for dynamic analysis\n');
fprintf(lsm, '%%HEADER.ANALYSIS.TOLERANCE\n');
fprintf(lsm, sprintf('%d',model.t));
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Damping ratio
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Provide type of damping consideration and according inputs\n');
% fprintf(lsm, 'First line: Damping type\n');
% fprintf(lsm, 'Obs.: Type = 0 -> critical damping ratio of 1st mode.          Following line: xi\n');
% fprintf(lsm, '      Type = 1 -> critical damping ratio of 1st and 2nd modes. Following line: xi_1  xi_2\n');
% fprintf(lsm, '      Type = 2 -> rayleigh coefficients.                       Following line: mass_coeff  stiff_coeff\n');
fprintf(lsm, '%%MODEL.PROPERTY.DAMPING\n');
fprintf(lsm, sprintf('%i\n',model.damping));
if isempty(model.xi)
    fprintf(lsm, '0');
    fprintf(lsm, '\n\n');
else
    switch model.damping
        case XI_1ST_MODE
            fprintf(lsm, sprintf('%d',model.xi(1)));
            fprintf(lsm, '\n\n');
        case XI_1ST_2ND_MODE
            fprintf(lsm, sprintf('%d  %d',model.xi(1),model.xi(2)));
            fprintf(lsm, '\n\n');
        case RAYLEIGH_COEFFS
            fprintf(lsm, sprintf('%d  %d',model.massDampCoeff,model.stiffDampCoeff));
            fprintf(lsm, '\n\n');
    end
end

%--------------------------------------------------------------------------
% Mass properties
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Provide mass consideration properties\n');
% fprintf(lsm, 'Mass type, Mass mu coefficient\n');
% fprintf(lsm, 'Obs.: Type = 0 -> Lumped mass\n');
% fprintf(lsm, '      Type = 1 -> Consistent mass\n');
% fprintf(lsm, '      Type = 2 -> Mixed mass\n');
fprintf(lsm, '%%MODEL.PROPERTY.MASS\n');
fprintf(lsm, sprintf('%i %d',model.mass_type,model.mass_mi));
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Maximum number of modes
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Provide maximum number of vibration modes to be computed\n');
fprintf(lsm, '%%MODEL.PROPERTY.MODES\n');
fprintf(lsm, sprintf('%i',model.n_modes));
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Nodal coordinates
if model.nnp > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide nodal coordinates\n');
%     fprintf(lsm, 'First line: Total number of nodes\n');
%     fprintf(lsm, 'Following lines: Node ID, coord_X [m], coord_Y [m], coord_Z [m]\n');
    fprintf(lsm, '%%NODE.COORD\n');
    fprintf(lsm, '%d\n', model.nnp);
    for n = 1:model.nnp
        fprintf(lsm, '%d     %d  %d  %d\n', n, model.nodes(n).coord(1),...
                                               model.nodes(n).coord(2),...
                                               model.nodes(n).coord(3));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Support conditions
if model.nnp > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify essential boundary conditions (support conditions)\n');
%     fprintf(lsm, 'First line: Total number of nodes\n');
%     fprintf(lsm, 'Following lines: Node ID, ebc_dX, ebc_dY, ebc_dZ, ebc_rX, ebc_rY, ebc_rZ\n');
%     fprintf(lsm, 'ebc = 0 --> Free degree-of-freedom\n');
%     fprintf(lsm, 'ebc = 1 --> Fixed degree-of-freedom\n');
%     fprintf(lsm, 'ebc = 2 --> Spring degree-of-freedom\n');
    fprintf(lsm, '%%NODE.SUPPORT\n');
    fprintf(lsm, '%d\n', model.nnp);
    for n = 1:model.nnp
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
                model.nodes(n).ebc(1), model.nodes(n).ebc(2),...
                model.nodes(n).ebc(3), model.nodes(n).ebc(4),...
                model.nodes(n).ebc(5), model.nodes(n).ebc(6));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Inclined supports
if model.nnp > 0
    % counts inclined supps
    count = 0;
    % initialize aux vector
    inclinedSupps = zeros(1,model.nnp);
    for n = 1:model.nnp
        if model.nodes(n).isInclinedSupp
            count = count + 1;
            inclinedSupps(count) = n;
        end
    end
    if count ~= 0
        inclinedSupps = inclinedSupps(1,1:count);
    else
        inclinedSupps = [];
    end
    
    if count > 0
%         fprintf(lsm, '-----------------------------------------------------------------------\n');
%         fprintf(lsm, 'Specify inclined support conditions \n');
%         fprintf(lsm, 'First line: Total number of inclined supports \n');
%         fprintf(lsm, 'Following lines: Node ID, dir_x, dir_y, dir_z, Vy_Flag, dir_vy_x, dir_vy_y, dir_vy_z\n');
%         fprintf(lsm, 'obs.: For 2D models, dir_z = 0\n');
%         fprintf(lsm, 'obs.2: Vy_Flag = 0 -> There is no vy direction vector ; Vy_Flag = 1 -> There is a vy direction vector\n');
        fprintf(lsm, '%%NODE.INCLINED.SUPPORT\n');
        fprintf(lsm, '%d\n', count);
        for n = inclinedSupps
            if isempty(model.nodes(n).inclSupp_vy)
                flag_vy = false;
                vy = [0,0,0];
            else
                flag_vy = true;
                vy = model.nodes(n).inclSupp_vy;
            end
            fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d  %d\n', model.nodes(n).id,...
                                                                model.nodes(n).inclSuppDir(1),...
                                                                model.nodes(n).inclSuppDir(2),...
                                                                model.nodes(n).inclSuppDir(3),...
                                                                flag_vy,...
                                                                vy(1),...
                                                                vy(2),...
                                                                vy(3));
        end
        fprintf(lsm, '\n');
    end
end

%--------------------------------------------------------------------------
% Spring Stiffness
nnsp = 0;
for n = 1:model.nnp
    if isempty(model.nodes(n).springStiff) == 0
        nnsp = nnsp + 1;
    end
end    

if nnsp ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide spring stiffness coefficients\n');
%     fprintf(lsm, 'First line: Total number of nodes with springs\n');
%     fprintf(lsm, 'Following lines: Node ID, k_dX, k_dY, k_dZ, k_rX, k_rY, k_rZ\n');
    fprintf(lsm, '%%NODE.SUPPORT.SPRING\n');
    fprintf(lsm, '%d\n', nnsp);
    for n = 1:model.nnp
        if isempty(model.nodes(n).springStiff) == 0
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
                model.nodes(n).springStiff(1), model.nodes(n).springStiff(2),...
                model.nodes(n).springStiff(3), model.nodes(n).springStiff(4),...
                model.nodes(n).springStiff(5), model.nodes(n).springStiff(6));
        end    
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Initial dynamic conditions
if model.nnp > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify initial conditions for dynamic analysis\n');
%     fprintf(lsm, 'First line: Total number of nodes\n');
%     fprintf(lsm, 'Following lines: Node ID, dX, dY, dZ, rX, rY, rZ, vdX, vdY, vdZ, vrX, vrY, vrZ\n');
    fprintf(lsm, '%%NODE.INITIAL.CONDITIONS\n');
    fprintf(lsm, '%d\n', model.nnp);
    for n = 1:model.nnp
        iniCon = full(model.nodes(n).initCond);
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d     %d  %d  %d  %d  %d  %d\n', n,...
                iniCon(1,1), iniCon(2,1),...
                iniCon(3,1), iniCon(4,1),...
                iniCon(5,1), iniCon(6,1),...
                iniCon(1,2), iniCon(2,2),...
                iniCon(3,2), iniCon(4,2),...
                iniCon(5,2), iniCon(6,2));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Dynamic time functions (load variation patterns over time)

if ~isempty(model.timeFcns)
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify time functions for dynamic loads\n');
%     fprintf(lsm, 'First line: Total number of lists of time functions\n');
%     fprintf(lsm, 'Following lines: Function list (load pattern) id, number of time functions on this list\n');
%     fprintf(lsm, 'Function list name \n');
%     fprintf(lsm, '                 Function id, type*, ti [s], tf [s]\n');
%     fprintf(lsm, '                 Function input**\n');
%     fprintf(lsm, 'obs.*: type = 0 -> CONSTANT ; type = 1 -> HARMONIC ; type = 2 -> SLOPE ; type = 3 -> TABLE\n' );
%     fprintf(lsm, 'obs.**:CONSTANT  -> weight factor\n');
%     fprintf(lsm, '       HARMONIC -> weight factor, w [rad/s], phi [rad]\n');
%     fprintf(lsm, '       SLOPE    -> weight factor\n');
%     fprintf(lsm, '       TABLE    -> weight factor, binary file\n');
    fprintf(lsm, '%%LOAD.CASE.TIME.FUNCTION\n');
    fprintf(lsm, '%d\n', length(model.timeFcns));
    
    for n = 1:length(model.timeFcns)
        ptr =  model.timeFcns{n};
        [~,nFcn] = ptr.goThrough();
        
        fprintf(lsm, '%i %i\n',n,nFcn);
        fprintf(lsm, '%s\n',model.strTimeFcns{n});
        
         while ~isempty(ptr)
            fprintf(lsm, '%i  %i  %d  %d\n',...
                    ptr.id, ptr.type, ptr.ti, ptr.tf);
            switch ptr.type
                case STATIC
                    fprintf(lsm, '%d\n',ptr.weightFctr);
                    
                case PERIODIC
                    fprintf(lsm, '%d  %d  %d\n',...
                        ptr.weightFctr,ptr.w,ptr.phi);
                    
                case SLOPE
                    fprintf(lsm, '%d\n',ptr.weightFctr);
                    
                case TABLE
                     [~,file_name,ext] = fileparts(ptr.src_file);
                     fprintf(lsm, '%d %s\n',ptr.weightFctr,strcat(file_name,ext));                     
            end
            ptr = ptr.next;
        end
    end
    fprintf(lsm, '\n');
else
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'No time functions were defined for this model\n');
    fprintf(lsm, '%%LOAD.CASE.TIME.FUNCTION\n');
    fprintf(lsm, '%d\n\n', 0); 
end

%--------------------------------------------------------------------------
% Dynamic nodal loads
nndl = zeros(1,model.nnp);
for n = 1:model.nnp
    if ~isempty(model.nodes(n).load.dynamic)
        nndl(n) = n;
    end
end
nndl = nndl(nndl>0);

if ~isempty(nndl)
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify applied dynamic nodal loads\n');
%     fprintf(lsm, 'First line: Total number of nodes with applied dynamic load\n');
%     fprintf(lsm, 'Following lines: Node ID, fX [kN], fY [kN], fZ [kN], mX [kNm], mY [kNm], mZ [kNm], Function list (load pattern) id\n');
    fprintf(lsm, '%%LOAD.CASE.NODAL.FORCE.DYNAMIC\n');
    fprintf(lsm, '%d\n', length(nndl));
    
    for n = nndl
        fcn_id = 0;
        for ii =1:length(model.timeFcns)
            if model.nodes(n).load.getFcn() == model.timeFcns{ii}
                fcn_id = ii;
                break
            end
        end
        
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d    %i\n', n,...
            model.nodes(n).load.dynamic(1), model.nodes(n).load.dynamic(2),...
            model.nodes(n).load.dynamic(3), model.nodes(n).load.dynamic(4),...
            model.nodes(n).load.dynamic(5), model.nodes(n).load.dynamic(6),...
            fcn_id);
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Concentrated Nodal Mass
nncm = zeros(1,model.nnp);
for n = 1:model.nnp
    if ~isempty(model.nodes(n).load.dynamic)
        nncm(n) = n;
    end
end
nncm = nncm(nncm>0);

if ~isempty(nncm)
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify concentrated nodal mass\n');
%     fprintf(lsm, 'First line: Total number of nodes with concentrated nodal mass\n');
%     fprintf(lsm, 'Following lines: Node ID, concentrated nodal mass [kg]  \n');
    fprintf(lsm, '%%LOAD.CASE.NODAL.MASS.DYNAMIC\n');
    fprintf(lsm, '%d\n', length(nncm));
    
    for n = nncm
        if ~isempty(model.nodes(n).displMass)
            fprintf(lsm, '%d     %d\n', n, model.nodes(n).displMass*1000);
        end
    end    
   fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Material properties
if model.nmat > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide materials\n');
%     fprintf(lsm, 'First line: Total number of materials\n');
%     fprintf(lsm, 'Following lines: Material ID, Elasticity [MPa], Poisson Ratio, Thermal exp. coeff. [/oC], Density [kg/m�]\n');
    fprintf(lsm, '%%MATERIAL.ISOTROPIC\n');
    fprintf(lsm, '%d\n', model.nmat);
    for m = 1:model.nmat
        fprintf(lsm, '%d     %d  %d  %d  %d\n', m,...
                1e-3 * model.materials(m).elasticity,...
                model.materials(m).poisson,...
                model.materials(m).thermExp,...
                1e3  * model.materials(m).density);
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Cross-section properties
if model.nsec > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide cross-sections\n');
%     fprintf(lsm, 'First line: Total number of cross-sections\n');
%     fprintf(lsm, 'Following lines: Cross-section ID, Area_x [cm2], Area_y [cm2], Area_z [cm2]\n');
%     fprintf(lsm, '                 Inertia_x [cm4], Inertia_y [cm4], Inertia_z [cm4]\n');
%     fprintf(lsm, '                 Height_y [cm], Height_z [cm]\n');
    fprintf(lsm, '%%SECTION.PROPERTY\n');
    fprintf(lsm, '%d\n', model.nsec);
    for s = 1:model.nsec
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d  %d  %d\n', s,...
                1e4*model.sections(s).area_x,...
                1e4*model.sections(s).area_y,...
                1e4*model.sections(s).area_z,...
                1e8*model.sections(s).inertia_x,...
                1e8*model.sections(s).inertia_y,...
                1e8*model.sections(s).inertia_z,...
                1e2*model.sections(s).height_y,...
                1e2*model.sections(s).height_z);
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Beam end liberation
% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Flag specification for different beam end liberation conditions\n');
% fprintf(lsm, '1 - Continuous end / Continuous end\n');
% fprintf(lsm, '2 - Hinged end / Continuous end\n');
% fprintf(lsm, '3 - Continuous end / Hinged end\n');
% fprintf(lsm, '4 - Hinged end / Hinged end\n');
% fprintf(lsm, '5 - Continuous end / Semi-rigid end\n');
% fprintf(lsm, '6 - Semi-rigid end / Continuous end\n');
% fprintf(lsm, '7 - Hinged end / Semi-rigid end\n');
% fprintf(lsm, '8 - Semi-rigid end / Hinged end\n');
% fprintf(lsm, '9 - Semi-rigid end / Semi-rigid end\n');
fprintf(lsm, '%%BEAM.END.LIBERATION\n');
fprintf(lsm, '9\n');
fprintf(lsm, '1     1  1  1  1  1  1  1  1  1  1  1  1\n');
fprintf(lsm, '2     1  1  1  0  0  0  1  1  1  1  1  1\n');
fprintf(lsm, '3     1  1  1  1  1  1  1  1  1  0  0  0\n');
fprintf(lsm, '4     1  1  1  0  0  0  1  1  1  0  0  0\n');
fprintf(lsm, '5     1  1  1  1  1  1  1  1  1  2  2  2\n');
fprintf(lsm, '6     1  1  1  2  2  2  1  1  1  1  1  1\n');
fprintf(lsm, '7     1  1  1  0  0  0  1  1  1  2  2  2\n');
fprintf(lsm, '8     1  1  1  2  2  2  1  1  1  0  0  0\n');
fprintf(lsm, '9     1  1  1  2  2  2  1  1  1  2  2  2\n');
fprintf(lsm, '\n');

%--------------------------------------------------------------------------
% Element information
nsrj = 0; % initialize counter of semi-rigid joints
if model.nel > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide elements\n');
%     fprintf(lsm, 'First line: Total number of elements\n');
%     fprintf(lsm, 'Following lines: Element ID, Element type*, Element mass consideration**, Material ID, Cross-section ID, End lib. flag,\n');
%     fprintf(lsm, '                 Init. Node ID, Final Node ID, vz_X, vz_Y, vz_Z\n');
%     fprintf(lsm, '*Element types: 0 = Navier, 1 = Timoshenko\n');
%     fprintf(lsm, '**Element mass consideration: 0 = No, 1 = Yes\n');
    fprintf(lsm, '%%ELEMENT.BEAM\n');
    fprintf(lsm, '%d\n', model.nel);
    
    for e = 1:model.nel
        if (model.elems(e).hingei == CONTINUOUS_END) && (model.elems(e).hingef == CONTINUOUS_END)
            lib = 1;
        elseif (model.elems(e).hingei == HINGED_END) && (model.elems(e).hingef == CONTINUOUS_END)
            lib = 2;
        elseif (model.elems(e).hingei == CONTINUOUS_END) && (model.elems(e).hingef == HINGED_END)
            lib = 3;
        elseif (model.elems(e).hingei == HINGED_END) && (model.elems(e).hingef == HINGED_END)
            lib = 4;
        elseif (model.elems(e).hingei == CONTINUOUS_END) && (model.elems(e).hingef == SEMIRIGID_END)
            lib = 5;
            nsrj = nsrj + 1;
        elseif (model.elems(e).hingei == SEMIRIGID_END) && (model.elems(e).hingef == CONTINUOUS_END)
            lib = 6;
            nsrj = nsrj + 1;
        elseif (model.elems(e).hingei == HINGED_END) && (model.elems(e).hingef == SEMIRIGID_END)
            lib = 7;
            nsrj = nsrj + 1;
        elseif (model.elems(e).hingei == SEMIRIGID_END) && (model.elems(e).hingef == HINGED_END)
            lib = 8;
            nsrj = nsrj + 1;
        elseif (model.elems(e).hingei == SEMIRIGID_END) && (model.elems(e).hingef == SEMIRIGID_END)
            lib = 9;
            nsrj = nsrj + 2;
        end
        
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n', e,...
                model.elems(e).type,model.elems(e).mass_consideration,...
                model.elems(e).material.id, model.elems(e).section.id, lib,...
                model.elems(e).nodes(1).id, model.elems(e).nodes(2).id,...
                model.elems(e).vz(1), model.elems(e).vz(2), model.elems(e).vz(3));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Element unsolved intersection information
intersections = getappdata(0,'intersections');
if ~isempty(intersections)
    nints = size(intersections,2);
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide element intersections that are not characterized as nodes\n');
%     fprintf(lsm, 'First line: Total number of intersections\n');
%     fprintf(lsm, 'Following lines: 1st -> Intersection ID, coord_X [m], coord_Y [m], coord_Z [m], number of elements\n');
%     fprintf(lsm, '                 2nd -> Element identifiers\n');
    fprintf(lsm, '%%ELEMENT.BEAM.INTERSECTIONS\n');
    fprintf(lsm, '%d\n', nints);
    
    for ni = 1:nints
        fprintf(lsm, '%d     %d  %d  %d  %d\n', ni,...
                intersections(ni).coord(1),...
                intersections(ni).coord(2),...
                intersections(ni).coord(3),...
                size(intersections(ni).elems,2));
        fprintf(lsm, '     ');
        for elems_int = intersections(ni).elems
            fprintf(lsm, '  %d', elems_int);
        end
        fprintf(lsm, '\n');
    end
    fprintf(lsm, '\n');
end
clear intersections

%--------------------------------------------------------------------------
% Semi-rigid joint information
if model.nel > 0 && nsrj > 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Provide semi-rigid joints information\n');
%     fprintf(lsm, 'First line: Total number of semi-rigid joints\n');
%     fprintf(lsm, 'Following lines: Element ID, Joint 1*, Joint 2*, J1_krx, J1_kry, J1_krz, J2_krx, J2_kry, J2_krz\n');
%     fprintf(lsm, '*Joint 1, Joint 2: 0 = Continuous/hinged joint, 1 = Semi-rigid joint\n');
    fprintf(lsm, '%%ELEMENT.BEAM.JOINTS\n');
    fprintf(lsm, '%d\n', nsrj);
    
    for e = 1:model.nel
        if model.elems(e).hingei == SEMIRIGID_END || model.elems(e).hingef == SEMIRIGID_END
            j1 = model.elems(e).hingei * (model.elems(e).hingei - 1) / 2;
            j2 = model.elems(e).hingef * (model.elems(e).hingef - 1) / 2;
            
            if ~isempty(model.elems(e).kri)
                kri = model.elems(e).kri;
            else
                kri = zeros(1,3);
            end
            
            if ~isempty(model.elems(e).krf)
                krf = model.elems(e).krf;
            else
                krf = zeros(1,3);
            end
            
            fprintf(lsm, '%i     %i  %i  %d  %d  %d  %d  %d  %d \n',...
                    e, j1, j2,...
                    kri(1), kri(2), kri(3),...
                    krf(1), krf(2), krf(3));
        end
    end
    fprintf(lsm, '\n');
end

% %--------------------------------------------------------------------------
% % Nodal prescribed displacements
% nnpd = 0;
% for n = 1:model.nnp
%     if isempty(model.nodes(n).prescDispl) == 0
%         nnpd = nnpd + 1;
%     end
% end
% 
% if nnpd ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify nodal prescribed displacements\n');
%     fprintf(lsm, 'First line: Total number of nodes with prescribed displacement\n');
%     fprintf(lsm, 'Following lines: Node ID, dX [mm], dY [mm], dZ [mm], rX [rad], rY [rad], rZ [rad]\n');
%     fprintf(lsm, '%%LOAD.CASE.NODAL.DISPLACEMENT\n');
%     fprintf(lsm, '%d\n', nnpd);
%     for n = 1:model.nnp
%         if isempty(model.nodes(n).prescDispl) == 0
%             fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
%                     1e3*model.nodes(n).prescDispl(1), 1e3*model.nodes(n).prescDispl(2),...
%                     1e3*model.nodes(n).prescDispl(3),     model.nodes(n).prescDispl(4),...
%                         model.nodes(n).prescDispl(5),     model.nodes(n).prescDispl(6));
%         end
%     end
%     fprintf(lsm, '\n');
% end

%--------------------------------------------------------------------------
% Load cases

% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Specify load cases\n');
% fprintf(lsm, 'First line: Total number of load cases\n');
% fprintf(lsm, 'Following lines: Load case labels\n');
fprintf(lsm, '%%LOAD\n');
fprintf(lsm, '%d\n', model.nlc);
for i = 1:model.nlc
    fprintf(lsm, '%s\n', char(model.strLc(i,:)));
end
fprintf(lsm, '\n');

%--------------------------------------------------------------------------
% Load case combinations
if isempty(model.loadComb) == 0
    ncomb = size(model.loadComb,2);
else
    ncomb = 0;
end    

if ncomb ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify load case combinations\n');
%     fprintf(lsm, 'First line: Total number of load case combinations\n');
%     fprintf(lsm, 'Following lines: Combination label\n');
%     fprintf(lsm, '                 Comb factors\n');
    fprintf(lsm, '%%LOAD.COMBINATION\n');
    fprintf(lsm, '%d\n',ncomb);
    for i = 1:ncomb
        fprintf(lsm, '%s\n',char(model.strComb(i,:)));
        for j = 1:model.nlc
            if j ~= model.nlc
                fprintf(lsm, '%f  ', model.loadComb(j,i));
            else
                fprintf(lsm, '%f\n', model.loadComb(j,i));
            end    
        end 
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Current load case

% fprintf(lsm, '-----------------------------------------------------------------------\n');
% fprintf(lsm, 'Specify current load case/combination\n');
% fprintf(lsm, 'Following line: Current load case/combination ID\n');
% fprintf(lsm, 'ATTENTION: To set a combination as the current load case,\n');
% fprintf(lsm, 'change the ID number to (comb ID + number of load cases)\n');
fprintf(lsm, '%%LOAD.CASE\n');
fprintf(lsm, '%d\n', 1);  % The first load case is saved as the current one
fprintf(lsm, '\n');
    
%--------------------------------------------------------------------------
% Nodal load cases
nnlc = 0;
for n = 1:model.nnp
    allLogic = all(model.nodes(n).nodalLoadCase == 0);
    if isempty(model.nodes(n).nodalLoadCase) == 0
        if ~all(allLogic == 1)
            nnlc = nnlc + 1;
        end    
    end
end

if nnlc ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify nodal load cases\n');
%     fprintf(lsm, 'First line: Total number of nodes with applied load and/or prescribed displacements\n');
%     fprintf(lsm, 'in at least one load case\n');
%     fprintf(lsm, 'Following lines: Node ID, fX [kN], fY [kN], fZ [kN], mX [kNm], mY [kNm], mZ [kNm],\n');
%     fprintf(lsm, 'dx[mm], dy[mm], dz[mm], rx[rad], ry[rad], rz[rad], number of load cases\n');
    fprintf(lsm, '%%LOAD.CASE.NODAL.CASES\n');
    fprintf(lsm, '%d\n', nnlc);
    fprintf(lsm, '\n');
    for n = 1:model.nnp
        allLogic = all(model.nodes(n).nodalLoadCase == 0);
        if isempty(model.nodes(n).nodalLoadCase) == 0
            if ~all(allLogic == 1)
                fprintf(lsm, '%d     ', n);
                for i = 1:size(model.nodes(n).nodalLoadCase,2)
                    if i == 1
                        fprintf(lsm, '%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d      %d\n', ...
                                model.nodes(n).nodalLoadCase(1,i), model.nodes(n).nodalLoadCase(2,i),...
                                model.nodes(n).nodalLoadCase(3,i), model.nodes(n).nodalLoadCase(4,i),...
                                model.nodes(n).nodalLoadCase(5,i), model.nodes(n).nodalLoadCase(6,i),...
                                1e3*model.nodes(n).nodalLoadCase(7,i), 1e3*model.nodes(n).nodalLoadCase(8,i),...
                                1e3*model.nodes(n).nodalLoadCase(9,i), model.nodes(n).nodalLoadCase(10,i),...
                                model.nodes(n).nodalLoadCase(11,i), model.nodes(n).nodalLoadCase(12,i),...
                                size(model.nodes(n).nodalLoadCase,2));
                    else
                        fprintf(lsm, '      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n', ...
                                model.nodes(n).nodalLoadCase(1,i), model.nodes(n).nodalLoadCase(2,i),...
                                model.nodes(n).nodalLoadCase(3,i), model.nodes(n).nodalLoadCase(4,i),...
                                model.nodes(n).nodalLoadCase(5,i), model.nodes(n).nodalLoadCase(6,i),...
                                1e3*model.nodes(n).nodalLoadCase(7,i), 1e3*model.nodes(n).nodalLoadCase(8,i),...
                                1e3*model.nodes(n).nodalLoadCase(9,i), model.nodes(n).nodalLoadCase(10,i),...
                                model.nodes(n).nodalLoadCase(11,i), model.nodes(n).nodalLoadCase(12,i));
                    end    
                end    
                fprintf(lsm, '\n');
            end    
        end
    end
    fprintf(lsm, '\n');
end

% %--------------------------------------------------------------------------
% % Nodal loads
% nnnl = 0;
% for n = 1:model.nnp
%     if isempty(model.nodes(n).load.static) == 0
%         nnnl = nnnl + 1;
%     end
% end
% 
% if nnnl ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify applied nodal loads\n');
%     fprintf(lsm, 'First line: Total number of nodes with applied load\n');
%     fprintf(lsm, 'Following lines: Node ID, fX [kN], fY [kN], fZ [kN], mX [kNm], mY [kNm], mZ [kNm]\n');
%     fprintf(lsm, '%%LOAD.CASE.NODAL.FORCE\n');
%     fprintf(lsm, '%d\n', nnnl);
%     for n = 1:model.nnp
%         if isempty(model.nodes(n).load.static) == 0
%             fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
%                     model.nodes(n).load.static(1), model.nodes(n).load.static(2),...
%                     model.nodes(n).load.static(3), model.nodes(n).load.static(4),...
%                     model.nodes(n).load.static(5), model.nodes(n).load.static(6));
%         end
%     end
%     fprintf(lsm, '\n');
% end

%--------------------------------------------------------------------------
% Element load cases
nelc = 0;
for e = 1:model.nel
    if isempty(model.elems(e).load.elemLoadCase) == 0
        nelc = nelc + 1;   
    end
end

if nelc ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify element load cases\n');
%     fprintf(lsm, 'First line: Total number of elements with distributed and/or thermal loads in\n');
%     fprintf(lsm, 'at least one load case\n');
%     fprintf(lsm, 'Following lines: Element ID, UnifDir, qX [kN/m], qY [kN/m], qZ [kN/m], LinearDir,\n');
%     fprintf(lsm, 'qX1 [kN/m], qY1 [kN/m], qZ1 [kN/m], qX2[kN/m], qY2[kN/m], qZ2[kN/m], dt_x[oC],\n'); 
%     fprintf(lsm, 'dt_y[oC], dt_z[oC], number of load cases\n');
%     fprintf(lsm, 'Obs.: Load direction (0->Global, 1->Local)\n');
    fprintf(lsm, '%%LOAD.CASE.BEAM.CASES\n');
    fprintf(lsm, '%d\n', nelc);
    fprintf(lsm, '\n');
    for e = 1:model.nel
        if isempty(model.elems(e).load.elemLoadCase) == 0
            fprintf(lsm, '%d     ', e);
            for i = 1:size(model.elems(e).load.elemLoadCase,2)
                if i == 1
                    fprintf(lsm, '%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d      %d\n', ...
                            model.elems(e).load.elemLoadCase(1,i), model.elems(e).load.elemLoadCase(2,i),...
                            model.elems(e).load.elemLoadCase(3,i), model.elems(e).load.elemLoadCase(4,i),...
                            model.elems(e).load.elemLoadCase(5,i), model.elems(e).load.elemLoadCase(6,i),...
                            model.elems(e).load.elemLoadCase(7,i), model.elems(e).load.elemLoadCase(8,i),...
                            model.elems(e).load.elemLoadCase(9,i), model.elems(e).load.elemLoadCase(10,i),...
                            model.elems(e).load.elemLoadCase(11,i), model.elems(e).load.elemLoadCase(12,i),...
                            model.elems(e).load.elemLoadCase(13,i), model.elems(e).load.elemLoadCase(14,i),...
                            size(model.elems(e).load.elemLoadCase,2));
                else
                    fprintf(lsm, '      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n', ...
                            model.elems(e).load.elemLoadCase(1,i), model.elems(e).load.elemLoadCase(2,i),...
                            model.elems(e).load.elemLoadCase(3,i), model.elems(e).load.elemLoadCase(4,i),...
                            model.elems(e).load.elemLoadCase(5,i), model.elems(e).load.elemLoadCase(6,i),...
                            model.elems(e).load.elemLoadCase(7,i), model.elems(e).load.elemLoadCase(8,i),...
                            model.elems(e).load.elemLoadCase(9,i), model.elems(e).load.elemLoadCase(10,i),...
                            model.elems(e).load.elemLoadCase(11,i), model.elems(e).load.elemLoadCase(12,i),...
                            model.elems(e).load.elemLoadCase(13,i), model.elems(e).load.elemLoadCase(14,i));
                end    
            end 
            fprintf(lsm, '\n');
        end
    end
    fprintf(lsm, '\n');
end

% %--------------------------------------------------------------------------
% % Element uniformly distributed loads
% neul = 0;
% for e = 1:model.nel
%     if isempty(model.elems(e).load.uniformGbl) == 0
%         neul = neul + 1;
%     end
% end
% 
% if neul ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify element uniformly distributed loads\n');
%     fprintf(lsm, 'First line: Total number of elements with uniformly distributed load\n');
%     fprintf(lsm, 'Following pair of lines:\n');
%     fprintf(lsm, '                         Element ID, Load direction (0->Global, 1->Local)\n');
%     fprintf(lsm, '                         Qx [kN/m], Qy [kN/m], Qz [kN/m]\n');
%     fprintf(lsm, '%%LOAD.CASE.BEAM.UNIFORM\n');
%     fprintf(lsm, '%d\n', neul);
%     for e = 1:model.nel
%         if isempty(model.elems(e).load.uniformGbl) == 0
%             fprintf(lsm, '%d     %d\n', e, model.elems(e).load.uniformDir);
%             
%             if model.elems(e).load.uniformDir == 0
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.uniformGbl(1),...
%                         model.elems(e).load.uniformGbl(2),...
%                         model.elems(e).load.uniformGbl(3));
%                 
%             elseif model.elems(e).load.uniformDir == 1
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.uniformLcl(1),...
%                         model.elems(e).load.uniformLcl(2),...
%                         model.elems(e).load.uniformLcl(3));
%             end
%         end
%     end
%     fprintf(lsm, '\n');
% end
% 
% %--------------------------------------------------------------------------
% % Element linearly distributed loads
% nell = 0;
% for e = 1:model.nel
%     if isempty(model.elems(e).load.linearGbl) == 0
%         nell = nell + 1;
%     end
% end
% 
% if nell ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify element linearly distributed loads\n');
%     fprintf(lsm, 'First line: Total number of elements with linearly distributed load\n');
%     fprintf(lsm, 'Following triple of lines:\n');
%     fprintf(lsm, '                         Element ID, Load direction (0->Global, 1->Local)\n');
%     fprintf(lsm, '                         Qx_init [kN/m], Qy_init [kN/m], Qz_init [kN/m]\n');
%     fprintf(lsm, '                         Qx_final [kN/m], Qy_final [kN/m], Qz_final [kN/m]\n');
%     fprintf(lsm, '%%LOAD.CASE.BEAM.LINEAR\n');
%     fprintf(lsm, '%d\n', nell);
%     for e = 1:model.nel
%         if isempty(model.elems(e).load.linearGbl) == 0
%             fprintf(lsm, '%d     %d\n', e, model.elems(e).load.linearDir);
%             
%             if model.elems(e).load.linearDir == 0
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.linearGbl(1),...
%                         model.elems(e).load.linearGbl(2),...
%                         model.elems(e).load.linearGbl(3));
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.linearGbl(4),...
%                         model.elems(e).load.linearGbl(5),...
%                         model.elems(e).load.linearGbl(6));
%                 
%             elseif model.elems(e).load.linearDir ==1
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.linearLcl(1),...
%                         model.elems(e).load.linearLcl(2),...
%                         model.elems(e).load.linearLcl(3));
%                 fprintf(lsm, '%d  %d  %d\n',...
%                         model.elems(e).load.linearLcl(4),...
%                         model.elems(e).load.linearLcl(5),...
%                         model.elems(e).load.linearLcl(6));
%             end
%         end
%     end
%     fprintf(lsm, '\n');
% end
% 
% %--------------------------------------------------------------------------
% % Element temperature variations
% netv = 0;
% for e = 1:model.nel
%     dtx = model.elems(e).load.tempVar_X;
%     dty = model.elems(e).load.tempVar_Y;
%     dtz = model.elems(e).load.tempVar_Z;
%     if (dtx ~= 0) || (dty ~= 0) || (dtz ~= 0)
%         netv = netv + 1;
%     end
% end
% 
% if netv ~= 0
%     fprintf(lsm, '-----------------------------------------------------------------------\n');
%     fprintf(lsm, 'Specify element thermal loads\n');
%     fprintf(lsm, 'First line: Total number of elements with thermal load\n');
%     fprintf(lsm, 'Following lines: Element ID, dt_x [oC], dt_y [oC], dt_z [oC]\n');
%     fprintf(lsm, 'dt_x       --> element temperature variation on its center of gravity axis\n');
%     fprintf(lsm, 'dt_y, dt_z --> element temperature gradient relative to local axes Y and Z\n');
%     fprintf(lsm, '               (bottomFaceTempVar - topFaceTempVar)\n');
%     fprintf(lsm, '%%LOAD.CASE.BEAM.TEMPERATURE\n');
%     fprintf(lsm, '%d\n', netv);
%     for e = 1:model.nel
%         dtx = model.elems(e).load.tempVar_X;
%         dty = model.elems(e).load.tempVar_Y;
%         dtz = model.elems(e).load.tempVar_Z;
%         if (dtx ~= 0) || (dty ~= 0) || (dtz ~= 0)
%             fprintf(lsm, '%d     %d  %d  %d\n', e, dtx, dty, dtz);
%         end
%     end
%     fprintf(lsm, '\n');
% end
% 
%--------------------------------------------------------------------------
fprintf(lsm, '%%END');
end