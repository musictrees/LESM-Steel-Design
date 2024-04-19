%% Auxiliary Modelling Function
%
% This file contains functions that are called to do modelling tasks that
% are not dependent of an Emouse object.
%
% Input:
% -> whichFunction: string that indicates which function is being called
% -> fctnArgIn: 'whichFunction' input arguments
%
%% ------------------------------------------------------------------------
% Works as a switch, calls other functions.
function fctnArgOut = auxModelFctn(whichFunction,fctnArgIn)
    fctnArgOut = [];
    switch whichFunction
        case 'isPointInNode'
            fctnArgOut = isPointInNode(fctnArgIn);
        case 'isPointInElem'
            fctnArgOut = isPointInElem(fctnArgIn);
        case 'isPointInLine3D'
            fctnArgOut = isPointInLine3D(fctnArgIn);
        case 'isElemEqual'
            fctnArgOut = isElemEqual(fctnArgIn);
        case 'areLinesCrossed3D'
            fctnArgOut = areLinesCrossed3D(fctnArgIn);
        case 'areElemsCrossed'
            fctnArgOut = areElemsCrossed(fctnArgIn);
        case 'getCrossIntSectPoints'
            fctnArgOut = getCrossIntSectPoints(fctnArgIn);
        case 'getCrossNodePoints'
            [crN,cE,eC] = getCrossNodePoints(fctnArgIn);
            fctnArgOut = {crN,cE,eC};
        case 'getCrossElemPoints'
            [crPts,newNodes] = getCrossElemPoints(fctnArgIn);
            fctnArgOut = {crPts,newNodes};
        case 'getNewNodes'
            fctnArgOut = getNewNodes(fctnArgIn);
        case 'getElemPointDisplAndStress'
            fctnArgOut = getElemPointDisplAndStress(fctnArgIn);
        case 'divideElement'
            divideElement(fctnArgIn)
        case 'deleteNodes'
            deleteNodes(fctnArgIn)
        case 'deleteElems'
            deleteElems(fctnArgIn)
    end
end
%% ------------------------------------------------------------------------
function whichNode = isPointInNode(coords)
    % Initialize variable to be returned
    whichNode = [];
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Get number of existing nodes
    nnp = getappdata(0,'nnp');
    
    for n = 1:nnp
        if coords(1) == nodes(n).coord(1) && coords(2) == nodes(n).coord(2)...
           && coords(3) == nodes(n).coord(3)
            whichNode = n;
            break
        end
    end
end

%--------------------------------------------------------------------------
% Checks if current position of cursor is inside any element
% Returns element id. If whichElem = 0, point is not inside any element.
function whichElems = isPointInElem(inputArg)
    % Check if input is cell to determine graphic tolerance
    if iscell(inputArg)
        coords = inputArg{1};
        tol = inputArg{2};
    else
        coords = inputArg;
        tol = 1e-15;
    end
    
    % Get number of elements
    nel = getappdata(0,'nel');
    
    % Initialize flag
    whichElems = zeros(1,nel);
    
    % Check if there are elements
    if nel ~= 0
        
        % Get vector of handles to elem objects
        elems = getappdata(0,'elems');
        
        for e = 1:nel
            % Get element end coordinates
            xe1 = elems(e).nodes(1).coord(1);
            ye1 = elems(e).nodes(1).coord(2);
            ze1 = elems(e).nodes(1).coord(3);
            xe2 = elems(e).nodes(2).coord(1);
            ye2 = elems(e).nodes(2).coord(2);
            ze2 = elems(e).nodes(2).coord(3);
            elemCoords = [xe1 ye1 ze1;
                          xe2 ye2 ze2];
            
            % Get flag for point inside element 'e'
            if coords(3) == ze1 && coords(3) == ze2
                inElemFlag = isPointInThisElem(coords,elemCoords,tol,'2D');
            else
                inElemFlag = isPointInThisElem(coords,elemCoords,tol,'3D');
            end
            
            % If point is inside an element 'e', update whichElems
            if inElemFlag == true
                whichElems(e) = e;
            end
        end
        whichElems = nonzeros(whichElems)';
    end
end

%--------------------------------------------------------------------------
function pointInFlag = isPointInLine3D(fctnArgIn)
    pointCoords = fctnArgIn{1};
    lineCoords = fctnArgIn{2};
    tol = fctnArgIn{3};
    if size(fctnArgIn,2) == 4
        clickFlag = fctnArgIn{4};
    else
        clickFlag = [];
    end
    if isempty(clickFlag)
        pointInFlag = isPointInThisElem(pointCoords,lineCoords,tol,'3D');
    else
        pointInFlag = isPointInThisElem(pointCoords,lineCoords,tol,'3D',clickFlag);
    end
end

%--------------------------------------------------------------------------
function whichElem = isElemEqual(inputArg)
    ni = inputArg{1}; % handle to initial node object
    nf = inputArg{2}; % handle to final node object
    
    % Initialize variable to be returned
    whichElem = 0;
    
    % Get number of elements and vector of handles to elem objects
    nel = getappdata(0,'nel');
    elems = getappdata(0,'elems');
    
    % Check if nodal coordinates are coincident with an existing element
    if nel ~= 0
        for e = 1:nel
            elemNodes = [elems(e).nodes(1).id elems(e).nodes(2).id];
            if (elemNodes(1) == ni.id && elemNodes(2) == nf.id) ||...
               (elemNodes(1) == nf.id && elemNodes(2) == ni.id)

                whichElem = e;
                break
            end
        end
    end
end

%--------------------------------------------------------------------------
function output = areLinesCrossed3D(fctnArgIn)
    line_1 = fctnArgIn{1}; % array of first line end coordinates
    line_2 = fctnArgIn{2}; % array of second line end coordinates
    tol = fctnArgIn{3}; % graphic tolerance
    
    flag = false; % Initialize flag to be returned
    
    line_3 = perpLine3D(line_1,line_2); % Returns array of ortho line
    if ~isempty(line_3)
        if norm(line_3(2,:)-line_3(1,:)) <= tol
            flag = true;
        end
    end
    
    output = {flag,line_3};
end

%--------------------------------------------------------------------------
% Works as a switch to call areElemsCrossed2D or areElemsCrossed3D
function crossPoint = areElemsCrossed(fctnArgIn)
    coords = fctnArgIn{1}; % new elem coordinates
    e = fctnArgIn{2};      % existing elem id
    
    sz = size(coords,2);
    
    switch sz
        case 2
            crossPoint = areElemsCrossed2D(coords,e);
            if ~isempty(crossPoint)
                crossPoint(3) = 0;  % z = 0
            end
        case 3
            tol = fctnArgIn{3};    % graphic tolerance
            crossPoint = areElemsCrossed3D(coords,e,tol);
    end
end

%--------------------------------------------------------------------------
% Check if an inserted element is crossing (or is collinear with) other 
% element. Returns corrdinates of crossing point.
function crossPoint = areElemsCrossed2D(coords,e)
    % Initialize variable to be returned
    crossPoint = [];
    
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Get existing elem end coordinates
    e1_nodei_x = elems(e).nodes(1).coord(1);
    e1_nodei_y = elems(e).nodes(1).coord(2);
    e1_nodef_x = elems(e).nodes(2).coord(1);
    e1_nodef_y = elems(e).nodes(2).coord(2);
    
    % Get new elem end coordinates
    e2_nodei_x = coords(1,1);
    e2_nodei_y = coords(1,2);
    e2_nodef_x = coords(2,1);
    e2_nodef_y = coords(2,2);
    
    % Define triangle points
    a = [e1_nodei_x, e1_nodei_y];
    b = [e1_nodef_x, e1_nodef_y];
    c = [e2_nodei_x, e2_nodei_y];
    d = [e2_nodef_x, e2_nodef_y];
    
    % Get signed areas
    orientABC = orient2D(a,b,c);
    orientABD = orient2D(a,b,d);
    orientCDA = orient2D(c,d,a);
    orientCDB = orient2D(c,d,b);
    
    % Check if elements cross
    if orientCDA > 0 && orientCDB > 0 % new element is over existing element
        return
    elseif orientCDA < 0 && orientCDB < 0 % new element is under existing element
        return
    elseif orientABC > 0 && orientABD > 0 % existing element is over new element
        return
    elseif orientABC < 0 && orientABD < 0 % existing element is under new element
        return
    elseif (orientABC == 0 && orientABD == 0) || (orientCDA == 0 && orientCDB == 0) % elements are collinear
        return
    elseif orientCDA == 0 && orientCDB ~= 0 % new element initial node touches existing element
        % crossPoint = a;
        return
    elseif orientCDA ~= 0 && orientCDB == 0 % new element final node touches existing element
        % crossPoint = b;
        return
    elseif orientABC == 0 && orientABD ~= 0 % existing element initial node touches new element
        % crossPoint = c;
        return
    elseif orientABC ~= 0 && orientABD == 0 % existing element final node touches new element
        % crossPoint = d;
        return
    else % elements cross
        crossPoint = c + (orientABC / (orientABC - orientABD)) * (d - c);
    end
end

%--------------------------------------------------------------------------
% Check if an inserted element is crossing (or is collinear with) other 
% element. Returns corrdinates of crossing point.
% Input:
% -> coords: new element end coordinates
% -> e: existing element id
function crossPoint = areElemsCrossed3D(coords,e,tol)
    % Initialize variable to be returned
    crossPoint = [];
    
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Get existing elem end coordinates
    elemCoords = [elems(e).nodes(1).coord;
                  elems(e).nodes(2).coord];
    
    % Check if elements are crossed
    fctnArgOut = areLinesCrossed3D({coords,elemCoords,tol});
    elemsAreCrossed = fctnArgOut{1};
    perpLine = fctnArgOut{2};
    if elemsAreCrossed == true
        crossPoint = perpLine(2,:);
    end
end
%--------------------------------------------------------------------------
function crossInterSections = getCrossIntSectPoints(fctnArgIn)
    % Check if a numeric tolerance was given or needs to be defined
    if iscell(fctnArgIn)
        coords = fctnArgIn{1};
        tol = fctnArgIn{2};
    else
        coords = fctnArgIn;
        tol = [];
        % Get axisWidth
        mdata = guidata(findobj('Tag','GUI_Main'));
        dfltUnits = get(mdata.axes_Canvas,'units');
        set(mdata.axes_Canvas,'units','normalized');
        limits = get(mdata.axes_Canvas,'Position');
        set(mdata.axes_Canvas,'units',dfltUnits);
        axisWidth = limits(3);
    end
    
    % Get vector of handles to elemIntersection objects and number of those
    intersections = getappdata(0,'intersections');
    nis = size(intersections,2);
    
    % Initialize matrix to be returned
    crossInterSections = zeros(nis,5);
    
    % Initialize auxiliar vector that indicates which lines of
    % crossInterSections will be deleted before its returned.
    delCrossInterSections = zeros(1,nis);
    
    % Check if model is 2D or 3D
    sz = size(coords,2);
    switch sz
        case 2
            planeOrSpatial = '2D';
            if isempty(tol)
                % Adjust axisWidth parameter based on axis limits
                aux = [mdata.axes_Canvas.XLim;mdata.axes_Canvas.YLim];
                scl = max([diff(aux(1,:)),diff(aux(2,:))])/2;
                axisWidth = 0.75*axisWidth * 10^(floor(log10(scl)));
                % Get numeric tolerance
                tol = axisWidth/18;
            end
        case 3
            planeOrSpatial = '3D';
            if isempty(tol)
                % Get numeric tolerance
                tol = axisWidth/20;
            end
    end
    
    count = 0;
    for n = 1:nis
        % Get intersection coordinates
        intSectCoords = intersections(n).coord;

        % Check if intersection 'n' is inside new element
        intSectInFlag = isPointInThisElem(intSectCoords,coords,tol,planeOrSpatial);
        
        % Check if new element does not cross intersection 'n'
        if intSectInFlag == false
            count = count + 1;
            delCrossInterSections(count) = n;
        else % new element crosses intersection 'n'
            r = norm(intersections(n).coord(1:sz) - coords(1,:));
            crossInterSections(n,:) = [n intersections(n).coord r];
        end
    end

    % Check if all lines of crossInterSections need to be deletd (new element
    % does not cross any existing intersection).
    if all(delCrossInterSections ~= 0)
        crossInterSections = [];
        return
    end
    
    % Delete all unnecessary lines from crossInterSections matrix.
    for dci = 1:size(delCrossInterSections,2)
        if delCrossInterSections(dci) == 0
            delCrossInterSections(dci:end) = [];
            crossInterSections(delCrossInterSections,:) = [];
            break
        end
    end
    
    % Reorder crossInterSections based on distance to elem node_i
    dist = crossInterSections(:,5);
    crossInterSections(:,5) = [];
    reorder = zeros(1,size(crossInterSections,1));
    for ci = 1:size(crossInterSections,1)
        [~,minIndex] = min(dist);
        reorder(ci) = minIndex;
        dist(minIndex) = 1.1 * max(dist);
    end
    crossInterSections = crossInterSections(reorder,:);
end

%--------------------------------------------------------------------------
% Check if new element crosses existing nodes.
% Returns: 
% -> Matrix indicating eventual cross points, where each line yields
% [node_id, crossPoint_x, crossPoint_y];
% -> Vector of collinear elements;
% -> Matrix that indicates new elements connectivity
function [crossNodes,collinearElems,elemConnect] = getCrossNodePoints(fctnArgIn)
    % Check if a numeric tolerance was given or needs to be defined
    if iscell(fctnArgIn)
        coords = fctnArgIn{1};
        tol = fctnArgIn{2};
    else
        coords = fctnArgIn;
        tol = [];
        % Get axisWidth
        mdata = guidata(findobj('Tag','GUI_Main'));
        dfltUnits = get(mdata.axes_Canvas,'units');
        set(mdata.axes_Canvas,'units','normalized');
        limits = get(mdata.axes_Canvas,'Position');
        set(mdata.axes_Canvas,'units',dfltUnits);
        axisWidth = limits(3);
    end
    
    % Get number of nodes
    nnp = getappdata(0,'nnp');
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Initialize matrix to be returned (collinearPoints)
    crossNodes = zeros(nnp,5);
    
    % Initialize auxiliar variable that indicates which elements are
    % collinear with the new one
    collinearElems = [];
    
    % Initialize elem connect
    elemConnect = [];
    
    % Initialize auxiliar vector that indicates which lines of crossNodes
    % will be deleted before its returned.
    delCrossNodes = zeros(1,nnp);
    
    % Check if model is 2D or 3D
    sz = size(coords,2);
    switch sz
        case 2
            planeOrSpatial = '2D';
            if isempty(tol)
                % Adjust axisWidth parameter based on axis limits
                aux = [mdata.axes_Canvas.XLim;mdata.axes_Canvas.YLim];
                scl = max([diff(aux(1,:)),diff(aux(2,:))])/2;
                axisWidth = 0.75*axisWidth * 10^(floor(log10(scl)));
                % Get numeric tolerance
                tol = axisWidth/18;
            end
        case 3
            planeOrSpatial = '3D';
            if isempty(tol)
                % Get numeric tolerance
                tol = axisWidth/20;
            end
    end
    
    count = 0;
    for n = 1:nnp
        % Get node coordinates
        nodeCoords = [nodes(n).coord(1), nodes(n).coord(2), nodes(n).coord(3)];

        % Check if node 'n' is inside new element
        nodeInFlag = isPointInThisElem(nodeCoords,coords,tol,planeOrSpatial);
        
        % If new element does not cross node 'n'
        if nodeInFlag == false
            count = count + 1;
            delCrossNodes(count) = n;
        else % new element crosses node 'n'
            r = norm(nodes(n).coord(1:sz) - coords(1,:));
            crossNodes(n,:) = [n nodes(n).coord r];
        end
    end

    % Check if all lines of crossNodes need to be deletd (new element
    % does not cross any existing node).
    if all(delCrossNodes ~= 0)
        crossNodes = [];
        return
    end
    
    % Delete all unnecessary lines from crossNodes matrix.
    for dcp = 1:size(delCrossNodes,2)
        if delCrossNodes(dcp) == 0
            delCrossNodes(dcp:end) = [];
            crossNodes(delCrossNodes,:) = [];
            break
        end
    end
    
    % Reorder crossNodes based on distance to elem node_i
    dist = crossNodes(:,5);
    crossNodes(:,5) = [];
    reorder = zeros(1,size(crossNodes,1));
    for cp = 1:size(crossNodes,1)
        [~,minIndex] = min(dist);
        reorder(cp) = minIndex;
        dist(minIndex) = 1.1 * max(dist);
    end
    crossNodes = crossNodes(reorder,:);
    
    % Check if there are collinear elements and set new elements connect
    collinearElems = zeros(1,getappdata(0,'nel'));
    countN = 0;
    countE = 0;
    elemConnect = zeros(size(crossNodes,1)-1,2);
    for n = 1:(size(crossNodes,1)-1)
        ni = nodes(crossNodes(n,1));
        nf = nodes(crossNodes(n+1,1));
        existentElemFlag = isElemEqual({ni,nf});
        if existentElemFlag == false
            countN = countN + 1;
            elemConnect(countN,:) = [ni.id nf.id];
        else
            countE = countE + 1;
            collinearElems(countE) = existentElemFlag;
        end
    end

    collinearElems = nonzeros(collinearElems)';

    for i = 1:size(elemConnect,1)
        if all(elemConnect(i,:) == 0)
            elemConnect(i:end,:) = [];
            break
        end
    end
end

%--------------------------------------------------------------------------
% Check if new element crosses other elements.
% Returns:
% -> Matrix indicating eventual cross points, where each line yields
% [element_id, crossPoint_x, crossPoint_y];
% -> Nodes that need to be created.
function [crossPoints,newNodes] = getCrossElemPoints(inputArg)
    coords = inputArg{1};
    crossNodePoints = inputArg{2};
    collinearElems = inputArg{3};
    if size(inputArg,2) == 4
        tol = inputArg{4};
    else
        tol = 1e-10;
    end
    
    % Get coords array size (2 == 2D, 3 == 3D)
    sz_coords = size(coords,2);
    
    % Get number of existing nodes and vector of handles to node objects
    nnp = getappdata(0,'nnp');
    nodes = getappdata(0,'nodes');
    
    % Get number of existing elements (not counting the new element yet)
    nel = getappdata(0,'nel');
    
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Initialize matrix to be returned (crossPoints)
    crossPoints = zeros(nel,5);
    
    % Initialize return variable
    newNodes = zeros(size(crossPoints,1) - size(collinearElems,2),3);
    
    % Initialize auxiliar vector that indicates which lines of crossPoints
    % will be deleted before its returned.
    delCrossPoints = zeros(1,nel);
    
    count = 0;
    for e = 1:nel
        % Get cross point between new element and element 'e'
        crPt = areElemsCrossed({coords,e,tol});
        
        % Check if new element does not cross element 'e'
        if isempty(crPt)
            count = count + 1;
            delCrossPoints(count) = e;
        else % new element crosses element 'e'
            existentNodeFlag = false;
            for n = 1:nnp
                if norm(nodes(n).coord-crPt) <= tol
                    existentNodeFlag = true;
                    break
                end
            end
            if existentNodeFlag == false
                r = norm(crPt(1:sz_coords) - coords(1,:));
                crossPoints(e,:) = [e crPt r];
            else
                count = count + 1;
                delCrossPoints(count) = e;
            end
        end
    end

    % Check if all lines of crossPoints need to be deletd (new element does
    % not cross any existing element).
    if all(delCrossPoints ~= 0)
        crossPoints = [];
    else
        % Delete all unnecessary lines from crossPoints matrix.
        for dcp = 1:size(delCrossPoints,2)
            if delCrossPoints(dcp) == 0
                delCrossPoints(dcp:end) = [];
                crossPoints(delCrossPoints,:) = [];
                break
            end
        end
    end
    
    % Reorder crossPoints based on distance to elem node_i
    if ~isempty(crossPoints)
        dist = crossPoints(:,5);
        crossPoints(:,5) = [];
        reorder = zeros(1,size(crossPoints,1));
        for cp = 1:size(crossPoints,1)
            [~,minIndex] = min(dist);
            reorder(cp) = minIndex;
            dist(minIndex) = 1.1 * max(dist);
        end
        crossPoints = crossPoints(reorder,:);
    end
    
    % Get nodes to be created (discard redundant nodes)
    countNN = 0;
    for cep = 1:size(crossPoints,1)
        e = crossPoints(cep,1);
        if all(collinearElems ~= e)
            redundantNodeFlag = false;
            e_node_i = elems(e).nodes(1).coord;
            e_node_f = elems(e).nodes(2).coord;
            for cnp = 1:size(crossNodePoints,1)
                % If new elem crosses existing elem and one of its nodes,
                % cross point is very close to an existing node.
                if (e_node_i(1) == crossNodePoints(cnp,2) && ...
                    e_node_i(2) == crossNodePoints(cnp,3) && ...
                    e_node_i(3) == crossNodePoints(cnp,4)) ||...
                   (e_node_f(1) == crossNodePoints(cnp,2) && ...
                    e_node_f(2) == crossNodePoints(cnp,3) && ...
                    e_node_f(3) == crossNodePoints(cnp,4))
                
                    redundantNodeFlag = true;
                    break
                end
            end
            if redundantNodeFlag == false
                countNN = countNN + 1;
                newNodes(countNN,:) = crossPoints(cep,2:end);
            end
        end
    end
    if countNN ~= 0
        newNodes = newNodes(1:countNN,:);
    else
        newNodes = [];
    end
end

%--------------------------------------------------------------------------
% Returns nodes to be created and drawn on intersections.
function newNodes = getNewNodes(inputArg)
    crossNodePoints = inputArg{1};
    collinearElems = inputArg{2};
    elemConnect = inputArg{3};
    if size(inputArg,2) == 4
        tol = inputArg{4};
    else
        tol = 1e-10;
    end

    % Get number of elements (not counting elems being created)
    nel  = getappdata(0,'nel');
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Initialize matrix to be returned (new node coordinates)
    newNodes = zeros(nel,3);
    
    % Initialize new nodes counter
    countNewNode = 0;
    
    % Loop through all new elements
    for e = 1:size(elemConnect,1)
        % Get new element end coordinates
        new_elem_coords = [nodes(elemConnect(e,1)).coord;
                           nodes(elemConnect(e,2)).coord];

        % Check if new element cross existing elements
        [~,newNodes_e] = getCrossElemPoints({new_elem_coords,crossNodePoints,collinearElems,tol});
        for i = 1:size(newNodes_e,1)
            if countNewNode == 0
                % Update counter
                countNewNode = countNewNode + 1;
                newNodes(countNewNode,:) = newNodes_e(i,:);
            elseif all(abs(newNodes(1:countNewNode,1)-newNodes_e(i,1)) > tol) || ...
                   all(abs(newNodes(1:countNewNode,2)-newNodes_e(i,2)) > tol) || ...
                   all(abs(newNodes(1:countNewNode,3)-newNodes_e(i,3)) > tol)
                % Update counter
                countNewNode = countNewNode + 1;
                newNodes(countNewNode,:) = newNodes_e(i,:);
            end
        end
    end
    
    newNodes = newNodes(1:countNewNode,:);
end

%--------------------------------------------------------------------------
function elemPointValues = getElemPointDisplAndStress(fctnArgIn)
    % Get function input arguments
    coords = fctnArgIn{1}; % click point on element
    e = fctnArgIn{2};      % element id

    % Get vector of handles to element objects
    elems = getappdata(0,'elems');
    
    % Get click point coordinate along element local X axis
    local_X = norm(coords - elems(e).nodes(1).coord(1:size(coords,2)));
    
    % Parametric coordinate along element local X axis
    t = local_X/elems(e).length;
    
    % Get closest element segment to click point
    whichSegment = (size(elems(e).intCoords,2) - 1) * t;
    
    % Parametric coordinate along 'whichSegment'
    if t == 1
        t_seg = 1;
    else
        t_seg = rem(whichSegment,1);
    end

    % Get closest internal points to click point
    whichPoint = whichSegment - t_seg + 1;
    intPts = [whichPoint, whichPoint + 1];
    
    % Linear interpolation of displacement values
    displ = elems(e).intDispl(:,intPts(1)) * (1 - t_seg) + ...
            elems(e).intDispl(:,intPts(2)) * t_seg;
    
    % Avoid residual values
    for i = 1:size(displ,1)
        if abs(displ(i)) <= 10^-12
            displ(i) = 0;
        end
    end
    
    % Linear interpolation of srtess values
    if size(elems(e).intStresses,2) > 2
        stress = elems(e).intStresses(:,intPts(1)) * (1 - t_seg) + ...
                 elems(e).intStresses(:,intPts(2)) * t_seg;
    elseif size(elems(e).intStresses,2) == 2
        stress = elems(e).intStresses(:,1) * (1 - t) +...
                 elems(e).intStresses(:,2) * t;     
    else
        stress = elems(e).intStresses(:,1);
    end
    
    % Avoid residual values
    for i = 1:size(stress,1)
        if abs(stress(i)) <= 10^-12
            stress(i) = 0;
        end
    end
    
    % Assemble array to be returned
    sz = max([size(stress,1), size(displ,1)]);
    elemPointValues = zeros(sz,3);
    elemPointValues(1,1) = local_X;
    elemPointValues(1:size(displ,1),2) = displ;
    elemPointValues(1:size(stress,1),3) = stress;
end

%--------------------------------------------------------------------------
function divideElement(inputArg)
    e = inputArg{1};  % element to be divided
    n = inputArg{2};  % new node id
    
    % Get handle to GUI_Main
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Get analysis model
    anm = get(mdata.popupmenu_Anm,'Value');
    
    % Get vector of handles to node objects
    nodes = getappdata(0,'nodes');
    
    % Get coordinates of the point where element will be divided
    x = nodes(n).coord(1);
    y = nodes(n).coord(2);
    z = nodes(n).coord(3);
     
    % Get vector of handles to elem objects
    elems = getappdata(0,'elems');
    
    % Get elem object properties
    type = elems(e).type;
    analisys = elems(e).anm;
    mat = elems(e).material;
    sec = elems(e).section;
    ni = elems(e).nodes(1);
    nf = elems(e).nodes(2);
    hi = elems(e).hingei;
    hf = elems(e).hingef;
    vz = elems(e).vz;
    kri = elems(e).kri;
    krf = elems(e).krf;
    load = elems(e).load;

    % Check if model is truss or frame
    if anm == 1 || anm == 4 
        newHinge = 0;
    else
        newHinge = 1;
    end

    % Check if there are element loads on elem to be divided
    load_1_elemLoadCase = load.elemLoadCase;
    load_2_elemLoadCase = load.elemLoadCase;
    linearGbl_1 = [0 0 0 0 0 0];
    linearGbl_2 = [0 0 0 0 0 0];
    if ~isempty(load.elemLoadCase)
        % Check if there is linear load
        % Other Lelem properties remain the same for both elems
        if all(all(load.elemLoadCase(6:11,1:end) == 0) == 1) == 0
            %�Get elem dimensions
            L = elems(e).length;
            xe1 = elems(e).nodes(1).coord(1);
            ye1 = elems(e).nodes(1).coord(2);
            ze1 = elems(e).nodes(1).coord(3);
            Lp = sqrt((x - xe1)^2 + (y - ye1)^2 + (z - ze1)^2);

            % Update load cases
            for nlc = 1:size(load.elemLoadCase,2)
                if all(load.elemLoadCase(6:11,nlc) == 0) == 0
                    % Get linear loads
                    qxi = load.elemLoadCase(6,nlc);
                    qxf = load.elemLoadCase(9,nlc);
                    dqx = qxf - qxi;
                    qyi = load.elemLoadCase(7,nlc);
                    qyf = load.elemLoadCase(10,nlc);
                    dqy = qyf - qyi;
                    qzf = load.elemLoadCase(8,nlc);
                    qzi = load.elemLoadCase(11,nlc);
                    dqz = qzf - qzi;

                    % Set new lin loads to load cases (elem_1)
                    load_1_elemLoadCase(6:8,nlc) = [qxi qyi qzi];
                    load_1_elemLoadCase(9:11,nlc) = [(qxi+dqx*Lp/L) (qyi+dqy*Lp/L) (qzi+dqz*Lp/L)];

                    % Set new lin loads to load cases (elem_2)
                    load_2_elemLoadCase(6:8,nlc) = [(qxi+dqx*Lp/L) (qyi+dqy*Lp/L) (qzi+dqz*Lp/L)];
                    load_2_elemLoadCase(9:11,nlc) = [qxf qyf qzf];
                end
            end

            % Update current load case
            if ~isempty(load.linearGbl)
                % Get current load case
                lc = get(mdata.popupmenu_LoadCase,'value');

                % Get linear loads
                qxi = load.linearGbl(1);
                qxf = load.linearGbl(4);
                dqx = qxf - qxi;
                qyi = load.linearGbl(2);
                qyf = load.linearGbl(5);
                dqy = qyf - qyi;
                qzi = load.linearGbl(3);
                qzf = load.linearGbl(6);
                dqz = qzf - qzi;

                % Set new linear loads (elem_1)
                linearGbl_1(1:3) = [qxi qyi qzi];
                linearGbl_1(4:6) = [(qxi+dqx*Lp/L) (qyi+dqy*Lp/L) (qzi+dqz*Lp/L)];
                load_1_elemLoadCase(5:11,lc) = [0; linearGbl_1'];

                % Set new linear loads (elem_2)
                linearGbl_2(1:3) = [(qxi+dqx*Lp/L) (qyi+dqy*Lp/L) (qzi+dqz*Lp/L)];
                linearGbl_2(4:6) = [qxf qyf qzf];
                load_2_elemLoadCase(5:11,lc) = [0; linearGbl_2'];
            end
        end
    end

    % Create two new elements
    % Elem_1
    dividedElem(1) = Elem(type,analisys,mat,sec,[ni nodes(n)],hi,newHinge,vz,kri,[]);
    dividedElem(1).load.elemLoadCase = load_1_elemLoadCase;
    dividedElem(1).load.uniformDir = load.uniformDir;
    dividedElem(1).load.uniformGbl = load.uniformGbl;
    dividedElem(1).load.uniformLcl = load.uniformLcl;
    dividedElem(1).load.linearDir = 0;
    dividedElem(1).load.linearGbl = [];
    dividedElem(1).load.linearLcl = [];
    dividedElem(1).load.setLinearLoad(linearGbl_1,0)
    dividedElem(1).load.tempVar_X = load.tempVar_X;
    dividedElem(1).load.tempVar_Y = load.tempVar_Y;
    dividedElem(1).load.tempVar_Z = load.tempVar_Z;

    % Elem_2
    nel = getappdata(0,'nel') + 1;
    dividedElem(2) = Elem(type,analisys,mat,sec,[nodes(n) nf],newHinge,hf,vz,[],krf);
    dividedElem(2).load.elemLoadCase = load_2_elemLoadCase;
    dividedElem(2).load.uniformDir = load.uniformDir;
    dividedElem(2).load.uniformGbl = load.uniformGbl;
    dividedElem(2).load.uniformLcl = load.uniformLcl;
    dividedElem(2).load.linearDir = 0;
    dividedElem(2).load.linearGbl = [];
    dividedElem(2).load.linearLcl = [];
    dividedElem(2).load.setLinearLoad(linearGbl_2,0)
    dividedElem(2).load.tempVar_X = load.tempVar_X;
    dividedElem(2).load.tempVar_Y = load.tempVar_Y;
    dividedElem(2).load.tempVar_Z = load.tempVar_Z;
    
    % Check if new elements are equal to an existing element
    elemEqualFlag_1 = isElemEqual({dividedElem(1).nodes(1),dividedElem(1).nodes(2)});
    elemEqualFlag_2 = isElemEqual({dividedElem(2).nodes(1),dividedElem(2).nodes(2)});
    
    % Delete original element
    elems(e) = [];
    
    % Update elems vector
    if elemEqualFlag_1 == false && elemEqualFlag_2 == false
        elems = horzcat(elems(1:e-1),dividedElem(1),elems(e:end),dividedElem(2));
    elseif elemEqualFlag_1 ~= false
        elems = horzcat(elems,dividedElem(2));
        nel = nel - 1;
    elseif elemEqualFlag_2 ~= false
        elems = horzcat(elems(1:e-1),dividedElem(1),elems(e:end));
        nel = nel - 1;
    end

    % Update draw and model object
    model = getappdata(0,'model');
    model.elems = elems;
    model.nel = nel;
    draw = getappdata(0,'draw');
    draw.mdl = model;
    
    % Set info panel data
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(4,:) = {'Elements',nel};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    
    if (strcmp(get(mdata.elemIDButton,'Checked'),'on') == 1) % Check if elements ID is on
        delete(findobj('tag','textElemID'))
        draw.elementID();
    end
    if (strcmp(get(mdata.orientationButton,'Checked'),'on') == 1) % Check if elements orientation is on
        delete(findobj('tag','drawElemOrient'))
        draw.elementOrientation();
    end
    
    % Return variables to root
    setappdata(0,'nel',nel)
    setappdata(0,'elems',elems)
    setappdata(0,'model',model)
    setappdata(0,'draw',draw)
    
    % Get graphic tolerance
    dfltUnits = get(mdata.axes_Canvas,'units');
    set(mdata.axes_Canvas,'units','normalized');
    limits = get(mdata.axes_Canvas,'Position');
    set(mdata.axes_Canvas,'units',dfltUnits);
    tol = limits(3)/22;
    
    % Check if divided element had more intersections, and allocate them in
    % the corresponding new element
    intersections = getappdata(0,'intersections');
    if ~isempty(intersections)
        nis = size(intersections,2);
        delIntersections = zeros(1,nis);
        for ni = 1:nis
            if intersections(ni).coord(1) == x && intersections(ni).coord(2) == y &&...
               intersections(ni).coord(3) == z
                delIntersections(ni) = ni;
            elseif ~all(intersections(ni).elems ~= e)
                whichElems = isPointInElem({intersections(ni).coord,tol});
                intersections(ni).elems = whichElems;
            end
        end
        intersections(nonzeros(delIntersections)') = [];
        % Return variable to root
        setappdata(0,'intersections',intersections)
    end
end

%--------------------------------------------------------------------------
% Auxiliary function to GUI_Main_KeyPressFcn.
% Delete a Node object from the list of nodes when 'delete' key is pressed.
function deleteNodes(del_n)
    model = getappdata(0,'model');
    draw = getappdata(0,'draw');
    nodes = getappdata(0,'nodes');
    nnp = getappdata(0,'nnp');
    elems = getappdata(0,'elems');
    nel = getappdata(0,'nel');
    mdata = guidata(findobj('Tag','GUI_Main'));

    % Get number of elements connected to the deleted node
    n_elem = 0;
    for e = 1:nel
        if (elems(e).nodes(1).id == nodes(del_n).id) || (elems(e).nodes(2).id == nodes(del_n).id)
            n_elem = n_elem + 1;
            del_elem(n_elem) = e; %#ok<AGROW>
        end
    end

    if n_elem ~= 0
        choice = questdlg('Deleting this node will automatically delete all elements connected to it.','Delete Node','OK','Cancel','Cancel');
    else
        choice = 'OK';
    end

    if strcmp(choice,'OK')
        % Initialize flags
        elemsNeedToBeRedrawn      = false;
        nodalLoadsNeedToBeRedrawn = false;
        elemLoadsNeedToBeRedrawn  = false;
        
        % Delete all elements connected to the deleted node
        for e = 1:n_elem
            % Get ID of deleted element
            elem_id = del_elem(e) - (e-1);
            
            % Check if deleted element had distributed or thermal loads
            if ~isempty(elems(elem_id).load.uniformGbl) || ~isempty(elems(elem_id).load.uniformLcl) ||...
               ~isempty(elems(elem_id).load.linearGbl) || ~isempty(elems(elem_id).load.linearLcl) ||...
               elems(elem_id).load.tempVar_X ~= 0 || elems(elem_id).load.tempVar_Y ~= 0 ||...
               elems(elem_id).load.tempVar_Z ~= 0
                elemLoadsNeedToBeRedrawn = true;
            end

            % Remove deleted element from vector of elements
            if elem_id == 1
                elems = elems(2:end);
            elseif elem_id == nel
                elems = elems(1:end-1);
            else
                ea = elems(1:elem_id-1);
                eb = elems(elem_id+1:end);
                elems = cat(2,ea,eb);
            end

            % Update number of elements
            nel = nel - 1;

            % Set model object properties
            model.elems = elems;
            model.nel = nel;

            % Check if element had any unresolved intersections
            intersections = getappdata(0,'intersections');
            if ~isempty(intersections)
                delIntSec = zeros(1,size(intersections,2));
                for nis = 1:size(intersections,2)
                    if ~all(intersections(nis).elems ~= del_elem(e))
                        intersections(nis).elems = nonzeros(intersections(nis).elems - del_elem(e))' + del_elem(e);
                        if size(intersections(nis).elems,2) <= 1
                            delIntSec(nis) = nis;
                        end
                    end
                end
                if ~all(delIntSec == 0)
                    intersections(nonzeros(delIntSec)') = [];
                end
                setappdata(0,'intersections',intersections)
            end
        end
        
        % Check if elements were deleted
        if n_elem ~= 0
            elemsNeedToBeRedrawn = true;
        end
        
        % Check if deleted node had loads or precribed displacements
        if ~isempty(nodes(del_n).load.static)  ||...
           ~isempty(nodes(del_n).load.dynamic) ||...
           ~isempty(nodes(del_n).prescDispl)   ||...
           ~isempty(nodes(del_n).initCond)     ||...
           ~isempty(nodes(del_n).displMass)    ||...
           ~isempty(nodes(del_n).rotMass)
            if ~all(nodes(del_n).load.static  == 0) ||...
               ~all(nodes(del_n).load.dynamic == 0) ||...
               ~all(nodes(del_n).prescDispl   == 0) ||...
               any(any(nodes(del_n).initCond))      ||...
               ~all(nodes(del_n).displMass    == 0) ||...
               ~all(nodes(del_n).rotMass      == 0)
                nodalLoadsNeedToBeRedrawn = true;
            end
        end

        % Get number of fixed and spring dofs related to the deleted node
        countFixed = 0;
        countSpring = 0;
        for i = 1:6
            if nodes(del_n).ebc(i) == 1
                countFixed = countFixed + 1;
            elseif nodes(del_n).ebc(i) == 2
                countSpring = countSpring + 1;
            end    
        end

        % Remove deleted node from vector of nodes
        if del_n == 1
            nodes = nodes(2:end);
        elseif del_n == nnp
            nodes = nodes(1:end-1);
        else
            na = nodes(1:del_n-1);
            nb = nodes(del_n+1:end);
            nodes = cat(2,na,nb);
        end

        % Update number of nodes
        nnp = nnp - 1;

        % Update nodes id
        for n = 1:nnp
            nodes(n).id = n;
        end

        % Update information panel in GUI_Main
        mouse = getappdata(0,'mouse');
        set(mdata.uitable_infoPanel,'Data',mouse.originalData)
        infoPanelData = get(mdata.uitable_infoPanel,'Data');
        infoPanelData(3,:) = {'Nodes',nnp};
        infoPanelData(4,:) = {'Elements',nel};
        anm = get(mdata.popupmenu_Anm,'Value');
        ndof = cell2mat(infoPanelData(5,2));
        nfreedof = cell2mat(infoPanelData(6,2));
        nfixeddof = cell2mat(infoPanelData(7,2));
        nspringdof = cell2mat(infoPanelData(8,2));
        if anm == 1
            ndof = ndof - 2;
            nfreedof = nfreedof - 2 + countFixed + countSpring;
        elseif anm == 2 || anm == 3 || anm == 4
            ndof = ndof - 3;
            nfreedof = nfreedof - 3 + countFixed + countSpring;
        else
            ndof = ndof - 6;
            nfreedof = nfreedof - 6 + countFixed + countSpring;
        end
        nfixeddof = nfixeddof - countFixed;
        nspringdof = nspringdof - countSpring;
        infoPanelData(5,:) = {'DOFs',ndof};
        infoPanelData(6,:) = {'Free DOFs',nfreedof};
        infoPanelData(7,:) = {'Fixed DOFs',nfixeddof};
        infoPanelData(8,:) = {'Springs',nspringdof};
        set(mdata.uitable_infoPanel,'Data',infoPanelData)
        set(mdata.uitable_infoPanelEditable,'Enable','off','Data',{})

        % Set model object properties
        model.nodes = nodes;
        model.nnp = nnp;

        % Enable/disable solve intersections pushbutton (toolbar)
        if size(getappdata(0,'intersections'),2) >= 1
            set(mdata.pushbutton_SolveIntSects,'enable','on')
        else
            set(mdata.pushbutton_SolveIntSects,'enable','off')
        end

        % Enable "Process Data" button in main GUI
        set(mdata.pushbutton_ProcessData,'Enable','on');

        % Enable model type option
        if nnp == 0
            set(mdata.popupmenu_Anm,'Enable','on');
            setappdata(0,'vis',0);
        end

        % Disable result buttons
        if get(mdata.popupmenu_Results,'value') ~= 1
            nodalLoadsNeedToBeRedrawn = true;
            elemLoadsNeedToBeRedrawn = true;
        end
        set(mdata.popupmenu_Results,'Enable','off','value',1)
        set(mdata.pushbutton_Textual,'Enable','off');
        set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
        set(mdata.text_Element,'string','Elements');
        set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
        set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
        set(mdata.pushbutton_DynamicResults,'enable','off');
        set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
        
        % Return variables to root
        setappdata(0,'resultType',0);
        setappdata(0,'elems',elems);
        setappdata(0,'nodes',nodes);
        setappdata(0,'nel',nel);
        setappdata(0,'nnp',nnp);
        setappdata(0,'model',model);
        draw.mdl = model;
        setappdata(0,'draw',draw);

        % Draw updated model
        if anm == 1 || anm == 2 || anm == 3 
            redraw(mdata,'Nodes',false)
        else
            redraw(mdata,'Nodes')
        end
        if elemsNeedToBeRedrawn == true
            redraw(mdata,'Elements')
        end
        if nodalLoadsNeedToBeRedrawn == true && elemLoadsNeedToBeRedrawn == false
            redraw(mdata,'Nodal Loads')
        elseif nodalLoadsNeedToBeRedrawn == false && elemLoadsNeedToBeRedrawn == true
            redraw(mdata,'Element Loads')
        elseif nodalLoadsNeedToBeRedrawn == true && elemLoadsNeedToBeRedrawn == true
            redraw(mdata,'Loads')
        end
        
        delete(findobj('tag', 'snapNode'));
        delete(findobj('tag', 'snapNode2'));
        delete(findobj('tag', 'snapElem'));
        delete(findobj('tag', 'snapElem2'));
        delete(findobj('tag', 'selectedNode'));
        
        mouse.originalData = {};
        mouse.selectedNode = 0;
        mouse.nodeSnap = [];
        mouse.whichNodeSnap = 0;
        if getappdata(0,'nnp') == 0
            mouse.sizeFlag = 0;
        end

        setappdata(0,'mouse',mouse);
    end
end

%--------------------------------------------------------------------------
% Auxiliary function to GUI_Main_KeyPressFcn.
% Delete an Element object from the list of elements when 'delete' key is 
% pressed.
function deleteElems(del_elem)
    model = getappdata(0,'model');
    elems = getappdata(0,'elems');
    nel = getappdata(0,'nel');
    mdata = guidata(findobj('Tag','GUI_Main'));
    
    % Check if deleted element had distributed or thermal loads
    if ~isempty(elems(del_elem).load.uniformGbl) || ~isempty(elems(del_elem).load.uniformLcl) ||...
       ~isempty(elems(del_elem).load.linearGbl) || ~isempty(elems(del_elem).load.linearLcl) ||...
       elems(del_elem).load.tempVar_X ~= 0 || elems(del_elem).load.tempVar_Y ~= 0 ||...
       elems(del_elem).load.tempVar_Z ~= 0
        loadsNeedToBeRedrawn = true;
    else
        loadsNeedToBeRedrawn = false;
    end

    % Remove deleted element from vector of elements
    if del_elem == 1
        elems = elems(2:end);
    elseif del_elem == nel
        elems = elems(1:end-1);
    else
        ea = elems(1:del_elem-1);
        eb = elems(del_elem+1:end);
        elems = cat(2,ea,eb);
    end

    % Update number of elements
    nel = nel - 1;

    % Set model object properties
    model.elems = elems;
    model.nel = nel;

    % Check if element had any unresolved intersections
    intersections = getappdata(0,'intersections');
    if ~isempty(intersections)
        delIntSec = zeros(1,size(intersections,2));
        for nis = 1:size(intersections,2)
            if ~all(intersections(nis).elems ~= del_elem)
                intersections(nis).elems = nonzeros(intersections(nis).elems - del_elem)' + del_elem;
                if size(intersections(nis).elems,2) <= 1
                    delIntSec(nis) = nis;
                end
            end
        end
        if ~all(delIntSec == 0)
            intersections(nonzeros(delIntSec)') = [];
        end
        setappdata(0,'intersections',intersections)
    end

    % Enable/disable solve intersections pushbutton (toolbar)
    if size(getappdata(0,'intersections'),2) >= 1
        set(mdata.pushbutton_SolveIntSects,'enable','on')
    else
        set(mdata.pushbutton_SolveIntSects,'enable','off')
    end

    % Enable "Process Data" button in main GUI
    set(mdata.pushbutton_ProcessData,'Enable','on');

    % Disable result buttons
    allLoadsNeedToBeRedrawn = false;
    if get(mdata.popupmenu_Results,'value') ~= 1
        allLoadsNeedToBeRedrawn = true;
    end
    set(mdata.popupmenu_Results,'Enable','off','value',1);
    set(mdata.pushbutton_Textual,'Enable','off');
    set(mdata.checkbox_Reactions,'Enable','off', 'Value', 0);
    set(mdata.text_Element,'string','Elements');
    set(mdata.popupmenu_ElementResults,'Visible','off','Enable','off','String',' ','Value',1,'Max',1);
    set(mdata.edit_ElementResults,'Visible','on','Enable','off','String','All');
    set(mdata.pushbutton_DynamicResults,'enable','off');
    set(mdata.pushbutton_PlayDynamicResults,'enable','off','UserData',false,'ForegroundColor',[0,0.7,0.2]);
    
    % Update information panel in GUI_Main
    mouse = getappdata(0,'mouse');
    set(mdata.uitable_infoPanel,'Data',mouse.originalData)
    infoPanelData = get(mdata.uitable_infoPanel,'Data');
    infoPanelData(4,:) = {'Elements',nel};
    set(mdata.uitable_infoPanel,'Data',infoPanelData)
    set(mdata.uitable_infoPanelEditable,'Enable','off','Data',{})

    % Return variables to root
    setappdata(0,'resultType',0);
    setappdata(0,'elems',elems);
    setappdata(0,'nel',nel);
    setappdata(0,'model',model);

    % Draw updated model
    redraw(mdata,'Elements')
    if allLoadsNeedToBeRedrawn == true
        redraw(mdata,'Loads')
    elseif loadsNeedToBeRedrawn == true
        redraw(mdata,'Element Loads')
    end
    delete(findobj('tag', 'snapNode'));
    delete(findobj('tag', 'snapNode2'));
    delete(findobj('tag', 'snapElem'));
    delete(findobj('tag', 'snapElem2'));

    % Reinitialize object for mouse events and save it in root (make them
    % acessible to all GUIs)
    mouse.originalData = {};
    mouse.selectedElem = 0;
    mouse.elemSnap = [];
    mouse.whichElemSnap = 0;
    setappdata(0,'mouse',mouse);
end

%% ------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO isPointInElem AND getCrossIntSectPoints.
% NOT CALLED BY auxModelFctn.
%
% Work as a switch between 2D and 3D functions
% Check if given point is inside a specific element.
% Returns flag. 0 = no; 1 = yes.
function pointIn = isPointInThisElem(pCoords,eCoords,tol,planeOrSpatial,clickFlag)
    if nargin < 4
        sz = size(eCoords,2);
        switch sz
            case 2
                planeOrSpatial = '2D';
            case 3
                planeOrSpatial = '3D';
                clickFlag = false;
        end
    elseif nargin < 5
        clickFlag = false;
    end

    switch planeOrSpatial
        case '2D'
            pointIn = isPointInThisElem2D(pCoords,eCoords,tol);
        case '3D'
            pointIn = isPointInThisElem3D(pCoords,eCoords,tol,clickFlag);
    end
end

% -------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO isPointInElem AND getCrossIntSectPoints.
% NOT CALLED BY auxModelFctn.
%
% Check if given point is inside a specific element.
% Returns flag. 0 = no; 1 = yes.
function pointIn = isPointInThisElem2D(pointCoords,elemCoords,tol)
    % Get numeric tolerance
    if nargin < 3
        tol = 1e-15;
    end

    % Initialize variable to be returned
    pointIn = false;
    
    % Get point coordinates
    point_x = pointCoords(1);
    point_y = pointCoords(2);
    
    % Get elem end coordinates
    e_nodei_x = elemCoords(1,1);
    e_nodei_y = elemCoords(1,2);
    e_nodef_x = elemCoords(2,1);
    e_nodef_y = elemCoords(2,2);
    
    % Define triangle points
    a = [point_x, point_y];
    b = [e_nodei_x, e_nodei_y];
    c = [e_nodef_x, e_nodef_y];
    
    % Get elem length
    L = norm(c-b);
    
    % Get signed area
    orientABC = orient2D(a,b,c);
    
    % Define min and max x and y coords
    minX = min([b(1) c(1)]);
    maxX = max([b(1) c(1)]);
    minY = min([b(2) c(2)]);
    maxY = max([b(2) c(2)]);
    
    distMinX = a(1) - minX;
    distMaxX = a(1) - maxX;
    distMinY = a(2) - minY;
    distMaxY = a(2) - maxY;
    
    % Check if point is inside elem (numeric tolerance)
    if (abs(orientABC)*2/L) < tol && distMinX > - tol &&...
        distMaxX < tol && distMinY > - tol && distMaxY < tol

        pointIn = true;
    end
    
end

%--------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO isPointInElem.
% NOT CALLED BY auxModelFctn.
%
% Check if given point is inside a specific element.
% Returns flag. 0 = no; 1 = yes.
function pointIn = isPointInThisElem3D(coords,elemCoords,tol,clickFlag)
    % Get numeric tolerance
    if nargin < 3
        tol = 1e-15;
        clickFlag = false;
    elseif nargin < 4
        clickFlag = false;
    end
    tol_diff = 1e-10;
    
    % Initialize variable to be returned
    pointIn = false;
    
    % Get point coordinates
    x = coords(1);
    y = coords(2);
    z = coords(3);

    % Get element end coordinates
    xe1 = elemCoords(1,1);
    ye1 = elemCoords(1,2);
    ze1 = elemCoords(1,3);
    xe2 = elemCoords(2,1);
    ye2 = elemCoords(2,2);
    ze2 = elemCoords(2,3);
    
    % Check if point is equal to one of the element's ends.
    if (abs(x-xe1) < tol_diff && abs(y-ye1) < tol_diff && abs(z-ze1) < tol_diff) ||...
       (abs(x-xe2) < tol_diff && abs(y-ye2) < tol_diff && abs(z-ze2) < tol_diff)
        pointIn = true;
        return
    end
    
    if clickFlag == true
        if abs(xe1-xe2) < tol_diff && abs(ye1-ye2) < tol_diff
            if norm([xe1 - x, ye1 - y]) <= tol
                pointIn = true;
                return
            end
        elseif abs(ze1-ze2) < tol_diff && abs(ye1-ye2) < tol_diff
            if norm([ye1 - y, ze1 - z]) <= tol
                pointIn = true;
                return
            end
        elseif abs(xe1-xe2) < tol_diff && abs(ze1-ze2) < tol_diff
            if norm([xe1 - x, ze1 - z]) <= tol
                pointIn = true;
                return
            end
        end
    end
    
    pLine = perpPointLine3D([xe1 ye1 ze1; xe2 ye2 ze2], [x y z]);
    
    if isempty(pLine) 
        return
    end
    
    if norm(pLine(1,:) - pLine(2,:)) <=  tol
        pointIn = true;
    end
    return
        

    % Get element angle with X axis (plane XY)
    alpha = atan((ye2-ye1)/(xe2-xe1));
    if isnan(alpha)
        if (ye2 - ye1) > 0
            alpha = pi * 0.5;
        else
            alpha = - pi * 0.5;
        end
    elseif (xe2-xe1) < 0
        alpha = alpha + pi;
    end

    % Get node_i-to-point angle with X axis (plane XY)
    alphaP = atan((y-ye1)/(x-xe1));
    if isnan(alphaP)
        if (y - ye1) > 0
            alphaP = pi * 0.5;
        else
            alphaP = - pi * 0.5;
        end
    elseif (x-xe1) < 0
        alphaP = alphaP + pi;
    end

    % Get element angle with plane XY (3D)
    beta = atan((ze2-ze1)/norm([(xe2-xe1) (ye2-ye1)]));
    if isnan(beta)
        if (ze2-ze1) >= 0
            beta = pi * 0.5;
        else
            beta = - pi * 0.5;
        end
    end

    % Get node_i-to-point angle with plane XY (3D)
    betaP = atan((z-ze1)/norm([(x-xe1) (y-ye1)]));
     if isnan(betaP)
        if (z-ze1) >= 0
            betaP = pi * 0.5;
        else
            betaP = - pi * 0.5;
        end
    end

    % Get angle between element and node_i-to-point (plane XY)
    theta = alphaP - alpha;

    % Get angle between element and node_i-to-point (3D)
    phi = betaP - beta;

    % Get element length
    L = norm([(xe2-xe1) (ye2-ye1) (ze2-ze1)]);

    % Get node_i-to-point length
    Lp = norm([(x-xe1) (y-ye1) (z-ze1)]);

    % Check if point is inside element
    if abs(theta) <= tol && abs(phi) <= tol && Lp <= L
        pointIn = true;
    end
end

%--------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO areLinesCrossed3D.
% NOT CALLED BY auxModelFctn.
%
% Compute shortest orthogonal line between two linear directions in 3D
function line_3 = perpLine3D(line_1,line_2)
    % Get line directions
    dir_1 = (line_1(2,:) - line_1(1,:));
    dir_2 = (line_2(2,:) - line_2(1,:));
    
    % Get line lengths
    length_1 = norm(line_1(2,:)-line_1(1,:));
    length_2 = norm(line_2(2,:)-line_2(1,:));
    
    % Get direction othogonal to both lines
    dir_3 = cross(dir_1,dir_2);
    
    % Assemble auxiliar 3x3 matrix
    A = [dir_1', dir_2', dir_3'];
    
    % Get scale factors
    dx = line_2(1,1) - line_1(1,1);
    dy = line_2(1,2) - line_1(1,2);
    dz = line_2(1,3) - line_1(1,3);
    
    % Check if A is ill conditioned (L1 and L2 are parallel)
    if abs(rcond(A)) > 10^-10
        pCoord = A \ [dx dy dz]';
        pCoord(2) = -pCoord(2);
    else
        line_3 = [];
        return
    end
    
    % Compute and assemble ortho line coordinates
    x1 = line_1(1,1) + dir_1(1)*pCoord(1);
    y1 = line_1(1,2) + dir_1(2)*pCoord(1);
    z1 = line_1(1,3) + dir_1(3)*pCoord(1);
    
    x2 = line_2(1,1) + dir_2(1)*pCoord(2);
    y2 = line_2(1,2) + dir_2(2)*pCoord(2);
    z2 = line_2(1,3) + dir_2(3)*pCoord(2);
    
    if norm([x1 y1 z1] - line_1(1,:)) <= length_1 &&...
       norm([x1 y1 z1] - line_1(2,:)) <= length_1 &&...
       norm([x2 y2 z2] - line_2(1,:)) <= length_2 &&...
       norm([x2 y2 z2] - line_2(2,:)) <= length_2
    
        line_3 = [x1 y1 z1;
                  x2 y2 z2];
    else
        line_3 = [];
    end
end

%--------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO areLinesCrossed3D.
% NOT CALLED BY auxModelFctn.
%
% Compute shortest orthogonal line between a line and a ponit in 3D
function pLine = perpPointLine3D(line,point)

    % Get vertexes - line: AB, point: C
    A = line(1,:);
    B = line(2,:);
    C = zeros(1,3); C(:) = point(:);
    
    % Compute directions
    AB = B-A;
    AC = C-A;
    
    % Compute parametric coords in AB
    tau = dot(AB,AC) / dot(AB,AB);
    
    % Check if projection is out of AB bounds
    if tau > 1.0 || tau < 0.0
        pLine = [];
        return;
    end
    
    % Compute perp line
    pLine = [ A + tau*AB; C ];
end

%--------------------------------------------------------------------------
% THIS IS AN AUXILIARY FUNCTION TO areElemsCrossed AND isPointInThisElem2D.
% NOT CALLED BY auxModelFctn.
%
% Calculate area of a triangle, respecting orientation (signed area).
function Area = orient2D(A,B,C)
    Area = 0.5 * det([A(1) A(2) 1;
                      B(1) B(2) 1;
                      C(1) C(2) 1]);
end
