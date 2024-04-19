%% Dynamic Results Dialog Callback Functions
% This file contains the callback functions associated with the "Dynamic
% Results" dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
function varargout = GUI_DynamicResults(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_DynamicResults_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_DynamicResults_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%--------------------------------------------------------------------------
% --- Executes just before GUI_DynamicResults is made visible.
function GUI_DynamicResults_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
include_constants;

% Move GUI to the center of the screen (ensures size 950px X 530px)
w = 950; h = 530;
set(handles.GUI_DynamicResults,'Units','pixels');
set(handles.GUI_DynamicResults,'Position',[-2*w, -2*h, w, h]);
set(handles.GUI_DynamicResults,'Units','normalized');
pos = get(handles.GUI_DynamicResults,'Position');
set(handles.GUI_DynamicResults,'Position',[(1-pos(3))/2, (1-pos(4))/2, pos(3), pos(4)]);

% Get model object from root
model = getappdata(0,'model');

% Set nodes popupmenu properties
string_nodes = num2str(1:model.nnp,'%d\n');
set(handles.popupmenu_Nodes,'string',string_nodes,'value',1,'Max',model.nnp,'enable','on')

% Determine string and maximum value of the dof popupmenu
switch model.anm.analysis_type
    case TRUSS2D_ANALYSIS
        string_dof = {'Dx';'Dy'};
        max_dof = 2;
    case FRAME2D_ANALYSIS
        string_dof = {'Dx';'Dy';'Rz'};
        max_dof = 3;
    case GRILLAGE_ANALYSIS
        string_dof = {'Rx';'Ry';'Dz'};
        max_dof = 3;
    case TRUSS3D_ANALYSIS
        string_dof = {'Dx';'Dy';'Dz'};
        max_dof = 3;
    case FRAME3D_ANALYSIS
        string_dof = {'Dx';'Dy';'Dz';'Rx';'Ry';'Rz'};
        max_dof = 6;
end

% Set dof popupmenu properties
set(handles.popupmenu_DOF,'string',string_dof,'value',1,'Max',max_dof,'enable','on')

% Set result popupmenu properties
set(handles.popupmenu_Result,'value',1,'enable','on')

% Set time limits on editable text boxes
set(handles.edit_Ti,'string','0');
set(handles.edit_Tf,'string',num2str(model.t));

% Check for results type
if model.results.type == DYNAMIC_NEWMARK_LINEAR || model.results.type == DYNAMIC_RK4_LINEAR || model.results.type == DYNAMIC_AM3_LINEAR || model.results.type == DYNAMIC_WILSON_LINEAR
    % Write tableData 
    tableData = {'All', true, 'Default', 'Default', 'Default'};
    
    % Set uitable_PlotModes properties
    set(handles.uitable_PlotModes,'Enable','off','Data',tableData,...
            'CellEditCallback',@uitable_PlotModes_CellEditCallback,...
            'ColumnEditable',[false false false false false],'ColumnFormat',{'char' 'logical' 'char' 'char' 'char'});
    
elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes == 1
    % Write tableData 
    tableData(1,:) = {'All', true, 'Default', 'Default', 'Default'};
    tableData(2,:) = {sprintf('%i (%.2e Hz)',1,model.W / (2 * pi)),   true, 'Default', 'Default', 'Default'};
    
    % Set uitable_PlotModes properties
    set(handles.uitable_PlotModes,'Enable','off','Data',tableData,...
            'CellEditCallback',@uitable_PlotModes_CellEditCallback,...
            'ColumnEditable',[false false false false false],'ColumnFormat',{'char' 'logical' 'char' 'char' 'char'});
    
elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1
    % Initialize tableData as a cell
    tableData = cell(model.n_modes + 1, 5);
    
    % Set first row of tableData 
    tableData(1,:) = {'All', true, 'Red', 1.00, '-'};
    
    % Auxiliar cell
    aux = {'Red' 'Green' 'Blue' 'Black' 'Magenta' 'Cyan'};
    counter = 0;
    
    % Get natural frequencies to be displayed on the uitable
    f = model.W / (2 * pi);
    
    % Set following rows of tableData 
    for i = 2:model.n_modes+1
        tableData(i,:) = {sprintf('%i (%.2e Hz)',i-1,f(i-1)), false, aux{i-6*counter}, 1.00, '-'};
        if rem(i,6) == 0
            counter = counter + 1;
        end
    end
    
    % Update auxilar cell
    aux{end+1} = 'New Color';
    
    % Set uitable_PlotModes properties
    set(handles.uitable_PlotModes,'Enable','on','Data',tableData,...
            'CellEditCallback',@uitable_PlotModes_CellEditCallback,...
            'ColumnEditable',[false true true true true],'ColumnFormat',{'char' 'logical' aux [] {'-' '--' ':' '-.'}});
end

% Initialize canvases
initCanvases(handles);

% Update result diagrams
popupmenu_DOF_Callback(handles.popupmenu_DOF, [], handles);

% Set callback for pressing down mouse buttons on this figure
set(handles.GUI_DynamicResults, 'WindowButtonDownFcn',  @buttonDownFcn);
set(handles.GUI_DynamicResults, 'WindowButtonMotionFcn',@moveCursorFcn);
set(handles.GUI_DynamicResults, 'WindowButtonUpFcn',    @buttonUpFcn);

% Centralize dialog on screen
pos = get(handles.GUI_DynamicResults,'Position');
set(handles.GUI_DynamicResults,'Position',[(1-pos(3))/2, (1-pos(4))/2, pos(3), pos(4)]);

% Choose default command line output for GUI_DynamicResults
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_PlotModes
function uitable_PlotModes_CellEditCallback(~, eventdata, ~)
% Get handle to this GUI
gui = guidata(findobj('tag','GUI_DynamicResults'));

% Get table data
tableData = get(gui.uitable_PlotModes,'Data');

% Get id of edited cell
if isempty(eventdata.Indices)
    return
elseif size(eventdata.Indices,1) > 1
    id = [eventdata.Indices(end,1),1];
else
    id = eventdata.Indices;
end

% TEMPORARY: CHECK IF SOLUTION IS ON FREQ DOMAIN
if get(gui.popupmenu_Domain,'Value') == 2
    if id(1) > 1
        return
    elseif id(2) == 2
        set(gui.popupmenu_Domain,'Value',1);
        popupmenu_Domain_Callback(gui.popupmenu_Domain,[],gui);
    end
end

% Verify if edited cell was a "display" logical checkbox
if id(2) == 2
    % Get edited cell value
    check = tableData{id(1),2};

    % Plot or delete diagram, according to check
    switch check
        case true  % Draw this diagram
            % Get model object from root
            model = getappdata(0,'model');
            
            % Get node and dof IDs
            n = get(gui.popupmenu_Nodes,'value');
            dof = model.ID(get(gui.popupmenu_DOF,'value'),n);
            
            % Get displacement, Veloc and acceleration vectors
            if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
                d = zeros(1,2);
                v = zeros(1,2);
                a = zeros(1,2);
            elseif id(1) == 1
                % Get results
                if get(gui.popupmenu_Result,'value') == 1 % Total results
                    dd = model.results.dynamicDispl(dof,:,:);
                    vv = model.results.dynamicVeloc(dof,:,:);
                    aa = model.results.dynamicAccel(dof,:,:);
                elseif get(gui.popupmenu_Result,'value') == 2 % Free vibration
                    dd = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
                    vv = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
                    aa = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);
                else % get(gui.popupmenu_Result,'value') == 3 % Forced vibration
                    dd = model.results.dynamicDisplForced(dof,:,:);
                    vv = model.results.dynamicVelocForced(dof,:,:);
                    aa = model.results.dynamicAccelForced(dof,:,:);
                end
                
                % Auxiliar matrices
                aux_d(:,:) = dd;
                aux_v(:,:) = vv;
                aux_a(:,:) = aa;
                
                % Set displacement, Veloc and acceleration vectors as sum
                % of all modes contributions
                d = sum(aux_d,2)';
                v = sum(aux_v,2)';
                a = sum(aux_a,2)';
                
                % Clear auxiliar variables
                clear dd vv aa aux_d aux_v aux_a
            else
                if get(gui.popupmenu_Result,'value') == 1 % Total results
                    d = model.results.dynamicDispl(dof,:,id(1)-1);
                    v = model.results.dynamicVeloc(dof,:,id(1)-1);
                    a = model.results.dynamicAccel(dof,:,id(1)-1);
                elseif get(gui.popupmenu_Result,'value') == 2 % Free vibration
                    d = model.results.dynamicDispl(dof,:,id(1)-1) - model.results.dynamicDisplForced(dof,:,id(1)-1);
                    v = model.results.dynamicVeloc(dof,:,id(1)-1) - model.results.dynamicVelocForced(dof,:,id(1)-1);
                    a = model.results.dynamicAccel(dof,:,id(1)-1) - model.results.dynamicAccelForced(dof,:,id(1)-1);
                else % get(gui.popupmenu_Result,'value') == 3 % Forced vibration
                    d = model.results.dynamicDisplForced(dof,:,id(1)-1);
                    v = model.results.dynamicVelocForced(dof,:,id(1)-1);
                    a = model.results.dynamicAccelForced(dof,:,id(1)-1);
                end
            end

            % Get time interval limits to be displayed
            ti = str2double(get(gui.edit_Ti,'string'));
            tf = str2double(get(gui.edit_Tf,'string'));
            
            % Get color
            c = tableData{id(1),3};
            if strcmp(c(1),'[')
                c = [str2double(c(2:5)),str2double(c(7:10)),str2double(c(12:15))];
            end
            
            % Get LineWidth
            lw = tableData{id(1),4};

            % Get LineStyle
            ls = tableData{id(1),5};
            
            % Draw displacement
            axes(gui.axes_Displ);
            grid on
            hold on
            plot(linspace(0,model.t,length(d)),d,'color',c,'LineWidth',lw,'LineStyle',ls,...
                 'Tag','DynamicDispl','UserData',id(1));
            xlim([ti,tf])
            
            % Draw Veloc
            axes(gui.axes_Veloc);
            grid on
            hold on
            plot(linspace(0,model.t,length(v)),v,'color',c,'LineWidth',lw,'LineStyle',ls,...
                 'Tag','DynamicVeloc','UserData',id(1));
            xlim([ti,tf])
            
            % Draw acceleration
            axes(gui.axes_Accel);
            grid on
            hold on
            plot(linspace(0,model.t,length(a)),a,'color',c,'LineWidth',lw,'LineStyle',ls,...
                 'Tag','DynamicAccel','UserData',id(1));
            xlim([ti,tf])
            
            % Check time steps interval for phase plane to be ploted
            if size(d,2) == 2
                ni = 1;
                nf = 2;
            else
                step = model.t / model.n_steps;
                ni = ti/step;
                if rem(ni,step) < (step/2)
                    ni = ni - rem(ni,step) + 1;
                else
                    ni = ni - rem(ni,step) + step + 1;
                end
                ni = round(ni,0);

                nf = tf/step;
                if rem(nf,step) < (step/2)
                    nf = nf - rem(nf,step) + 1;
                else
                    nf = nf - rem(nf,step) + step + 1;
                end
                nf = round(nf,0);
            end
            
            % Compute new limits for phase portrait
            maxD = max(d(ni:nf));
            minD = min(d(ni:nf));
            
            maxV = max(v(ni:nf));
            minV = min(v(ni:nf));

            maxA = max(a(ni:nf));
            minA = min(a(ni:nf));

            plots = findobj('Tag','DynamicPPlan');
            if isempty(plots)
                newXLim = ([minD maxD]);
                newYLim = ([minV maxV]);
                newZLim = ([minA maxA]);
            else
                origX = get(gui.axes_PhasePlane_1,'XLim');
                origY = get(gui.axes_PhasePlane_1,'YLim');
                origZ = get(gui.axes_PhasePlane_1,'ZLim');
                
                newXLim = ([min([origX(1), minD]) max([origX(2), maxD])]);
                newYLim = ([min([origY(1), minV]) max([origY(2), maxV])]);
                newZLim = ([min([origZ(1), minA]) max([origZ(2), maxA])]);
            end
            
            % Check equal limits to avoid erros
            if newXLim(1) == newXLim(2)
                newXLim(1) = newXLim(1) - 1;
                newXLim(2) = newXLim(2) + 1;
            end
            if newYLim(1) == newYLim(2)
                newYLim(1) = newYLim(1) - 1;
                newYLim(2) = newYLim(2) + 1;
            end
            if newZLim(1) == newZLim(2)
                newZLim(1) = newZLim(1) - 1;
                newZLim(2) = newZLim(2) + 1;
            end
            
            % Draw (Veloc x displacement) phase plane diagram
            auxPos = get(gui.axes_PhasePlane_1,'UserData');
            axes(gui.axes_PhasePlane_1);
            plot3(d(ni:nf),v(ni:nf),a(ni:nf),'color',c,'LineWidth',lw,'LineStyle',ls,...
                 'Tag','DynamicPPlan','UserData',id(1));
            hold on
            xlim(newXLim);
            ylim(newYLim);
            zlim(newZLim);
            grid on
            set(gui.axes_PhasePlane_1,'UserData',auxPos);

        case false % Delete this diagram
            delete(findobj('Tag','DynamicDispl','UserData',id(1)));
            delete(findobj('Tag','DynamicVeloc','UserData',id(1)));
            delete(findobj('Tag','DynamicAccel','UserData',id(1)));
            delete(findobj('Tag','DynamicPPlan','UserData',id(1)));
            
            plots = findobj('Tag','DynamicPPlan');
            if isempty(plots)
                maxX =  1;
                minX = -1;
                maxY =  1;
                minY = -1;
                maxZ =  1;
                minZ = -1;
            else
                maxX = max(plots(1).XData);
                minX = min(plots(1).XData);
                maxY = max(plots(1).YData);
                minY = min(plots(1).YData);
                maxZ = max(plots(1).ZData);
                minZ = min(plots(1).ZData);
                for i = 2:length(plots)
                    maxX = max([maxX,max(plots(i).XData)]);
                    minX = min([minX,min(plots(i).XData)]);
                    maxY = max([maxY,max(plots(i).YData)]);
                    minY = min([minY,min(plots(i).YData)]);
                    maxZ = max([maxZ,max(plots(i).ZData)]);
                    minZ = min([minZ,min(plots(i).ZData)]);
                end
            end
            
            % Chech equal limits to avoid erros
            newXLim = ([minX  maxX]);
            if all(newXLim == 0)
                newXLim = [-1 1];
            end
            newYLim = ([minY  maxY]);
            if all(newYLim == 0)
                newYLim = [-1 1];
            end
            newZLim = ([minZ  maxZ]);
            if all(newZLim == 0)
                newZLim = [-1 1];
            end
            
            if newXLim(1) == newXLim(2)
                newXLim(1) = newXLim(1) - 1;
                newXLim(2) = newXLim(2) + 1;
            end
            if newYLim(1) == newYLim(2)
                newYLim(1) = newYLim(1) - 1;
                newYLim(2) = newYLim(2) + 1;
            end
            if newZLim(1) == newZLim(2)
                newZLim(1) = newZLim(1) - 1;
                newZLim(2) = newZLim(2) + 1;
            end
            
            auxPos = get(gui.axes_PhasePlane_1,'UserData');
            axes(gui.axes_PhasePlane_1);
            hold on
            xlim(newXLim);
            ylim(newYLim);
            zlim(newZLim);
            grid on
            set(gui.axes_PhasePlane_1,'UserData',auxPos);
    end
    
    % Reset axes y-limits
    for ax = [ gui.axes_Displ, gui.axes_Veloc, gui.axes_Accel ]
        % Get list of handles to graphic objects drawn on ax
        plots = ax.Children;
        if ~isempty(plots)
            ymax = max(plots(1).YData);
            ymin = min(plots(1).YData);
            for i = 2:length(plots)
                ymax = max([ymax,max(plots(i).YData)]);
                ymin = min([ymin,min(plots(i).YData)]);
            end
            if ymax == ymin
                ylim(ax, ylim);
            else
                Ly = (ymax - ymin) * 1.14;
                My = (ymax + ymin) * 0.5;
                ylim(ax, My + Ly*[-0.5,0.5]);
            end
        end
    end

elseif id(2) == 3 % User changed color
    % Check if user is trying to enter an rgb vector
    str = tableData{id(1),id(2)};
    if strcmp(str,'New Color')
        % Get terms of rgb color vecotr entered by user
        %answer = inputdlg({'Red','Green','Blue'},'RGB',[1 30],{'0','0','0'});
        clr = uisetcolor;
        
        % Verify if there was an input 
        if ~clr
            % Get diagram color to reset popupmenu on table
            displ = findobj('Tag','DynamicDispl','UserData',id(1));
            c = get(displ,'color');
            
            % Check which color it is
            if all(c == [1 0 0])
                tableData(id(1),3) = {'Red'};
            elseif all(c == [0 1 0])
                tableData(id(1),3) = {'Green'};
            elseif all(c == [0 0 1])
                tableData(id(1),3) = {'Blue'};
            elseif all(c == [0 0 0])
                tableData(id(1),3) = {'Black'};
            elseif all(c == [1 0 1])
                tableData(id(1),3) = {'Magenta'};
            elseif all(c == [0 1 1])
                tableData(id(1),3) = {'Cyan'};
            else
                tableData(id(1),3) = {sprintf('[%.2f,%.2f,%.2f]',c(1),c(2),c(3))};
            end
            % Reset popupmenu on table
            set(gui.uitable_PlotModes,'Data',tableData)
            return
        end
        
        % Update table data
        clr_txt = sprintf('[%.2f,%.2f,%.2f]',clr(1),clr(2),clr(3));
        tableData(id(1),3) = {clr_txt};
        
        % Update color popupmenu inside table
        formats = get(gui.uitable_PlotModes,'ColumnFormat');
        popup = formats{3};
        popup{end} = clr_txt;
        popup{end + 1} = 'New Color';
        formats{3} = popup;
        
        set(gui.uitable_PlotModes,'Data',tableData,'ColumnFormat',formats)
    elseif strcmp(str(1),'[')
        % Get selected color
        clr = [str2double(str(2:5)),str2double(str(7:10)),str2double(str(12:15))];
    else
        % Get selected color
        clr = tableData{id(1),id(2)};
    end
    
    % Update diagrams colors
    displ = findobj('Tag','DynamicDispl','UserData',id(1));
    veloc = findobj('Tag','DynamicVeloc','UserData',id(1));
    accel = findobj('Tag','DynamicAccel','UserData',id(1));
    pplan = findobj('Tag','DynamicPPlan','UserData',id(1));
    set(displ,'color',clr);
    set(veloc,'color',clr);
    set(accel,'color',clr);
    set(pplan,'color',clr);
    
elseif id(2) == 4 % User changed LineWidth
    % Get input
    width = tableData{id(1),id(2)};
    
    % Get handle to diagrams' plots
    displ = findobj('Tag','DynamicDispl','UserData',id(1));
    veloc = findobj('Tag','DynamicVeloc','UserData',id(1));
    accel = findobj('Tag','DynamicAccel','UserData',id(1));
    pplan = findobj('Tag','DynamicPPlan','UserData',id(1));
    
    % Check validity
    if isnan(width) || width > 5 || width <= 0
        % Reset tableData
        tableData(id(1),id(2)) = {get(displ,'LineWidth')};
        set(gui.uitable_PlotModes,'Data',tableData);
        return
    end
    
    % Update diagrams LineWidth
    set(displ,'LineWidth',width);
    set(veloc,'LineWidth',width);
    set(accel,'LineWidth',width);
    set(pplan,'LineWidth',width);
    
elseif id(2) == 5 % User changed LineStyle
    % Get input
    style = tableData{id(1),id(2)};
    
    % Update diagrams LineStyle
    displ = findobj('Tag','DynamicDispl','UserData',id(1));
    veloc = findobj('Tag','DynamicVeloc','UserData',id(1));
    accel = findobj('Tag','DynamicAccel','UserData',id(1));
    pplan = findobj('Tag','DynamicPPlan','UserData',id(1));
    set(displ,'LineStyle',style);
    set(veloc,'LineStyle',style);
    set(accel,'LineStyle',style);
    set(pplan,'LineStyle',style);
end

%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = GUI_DynamicResults_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Nodes.
function popupmenu_Nodes_Callback(hObject, ~, handles)
include_constants;

% Delete any previous snap points drawing
delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Get model object from root
model = getappdata(0,'model');

% Get node and dof IDs
n = get(hObject,'value');
dof = model.ID(get(handles.popupmenu_DOF,'value'),n);

% TEMPORARY: CHECK IF SOLUTION IS ON FREQ DOMAIN
if get(handles.popupmenu_Domain,'Value') == 2
    set(handles.popupmenu_Domain,'Value',1);
    popupmenu_Domain_Callback(handles.popupmenu_Domain,[],handles);
end

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);    
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:); 
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = vertcat(sum(aux_d,2)',aux_d');
    v = vertcat(sum(aux_v,2)',aux_v');
    a = vertcat(sum(aux_a,2)',aux_a');
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
    
elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    % Reorganize results. Obs: First row is sum of all modes
    d = zeros(model.n_modes+1,2);
    v = zeros(model.n_modes+1,2);
    a = zeros(model.n_modes+1,2);
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
else
    % Only diagrams for all modes together can be plotted
    modes = true;
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
end

% Get time interval limits to be displayed
ti = str2double(get(handles.edit_Ti,'string'));
tf = str2double(get(handles.edit_Tf,'string'));

% Determine if dof is displ or rot
rotationalDOF = false;
switch model.anm.analysis_type
    case FRAME2D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') == 3
            rotationalDOF = true;
        end
    case GRILLAGE_ANALYSIS
        if get(handles.popupmenu_DOF,'value') <= 2
            rotationalDOF = true;
        end
    case FRAME3D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') >= 4
            rotationalDOF = true;
        end
end

updateCanvases(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF);

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_DOF.
function popupmenu_DOF_Callback(hObject, ~, handles)
include_constants;

% Delete any previous snap points drawing
delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Get model object from root
model = getappdata(0,'model');

% Get node and dof IDs
n = get(handles.popupmenu_Nodes,'value');
dof = model.ID(get(hObject,'value'),n);

% TEMPORARY: CHECK IF SOLUTION IS ON FREQ DOMAIN
if get(handles.popupmenu_Domain,'Value') == 2
    set(handles.popupmenu_Domain,'Value',1);
    popupmenu_Domain_Callback(handles.popupmenu_Domain,[],handles);
end

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);    
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:);    
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = vertcat(sum(aux_d,2)',aux_d');
    v = vertcat(sum(aux_v,2)',aux_v');
    a = vertcat(sum(aux_a,2)',aux_a');
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);

elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    % Reorganize results. Obs: First row is sum of all modes
    d = zeros(model.n_modes+1,2);
    v = zeros(model.n_modes+1,2);
    a = zeros(model.n_modes+1,2);
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
else
    % Only diagrams for all modes together can be plotted
    modes = true;
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
end

% Get time interval limits to be displayed
ti = str2double(get(handles.edit_Ti,'string'));
tf = str2double(get(handles.edit_Tf,'string'));

% Determine if dof is displ or rot
rotationalDOF = false;
switch model.anm.analysis_type
    case FRAME2D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') == 3
            rotationalDOF = true;
        end
    case GRILLAGE_ANALYSIS
        if get(handles.popupmenu_DOF,'value') <= 2
            rotationalDOF = true;
        end
    case FRAME3D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') >= 4
            rotationalDOF = true;
        end
end

updateCanvases(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF);

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Result.
function popupmenu_Result_Callback( ~, ~, handles)
popupmenu_DOF_Callback(handles.popupmenu_DOF, [], handles);

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Domain.
function popupmenu_Domain_Callback(hObject, ~, handles)
include_constants;

% Delete any previous snap points drawing
delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Get model object from root
model = getappdata(0,'model');

% Check if selection of domain was time
if get(hObject,'Value') == 1
    set(handles.text_Ti,'String','Ti [s]');
    set(handles.text_Tf,'String','Tf [s]');
    set(handles.edit_Ti,'string','0');
    set(handles.edit_Tf,'string',num2str(model.t));
    popupmenu_DOF_Callback(handles.popupmenu_DOF, [], handles);
    return
end

% Set text boxes for frequency
set(handles.text_Ti,'String','Fi [Hz]');
set(handles.text_Tf,'String','Ff [Hz]');
set(handles.edit_Ti,'string','0');
set(handles.edit_Tf,'string',num2str(floor((model.n_steps/model.t)/2)));

% Get node and dof IDs
n = get(handles.popupmenu_Nodes,'value');
dof = model.ID(get(handles.popupmenu_DOF,'value'),n);

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);    
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:);    
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = sum(aux_d,2)';
    v = sum(aux_v,2)';
    a = sum(aux_a,2)';
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData{1,3};
    
    % Get flags for LineWidth
    width = tableData{1,4};
    
    % Get flags for LineStyle
    style = tableData{1,5};

else
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData{1,3};
    
    % Get flags for LineWidth
    width = tableData{1,4};
    
    % Get flags for LineStyle
    style = tableData{1,5};
end

% Bug correction
if strcmp(width,'Default')
    width = 1.05;
end
if strcmp(style,'Default')
    style = '-';
end

% Get frequency interval limits to be displayed
fi = str2double(get(handles.edit_Ti,'string'));
ff = str2double(get(handles.edit_Tf,'string'));

% Compute number of points for ffr
npts = 2^nextpow2(size(d,2));

% Compute sample frequency
fs = model.n_steps/model.t;

% Compute fft of solution on time domain for displacement
Y = fft(d,npts);
D2 = abs(Y/size(d,2));
D1 = 2 * D2(1:npts/2+1);
fD = 0:fs/npts:fs*(1/2 - 1/npts);

maxD = max(max(D1));
minD = min(min(D1));
if minD == maxD
    maxD = maxD + 1;
    minD = minD - 1;
elseif isempty(minD) || isempty(maxD)
    maxD =  1;
    minD = -1;
else
    medD = (maxD + minD) * 0.5;
    intD = maxD - minD;
    maxD = medD + 0.6*intD;
    minD = medD - 0.6*intD;
end

% Draw (displacement x freq) diagram
axes(handles.axes_Displ);
delete(findobj('Tag','DynamicDispl'))
grid on
hold on
xlabel('f [Hz]');
ylabel('');
if strcmp(clr,'Default')
    c = [0.75,0,0];
else
    c = clr;
end
plot(fD,D1(1:npts/2),'color',c,'LineWidth',width,'LineStyle',style,...
     'Tag','DynamicDispl','UserData',1);
xlim([fi, ff]);
ylim([minD, maxD]);

% Compute fft of solution on time domain for velocity
Y = fft(v,npts);
V2 = abs(Y/size(v,2));
V1 = 2 * V2(1:npts/2+1);
fV = 0:fs/npts:fs*(1/2 - 1/npts);

maxV = max(max(V1));
minV = min(min(V1));
if minV == maxV
    maxV = maxV + 1;
    minV = minV - 1;
elseif isempty(minV) || isempty(maxV)
    maxV =  1;
    minV = -1;
else
    medV = (maxV + minV) * 0.5;
    intV = maxV - minV;
    maxV = medV + 0.6*intV;
    minV = medV - 0.6*intV;
end

% Draw (velocity x freq) diagram
axes(handles.axes_Veloc);
delete(findobj('Tag','DynamicVeloc'))
grid on
hold on
xlabel('f [Hz]');
ylabel('');
if strcmp(clr,'Default')
    c = [0,0.6,0];
else
    c = clr;
end
plot(fV,V1(1:npts/2),'color',c,'LineWidth',width,'LineStyle',style,...
     'Tag','DynamicVeloc','UserData',1);
xlim([fi, ff]);
ylim([minV, maxV]);

% Compute fft of solution on time domain for acceleration
Y = fft(a,npts);
A2 = abs(Y/size(a,2));
A1 = 2 * A2(1:npts/2+1);
fA = 0:fs/npts:fs*(1/2 - 1/npts);

maxA = max(max(A1));
minA = min(min(A1));
if minA == maxA
    maxA = maxA + 1;
    minA = minA - 1;
elseif isempty(minA) || isempty(maxA)
    maxA =  1;
    minA = -1;
else
    medA = (maxA + minA) * 0.5;
    intA = maxA - minA;
    maxA = medA + 0.6*intA;
    minA = medA - 0.6*intA;
end

% Draw (acceleration x freq) diagram
axes(handles.axes_Accel);
delete(findobj('Tag','DynamicAccel'))
grid on
hold on
xlabel('f [Hz]');
ylabel('');
if strcmp(clr,'Default')
    c = [0,0,0.8];
else
    c = clr;
end
plot(fA,A1(1:npts/2),'color',c,'LineWidth',width,'LineStyle',style,...
     'Tag','DynamicAccel','UserData',1);
xlim([fi, ff]);
ylim([minA, maxA]);


%--------------------------------------------------------------------------

% --- Executes after string edition on editable textbox edit_Ti.
function edit_Ti_Callback(hObject, ~, handles)
model = getappdata(0,'model');

% Get current diagram initial time being displayed on canvas
current_T = get(handles.axes_Displ,'XLim');
current_Ti = current_T(1);
current_Tf = current_T(2);

% Get user input
txt_Ti = str2double(get(hObject,'string'));

% Check if input is valid
if isnan(txt_Ti) || txt_Ti >= current_Tf
    set(hObject,'string',num2str(current_Ti));
    return
elseif txt_Ti < 0
    txt_Ti = 0;
    set(hObject,'string','0');
end

% Delete any previous snap points drawing
delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Get plots data
plotsDispl = findobj('Tag','DynamicDispl');
plotsVeloc = findobj('Tag','DynamicVeloc');
plotsAccel = findobj('Tag','DynamicAccel');

% Compute step interval being shown
if length(plotsDispl(1).YData) == 2
    ni = 1;
    nf = 2;
elseif get(handles.popupmenu_Domain,'Value') == 2 % Freq
    disp_xlim = plotsDispl(1).XData;
    [~,ni] = min(abs(disp_xlim-txt_Ti));
    [~,nf] = min(abs(disp_xlim-current_Tf));
else
    step = model.t / model.n_steps;
    ni = txt_Ti/step;
    if rem(ni,step) < (step/2)
        ni = ni - rem(ni,step) + 1;
    else
        ni = ni - rem(ni,step) + step + 1;
    end
    ni = round(ni,0);
    
    nf = current_Tf/step;
    if rem(nf,step) < (step/2)
        nf = nf - rem(nf,step) + 1;
    else
        nf = nf - rem(nf,step) + step + 1;
    end
    nf = round(nf,0);
end

% Compute max and min YData
for i = 1:length(plotsDispl)
    maxYD = max(plotsDispl(i).YData(ni:nf));
    minYD = min(plotsDispl(i).YData(ni:nf));
    if minYD == maxYD
        maxYD = maxYD + 1;
        minYD = minYD - 1;
    elseif isempty(minYD) || isempty(maxYD)
        maxYD =  1;
        minYD = -1;
    else
        medYD = (maxYD + minYD) * 0.5;
        intD = maxYD - minYD;
        maxYD = medYD + 0.57*intD;
        minYD = medYD - 0.57*intD;
    end
end

for i = 1:length(plotsVeloc)
    maxYV = max(plotsVeloc(i).YData(ni:nf));
    minYV = min(plotsVeloc(i).YData(ni:nf));
    if minYV == maxYV
        maxYV = maxYV + 1;
        minYV = minYV - 1;
    elseif isempty(minYV) || isempty(maxYV)
        maxYV =  1;
        minYV = -1;
    else
        medYV = (maxYV + minYV) * 0.5;
        intV = maxYV - minYV;
        maxYV = medYV + 0.57*intV;
        minYV = medYV - 0.57*intV;
    end
end

for i = 1:length(plotsAccel)
    maxYA = max(plotsAccel(i).YData(ni:nf));
    minYA = min(plotsAccel(i).YData(ni:nf));
    if minYA == maxYA
        maxYA = maxYA + 1;
        minYA = minYA - 1;
    elseif isempty(minYA) || isempty(maxYA)
        maxYA =  1;
        minYA = -1;
    else
        medYA = (maxYA + minYA) * 0.5;
        intA = maxYA - minYA;
        maxYA = medYA + 0.57*intA;
        minYA = medYA - 0.57*intA;
    end
end

% Reset axes limits
axes(handles.axes_Displ);
xlim([txt_Ti, current_Tf]);
ylim([minYD, maxYD]);
axes(handles.axes_Veloc);
xlim([txt_Ti, current_Tf]);
ylim([minYV, maxYV]);
axes(handles.axes_Accel);
xlim([txt_Ti, current_Tf]);
ylim([minYA, maxYA]);

% TEMPORARY: CHECK IF DOMAIN IS FREQUENCY
if get(handles.popupmenu_Domain,'Value') == 2
    return
end

% UPDATE PHASE PLANE-------------------------------------------------------
include_constants;

% Get model object from root
model = getappdata(0,'model');

% Get node and dof IDs
n = get(handles.popupmenu_Nodes,'value');
dof = model.ID(get(handles.popupmenu_DOF,'value'),n);

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:);
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = vertcat(sum(aux_d,2)',aux_d');
    v = vertcat(sum(aux_v,2)',aux_v');
    a = vertcat(sum(aux_a,2)',aux_a');
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);

elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    % Reorganize results. Obs: First row is sum of all modes
    d = zeros(model.n_modes+1,2);
    v = zeros(model.n_modes+1,2);
    a = zeros(model.n_modes+1,2);
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
else
    % Only diagrams for all modes together can be plotted
    modes = true;
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
end

% Get initial and final time
ti = txt_Ti;
tf = current_Tf;

% Determine if dof is displ or rot
rotationalDOF = false;
switch model.anm.analysis_type
    case FRAME2D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') == 3
            rotationalDOF = true;
        end
    case GRILLAGE_ANALYSIS
        if get(handles.popupmenu_DOF,'value') <= 2
            rotationalDOF = true;
        end
    case FRAME3D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') >= 4
            rotationalDOF = true;
        end
end

updatePhasePortrait(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF);

%--------------------------------------------------------------------------

% --- Executes after string edition on editable textbox edit_Tf.
function edit_Tf_Callback(hObject, ~, handles)
% Get model object from root
model = getappdata(0,'model');

% Get current diagram initial time being displayed on canvas
current_T = get(handles.axes_Displ,'XLim');
current_Ti = current_T(1);
current_Tf = current_T(2);

% Get user input
txt_Tf = str2double(get(hObject,'string'));

% Check if input is valid
if isnan(txt_Tf) || txt_Tf <= current_Ti
    set(hObject,'string',num2str(current_Tf));
    return
elseif txt_Tf > model.t && get(handles.popupmenu_Domain,'Value') == 1
    txt_Tf = model.t;
    set(hObject,'string',num2str(model.t));
elseif txt_Tf > floor((model.n_steps/model.t)/2) && get(handles.popupmenu_Domain,'Value') == 2
    txt_Tf = floor((model.n_steps/model.t)/2);
    set(hObject,'string',num2str(floor(txt_Tf)));
end

% Delete any previous snap points drawing
delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Get plots data
plotsDispl = findobj('Tag','DynamicDispl');
plotsVeloc = findobj('Tag','DynamicVeloc');
plotsAccel = findobj('Tag','DynamicAccel');

% Compute step interval being shown
if length(plotsDispl(1).YData) == 2
    ni = 1;
    nf = 2;
elseif get(handles.popupmenu_Domain,'Value') == 2 % Freq
    disp_xlim = plotsDispl(1).XData;
    [~,ni] = min(abs(disp_xlim-current_Ti));
    [~,nf] = min(abs(disp_xlim-txt_Tf));
else
    step = model.t / model.n_steps;
    ni = current_Ti/step;
    if rem(ni,step) < (step/2)
        ni = ni - rem(ni,step) + 1;
    else
        ni = ni - rem(ni,step) + step + 1;
    end
    ni = round(ni,0);
    
    nf = txt_Tf/step;
    if rem(nf,step) < (step/2)
        nf = nf - rem(nf,step) + 1;
    else
        nf = nf - rem(nf,step) + step + 1;
    end
    nf = round(nf,0);
end

% Compute max and min YData
for i = 1:length(plotsDispl)
    maxYD = max(plotsDispl(i).YData(ni:nf));
    minYD = min(plotsDispl(i).YData(ni:nf));
    if minYD == maxYD
        maxYD = maxYD + 1;
        minYD = minYD - 1;
    elseif isempty(minYD) || isempty(maxYD)
        maxYD =  1;
        minYD = -1;
    else
        medYD = (maxYD + minYD) * 0.5;
        intD = maxYD - minYD;
        maxYD = medYD + 0.57*intD;
        minYD = medYD - 0.57*intD;
    end
end

for i = 1:length(plotsVeloc)
    maxYV = max(plotsVeloc(i).YData(ni:nf));
    minYV = min(plotsVeloc(i).YData(ni:nf));
    if minYV == maxYV
        maxYV = maxYV + 1;
        minYV = minYV - 1;
    elseif isempty(minYV) || isempty(maxYV)
        maxYV =  1;
        minYV = -1;
    else
        medYV = (maxYV + minYV) * 0.5;
        intV = maxYV - minYV;
        maxYV = medYV + 0.57*intV;
        minYV = medYV - 0.57*intV;
    end
end

for i = 1:length(plotsAccel)
    maxYA = max(plotsAccel(i).YData(ni:nf));
    minYA = min(plotsAccel(i).YData(ni:nf));
    if minYA == maxYA
        maxYA = maxYA + 1;
        minYA = minYA - 1;
    elseif isempty(minYA) || isempty(maxYA)
        maxYA =  1;
        minYA = -1;
    else
        medYA = (maxYA + minYA) * 0.5;
        intA = maxYA - minYA;
        maxYA = medYA + 0.57*intA;
        minYA = medYA - 0.57*intA;
    end
end

% Reset axes limits
axes(handles.axes_Displ);
xlim([current_Ti, txt_Tf]);
ylim([minYD, maxYD]);
axes(handles.axes_Veloc);
xlim([current_Ti, txt_Tf]);
ylim([minYV, maxYV]);
axes(handles.axes_Accel);
xlim([current_Ti, txt_Tf]);
ylim([minYA, maxYA]);

% TEMPORARY: CHECK IF DOMAIN IS FREQUENCY
if get(handles.popupmenu_Domain,'Value') == 2
    return
end

% UPDATE PHASE PLANE-------------------------------------------------------
include_constants;

% Get node and dof IDs
n = get(handles.popupmenu_Nodes,'value');
dof = model.ID(get(handles.popupmenu_DOF,'value'),n);

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:);
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = vertcat(sum(aux_d,2)',aux_d');
    v = vertcat(sum(aux_v,2)',aux_v');
    a = vertcat(sum(aux_a,2)',aux_a');
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);

elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    % Reorganize results. Obs: First row is sum of all modes
    d = zeros(model.n_modes+1,2);
    v = zeros(model.n_modes+1,2);
    a = zeros(model.n_modes+1,2);
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
else
    % Only diagrams for all modes together can be plotted
    modes = true;
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for colors
    clr = tableData(:,3);
    
    % Get flags for LineWidth
    width = tableData(:,4);
    
    % Get flags for LineStyle
    style = tableData(:,5);
end

% Get initial and final time
ti = current_Ti;
tf = txt_Tf;

% Determine if dof is displ or rot
rotationalDOF = false;
switch model.anm.analysis_type
    case FRAME2D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') == 3
            rotationalDOF = true;
        end
    case GRILLAGE_ANALYSIS
        if get(handles.popupmenu_DOF,'value') <= 2
            rotationalDOF = true;
        end
    case FRAME3D_ANALYSIS
        if get(handles.popupmenu_DOF,'value') >= 4
            rotationalDOF = true;
        end
end

updatePhasePortrait(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF);

%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_Phase.
function popupmenu_Phase_Callback(hObject, ~, handles)
cam = getappdata(0,'phaseCam');
if ~isempty(cam)
    setappdata(0,'phaseCam',[]);
end

val = get(hObject,'Value');
axes(handles.axes_PhasePlane_1);
if val == 1 % 3D
    view(3);
elseif val == 2 % Displ x Veloc
    view(0,90);
elseif val == 3 % Displ x Accel
    view(0,0);
else % Veloc x Accel
    view(90,0);
end

%--------------------------------------------------------------------------
% --- Executes by pressing pushbutton_ExportValues
function pushbutton_ExportValues_Callback(hObject, ~, handles)
include_constants;

% Disable button while tasks are being performed
set(hObject,'Enable','Off');

filterspec = {'*.txt';'*.xlsx'};
DialogTitle = 'LESM - Export Result Values';
DefaultName = 'untitled';
[filename,pathname,FilterIndex] = uiputfile(filterspec,DialogTitle,DefaultName);

if ~filename
    % Enable button for future use
    set(hObject,'Enable','On');
    return
end
fullname = strcat(pathname,filename);

model = getappdata(0,'model');

% Get node and dof IDs
n = get(handles.popupmenu_Nodes,'value');
dofVal = get(handles.popupmenu_DOF,'value');
dof = model.ID(dofVal,n);

dofString = get(handles.popupmenu_DOF,'string');
dofLabel = char(dofString(dofVal,:));

% Get displacement, Veloc and acceleration vectors
if dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    d = zeros(1,2);
    v = zeros(1,2);
    a = zeros(1,2);
elseif get(handles.popupmenu_Result,'value') == 1 % Total results
    d = model.results.dynamicDispl(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:);
elseif get(handles.popupmenu_Result,'value') == 2 % Free vibration
    d = model.results.dynamicDispl(dof,:,:) - model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVeloc(dof,:,:) - model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccel(dof,:,:) - model.results.dynamicAccelForced(dof,:,:);
else % get(handles.popupmenu_Result,'value') == 3 % Forced vibration
    d = model.results.dynamicDisplForced(dof,:,:);
    v = model.results.dynamicVelocForced(dof,:,:);
    a = model.results.dynamicAccelForced(dof,:,:);
end

% Check if analysis was uncoupled (results separated by mode contribution)
if model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && (dof <= (model.neq - model.neqfixed) || model.nodes(n).isInclinedSupp)
    % Initialize auxiliar matrices
    aux_d(:,:) = d;
    aux_v(:,:) = v;
    aux_a(:,:) = a;
    
    % Concatenate results. Obs: First row is sum of all modes
    d = vertcat(sum(aux_d,2)',aux_d');
    v = vertcat(sum(aux_v,2)',aux_v');
    a = vertcat(sum(aux_a,2)',aux_a');
    
    % Clear auxiliar variables
    clear aux_d
    clear aux_v
    clear aux_a
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    if all(modes == false)
        modes = true;
    end
    
elseif model.results.type == DYNAMIC_MODALSUP_LINEAR && model.n_modes > 1 && dof > (model.neq - model.neqfixed) && ~model.nodes(n).isInclinedSupp
    % Reorganize results. Obs: First row is sum of all modes
    d = zeros(model.n_modes+1,model.n_steps + 1);
    v = zeros(model.n_modes+1,model.n_steps + 1);
    a = zeros(model.n_modes+1,model.n_steps + 1);
    
    % Get table data
    tableData = get(handles.uitable_PlotModes,'Data');
    
    % Get flags for which diagrams should be plotted
    modes = zeros(1,model.n_modes+1);
    for i = 1:model.n_modes+1
        modes(i) = tableData{i,2};
    end
    if all(modes == false)
        modes = true;
    end
else
    % Only diagrams for all modes together can be plotted
    modes = true;
end
modes = find(modes);

% Compute number of points for ffr
npts = 2^nextpow2(size(d,2));

% Compute sample frequency
fs = model.n_steps/model.t;

% Compute range of frequency
ff = 0:fs/npts:fs*(1/2 - 1/npts);

% Compute fft of solution on time domain for displacement
Y = fft(d(1,:),npts);
D2 = abs(Y/size(d,2));
D1 = 2 * D2(1:npts/2);

% Compute fft of solution on time domain for velocity
Y = fft(v(1,:),npts);
V2 = abs(Y/size(v,2));
V1 = 2 * V2(1:npts/2);

% Compute fft of solution on time domain for acceleration
Y = fft(a(1,:),npts);
A2 = abs(Y/size(a,2));
A1 = 2 * A2(1:npts/2);

if FilterIndex == 1 % .txt
    % Create folder for text files to be saved in
    status = mkdir(strcat(pathname,filename(1:end-4)));
    if ~status
        set(hObject,'Enable','On');
        msgbox('Failed to create new folder for text files.','Error','error');
        return
    end
    
    % Adjust fullname to created folder path
    fullname = strcat(fullname(1:end-4),'\',filename);
    
    units = 'm';
    switch model.anm.analysis_type
        case FRAME2D_ANALYSIS
            if get(handles.popupmenu_DOF,'value') == 3
                units = 'rad';
            end
        case GRILLAGE_ANALYSIS
            if get(handles.popupmenu_DOF,'value') <= 2
                units = 'rad';
            end
        case FRAME3D_ANALYSIS
            if get(handles.popupmenu_DOF,'value') >= 4
                units = 'rad';
            end
    end
    
    % File names
    files = cell(length(modes)+1,1);
    for i = 1:length(modes)
        if (modes(i) == 1) % all modes
            files(i) = {horzcat(fullname(1:end-4),'_AllModes',fullname(end-3:end))};
        else
            files(i) = {horzcat(fullname(1:end-4),sprintf('_Mode_%i',modes(i)-1),fullname(end-3:end))};
        end
    end
    files(length(modes)+1) = {horzcat(fullname(1:end-4),'_fft',fullname(end-3:end))};
    
    % Print files
    h = waitbar(0,'Exporting data to txt...','WindowStyle','modal');
    tTime = linspace(0,model.t,model.n_steps+1)';
    counter = 0;
    
    for m = modes
        counter = counter + 1;
        waitbar(counter/(length(modes)+2));
        
        % Open file
        txt = fopen(files{counter},'wt');
        
        fprintf(txt, '\n=========================================================\n');
        fprintf(txt, ' LESM - Linear Elements Structure Model analysis program\n');
        fprintf(txt, '                 DYNAMIC ANALYSIS RESULTS\n');
        fprintf(txt, '\n');
        fprintf(txt, sprintf(' Node: %i\n',n));
        fprintf(txt, sprintf(' DOF : %s\n',dofLabel));
        if m == 1
            fprintf(txt, ' All modes\n');
        else
            fprintf(txt, sprintf(' Mode: %i (%.5e Hz)\n',m-1,model.W(m-1)/(2*pi)));
        end
        fprintf(txt, '=========================================================\n');
        fprintf(txt, sprintf('   TIME [S] | DISPL [%s] | VELOC[%s/s] | ACCEL[%s/s2]\n',units,units,units));
        for i = 1:model.n_steps+1
            fprintf(txt, sprintf(' %.5e  %.5e  %.5e  %.5e\n',tTime(i),d(m,i),v(m,i),a(m,i)));
        end

        % Close file
        fclose(txt);
    end
    
    counter = counter + 1;
    waitbar(counter/(length(modes)+2));
    
    % Open file
    txt = fopen(files{end},'wt');
    
    fprintf(txt, '\n=========================================================\n');
    fprintf(txt, ' LESM - Linear Elements Structure Model analysis program\n');
    fprintf(txt, '                 DYNAMIC ANALYSIS RESULTS\n');
    fprintf(txt, '\n');
    fprintf(txt, sprintf(' Node: %i\n',n));
    fprintf(txt, sprintf(' DOF : %s\n',dofLabel));
    fprintf(txt, ' FFT (solution on frequency domain)\n');
    fprintf(txt, '=========================================================\n');
    fprintf(txt, sprintf('   FREQ [Hz] | DISPL [%s] | VELOC[%s/s] | ACCEL[%s/s2]\n',units,units,units));
    for i = 1:npts/2
        fprintf(txt, sprintf(' %.5e  %.5e  %.5e  %.5e\n',ff(i),D1(i),V1(i),A1(i)));
    end
    
    % Close file
    fclose(txt);
    
    % Close waitbar
    waitbar(1);
    close(h);
    
    % Enable button for future use
     set(hObject,'Enable','On');
    
else % FilterIndex == 2  % .xlsx
    
    tNode = cell(model.n_steps + 1,1);
    tNode(1) = {n};
    
    tDof = tNode;
    tDof(1) = {dofLabel};
    
    tMode  = cell(model.n_steps + 1,1);
    
    tTime = linspace(0,model.t,model.n_steps+1)';
    
    h = waitbar(0,'Exporting data to spreadsheet...','WindowStyle','modal');
    
    counter = 0;
    for m = modes
        counter = counter + 1;
        waitbar(counter/(length(modes)+2));
        
        if m == 1
            tMode(1) = {'All'};
            tMode(2) = {[]};
        else
            tMode(1) = {m-1};
            tMode(2) = {sprintf('%.4e Hz',model.W(m-1)/(2*pi))};
        end
        tDispl = d(m,:)';
        tVeloc = v(m,:)';
        tAccel = a(m,:)';
        
        T = table(tNode,tDof,tMode,tTime,tDispl,tVeloc,tAccel,'VariableNames',{'Node';'DOF';'Mode';'Time';'Displ';'Veloc';'Accel'});
        try
            writetable(T,fullname,'Sheet',counter,'Range','B2');
        catch 
            close(h);
            set(hObject,'Enable','On');
            msgbox(sprintf('Failed to write mode %i data because Excel is running.',m-1),'Error','error');
            return
        end
    end
    
    counter = counter + 1;
    waitbar(counter/(length(modes)+2));
    
    tNode = cell(npts/2,1);
    tNode(1) = {n};
    
    tDof = tNode;
    tDof(1) = {dofLabel};
    
    tDispl = D1';
    tVeloc = V1';
    tAccel = A1';
    tFreq  = ff';
    
    T = table(tNode,tDof,tFreq,tDispl,tVeloc,tAccel,'VariableNames',{'Node';'DOF';'Freq_Hz';'Displ';'Veloc';'Accel'});
    warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;
    try
        writetable(T,fullname,'Sheet',counter,'Range','B2');
    catch
        msgbox('Failed to write fft data because Excel is running.','Error','error');
    end
    
    waitbar(1);
    close(h);
    
    % Enable button for future use
    set(hObject,'Enable','On');
end

%--------------------------------------------------------------------------
% --- Executes by pressing pushbutton_ExportImages
function pushbutton_ExportImages_Callback(hObject, ~, handles)
% Disable button while tasks are being performed
set(hObject,'Enable','Off');

% Check if anything is drawn on axes
tableData = get(handles.uitable_PlotModes,'Data');
aux = zeros(1,size(tableData,1));
for i = 1:size(tableData,1)
    aux(i) = tableData{i,2};
end
if all(aux == false)
    set(hObject,'Enable','On');
    msgbox('No diagrams are plotted on the canvases!','Warning','warn');
    return
end
clear aux;

% Get handle to GUI_Main
mainDlg = findobj('Tag','GUI_Main');

% Get GUI_Main's name
name = get(mainDlg,'Name');

% Get current model name
name = erase(name,'LESM - Linear Elements Structure Model');
if isempty(name)
    name = 'untitled';
else
    name = erase(name,'-');
end

% Specifiy possible export formats
filterspec = {'*.jpeg';'*.png';'*.eps';'*.emf';'*.svg'};

% Open dialog for user to save images
DialogTitle = 'LESM - Save figure';
[filename,pathname,indx] = uiputfile(filterspec,DialogTitle,name);

% If user provided a valid name
if filename ~= 0
    
    % Get node and dof IDs
    n = get(handles.popupmenu_Nodes,'value');
    dofVal = get(handles.popupmenu_DOF,'value');
    dofString = get(handles.popupmenu_DOF,'string');
    dofLabel = char(dofString(dofVal,:));
    
    % Get time interval limits being displayed
    ti = str2double(get(handles.edit_Ti,'string'));
    tf = str2double(get(handles.edit_Tf,'string'));
    
    % Get result type (total, free, forced)
    if get(handles.popupmenu_Result,'value') == 1 % total
        resLabel = 'Total Response';
    elseif get(handles.popupmenu_Result,'value') == 2 % free
        resLabel = 'Free Vibration';
    else % forced
        resLabel = 'Forced Vibration';
    end
    
    % Get solution domain
    domVal = get(handles.popupmenu_Domain,'Value');
    
    % Compute diagrams titles
    if domVal == 1  % Time domain
        howManyAxes = 4;
        axTitle = sprintf('LESM - Dynamic Analysis - %s - Model: %s - Node: %i  DOF: %s',resLabel,name,n,dofLabel);
        axTitles = { axTitle ;
                     axTitle ;
                     axTitle ;
                     sprintf('LESM - Dynamic Analysis - %s - Model: %s - Node: %i  DOF: %s - ti: %05.2f [s]  tf: %05.2f [s]',resLabel,name,n,dofLabel,ti,tf) };

        % Determine axes legends
        axLegend = cell(1,size(tableData,1));
        count = 0;
        for i = 1:size(tableData,1)
            if tableData{i,2}
                count = count + 1;
                if strcmp(tableData{i,1},'All')
                    axLegend(count) = {'All Modes'};
                else
                    axLegend(count) = {sprintf('Mode %s',tableData{i,1})};
                end
            end
        end
        aux = axLegend(1:count);
        axLegend = aux;
        clear aux;
        
    else % domVal == 2 % Freq domain
        howManyAxes = 3;
        axTitles = { sprintf('LESM - Dynamic Analysis - %s - Model: %s - Node: %i  DOF: %s - Displacement Frequency Spectrum',resLabel,name,n,dofLabel)     ;
                     sprintf('LESM - Dynamic Analysis - %s - Model: %s - Node: %i  DOF: %s - Velocity Frequency Spectrum',resLabel,name,n,dofLabel)         ;
                     sprintf('LESM - Dynamic Analysis - %s - Model: %s - Node: %i  DOF: %s - Acceleration Frequency Spectrum',resLabel,name,n,dofLabel)     };
        axLegend = {'All modes'};
    end
    
    % Compute auxiliar index value, ammount of chars of file format
    if indx == 1
        aux_indx = 5;
    else
        aux_indx = 4;
    end
    
    % Auxiliary formats cell
    formats = {'jpeg';'png';'eps';'meta';'svg'};
    
    % Axes handles array
    axs = [ handles.axes_Displ ;
            handles.axes_Veloc ;
            handles.axes_Accel ;
            handles.axes_PhasePlane_1 ];
    
    % New filename specifiers array
    names = [ '_Displ' ;
              '_Veloc' ;
              '_Accel' ;
              '_Phase' ];
          
    % Check if there is already a folder with this name
    list = dir(pathname);
    aux = [];
    for i = 1:length(list)
        if strcmp(filename(1:end-aux_indx),list(i).name)
            % if a folder with this name already exists, use clock to create a new foldername
            c = clock;
            aux = sprintf('_%02.0f-%02.0f-%04.0f_%02.0f-%02.0f-%02.0f',c(2),c(3),c(1),c(4),c(5),c(6));
            break
        end
    end
    
    % Create folder for image files to be saved in
    status = mkdir(strcat(pathname,filename(1:end-aux_indx),aux));
    if ~status
        set(hObject,'Enable','On');
        msgbox('Failed to create new folder for image files.','Error','error');
        return
    end
    
    % Adjust pathname to created folder path
    pathname = strcat(pathname,filename(1:end-aux_indx),aux,'\');
    
    % Create waitbar
    h = waitbar(0,'Exporting images...','WindowStyle','modal');
    waitbar(1/(howManyAxes+1));
    
    % Get screensize from root
    screensize = get(groot,'Screensize');
    
    % Save each axes on a different file
    for i = 1:howManyAxes
        f_new = figure;
        set(f_new,'Visible','off','Position',screensize);
        ax_new = copyobj(axs(i),f_new);
        title(axTitles{i},'Interpreter','none');
        legend(axLegend,'Location','northeast');
        
        fullname = strcat(pathname,filename(1:end-aux_indx),names(i,:),filename(end-aux_indx+1:end));
        saveas(ax_new,fullname,formats{indx})
        
        close(f_new)
        
        waitbar((i+1)/(howManyAxes+1));
    end
    
    % Close waitbar
    close(h);
end

% Enable button for future use
set(hObject,'Enable','On');

%--------------------------------------------------------------------------
function initCanvases(handles)
% Get handle to model
model = getappdata(0,'model');

% Set axes_Displ
axes(handles.axes_Displ);
xlim([0 model.t]);
grid on
hold on

% Set axes_Veloc
axes(handles.axes_Veloc);
xlim([0 model.t]);
grid on
hold on

% Set axes_Accel
axes(handles.axes_Accel);
xlim([0 model.t]);
grid on
hold on

% Set axes_PhasePlane_1
axes(handles.axes_PhasePlane_1);
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
grid on
handles.axes_PhasePlane_1.Clipping = 'off';

% Compute axes limits within whole dialog
% axes_Displ
axDisplPos = get(handles.axes_Displ,'Position');     % units = normalized
pnDisplPos = get(handles.uipanel_Displ,'Position');  % units = normalized
displX = pnDisplPos(1) + axDisplPos(1) * pnDisplPos(3);
displY = pnDisplPos(2) + axDisplPos(2) * pnDisplPos(4);
axDisplCorners = [displX, displY;...
                  displX + axDisplPos(3) * pnDisplPos(3), displY;...
                  displX + axDisplPos(3) * pnDisplPos(3), displY + axDisplPos(4) * pnDisplPos(4);...
                  displX, displY + axDisplPos(4) * pnDisplPos(4)];
set(handles.axes_Displ,'UserData',axDisplCorners);

% axes_Veloc
axVelocPos = get(handles.axes_Veloc,'Position');     % units = normalized
pnVelocPos = get(handles.uipanel_Veloc,'Position');  % units = normalized
velocX = pnVelocPos(1) + axVelocPos(1) * pnVelocPos(3);
velocY = pnVelocPos(2) + axVelocPos(2) * pnVelocPos(4);
axVelocCorners = [velocX, velocY;...
                  velocX + axVelocPos(3) * pnVelocPos(3), velocY;...
                  velocX + axVelocPos(3) * pnVelocPos(3), velocY + axVelocPos(4) * pnVelocPos(4);...
                  velocX, velocY + axVelocPos(4) * pnVelocPos(4)];
set(handles.axes_Veloc,'UserData',axVelocCorners);

% axes_Accel
axAccelPos = get(handles.axes_Accel,'Position');     % units = normalized
pnAccelPos = get(handles.uipanel_Accel,'Position');  % units = normalized
accelX = pnAccelPos(1) + axAccelPos(1) * pnAccelPos(3);
accelY = pnAccelPos(2) + axAccelPos(2) * pnAccelPos(4);
axAccelCorners = [accelX, accelY;...
                  accelX + axAccelPos(3) * pnAccelPos(3), accelY;...
                  accelX + axAccelPos(3) * pnAccelPos(3), accelY + axAccelPos(4) * pnAccelPos(4);...
                  accelX, accelY + axAccelPos(4) * pnAccelPos(4)];
set(handles.axes_Accel,'UserData',axAccelCorners);

% axes_PhasePlane_1
axPhasePos = get(handles.axes_PhasePlane_1,'Position');% units = normalized
pnPhasePos = get(handles.uipanel_Phase,'Position');    % units = normalized
phaseX = pnPhasePos(1) + axPhasePos(1) * pnPhasePos(3);
phaseY = pnPhasePos(2) + axPhasePos(2) * pnPhasePos(4);
axPhaseCorners = [phaseX, phaseY;...
                  phaseX + axPhasePos(3) * pnPhasePos(3), phaseY;...
                  phaseX + axPhasePos(3) * pnPhasePos(3), phaseY + axPhasePos(4) * pnPhasePos(4);...
                  phaseX, phaseY + axPhasePos(4) * pnPhasePos(4)];
set(handles.axes_PhasePlane_1,'UserData',axPhaseCorners);

%--------------------------------------------------------------------------
function updateCanvases(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF)
% Get handle to model
model = getappdata(0,'model');

% Check time steps interval to be ploted
if size(d,2) == 2
    ni = 1;
    nf = 2;
else
    step = model.t / model.n_steps;
    ni = ti/step;
    if rem(ni,step) < (step/2)
        ni = ni - rem(ni,step) + 1;
    else
        ni = ni - rem(ni,step) + step + 1;
    end
    ni = round(ni,0);
    
    nf = tf/step;
    if rem(nf,step) < (step/2)
        nf = nf - rem(nf,step) + 1;
    else
        nf = nf - rem(nf,step) + step + 1;
    end
    nf = round(nf,0);
end

% Determine limits
index = find(modes);

maxD = max(max(d(index,ni:nf)));
minD = min(min(d(index,ni:nf)));
if minD == maxD
    maxD = maxD + 1;
    minD = minD - 1;
elseif isempty(minD) || isempty(maxD) || isinf(minD) || isinf(maxD)
    maxD =  1;
    minD = -1;
else
    medD = (maxD + minD) * 0.5;
    intD = maxD - minD;
    maxD = medD + 0.57*intD;
    minD = medD - 0.57*intD;
end

maxV = max(max(v(index,ni:nf)));
minV = min(min(v(index,ni:nf)));
if minV == maxV
    maxV = maxV + 1;
    minV = minV - 1;
elseif isempty(minV) || isempty(maxV) || isinf(minV) || isinf(maxV)
    maxV =  1;
    minV = -1;
else
    medV = (maxV + minV) * 0.5;
    intV = maxV - minV;
    maxV = medV + 0.57*intV;
    minV = medV - 0.57*intV;
end

maxA = max(max(a(index,ni:nf)));
minA = min(min(a(index,ni:nf)));
if minA == maxA
    maxA = maxA + 1;
    minA = minA - 1;
elseif isempty(minA) || isempty(maxA) || isinf(minA) || isinf(maxA)
    maxA =  1;
    minA = -1;
else
    medA = (maxA + minA) * 0.5;
    intA = maxA - minA;
    maxA = medA + 0.57*intA;
    minA = medA - 0.57*intA;
end

% Determine string of ylabel for (displacement x time) diagram
if ~rotationalDOF
    str_ylabel = 'displ [m]';
else
    str_ylabel = 'rot [rad]';
end

% Initialize counter
counter = 1;

% Draw (displacement x time) diagram
axes(handles.axes_Displ);
delete(findobj('Tag','DynamicDispl'))
xlabel('t [s]');
ylabel(str_ylabel);
while counter <= length(modes)
    if modes(counter)
        c = clr{counter};
        if strcmp(c,'Default')
            c = [0.75,0,0];
        elseif strcmp(c(1),'[')
            c = [str2double(c(2:5)),str2double(c(7:10)),str2double(c(12:15))];
        end
        lw = width{counter};
        if strcmp(lw,'Default')
            lw = 1.05;
        end
        ls = style{counter};
        if strcmp(ls,'Default')
            ls = '-';
        end
        plot(linspace(0,model.t,length(d(counter,:))),d(counter,:),'color',c,'LineWidth',lw,'LineStyle',ls,...
             'Tag','DynamicDispl','UserData',counter);
    end
    counter  = counter + 1;
end
xlim([ti, tf]);
ylim([minD, maxD]);

% Determine string of ylabel for (Veloc x time) diagram
if ~rotationalDOF
    str_ylabel = 'veloc [m/s]';
else
    str_ylabel = 'veloc [rad/s]';
end

% Initialize counter
counter = 1;

% Draw (Veloc x time) diagram
axes(handles.axes_Veloc);
delete(findobj('Tag','DynamicVeloc'))
xlabel('t [s]');
ylabel(str_ylabel);
while counter <= length(modes)
    if modes(counter)
        c = clr{counter};
        if strcmp(c,'Default')
            c = [0,0.6,0];
        elseif strcmp(c(1),'[')
            c = [str2double(c(2:5)),str2double(c(7:10)),str2double(c(12:15))];
        end
        lw = width{counter};
        if strcmp(lw,'Default')
            lw = 1.05;
        end
        ls = style{counter};
        if strcmp(ls,'Default')
            ls = '-';
        end
        plot(linspace(0,model.t,length(v(counter,:))),v(counter,:),'color',c,'LineWidth',lw,'LineStyle',ls,...
             'Tag','DynamicVeloc','UserData',counter);
    end
    counter  = counter + 1;
end
xlim([ti, tf]);
ylim([minV, maxV]);

% Determine string of ylabel for (acceleration x time) diagram
if ~rotationalDOF
    str_ylabel = 'accel [m/s^2]';
else
    str_ylabel = 'accel [rad/s^2]';
end

% Initialize counter
counter = 1;

% Draw (acceleration x time) diagram
axes(handles.axes_Accel);
delete(findobj('Tag','DynamicAccel'))
xlabel('t [s]');
ylabel(str_ylabel);
while counter <= length(modes)
    if modes(counter)
        c = clr{counter};
        if strcmp(c,'Default')
            c = [0,0,0.8];
        elseif strcmp(c(1),'[')
            c = [str2double(c(2:5)),str2double(c(7:10)),str2double(c(12:15))];
        end
        lw = width{counter};
        if strcmp(lw,'Default')
            lw = 1.05;
        end
        ls = style{counter};
        if strcmp(ls,'Default')
            ls = '-';
        end
        plot(linspace(0,model.t,length(a(counter,:))),a(counter,:),'color',c,'LineWidth',lw,'LineStyle',ls,...
             'Tag','DynamicAccel','UserData',counter);
    end
    counter  = counter + 1;
end
xlim([ti, tf]);
ylim([minA, maxA]);

updatePhasePortrait(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF);

%--------------------------------------------------------------------------
function updatePhasePortrait(handles,d,v,a,modes,clr,width,style,ti,tf,rotationalDOF)
% Get handle to model
model = getappdata(0,'model');

% Check time steps interval to be ploted
if size(d,2) == 2
    ni = 1;
    nf = 2;
else
    step = model.t / model.n_steps;
    ni = ti/step;
    if rem(ni,step) < (step/2)
        ni = ni - rem(ni,step) + 1;
    else
        ni = ni - rem(ni,step) + step + 1;
    end
    ni = round(ni,0);
    
    nf = tf/step;
    if rem(nf,step) < (step/2)
        nf = nf - rem(nf,step) + 1;
    else
        nf = nf - rem(nf,step) + step + 1;
    end
    nf = round(nf,0);
end

% Determine string of xlabel and ylabel for (Veloc x displacement) phase plane diagram
if ~rotationalDOF
    str_xlabel = 'displ [m]';
    str_ylabel = 'veloc [m/s]';
    str_zlabel = 'accel [m/s^2]';
else
    str_xlabel = 'rot [rad]';
    str_ylabel = 'veloc [rad/s]';
    str_zlabel = 'accel [rad/s^2]';
end

% Determine limits
index = find(modes);

maxD = max(max(d(index,ni:nf)));
minD = min(min(d(index,ni:nf)));
if minD == maxD
    maxD = maxD + 1;
    minD = minD - 1;
elseif isempty(minD) || isempty(maxD)
    maxD =  1;
    minD = -1;
end
maxV = max(max(v(index,ni:nf)));
minV = min(min(v(index,ni:nf)));
if minV == maxV
    maxV = maxV + 1;
    minV = minV - 1;
elseif isempty(minV) || isempty(maxV)
    maxV =  1;
    minV = -1;
end
maxA = max(max(a(index,ni:nf)));
minA = min(min(a(index,ni:nf)));
if minA == maxA
    maxA = maxA + 1;
    minA = minA - 1;
elseif isempty(minA) || isempty(maxA)
    maxA =  1;
    minA = -1;
end

% Initialize counter
counter = 1;

% Draw phase portrait diagram
auxPos = get(handles.axes_PhasePlane_1,'UserData');
axes(handles.axes_PhasePlane_1);
delete(findobj('Tag','DynamicPPlan'))
while counter <= length(modes)
    if modes(counter)
        c = clr{counter};
        if strcmp(c,'Default')
            c = [0.6,0.2,0.6];
        elseif strcmp(c(1),'[')
            c = [str2double(c(2:5)),str2double(c(7:10)),str2double(c(12:15))];
        end
        lw = width{counter};
        if strcmp(lw,'Default')
            lw = 1.05;
        end
        ls = style{counter};
        if strcmp(ls,'Default')
            ls = '-';
        end
        plot3(d(counter,ni:nf),v(counter,ni:nf),a(counter,ni:nf),'color',c,'LineWidth',lw,'LineStyle',ls,...
             'Tag','DynamicPPlan','UserData',counter);
        hold on
    end
    counter  = counter + 1;
end
xlim([minD maxD]);
ylim([minV maxV]);
zlim([minA maxA]);
xlabel(str_xlabel);
ylabel(str_ylabel);
zlabel(str_zlabel);
grid on
set(handles.axes_PhasePlane_1,'UserData',auxPos);

%--------------------------------------------------------------------------
% --- Executes when cursor moves on a canvas
function moveCursorFcn(~, ~, ~)
cam = getappdata(0,'phaseCam');
if isempty(cam)
    return
end

% Get dialog
dlg = findobj('Tag','GUI_DynamicResults');

% Get handle to GUI_DynamicResults
gui = guidata(dlg);

% Get cooridnates within dialog
dlgCoords = get(dlg,'CurrentPoint');
if isempty(dlgCoords)
    return
end

% Get which axes was clicked on
axPhasePos = get(gui.axes_PhasePlane_1,'UserData');
if inpolygon(dlgCoords(1),dlgCoords(2),axPhasePos(:,1),axPhasePos(:,2))
    % Catches the actual coordinates on window.
    set(dlg,'Units','pixels');
    wpt = get(dlg, 'CurrentPoint');
    x = wpt(1);
    y = wpt(2);
    set(dlg,'Units','normalized');

    % Calculate the displacements of azimute and elevation.
    dAz = cam.rotIniX - x;
    dEl = cam.rotIniY - y;
    if dAz > -10 && dAz < 10
        dAz = 0;
    end

%     aux_max = max(abs([dAz,dEl]));
%     if aux_max == abs(dAz)
%         dEl = 0;
%     else
%         dAz = 0;
%     end

    % Incorporate the displacements calculated before in the
    % inicial azimute and elevation.
    az = cam.rotIniAz + dAz/2;
    el = cam.rotIniEl + dEl/2;

    % Checking the elevation values and changing them when its
    % necessary.
    if (el > 90)
        el = 90;
    elseif (el < -90)
        el = -90;
    end

    % Sets the new view.
    view(gui.axes_PhasePlane_1,[az,el]);
end

%--------------------------------------------------------------------------
% --- Executes when button is released
function buttonUpFcn(~, ~, ~)
cam = getappdata(0,'phaseCam');
if ~isempty(cam)
    setappdata(0,'phaseCam',[]);
    return
end

%--------------------------------------------------------------------------
% --- Executes when button is pressed down on a canvas
function buttonDownFcn(~, ~, ~)
% Get dialog
dlg = findobj('Tag','GUI_DynamicResults');

% Get handle to GUI_DynamicResults
gui = guidata(dlg);

% Delete any previous snap points drawing
if ~strcmp(get(dlg, 'SelectionType'),'alt')
    delete(findobj('Tag','dynResSnapPointDispl'))
    delete(findobj('Tag','dynResSnapPointVeloc'))
    delete(findobj('Tag','dynResSnapPointAccel'))
    delete(findobj('Tag','dynResSnapPointPPlan'))
end

% Get cooridnates within dialog
dlgCoords = get(dlg,'CurrentPoint');
if isempty(dlgCoords)
    return
end

% Get which axes was clicked on
thisAxes = [];
otherAxes = [];
phaseAxes = [];
axDisplPos = get(gui.axes_Displ,'UserData');
if inpolygon(dlgCoords(1),dlgCoords(2),axDisplPos(:,1),axDisplPos(:,2))
    thisAxes = gui.axes_Displ;
    otherAxes = [ gui.axes_Veloc ; gui.axes_Accel ];
else
    axVelocPos = get(gui.axes_Veloc,'UserData');
    if inpolygon(dlgCoords(1),dlgCoords(2),axVelocPos(:,1),axVelocPos(:,2))
        thisAxes = gui.axes_Veloc;
        otherAxes = [ gui.axes_Displ ; gui.axes_Accel ];
    else
        axAccelPos = get(gui.axes_Accel,'UserData');
        if inpolygon(dlgCoords(1),dlgCoords(2),axAccelPos(:,1),axAccelPos(:,2))
            thisAxes = gui.axes_Accel;
            otherAxes = [ gui.axes_Displ ; gui.axes_Veloc ];
        else
            axPhasePos = get(gui.axes_PhasePlane_1,'UserData');
            if inpolygon(dlgCoords(1),dlgCoords(2),axPhasePos(:,1),axPhasePos(:,2))
                phaseAxes = gui.axes_PhasePlane_1;
            end
        end
    end
end

if isempty(thisAxes) && isempty(phaseAxes)
    delete(findobj('Tag','dynResSnapPointDispl'))
    delete(findobj('Tag','dynResSnapPointVeloc'))
    delete(findobj('Tag','dynResSnapPointAccel'))
    delete(findobj('Tag','dynResSnapPointPPlan'))
    return
elseif isempty(thisAxes) && ~isempty(phaseAxes)
    cam = getappdata(0,'phaseCam');
    if ~isempty(cam)
        setappdata(0,'phaseCam',[]);
        return
    end
    
    if ~strcmp(get(dlg, 'SelectionType'),'alt')
        setappdata(0,'phaseCam',[]);
        return
    end
    
    set(gui.popupmenu_Phase,'value',1);
    
    % Catches the inicial coordinates on window.
    set(dlg,'Units','pixels');
    crd = get(dlg, 'CurrentPoint');
    cx = crd(1);
    cy = crd(2);
    [az,el] = view(phaseAxes);
    set(dlg,'Units','normalized');
    
    % Incorporate the values on the variables.
    cam.rotIniAz = az;
    cam.rotIniEl = el;
    cam.rotIniX = cx;
    cam.rotIniY = cy;
    setappdata(0,'phaseCam',cam);
    return
end

delete(findobj('Tag','dynResSnapPointDispl'))
delete(findobj('Tag','dynResSnapPointVeloc'))
delete(findobj('Tag','dynResSnapPointAccel'))
delete(findobj('Tag','dynResSnapPointPPlan'))

% Check if frequency is being displayed
if get(gui.popupmenu_Domain,'Value') == 2
    return
end

% Get current cursor coordinates
coords = get(thisAxes,'CurrentPoint');
coords = coords(1,1:2);

% Get all plots on this canvas
name = get(thisAxes,'Tag');
plots = findobj('Tag',sprintf('Dynamic%s',name(6:end)));

% Check if there aren't any plots
if isempty(plots)
    return
end

% Check if coords are within solution interval
X = get(thisAxes,'XLim');
if coords(1) < X(1) || coords(1) > X(2)
    return
end

% Check if there are only two plotting points (null solution)
if length(plots(1).XData) <= 2
    return
end

% Get names of other axes
axNames = [ 'Displ' ; 'Veloc' ; 'Accel'];
if strcmp(axNames(1,:),name(6:end))
    axNames(1,:) = [];
elseif strcmp(axNames(2,:),name(6:end))
    axNames(2,:) = [];
elseif strcmp(axNames(3,:),name(6:end))
    axNames(3,:) = [];
else
    return
end

% Compute time step
step = plots(1).XData(2) - plots(1).XData(1);

% Get index of x coord on XData
i_1 = find(abs(plots(1).XData - coords(1)) <= step/2);

% Check if index is out of matrix
if isempty(i_1)
    return
end
i_1 = i_1(1);

% Get closer step point to current point (x)
dx = plots(1).XData(i_1);

% Get y coordinate
aux_y = zeros(length(plots),2);
for j = 1:length(plots)
    aux_y(j,:) = [plots(j).YData(i_1), j];
end
[~,indx] = min(abs(aux_y(:,1)-coords(2)));
y = aux_y(indx(1),1);
i = aux_y(indx(1),2);
pt = [dx, y];

% Draw visual feedback
% Selected axes
axes(thisAxes);

axYLim = get(thisAxes,'YLim');
plot([pt(1), pt(1)],axYLim,'Color',[0,0,0,0.5],'LineWidth',0.5,...
     'LineStyle','--','Tag',sprintf('dynResSnapPoint%s',name(6:end)));
 
scatter(pt(1), pt(2), 25,...
        'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],...
        'Tag',sprintf('dynResSnapPoint%s',name(6:end)));

% Not selected axes 1
axes(otherAxes(1));

plots_ax1 = findobj('Tag',sprintf('Dynamic%s',axNames(1,:)));
pt_ax1 = [dx plots_ax1(i).YData(i_1)];

axYLim = get(otherAxes(1),'YLim');
plot([pt_ax1(1), pt_ax1(1)],axYLim,'Color',[0,0,0,0.5],'LineWidth',0.5,...
    'LineStyle','--','Tag',sprintf('dynResSnapPoint%s',axNames(1,:)));

scatter(pt_ax1(1), pt_ax1(2), 25,...
        'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],...
        'Tag',sprintf('dynResSnapPoint%s',axNames(1,:)));

% Not selected axes 2
axes(otherAxes(2));

plots_ax2 = findobj('Tag',sprintf('Dynamic%s',axNames(2,:)));
pt_ax2 = [dx plots_ax2(i).YData(i_1)];

axYLim = get(otherAxes(2),'YLim');
plot([pt_ax2(1), pt_ax2(1)],axYLim,'Color',[0,0,0,0.5],'LineWidth',0.5,...
    'LineStyle','--','Tag',sprintf('dynResSnapPoint%s',axNames(2,:)));

scatter(pt_ax2(1), pt_ax2(2), 25,...
        'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],...
        'Tag',sprintf('dynResSnapPoint%s',axNames(2,:)));
    
% Draw visual feedback on phase portrait axes
displ = findobj('Tag','DynamicDispl');
veloc = findobj('Tag','DynamicVeloc');
accel = findobj('Tag','DynamicAccel');

axes(gui.axes_PhasePlane_1);
gui.axes_PhasePlane_1.Clipping = 'off';

axXLim = get(gui.axes_PhasePlane_1,'XLim');
axYLim = get(gui.axes_PhasePlane_1,'YLim');
axZLim = get(gui.axes_PhasePlane_1,'ZLim');

plot3(axXLim,[veloc(i).YData(i_1), veloc(i).YData(i_1)],[accel(i).YData(i_1), accel(i).YData(i_1)],'Color',[0,0,0,0.5],...
    'LineWidth',0.5,'LineStyle','--','Tag','dynResSnapPointPPlan');
plot3([displ(i).YData(i_1), displ(i).YData(i_1)],axYLim,[accel(i).YData(i_1), accel(i).YData(i_1)],'Color',[0,0,0,0.5],...
    'LineWidth',0.5,'LineStyle','--','Tag','dynResSnapPointPPlan');
plot3([displ(i).YData(i_1), displ(i).YData(i_1)],[veloc(i).YData(i_1), veloc(i).YData(i_1)],axZLim,'Color',[0,0,0,0.5],...
    'LineWidth',0.5,'LineStyle','--','Tag','dynResSnapPointPPlan');

scatter3(displ(i).YData(i_1), veloc(i).YData(i_1),accel(i).YData(i_1), 25,...
        'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 0],...
        'Tag','dynResSnapPointPPlan');
    
%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function popupmenu_DOF_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Nodes_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Result_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Domain_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Phase_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
