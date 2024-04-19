%% Edit Load Case Dialog Callback Functions
% This file contains the callback functions associated with the edit load
% cases dialog of the graphical version of the LESM program.
% Common input arguments for all callback functions:
%  hObject: handle to interface object related to the function
%  eventdata: reserved - to be defined in a future version of MATLAB
%  handles: structure with handles and user data
%
%% ------------------------------------------------------------------------
% GUI initialization function.
function varargout = GUI_EditLoadCase(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_EditLoadCase_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_EditLoadCase_OutputFcn, ...
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

%--------------------------------------------------------------------------
% Executes just before GUI_EditLoadCase is made visible.
% Sets GUI initial properties.
function GUI_EditLoadCase_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% Choose default command line output for GUI_EditLoadCase
handles.output = hObject;

% Move GUI to the center of the screen
if getappdata(0,'move') == 1
    movegui(hObject,'center');
end

% Make GUI a modal window
set(hObject,'WindowStyle','modal');

% Get model object
model = getappdata(0,'model');
  
% Create a list of load cases
lc = model.strLc;
if ~isempty(lc)
    set(handles.listbox_LoadCase,'string',lc,'Value',1,'Max',size(lc,1))
end    

% Enable/Disable delete load case button
nlc = model.nlc;
if nlc >= 2
    set(handles.pushbutton_Delete,'Enable','on')
else
    set(handles.pushbutton_Delete,'Enable','off')
end    

% Create a list of combinations
combs = model.strComb;
if ~isempty(combs)
    set(handles.listbox_Combinations,'string',combs,'Value',1,'Max',size(combs,1))
end    

% Enable/Disable delete combinations pushbutton, listbox and uitable
ncomb = model.ncomb;
if ncomb >= 1
    set(handles.pushbutton_DeleteComb,'Enable','on')
    set(handles.listbox_Combinations,'Enable','on')
    set(handles.uitable_Combinations,'Enable','on')
else
    set(handles.pushbutton_DeleteComb,'Enable','off')
    set(handles.listbox_Combinations,'Enable','off')
    set(handles.uitable_Combinations,'Enable','off')
end

% Set callback for cell edition in uitable_Combinations
set(handles.uitable_Combinations,'CellEditCallback',@uitable_Combinations_CellEditCallback)

% Write first combination's info on uitable
if ncomb >= 1
    loadComb = model.loadComb(:,1);
    set(handles.listbox_Combinations,'Value',1)
    count = 0;
    for i = 1:size(loadComb,1)
        if model.loadCombID(i,1) ~= 0
            count = count + 1;
            nextData = {char(lc(i)),loadComb(i)};
            if count == 1
                set(handles.uitable_Combinations,'Data',nextData);
            else
                tableData = get(handles.uitable_Combinations,'Data');
                tableData(size(tableData,1)+1,:) = nextData;
                set(handles.uitable_Combinations,'Data',tableData);
            end
        end
    end  
else
    set(handles.uitable_Combinations,'Enable','off','Data',{})
end

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = GUI_EditLoadCase_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_New.
function pushbutton_New_Callback(~, ~, handles) %#ok<DEFNU>
model = getappdata(0,'model');
nlc = model.nlc;
lc = model.strLc;
newLoadCase = char(inputdlg('Enter new load case name','New Load Case',[1 50]));

if ~isempty(newLoadCase)
    % Check if chosen name is already in use, or is empty
    if strcmp(newLoadCase,get(handles.listbox_LoadCase,'String')) == 1
        msgbox('There is already a load case with this name, please choose another name.', 'Error','error');
        return
    elseif strcmp(newLoadCase,' ')
        return
    elseif strcmp(newLoadCase,'  ')
        return
    elseif strcmp(newLoadCase,'   ')
        return
    elseif strcmp(newLoadCase,'    ')
        return
    elseif strcmp(newLoadCase,'     ')
        return
    elseif strcmp(newLoadCase,'      ')
        return
    end
    
    for i = 1:size(lc,1)
        if strcmp(newLoadCase,lc(i)) == 1
            msgbox('There is already a load case with this name, please choose another name.', 'Error','error');
            return  
        end    
    end
    
    % Write new load case name on listbox and update model object
    lc = char(lc);
    set(handles.listbox_LoadCase,'string',{lc,newLoadCase},'Max',get(handles.listbox_LoadCase,'Max')+1)
    set(handles.pushbutton_Delete,'Enable','on')
    model.nlc = nlc + 1;
    model.strLc = get(handles.listbox_LoadCase,'string');
    
    % Write new load case in the popupmenu (GUI_Main) and update model.loadComb
    ncomb = model.ncomb;
    comb = model.strComb;
    mdata = guidata(findobj('Tag','GUI_Main'));
    currentLc = get(mdata.popupmenu_LoadCase,'Value');
    if ncomb == 0
        set(mdata.popupmenu_LoadCase,'String',get(handles.listbox_LoadCase,'String'),...
            'Max',get(handles.listbox_LoadCase,'Max'));
    else
        if currentLc <= nlc
            set(mdata.popupmenu_LoadCase,'String',{char(get(handles.listbox_LoadCase,'String')),char(comb)},...
               'Max',get(handles.listbox_LoadCase,'Max')+ncomb);
        else
            set(mdata.popupmenu_LoadCase,'String',{char(get(handles.listbox_LoadCase,'String')),char(comb)},...
               'Max',get(handles.listbox_LoadCase,'Max')+ncomb,'Value',currentLc+1);
            setappdata(0,'currentLc',currentLc+1);
        end
        model.loadComb(size(model.loadComb,1)+1,:) = zeros(1,ncomb);
        model.loadCombID(size(model.loadCombID,1)+1,:) = zeros(1,ncomb);
    end
    
    % Save model in root
    setappdata(0,'model',model)
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Delete.
function pushbutton_Delete_Callback(hObject, ~, handles) %#ok<DEFNU>
mdata = guidata(findobj('Tag','GUI_Main'));
nodes = getappdata(0,'nodes');
elems = getappdata(0,'elems');
model = getappdata(0,'model');
nlc = model.nlc;
ncomb = model.ncomb;
lc = model.strLc;
combs = model.strComb;

% Load cases to be deleted
delLc = get(handles.listbox_LoadCase,'Value');

% Check if all load cases are selected
if size(delLc,2) == size(lc,1)
    msgbox('Cannot delete all load cases.', 'Error','error');
    return
end

% Initialize redraw flag (info if the model needs to be redrawn in canvas)
redrawFlag = 0;

% Delete selected load cases from load case names vector
lc(delLc) = [];

% Set new value and max propeties in listbox_LoadCase
% (prevents MATLAB warnings)
if min(delLc) ~= 1
    value = min(delLc) - 1;
else
    value = min(delLc);
end
maxLc = get(handles.listbox_LoadCase,'Max');
set(handles.listbox_LoadCase,'string',lc,'Value',value,'Max',maxLc-size(delLc,2))

% Get current popupmenu value
currentVal = get(mdata.popupmenu_LoadCase,'Value');

% Delete selected load cases from popupmenu in GUI_Main
if ncomb == 0
    if all(delLc ~= currentVal) == 1 % none of the deleted load cases are the current one
        auxValue = 0;
        for i = 1:size(delLc,2)
            if delLc(i) < currentVal
                auxValue = auxValue + 1;
            end
        end
        set(mdata.popupmenu_LoadCase,'String',get(handles.listbox_LoadCase,'String'),'Max',maxLc-size(delLc,2),'Value',currentVal - auxValue);
    else % the current load case is among those being deleted
        set(mdata.popupmenu_LoadCase,'String',get(handles.listbox_LoadCase,'String'),'Max',maxLc-size(delLc,2),'Value',value);
        redrawFlag = 1;
    end
else % there are combinations
    if all(delLc ~= currentVal) == 1 && currentVal <= nlc  % none of the deleted load cases are the current one
        auxValue = 0;
        for i = 1:size(delLc,2)
            if delLc(i) < currentVal
                auxValue = auxValue + 1;
            end
        end
        set(mdata.popupmenu_LoadCase,'String',{char(get(handles.listbox_LoadCase,'String')),char(combs)},...
            'Max',maxLc-size(delLc,2)+ncomb,'Value',currentVal - auxValue);
    elseif all(delLc ~= currentVal) == 1 && currentVal > nlc % current load case is a combination
        set(mdata.popupmenu_LoadCase,'String',{char(get(handles.listbox_LoadCase,'String')),char(combs)},...
            'Max',maxLc-size(delLc,2)+ncomb,'Value',currentVal - size(delLc,2));
    elseif all(delLc ~= currentVal) == 0 % the current load case is among those being deleted
        set(mdata.popupmenu_LoadCase,'String',{char(get(handles.listbox_LoadCase,'String')),char(combs)},...
            'Max',maxLc-size(delLc,2)+ncomb,'Value',value);
        redrawFlag = 1;
    end    
end    

% Check if there is only one load case remaining. If so, disable delete
% load case pushbutton.
if get(handles.listbox_LoadCase,'Max') == 1
    set(hObject,'Enable','off')
end    

% Delete selected load cases from load case matrix in each node and elem
auxDelLc = fliplr(delLc);
for dlc = 1:size(delLc,2)  
    for n = 1:model.nnp
        if auxDelLc(dlc) <= size(nodes(n).nodalLoadCase,2)
            nodes(n).nodalLoadCase(:,auxDelLc(dlc)) = [];
        end
    end

    for e = 1:model.nel
        if auxDelLc(dlc) <= size(elems(e).load.elemLoadCase,2)
            elems(e).load.elemLoadCase(:,auxDelLc(dlc)) = [];
        end
    end   
end

% Delete selected load cases from loadComb matrix
if ncomb ~= 0
    model.loadComb(delLc,:) = [];
    model.loadCombID(delLc,:) = [];
end

% Delete load combinations that no longer have any load cases
nc = ncomb; % initialize auxiliary variables
delLComb = 0;

for i = 1:nc
    if all(model.loadCombID(:,i) == 0) == 1 % if all terms in a column of model.loadCombID are zero
        ncomb = ncomb - 1;                  % all load cases of the respective combination have
                                            % been deleted.
        
        if get(handles.listbox_Combinations,'Value') == i % the now empty load comb is
            if i ~= 1                                     % the selected one.
                value = i -1;
            else
                value = i;
            end
        end
        
        % Update combinations listbox
        if ncomb ~= 0 % not all load combinations are empty
            strComb = get(handles.listbox_Combinations,'String');
            strComb(i,:) = []; % delete empty load comb name from string
            set(handles.listbox_Combinations,'String',char(strComb),'Max',ncomb,'Value',value)
            listbox_Combinations_Callback([], [], handles);
        else % all load combinations became empty (were deleted)
            set(handles.listbox_Combinations,'String',' ','Enable','off','Value',1,'Max',1)
            set(handles.pushbutton_DeleteComb,'Enable','off')
            data = {}; % clear uitable
            set(handles.uitable_Combinations,'Enable','off','Data',data)
        end
        
        % Update load cases popupmenu (GUI_Main)
        strPopupmenu = get(mdata.popupmenu_LoadCase,'String');
        strPopupmenu(nlc - size(delLc,2) + i - delLComb,:) = []; % delete empty load comb name from string
        set(mdata.popupmenu_LoadCase,'String',char(strPopupmenu),'Max',nlc - size(delLc,2) + ncomb)
        
        % Check if deleted combination is the current load case
        if get(mdata.popupmenu_LoadCase,'Value') == i + nlc
            set(mdata.popupmenu_LoadCase,'Value',1)
            redrawFlag = 1;
        end
        
        % Update counter
        delLComb = delLComb + 1;
    end
end

% Update uitable_Combinations
nlc = nlc - size(delLc,2);
strLc = get(handles.listbox_LoadCase,'string');
if ncomb ~= 0
    comb = get(handles.listbox_Combinations,'value');
    tableData = cell(nlc,2);
    count = 0;
    for i = 1:nlc
        if model.loadCombID(i,comb) ~= 0
            count = count + 1;
            tableData(count,:) = {char(strLc(i,:)),model.loadComb(i,comb)};
        else
            tableData(end,:) = [];
        end
    end
    if ~isempty(tableData)
        set(handles.uitable_Combinations,'Data',tableData);
    else
        set(handles.uitable_Combinations,'Data',{});
    end
end

% Delete empty load combs from model.loadComb
indexDelComb = zeros(delLComb,1);
count = 0;
for j = 1:size(model.loadComb,2)
    if all(model.loadCombID(:,j) == 0) == 1
        count = count + 1;
        indexDelComb(count) = j;
    end
end
model.loadComb(:,indexDelComb) = [];
model.loadCombID(:,indexDelComb) = [];

% Check if model needs to be redrawn (current load case was deleted)
if redrawFlag == 1    
    % Make GUI a normal window to avoid warning sound
    gui = findobj('Tag','GUI_EditLoadCase');
    set(gui,'WindowStyle','normal');
    
    setappdata(0,'currentLc',0)
    GUI_Main('popupmenu_LoadCase_Callback',mdata.popupmenu_LoadCase, [], mdata)
    
    % Make GUI a modal window again
    set(gui,'WindowStyle','modal');
end

% Update model object
model.nlc = nlc;
model.strLc = strLc;
model.ncomb = ncomb;
model.strComb = get(handles.listbox_Combinations,'String');
model.nodes = nodes;
model.elems = elems;

% Save variables in root
setappdata(0','currentLc',get(mdata.popupmenu_LoadCase,'value'))
setappdata(0,'nodes',nodes)
setappdata(0,'elems',elems)
setappdata(0,'model',model)

%--------------------------------------------------------------------------
% --- Executes on selection change in listbox_LoadCase.
function listbox_LoadCase_Callback(~, ~, ~) %#ok<DEFNU>

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function listbox_LoadCase_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Combination.
function pushbutton_Combination_Callback(~, ~, handles) %#ok<DEFNU>
model = getappdata(0,'model');
ncomb = model.ncomb;
combs = model.strComb;
newComb = char(inputdlg('Enter new combination name','New Combination',[1 50]));

% Add combination
if ~isempty(newComb)
    % Check if chosen name is already in use, or is empty
    if strcmp(newComb,get(handles.listbox_Combinations,'String')) == 1
        msgbox('There is already a load case with this name, please choose another name.', 'Error','error');
        return
    elseif strcmp(newComb,' ')
        return
    elseif strcmp(newComb,'  ')
        return
    elseif strcmp(newComb,'   ')
        return
    elseif strcmp(newComb,'    ')
        return
    elseif strcmp(newComb,'     ')
        return
    elseif strcmp(newComb,'      ')
        return
    end
    
    for i = 1:size(combs,1)
        if strcmp(newComb,combs(i,:)) == 1
            msgbox('There is already a combination with this name, please choose another name.', 'Error','error');
            return  
        end    
    end
    
    % Write name on combinations listbox
    if ncomb == 0   % there were previously no combinations
        set(handles.listbox_Combinations,'Enable','on','String',newComb,'Max',1)
    else   % there was already at least one combination
        combs = char(combs);
        set(handles.listbox_Combinations,'Enable','on','String',{combs,newComb},'Max',get(handles.listbox_Combinations,'Max')+1)
    end
    ncomb = ncomb + 1;
    set(handles.listbox_Combinations,'Value',ncomb)
    
    % Enable delete pushbutton and uitable_Combinations
    set(handles.pushbutton_DeleteComb,'Enable','on')
    set(handles.uitable_Combinations,'Enable','on')
    
    % Write new comb info on uitable
    lc = (get(handles.listbox_LoadCase,'Value'))'; % selected load cases that form combination
    strLc = get(handles.listbox_LoadCase,'String'); 
    loadCases = strLc(lc,:);  % names of selected load cases
    combFactors = ones(size(loadCases,1),1);  % default comb factor = 1
    
    for i = 1:size(lc,1)
        if i == 1 % Set first line of uitable
            set(handles.uitable_Combinations,'Data',{char(loadCases(i,:)),combFactors(i)}) 
        else % set the other lines
            data = get(handles.uitable_Combinations,'Data');
            newData = {char(loadCases(i,:)),combFactors(i)};
            data(size(data,1)+1,:) = newData;
            set(handles.uitable_Combinations,'Data',data) 
        end    
    end 
    
    % Update model.loadComb matrix
    nlc = model.nlc;
    loadComb = zeros(nlc,1);
    loadComb(lc) = ones(size(lc,1),1);
    if isempty(model.loadComb)
        model.loadComb = zeros(nlc,ncomb); 
    end
    model.loadComb(:,ncomb) = loadComb;
    
    % Update model.loadCombID matrix
    if isempty(model.loadCombID)
        model.loadCombID = loadComb;
    else
        model.loadCombID(:,ncomb) = loadComb;
    end
    
    % Update model object
    model.ncomb = ncomb;
    model.strComb = get(handles.listbox_Combinations,'String');
    setappdata(0,'model',model)
    
    % Writes new comb name in the popupmenu (GUI_Main)
    mdata = guidata(findobj('Tag','GUI_Main'));
    str = char(get(mdata.popupmenu_LoadCase,'String'));
    set(mdata.popupmenu_LoadCase,'String',{str,newComb},'Max',get(handles.listbox_Combinations,'Max')+nlc)  
end

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_DeleteComb.
function pushbutton_DeleteComb_Callback(hObject, ~, handles) %#ok<DEFNU>
mdata = guidata(findobj('Tag','GUI_Main'));
model = getappdata(0,'model');
nlc = model.nlc;
lc = model.strLc;
ncomb = model.ncomb;
combs = model.strComb;

% Initialize redraw flag (info if the model needs to be redrawn in canvas)
redrawFlag = 0;

% Load case combinations to be deleted
delLcomb = get(handles.listbox_Combinations,'Value');

% Delete selected combinations from combinations' names vector
combs(delLcomb,:) = [];

% Set new value property in listbox_Combinations (prevents MATLAB warnings)
if min(delLcomb) ~= 1
    value = min(delLcomb) - 1;
else
    value = min(delLcomb);
end

maxLcomb = get(handles.listbox_Combinations,'Max');
ncomb = ncomb - size(delLcomb,2);

% Check if all combinations were deleted
if ncomb == 0
    set(handles.listbox_Combinations,'Enable','off','string',' ','Value',1,'Max',1)
    set(hObject,'Enable','off') % disable delete combinations button
    model.loadComb = []; % clean loadComb matrix
    model.loadCombID = []; % clean loadCombID matrix
    set(handles.uitable_Combinations,'Enable','off','Data',{}) % clear uitable
    if get(mdata.popupmenu_LoadCase,'Value') > nlc % current load case is a combination
        set(mdata.popupmenu_LoadCase,'Value',1);
        redrawFlag = 1;
    end
    set(mdata.popupmenu_LoadCase,'String',lc,'Max',nlc);
else  % there are combinations left
    set(handles.listbox_Combinations,'string',char(combs),'Value',value,'Max',maxLcomb-size(delLcomb,2))
    model.loadComb(:,delLcomb) = []; % delete load combinations from loadComb matrix
    model.loadCombID(:,delLcomb) = []; % delete load combinations from loadCombID matrix
    
    % Write selected combination info on the uitable
    count = 0;
    loadComb = model.loadComb(:,value);
    for i = 1:size(loadComb,1)
        if model.loadCombID(i,value) ~= 0
            count = count + 1;
            newData = {char(lc(i)),loadComb(i)};
            if count == 1  % first line of the uitable
                set(handles.uitable_Combinations,'Data',newData);
            else
                data = get(handles.uitable_Combinations,'Data');
                data(size(data,1)+1,:) = newData;
                set(handles.uitable_Combinations,'Data',data);
            end
        end
    end
    set(mdata.popupmenu_LoadCase,'String',{char(lc),char(combs)},'Max',maxLcomb-size(delLcomb,2)+nlc);
    if get(mdata.popupmenu_LoadCase,'Value') > nlc % current load case is a combination
        set(mdata.popupmenu_LoadCase,'Value',value + nlc);
        redrawFlag = 1;
    end
end

% Draw and update model
if redrawFlag == 1
    % Make GUI a normal window to avoid warning sound
    gui = findobj('Tag','GUI_EditLoadCase');
    set(gui,'WindowStyle','normal');
    
    setappdata(0,'currentLc',0);
    GUI_Main('popupmenu_LoadCase_Callback',mdata.popupmenu_LoadCase,[],mdata);
    
    % Make GUI a modal window again
    set(gui,'WindowStyle','modal');
end

% Update model object
model.ncomb = ncomb;
model.strComb = get(handles.listbox_Combinations,'String');

% Save variables in root
setappdata(0','currentLc',get(mdata.popupmenu_LoadCase,'value'));
setappdata(0,'model',model);

%--------------------------------------------------------------------------
% --- Executes on selection change in listbox_Combinations.
function listbox_Combinations_Callback(~, ~, handles)
vComb = get(handles.listbox_Combinations,'Value');
if size(vComb,2) >= 2 % if multiple combs, only consider the first one.
    vComb = vComb(1);
end
model = getappdata(0,'model');
lc = model.strLc;
loadComb = model.loadComb(:,vComb);

% Write selected combination info on the uitable
count = 0;
for i = 1:size(loadComb,1)
    if model.loadCombID(i,vComb) ~= 0
        count = count + 1;
        newData = {char(lc(i)),loadComb(i)};
        if count == 1
            set(handles.uitable_Combinations,'Data',newData);
        else
            data = get(handles.uitable_Combinations,'Data');
            data(size(data,1)+1,:) = newData;
            set(handles.uitable_Combinations,'Data',data);
        end
    end
end

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function listbox_Combinations_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% Executes when cell is edited in uitable_Combinations
function uitable_Combinations_CellEditCallback(~,~,~)
mdata = guidata(findobj('tag','GUI_EditLoadCase'));
set(mdata.pushbutton_Apply,'Enable','on')

%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_Apply.
function pushbutton_Apply_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% Get model object
model = getappdata(0,'model');

% Get uitable data and selected load comb
tableData = get(handles.uitable_Combinations,'Data');
whichComb = get(handles.listbox_Combinations,'value');

% Set uitable data (user input) to model (COMB FACTORS)
count = 0;
combFactorsChanged = false;
for i = 1:size(model.loadComb,1)
    if model.loadCombID(i,whichComb) ~= 0
        count = count + 1;
        if ~isnan(cell2mat(tableData(count,2)))
            if cell2mat(tableData(count,2)) ~= model.loadComb(i,whichComb)
                model.loadComb(i,whichComb) = cell2mat(tableData(count,2));
                combFactorsChanged = true;
            end
        end
    end
end

% Save model in root
setappdata(0,'model',model)

% Update uitable
listbox_Combinations_Callback([], [], handles)

% Check if combination that changed is current load case
if combFactorsChanged == true
    mdata = guidata(findobj('tag','GUI_Main'));
    currentLc = get(mdata.popupmenu_LoadCase,'value');
    if currentLc == (whichComb + model.nlc)
        setappdata(0,'currentLc',0)
        GUI_Main('popupmenu_LoadCase_Callback',mdata.popupmenu_LoadCase,eventdata,mdata)
    end
end

% Disable apply pushbutton
set(hObject,'enable','off')

%--------------------------------------------------------------------------
% --- Executes when user attempts to close GUI_EditLoadCase.
function GUI_EditLoadCase_CloseRequestFcn(hObject, ~, ~) %#ok<DEFNU>
delete(hObject);

%--------------------------------------------------------------------------
% --- Executes on key press with focus on GUI_EditLoadCase and none of its controls.
function GUI_EditLoadCase_KeyPressFcn(~,~,~) %#ok<DEFNU>
