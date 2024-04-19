function varargout = GUI_About(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_About_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_About_OutputFcn, ...
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

function GUI_About_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;

% Move to center
movegui(hObject,'center');

% Set logos
axes(handles.logo_puc);
imshow('logo_puc_about.png');
axis image;

axes(handles.logo_tecgraf);
imshow('logo_tecgraf_about.png');
axis image;

guidata(hObject, handles);

function varargout = GUI_About_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;
