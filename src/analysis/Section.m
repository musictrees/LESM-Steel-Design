%% Cross-Section Class
%
%% Description
%
% This is a handle class for the definition of an element cross-section.
%
% All cross-sections in LESM are considered to be of a generic type,
% which means that their shapes are not specified, only their geometric
% properties are provided, such as area, moment of inertia and height.
%
classdef Section < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id        = 0;   % identification number
        area_x    = 0;   % area relative to local x-axis (full area)
        area_y    = 0;   % area relative to local y-axis (effective shear area)
        area_z    = 0;   % area relative to local z-axis (effective shear area)
        inertia_x = 0;   % moment of inertia relative to local x-axis (torsion inertia)
        inertia_y = 0;   % moment of inertia relative to local y-axis (bending inertia)
        inertia_z = 0;   % moment of inertia relative to local z-axis (bending inertia)
        height_y  = 0;   % height relative to local y-axis
        height_z  = 0;   % height relative to local z-axis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function section = Section(id,Ax,Ay,Az,Ix,Iy,Iz,Hy,Hz)
            if (nargin > 0)
                section.id        = id;
                section.area_x    = Ax;
                section.area_y    = Ay;
                section.area_z    = Az;
                section.inertia_x = Ix;
                section.inertia_y = Iy;
                section.inertia_z = Iz;
                section.height_y  = Hy;
                section.height_z  = Hz;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a Section object.
        function clean(section)
            section.id        = 0;
            section.area_x    = 0;
            section.area_y    = 0;
            section.area_z    = 0;
            section.inertia_x = 0;
            section.inertia_y = 0;
            section.inertia_z = 0;
            section.height_y  = 0;
            section.height_z  = 0;
        end
    end
end