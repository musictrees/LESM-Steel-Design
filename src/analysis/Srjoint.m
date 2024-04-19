%% Srjoint (Semi-Rigid Joint) Class
%
%% Description
%
% This is a handle class for the definition of a semi rigid joint.
%
% This class generically handles a three-dimensional behavior of a semi-
% rigid joint. An object of the <anm.html *Anm*> class is responsible to
% "project" this generic 3D behavior to a specific model behavior, such
% as 2D frame, 2D truss, grillage, 3D truss or 3D frame model.
%
% Please refer to <elem.html *Elem*> class for definition of local
% coordinate system, since a semi-rigid joint uses the local coordinate
% system of its parents element.
%
classdef Srjoint < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        anm  = [];   % handle to an object of the Anm class
        elem = [];   % handle to parent element
        node = [];   % handles to node to which semi-rigid joint is connected
        krx  = [];   % rotational stiffness around local x axis
        kry  = [];   % rotational stiffness around local y axis
        krz  = [];   % rotational stiffness around local z axis
        rot  = [];   % joint rotation transformation matrix
        eqs  = [];   % equation numbers of additional joint rot. d.o.f.
        glj  = [];   % gather vector (stores element d.o.f. eqn. numbers)
        kjl  = [];   % stiffness matrix in local system
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function srj = Srjoint(anm,elem,node,krx,kry,krz)
            if (nargin > 0)
                srj.anm  = anm;
                srj.elem = elem;
                srj.node = node;
                srj.krx  = krx;
                srj.kry  = kry;
                srj.krz  = krz;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Sets semi-rigid joint object properties, accordingly to
        % respective element.
        % Input:
        %  elem: handle to object of the Elem class
        function setElemToSrj(srj,elem)
            srj.anm  = elem.anm;
            srj.elem = elem;
            srj.rot  = elem.anm.gblToLocSrjointRotMtx(srj);
        end
        
        %------------------------------------------------------------------
        % Computes joint stiffness matrix in global system.
        % Output:
        %  kjg: joint stiffness matrix in global system
        function kjg = gblStiffMtx(srj)
            % Compute and store joint stiffness matrix in local system
            srj.kjl = srj.anm.jointLocStiffMtx(srj);
            
            % Transform joint stiffness matrix from local to global system
            kjg = srj.rot' * srj.kjl * srj.rot;
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of an Srjoint object.
        function clean(srj)
            srj.anm  = [];
            srj.elem = [];
            srj.node = [];
            srj.krx  = [];
            srj.kry  = [];
            srj.krz  = [];
            srj.rot  = [];
            srj.glj  = [];
            srj.kjl  = [];
        end
    end
end