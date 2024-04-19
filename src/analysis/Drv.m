%% Drv (Driver) Class
%
%% Description
%
% This is a handle super-class for the definition of an analysis driver.
%
% An analysis driver object is responsible for running the functions that
% drives the structural analysis.
%
% This super-class defines abstract methods that must be implemented in
% the derived *sub-classes*:
%
% * <drv_led.html linear-elastic static analysis>.
% * <drv_les.html linear-elastic dynamic analysis>.
%
%% History
%
% Updated to version 2.0: August 2018 by PCLopes
%    When spring supports and load cases were added to LESM.
%    Internal displacements and stresses are now computed by using vector
%    and matrix operations, instead of a loop through all internal points
%    on elements.
%    Created properties neqspring and gblSprStiff.
%    Created properties nlc, ncomb, strLc, strComb, loadComb and loadCombID.
%    Created method elemIntStress.
%    Modified methods process and procces_GUI to call Anm method 
%    setupSpringStiff.
%    Modified methods assembleDOFNum and solveEqnSystem to account for 
%    spring degrees-of-freedom.
%
% Modified on 05-dec-2018 by LFMartha
%    Version 2.1: When created semi-rigid joint class.
%    Created properties srjoints and njoints.
%    Created methods computeNumSrjoints and createsSrjoints.
%    Modified method assembleGle to call abstract methods assembleGle and
%    assembleGlj of Anm object, and renamed it to assembleGatherVectors.
%    Created method assembleSrjointMtx.
%
% Modified on 21-feb-2019 by PCLopes
%    Changed properties strLc and strComb to cell arrays, to avoid bugs on
%    the GUI_EditLoadCase.
%    Modified method dimKFD so that the dimensioning of the stiffness
%    matrix, the forcing vector and the displacement vector are sparse.
%    Modified method solveEqnSystem to account for a sparse stiffness
%    matrix when using the rcond built-in function, to check if structure
%    is unstable. Also added a line to call an anm method to rotate
%    displacements on local inclined support axis to global axis.
%
% Modified on 29-apr-2019 by PCLopes    #####  ATTENTION  #####
%    MAJOR CHANGES due to program restructuring. All properties related to
%    model information were moved to a new class, called "Model". From now
%    on, this is an abstact superclass, idealized to handle solely the
%    process of different types of analysis.
%    Properties now consist on an identifier as to what type of analysis is
%    being considered, a flag for graphical or non-graphical, a handle to a
%    Model object and a handle to an Algorithm object.
%    In essence, the behavior of the analysis is contemplated by this
%    superclass, while all the model info migrated to the model class. This
%    is to achieve code modularization and organization while dealing with
%    multiple analysis types. For now, linear elastic static and linear
%    elastic dynamic are considered, but a path is being laid out for more
%    possibilities to come.
%
classdef Drv < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        analysis  =    0;     % flag for type of analysis
        graphical = true;     % flag for graphical or non-graphical
        model     =   [];     % handle to a model object
        solver    =   [];     % handle to an algorithm object
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function drv = Drv(analysis_flag,graphical_flag,model)
            drv.analysis  =  analysis_flag;
            drv.graphical = graphical_flag;
            drv.model     =          model;
            drv.setSolver;
        end
    end
    
    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % This method is called upon Drv object creation.
        % Creates a Solver object, according to analysis type, and
        % store a handle to it as a "solver" property of the Drv object.
        function setSolver(drv)
            include_constants;
            switch drv.analysis
                case STATIC_LINEAR
                    drv.solver = Solver_LES(drv);
                case DYNAMIC_NEWMARK_LINEAR
                    drv.solver = Solver_LED_Newmark(drv);
                case DYNAMIC_MODALSUP_LINEAR
                    drv.solver = Solver_LED_ModalSup(drv);
                case DYNAMIC_RK4_LINEAR
                    drv.solver = Solver_LED_RK4(drv);
                case DYNAMIC_AM3_LINEAR
                    drv.solver = Solver_LED_AM3(drv);
                case DYNAMIC_WILSON_LINEAR
                    drv.solver = Solver_LED_Wilson(drv);
            end
        end
    end
    
    %% Abstract methods
    methods (Abstract)
        %------------------------------------------------------------------
        % Dimensions and initializes global matrices and vectors needed for
        % the current type of analysis
        dimMtxVctrs(drv)
        
        %------------------------------------------------------------------
        % Assembles global matrices
        gblMtx(drv)
        
        %------------------------------------------------------------------
        % Computes element internal forces in several cross-section
        % positions.
        elemIntForces(drv)
        
        %------------------------------------------------------------------
        % Computes element internal displacements in several cross-section
        % positions.
        elemIntDispl(drv)
        
        %------------------------------------------------------------------
        % Processes current model data according to the current analysis
        % type
        status = process(drv)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of the driver object.
        function clean(drv)
            drv.analysis  =  0;
            drv.graphical = [];
            drv.model     = [];
            drv.solver    = [];
        end
    end
end