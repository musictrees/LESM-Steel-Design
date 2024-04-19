%% Global Constants
% This file contains variables that have a global scope and store flags to
% help understanding the program.
% It is necessary to include this file in every function that uses a global
% constant.
%
%% Types of analysis models
global TRUSS2D_ANALYSIS    % 2D truss analysis
global FRAME2D_ANALYSIS    % 2D frame analysis
global GRILLAGE_ANALYSIS   % grillage analysis
global TRUSS3D_ANALYSIS    % 3D truss analysis
global FRAME3D_ANALYSIS    % 3D frame analysis

TRUSS2D_ANALYSIS  = 0;
FRAME2D_ANALYSIS  = 1;
GRILLAGE_ANALYSIS = 2;
TRUSS3D_ANALYSIS  = 3;
FRAME3D_ANALYSIS  = 4;

%% Types of elements
global MEMBER_NAVIER       % Navier (Euler-Bernoulli) beam element
global MEMBER_TIMOSHENKO   % Timoshenko beam element

MEMBER_NAVIER     = 0;
MEMBER_TIMOSHENKO = 1;

%% Types of mass matrices
global LUMPED_MASS       % lumped mass matrices
global CONSISTENT_MASS   % displacement function consistent mass matrices
global MIXED_MASS        % half lumped mass, half consistent mass

LUMPED_MASS     = 0;
CONSISTENT_MASS = 1;
MIXED_MASS      = 2;

%% Types of continuity conditions
global HINGED_END       % hinged element end
global CONTINUOUS_END   % continuous element end
global SEMIRIGID_END    % semi rigid joint

HINGED_END     = 0;
CONTINUOUS_END = 1;
SEMIRIGID_END  = 2;

%% Types of load directions
global GLOBAL_LOAD   % element load applied in global system
global LOCAL_LOAD    % element load applied in local system

GLOBAL_LOAD = 0;
LOCAL_LOAD  = 1;

%% Types of essential boundary conditions
global FICTFIXED_DOF   % fictitious rotation constraints
global FREE_DOF        % free degree of freedom
global FIXED_DOF       % fixed degree of freedom
global SPRING_DOF      % degree of freedom partially constrained by spring

FICTFIXED_DOF = -1;
FREE_DOF      =  0;
FIXED_DOF     =  1;
SPRING_DOF    =  2;

%% Types of analysis
global LE_STATIC_ANALYSIS    % linear-elastic static analysis
global LE_DYNAMIC_ANALYSIS   % linear-elastic dynamic analysis

LE_STATIC_ANALYSIS  = 1;
LE_DYNAMIC_ANALYSIS = 2;

%% Types of dynamic analysis responses
global MODAL_ANALYSIS         % modal analysis only
global TRANS_MODAL_ANALYSIS   % transient + modal analysis

MODAL_ANALYSIS       = 0;
TRANS_MODAL_ANALYSIS = 1;

%% Types of solvers
global STATIC_LINEAR             % linear-elastic static solver
global DYNAMIC_NEWMARK_LINEAR    % linear-elastic dynamic numerical solver
global DYNAMIC_MODALSUP_LINEAR   % linear-elastic dynamic semi-analytical uncoupled solver
global DYNAMIC_RK4_LINEAR        % linear-elastic dynamic rk4 coupled solver
global DYNAMIC_AM3_LINEAR        % linear-elastic dynamic Adams-Moulton coupled solver
global DYNAMIC_WILSON_LINEAR     % linear-elastic dynamic Wilson Theta coupled solver

STATIC_LINEAR           =  0;
DYNAMIC_NEWMARK_LINEAR  =  1;
DYNAMIC_MODALSUP_LINEAR =  2;
DYNAMIC_RK4_LINEAR      =  3;
DYNAMIC_AM3_LINEAR      =  4;
DYNAMIC_WILSON_LINEAR   =  5;

%% Types of loads
global STATIC     % static load
global PERIODIC   % dynamic periodic load
global SLOPE      % dynamic slope load
global TABLE      % dynamic time table load

STATIC   =  0;
PERIODIC =  1;
SLOPE    =  2;
TABLE    =  3;

%% Types of damping
global XI_1ST_MODE       % critical damping ratio of 1st mode
global XI_1ST_2ND_MODE   % critical damping ratio of 1st and 2nd modes
global RAYLEIGH_COEFFS   % rayleigh coeffs provided

XI_1ST_MODE     =  0;
XI_1ST_2ND_MODE =  1;
RAYLEIGH_COEFFS =  2;
