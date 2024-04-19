classdef Steel < handle & Material
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        A36 ={"Steel-ASTM A36",200*1e9, 0.32, 75*1e9, 12*1e-6, 7850, 250*1e6, 400*1e6; }; %Pa,ad,Pa,nao sei,kg/m^3,Pa,Pa
         A57 ={"Steel-ASTM A570",205*1e9, 0.32, 75*1e9, 12*1e-6, 7850, 345*1e6, 400*1e6; };
        
        steelList={Steel.A36,Steel.A57};
    end
    
    methods
        function steel=Steel(id,e,v,te,rho,leak,yld)
            if (nargin > 0)
                steel.id = id;
                steel.elasticity = e;
                steel.poisson = v;
                steel.shear = e / (2 * (1 + v));
                steel.thermExp = te;
                steel.density = rho;
                steel.leakage=leak;
                steel.yield=yld;
               
            end
        end
    end
end

