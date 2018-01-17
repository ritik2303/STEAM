classdef MODDIMS
    properties ( Access = protected )
        n_other_ios = 0;
        n_internal_unks = 0;
        n_explicit_outs = 0;
        n_implicit_eqns = 0;
    end
    
    methods ( Access = public )
        function fqei = MODDIMS(MOD)   % Pass in a ModSpec model to gather data
            % Setting up values of sizing parameters
            fqei.n_other_ios = MOD.n_other_ios;
            fqei.n_internal_unks = MOD.n_internal_unks;
            fqei.n_explicit_outs = MOD.n_explicit_outs;
            fqei.n_implicit_eqns = MOD.n_implicit_eqns;
        end
    end
end
