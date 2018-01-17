classdef STEAM < ModSpecModel
    properties (Access = protected)
        m_interp = struct();
        i_method = ''; % Interpolation Method
        d_method = ''; % Discretization Method
        d_order  = 0;
        d_bounds = cell(0);
    end
    
    properties (SetAccess = protected, GetAccess = public)
        n_in_dims = 0;
    end

    methods (Access = protected)
        tabulateModSpec(oMOD, iMOD, i_method, d_method, d_order, d_bounds);
    end

    methods (Access = public)
        function Obj = STEAM(MOD, i_method, d_method, d_order, d_bounds)
        % Class constuctor
        % Takes in parameters relevant to the tabulation of a model and creates
        % a STEAM approximation for it.
            Obj           = Obj@ModSpecModel(MOD);
            Obj.n_in_dims = MOD.n_other_ios + MOD.n_internal_unks;
            Obj.tabulateModSpec(MOD, i_method, d_method, d_order, d_bounds);
        end

        function plotChebCoeffs(MOD, pc_index)
        % Plot the Chebyshev coefficients for the interpolant: Useful for
        % getting insights about the Chebyshev series coefficients.
            MOD.m_interp.plotChebCoeffs(pc_index);
        end

        % Functions described in auxiliary files.
        [fe, qe, fi, qi] = fqei(mod, x_in, flag)
        [fqei_vals, J]   = fqeiJ(mod, vecX, vecY, vecU, flag, ~);
    end
end

