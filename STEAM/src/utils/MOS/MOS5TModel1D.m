classdef MOS5TModel1D < MOS5TModel2D
% Class MOS5TModel1D describes a diode-connected MOS5T Transistor. This diode
% has just 2 ports. The gate and the drain are connected, and similarly, the
% bulk and source are connected.  This model is derived from MOSModel2D, so the
% internal unknowns have been removed.
    methods (Access = public)
        function mod = MOS5TModel1D(m_handle, model_identifier, data_pts)
            mod = mod@MOS5TModel2D(m_handle, model_identifier);

            % Handling NIL Manually
            mod.NIL                    = NIL({'d', 's'});
            mode.explicit_output_names = {'ids'};
            mod.other_IO_names         = {'vds'};
            mod.IO_names               = {mod.other_IO_names{:}, mod.explicit_output_names{:}};

            % Setting up the number of unknowns for the model
            mod.n_explicit_outs = 1;
            mod.n_other_ios     = 1;
            mod.J_def           = FQJAC(mod);
        end

        function [fqei_vals, J] = fqeiJ(mod, x_in, maybe_flag, ~, flag, ~)
        % Alternate calling syntax: fqeiJ(mod, vecX, vecY, vecU, elflag, elmodel)
            if (nargin == 3)
                flag = maybe_flag;
            end
            [fqei_vals, J_def] = fqeiJ@MOS5TModel2D(mod, [x_in; x_in], flag);
            
            % Allocating the raw function values into appropriate places
            fqei_vals.qe = [fqei_vals.qe(1,1) + fqei_vals.qe(2,1)];
            fqei_vals.fe = [fqei_vals.fe(1,1) + fqei_vals.fe(2,1)];
            fqei_vals.fi = zeros(mod.n_implicit_eqns, 1);
            fqei_vals.qi = zeros(mod.n_implicit_eqns, 1);

            % Allocating the raw derivatives in appropriate places and using the chain rule
            % Let the common voltage be called Vc. We want d(I_total)/dVc.
            % We have d(Ids)/dvds, d(Igs)/dvds, d(Ids)/dvgs  and d(Igs)/dVgs, and therefore, 
            % d(I_total)/dvds = d(Ids)/dvds + d(Igs)/dvds
            % d(I_total)/dvgs = d(Ids)/dvgs + d(Igs)/dvgs
            % d(I_total)/dVc can be obtained with chain rule and the partial derivatives.
            % In this case, as Vc = vds = vgs, the above quantity is the sum of all the 4
            J = mod.J_def;
            J.Jfe.dfe_dvecX = sum(J_def.Jfe.dfe_dvecX(:));
            J.Jqe.dqe_dvecX = sum(J_def.Jqe.dqe_dvecX(:));
        end
    end
end
