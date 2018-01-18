classdef MOS5TModel < MOSModel
% MOS5TMODEL
% A 3D model of a MOS5T transistor
% Model has 4 terminal nodes, D, G, S, and B.

    properties (Access = public)
        
    end

    methods (Access = public)
        function mod = MOS5TModel(varargin) 
            mod = mod@MOSModel(varargin{:});

            % Handling the NIL Manually
            mod.NIL                   = NIL({'d', 'g', 's', 'b'});
            mod.explicit_output_names = {'idb', 'igb', 'isb'};
            mod.internal_unk_names    = {};
            mod.other_IO_names        = {'vdb', 'vgb', 'vsb'};
            mod.IO_names              = {mod.other_IO_names{:}, mod.explicit_output_names{:}};

            % Set up the number of variables for the model
            mod.n_explicit_outs = 3;
            mod.n_other_ios     = 3;
            mod.n_internal_unks = 0;
            mod.n_implicit_eqns = 0;

            % Default Jacobian container
            mod.J_def = FQJAC(mod);
        end

        function [fqei_vals, J] = fqeiJ(mod, x_in, maybe_flag, ~, flag, ~)
        % Alternate calling syntax: fqeiJ(mod, vecX, vecY, vecU, elflag, elmodel)
            if (nargin == 3)
                flag = maybe_flag; 
            end
            VecX = x_in;
            VecY = [x_in(1); x_in(3)];
            [fqei_vals, J_def] = mod.base_model.fqeiJ(VecX, VecY, [], flag, mod.base_model);
            fqei_vals.qe = [fqei_vals.qi(1,1); fqei_vals.qe(2, 1); fqei_vals.qi(2,1)];
            fqei_vals.fe = [fqei_vals.fi(1,1); fqei_vals.fe(2, 1); fqei_vals.fi(2,1)];
            fqei_vals.fi = zeros(mod.n_implicit_eqns, 1);
            fqei_vals.qi = zeros(mod.n_implicit_eqns, 1);

            if (nargout > 1)
                % Allocating the raw derivative values in appropriate places and using DISTRIBUTION OF DERIVATIVE
                J = mod.J_def;
                J.Jqe.dqe_dvecX = [...
                    [J_def.Jqi.dqi_dvecY(1,1), J_def.Jqe.dqe_dvecX(1,2), ...
                        J_def.Jqi.dqi_dvecY(1,2)]; ...
                    [J_def.Jqe.dqe_dvecY(2,1), J_def.Jqe.dqe_dvecX(2,2), ...
                        J_def.Jqe.dqe_dvecX(2,3)]; ...
                    [J_def.Jqi.dqi_dvecY(2,1), J_def.Jqe.dqe_dvecX(3,2), ...
                        J_def.Jqi.dqi_dvecY(2,2)] ];
                J.Jfe.dfe_dvecX = [...
                    [J_def.Jfi.dfi_dvecY(1,1), J_def.Jfe.dfe_dvecX(1,2), ...
                        J_def.Jfi.dfi_dvecY(1,2)]; ...
                    [J_def.Jfe.dfe_dvecY(2,1), J_def.Jfe.dfe_dvecX(2,2), ...
                        J_def.Jfe.dfe_dvecX(2,3)]; ...
                    [J_def.Jfi.dfi_dvecY(2,1), J_def.Jfe.dfe_dvecX(3,2), ...
                        J_def.Jfi.dfi_dvecY(2,2)] ];
            end
        end
    % end methods
    end
    
% end classdef
end
