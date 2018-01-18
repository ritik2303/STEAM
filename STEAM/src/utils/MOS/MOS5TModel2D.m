classdef MOS5TModel2D < MOSModel2D
    methods (Access = public)
        function mod = MOS5TModel2D(varargin)
        % function mod = MOS5TModel2D(varargin)
        % Class constructor
            mod = mod@MOSModel2D(varargin{:});
        end

        function [fqei_vals, J] = fqeiJ(mod, x_in, maybe_flag, ~, flag, ~)
        % Alternate calling syntax: fqeiJ(mod, vecX, vecY, vecU, elflag, elmodel)
            if (nargin == 3)
                flag = maybe_flag; 
            end
            VecX = [x_in; 0];
            VecY = [x_in(1); 0];

            [fqei_vals, J_def] = mod.base_model.fqeiJ(VecX, VecY, [], flag, mod.base_model);
            
            % Allocating the raw function values into appropriate places
            fqei_vals.qe = [fqei_vals.qi(1,1) + fqei_vals.qe(1, 1); fqei_vals.qe(2,1)]; % Functional cancellation of
                                                                                        % resistor current/charge terms
            fqei_vals.fe = [fqei_vals.fi(1,1) + fqei_vals.fe(1, 1); fqei_vals.fe(2,1)];
            fqei_vals.fi = zeros(mod.n_implicit_eqns, 1);
            fqei_vals.qi = zeros(mod.n_implicit_eqns, 1);

            if (nargout > 1)
                % Allocating the raw derivative values in appropriate places and using DISTRIBUTION OF DERIVATIVE
                J = mod.J_def;
                J.Jqe.dqe_dvecX = [ [J_def.Jqi.dqi_dvecY(1,1) + J_def.Jqe.dqe_dvecY(1,1), J_def.Jqi.dqi_dvecX(1,2) + J_def.Jqe.dqe_dvecX(1,2)]; ...
                                    [J_def.Jqe.dqe_dvecY(2,1), J_def.Jqe.dqe_dvecX(2,2)] ];
                J.Jfe.dfe_dvecX = [ [J_def.Jfi.dfi_dvecY(1,1) + J_def.Jfe.dfe_dvecY(1,1), J_def.Jfi.dfi_dvecX(1,2) + J_def.Jfe.dfe_dvecX(1,2)]; ...
                                    [J_def.Jfe.dfe_dvecY(2,1), J_def.Jfe.dfe_dvecX(2,2)] ];
            end
        end

        function val = getIds(Obj, x_in)
            [~, ~, fi, ~] = Obj.fqei(x_in, flag);
            val = fi(1, 1); % TODO: This needs to be changed for vectorization
        end
    end
end
