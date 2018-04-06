classdef MOSModel < ModSpecModel
    properties (Access = protected)
        base_model = struct();
    end

    methods (Access = protected)    % Internal functions

    end

    methods (Abstract, Access = public) % Abstract methods

    end

    methods (Access = public)
        function mod = MOSModel(m_handle, model_identifier)
            base_model = setMOSFETParms(m_handle(model_identifier), model_identifier);
            mod = mod@ModSpecModel(base_model);
            mod.base_model = base_model;
        end

        function [fe, qe, fi, qi] = fqei(mod, x_in, flag)
            if (nargin < 3)
                flag = mod.flag;
            end
            fqei_vals = mod.fqeiJ(x_in, flag);

            % The results from the original model need to be resized
            qe = fqei_vals.qe;
            fe = fqei_vals.fe;
            fi = fqei_vals.fi;
            qi = fqei_vals.qi;

            if (nargout == 1)
                fe = cat(1, qe, qi, fe, fi);
            end
        end

        % Default implementation of fqeiJ
        function [fqei, J] = fqeiJ(mod, x_in, maybe_flag, ~, flag, ~)
            if (nargin == 3)
                flag = maybe_flag;
            end
            [fqei, J] = mod.base_model.fqeiJ(x_in, [], [], flag, ...
                mod.base_model);
        end
    end
end
