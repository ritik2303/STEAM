classdef MOSModel2D < MOSModel
    properties (Access = protected)
        sampling_args = [];
    end

    methods (Access = public)
        function mod = MOSModel2D(m_handle, model_identifier, data_pts)
            mod = mod@MOSModel(m_handle, model_identifier);
            if (nargin > 2)
                mod.initialize(data_pts);
            end

            % Handling the NIL Manually
            mod.NIL                   = NIL({'d', 'g', 's'});
            mod.explicit_output_names = {'ids', 'igs'};
            mod.internal_unk_names    = {};
            mod.other_IO_names        = {'vds', 'vgs'};
            mod.IO_names              = {mod.other_IO_names{:}, mod.explicit_output_names{:}};

            % Set up the number of variables
            mod.n_explicit_outs = 2;
            mod.n_other_ios     = 2;
            mod.n_internal_unks = 0;
            mod.n_implicit_eqns = 0;

            % Set up the Jacobian container
            mod.J_def = FQJAC(mod);
        end

        function n_dims = getModelDims(Obj)
            n_dims = Obj.sampling_args.getNDims();
        end
    end
end
