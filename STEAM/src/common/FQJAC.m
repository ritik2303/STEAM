classdef FQJAC < MODDIMS
    properties ( Access = public )
        Jfe = struct();
        Jfi = struct();
        Jqe = struct();
        Jqi = struct();
    end

    methods ( Access = public )
        function J = FQJAC(MOD)
            J = J@MODDIMS(MOD);

            J.Jfe.dfe_dvecX = sparse(J.n_explicit_outs, J.n_other_ios);
            J.Jfe.dfe_dvecY = sparse(J.n_explicit_outs, J.n_internal_unks);
            J.Jfe.dfe_dvecU = sparse(J.n_explicit_outs, 0);
            J.Jfe.dfe_dvecLim = sparse(J.n_explicit_outs, 0);

            J.Jfi.dfi_dvecX = sparse(J.n_implicit_eqns, J.n_other_ios);
            J.Jfi.dfi_dvecY = sparse(J.n_implicit_eqns, J.n_internal_unks);
            J.Jfi.dfi_dvecU = sparse(J.n_implicit_eqns, 0);
            J.Jfi.dfi_dvecLim = sparse(J.n_implicit_eqns, 0);

            J.Jqe.dqe_dvecX = sparse(J.n_explicit_outs, J.n_other_ios);
            J.Jqe.dqe_dvecY = sparse(J.n_explicit_outs, J.n_internal_unks);
            J.Jqe.dqe_dvecU = sparse(J.n_explicit_outs, 0);
            J.Jqe.dqe_dvecLim = sparse(J.n_explicit_outs, 0);

            J.Jqi.dqi_dvecX = sparse(J.n_implicit_eqns, J.n_other_ios);
            J.Jqi.dqi_dvecY = sparse(J.n_implicit_eqns, J.n_internal_unks);
            J.Jqi.dqi_dvecU = sparse(J.n_implicit_eqns, 0);
            J.Jqi.dqi_dvecLim = sparse(J.n_implicit_eqns, 0);
        end
    end
end
