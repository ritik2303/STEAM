function [fqei_vals, J] = fqeiJ(MOD, vecX, vecY, vecU, flag, ~)
    if  (nargin > 3)
        x_in = [vecX; vecY; vecU];
    else
        x_in = vecX;
        flag = vecY;
    end
    % Recall that in file MOSModel2D.m, the calculated values are stored in the
    % order: [qe; qi; fe; fi]

    % First we calculate everything using the interpolant
    [all_vals, all_ders] = MOD.m_interp.computeWithDer( x_in );

    % Note that creating new structs and putting data into them might take up a
    % lot of time (Just saying that there might be faster ways of doing whatever
    % we are trying to do here).

    % Assigning values to fqei
    fqei_vals.qe = all_vals(1:MOD.n_explicit_outs, 1);
    fqei_vals.qi = all_vals(MOD.n_explicit_outs+1:MOD.n_explicit_outs+MOD.n_implicit_eqns, 1);
    fqei_vals.fe = all_vals(MOD.n_explicit_outs+MOD.n_implicit_eqns+1:2*MOD.n_explicit_outs+MOD.n_implicit_eqns, 1);
    fqei_vals.fi = all_vals(2*MOD.n_explicit_outs+MOD.n_implicit_eqns+1:end, 1);

    if (nargout < 2)
        return;
    end

    J = struct();
    if (flag.J)
        J.Jqe.dqe_dvecX = all_ders(1:MOD.n_explicit_outs, 1:MOD.n_other_ios); 
        J.Jqe.dqe_dvecY = all_ders(1:MOD.n_explicit_outs, 1+MOD.n_other_ios:MOD.n_other_ios+MOD.n_implicit_eqns);
        J.Jqe.dqe_dvecLim = sparse(MOD.n_explicit_outs, 0);
        J.Jqe.dqe_dvecU = sparse(MOD.n_explicit_outs, 0);

        J.Jqi.dqi_dvecX = all_ders(MOD.n_explicit_outs+1:MOD.n_explicit_outs+MOD.n_implicit_eqns, 1:MOD.n_other_ios);
        J.Jqi.dqi_dvecY = all_ders(MOD.n_explicit_outs+1:MOD.n_explicit_outs+MOD.n_implicit_eqns, 1+MOD.n_other_ios:MOD.n_other_ios+MOD.n_implicit_eqns);
        J.Jqi.dqi_dvecLim = sparse(MOD.n_implicit_eqns, 0);
        J.Jqi.dqi_dvecU = sparse(MOD.n_implicit_eqns, 0);

        J.Jfe.dfe_dvecX = all_ders(MOD.n_explicit_outs+MOD.n_implicit_eqns+1:2*MOD.n_explicit_outs+MOD.n_implicit_eqns, 1:MOD.n_other_ios);
        J.Jfe.dfe_dvecY = all_ders(MOD.n_explicit_outs+MOD.n_implicit_eqns+1:2*MOD.n_explicit_outs+MOD.n_implicit_eqns, 1+MOD.n_other_ios:MOD.n_other_ios+MOD.n_implicit_eqns);
        J.Jfe.dfe_dvecLim = sparse(MOD.n_explicit_outs, 0);
        J.Jfe.dfe_dvecU = sparse(MOD.n_explicit_outs, 0);

        J.Jfi.dfi_dvecX = all_ders(2*MOD.n_explicit_outs+MOD.n_implicit_eqns+1:end, 1:MOD.n_other_ios);
        J.Jfi.dfi_dvecY = all_ders(2*MOD.n_explicit_outs+MOD.n_implicit_eqns+1:end, 1+MOD.n_other_ios:MOD.n_other_ios+MOD.n_implicit_eqns);
        J.Jfi.dfi_dvecLim = sparse(MOD.n_implicit_eqns, 0);
        J.Jfi.dfi_dvecU = sparse(MOD.n_implicit_eqns, 0);
    end
end
