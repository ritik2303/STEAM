function [fqout, dfqout] = fqei_dfqeiXYU(mod, x_in, forq, eori, XorY)
    if nargout > 1
        derivative_required = 1;
        if (nargin < 5)
            error('ERROR: Derivative required, but input missing');
        end
    else
        derivative_required = 0;
    end

    mod.flag.J = derivative_required;

    [fqei_vals, J] = mod.fqeiJ(x_in, mod.flag);

    if ((1 == strcmp(forq, 'f'))  &&  (1==strcmp(eori,'e')))
        fqout = fqei_vals.fe;
        if (derivative_required)
            dfqout_dX = J.Jfe.dfe_dvecX;
            dfqout_dY = J.Jfe.dfe_dvecY;
        end
    elseif ((1 == strcmp(forq, 'f'))  &&  (1==strcmp(eori,'i')))
        fqout = fqei_vals.fi;
        if (derivative_required)
            dfqout_dX = J.Jfi.dfi_dvecX;
            dfqout_dY = J.Jfi.dfi_dvecY;
        end
    elseif ((1 == strcmp(forq, 'q'))  &&  (1==strcmp(eori,'e')))
        fqout = fqei_vals.qe;
        if (derivative_required)
            dfqout_dX = J.Jqe.dqe_dvecX;
            dfqout_dY = J.Jqe.dqe_dvecY;
        end
    elseif ((1 == strcmp(forq, 'q'))  &&  (1==strcmp(eori,'i')))
        fqout = fqei_vals.qi;
        if (derivative_required)
            dfqout_dX = J.Jqi.dqi_dvecX;
            dfqout_dY = J.Jqi.dqi_dvecY;
        end
    else
        fprintf(2,'Invalid input for f or q\n');
        fqout = [];
        dfqout = [];
    end

    if (derivative_required)
        if 1 == strcmp(XorY, 'X')
            dfqout = dfqout_dX;
        elseif 1 == strcmp(XorY, 'Y')
            dfqout = dfqout_dY;
        else
            fprintf(2,'Invalid input for X or Y\n');
            dfqout = [];
        end
    end
end
