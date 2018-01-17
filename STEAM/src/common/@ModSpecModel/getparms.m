function parms = getparms(MOD, parm_names, ~)
    if ( nargin < 3 )
        parms = MOD.parm_vals;
    else
        pidx = find(strcmp(parm_names, MOD.parm_names));
        parms = MOD.parm_vals{pidx};
    end
end
