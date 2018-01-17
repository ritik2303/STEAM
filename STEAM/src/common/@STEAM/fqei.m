function [fe, qe, fi, qi] = fqei(MOD, x_in, flag)
    [fqei_vals, ~] = MOD.fqeiJ(x_in, flag);
    fe = fqei_vals.fe; 
    qe = fqei_vals.qe;  
    fi = fqei_vals.fi;  
    qi = fqei_vals.qi;  
end
