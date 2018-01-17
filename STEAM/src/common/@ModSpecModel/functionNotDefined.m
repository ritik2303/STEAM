function out = functionNotDefined(MOD, eori)
    if (nargin > 1)
        if (strcmp(eori,'e'))
            out = sparse(MOD.n_explicit_outs, 0);
        else
            out = sparse(MOD.n_implicit_eqns, 0);
        end
    else
        error('ERROR: Undefined Function called. Try dbstop at this line for Stack Trace');
    end
end
