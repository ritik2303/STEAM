function f_handle = getEvalHandle(Obj)
    if ( isa(Obj, 'function_handle') )
        f_handle = Obj;
    elseif ( isa(Obj, 'ModSpecModel') )
        f_handle = @(x_in) full(Obj.fe(x_in, [], []));
    elseif ( isa(Obj, 'Interpolant') | isa(Obj, 'PiecewiseInterpolant') )
        f_handle = @(x_in) Obj.computeWithDer( x_in );
    else
        error('ERROR: Invalid object type for extracting function handle\n');
    end
end
