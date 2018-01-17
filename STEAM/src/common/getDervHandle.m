function f_handle = getDervHandle(Obj)
    if ( isa(Obj, 'function_handle') )
        f_handle = Obj;
    elseif ( isa(Obj, 'ModSpecModel') )
        f_handle = @(x_in) Obj.dfe_dvecX(x_in, [], [], 0);    % Eveything except x_in is dummy
    elseif ( isa(Obj, 'Interpolant') | isa(Obj, 'PiecewiseInterpolant') )
        f_handle = @(x_in) Obj.firstDerivative( x_in );
    else
        error('ERROR: Invalid object type for extracting function handle\n');
    end
end
