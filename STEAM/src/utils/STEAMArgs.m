classdef STEAMArgs < handle
% Author: Archit Gupta (November 06, 2016), Bumped up to a class on March 17, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This class takes in an optional "method" argument and generates an object that
% can be used for running circuits with transistor models derived with
% compressed sensing (or otherwise)
%
% INPUTS:
%   method: Currently, the following methods are supported
%       1. spline
%       2. passive [TODO]
%       3. Lagrange
%       4. Bli
%
%   Of the methods mentioned above, only spline and passive splines (which are just another version of spline
%   interpolation) support compressed sensing. Please check the branch you are on. As of today, the method would only
%   work on BRANCH: ag-COSINES-with-CS
%
% OUTPUTS:
%   method_args (class object) with the following fields/properties:
%       interpMethod: (method used as default) 
%       discretizatonMethod: c
%       stepSize: 10mV
%       compression_factor: 4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        V_MIN = -1.5;
        V_MAX = 1.5;
    end

    properties ( SetAccess = private, GetAccess = public)
        n_dims             = 1; % Model Dimensions (input)
        i_method           = 'spline'; % Interpolation Method
        d_method           = 'chebyshev'; % Discretization Methods
        order              = 5;
        compression_factor = 4;
        bounds             = [-1.5; 1.5];
    end

    methods ( Access = private )
        function o_qty = resizeQty(Obj, i_qty) % [TODO].. Need Superclass for any object with n_dims, PUT THIS THERE
            if ( isscalar(i_qty) )
                fprintf(2, 'WARNING: Scalar quantity provided, expanding to appropriate size\n');
                o_qty = i_qty * ones(1, Obj.n_dims);
            else
                if (prod( size(i_qty) == [1 Obj.n_dims] ))
                    o_qty = i_qty;
                else
                    error('Mismatch in supplied order and model dimensions');
                end
            end
        end
    end

    methods ( Access = public )
        function Obj = STEAMArgs(n_dims, i_method, d_method, n_pieces_or_bounds)
            Obj.n_dims   = n_dims;
            Obj.i_method = i_method;
            Obj.d_method = d_method;
            Obj.order    = Obj.resizeQty(Obj.order);
            if ( iscell(n_pieces_or_bounds) )
                Obj.bounds = n_pieces_or_bounds;
            else
                Obj.bounds = cell(1, Obj.n_dims);
                n_pieces   = Obj.resizeQty(n_pieces_or_bounds);
                for d_i = 1 : Obj.n_dims
                    Obj.bounds{d_i} = linspace(Obj.V_MIN, Obj.V_MAX, 1 + n_pieces(d_i))';
                end
            end
        end

        function setOrder(Obj, n_order)
            Obj.order = Obj.resizeQty(n_order);
        end
    end
end
