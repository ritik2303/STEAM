function tabulateModSpec(oMOD, iMOD, i_method, d_method, d_order, d_bounds, ...
    compression_factor)
% function tabulateModSpec(oMOD, iMOD, i_method, d_method, d_order, d_bounds, ...
%     compression_factor)
% Author: Archit Gupta (March 12, 2016)
% Update: March 11, 2017 (Adding new interpolant classes: (piecewise) Lagrange
% and piecewise BLI to the available options of interpolants. Also, Chebyshev
% sampling points have now been added
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             Tabulating an input ModSpec Object: Version 3.0
% This wrapper is specific for 2D models (only 3 or less terminal device can be
% tabulated using this wrapper. Here, the oMOD object over-writes the fqei
% functions, as well as the function handles for the derivates (trying to
% eliminate the use of vecvalder objects) in order to speed up the entire
% simulation process
%
% INPUTS:
%   iMOD: Input ModSpec object (should have atmost 3 terminals and no internal
%   unknowns)
%   interpMethod: Select from the available interpolation methods. This could be
%       1. spline - for bicubic coefficient based spline interpolation
%       2. passive - for bicubic coeffcient based passive spline interpolation
%       3. cosine - for a DCT based interpolation with linear extrapolation
%       4. lagrange - Piecwise lagrange interpolant
%       5. bli - Piecewise barycentric lagrange interpolant
%   d_method: Select the discretization method, i.e. which points do
%   we evaluate the iMOD on, in order to construct oMOD's table. Available
%   d_methods are:
%       1. 'u': Uniform Discretization
%       2. 'c': Chebyshev Points
%       3. 's': Compressed Sensing: The input space is discretized with a
%       evaluating the device functions at each of the sample points, a few
%       sample points are randomly chosen (ratio given by compression factor)
%   discretizationStep: This decides the degree of refineness in gridding the
%   voltage/input axes. For uniform discretization, this is the length of the
%   interval between 2 sampling points (e.g. 0.01 selects a 10mV discretization
%   step)
%   compression_factor: In case the d_method is chosen to be 's',
%   i.e. compressed sensing, this determines how fewer points will be sampled
%   relative to the uniform sampling case.
%
% OUTPUTS:
%   oMOD: Output ModSpec object. We will add a the following fields to the base
%   ModSpec class (name and functionality specified here):
%   1. inputSweep - sweep values for the inputs (need to know what the
%   inputs and what are the outputs)
%   2. lookUpTable - a 3 dimensional array (dimensionalty = #inputs) that
%   specifies the values of fe and qe for the inputSweep (fi/qi not supported)
%   3. @generate_LUT - function handle to regenrate the lookUpTable. This
%   can take in (optional) arguments to set the d_method in the
%   future.
%   4. In addition to these, we will add "interpMethod" and
%   "d_method" as parameters (use setparms/getparms).
%
%   Interpolation methods that are currently supported are: SPLINE, COSINES
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data_dir = ...
        '/home/ag/0_Unison/2_Academic/Spring2016/2_EE_290A_Numerical_Modelling_2/2_Project_Ideas/table_lookup_data/';
    defaultDiscretizationStep = 0.05;
    %checkInputModSpecSanity = check_ModSpec(iMOD);
    checkInputModSpecSanity = 1;
    if (~checkInputModSpecSanity)
        fprintf(2,'ERROR: Cannot verify input ModSpec sanity\n');
        oMOD = struct();
    end
    fprintf(2, 'Running tabulate-ModSpec Version 03...\n');

    % Beside the optional argument "compression factor", EVERYTHING has to be supplied % for the script to work. Have so
    % many possibilities in what you can provide makes confusing (and ambiguous) to use the script. We can have a
    % separate helper script to generate all the default arguments later on [TODO]

    narginchk(5, 6);

    % Setting up new parameter names and values
    oMOD.i_method = i_method;
    oMOD.d_method = d_method;
    oMOD.d_order = d_order;
    oMOD.d_bounds = d_bounds;

    % TODO: Current identifier works only when scalar values are provided (and are then expanded in replicated vectors
    % of appropriate size... We should make a hash out of the supplied data)
    str_order = strcat('_', sprintf('%-1d', oMOD.d_order));
    n_pieces = zeros(size( oMOD.d_bounds ));
    for np_i = 1 : size( oMOD.d_bounds, 2 )
        n_pieces(np_i) = size( oMOD.d_bounds{np_i}, 1 ) - 1;
    end
    str_pieces = strcat('_', sprintf('%-1d', n_pieces));
    oMOD.identifier = strcat('st_octave_', oMOD.i_method, '_', oMOD.d_method, str_pieces, str_order, '_', oMOD.identifier);

    switch nargin
        case 5
            if (d_method == 's')
                compression_factor = 4;
                fprintf(2, ['WARNING: No compression factor supplied. Using default', ...
                    '(%d)\n'], compression_factor);
            end
        case 6
            % DO NOTHING
        otherwise
            fprintf(2, 'Too many input arguments\n')
            oMOD = struct();
            return;
    end

    % Selecting an appropriate class handle for the given interpolation method

    i_handle = [];
    switch (oMOD.i_method)
    case 'spline'
        i_handle = @Spline2D;
    case 'passive'
        i_handle = @PSpline2D;
    case 'lagrange'
        i_handle = @PiecewiseLagrange2;
    case 'bli'
        i_handle = @PiecewiseBLI2;
    otherwise
        error('ERROR: Invalid interpolant requested. Aborting!\n');
    end

    fprintf(2, 'Discretization order is: ');
    disp(oMOD.d_order);
    
    global_FQ_varname = strcat('FEQE_', oMOD.identifier);
    base_FQ_filename = strcat(global_FQ_varname, '.mat');
    FQ_filename = strcat(data_dir, base_FQ_filename);
    eval(['global ',  global_FQ_varname, ';']);
    interp_exists = eval(['isa(', global_FQ_varname, ', ''Interpolant'') || isa(', global_FQ_varname, ...
        ', ''PiecewiseInterpolant'')']);
    if (interp_exists)
        fprintf(2, ['Found another model with the same Look Up Table %s,', ...
            'Skipping reading from memory\n'], base_FQ_filename);
        eval(['m_interp = ', global_FQ_varname, ';']);
    else
        if (exist(FQ_filename, 'file'))
            fprintf(2, 'Found Interpolation Object %s\n', base_FQ_filename);
            m_interp = i_handle();
            m_interp.load(FQ_filename);
            %fprintf(2, 'Table was computed in %ds\n', interpolant_eval_time);
        else
            tic;
            f_han = @(x_in) iMOD.fqei(x_in);
            m_interp = i_handle(f_han, oMOD.n_in_dims, d_bounds, d_order, d_method);
            interpolant_eval_time = toc;
            fprintf('Fit %s interpolant in %d s. Saving Interpolant Object\n', ...
                oMOD.i_method, interpolant_eval_time);
            % TODO: Add functionality to save the time required to generate this
            m_interp.save(FQ_filename);
        end
        eval([global_FQ_varname, ' = m_interp;']);
    end
    oMOD.m_interp = m_interp;
end
