% Comparing Barycentric Lagrange Interpolant with Splines for a set of smooth
% functions, including a detailed Chebyshev series analysis.

disp('Comparing BLI and Splines over smooth functions.');
% Start with a smooth test function to interpolate
[fun,dfun_dx]  = getTestFHandle(1, 1, 'c');

disp('Generating arguments for the interpolant.');
args           = defaultInterpolantArgs();

disp('Instantiating low order interpolants.')
args{3}(:)     = 3;
bli_args       = args;
bli_args{4}    = 'c';
lo_bl_interpolant = BLI(fun, bli_args{:});
lo_sp_interpolant = Spline1D(fun, args{:});

disp('Instantiating high order interpolants.')
args{3}(:)     = 7;
bli_args       = args;
bli_args{4}    = 'c';
hi_bl_interpolant = BLI(fun, bli_args{:});
hi_sp_interpolant = Spline1D(fun, args{:});

n_test_pts     = 1000;
bounds         = args{2};

% Compare the interpolants at low order
[er, sp, plt] = compareIObjs(n_test_pts, bounds, fun, lo_bl_interpolant, 'LO-BLI', lo_sp_interpolant, 'LO-SPLINE');
set(gca, 'YScale', 'log');
compareIDers(n_test_pts, bounds, dfun_dx, lo_bl_interpolant, 'LO-BLI', lo_sp_interpolant, 'LO-SPLINE');

% Compare the interpolants at high order - Lots of sample points
[er, sp, plt] = compareIObjs(n_test_pts, bounds, fun, hi_bl_interpolant, 'HI-BLI', hi_sp_interpolant, 'HI-SPLINE');
set(gca, 'YScale', 'log');
compareIDers(n_test_pts, bounds, dfun_dx, hi_bl_interpolant, 'HI-BLI', hi_sp_interpolant, 'HI-SPLINE');

% Plot the Chebyshev Series for the functions
% This can be obtained directly from BLI - Do not need to explicitly fit a
% Chebyshev Series to the data.
hi_bl_interpolant.plotChebCoeffs();