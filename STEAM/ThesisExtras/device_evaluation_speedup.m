%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Test script for comparing the speedup and accuracy with STEAM.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This compares a few selected compact models (BSIM) with their STEAM
% approximations for speedup and accuracy in evaluating the models alone. These
% are raw speedups that can be obtained in evaluating the fqei functions for the
% device and the error in computing these.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transistor model parameters
m_handle = @bsim3;
model_id = 'VA_BSIM3_NMOS';

% Trying a 2D Model
n_dims   = 2;
model    = MOS5TModel2D(m_handle, model_id);

% STEAM (discretization/interpolation) parameters
d_method = 'u';
d_bounds = cell( 1, n_dims );

% Interpolation parameters (Which interpolant to use post-discretization)
n_pieces = randi([1, 1], 1, n_dims);
i_method = 'spline';

for d_i = 1 : size(d_bounds, n_dims)
    d_bounds{1, d_i} = linspace(-1, 1, n_pieces(d_i)+1)';
end

% Draw a comparison between the original model and the STEAM generated
% model to see the differences between the two
n_test_pts     = 50;

d_orders_to_try = 2:5;
n_orders_to_try = length(d_orders_to_try);
steam_model_errors  = zeros(n_orders_to_try, 1);
steam_model_speedup = zeros(n_orders_to_try, 1);

% Create a new model for each different polynomial order and compare this
% model with the baseline
fprintf(2, 'Testing STEAM generated model(s) against original!\n')
for dor = 1:n_orders_to_try
    % Initializing a STEAM model
    d_order     = d_orders_to_try(dor) * ones( 1, n_dims);
    steam_model = STEAM(model, i_method, d_method, d_order, d_bounds);

    % Clear the function cache to remove any JIT effects that we might have
    % in the code. You might really dislike this, and I agree with you, but
    % speedup comparison in MATLAB is very strange!
    clear functions;

    % Get the error and speedup for this STEAM approximation
    [err, sp] = compareIObjs(n_test_pts, d_bounds, model, steam_model, 'STEAM');
    steam_model_errors(dor)  = err;
    steam_model_speedup(dor) = sp;
end

% TODO: Currently the x-axis labels are model ORDERs. The sample points are
% related to this as 2^order + 1. The discretization step size is given by
% <domain size>/<sample points>. We should be plotting the discretization
% step size as it is the most intuitive measure.

% Produce a bar graph for the modeling errors and speedup values.
figure, bar(steam_model_errors);
set(gca, 'FontSize', 28);
set(gca, 'YScale', 'log')
xlabel('Sample Points');
ylabel('Mean relative Error');
grid on;

figure, bar(steam_model_speedup);
xlabel('Sample Points');
ylabel('Speedup');
grid on;