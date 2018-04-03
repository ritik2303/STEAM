n_discretization_steps = 10;
discretization_steps = round(logspace(-0.2, -2.1, n_discretization_steps), 3);
error_threshold = 1e-3;

required_fraction_of_singular_values = svd_compressibility( ...
    discretization_steps, error_threshold);
required_fraction_of_dct_coeffs = dct_compressibility( ...
    discretization_steps, error_threshold);

plot_line_width = 2.5;

figure(), stem(required_fraction_of_dct_coeffs, 'linestyle' , ...
    'none', 'LineWidth', plot_line_width);
hold on;
stem(required_fraction_of_singular_values, 'linestyle' , ...
    'none', 'LineWidth', plot_line_width);
legend('DCT', 'SVD');
set(gca, 'XTickLabel', int2str(1000*discretization_steps'));
xlabel('discretization step (mV)');
ylabel('Fraction of Coefficients');
grid on;
set(gca, 'FontSize', 28);

function n_vals = svd_compressibility(discretization_steps, error_threshold)
    n_discretization_steps = length(discretization_steps);
    required_fraction_of_singular_values = zeros(n_discretization_steps, 1);
    for d_step_index = 1:n_discretization_steps
        [required_fraction_of_singular_values(d_step_index, 1), ~] = ...
            BSIM_2D_SVD(discretization_steps(d_step_index), error_threshold, 1);
    end
    n_vals = required_fraction_of_singular_values;
end

function n_vals = dct_compressibility(discretization_steps, error_threshold)
    n_discretization_steps = length(discretization_steps);
    required_fraction_of_dct_coeffs = zeros(n_discretization_steps, 1);
    for d_step_index = 1:n_discretization_steps
        required_fraction_of_dct_coeffs(d_step_index, 1) = ...
            BSIM_2D_DCT(discretization_steps(d_step_index), error_threshold);
    end
    n_vals = required_fraction_of_dct_coeffs;
end