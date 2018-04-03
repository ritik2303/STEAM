function required_dct_coeffs = BSIMdct(discretization_step, error_thresh)
% function required_dct_coeffs = BSIMdct(discretization_step, error_thresh)
% Author: Archit Gupta (Nov. 8, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   In this experiment, we look at the DCT of the Ids data for BSIM. For a given
%   discretization step size, we construct the Ids matrix and compute the
%   fraction of dct coefficients required to reconstruct the matrix for a given
%   error
%   threshold
%   INPUTS:
%       discretization_step: discretization step size to be used for vds and vgs
%       error_thresh: Limit on Max relative error for the reconstruction to be
%       considered exact
%
%   OUTPUTS:
%       required_dct_coeffs: #DCT-coefficients required for the max
%       reconstruction error to be less than error_thresh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    args = BSIMSamplingArgs();  
    run_recovery = 0;
    doPlots = 0;
    ERR_MODE = 'MEAN_{REL}';

    % Parameters for a 41x41 data matrix
    args.VMAX = 1.6;
    ags.VMIN = -1.55;
    if (nargin < 1)
        doPlots = 1;
        discretization_step = 0.08;
        compression_factor = 4;
        error_thresh = 1e-4;
    end

    [M, vds, vgs] = loadBSIMImage(discretization_step, args);
    n_entries_M = prod(size(M));
    size_M = size(M);
    
    % Singular Value decomposition of BSIM data
    dct_coeffs  = dct2(M);
    dct_vec = reshape(dct_coeffs, [], 1);
    [dct_sorted, coeff_index] = sort(abs(dct_vec), 1, 'descend');

    plot_line_width = 2.5;
    % Plot for the singular values of M
    if (doPlots)
        figure();
        plot_handle = stem(abs(dct_vec), 'linestyle', 'none', ...
            'LineWidth', plot_line_width);
        grid on;
        set(gca, 'yscale', 'log');
        set(gca, 'FontSize', 28);
        xlabel('index');
        ylabel('DCT Coefficient');
    end

    % From the SVD graph, we can see that the singular values drop dramatically.
    % We will now do a reconstruction by using k < N singular values, where N is
    % the size of the matrix

    MAX_DCT_COEFFS = n_entries_M;
    MIN_DCT_COEFFS = ceil(0.5*sqrt(MAX_DCT_COEFFS));
    SKIP_DCT_COEFFS = ceil(0.5*sqrt(MAX_DCT_COEFFS));
    partial_dct_matrix = zeros(size_M);

    %dct_coeffs_to_try = cat(1, 0, [MIN_DCT_COEFFS:SKIP_DCT_COEFFS: ...
    %    MAX_DCT_COEFFS-SKIP_DCT_COEFFS]', [MAX_DCT_COEFFS-SKIP_DCT_COEFFS+1 ...
    %    :1:MAX_DCT_COEFFS]');
    dct_coeffs_to_try = cat(1, 0, [MIN_DCT_COEFFS:SKIP_DCT_COEFFS: ...
        MAX_DCT_COEFFS]');
    n_dct_coeffs_to_try = length(dct_coeffs_to_try);
    
    index_to_observe = ceil(length(dct_coeffs_to_try)/2)
    partial_reconstruction_error = zeros(n_dct_coeffs_to_try-1, 1);
    voltage_to_exclude_for_err = 100e-3;
    n_indices_to_exclude = ceil(voltage_to_exclude_for_err/discretization_step);
    error_calculation_range_x = [n_indices_to_exclude:size_M(1)-n_indices_to_exclude];
    error_calculation_range_y = [n_indices_to_exclude:size_M(2)-n_indices_to_exclude];

    M_to_check = M(error_calculation_range_x, error_calculation_range_y);
    for dct_coeff_index = 2:n_dct_coeffs_to_try
        recovered_indices = dct_coeffs_to_try(dct_coeff_index-1)+1:dct_coeffs_to_try(dct_coeff_index);
        partial_dct_matrix(coeff_index(recovered_indices)) = dct_coeffs(coeff_index(recovered_indices));
        M_recon = idct2(partial_dct_matrix);
        M_recon_to_check = M_recon(error_calculation_range_x, error_calculation_range_y);
        partial_reconstruction_error(dct_coeff_index-1,1) = find_error(...
            M_to_check-M_recon_to_check, M_to_check, error_thresh, ERR_MODE);
        if (doPlots && (dct_coeff_index == index_to_observe))
            figure(), surf(M_recon_to_check);
            freezeColors;
            hold on;
            surf(M_to_check);
            colormap('spring');
            xlabel('v_{ds}');
            ylabel('v_{gs}');
            zlabel('I_{ds}');
            set(gca, 'FontSize', 28);
            legend('Reconstructed', 'Original');
            title('I_{ds} Matrix: Reconstructed and Original');
        end
    end

    % See where we can do almost as good as numerical precision (that is
    % reconstruction using all the singular values)
    dct_coeffs_to_try = dct_coeffs_to_try(2:end);
    n_dct_coeffs_to_try = length(dct_coeffs_to_try);
    required_dct_coeffs = find(partial_reconstruction_error < ...
        error_thresh, 1, 'first')/n_dct_coeffs_to_try;

    % Plotting reconstruction error
    if (doPlots)
        figure();
        plot_handle = stem(dct_coeffs_to_try, 100*partial_reconstruction_error, 'linestyle', ...
            'none', 'LineWidth', plot_line_width);
        xlabel('# DCT Coefficients');
        ylabel(strcat(ERR_MODE, 'reconstruction Error (%)'));
        set(gca, 'yscale', 'log');
        set(gca, 'FontSize', 28);
        grid on;
    end

    % Instead of constructing the indicator matrix, we can just pick k out of
    % N^2 entries of the M matrix for sampling

    if (run_recovery)
        sample_indices = sort(randsample(n_entries_M, ceil(n_entries_M/compression_factor)));
        M_sampled = nan(size_M);
        M_sampled(sample_indices) = M(sample_indices);
        
        if (doPlots)
            figure(), surf(vds, vgs, M_sampled');
            hold on;
            freezeColors;
            %surf(vds, vgs, M');
            %colormap('spring');
            %legend('Sampled', 'Original');
            xlabel('v_{ds}');
            ylabel('v_{gs}');
            zlabel('I_{ds}');
            set(gca, 'FontSize', 28);
            title('I_{ds} Matrix: Sampled and Original');
        end

        % Using CVX 
        observations = M(sample_indices);
        mu = 2e-5; % Smooting parameter: We might have to play around with this a
                    % little in order to get some accuracy

        tic;
            M_recon = solver_sNuclearBP( {size_M(1), size_M(2), sample_indices}, observations, mu);
        solver_time =  toc;

        % This package doesn't work... We should try something else (CVX maybe)
        % indicator_matrix = (rand(size(M)) < 1/compression_factor);
        %[M_recon, success] = MatrixCompletion(M.*indicator_matrix, ...
        %    indicator_matrix, 10000, 'nuclear', 1e-14, 1e-10, 1);

        figure, surf(vds, vgs, M_recon);
        hold on;
        freezeColors;
        surf(vds, vgs, M');
        colormap('spring');
        legend('Reconstructed', 'Original');
        xlabel('v_{ds}');
        ylabel('v_{gs}');
        zlabel('I_{ds}');
        set(gca, 'FontSize', 28);
        title('I_{ds} Matrix: Reconstructed and Original');

        error = find_error(M-M_recon, M, error_thresh, ERR_MODE)
    end

end
