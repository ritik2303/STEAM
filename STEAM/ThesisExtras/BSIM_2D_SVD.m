function [required_singular_values, reconstruction_error]  = BSIM_2D_SVD(...
    discretization_step, error_thresh, compression_factor)
% function [required_singular_values, reconstruction_error] = BSIM_2D_SVD(...
%   discretization_step, compression_factor)
% Author: Archit Gupta (Nov. 8, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   In this experiment, we look at the SVD of the Ids data for BSIM. For a given
%   discretization step size, we construct the Ids matrix and compute the
%   fraction of singular values required to reconstruct the matrix for a given
%   error
%   threshold
%   INPUTS:
%       discretization_step: discretization step size to be used for vds and vgs
%       error_thresh: Limit on Max relative error for the reconstruction to be
%       considered exact
%
%   OUTPUTS:
%       required_singular_values: #Singular-Values required for the max
%       reconstruction error to be less than error_thresh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    args = BSIMSamplingArgs();  
    run_recovery = 0;
    doPlots = 0;
    ERR_MODE = 'MAX_{REL}';
    reconstruction_error = 'INVALID';

    % Parameters for a 41x41 data matrix
    args.VMAX = 1.6;
    ags.VMIN = -1.55;
    if (nargin < 1)
        doPlots = 1;
        run_recovery = 1;
        discretization_step = 0.08;
        error_thresh = 1e-6;
    end

    if (nargin < 3 && run_recovery)
        compression_factor = 4;
    end

    [M, vds, vgs] = loadBSIMImage(discretization_step, args);
    n_entries_M = prod(size(M));
    size_M = size(M);
    
    % Singular Value decomposition of BSIM data
    [U, S, V] = svd(M);
    singular_values = diag(S);
    
    plot_line_width = 2.5;
    % Plot for the singular values of M
    if (doPlots)
        figure();
        plot_handle = stem(abs(singular_values), 'linestyle', 'none', ...
            'LineWidth', plot_line_width);
        grid on;
        set(gca, 'yscale', 'log');
        set(gca, 'FontSize', 28);
        xlabel('index');
        ylabel('Singular Value');
    end

    % From the SVD graph, we can see that the singular values drop dramatically.
    % We will now do a reconstruction by using k < N singular values, where N is
    % the size of the matrix

    %MAX_SINGULAR_VALUES = ceil(min(size_M)/2);
    MAX_SINGULAR_VALUES = min(size_M);
    reconstruction_error = zeros(MAX_SINGULAR_VALUES, 1);
    for n_singular_values = 1:MAX_SINGULAR_VALUES
        M_recon = U(:,1:n_singular_values)*S(1:n_singular_values,1:n_singular_values)*V(:,1:n_singular_values)'; 
        reconstruction_error(n_singular_values,1) = find_error(...
            M-M_recon, M, error_thresh, ERR_MODE);
    end

    % See where we can do almost as good as numerical precision (that is
    % reconstruction using all the singular values)
    required_singular_values = find(reconstruction_error < ...
        error_thresh, 1, 'first')/MAX_SINGULAR_VALUES;

    % Plotting reconstruction error
    if (doPlots)
        figure();
        plot_handle = stem(100*reconstruction_error, 'linestyle', ...
            'none', 'LineWidth', plot_line_width);
        xlabel('# Singular Values');
        ylabel('Max reconstruction Error (%)');
        set(gca, 'YScale', 'log');
        set(gca, 'FontSize', 28);
        grid on;
    end
end