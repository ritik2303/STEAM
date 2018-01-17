function x_vec = reconstructWithCS(b, A, dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 COMPRESSED SENSING SPECIFIC FUNCTIONS                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds an x_vec s.t. its dct has the min l1-norm and satisfies 
%                           b = A*idct2(x_dct)
% where x_dct is the discrete cosine transform of x_vec. Here, b corresponds to
% the measurements of x_vec (nx1) and A is a fat matrix that is formed by
% selecting k (<n) rows of the identity matrix

    % In some cases the norm of b might be very small, it might be a good idea
    % to detect that case and stop wasting energy on it. The threshold is based
    % on observed values. The smallest values that we encounter are in qe/qi,
    % which are of the order ~ 1e-16

    b_max_val = max(abs((b)));
    if ( b_max_val < 1e-20)
        [M, N] = size(A);
        x_vec = zeros(N,1);
        return;
    end

    % Function handles to be passed to the optimization routine for conjugate
    % gradient method
    coeff_to_b = @(coeff) A*reshape(idct2(reshape(coeff, dims)), [], 1);
    b_to_coeff = @(b) reshape(dct2(reshape(A'*b, dims)), [], 1);

    x_dct_init_guess = b_to_coeff(b);

    x_dct_l1_min = l1qc_logbarrier(x_dct_init_guess, coeff_to_b, b_to_coeff, ...
        b, 1e-8*b_max_val, 1e-10, 10, 1e-10, 500);
    x_vec = idct2(reshape(x_dct_l1_min, dims));
end

