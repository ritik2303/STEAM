function reconstruction_error = image_reconstruction(compression_factor)
% function reconstruction_error = image_reconstruction(compression_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this experiment, I will pick one of the sceneries from our Sandia trip,
% look at some of its pixels and try to reconstruct the entire image using these
% pixels alone. Intutively, this should be almost impossible to do because if
% you can't look at the entire image, pixel by pixel, it should not be possible
% to tell what these pixels are either
%
% RESULTS MATCH EXPECTATIONS: If the measurement matrix is a subset of columns
% of the identity matrix, DCT works and DWT doesn't
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin < 1)
        compression_factor = 2;
    end

    image_name = strcat('scene_', num2str(ceil(3*rand)));

    % Converting to grayscale to get the intensity at each pixel

    image = imread(image_name, 'png');
    image_gray = rgb2gray(image);    % Image size should be 332x256
    [M, N] = size(image_gray);

    patch_size = 256;    % Starting with a small patch
    x_patch_init = max(ceil(rand*(M-patch_size)), 1);
    y_patch_init = max(ceil(rand*(N-patch_size)), 1);   % When using patch_size
                                                        % = 256, the calculated
                                                        % value can be 0, hence
                                                        % the max() with 1
    I = double(image_gray(x_patch_init:x_patch_init+patch_size-1, ...
        y_patch_init:y_patch_init+patch_size-1));
    
    nPixels = prod(size(I));
    nPixelSamples = ceil(nPixels/compression_factor);

    sample_indices = randperm(nPixels, nPixelSamples)';
    %sample_indices = [1:ceil(nPixels/2)]';
    E = speye(nPixels);

    %%%%% INTERESTING: How do you generate an identity matrix with some (maybe
    % randomly chosen) rows zeroed out? 
    % ANSWER: E(sample_indices,:)'*E(sample_indices,:) does the trick
    % First note that all the non-diagonal entries are 0! This is because entry
    % (i,j) is obtained by multiplying the i^th column of E with the jth column
    % of E. If (i != j), this product is 0.
    % Next, if the k^th index is a part of sample indices, then the entry (k,k)
    % will be obtained by multiplying the k^th column of E with itself. Because
    % the k^th row if identity matrix has been selected in sample_indices, there
    % will be atleast 1 non-zero term in the k^th column of E. Similiar argument
    % can be made for indices that are not present in sample_indices

    A = E(sample_indices,:);
    b = A*reshape(I, [], 1);                    % Meaurement
    pixel_samples = reshape(A'*b, size(I));     % Visualization

    % Discrete Cosine Transform
    [C_dct] = dct2(I);
    sparsity_DCT = nnz(C_dct > 1e-4);

    % Wavelet Transform
    dwtmode('per');
    nLevels = 6;
    wavelet_basis = 'haar';

    % Generating the book-keeping vector required for wavelet reconstruction
    [C_dwt, S] = wavedec2(I, nLevels, wavelet_basis);

    sparsity_DWT = nnz(C_dwt > 1e-4);
    basis = 1;
    
    if (1 == basis) 
        coeff_to_b = @(coeff) A*reshape(idct2(reshape(coeff, size(I))), ...
            [], 1);
        b_to_coeff = @(b) reshape(dct2(reshape(A'*b, size(I))), [], 1);
    else
        coeff_to_b = @(coeff) A*reshape(waverec2(coeff, S, ...
            wavelet_basis), [], 1);
        b_to_coeff = @(b) reshape(wavedec2(reshape(A'*b, size(I)),...
            nLevels, wavelet_basis), [], 1);
    end

    % From here on, Ir is the reconstructed image
    Ir_dct_init_guess = b_to_coeff(b);

    algorithm = '[l1-magic] L1-Min';
    Ir_dct_l1min = l1qc_logbarrier(Ir_dct_init_guess, coeff_to_b, b_to_coeff, b, 1e-8, ...
        1e-6, 10, 1e-10, 500);
    Ir_dct_tvmin = tvqc_logbarrier(Ir_dct_init_guess, coeff_to_b, b_to_coeff, b, 1e-8, ...
        1e-6, 10, 1e-10, 500);

    Ir_image_l1 = idct2(reshape(Ir_dct_l1min, size(I)));
    Ir_image_tv = idct2(reshape(Ir_dct_tvmin, size(I)));

    figure(), imshow(uint8(I));
    title('Original Image');
    figure(), imshow(uint8(pixel_samples));
    title('Observed Image');
    figure(), imshow(uint8(Ir_image_l1));
    title('Reconstructed Image - L1 Min');
    figure(), imshow(uint8(Ir_image_tv));
    title('Reconstructed Image - TV Min');

    fprintf('DCT Sparsity = %d\n', sparsity_DCT);
    fprintf('DWT Sparsity = %d\n', sparsity_DWT);
end
