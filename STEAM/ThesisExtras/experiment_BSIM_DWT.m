function reconstruction_error = experiment_BSIM_DWT(discretization_step, ...
    compression_factor, sampling_args, error_args)
% function reconstruction_error = experiment_BSIM_DWT(discretization_step, ...
%     compression_factor, sampling_args, error_args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In the previous experiment, with the discrete cosine transform, we saw that
% BSIM data that had been sampled along a line on the vds, vgs plane had sparse
% DCT coefficients. However, the recovery was very jagged! Some of the high
% frequency coefficients, which had a value close to 0 for the original function
% got significantly high values after the L1-recovery, making the recovered
% time-domain waveform very erroneous. We will repeat the same experiment with
% the wavelet basis to see its results. One of the probable reasons for the
% poor performance of the DCT could be the weakness of decay of the
% coefficients. While it seems that the coefficients are sparse to the eye, a
% plot of the absolute value of the coefficients on a log scale shows that the
% decay is indeed not very fast (for a required amount of reconstruction
% accuracy)
%
% INPUTS:
%   discretization_step: The discretization step size to be used for BSIM vds.
%   compression_factor: The sensing matrix has a size cnxn, where c is a 
%   fraction. compression factor specifies 1/c 
%   sampling_args: STRUCT that specifies the variables l and m. vds, vgs for the 
%   transistor are sampled as l*vds + m*vgs = 1; 
%   error_args: Arguments to be used for error calculation. See help find_error 
%   for more details 
% OUTPUTS: 
%   reconstruction_error: Single number denoting the error (MAX_ABS, MEAN_ABS, 
%   MAX_REL, MEAN_REL etc.) depending on the supplied error args 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    default_error_args.threshold = 1e-4;
    default_error_args.ERR_MODE = 'MAX_{REL}';

    doPlots = (nargout == 0);

    switch(nargin)
    case 0
        discretization_step = 0.01;
        compression_factor = 2;
        sampling_args = BSIMSamplingArgs();
        v2struct(default_error_args);
    case 2
        sampling_args = default_BSIM_sampling_args();
        v2struct(default_error_args);
    case 3
        v2struct(default_error_args);
    case 4
        v2struct(error_args);
    otherwise
        fprintf(2, ['ABORTING: Provide atleast 2 arguments or 4 or NONE to', ... 
            'use defaults\n']);
        return;
    end

    [x, vds, vgs] = loadBSIMSlice(discretization_step, sampling_args);
    n_samples = length(x);

    % Wavelet transform tends to be very painful if the length of x is an odd
    % number. This is because the outputs of DWT are two vectors, CA and CD,
    % which are equal in size.

    % In order to  be able to use tvqc, we need x to be an even perfect square!
    root = floor(sqrt(n_samples));
    n_samples = (root-rem(root,2))^2;
    x = x(1:n_samples,:);
    vds = vds(1:n_samples,:);
    vgs = vgs(1:n_samples,:);

    wavelet_filter_name = 'haar';
    %wavelet_filter_name = 'db2';
    %wavelet_filter_name = 'db4';
    transform = strcat('''', wavelet_filter_name, ''' wavelet');
    [x_ca_dwt, x_cd_dwt] = m_DWT(x, wavelet_filter_name);
    % ^^              ^^
    % approximation   detail
    % coefficient     coefficients
    %
    % The approximation coefficient is obtained by convolving the input signal
    % with a low pass filter and thereafter, downsampling by a factor of 2
    % (considering only the even entries). The detail coefficient on the other
    % hand is obtained by convolving the result with a high-pass filter and
    % downsampling by a factor of 2

    % Checking sanity of Wavelet Reconstruction
    x_trivial = idwt(x_ca_dwt, x_cd_dwt, wavelet_filter_name);
    fprintf(2, '|x_reconstructed| = %d, |x| = %d for trivial recovery\n', ... 
        length(x_trivial), n_samples); 
    x_trivial = x_trivial(1:n_samples,:);

    % Reconstruction using only a few coefficients: We have seen that the
    % coefficients CA and CD are mostly 0s. In this case, we will use some 'k'
    % coefficients from CA and CD and set eveything else to exactly 0.  We will
    % compare this reconstruction with the reconstruction using a small number
    % of DCT coefficients

    % INTERESTING INSIGHT: (TODO) Because the wavelets are localized in "space"
    % reconstruction with the "first few" wavelet coefficients is rather
    % pointless. This makes sense only in case of DFT or DCT, where we know that
    % for a smooth/continuous function, the coefficients decay. In case of the
    % Haar Wavelet, CA basically gives a local average, i.e. x(i) + x(i+1) and
    % CD give the local deltas, i.e. x(i) - x(i+1), scaled by sqrt(2). This is
    % the reason why the plots that we get here look "correct" for a wavefirn
    % that flattens out to 0 after a point but  wrong for everything else. CAs
    % will be zero only in the region that is close to 0, whereas CDs will be
    % close to 0 where the waveform is flat.

    N = n_samples;
    nX_dwt_MAX = floor(N/2);
    nX_dwt_MIN = floor(N/10);
    nX_dwt = floor((nX_dwt_MAX+nX_dwt_MIN)/2); 
    % Note that the number of coefficients being used here is half of what we
    % had used for DCT. This is because we will be taking nX_dwt coefficients
    % from both the approximation and the detail coefficients

    [TOTAL_COEFFICIENTS NDIMS] = size(x_ca_dwt);
    x_ca_partial_recovery = cat(1, x_ca_dwt(1:nX_dwt,:), zeros( ...
        TOTAL_COEFFICIENTS - nX_dwt, NDIMS));
    x_cd_partial_recovery = cat(1, x_cd_dwt(1:nX_dwt,:), zeros( ...
        TOTAL_COEFFICIENTS - nX_dwt, NDIMS));

    dropped_ca_coefficients = x_ca_dwt(nX_dwt+1:end,:);
    dropped_cd_coefficients = x_cd_dwt(nX_dwt+1:end,:);

    x_partial_recovery = idwt(x_ca_partial_recovery, x_cd_partial_recovery, ...
        wavelet_filter_name);
    fprintf(2, '|x_reconstructed| = %d, |x| = %d for partial recovery\n', ... 
        length(x_partial_recovery), N); 
    x_partial_recovery = x_partial_recovery(1:N,:);

    plot_line_width = 2.5;
    if (doPlots)
        figure();
        title('DWT coefficients for sample data');
        subplot(2,1,1), semilogy(abs(x_ca_dwt));
        xlim([1 ceil(n_samples/2)]);
        xlabel('Coefficient index (i)', 'FontSize', 28);
        ylabel('Approximation Coefficient ca(i)', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        grid on;
        subplot(2,1,2), semilogy(abs(x_cd_dwt));
        xlim([1 ceil(n_samples/2)]);
        xlabel('Coefficient index (i)', 'FontSize', 28);
        ylabel('Detail Coefficient cd(i)', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        grid on;

        figure, plot(vds, x, 'Color', 'blue', 'LineWidth', plot_line_width/2);
        hold on;
        plot(vds, x_trivial, 'Color', 'red', 'LineWidth', plot_line_width/2);
        hold on;
        plot(vds, x_partial_recovery, 'Color', 'black', 'LineWidth', ...
            plot_line_width/2);
        legend('Original Signal', 'DWTinv(DWT(x))', 'DWTinv(TruncatedDWT(x))');
        xlabel('Sample index (k)', 'FontSize', 28);
        ylabel('Sample value x(k)', 'FontSize', 28);
        title(strcat('Stability of ''', transform,''': Truncated to first ''', ...
            num2str(2*nX_dwt), ''' coefficients'));
        grid on;
        set(gca, 'FontSize', 28);
    end

    % Initializing variables that need to be looped over for error analysis
    max_recovery_error = zeros(nX_dwt_MAX-nX_dwt_MIN+1,1);
    mean_recovery_error = zeros(nX_dwt_MAX-nX_dwt_MIN+1,1);

    x_ca_partial_recovery = cat(1, x_ca_dwt(1:nX_dwt_MIN,:), zeros( ...
        TOTAL_COEFFICIENTS - nX_dwt_MIN, NDIMS));
    x_cd_partial_recovery = cat(1, x_cd_dwt(1:nX_dwt_MIN,:), zeros( ...
        TOTAL_COEFFICIENTS - nX_dwt_MIN, NDIMS));

    for nX_dwt = nX_dwt_MIN:nX_dwt_MAX
        x_ca_partial_recovery(nX_dwt,:) = x_ca_dwt(nX_dwt,:);
        x_cd_partial_recovery(nX_dwt,:) = x_cd_dwt(nX_dwt,:);
        x_partial_recovery = idwt(x_ca_partial_recovery, ...
            x_cd_partial_recovery, wavelet_filter_name);
        x_partial_recovery = x_partial_recovery(1:N,:); % Reconstruction always
        % has an even size, so this creates some trouble while plotting
        max_recovery_error(1+nX_dwt-nX_dwt_MIN,1) = find_error( ...
            x-x_partial_recovery, x, 1e-4, 'MAX_{REL}');
        mean_recovery_error(1+nX_dwt-nX_dwt_MIN,1) = find_error( ...
            x-x_partial_recovery, x, 1e-4, 'MEAN_{REL}');
    end

    if (doPlots)
        figure, semilogy(2*[nX_dwt_MIN:nX_dwt_MAX]', max_recovery_error, ...
            'Color', 'blue', 'Marker', 'o', 'MarkerEdgeColor', 'blue', ...
            'LineWidth', plot_line_width);
        grid on;
        xlabel('#DWT Coefficients', 'FontSize', 28);
        ylabel('Max reconstruction Error', 'FontSize', 28);
        title(transform);
        set(gca, 'FontSize', 28);
        figure, semilogy(2*[nX_dwt_MIN:nX_dwt_MAX]', mean_recovery_error, ...
            'Color', 'blue', 'Marker', 'o', 'MarkerEdgeColor', 'blue', ...
            'LineWidth', plot_line_width);
        grid on;
        xlabel('#DWT Coefficients', 'FontSize', 28);
        ylabel('Mean reconstruction Error', 'FontSize', 28);
        title(transform);
        set(gca, 'FontSize', 28);
    end

    % Rather trivial reconstruction scheme. Instead of using linear measurements
    % of the vector x, we use linear measurements of the DWT coefficients and
    % use a few of these to reconstruct all the DWT coefficients, and therefore,
    % the original vector
    
    wavelet_basis = 'haar';
    [CACD, S] = wavedec(x, 1, wavelet_basis);

    N = length(CACD);
    k = floor(N/compression_factor);

    A = randn(k, N);
    A = orth(A')';

    bT = A*CACD;
    cacdT_init_guess = A'*bT;
    cacdT_l1_min = tvqc_logbarrier(cacdT_init_guess, A, [], bT, 1e-8, 1e-6);
    %cacdT_l1_min = l1qc_logbarrier(cacdT_init_guess, A, [], bT, 1e-8, 1e-6);
    xT_l1_min = waverec(cacdT_l1_min, S, wavelet_basis);

    % Supplying compressed sensing to the data at hand. We know that that the
    % waveform that we have is sparse in the DWT basis. If we have a measurement
    % matrix M, s.t, M is incoherent with the DWT basis, then M_*IDWT*x_dwt,
    % where M_ is a matrix formed by a few rows of M, can be used to reconstruct
    % the entire vector x_dwt. More importantly, IDWT*x_dwt represents the
    % actual measurement "x", that we can obtain.

    % Recovering DWT coefficients using the subsampled (more like linear
    % measurements) of the waveform - One of the technical details in doing this
    % for the DWT is the handling of coefficients. DWT typically return the
    % approximation coefficients and the detail coefficients separately. We can
    % bundle them together into [ca, cd], which is always an even number. We
    % need to solve for the optimization problem:
    %
    %               min_x ||x||_1 s.t. A*idwt(cacd) = A*x, 
    %
    % where cacd is the concactenation of ca and cd.
    
    nObservations = floor(n_samples/compression_factor);
    disp('Creating non-trivial measurement matrix...');
    %A = randn(nObservations, n_samples);
    %A = orth(A')';  % All the papers on compressed sensing start with an
                    % orthogonal measurement matrix (in general an invertible
                    % matrix) because if A is not invertible, it is
                    % theoritically impossible to determine a vector x from the
                    % product Ax, as anything from the null space could be added
                    % to it. Can sparsity get rid of this too? TODO: Think
    E = eye(n_samples);   % Identity matric 
    A = E(randperm(n_samples, nObservations)',:);
    bNT = A*x;

    % TODO: Come up with a better intial guess
    cacdNT_init_guess = wavedec(A'*bNT, 1, wavelet_basis);
    cacdNT_to_y = @(cacd) A*waverec(cacd, S, wavelet_basis);
    y_to_cacdNT = @(y) wavedec(A'*y, 1, wavelet_basis);

    %cacdNT_l1_min = l1eq_pd(cacdNT_init_guess, cacd_to_y, y_to_cacd, b, 1e-5);
    %cacdNT_l1_min = l1qc_logbarrier(cacdNT_init_guess, cacdNT_to_y, ...
    %    y_to_cacdNT, bNT, 1e-8, 1e-6);
    cacdNT_l1_min = tvqc_logbarrier(cacdNT_init_guess, cacdNT_to_y, ...
        y_to_cacdNT, bNT, 1e-8, 1e-6);
    x_l1_min = waverec(cacdNT_l1_min, S, wavelet_basis);

    trivial_reconstruction_error = find_error(xT_l1_min-x, x, 1e-1, ...
        'MEAN_{REL}')
    non_trivial_reconstruction_error = find_error(x_l1_min-x, x, 1e-1, ...
        'MEAN_{REL}')

    figure, plot(x_l1_min-x);

    if (doPlots)
        figure();
        title('Comapring  DWT coefficients from TV Minimization');
        semilogy(abs(CACD), 'LineWidth', plot_line_width/2);
        hold on;
        semilogy(abs(cacdT_l1_min), 'LineWidth', plot_line_width/2);
        hold on;
        semilogy(abs(cacdNT_l1_min), 'LineWidth', plot_line_width/2);
        xlabel('Coefficient index (i)', 'FontSize', 28);
        ylabel('DWT Coefficient ca,cd(i)', 'FontSize', 28);
        legend('Original', 'Trivial TV-min', 'Non-Trivial TV-Min');
        set(gca, 'FontSize', 28);
        grid on;

        figure, plot(vds, x, 'Color', 'blue', 'LineWidth', plot_line_width/2);
        hold on;
        plot(vds, xT_l1_min, 'Color', 'red', 'LineWidth', plot_line_width/2);
        hold on;
        plot(vds, x_l1_min, 'Color', 'black', 'LineWidth', plot_line_width/2);
        legend('Original Signal', 'Trivial TV-Min', 'Non-Trivial TV-Min');
        xlabel('Sample index (k)', 'FontSize', 28);
        ylabel('Sample value x(k)', 'FontSize', 28);
        title('Reconstruction x(k) with TV Minimization');
        grid on;
        set(gca, 'FontSize', 28);
    end
