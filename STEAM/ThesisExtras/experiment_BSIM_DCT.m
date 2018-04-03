function reconstruction_error = experiment_BSIM_DCT(discretization_step, ...
    error_thresh, compression_factor)
% function reconstruction_error = experiment_BSIM_DCT(discretization_step, ...
%     error_thresh, compression_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Previously, we saw that the reconstruction of BSIM data along some line in the
% vgs, vds plane was a little troublesome. It appeared that only the first few
% coefficients were non-zero. However, if we used only the first 40 coefficients
% to recover the BSIM data, there was a large mismatch at the initial few points
% and then the recovery with first 40 DCT coefficients seemed to do very well.
%
% Now, we have implemented several versions of DCT and we will be testing them
% with CS to see which transforms work and which don't and MORE IMPORTANTLY,
% WHY!
%
% INPUTS:
%   discretization_step: The discretization step size to be used for BSIM vds.
%   compression_factor: The sensing matrix has a size cnxn, where c is a 
%   fraction. compression factor specifies 1/c 
%   error_thresh: Arguments to be used for error calculation. See help find_error 
%   for more details 
% OUTPUTS: 
%   reconstruction_error: Single number denoting the error (MAX_ABS, MEAN_ABS, 
%   MAX_REL, MEAN_REL etc.) depending on the supplied error args 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    reconstruction_error = 'INVALID';
    args = BSIMSamplingArgs();
    ERR_MODE = 'MAX_{REL}';
    doPlots = 0;

    if (nargin < 1)
        discretization_step = 0.008;
        compression_factor = 10;
        error_thresh = 1e-4;
        doPlots = 1;
    end

    [x, vds, vgs] = loadBSIMSlice(discretization_step, args);
    n_samples = length(x);
    size_x = size(x);
    
    % In order to  be able to use tvqc, we need x to be an even perfect square!
    root = floor(sqrt(n_samples));
    n_samples = (root-rem(root,2))^2;
    x = x(1:n_samples,:);
    vds = vds(1:n_samples,:);
    vgs = vgs(1:n_samples,:);

    % Basis transformation - We have 4 different DCT implementations, namely,
    % DCT1, DCT2, DCT3 and DCT4. Let's test them out one by one. Besides, we are
    % checking the sanctity of the recovery using:
    % 1. All the DCT coefficients
    % 2. Some n_x_dct DCT coefficients. It might be a good idea to vary n_x_dct
    % over some large range (TODO: What  if we pick coefficients at random? Does
    % this have the same property as the DFT? It is essentially a DFT of some
    % longer sequence)

    %{
    [x_cos, dctmat] = m_DCT1(x);
    recovery_mat = dctmat;  % DCT1 is symmetric and orthogonal
    transform = 'DCT - I';
    %}

    [x_cos, dctmat] = m_DCT2(x);
    recovery_mat = dctmat'; % DCT2 is orthogonal. Infact inv(DCT2) = DCT3 
    transform = 'DCT - II';

    %{
    [x_cos, dctmat] = m_DCT3(x);
    recovery_mat = dctmat'; % DCT3 is orthogonal. inv(DCT3) = DCT2
    transform = 'DCT - III';

    [x_cos, dctmat] = m_DCT4(x);
    recovery_mat = dctmat;  % DCT 4 is symmetric and orthogonal
    transform = 'DCT - IV';
    %}
    
    % Basis transformation - We have 4 different DST implementations, namely,
    % DST1, DST2, DST3 and DST4. Let's test them out one by one. Besides, we are
    % checking the sanctity of the recovery using:
    % 1. All the DST coefficients
    % 2. Some nX_dst DST coefficients. It might be a good idea to vary nX_dst
    % over some large range (TODO: What  if we pick coefficients at random? Does
    % this have the same property as the DFT? It is essentially a DFT of some
    % longer sequence)

    %{
    [x_cos, dctmat] = m_DST1(x);
    recovery_mat = dctmat;  % DST1 is symmetric and orthogonal
    transform = 'DST - I';

    [x_cos, dctmat] = m_DST2(x);
    recovery_mat = dctmat'; % DST2 is orthogonal. Infact inv(DST2) = DST3 
    transform = 'DST - II';

    [x_cos, dctmat] = m_DST3(x);
    recovery_mat = dctmat'; % DST3 is orthogonal. inv(DST3) = DST2
    transform = 'DST - III';

    [x_cos, dctmat] = m_DST4(x);
    recovery_mat = dctmat;  % DST 4 is symmetric and orthogonal
    transform = 'DST - IV';
    %}
    
    % Checking sanity of transform
    x_trivial = recovery_mat * x_cos;

    % Using first n_x_dct DCT coefficients for recovery: We will vary n_x_dct in
    % some range and plot the error a function of n_x_dct. A few typical cases
    % will also be plotted to see where the error is!

    n_x_dct_max = n_samples;
    n_x_dct_min = floor(sqrt(0.5*n_samples));
    n_x_dct_skip = floor(sqrt(0.5*n_samples));

    n_x_dct_typical = floor((n_x_dct_max+n_x_dct_min)/2);

    x_partial_recovery = recovery_mat(:,1:n_x_dct_typical) * ...
        x_cos(1:n_x_dct_typical,1);

    dropped_coeffs = x_cos(n_x_dct_typical+1:end);

    plot_line_width = 2.5;
    if (doPlots)
        % For visualizing ourselves
        %figure, plot(vds, x, 'Color', 'blue');
        %hold on;
        %plot(vds, x_trivial, 'Color', 'red');
        %hold on;
        %plot(vds, x_partial_recovery, 'Color', 'black');

        % For publishing (Thicker lines)
        figure, plot(vds, x, '-.', 'Color', 'blue', 'Marker', 'o', ...
            'MarkerEdgeColor', 'blue', 'LineWidth', plot_line_width);
        hold on;
        plot(vds, x_trivial, 'Color', 'red', 'Marker', 's', ...
            'MarkerEdgeColor', 'red', 'LineWidth', plot_line_width);
        hold on;
        plot(vds, x_partial_recovery, 'Color', 'black', 'Marker', 'd', ...
            'MarkerEdgeColor', 'black', 'LineWidth', plot_line_width);
        legend('Original Signal', 'DCTinv(DCT(x))', 'DCTinv(TruncatedDCT(x))');
        xlabel('Sample index (k)', 'FontSize', 28);
        ylabel('Sample value x(k)', 'FontSize', 28);
        title(strcat('Stability of ''', transform,''': Truncated to first ''', ...
            num2str(n_x_dct_typical), ''' coefficients'));
        grid on;
        set(gca, 'FontSize', 28);
    end

    % Varying n_x_dct and plotting the MAX/MEAN reconstruction error. To do this,
    % we first reconstruct the waveform the minimal number of DCT coefficients
    % and add the deltas corresponding to the additional coefficients later on.

    dct_coeffs_to_try =  cat(1, 0, [n_x_dct_min:n_x_dct_skip:n_x_dct_max]');
    n_dct_coeffs_to_try = length(dct_coeffs_to_try);
    max_recovery_error = zeros(n_dct_coeffs_to_try-1, 1);
    mean_recovery_error = zeros(n_dct_coeffs_to_try-1, 1);
    dropped_energy = zeros(n_dct_coeffs_to_try-1, 1);
    TOTAL_ENERGY = sumsqr(x_cos);

    % What a mistake: Need to re-initialize this mistake since the variable had
    % been previously used
    x_partial_recovery = 0;
    % The first entry of n_dct_coeffs_to_try is 0
    for n_x_dct_index = 2:n_dct_coeffs_to_try
        dropped_coeffs = x_cos(dct_coeffs_to_try(n_x_dct_index)+1:end);
        dropped_energy(n_x_dct_index-1, 1) = sqrt(sumsqr(dropped_coeffs)/ ...
            TOTAL_ENERGY);
        recovered_indices = dct_coeffs_to_try(n_x_dct_index-1)+1:...
            dct_coeffs_to_try(n_x_dct_index);
        recovered_delta =  recovery_mat(:,recovered_indices)*x_cos(recovered_indices,1);
        x_partial_recovery = x_partial_recovery + recovered_delta;
        max_recovery_error(n_x_dct_index-1,1) = find_error( ...
            x-x_partial_recovery, x, error_thresh, 'MAX_{REL}');
        mean_recovery_error(n_x_dct_index-1,1) = find_error( ...
            x-x_partial_recovery, x, error_thresh, 'MEAN_{REL}');
    end

    % Removing the 1 from the beginning
    dct_coeffs_to_try = dct_coeffs_to_try(2:end);
    if (doPlots)
        figure, semilogy(100*dct_coeffs_to_try/n_samples, ...
        max_recovery_error, 'Color', 'blue', 'Marker', 'o', ...
        'MarkerEdgeColor', 'blue', 'LineWidth', plot_line_width);
        grid on;
        xlabel('(%) DCT Coefficients', 'FontSize', 28);
        ylabel('Max reconstruction Error', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        %set(gca, 'xdir', 'reverse');
        title(transform);

        figure, semilogy(100*dct_coeffs_to_try/n_samples, ...
            mean_recovery_error, 'Color', 'blue', 'Marker', ...
            'o', 'MarkerEdgeColor', 'blue', 'LineWidth', plot_line_width);
        grid on;
        xlabel('(%) DCT Coefficients', 'FontSize', 28);
        ylabel('Mean reconstruction Error', 'FontSize', 28);
        set(gca, 'FontSize', 28);
        %set(gca, 'xdir', 'reverse');
        title(transform);

        figure, semilogy(100*dct_coeffs_to_try/n_samples, ...
            100*dropped_energy, 'LineWidth', plot_line_width);
        xlabel('(%) DCT Coefficients', 'FontSize', 28);
        ylabel('Dropped Energy (%)', 'FontSize', 28);
        grid on;
        title(transform);
        set(gca, 'FontSize', 28);
        %set(gca, 'xdir', 'reverse');
    end

     % Observations - TODO: One idea could could be to use a fat "eye" matrix.
     % This would help reduce the burden of storing the sensing matrix
     % significantly. This is discussed in the paper titled Basis pursuit by
     % Candes
     nObservations = floor(n_samples/compression_factor);
     disp('Creating non-trivial measurement matrix ...');
     E = eye(n_samples);   % Identity matrix 
     A = E(randperm(n_samples, nObservations)',:);
     y = A*x;        % x refers to the time domain waveform

     % Recovering the DCT coefficients using the sub-sampled time-domain data 
     fp = dctmat'*A'*y;  % Initial  guess for the DCT coefficients
     %fp = l1eq_pd(fp, A*dctmat, [], y, 1e-8);
     %fp = l1decode_pd(fp, A*dctmat, [], y, 1e-8, 50, 1e-8, 200);
     %fp = tveq_newton(fp, A*dctmat, [], y, 1e-8);

    for n_iters = 0:50
        fprintf(2, 'ITERATION: %d\n', n_iters);
        %fp = tvqc_logbarrier(fp, A*dctmat, [], y, 10^(6-n_iters/10), 10^(-2), 2, 10^(-1-n_iters/10), 500);
        fp = l1qc_logbarrier(fp, A*dctmat, [], y, 10^(-1-n_iters/10), 10^(-2-n_iters/10), 2, 10^(-5-n_iters/10), 500);
    end

     xp = dctmat*fp;
     plot_time_domain_waveforms = 1;
     epsilon = 1e-18;
     if (doPlots)
         figure();
         stem(log10(abs(x_cos)+epsilon), '-.', 'Color', 'blue', 'Marker', 'o', ...
             'MarkerEdgeColor', 'blue', 'LineWidth', plot_line_width, ...
             'linestyle', 'none');
         hold on;
         stem(log10(abs(fp)+epsilon), 'Color', 'red', 'Marker', 's', 'MarkerEdgeColor', ...
             'red', 'LineWidth', plot_line_width, 'linestyle', 'none');
         hold on;
         legend('Original Signal', 'Recovery-DCT');
         set(gca, 'FontSize', 26);
         xlabel('Frequency Index (w)');
         ylabel('log(|X_{p}(w)|)');
         title('DCT Coeffecients');
         xlim([0 n_samples]);
         set(gca, 'FontSize', 26);

         % The waveforms are usually so close that it is very hard to see any
         % difference whatsoever. We will plot the error to see this difference

         figure();
         plot(fp-x_cos, 'Color', 'blue', 'Marker', 'o', ...
             'MarkerEdgeColor', 'blue', 'LineWidth', plot_line_width);
         hold on;
         xlabel('Frequency Index (w)');
         ylabel('\Delta{X_{p}}(w)');
         title('Error in recovery of DCT coeffecients');
         xlim([0 n_samples]);
         set(gca, 'FontSize', 26);
         grid on;
         
         if (plot_time_domain_waveforms == 1)
             figure();
             scatter(A*vds, A*x, 's', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', ...
             'blue', 'LineWidth', plot_line_width);
             hold on;
             plot(vds, xp, 'Color', 'red');
             legend('Measured Signal', 'Recovery');
             xlabel('vds (in Volts)');
             ylabel('X_{p}(k)');
             title('Measurement vs. Recovery');
             xlim([vds(1,1) vds(n_samples,1)]);
             set(gca, 'FontSize', 26);
             grid on;

             figure();
             plot(vds, x, '-.', 'Color', 'blue', 'Marker', 'o', 'MarkerEdgeColor', ...
             'blue', 'LineWidth', plot_line_width);
             hold on;
             plot(vds, xp, 'Color', 'red', 'Marker', 's', 'MarkerEdgeColor', ...
                 'red', 'LineWidth', plot_line_width);
             legend('Original Signal', 'Recovery');
             xlabel('vds (in Volts)');
             ylabel('X_{p}(k)');
             title('Recovered time-domain waveform');
             xlim([vds(1,1) vds(n_samples,1)]);
             set(gca, 'FontSize', 26);
             grid on;

             figure();
             plot(vds, xp-x, 'Color', 'blue', 'Marker', 'o', 'MarkerEdgeColor', ...
                 'blue', 'LineWidth', plot_line_width);
             hold on;
             xlabel('vds (in Volts)');
             ylabel('\Delta{X_{p}}(k)');
             title('Error in time-domain recovery');
             xlim([vds(1,1) vds(n_samples,1)]);
             set(gca, 'FontSize', 26);
             grid on;
         end
     end
     reconstruction_error = find_error(xp-x, x, error_thresh, ERR_MODE);
end