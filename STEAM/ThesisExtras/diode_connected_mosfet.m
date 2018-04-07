function diode_connected_mosfet(m_handle, m_identifier, m_generator, do_plots)
    if (nargin < 4)
        do_plots = 0;
        if (nargin == 0)
            do_plots = 1;
            m_generator = @MOS5TModel1D;

            m_handle = @bsim3;
            m_identifier = 'VA_BSIM3_NMOS';
            % In order to look at PSP, uncomment the following
            %{
                m_handle = @psp101;
                m_identifier = 'VA_PSP_NMOS';
            %}
        end
    end

    MOD = m_generator(m_handle, m_identifier);
    f_handle = getEvalHandle(MOD);
    i_type_bli = 'chebyshev';
    i_type_spl = 'uniform';
    n_test_pts = 10000;

    % VMIN = -0.33625;    % Brings out the discontinuity in the second derivative
    % VMAX = -0.33620;

    VMIN = -1;
    VMAX = 1; % Look at this to fully appreciate the problem

    %VMIN = -0.3;
    %VMAX = 1;

    % Single Polynomial comparison
    bounds_single_piece = {[VMIN; VMAX]};
    orders_to_try = [11 : 11];
    sp_speedup_spline = zeros(length(orders_to_try), 1);
    sp_speedup_bli = zeros(length(orders_to_try), 1);
    sp_err_spline = zeros(length(orders_to_try), 1);
    sp_err_bli = zeros(length(orders_to_try), 1);

    % This takes some non-negligible amount of time, and more importantly, we
    % would like to compute these results with Octave but plot them using MATLAB
    % becasue there, the plots look significantly better

    filename = [m_identifier, '_singlePPoly.mat'];
    if ( exist(filename) )
        fprintf(2, 'Found experimental data in %s. Skipping evaluation...\n', filename);
        load(filename, ...
            'sp_err_bli', ...
            'sp_err_spline', ...
            'sp_speedup_bli', ...
            'sp_speedup_spline', ...
            '-mat');
        m__outputSinglePPolyData(orders_to_try, sp_err_spline, sp_err_bli, ...
            sp_speedup_spline, sp_speedup_bli, do_plots);
        return;
    end

    for iter = 1:length(orders_to_try)
        order_single_piece = orders_to_try(iter);
        bli_obj = BLI(f_handle, 1, bounds_single_piece, order_single_piece, i_type_bli);
        spl_obj = Spline1D(f_handle, 1, bounds_single_piece, order_single_piece, i_type_spl);
        [est_sp_err, est_sp_speedup, err_plt] = compareIObjs(n_test_pts, bounds_single_piece, MOD, bli_obj, 'BLI', spl_obj, 'SPLINE');
        sp_speedup_spline(iter, 1) = est_sp_speedup(2);
        sp_speedup_bli(iter, 1) = est_sp_speedup(1);
        sp_err_spline(iter, 1) = est_sp_err(2);
        sp_err_bli(iter, 1) = est_sp_err(1);
    end

    % Plot the Chebyshev Series for the model
    bli_obj.plotChebCoeffs();
    m__outputSinglePPolyData(orders_to_try, sp_err_spline, sp_err_bli, ...
        sp_speedup_spline, sp_speedup_bli, do_plots);

    fprintf(2, 'Saving experimental data in %s.\n', filename);
    save(filename, ...
        'sp_err_bli', ...
        'sp_err_spline', ...
        'sp_speedup_bli', ...
        'sp_speedup_spline', ...
        '-v7');
end

function m__outputSinglePPolyData(orders_to_try, sp_err_spline, sp_err_bli, ...
    sp_speedup_spline, sp_speedup_bli, do_plots)

    % Plotting the results obtained from the single piece experiment to look at the variation in the error

    if (do_plots)
        n_pts_to_plot = 2.^orders_to_try(:) + 1;
        figure();
        H = bar(100*cat(2, sp_err_bli, sp_err_spline), 0.95);
        xlabel('# Sample Points');
        ylabel('Error (%)');
        legend('BLI', 'SPLINE');
        set(gca, 'XTickLabel', int2str(n_pts_to_plot));
        set(gca, 'FontSize', 28);
        set(gca, 'YScale', 'log');
        grid on;

        figure();
        H = bar(cat(2, sp_speedup_bli, sp_speedup_spline), 0.95);
        set(H, 'BaseValue', 1);
        xlabel('# Sample Points');
        ylabel('Speedup');
        legend('BLI', 'SPLINE');
        set(gca, 'XTickLabel', int2str(n_pts_to_plot));
        set(gca, 'FontSize', 28);
        %{
        set(gca, 'YScale', 'log');
        %}
        grid on;
    end

    fprintf(2, 'Sample Points Used: \n');
    disp(n_pts_to_plot);

    fprintf(2, 'Error in BLI: \n');
    disp(sp_err_bli);

    fprintf(2, 'Error in SPLINE: \n');
    disp(sp_err_spline);
end