function bestSplit(m_handle, m_identifier, m_generator, do_plots)
    if (nargin < 4)
        do_plots = 0;
        if (nargin == 0)
            % Default arguments
            m_generator = @MOS5TModel1D;
            m_identifier = 'VA_BSIM3_NMOS';
            m_handle = @bsim3;
            do_plots = 1;
        end
    end

    MOD = m_generator(m_handle, m_identifier);
    f_handle = getEvalHandle(MOD);
    i_type_bli = 'chebyshev';
    i_type_spl = 'uniform';
    n_test_pts = 1000;

    VMIN = -1;
    VMAX = 1;

    % Best split (hand-fit)
    % Best Split for BSIM
    bounds_single_piece = {[VMIN; -0.33625; -0.33620; -0.01; 0.0; 0.01; 0.5; VMAX]};

    % Best split for PSP
    %{
    bounds_single_piece = {[VMIN; -0.87520; -0.87515; -0.122; -0.121; -0.105; -0.095; 0.084; ...
        0.088; 0.1105; 0.1115; 0.18; 0.21; VMAX]};
    %}

    orders_to_try = [2 : 5];
    sp_speedup_spline = zeros(length(orders_to_try), 1);
    sp_speedup_bli = zeros(length(orders_to_try), 1);
    sp_err_spline = zeros(length(orders_to_try), 1);
    sp_err_bli = zeros(length(orders_to_try), 1);

    for iter = 1:length(orders_to_try)
        order_single_piece = orders_to_try(iter);
        spl_order = ceil( log10(length(bounds_single_piece{1})) / log10(2) ) + ...
            order_single_piece;
        bli_obj = PiecewiseBLI(f_handle, 1, bounds_single_piece, order_single_piece, i_type_bli);
        spl_obj = Spline1D(f_handle, 1, {[VMIN; VMAX]}, spl_order, i_type_spl);
        [est_sp_err, est_sp_speedup] = compareIObjs(n_test_pts, bounds_single_piece, MOD, bli_obj, 'BLI', spl_obj, 'SPLINE');
        sp_speedup_spline(iter, 1) = est_sp_speedup(2);
        sp_speedup_bli(iter, 1) = est_sp_speedup(1);
        sp_err_spline(iter, 1) = est_sp_err(2);
        sp_err_bli(iter, 1) = est_sp_err(1);
    end

    % Plotting the results obtained from the single piece experiment to look at the variation in the error
    n_pts_to_plot = length(bounds_single_piece{1}) * (2.^orders_to_try(:) + 1);

    if (do_plots)
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
        set(gca, 'YScale', 'log');
        grid on;

    end

    fprintf(2, 'Sample Points Used: \n');
    disp(n_pts_to_plot);

    fprintf(2, 'Error in BLI: \n');
    disp(sp_err_bli);

    fprintf(2, 'Error in SPLINE: \n');
    disp(sp_err_spline);
end
