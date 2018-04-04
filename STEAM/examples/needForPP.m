function needForPP(m_handle, m_identifier, m_generator, do_plots)
% function needForPP(m_handle, m_identifier, m_generator, do_plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A quick demonstration of the regions where compact model evaluation differs
% from a polynomial interpolant which was actually fit with some of the sampled
% data from the compact model evaluation.
%
% INPUT(s):
%   m_handle: Function handle for the model
%   m_identifier: Parameter string for identifying the parameters of the model
%   m_generator: A model generator function that can be used to take the model
%       parameter string, the model handle and create an appropriate executable
%       model. This is required because models like BSIM and MVS have 2 trivial
%       internal unknowns which we do not want to include in polynomial
%       interpolation.
%   do_plot: Should plots be done (true) or not (false).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    n_test_pts = 10000;

    VMIN = -1;
    VMAX = 1;

    % Uniformly Spaced Pieces
    bounds_single_piece = [VMIN; VMAX];
    % A. Overall analysis
    %{
    n_pieces_to_try = [4, 8, 16, 32];
    orders_to_try   = [4:2:8];
    %}

    % B. Variation with order alone
    %{
    n_pieces_to_try = [4];
    orders_to_try   = [2:2:8];
    %}

    % C. Variation with number of pieces alone
    n_pieces_to_try = [4, 8, 16, 32];
    orders_to_try   = [6];

    usp_speedup_spline = zeros(length(n_pieces_to_try), length(orders_to_try));
    usp_speedup_bli = zeros(length(n_pieces_to_try), length(orders_to_try));
    usp_err_spline = zeros(length(n_pieces_to_try), length(orders_to_try));
    usp_err_bli = zeros(length(n_pieces_to_try), length(orders_to_try));
    for p_iter = 1 : length(n_pieces_to_try)
        bounds = {linspace(VMIN, VMAX, 1+n_pieces_to_try(p_iter))'};
        for o_iter = 1 : length(orders_to_try)
            order_single_piece = orders_to_try(o_iter);
            spl_order = round( log10(n_pieces_to_try(p_iter))/log10(2) ) + order_single_piece;
            bli_obj = PiecewiseBLI(f_handle, 1, bounds, order_single_piece, i_type_bli);
            spl_obj = Spline1D(f_handle, 1, {[VMIN; VMAX]},  spl_order, i_type_spl);
            [est_err, est_speedup] = compareIObjs(n_test_pts, {bounds_single_piece}, MOD, bli_obj, 'BLI', spl_obj, 'SPLINE');
            usp_speedup_spline(p_iter, o_iter) = est_speedup(2);
            usp_speedup_bli(p_iter, o_iter) = est_speedup(1);
            usp_err_spline(p_iter, o_iter) = est_err(2);
            usp_err_bli(p_iter, o_iter) = est_err(1);
        end
    end

    n_pts_to_plot = [2.^orders_to_try(:) + 1] * n_pieces_to_try;
    order_labels = cellstr(int2str( 2.^orders_to_try(:) + 1 ));
    piece_labels = cellstr(int2str( n_pieces_to_try(:) ));
    [id1, id2] = meshgrid( grp2idx(order_labels), grp2idx(piece_labels) );
    data_labels = strcat(int2str(reshape(n_pts_to_plot', [], 1)), '(', order_labels(id1(:)), ',', piece_labels(id2(:)), ')');
    if (do_plots)
        figure();
        H_usp_err = bar(100*cat(2, reshape(usp_err_bli, [], 1), reshape(usp_err_spline, [], 1)), 0.95);
        xlabel('#Points ( pts pp, #pieces) ');
        ylabel('Error (%)');
        set(gca, 'XTickLabel', data_labels);
        set(gca, 'FontSize', 28);
        set(gca, 'YScale', 'log');
        legend('BLI', 'SPLINE');
        grid on;

        figure();
        H_usp_speedup = bar(cat(2, reshape(usp_speedup_bli, [], 1), reshape(usp_speedup_spline, [], 1)), 0.95);
        set(H_usp_speedup, 'BaseValue', 1); 
        xlabel('#Points ( pts pp, #pieces) ');
        ylabel('Speedup');
        set(gca, 'XTickLabel', data_labels);
        set(gca, 'FontSize', 28);
        %{
        set(gca, 'YScale', 'log');
        %}
        legend('BLI', 'SPLINE');
        grid on;
    end

    fprintf(2, 'Sample Points Used: \n');
    disp(n_pts_to_plot);

    fprintf(2, 'Error in BLI: \n');
    disp(usp_err_bli);

    fprintf(2, 'Error in SPLINE: \n');
    disp(usp_err_spline);
end
