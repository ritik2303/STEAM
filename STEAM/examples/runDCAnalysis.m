function [speedup, estimation_error, base_solution] = runDCAnalysis(...
    model, parm_string, n_pieces, order)
% function [speedup, estimation_error, base_solution] = runDCAnalysis(...
%     model, parm_string, n_pieces, order)
%
% Author: Archit Gupta (March 12, 2016).
% UPDATE (March 17, 2017): The script can now compare table-based model with
% BLI/Lagrange/DCT interpolants for interpolation with both the original
% model and amongst themselves
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script lets you run dc analysis on any dae which has a voltage source vin
% in the circuit (pre-EqnEngine stage). The simulation parameters, like start,
% step and stop values for voltage sweep should be supplied by the dae in a
% struct called simArgs.
%
% see daeMOSInverter for more details
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin == 0)
        model = @bsim3;
        parm_string = 'VA_BSIM3';
        n_pieces = 4;
        order = 4;
    else
        narginchk(4,4);
    end

    global data_dir;
    subCkt = [];

    % Arguments for running STEAM with compressed sensing
    method_bli = STEAMArgs(2, 'bli', 'chebyshev', n_pieces);
    method_bli.setOrder(order);

    spl_order = round( log10(n_pieces)/log10(2) ) + order;
    method_spl = STEAMArgs(2, 'spline', 'uniform', 1);
    method_spl.setOrder(spl_order);

    subCkt = getSubCkt(parm_string);

    % Running steam circuit (Maybe be tabulated/2-dimensional etc.)
    [steam_dae, steam_outs] = daeMOSDiffpair(model, parm_string, ...
        subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSDiffpair(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, simArgs] = daeMOSDiffpair(model, parm_string, subCkt);
    daeIdentifier = 'diffpair';

    %{
    [steam_dae, steam_outs] = daeMOSInverter(model, parm_string, ...
        subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSInverter(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, simArgs] = daeMOSInverter(model, parm_string, subCkt);
    daeIdentifier = 'inverter';
    %}

    v2struct(simArgs);
    dae_inputs = base_dae.uQSS(base_dae);
    vdd_val = dae_inputs(1);
    n_vdd_steps = base_dae.nunks(base_dae);
    vdd_step_range = [ linspace(0.01, vdd_val, n_vdd_steps) ];
    vdd_str = 'vdd:::E';

    % clear functions;
    % Running DC sweep on the tabulated model from STEAM
    tic;
        steamSweep = dcsweep(steam_dae, xinit, 'vin:::E', v_start:v_step:v_stop, vdd_str, ...
            vdd_step_range);
    tabulated_eval_time = toc;

    % clear functions;
    % Running DC sweep on the tabulated model from ALACARTE
    tic;
        bliSweep = dcsweep(bli_dae, xinit, 'vin:::E', v_start:v_step:v_stop, vdd_str, ...
            vdd_step_range);
    bli_eval_time = toc;

    % clear functions;
    % Running DC sweep on the base model
    base_filename = strcat(daeIdentifier, '_', parm_string,'_dc.mat');
    if (exist(data_dir,'file'))
        analysis_filename = strcat(data_dir, base_filename);
    else
        analysis_filename = base_filename;
    end

    if (exist(analysis_filename,'file'))
        fprintf(2, 'Found previous simulation results in %s, SKIPPING DC\n', ...
            base_filename);
        load(analysis_filename);
    else
        tic;
            baseSweep = dcsweep(base_dae, xinit, 'vin:::E', v_start:v_step:...
                v_stop, vdd_str, vdd_step_range);
        base_eval_time = toc;
        save(analysis_filename, 'baseSweep', 'base_eval_time');
    end 
    speedup.STEAM = base_eval_time/tabulated_eval_time;
    speedup.BLI = base_eval_time/bli_eval_time;
    speedup.STEAM_ITERS = baseSweep.total_iters / steamSweep.total_iters;
    speedup.BLI_ITERS = baseSweep.total_iters / bliSweep.total_iters;

    % Calculating the error in V_{OUT}
    [bliTpts, bliVals] = bliSweep.getSolution(bliSweep);
    [steamTpts, steamVals] = steamSweep.getSolution(steamSweep);
    [baseTpts, base_vals] = baseSweep.getSolution(baseSweep);

    bliOpIndices = bli_outs.OutputIndices(bli_outs);
    steamOpIndices = steam_outs.OutputIndices(steam_outs);
    base_op_indices = base_outs.OutputIndices(base_outs);

    % Check if voltage steps are increasing/decreasing: Same for STEAM and CS,
    % they come from the same dae generation script
    if (steamTpts(2)-steamTpts(1)>0)
        steam_vals_at_base_pts = spline(steamTpts, steamVals(steamOpIndices,:), ...
            baseTpts);
        bli_vals_at_base_pts = spline(bliTpts, bliVals(bliOpIndices,:), baseTpts);
    else
        steam_vals_at_base_pts = spline(-steamTpts, steamVals(steamOpIndices,:), ...
            -baseTpts);
        bli_vals_at_base_pts = spline(-bliTpts, bliVals(bliOpIndices,:), -baseTpts);
    end
    steam_estimation_error = steam_vals_at_base_pts - base_vals(base_op_indices,:);
    bli_estimation_error = bli_vals_at_base_pts - base_vals(base_op_indices,:);
    base_solution = base_vals(base_op_indices,:);
    
    e_args = ErrorArgs();
    e_args.setDefaults();
    e_args.setEps(1e-8);

    estimation_error.STEAM = find_error(steam_estimation_error, base_vals(base_op_indices,:), e_args.eps(), ...
        e_args.errorMode());
    estimation_error.BLI = find_error(bli_estimation_error, base_vals(base_op_indices,:), e_args.eps(), ...
        e_args.errorMode());

    if (nargout < 2)
        % Plotting the two simulations in a single graph
        if (strfind(parm_string, 'BSIM'))
            BSIM_or_MVS = 'BSIM3.2.4';
        else
            BSIM_or_MVS = 'MVS1.0.1';
        end
        plot_varargin = {'LineWidth', 2.6};
        steam_time_units = steam_dae.time_units;
        [figh, onames, colindex] = transientPlot(steam_dae, steamTpts, steamVals, ...
            steam_time_units, 'stateoutputs', steam_outs, 'lgndprefix', 'STEAM-', ...
            'linestyle', '--', 'plotvarargin', plot_varargin);

        bli_time_units = bli_dae.time_units;
        [figh, onames, colindex] = transientPlot(bli_dae, bliTpts, bliVals, ...
            bli_time_units, 'stateoutputs', bli_outs, 'lgndprefix', 'ALACARTE-', ...
            'linestyle', '-.', ...
            'plotvarargin', plot_varargin, 'fighandle', figh, 'legends', ...
            onames, 'clrindex', colindex);

        base_time_units = base_dae.time_units;
        [figh, ~, ~] = transientPlot(base_dae, baseTpts, base_vals, ...
            base_time_units, 'stateoutputs', base_outs, 'lgndprefix', BSIM_or_MVS, ...
            'linestyle', 'o-', 'plotvarargin', plot_varargin, 'fighandle', figh, ...
            'legends', onames, 'clrindex', colindex);

        % Correct the legends etc. in the plot
        figure(figh);
        set(gca, 'FontSize', 32);
        xlabel('V_{IN} (volts)', 'FontSize', 28);
        ylabel('V_{OUT} (volts)', 'FontSize', 28);

        err_plots = figure;
        plot(baseTpts, steam_estimation_error, plot_varargin{:}, 'Marker', 's');
        hold on;
        plot(baseTpts, bli_estimation_error, plot_varargin{:}, 'Marker', 'o');
        xlabel('V_{IN} (volts)');
        ylabel('\Delta{V}_{OUT} (volts)');
        set(gca, 'FontSize', 28);

        bli_plot_op_names = base_outs.OutputNames(base_outs);
        steam_plot_op_names = base_outs.OutputNames(base_outs);
        for index = 1:length(bli_plot_op_names)
            bli_plot_op_names{index} = strcat('ALACARTE \Delta{',escape_special_characters(...
                bli_plot_op_names{index}),'}');
            steam_plot_op_names{index} = strcat('STEAM \Delta{',escape_special_characters(...
                steam_plot_op_names{index}),'}');
        end
        grid on;
        legend(steam_plot_op_names{:}, bli_plot_op_names{:});
        set(gca, 'FontSize', 32);
    end
end