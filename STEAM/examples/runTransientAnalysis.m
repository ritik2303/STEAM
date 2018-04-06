function [speedup, estimation_error, base_solution] = runTransientAnalysis(model, parm_string, n_pieces, order)
% function [speedup, estimation_error, base_solution] = runTransientAnalysis(model, parm_string, n_pieces, order)
% Author: Archit Gupta (March 23, 2016)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script lets you run a transient analysis on a dae. The transient values
% for the all the inputs must be setup beforehand. See daeMOSInverter for more
% details
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin == 0)
        model = @bsim3;
        parm_string = 'VA_BSIM3';
        n_pieces = 2;
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

    if (strfind( parm_string, 'BSIM' ))
        subCkt = @BSIM_with_RsRd;
    elseif (strfind( parm_string, 'MVS' ))
        subCkt = @MVS_with_RsRd; 
    elseif (strfind( parm_string, 'PSP' ))
        subCkt = @PSP_subckt;
    end
    
    % Ring-oscillator or inverter chain
    n_inverters  = 1;
    %{
    daeGenerator = @daeMOSInverterChain; 
    [steam_dae, steam_outs] = daeGenerator(n_inverters, model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeGenerator(n_inverters, model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, sim_args] = daeGenerator(n_inverters, model, parm_string, subCkt);
    %}

    % Any other circuit
    daeGenerator = @daeMOSSourceFollower;
    [steam_dae, steam_outs] = daeGenerator(model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeGenerator(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, sim_args] = daeGenerator(model, parm_string, subCkt);

    v2struct(sim_args);  % Creates tstart, tstep and tstop;

    % Transient Analysis on base circuit's dae

    % clear functions;
    base_filename = strcat(daeIdentifier, '_',  parm_string,'_tran.mat');
    if (exist(data_dir,'file'))
        analysis_filename = strcat(data_dir, base_filename);
    else
        analysis_filename = base_filename;
    end

    if (exist(analysis_filename))
        fprintf(2, 'Found previous simulation results in %s, SKIPPING TRAN\n', ...
            base_filename);
        load(analysis_filename);
    else
        fprintf(2, 'Running TRAN for base model\n');
        tic;
           base_lms_obj = dot_transient(base_dae, xinit, tstart, tstep, tstop);
        base_eval_time = toc;

        fprintf(2, 'Runnning TRAN for STEAM model\n');
        tic;
            steam_lms_obj = dot_transient(steam_dae, xinit, tstart, tstep, tstop);
        steam_eval_time = toc;
              
        fprintf(2, 'Runnning TRAN for BLI model\n');
        tic;
            bli_lms_obj = dot_transient(bli_dae, xinit, tstart, tstep, tstop);
        bli_eval_time = toc;

        %{
        fprintf(2, 'Saving TRAN results in %s\n', analysis_filename);
        % No need to have the DAE's in the LMS objects. It blows up the size unnecessarily
        base_lms_obj.DAE = [];
        base_lms_obj.AFobj = [];
        bli_lms_obj.DAE = [];
        bli_lms_obj.AFobj = [];
        steam_lms_obj.DAE = [];
        steam_lms_obj.AFobj = [];
        save(analysis_filename, 'base_lms_obj', 'base_eval_time', 'steam_lms_obj', 'steam_eval_time', 'bli_lms_obj', 'bli_eval_time', '-v7.3');
        %}
    end

    speedup.STEAM = base_eval_time/steam_eval_time;
    speedup.BLI = base_eval_time/bli_eval_time

    % Calculating the error in transient
    [steam_tpts, steam_unk_vals] = steam_lms_obj.getSolution(steam_lms_obj);
    [bli_tpts, bli_unk_vals] = bli_lms_obj.getSolution(bli_lms_obj);
    [base_tpts, base_unk_vals] = base_lms_obj.getSolution(base_lms_obj);
    steam_vals = steam_dae.output_matrix * steam_unk_vals;
    bli_vals = bli_dae.output_matrix * bli_unk_vals;
    base_vals = base_dae.output_matrix * base_unk_vals;

    steam_op_indices = steam_outs.OutputIndices(steam_outs);
    bli_op_indices = bli_outs.OutputIndices(bli_outs);
    base_op_indices = base_outs.OutputIndices(base_outs);

    % When we are dealing with very high accuracy, this interpolation error
    % might start kicking in and corrupting the accuracy that we get with
    % ALACARTE.
    steam_vals_at_base_pts = spline(steam_tpts', steam_vals(steam_op_indices,:), ...
        base_tpts');
    bli_vals_at_base_pts = spline(bli_tpts', bli_vals(bli_op_indices,:), base_tpts');
    base_solution = base_vals(base_op_indices,:);


    e_args = ErrorArgs();
    e_args.setDefaults();
    e_args.setEps(1e-8);

    steam_est_err_wf = steam_vals_at_base_pts - base_solution;
    bli_est_err_wf = bli_vals_at_base_pts - base_solution;
    estimation_error.STEAM = find_error(steam_est_err_wf, base_solution, e_args.eps(), ...
        e_args.errorMode());
    estimation_error.BLI = find_error(bli_est_err_wf, base_solution, e_args.eps(), ...
        e_args.errorMode())

    if (nargout < 2)
        % Plotting the two simulations in a single graph
        if (findstr(parm_string, 'BSIM'))
            BSIM_or_MVS = 'BSIM3.2.4';
        else
            BSIM_or_MVS = 'MVS1.0.1';
        end
        plot_varargin = {'LineWidth', 2.6};
        fprintf(2, 'Plotting TRAN simulation results\n');
        steam_time_units = steam_dae.time_units;
        [figh, onames, colindex] = transientPlot(steam_dae, steam_tpts, steam_unk_vals, ...
            steam_time_units, 'stateoutputs', steam_outs, 'lgndprefix', 'STEAM', ...
            'linestyle', '-.', 'plotvarargin', plot_varargin);

        bli_time_units = base_dae.time_units;
        [figh, onames, colindex] = transientPlot(bli_dae, bli_tpts, bli_unk_vals, ...
            bli_time_units, 'stateoutputs', bli_outs, 'lgndprefix', 'BLI', ...
            'linestyle', 'o-', 'plotvarargin', plot_varargin, 'fighandle', figh, ...
            'legends', onames, 'clrindex', colindex);

        base_time_units = base_dae.time_units;
        [figh, onames, colindex] = transientPlot(base_dae, base_tpts, base_unk_vals, ...
            base_time_units, 'stateoutputs', base_outs, 'lgndprefix', BSIM_or_MVS, ...
            'linestyle', 's-', 'plotvarargin', plot_varargin, 'fighandle', figh, ...
            'legends', onames, 'clrindex', colindex);

        % Correct the stuff in the plot
        figure(figh);
        set(gca, 'FontSize', 32);
        xlabel('TIME (in s)', 'FontSize', 28);
        ylabel('OUTPUT(s)', 'FontSize', 28);

        % Plotting the error in V_{OUT}
        fprintf(2, 'Plotting TRAN simulation error(s)\n');
        err_plots = figure; hold on;
        plot(base_tpts, steam_est_err_wf, plot_varargin{:}, 'linestyle', '-.');
        plot(base_tpts, bli_est_err_wf, plot_varargin{:}, 'linestyle', ':');
        xlabel('TIME (in s)', 'FontSize', 28);
        ylabel('\Delta{V}_{OUT} (volts)', 'FontSize', 28);

        plot_op_names = base_outs.OutputNames(base_outs);
        steam_plot_op_names = base_outs.OutputNames(base_outs);
        bli_plot_op_names = base_outs.OutputNames(base_outs);
        for index = 1:length(plot_op_names)
            steam_plot_op_names{index} = strcat('STEAM - \Delta{',escape_special_characters(...
                plot_op_names{index}),'}');
            bli_plot_op_names{index} = strcat('BLI - \Delta{',escape_special_characters(...
                plot_op_names{index}),'}');
        end
        grid on;
        legend(steam_plot_op_names{:}, bli_plot_op_names{:}, 'Location', 'best');
        set(gca, 'FontSize', 32);
    end
end