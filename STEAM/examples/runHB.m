function [speedup, estimation_error, freq_error] = runHB(...
    model, parm_string, n_pieces, order)
% function [speedup, estimation_error, base_solution] = runHB(...
%     model, parm_string, n_pieces, order)
% Author: Archit Gupta (March 26, 2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs a single-tone Harmonic balance on the following circuits:
%   1. 3-stage ring oscillator
%       - done
%
%   2. Gilbert Cell
%       - ongoing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin == 0)
        model = @bsim3;
        parm_string = 'VA_BSIM3';
        n_pieces = 4;
        order = 4;
    else
        narginchk(4,4);
    end

    global STEAM_DATA_DIR;

    % Arguments for running STEAM with compressed sensing
    method_bli = STEAMArgs(2, 'bli', 'chebyshev', n_pieces);
    method_bli.setOrder(order);

    spl_order = round( log10(n_pieces)/log10(2) ) + order;
    method_spl = STEAMArgs(2, 'spline', 'uniform', 1);
    method_spl.setOrder(spl_order);

    subCkt = getSubCkt(parm_string);
    % 3-Stage ring-oscillator
    is_osc = 1;
    %f0 = 5.9106e+06;
    f0 = 1.0406e+07;
    [steam_dae, steam_outs] = daeMOSRingosc(3, model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSRingosc(3, model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, sim_args] = daeMOSRingosc(3, model, parm_string, subCkt);
    daeIdentifier = '3_Stage_RingOscillator';
    %{
    [~, steam_tran_xinit] = ringOscInitializer(3);
    %}
    xinit = rand(steam_dae.nunks(steam_dae), 1);

    % Gilbert Cell
    %{
    is_osc = 0;
    f0 = 1e6;
    [steam_dae, steam_outs] = daeMOSGilbertCell(model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSGilbertCell(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, sim_args] = daeMOSGilbertCell(model, parm_string, subCkt);
    steam_tran_xinit = zeros(steam_dae.nunks(steam_dae), 1);
    bli_tran_xinit = zeros(bli_dae.nunks(bli_dae), 1);
    base_tran_xinit = zeros(base_dae.nunks(base_dae), 1);
    daeIdentifier = 'GilbertCell';
    %}

    % Estimate of base frequency (f0) for the oscillator
    M = 27;  % Number of harmonics for single-tone HB
    n_cycles = 4;
    n_pts_per_cycle = 30;
    n_cycles_to_skip = 2;
    do_plots_for_initguess = 1;

    v2struct(sim_args);  
    base_filename = strcat(daeIdentifier, '_',  parm_string,'_1T_', num2str(M), '_HARMS_HB.mat');
    initguess_base_filename = strcat(daeIdentifier, '_',  parm_string,'_INITGUESS_1T_', num2str(M), '_HARMS_HB.mat');

    if (exist(STEAM_DATA_DIR, 'dir'))
        analysis_filename = strcat(STEAM_DATA_DIR, base_filename);
        initguess_filename = strcat(STEAM_DATA_DIR, initguess_base_filename);
    else
        analysis_filename = base_filename;
        initguess_filename = initguess_base_filename;
    end

    if (exist(initguess_filename))
        fprintf(2, 'Found previous inital guess in %s, SKIPPING TRAN\n', ...
            initguess_base_filename);
        load(initguess_filename, '-mat');
    else
        % Generating an initial guess for Harmonic Balance
        base_xinitguess_Nn = HB_initguess(base_dae, f0, M, [], xinit, ...
            'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
            do_plots_for_initguess);

        bli_xinitguess_Nn = HB_initguess(bli_dae, f0, M, [], xinit, ...
            'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
            do_plots_for_initguess);

        % For some reason this transient doesn't work. Rest of HB works fine with STEAM too
        steam_xinitguess_Nn = HB_initguess(steam_dae, f0, M, [], xinit, ...
            'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
            do_plots_for_initguess);
        fprintf(2, 'Saving HB initial guess results from TRAN in %s\n', initguess_filename);
        save(initguess_filename, 'base_xinitguess_Nn', 'bli_xinitguess_Nn', 'steam_xinitguess_Nn', '-v7');
    end

    if (exist(analysis_filename))
        fprintf(2, 'Found previous simulation results in %s, SKIPPING HB\n', ...
            base_filename);
        load(analysis_filename, '-mat');
    else
        fprintf(2, 'Running HB with baseline model\n');
        base_hb_obj = HB(base_dae, is_osc);
        tic;
            base_hb_obj = base_hb_obj.solve(base_hb_obj, base_xinitguess_Nn, M, f0);
        base_eval_time = toc;

        fprintf(2, 'Running HB with BLI model\n');
        bli_hb_obj = HB(bli_dae, is_osc);
        tic;
            bli_hb_obj = bli_hb_obj.solve(bli_hb_obj, bli_xinitguess_Nn, M, f0);
        bli_eval_time = toc;

        fprintf(2, 'Running HB with STEAM model\n');
        steam_hb_obj = HB(steam_dae, is_osc);
        tic;
            steam_hb_obj = steam_hb_obj.solve(steam_hb_obj, steam_xinitguess_Nn, M, f0);
        steam_eval_time = toc;

        % Getting the solutions
        base_sol_struct = base_hb_obj.getsolution(base_hb_obj);
        bli_sol_struct = bli_hb_obj.getsolution(bli_hb_obj);
        steam_sol_struct = steam_hb_obj.getsolution(steam_hb_obj);

        fprintf(2, 'Saving HB results in %s\n', analysis_filename);
        save(analysis_filename, 'base_sol_struct', 'base_eval_time', 'bli_sol_struct', 'bli_eval_time', ...
            'steam_sol_struct', 'steam_eval_time', '-v7');
    end

    % Error calculation
    base_freq = base_sol_struct.frequency;
    bli_freq = bli_sol_struct.frequency;
    steam_freq = steam_sol_struct.frequency;

    freq_error = struct();
    freq_error.STEAM = abs(steam_freq - base_freq)/base_freq;
    freq_error.BLI = abs(bli_freq - base_freq)/base_freq

    e_args = ErrorArgs();
    e_args.setMode('MEAN_{REL}');
    e_args.setEps(1e-8);

    estimation_error = struct();
    estimation_error.STEAM = find_error(steam_sol_struct.X_twoD - base_sol_struct.X_twoD, base_sol_struct.X_twoD, ...
        e_args.eps(), e_args.errorMode());
    estimation_error.BLI = find_error(bli_sol_struct.X_twoD - base_sol_struct.X_twoD, base_sol_struct.X_twoD, ...
        e_args.eps(), e_args.errorMode())

    speedup = struct();
    speedup.STEAM = base_eval_time/steam_eval_time;
    speedup.BLI = base_eval_time/bli_eval_time

    % Plotting
    base_str_for_plot = ''; 
    steam_str_for_plot = 'STEAM ';
    bli_str_for_plot = 'BLI ';
    if (strfind( parm_string, 'BSIM' ))
        base_str_for_plot = 'BSIM3v3 ';
    elseif (strfind( parm_string, 'MVS' ))
        base_str_for_plot = 'MVS101 ';
    elseif (strfind( parm_string, 'PSP' ))
        base_str_for_plot = 'PSP1.03 ';
    end

    base_sol_struct.X_twoD = base_dae.output_matrix * base_sol_struct.X_twoD;
    base_sol_struct.x_twoD = base_dae.output_matrix * base_sol_struct.x_twoD;

    steam_sol_struct.X_twoD = steam_dae.output_matrix * steam_sol_struct.X_twoD;
    steam_sol_struct.x_twoD = steam_dae.output_matrix * steam_sol_struct.x_twoD;

    bli_sol_struct.X_twoD = bli_dae.output_matrix * bli_sol_struct.X_twoD;
    bli_sol_struct.x_twoD = bli_dae.output_matrix * bli_sol_struct.x_twoD;

    fharms = [0, 1:M, -M:-1];
    fprintf(2, 'Plotting HB solutions\n');
    [figh, onames, colindex]= plot_FD_TD(base_str_for_plot, base_outs.outputnames, base_sol_struct.f0, fharms, ...
        base_sol_struct.X_twoD(base_outs.outputindices, :), base_sol_struct.tpts, ...
        base_sol_struct.x_twoD(base_outs.outputindices, :), base_dae.time_units);
    %{ These lead to too many unnecessary plots
    [figh, onames, colindex] = plot_FD_TD(base_str_for_plot, steam_outs.outputnames, steam_sol_struct.f0, fharms, ...
        steam_sol_struct.X_twoD(steam_outs.outputindices, :), steam_sol_struct.tpts, ...
        steam_sol_struct.x_twoD(steam_outs.outputindices, :), steam_dae.time_units, figh, onames, colindex);
    [figh, onames, colindex] = plot_FD_TD(base_str_for_plot, bli_outs.outputnames, bli_sol_struct.f0, fharms, ...
        bli_sol_struct.X_twoD(bli_outs.outputindices, :), bli_sol_struct.tpts, ...
        bli_sol_struct.x_twoD(bli_outs.outputindices, :), bli_dae.time_units, figh, onames, colindex);
    %}


    % Plotting the error (in both the Time domain and the frequency domain)
    fprintf(2, 'Plotting HB Solution Error(s)\n');
    bli_err_struct = bli_sol_struct;
    bli_err_struct.X_twoD = (base_sol_struct.X_twoD - bli_sol_struct.X_twoD) ./ base_sol_struct.X_twoD;
    bli_err_struct.x_twoD = base_sol_struct.x_twoD - bli_sol_struct.x_twoD;
    plot_FD_TD(['\Delta_{', bli_str_for_plot, '}'], bli_outs.outputnames, ...
        bli_err_struct.f0, fharms, bli_err_struct.X_twoD(bli_outs.outputindices, :), bli_err_struct.tpts, ...
        bli_err_struct.x_twoD(bli_outs.outputindices, :), bli_dae.time_units); 

    steam_err_struct = steam_sol_struct;
    steam_err_struct.X_twoD = (base_sol_struct.X_twoD - steam_sol_struct.X_twoD) ./ base_sol_struct.X_twoD;
    steam_err_struct.x_twoD = base_sol_struct.x_twoD - steam_sol_struct.x_twoD;
    [e_figh, e_onames, colindex] = plot_FD_TD(['\Delta_{', steam_str_for_plot, '}'], steam_outs.outputnames, ...
        steam_err_struct.f0, fharms, steam_err_struct.X_twoD(steam_outs.outputindices, :), steam_err_struct.tpts, ...
        steam_err_struct.x_twoD(steam_outs.outputindices, :), steam_dae.time_units);
end
