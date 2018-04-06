function runHBSweep(model, parm_string, n_pieces, order)
% function runHBSweep(model, parm_string, n_pieces, order)
% Author: Archit Gupta (March 26, 2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script runs a DC sweep on a Differential Amplifier circuit. For each DC
% operating point, it then performs Harmonic Balance and plots the DC, first,
% and second harmonics of the transfer function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargin < 2)
        model = @bsim3;
        parm_string = 'VA_BSIM3';
        %{
        model = @SH_MOS_ModSpec;
        parm_string = 'SH';
        %}
    else
        narginchk(2,4);
    end

    run_STEAM = (nargin > 2);

    global STEAM_DATA_DIR;

    if (run_STEAM)
        % Arguments for running STEAM
        method_bli = STEAMArgs(2, 'bli', 'chebyshev', n_pieces);
        method_bli.setOrder(order);

        spl_order = round( log10(n_pieces)/log10(2) ) + order;
        method_spl = STEAMArgs(2, 'spline', 'uniform', 1);
        method_spl.setOrder(spl_order);
    end

    sub_ckt = getSubCkt(parm_string);

    % Differential Pair
    is_osc = 0;
    [base_dae, base_outs, sim_args] = daeMOSDiffpair(model, parm_string, sub_ckt);
    base_tran_xinit = zeros(base_dae.nunks(base_dae), 1);

    if (run_STEAM)
        [steam_dae, steam_outs] = daeMOSDiffpair(model, parm_string, sub_ckt, method_spl);
        [bli_dae, bli_outs] = daeMOSDiffpair(model, parm_string, sub_ckt, method_bli);
        steam_tran_xinit = zeros(steam_dae.nunks(steam_dae), 1);
        bli_tran_xinit = zeros(bli_dae.nunks(bli_dae), 1);
    end
    daeIdentifier = 'DiffPair';

    % Estimate of base frequency (f0) for the oscillator
    M                      = 27;  % Number of harmonics for single-tone HB
    n_cycles               = 4;
    n_pts_per_cycle        = 30;
    n_cycles_to_skip       = 2;
    do_plots_for_initguess = 0;

    v2struct(sim_args); 
    if (run_STEAM)
        base_filename = strcat(daeIdentifier, '_',  parm_string,'_SWEEP_1T_', num2str(M), '_HARMS_HB.mat');
        initguess_base_filename = strcat(daeIdentifier, '_', ...
            parm_string,'_INITGUESS_SWEEP_1T_', num2str(M), '_HARMS_HB.mat');
    else
        base_filename = strcat(daeIdentifier, '_',  parm_string,'_SOLO_SWEEP_1T_', num2str(M), '_HARMS_HB.mat');
        initguess_base_filename = strcat(daeIdentifier, '_', ...
            parm_string,'_INITGUESS_SOLO_SWEEP_1T_', num2str(M), '_HARMS_HB.mat');
    end

    if (exist(STEAM_DATA_DIR, 'dir'))
        analysis_filename  = strcat(STEAM_DATA_DIR, base_filename);
        initguess_filename = strcat(STEAM_DATA_DIR, initguess_base_filename);
    else
        analysis_filename  = base_filename;
        initguess_filename = initguess_base_filename;
    end

    if (exist(initguess_filename))
        fprintf(2, 'Found previous inital guess in %s, SKIPPING TRAN\n', ...
            initguess_base_filename);
        load(initguess_filename, '-mat');
    else
        % Generating an initial guess for Harmonic Balance
        base_xinitguess_Nn = HB_initguess(base_dae, f0, M, [], base_tran_xinit, ...
            'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
            do_plots_for_initguess);

        if (run_STEAM)
            bli_xinitguess_Nn = HB_initguess(bli_dae, f0, M, [], bli_tran_xinit, ...
                'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
                do_plots_for_initguess);

            % steam_xinitguess_Nn = bli_xinitguess_Nn;
            % For some reason this transient doesn't work. Rest of HB works fine with STEAM too
            steam_xinitguess_Nn = HB_initguess(steam_dae, f0, M, [], steam_tran_xinit, ...
                'Transient', n_cycles, n_pts_per_cycle, n_cycles_to_skip, M, ...
                do_plots_for_initguess);

            fprintf(2, 'Saving HB initial guess results from TRAN in %s\n', initguess_filename);
            save(initguess_filename, 'base_xinitguess_Nn', 'bli_xinitguess_Nn', 'steam_xinitguess_Nn', '-v7');
        else
            fprintf(2, 'Saving HB initial guess results from TRAN in %s\n', initguess_filename);
            save(initguess_filename, 'base_xinitguess_Nn', '-v7.3');
        end
    end

    if (exist(analysis_filename))
        fprintf(2, 'Found previous simulation results in %s, SKIPPING HB\n', ...
            base_filename);
        load(analysis_filename, '-mat');
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run a DC sweep on the input vin
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qss_sweep    = dcsweep(base_dae, rand(base_dae.nunks(base_dae), 1), ...
                        'vin:::E', v_start:v_step:v_stop);
        [vins, sols] = qss_sweep.getSolution(qss_sweep);
        n_ins        = length(vins);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use the DC Sweep results to set up HB
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (run_STEAM)
            hb_sols   = cell(n_ins, 3);
        else
            hb_sols   = cell(n_ins, 1);
        end

        sinFDFunc = @(f, args) [args.offset, i*args.A/2, -i*args.A/2];

        % Need to automate the process of getting these parameters
        in_args.f = f0;
        in_args.A = 0.6e-3;
        in_args.Offset = 0.7;

        % Running Harmonic balance
        for v_in = 1 : n_ins
            % Set up the corresponding input
            in_args.offset = vins(v_in);
            base_dae       = base_dae.set_uHB('vin:::E', sinFDFunc, in_args, ...
                                base_dae);
            if (run_STEAM)
                steam_dae      = steam_dae.set_uHB('vin:::E', sinFDFunc, ...
                                    in_args, steam_dae);
                bli_dae        = bli_dae.set_uHB('vin:::E', sinFDFunc, in_args, ...
                                    bli_dae);
            end

            if (v_in ~= 1)
                % TODO: Use the previous solution
            end

            base_hb_obj = HB(base_dae, is_osc);
            tic;
                base_hb_obj = base_hb_obj.solve(base_hb_obj, base_xinitguess_Nn, M, f0);
            base_eval_time = toc;

            % Get the Harmonic Balance Solution at this Point
            hb_sols{v_in, 1}  = base_hb_obj.getsolution(base_hb_obj);

            if (run_STEAM)
                bli_hb_obj = HB(bli_dae, is_osc);
                tic;
                    bli_hb_obj = bli_hb_obj.solve(bli_hb_obj, bli_xinitguess_Nn, M, f0);
                bli_eval_time = toc;

                steam_hb_obj = HB(steam_dae, is_osc);
                tic;
                    steam_hb_obj = steam_hb_obj.solve(steam_hb_obj, steam_xinitguess_Nn, M, f0);
                steam_eval_time = toc;

                hb_sols{v_in, 2}  = bli_hb_obj.getsolution(bli_hb_obj);
                hb_sols{v_in, 3}  = steam_hb_obj.getsolution(steam_hb_obj);
            end
        end

        if (run_STEAM)
            save(analysis_filename, 'hb_sols', 'base_eval_time', 'bli_eval_time', ...
                'steam_eval_time', 'vins', 'n_ins', '-v7');
        else
            save(analysis_filename, 'hb_sols', 'base_eval_time', 'vins', ...
                'n_ins', '-v7');
        end
    end

    op_dims   = size(base_dae.output_matrix * hb_sols{1,1}.X_twoD);
    base_ops  = zeros([n_ins op_dims]);
    if (run_STEAM)
        bli_ops   = zeros([n_ins op_dims]);
        steam_ops = zeros([n_ins op_dims]);
    end

    % Rearrange the data in ops arrays
    for v_in = 1 : n_ins
        base_ops(v_in, :, :)   = abs(base_dae.output_matrix * hb_sols{v_in, 1}.X_twoD);
        if (run_STEAM)
            bli_ops(v_in, :, :)    = abs(bli_dae.output_matrix * hb_sols{v_in, 2}.X_twoD);
            steam_ops(v_in, :, :) = abs(steam_dae.output_matrix * hb_sols{v_in, 3}.X_twoD);
        end
    end

    plotHarmonicsForOutput(vins, base_ops, 1);
    if (run_STEAM)
        plotHarmonicsForOutput(vins, bli_ops, 1);
        plotHarmonicsForOutput(vins, steam_ops, 1);
    end
end

function plotHarmonicsForOutput(vins, sols, op_num)
    figure();
    subplot(3, 1, 1), scatter(vins, sols(:, op_num, 1), 'MarkerFaceColor', 'b');
    title('DC');
    set(gca, 'FontSize', 28);
    grid on;
    subplot(3, 1, 2), scatter(vins, sols(:, op_num, 2), 'MarkerFaceColor', 'b');
    title('1^{st} Harmonic');
    set(gca, 'FontSize', 28);
    grid on;
    subplot(3, 1, 3), scatter(vins, sols(:, op_num, 3), 'MarkerFaceColor', 'b');
    title('2^{nd} Harmonic');
    set(gca, 'FontSize', 28);
    grid on;
end