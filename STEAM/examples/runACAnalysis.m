function [speedup, estimation_error, baseSolution] = runACAnalysis( ...
    model, parm_string, n_pieces, order)
% function [speedup, estimation_error, baseSolution] = runACAnalysis( ...
%   model, parm_string, n_pieces, order)
% Author: Archit Gupta, April 20, 2016
% Updated: March 17, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script lets you run AC analysis on a _dae.
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

    if (strfind(parm_string,'BSIM'))
        subCkt = @BSIM_with_RsRd;
    elseif (strfind(parm_string, 'MVS'))
        subCkt = @MVS_with_RsRd; 
    end
    
    % Arguments for running STEAM with compressed sensing
    method_bli = STEAMArgs(2, 'bli', 'chebyshev', n_pieces);
    method_bli.setOrder(order);

    spl_order = round( log10(n_pieces)/log10(2) ) + order;
    method_spl = STEAMArgs(2, 'spline', 'uniform', 1);
    method_spl.setOrder(spl_order);

    % Running steam circuit (Maybe be tabulated/2-dimensional etc.)
    %{
    [steam_dae, steam_outs, sim_args] = daeMOSGilbertCell(model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSGilbertCell(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs] = daeMOSGilbertCell(model, parm_string, subCkt);
    daeIdentifier = 'GilbertCell';
    %}

    [steam_dae, steam_outs] = daeMOSInverter(model, parm_string, subCkt, method_spl);
    [bli_dae, bli_outs] = daeMOSInverter(model, parm_string, subCkt, method_bli);
    [base_dae, base_outs, sim_args] = daeMOSInverter(model, parm_string, subCkt);
    daeIdentifier = 'inverter';

    v2struct(sim_args);
    fstart = 1;
    fstop = 1e9;
    nsteps = 5;
    vstep = (VIN_MAX-VIN_MIN)/n_ops;

    steam_ac_time = zeros(n_ops,1);
    steam_freqs = cell(n_ops,1);
    steam_vals = cell(n_ops,1);

    base_ac_time = zeros(n_ops,1);
    base_freqs = cell(n_ops,1);
    base_vals = cell(n_ops,1);

    bli_ac_time = zeros(n_ops,1);
    bli_freqs = cell(n_ops,1);
    bli_vals = cell(n_ops,1);

    sweeptype = 'DEC';

    % Function for setting uLTISSS
    Ufargs.string = 'NO ARGUMENTS REQUIRED';
    Uffunc = @(f, args) ones(size(f)); 
    steam_dae = steam_dae.set_uLTISSS('vin:::E', Uffunc, Ufargs, steam_dae);
    bli_dae = bli_dae.set_uLTISSS('vin:::E', Uffunc, Ufargs, bli_dae);
    base_dae = base_dae.set_uLTISSS('vin:::E', Uffunc, Ufargs, base_dae);

    VIN = [VIN_MIN:vstep:VIN_MAX]';
    VIN = VIN(1:n_ops,1);    % Handling some exception cases, like n_ops=1

    % AC ANALYSIS FOR BASE MODEL
    base_filename = strcat(daeIdentifier, '_',  parm_string,'_ac.mat');
    if (exist(data_dir,'file'))
        analysis_filename = strcat(data_dir, base_filename);
    else
        analysis_filename = base_filename;
    end

    if (exist(analysis_filename))
        fprintf(2, 'Found previous simulation results in %s, SKIPPING AC\n', ...
            base_filename);
        load(analysis_filename);
    else
        dcop = xinit;
        for op_index = 1:n_ops
            %clear functions;
            base_dae = base_dae.set_uQSS('vin:::E', VIN(op_index,1), base_dae);
            vin = base_dae.uQSS(base_dae);
            tic;
                base_opt_pt = op(base_dae, dcop);
                dcop = base_opt_pt.getSolution(base_opt_pt);
                AC = LTISSS(base_dae, dcop, vin);
                AC = AC.solve(fstart, fstop, nsteps, sweeptype, AC);
            base_ac_time(op_index,1) = toc;
            [base_freqs{op_index}, base_vals{op_index}] = AC.getSolution(AC);
            %AC.plot(AC, base_outs);
        end
        save(analysis_filename, 'base_freqs', 'base_vals', 'base_ac_time');
    end

    
    % AC ANALYSIS FOR STEAM MODEL
    dcop = xinit;
    for op_index = 1:n_ops
        %clear functions;
        steam_dae = steam_dae.set_uQSS('vin:::E', VIN(op_index,1), steam_dae);
        vin = steam_dae.uQSS(steam_dae);
        tic;
            steam_opt_pt = op(steam_dae, dcop);
            dcop = steam_opt_pt.getSolution(steam_opt_pt);
            AC = LTISSS(steam_dae, dcop, vin);
            AC = AC.solve(fstart, fstop, nsteps, sweeptype, AC);
        steam_ac_time(op_index,1) = toc;
        [steam_freqs{op_index}, steam_vals{op_index}] = AC.getSolution(AC);
        %AC.plot(AC, steam_outs);
    end

    % AC ANALYSIS FOR BLI MODEL
    dcop = xinit;
    for op_index = 1:n_ops
        %clear functions;
        bli_dae = bli_dae.set_uQSS('vin:::E', VIN(op_index,1), bli_dae);
        vin = bli_dae.uQSS(bli_dae);
        tic;
            bli_opt_pt = op(bli_dae, dcop);
            dcop = bli_opt_pt.getSolution(bli_opt_pt);
            AC = LTISSS(bli_dae, dcop, vin);
            AC = AC.solve(fstart, fstop, nsteps, sweeptype, AC);
        bli_ac_time(op_index,1) = toc;
        [bli_freqs{op_index}, bli_vals{op_index}] = AC.getSolution(AC);
        %AC.plot(AC, bli_outs);
    end
    speedup.STEAM = mean(base_ac_time)/mean(steam_ac_time);
    speedup.BLI = mean(base_ac_time)/mean(bli_ac_time);
    
    % Plotting the Error (and also the phase-magnitude responses together)
    if (strfind(parm_string, 'BSIM'))
        BSIM_or_MVS = 'BSIM3.2.4: ';
    else
        BSIM_or_MVS = 'MVS1.0.1: ';
    end

    steam_op_indices = steam_outs.OutputIndices(steam_outs);
    bli_op_indices = bli_outs.OutputIndices(bli_outs);
    base_op_indices = base_outs.OutputIndices(base_outs);

    plot_op_names = base_outs.OutputNames(base_outs);

    steam_plot_op_names = base_outs.OutputNames(base_outs);
    bli_plot_op_names = base_outs.OutputNames(base_outs);
    base_plot_op_names = base_outs.OutputNames(base_outs);
    for index = 1:length(plot_op_names)
        steam_error_plot_op_names{index} = strcat('\Delta{', escape_special_characters(...
            plot_op_names{index}),'-STEAM}');
        bli_error_plot_op_names{index} = strcat('\Delta{', escape_special_characters(...
            plot_op_names{index}),'-BLI}');

        steam_plot_op_names{index} = strcat('STEAM: ', escape_special_characters(...
            plot_op_names{index}));
        bli_plot_op_names{index} = strcat('BLI: ', escape_special_characters(...
            plot_op_names{index}));
        base_plot_op_names{index} = strcat(BSIM_or_MVS, escape_special_characters(...
            plot_op_names{index}));
    end
    
    e_args = ErrorArgs();
    e_args.setDefaults();
    op_pt_error = zeros(n_ops,2);
    for op_index = 1:n_ops
        steam_vals{op_index} = squeeze(steam_vals{op_index});
        steam_vals{op_index} = steam_vals{op_index}(steam_op_indices,:); 

        bli_vals{op_index} = squeeze(bli_vals{op_index});
        bli_vals{op_index} = bli_vals{op_index}(bli_op_indices,:); 

        base_vals{op_index} = squeeze(base_vals{op_index});
        base_vals{op_index} = base_vals{op_index}(base_op_indices,:); 

        steam_vals_at_base_pts = spline(steam_freqs{op_index}, steam_vals{op_index}, ...
            base_freqs{op_index});
        bli_vals_at_base_pts = spline(bli_freqs{op_index}, bli_vals{op_index}, ...
            base_freqs{op_index});

        op_pt_error(op_index, 1) = find_error(abs(steam_vals_at_base_pts - base_vals{op_index}), ...
            abs(base_vals{op_index}), e_args.eps(), e_args.errorMode());
        op_pt_error(op_index, 2) = find_error(abs(bli_vals_at_base_pts - base_vals{op_index}), ...
            base_vals{op_index}, e_args.eps(), e_args.errorMode());

        steam_phase_estimation_error = unwrap(angle(steam_vals_at_base_pts), [], 1)- ...
            unwrap(angle(base_vals{op_index}), [], 1);
        steam_magnitude_estimation_error = abs(abs(steam_vals_at_base_pts)- ...
            abs(base_vals{op_index}))./abs(base_vals{op_index});

        bli_phase_estimation_error = unwrap(angle(bli_vals_at_base_pts), [], 1)- ...
            unwrap(angle(base_vals{op_index}), [], 1);
        bli_magnitude_estimation_error = abs(abs(bli_vals_at_base_pts)- ...
            abs(base_vals{op_index}))./abs(base_vals{op_index});

        if (nargout < 2)
            newMagFig = figure; hold on;
            % Plotting overlaid phase/magnitude plots
            plot(steam_freqs{op_index}, abs(steam_vals{op_index})', 'LineWidth', 2.6, ...
                'Marker', 's');
            plot(bli_freqs{op_index}, abs(bli_vals{op_index})', 'LineWidth', 2.6, ...
                'Marker', 'd');
            plot(base_freqs{op_index}, abs(base_vals{op_index})', 'LineWidth', 2.6, ...
                'Marker', 'o');
            grid on;
            legend(steam_plot_op_names{:}, bli_plot_op_names{:}, base_plot_op_names{:});
            set(gca, 'FontSize', 32);
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
            xlabel('Frequency (Hz)', 'FontSize', 28);
            ylabel('Magnitude (dB)', 'FontSize', 28);

            newPhaseFig = figure; hold on;
            plot(steam_freqs{op_index}, angle(steam_vals{op_index}), 'LineWidth', ...
                2.6, 'Marker', 's');
            plot(bli_freqs{op_index}, angle(bli_vals{op_index}), 'LineWidth', ...
                2.6, 'Marker', 'd');
            plot(base_freqs{op_index}, angle(base_vals{op_index}), 'LineWidth', ...
                2.6, 'Marker', 'o');
            grid on;
            legend(steam_plot_op_names{:}, bli_plot_op_names{:}, base_plot_op_names{:});
            set(gca, 'FontSize', 32);
            set(gca, 'XScale', 'log');
            xlabel('Frequency (Hz)', 'FontSize', 28);
            ylabel('Phase (radians)', 'FontSize', 28);
            
            newErrMagFig = figure; hold on;
            plot(steam_freqs{op_index}, steam_phase_estimation_error, ...
                'LineWidth', 2.6, 'Marker', 's');
            plot(bli_freqs{op_index}, bli_phase_estimation_error, ...
                'LineWidth', 2.6, 'Marker', 'd');
            grid on;
            legend(steam_error_plot_op_names{:}, bli_error_plot_op_names{:});
            set(gca, 'FontSize', 32);
            set(gca, 'XScale', 'log');
            xlabel('Frequency (Hz)', 'FontSize', 28);
            ylabel('\Delta{Phase} (radians)', 'FontSize', 28);

            newErrPhaseFig = figure; hold on;
            plot(steam_freqs{op_index}, steam_magnitude_estimation_error, ...
                'LineWidth', 2.6, 'Marker', 's');
            plot(bli_freqs{op_index}, bli_magnitude_estimation_error, ...
                'LineWidth', 2.6, 'Marker', 'd');
            grid on;
            legend(steam_error_plot_op_names{:}, bli_error_plot_op_names{:});
            set(gca, 'FontSize', 32);
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
            xlabel('Frequency (Hz)', 'FontSize', 28);
            ylabel('Relative Error in Magnitude (\Delta{Mag}/Mag)', 'FontSize', 28);
        end
    end

    estimation_error.STEAM = mean(op_pt_error(:,1));
    estimation_error.BLI = mean(op_pt_error(:,2));

    baseSolution = base_vals;
end
