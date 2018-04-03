function [dae, outputs, sim_args] = daeMOSInverterChain(n_elements, model, ...
    parm_string, subCkt, method)
% function [dae, outputs, sim_args] = daeMOSInverterChain(n_elements, model, ...
%     parm_string, subCkt, method)
% An 'N' stage inverter chain 
% Author: Archit Gupta, 2017/03/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   The Circuit:
%       An 'N' stage inverter chain. The inverters for this chain are
%       constructed using BSIM/MVS transistors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cktdata.cktname = ['MOS Inverter Chain (', num2str(n_elements), '-Stage)'];
    do_tabulate = 0;

    if (nargin > 4)
        do_tabulate = 1;
        cktdata.cktname = ['Tabulated', cktdata.cktname];
    end

    vdd_dc = 1.0;
    c_l = 1e-12;

    %===========================================================================
    % subcircuit: MOSinverter
    %---------------------------------------------------------------------------

    inverter = [];
    if (do_tabulate)
        inverter = MOSInverter(model, parm_string, subCkt, method);
    else
        inverter = MOSInverter(model, parm_string, subCkt);
    end

    inverter = add_element(inverter, capModSpec(), 'c_l', ...
       {'out', 'gnd'}, {{'C', c_l}}, {});

    %===========================================================================

    cktdata.groundnodename = 'gnd';

    cktdata.nodenames = cell(1, n_elements+2);
    cktdata.nodenames{1} = 'vdd';

    cktdata.nodenames{2} = 'in';
    cktdata.nodenames{end} = 'out';
    for n_inv = 1:n_elements-1  % Label the output terminals of all but the last
                                % inverters in the chain
        % Inverter i
        cktdata.nodenames{2+n_inv} = ['inv', num2str(n_inv)];

    end

    cktdata = add_element(cktdata, vsrcModSpec(), 'Vdd', ...
    {'vdd', 'gnd'}, {}, {{'E', {'DC', vdd_dc}}});
    cktdata = add_element(cktdata, vsrcModSpec(), 'vin', ...
    {'in', 'gnd'}, {}, {{'E', {'DC', 0}}});

    for n_inv = 1:n_elements
        cktdata = add_subcircuit(cktdata, inverter, ['X', num2str(n_inv)], ...
            {'vdd', cktdata.nodenames{1+n_inv}, cktdata.nodenames{2+n_inv}});
    end

    dae = MNA_EqnEngine(cktdata);
    outputs = StateOutputs(dae);
    outputs = outputs.DeleteAll(outputs);
    outputs = outputs.Add({'e_in', 'e_out'}, outputs);

    % Setting up simulation arguments for the device model / parameters
    
    sim_args.tstart = 0;
    n_dc_steps = 100;
    sim_args.v_start = 0;
    sim_args.v_stop = vdd_dc;
    sim_args.v_step = (sim_args.v_stop - sim_args.v_start)/n_dc_steps;

    v_in_offset = vdd_dc/2;
    sin_args.A = .2;
    sin_args.f = 1e4;
    sin_args.Offset = v_in_offset;

    % Transient Arguments
    sim_args.tstart = 0;
    sim_args.tstep = 0.01/sin_args.f;
    sim_args.tstop = 4/sin_args.f;

    if (findstr(parm_string, 'BSIM'))
    elseif (findstr(parm_string, 'MVS'))
    elseif (findstr(parm_string, 'PSP'))
    else
        error('ERROR: circuit parameters not found for the given MOS model');
    end

    sim_args.VDD = vdd_dc;
    vdd_args.VDD = vdd_dc;

    % Transient function descriptions 
    vddFunc = @(t, a) a.VDD;
    sinFunc = @(t,args) args.Offset+args.A*sin(args.f*2*pi*t);

    % Setting up the transient functions
    dae = dae.set_utransient('Vdd:::E', vddFunc, vdd_args, dae);
    dae = dae.set_utransient('vin:::E', sinFunc, sin_args, dae);

    % Setting up DC inputs
    QSS = dae.utransient(sim_args.tstart, dae);
    dae = dae.set_uQSS(QSS, dae);

    % Setting up inputs for HB
    dae = dae.set_uHB('Vdd:::E', vddFunc, vdd_args, dae);

    if (nargout > 2)
        % Checking if we have already computed an initial guess
        sim_args.daeIdentifier = [num2str(n_elements), '_', parm_string, '_inverter_chain'];
        global STEAM_DATA_DIR;
        base_initguess_filename = [sim_args.daeIdentifier, '_initguess.mat'];
        initguess_filename = [STEAM_DATA_DIR, base_initguess_filename];

        if (exist(initguess_filename, 'file'))
            fprintf(2, 'Found initial guess for DAE: %s in %s\n', sim_args.daeIdentifier, base_initguess_filename);
            load(initguess_filename, '-mat');
        else
            initguess = zeros(dae.nunks(dae), 1);
            % Generating an initial guess using voltage stepping
            n_vdd_steps = dae.nunks(dae);
            for vdd_val = [ linspace( 0.1, vdd_dc, n_vdd_steps)];
                stepping_dae = dae.set_uQSS('Vdd:::E', vdd_val, dae); 
                op_pt = op(stepping_dae, initguess);
                initguess = op_pt.getSolution(op_pt);
            end
            fprintf(2, 'Saving generated initial-guess in %s for later use\n', base_initguess_filename);
            save(initguess_filename, 'initguess', '-v7');
        end
        sim_args.xinit = initguess;
    end
end
