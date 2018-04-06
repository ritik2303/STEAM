function [dae, outputs, sim_args] = daeMOSSourceFollower(model, parm_string, ...
    subCkt, method)
% function [dae, outputs, sim_args] = daeMOSSourceFollower(model, parm_string, ...
%     subCkt, method)
% Author: Archit Gupta
% Date: April 06, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAE for a source follower circuit. Simple and straightforward source follower
% used to demonstrate the harmonic discontinuity on running Harmonic Balance with
% a crude model like Shichman Hodges.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cktdata.cktname = 'Source-Follower';
    cktdata.nodenames = {'vdd', 'd', 'g', 's'};
    cktdata.groundnodename = 'gnd';

    % Circuit parameters (just the resistance being used) 
    Rs  = 1.0e3;

    vinElem.name    = 'vin';
    vinElem.model   = vsrcModSpec('vin');
    vinElem.nodes   = {'g', 'gnd'};
    vinElem.parms   = {};

    vddElem.name    = 'vdd';
    vddElem.model   = vsrcModSpec('vdd');
    vddElem.nodes   = {'vdd', 'd'};
    vddElem.parms   = {};

    rsElem.name     = 'rs';
    rsElem.model    = resModSpec('rs');
    rsElem.nodes    = {'s', 'gnd'};
    rsElem.parms    = {Rs};

    cktdata.elements = {vinElem, vddElem, rsElem};

    mos_parm_string = strcat(parm_string, '_NMOS');
    if (nargin > 3)
        mos_sub_ckt     = subCkt(model, mos_parm_string, method);
    else
        mos_sub_ckt     = subCkt(model, mos_parm_string);
    end
    cktdata         = add_subcircuit(cktdata, mos_sub_ckt, 'MOS', {'d', 'g', 's'});
    cktdata         = add_output(cktdata, 'e(d)');
    cktdata         = add_output(cktdata, 'e(g)');
    cktdata         = add_output(cktdata, 'e(s)');

    dae             = MNA_EqnEngine(cktdata);

    % Create the simulation arguments
    VDD = 1.0;
    VInOffset = 0.5;
    sinAmp = 1e-3;

    % DC Sweep Arguments 
    sim_args.v_start = 0;
    sim_args.v_step = 0.02;
    sim_args.v_stop = VDD;

    % AC Analysis Arguments
    sim_args.VIN_MIN = VInOffset-0.1;
    sim_args.VIN_MAX = VInOffset+0.1;
    sim_args.n_ops = 1;

    simArgs.fstart = 1;
    simArgs.fstop = 1e9;
    simArgs.nsteps = 5; % 5 points per decade

    % Transient analysis Arguments
    sim_args.tstart = 0.0;
    sim_args.tstep  = 1.0e-5;
    sim_args.tstop  = 2.0e-3;

    % Harmonic Balance arguments
    % Same as the sin_args, the function's parameters can be reused

    % Different Transient Inputs to try
    vdd_args.Val    = VDD;
    sin_args.f      = 500;
    sin_args.A      = 0.6;
    sin_args.Offset = VInOffset;

    % Transient function description
    sinFunc     = @(t,args) args.Offset+args.A*sin(args.f*2*pi*t);
    constFunc   = @(t,a) a.Val;

    % Harmonic Balance function description
    sinFDFunc = @(f, args) [args.Offset, i*args.A/2, -i*args.A/2];

    inFunc      = sinFunc;
    in_args     = sin_args;

    % Principal frequency (needed for Harmonic Balance)
    sim_args.f0 = in_args.f;

    % Setting up function handles for transient simulation
    dae = dae.set_utransient('vin:::E', inFunc, in_args, dae);
    dae = dae.set_utransient('vdd:::E', constFunc, vdd_args, dae);

    % Setting up DC inputs
    QSS = dae.utransient(sim_args.tstart, dae);
    dae = dae.set_uQSS(QSS, dae);

    % Setting up Single-Tone Harmonic Balance inputs
    dae = dae.set_uHB('vdd:::E', constFunc, vdd_args, dae); % DC - Works
    dae = dae.set_uHB('vin:::E', sinFDFunc, sin_args, dae); 

    % TODO: Just make a separate function out of this, it is used in too many places.
    if (nargout > 2)
        % Checking if we have already computed an initial guess
        sim_args.daeIdentifier = [parm_string, '_source_follower'];
        global STEAM_DATA_DIR;
        base_initguess_filename = [sim_args.daeIdentifier, '_initguess.mat'];
        initguess_filename = [STEAM_DATA_DIR, base_initguess_filename];

        if (exist(initguess_filename, 'file'))
            fprintf(2, 'Found initial guess for DAE: %s in %s\n', sim_args.daeIdentifier, base_initguess_filename);
            load(initguess_filename);
        else
            initguess = randn(dae.nunks(dae), 1);
            % Generating an initial guess using voltage stepping
            n_vdd_steps = dae.nunks(dae);
            for vdd_val = [ linspace( 0.5*VDD, VDD, n_vdd_steps)];
                stepping_dae = dae.set_uQSS('vdd:::E', vdd_val, dae); 
                op_pt = op(stepping_dae, initguess);
                initguess = op_pt.getSolution(op_pt);
            end
            fprintf(2, 'Saving generated initial-guess in %s for later use\n', base_initguess_filename);
            save(initguess_filename, 'initguess');
        end
        sim_args.xinit = initguess;
    end

    outputs = StateOutputs(dae);