function [dae, outputs, sim_args] = daeMOSDiffpair(varargin)
% function [dae, outputs, sim_args] = daeMOSDiffpair(varargin)
% An ideal-ish (I don't know what this means. Please see help MVSdiffpair_ckt to
% know more). The diffpair used in this circuit is not symmetric w.r.t the
% input. 
% > There are 2 MOSFETS (both N-Type). MOSFET-L has Vin connected to its
% gate, whereas the gate for MOSFET-R has been grounded. 
% > The sources of these FETs are connected to a node nS, which is connected
% to an Ideal current source


    cktdata.cktname = 'MOS diffpair';
    do_tabulate = 0;
    if (nargin > 3)
        do_tabulate = 1;
        cktdata.cktname = 'STEAM MOS diffpair';
    end    

    %=========================================================================%
    %   Setting up the circuit Netlist
    %=========================================================================%
    
    cktdata.nodenames = {'vdd', 'DL', 'DR', 'S', 'in', 'CM'};
    cktdata.groundnodename = 'gnd';

    % Circuit Parameter Values
    
    RLe = 50000;
    RRe = 50000;
    RSS = 100000;
    CLe = 1e-6;
    CRe = 1e-6;

    %vdd_elem
    vdd_elem.name = 'vdd';
    vdd_elem.model = vsrcModSpec('vdd'); % get the ModSpec model
    vdd_elem.nodes = {'vdd', 'gnd'}; % positive and negative nodes of the vsrc
    vdd_elem.parms = {};

    cktdata.elements = {vdd_elem};

    %vcm_elem -- For setting the bias point
    vcm_elem.name = 'vcm';
    vcm_elem.model = vsrcModSpec('VCM'); % get the ModSpec model
    vcm_elem.nodes = {'CM', 'gnd'}; % positive and negative nodes of the vsrc
    vcm_elem.parms = {};

    cktdata.elements = {cktdata.elements{:}, vcm_elem};

    %vin_elem
    vin_elem.name = 'vin';
    vin_elem.model = vsrcModSpec('vin'); % get the ModSpec model
    vin_elem.nodes = {'in', 'gnd'}; % positive and negative nodes of the vsrc
    vin_elem.parms = {};

    cktdata.elements = {cktdata.elements{:}, vin_elem};

    %iss_elem
    iss_elem.name = 'IS';
    iss_elem.model = isrcModSpec('IS'); % get the ModSpec model
    iss_elem.nodes = {'S', 'gnd'}; % positive and negative nodes of the ssrc
    iss_elem.parms = {};

    cktdata.elements = {cktdata.elements{:}, iss_elem};

    %res_ss_elem
    res_ss_elem.name = 'RSS';
    res_ss_elem.model = resModSpec('RSS'); % get the ModSpec model
    res_ss_elem.nodes = {'S', 'gnd'}; % resistor nodes
    res_ss_elem.parms = {RSS}; % or use feval(res_l_elem.model.defaultparms,

    cktdata.elements = {cktdata.elements{:}, res_ss_elem};

    %res_l_elem
    res_l_elem.name = 'R1';
    res_l_elem.model = resModSpec('R1'); % get the ModSpec model
    res_l_elem.nodes = {'vdd', 'DL'}; % resistor nodes
    res_l_elem.parms = {RLe}; % or use feval(res_l_elem.model.defaultparms,

    cktdata.elements = {cktdata.elements{:}, res_l_elem};

    %res_r_elem
    res_r_elem.name = 'R2';
    res_r_elem.model = resModSpec('R2'); % get the ModSpec model
    res_r_elem.nodes = {'vdd', 'DR'}; % resistor nodes
    res_r_elem.parms = {RRe}; % or use feval(res_r_elem.model.defaultparms,

    cktdata.elements = {cktdata.elements{:}, res_r_elem};

    cktdata = add_subcircuit(cktdata, MOSDiffpair(varargin{:}), 'DP', ...
        {'DL', 'DR', 'in', 'CM', 'S'});
    cktdata = add_output(cktdata, 'e(DR)', 'e(DL)', 'vout');
    cktdata = add_output(cktdata, 'e(S)');
    cktdata = add_output(cktdata, 'e(in)', 'e(CM)', 'vin');

    dae = MNA_EqnEngine(cktdata);

    % Setting up transient/DC/AC Parameters

    VDD = 1.0;
    VCM = 0.7;
    IS = 0;
    VInOffset = 0;
    sinAmp = 0;

    % DC Sweep Arguments 
    sim_args.v_start = 0;
    sim_args.v_step = 0.01;
    sim_args.v_stop = VDD;

    % AC Analysis Arguments
    sim_args.VIN_MIN = VCM-0.1;
    sim_args.VIN_MAX = VCM+0.1;
    sim_args.n_ops = 1;

    simArgs.fstart = 1;
    simArgs.fstop = 1e9;
    simArgs.nsteps = 5; % 5 points per decade

    % Finding if the model is BSIM or MVS
    if (findstr(varargin{2}, 'BSIM'))
        IS = 2.4e-5;  
        VInOffset = VCM;
        simAmp = 0.1;

        % Transient Arguments
        sim_args.tstart = 0;
        sim_args.tstep = 6e-5;
        sim_args.tstop = 8e-3;
    elseif(findstr(varargin{2}, 'PSP'))
        IS = 8e-6;  
        VInOffset = VCM;
        simAmp = 0.1;

        % Transient Arguments
        sim_args.tstart = 0;
        sim_args.tstep = 6e-5;
        sim_args.tstop = 8e-3;
    elseif(findstr(varargin{2}, 'MVS'))
        IS = 1e-6;
        VInOffset = 0.4;
        sinAmp = 0.05;

        % Transient Arguments
        sim_args.tstart = 0;
        sim_args.tstep = 6e-5;
        sim_args.tstop = 8e-3;

        % AC Analysis Arguments
        sim_args.VIN_MIN = 0.4;
        sim_args.VIN_MAX = 0.6;
        sim_args.n_ops = 2;
    else
        fprintf(2, 'Circuit parameters not found for the given MOS model\n');
    end

    % Different Transient Inputs to try
    sin_args.f = 500;
    sin_args.A = 0.6;
    sin_args.Offset = VInOffset;

    step_args.A = 1;
    step_args.scale = 10*sim_args.tstep;
    step_args.Offset = 0;
    step_args.Onset = (sim_args.tstart+sim_args.tstop)/2;
    step_args.smoothing = 0.01;

    vdd_args.Val = VDD;
    is_args.Val = IS;
    vcm_args.Val = VCM;

    % Transient function description
    sinFunc = @(t,args) args.Offset+args.A*sin(args.f*2*pi*t);
    stepFunc = @(t,args) args.Offset+args.A*smoothstep((t-args.Onset)/...
        args.scale, args.smoothing);
    constFunc = @(t,a) a.Val;
    
    inFunc = sinFunc; in_args = sin_args;
    %inFunc = stepFunc; in_args = step_args;

    % Setting up transient inputs
    dae = dae.set_utransient('vin:::E', inFunc, in_args, dae);
    dae = dae.set_utransient('vcm:::E', constFunc, vcm_args, dae);
    dae = dae.set_utransient('vdd:::E', constFunc, vdd_args, dae);
    dae = dae.set_utransient('IS:::I', constFunc, is_args, dae);

    % Setting up DC inputs
    QSS = dae.utransient(sim_args.tstart, dae);
    dae = dae.set_uQSS(QSS, dae);

    if (nargout > 2)
        % Checking if we have already computed an initial guess
        sim_args.daeIdentifier = [varargin{2}, '_diffpair'];
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
    %{
    outputs = outputs.DeleteAll(outputs);
    %outputs = outputs.Add({'e_DL', 'e_DR', 'e_S', 'e_in', 'e_CM'}, outputs);
    %outputs = outputs.Add({'e_DL', 'e_DR', 'e_S'}, outputs); % For AC Analysis
    %}
end 
