function [dae, outputs, sim_args] = daeMOSGilbertCell(varargin)
% function [dae, outputs, sim_args] = daeMOSGilbertCell(model, parm_string, ...
%     method, subCkt)
% A make-shift Gilber cell circuit that acts as a mixer.
% Author: Archit Gupta (April 04, 2017)
% Motivation: Testing 2-tone (and later 3-tone) Harmonic balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP, INSTRUCTIONS [TODO]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cktdata.cktname = 'MOS Gilbert Cell';
    do_tabulate = 0;
    if (nargin > 3)
        do_tabulate = 1;
        cktdata.cktname = 'STEAM MOS Gilbert Cell';
    end

    cktdata.nodenames = { 'vdd', 'vss', 'S', 'v__lo_plus', 'v__lo_minus', 'v__ps', 'v__ns', 'v__s_plus', ...
        'v__s_minus', 'v__o_plus', 'v__o_minus' };
    cktdata.groundnodename = 'gnd';

    % Parameter declarations for the circuit
    RSS = 1e7;
    RL = 8e3;
    RR = 8e3;

    % VDD Element
    vdd_elem.name = 'Vdd';
    vdd_elem.model = vsrcModSpec('Vdd');
    vdd_elem.nodes = {'vdd', 'gnd'};
    vdd_elem.parms = {};

    % VSS Element
    vss_elem.name = 'vss';
    vss_elem.model = vsrcModSpec('VSS');
    vss_elem.nodes = {'vss', 'gnd'};
    vss_elem.parms = {};

    % VLO Element(s) - Local oscillator
    vlo_plus_elem.name = 'vlo_p';
    vlo_plus_elem.model = vsrcModSpec('VLOP');
    vlo_plus_elem.nodes = {'v__lo_plus', 'gnd'};
    vlo_plus_elem.parms = {};

    vlo_minus_elem.name = 'vlo_m';
    vlo_minus_elem.model = vsrcModSpec('VLOM');
    vlo_minus_elem.nodes = {'v__lo_minus', 'gnd'};
    vlo_minus_elem.parms = {};

    % Vs - The signal to be mixed
    vin_elem.name = 'vin';
    vin_elem.model = vsrcModSpec('VSP');
    vin_elem.nodes = {'v__s_plus', 'gnd'};
    vin_elem.parms = {};

    vin_minus_elem.name = 'vin_m';
    vin_minus_elem.model = vsrcModSpec('VSM');
    vin_minus_elem.nodes = {'v__s_minus', 'gnd'};
    vin_minus_elem.parms = {};

    % ISS - Current source
    iss_elem.name = 'IS';
    iss_elem.model = isrcModSpec('IS');
    iss_elem.nodes = {'S', 'vss'};
    iss_elem.parms = {};

    % Output resistances and capaticances
    rss_elem.name = 'RSS'; 
    rss_elem.model = resModSpec('RSS');
    rss_elem.nodes = {'S', 'vss'};
    rss_elem.parms = {RSS};

    rl_elem.name =  'RL';
    rl_elem.model = resModSpec('RL');
    rl_elem.nodes = {'vdd', 'v__o_plus'};
    rl_elem.parms = {RL};

    rr_elem.name = 'RR';
    rr_elem.model = resModSpec('RR');
    rr_elem.nodes = {'vdd', 'v__o_minus'};
    rr_elem.parms = {RR};

    % Filtering Circuit for down-conversion

    % This combination leaves a lot of ripples
    %{
    RFILTER = 2e6;
    CFILTER = 25e-12;    % Cut-off frequency: 20kHz
    %}
    RFILTER = 2e6;
    CFILTER = 50e-12;    % Cut-off frequency: 20kHz
    r_filter_elem.name = 'RF';
    r_filter_elem.model = resModSpec('RF');
    r_filter_elem.nodes = {'v__o_plus', 'v__o_minus'};
    r_filter_elem.parms = {RFILTER};
    c_filter_elem.name = 'CF';
    c_filter_elem.model = capModSpec('CF');
    c_filter_elem.nodes = {'v__o_plus', 'v__o_minus'};
    c_filter_elem.parms = {CFILTER};

    % Adding the voltage/current sources and resistances/capacitances to the list of circuit elements
    cktdata.elements = {vdd_elem, vss_elem, vin_elem, vin_minus_elem, vlo_plus_elem, vlo_minus_elem, iss_elem, ...
        rss_elem, rl_elem, rr_elem, r_filter_elem, c_filter_elem};

    % Adding diffpair subcircuits - Here is a list of nodenames for reference
    % { 'vdd', 'vss', 'S', 'v__lo_plus', 'v__lo_minus', 'v__ps', 'v__ns', 'v__s_plus', ...
    %     'v__s_minus', 'v__o_plus', 'v__o_minus' };

    % First Stage
    cktdata = add_subcircuit(cktdata, MOSDiffpair(varargin{:}), 'DP_STG_01', ...
        {'v__ps', 'v__ns', 'v__s_plus', 'v__s_minus', 'S'});

    % Second Stage: Positive Leg
    cktdata = add_subcircuit(cktdata, MOSDiffpair(varargin{:}), 'DPP_STG_02', ...
        {'v__o_plus', 'v__o_minus', 'v__lo_plus', 'v__lo_minus', 'v__ps'});

    % Second Stage: Negative Leg
    cktdata = add_subcircuit(cktdata, MOSDiffpair(varargin{:}), 'DPN_STG_02', ...
        {'v__o_minus', 'v__o_plus', 'v__lo_plus', 'v__lo_minus', 'v__ns'});

    cktdata = add_output(cktdata, 'e(v__o_plus)', 'e(v__o_minus)', 'vout');
    %{
    cktdata = add_output(cktdata, 'e(v__s_plus)', 'e(v__s_minus)', 'vin');
    cktdata = add_output(cktdata, 'e(v__lo_plus)', 'e(v__lo_minus)', 'vlo');
    cktdata = add_output(cktdata, 'e(v__o_plus)');
    cktdata = add_output(cktdata, 'e(v__ps)');
    cktdata = add_output(cktdata, 'e(v__ns)');
    %}
    cktdata = add_output(cktdata, 'e(S)');

    if (strfind(varargin{2}, 'PSP'))
        VCM = 0.4;
        VDD = 1;
    else
        VCM = 0.5;
        VDD = 1.5;
    end
    VSS = 0.0;
    v_lo_offset = 2*VCM;
    vin_offset = VCM;

    p_lo_args.f = 1e6;
    p_lo_args.A = 4e-1;
    p_lo_args.Offset = v_lo_offset + 6.557406991565868e-2;
    m_lo_args.f = p_lo_args.f;
    m_lo_args.A = -p_lo_args.A;
    m_lo_args.Offset = v_lo_offset - 9.571167857418955e-2;

    % Single Tone Function
    p_signal_args.fs = 500000; % Determines the number of transient steps
    p_signal_args.fc = p_lo_args.f;
    p_signal_args.A = 0;
    p_signal_args.Offset = vin_offset + 6.557406991565868e-2 ; % This is our signal (Random Numbers)
    m_signal_args.fs = 500000;
    m_signal_args.fc = p_lo_args.f;
    m_signal_args.A = 0;
    m_signal_args.Offset = vin_offset - 9.571167857418955e-2;

    % Two Tone Function
    %{
    p_signal_args.fs = 5000;
    p_signal_args.A = 0.5e-4;
    p_signal_args.Offset = vin_offset;
    m_signal_args.fs = 5000;
    m_signal_args.A = -0.5e-4;
    m_signal_args.Offset = vin_offset;

    % Using a modulation freuency
    p_signal_args.fc = p_lo_args.f;
    m_signal_args.fc = m_lo_args.f;
    %}

    % Transient Arguments
    sim_args.tstart = 0;
    sim_args.tstep = 0.02/p_lo_args.f;
    sim_args.tstop = sim_args.tstart + 1/p_signal_args.fs;

    % Finding if the model is BSIM or MVS
    if (findstr(varargin{2}, 'BSIM'))
        IS = 1.6e-6;  
        simAmp = 0.1;

        % DC Sweep Arguments
        sim_args.v_start = vin_offset + VDD/4.0;
        sim_args.v_step = -0.01;
        sim_args.v_stop = vin_offset - VDD/4.0;

        % AC Analysis Arguments
        sim_args.VIN_MIN = 0.4;
        sim_args.VIN_MAX = 0.6;
        sim_args.n_ops = 4;
    elseif (findstr(varargin{2}, 'PSP')) 
        IS = 1.6e-6;  
        simAmp = 0.1;

        % DC Sweep Arguments
        sim_args.v_start = vin_offset - VDD/4.0;
        sim_args.v_step = 0.01;
        sim_args.v_stop = vin_offset + VDD/4.0;

        % AC Analysis Arguments
        sim_args.VIN_MIN = 0.4;
        sim_args.VIN_MAX = 0.6;
        sim_args.n_ops = 4;
    elseif (findstr(varargin{2}, 'MVS'))
        IS = 1e-6;  % If this doesn't work, please look at daeMOSDiffpair, there is a long rant there
        sinAmp = 0.05;

        % AC Analysis Arguments
        sim_args.VIN_MIN = 0.4;
        sim_args.VIN_MAX = 0.6;
        sim_args.n_ops = 2;
    else
        error('Circuit parameters not found for the given MOS model');
    end

    dae = MNA_EqnEngine(cktdata);

    vdd_args.Val = VDD;
    vss_args.Val = VSS;
    iss_args.Val = IS;

    % Transient function description
    sinFunc = @(t,args) args.Offset+args.A*sin(args.f*2*pi*t);
    signalFunc = @(t, args) args.Offset + args.A*sin(args.fc*2*pi*t).*sin(args.fs*2*pi*t);
    stepFunc = @(t,args) args.Offset+args.A*smoothstep((t-args.Onset)/...
        args.scale, args.smoothing);
    constFunc = @(t,a) a.Val;

    % Setting up transient inputs
    dae = dae.set_utransient('Vdd:::E', constFunc, vdd_args, dae);
    dae = dae.set_utransient('vss:::E', constFunc, vss_args, dae);
    dae = dae.set_utransient('IS:::I', constFunc, iss_args, dae);
    dae = dae.set_utransient('vin:::E', signalFunc, p_signal_args, dae); % Named differently for the convenience of QSS
    dae = dae.set_utransient('vin_m:::E', signalFunc, m_signal_args, dae);
    dae = dae.set_utransient('vlo_p:::E', sinFunc, p_lo_args, dae);
    dae = dae.set_utransient('vlo_m:::E', sinFunc, m_lo_args, dae);

    % Setting up inputs for HB
    vin_p_args.Val =  p_signal_args.Offset;
    vin_m_args.Val = m_signal_args.Offset;
    v_lo_m_args.Val = m_lo_args.Offset;

    sinFDFunc = @(f, args) [args.Offset, i*args.A/2, -i*args.A/2];

    % Running single-tone
    dae = dae.set_uHB('Vdd:::E', constFunc, vdd_args, dae); % DC - Works
    dae = dae.set_uHB('vss:::E', constFunc, vss_args, dae); % DC - Works
    dae = dae.set_uHB('IS:::I', constFunc, iss_args, dae); % DC - Works
    dae = dae.set_uHB('vin:::E', constFunc, vin_p_args, dae); 
    dae = dae.set_uHB('vin_m:::E', constFunc, vin_m_args, dae);
    dae = dae.set_uHB('vlo_p:::E', sinFDFunc, p_lo_args, dae);
    dae = dae.set_uHB('vlo_m:::E', constFunc, v_lo_m_args, dae);

    % Running multi-tone

    % Setting up DC inputs
    QSS = dae.utransient(sim_args.tstart, dae);
    dae = dae.set_uQSS(QSS, dae);

    if (nargout > 2)
        % Checking if we have already computed an initial guess
        sim_args.daeIdentifier = [varargin{2}, '_gilbert_cell'];
        global STEAM_DATA_DIR;
        base_initguess_filename = [sim_args.daeIdentifier, '_initguess.mat'];
        initguess_filename = [STEAM_DATA_DIR, base_initguess_filename];

        if (exist(initguess_filename, 'file'))
            fprintf(2, 'Found initial guess for DAE: %s in %s\n', sim_args.daeIdentifier, base_initguess_filename);
            load(initguess_filename, '-mat');
        else
            initguess = zeros(dae.nunks(dae), 1);
            % Generating an initial guess using voltage stepping
            n_vdd_steps = 5;
            for vdd_val = [ linspace( 0.1, VDD, n_vdd_steps)];
                stepping_dae = dae.set_uQSS('Vdd:::E', vdd_val, dae); 
                % This is for DC Analysis only. Comment out for Transient
                stepping_dae = stepping_dae.set_uQSS('vin:::E', ...
                    sim_args.v_start, stepping_dae); 
                op_pt = op(stepping_dae, initguess);
                initguess = op_pt.getSolution(op_pt);
            end
            fprintf(2, 'Saving generated initial-guess in %s for later use\n', base_initguess_filename);
            save(initguess_filename, 'initguess', '-v7');
        end
        sim_args.xinit = initguess;
    end

    outputs = StateOutputs(dae);
    %{
    outputs = outputs.DeleteAll(outputs);
    outputs = outputs.Add({'vout'}, outputs); % For AC Analysis
    outputs = outputs.Add({'e_v__o_plus', 'e_v__lo_plus', 'e_v__s_plus', 'e_S', 'e_v__ps'}, outputs); % For DC (also helpful for TRAN)
    %}
end
