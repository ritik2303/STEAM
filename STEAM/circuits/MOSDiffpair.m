function subcktdata = MOSDiffpair(model, parm_string, subCkt, method)
% function subcktdata = MOSDiffpair(model, parm_string, subCkt, method)
    subckt.cktname = 'MOS Diffpair';
    do_tabulate = 0;

    if (nargin > 3)
        do_tabulate = 1;
        subcktdata.cktname =  'STEAM MOS Diffpair';
    end

    nmos_parm_string = strcat(parm_string, '_NMOS');

    % node (names)
    subcktdata.nodenames = {'v__o_plus', 'v__o_minus', 'v__i_plus', 'v__i_minus', 'v__s'};
    subcktdata.terminalnames = subcktdata.nodenames;
    subcktdata.groundnodename = 'gnd';
    subcktdata.elements = {};

    % list of elements 
    if (do_tabulate)
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string, method), 'P', {'v__o_plus', 'v__i_plus', 'v__s'});
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string, method), 'N', {'v__o_minus', 'v__i_minus', 'v__s'});
    else
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string), 'P', {'v__o_plus', 'v__i_plus', 'v__s'});
        subcktdata = add_subcircuit(subcktdata, subCkt(model, ...
            nmos_parm_string), 'N', {'v__o_minus', 'v__i_minus', 'v__s'});
    end

