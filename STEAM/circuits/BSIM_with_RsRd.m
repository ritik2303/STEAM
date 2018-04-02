function subcktnetlist = BSIM_with_RsRd(model, parm_string, method)
% function subcktnetlist = BSIM_with_RsRd(model, parm_string, method)
% Author: Archit Gupta
% This script takes a BSIM MOS model, tabulates it and returns a MOS model that
% has only 3  terminals. However, in order to compensate for the terminal
% resistances, it adds resistances Rs and Rd at the Drain and Source of the
% MOSFET. 
    
    fprintf(2, 'Creating BSIM with terminal resistances ...\n');
    pbaseStr = '';
    index = findstr(parm_string, 'VA');
    if (length(index) ~= 0)
        pbaseStr = 'parm_';
    end

    prshStr = strcat(pbaseStr,'rsh');
    pnrsStr = strcat(pbaseStr,'nrd');
    pnrdStr = strcat(pbaseStr,'nrs');

    subcktnetlist.cktname = 'BSIM_with_RsRd';
    subcktnetlist.nodenames = {'d', 'g', 's', 'di', 'si'};
    subcktnetlist.groundnodename = 'gnd';
    
    % Adding the core MOS Model
    MOSElem.name = parm_string;
    if (nargin > 2)
        MOSElem.model = STEAM( MOS5TModel2D(model, parm_string), method.i_method, method.d_method, ...
            method.order, method.bounds );
    else
        MOSElem.model = MOS5TModel2D(model, parm_string);
    end
    % [TODO] Compression Factor has been ignored here
    MOSElem.nodes = {'di', 'g', 'si'};
    MOSElem.parms = MOSElem.model.getparms();

    subcktnetlist.elements = {MOSElem};

    sheetResistance = MOSElem.model.getparms(prshStr, MOSElem.model);
    numSquaresDrain = MOSElem.model.getparms(pnrdStr, MOSElem.model);
    numSquaresSource = MOSElem.model.getparms(pnrsStr, MOSElem.model);

    RsVal = sheetResistance*numSquaresSource;
    RdVal = sheetResistance*numSquaresDrain;

    subcktnetlist = add_element(subcktnetlist, resModSpec('Rs'), strcat( ...
        parm_string,'Rs'), {'s', 'si'}, {{'R', RsVal}});
    subcktnetlist = add_element(subcktnetlist, resModSpec('Rd'), strcat( ...
        parm_string,'Rd'), {'d', 'di'}, {{'R', RdVal}});

    subcktnetlist.terminalnames = {'d', 'g', 's'};
end
