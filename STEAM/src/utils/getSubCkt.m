function subCkt = getSubCkt(parm_string)
% function subCkt = getSubCkt(parm_string)
% Author: Archit Gupta
% Date: April 04, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility function, returns a function handle for creating a sub-circuit out
% of a transistor model. Useful when a 5-Variable transistors has to be replaced
% with a 3-Variable transistor and 2 external resistances.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (strfind( parm_string, 'BSIM' ))
        subCkt = @BSIM_with_RsRd;
    elseif (strfind( parm_string, 'MVS' ))
        subCkt = @MVS_with_RsRd; 
    elseif (strfind( parm_string, 'PSP' ))
        subCkt = @PSP_subckt; 
    end
end