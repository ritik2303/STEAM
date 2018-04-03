function args = BSIMSamplingArgs()
% function args = BSIMSamplingArgs()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initializes a struct with the default sampling arguments for compressive
% sensing experiments.
%
% See: loadBSIMSlice, loadBSIMImage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    args.l = 1;
    args.m = 1;
    args.vds_offset = 2;
    args.vgs_offset = 1;
    args.VMIN = -1.6;
    args.VMAX = 1.6;

end