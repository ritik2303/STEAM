function [vals_FD_twoD, fharms, unif_vals_TD_twoD, unif_tpts] = FD_from_periodic_TD(tpts, vals_TD_twoD, M)
% compute Fourier components from non-uniformly-spaced periodic TD data for 1 cycle
%
% input args:
% 	tpts: 		time points at which vals is sampled. 
%			Note: the period T is taken to be tpts(end)-tpts(1).
%	vals_TD_twoD: 	samples of the data in twoD format (n x ntpts). Note: periodicity of the data is NOT CHECKED.
%	M:    		number of positive harmonics desired for the Fourier analysis. 
%	      		(The total number of F. coeffs will be N=2*M+1.)
% outputs:
%	vals_FD_twoD: 		the Fourier coeffs in twoD format (n x N) in standard FFT order (0,1,...,M,-M,...,-1)
%	fharms: 		the array [0, 1, ..., M, -M, ..., -1]
%	unif_vals_TD_twoD: 	vals_TD_twoD resampled on the uniform grid tpts(1) + (0:N-1)/N*T
%	unif_tpts:		the array tpts(1) + (0:N-1)/N*T.
%	
% Note: uses interp1 to resample onto a uniform grid, so numerical dynamic range poor
% Author: Jaijeet Roychowdhury <jr@berkeley.edu>, 2012/06/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2012 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	N = 2*M+1;
	fharms = [0:M, -M:-1];
	T = tpts(end) - tpts(1);
	unif_tpts = tpts(1) + (0:N-1)/N*T;

	% compute Fourier components of the eqn PPVs
	unif_vals_TD_twoD = interp1(tpts, vals_TD_twoD.', unif_tpts.').';

	vals_FD_twoD = 1/N*fft(unif_vals_TD_twoD.').';
end 
% end function FD_from_TD
