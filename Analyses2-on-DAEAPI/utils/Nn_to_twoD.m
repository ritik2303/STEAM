function X = Nn_to_twoD(vec_Nn, n, N)
%function X = Nn_to_twoD(vec_Nn, n, N)
% reshapes a vector in Nn format to twoD format (ie, a matrix with n rows and N cols)
% - the first row of X are the first N entries of vec_nN
% - the second row of X are the next n entries of vec_nN
% - and so on.
% 
% Recall:
% 1. nN format: N block-vectors (of F coeffs in standard order, or uniformly space time 
%    samples) stacked one after the other. Each block-vector is of size 
%    n (in the order of the DAE's unknowns)
% 2. Nn format: n block-vectors (in the order of the DAE's unknowns) stacked one after the other. 
%    Each block-vector is of size N (the F coeffs in standard order, or uniformly
%    spaced time samples) of the corresponding DAE unknown.
%
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/04
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2012 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	X = reshape(vec_Nn, N, n);
	X = X.';
	%%below: the long version, works fine
	%for i=1:n
	%	startidx_Nn = (i-1)*N+1;
	%	endidx_Nn = i*N;
	%	X(i,:) = vec_Nn(startidx_Nn:endidx_Nn).';
	%end
end
