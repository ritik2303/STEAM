function X = nN_to_twoD(vec_nN, n, N)
%function X = nN_to_twoD(vec_nN, n, N)
% reshapes a vector in nN format to twoD format (ie, a matrix with n rows and N cols)
% - the first column of X are the first n entries of vec_nN
% - the second column of X are the next n entries of vec_nN
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
%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/01
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2012 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	X = reshape(vec_nN, n, N);
	%%below: the long version, works fine
	%for j=1:N
	%	startidx_nN = (j-1)*n+1;
	%	endidx_nN = j*n;
	%	X(:,j) = vec_nN(startidx_nN:endidx_nN);
	%end
end
