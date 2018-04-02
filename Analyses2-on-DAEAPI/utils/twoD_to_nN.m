function out = twoD_to_nN(X)
%function out = twoD_to_nN(X)
% reshapes twoD format (ie, a matrix with n rows and N cols) to a vector in nN format
% - the first column of X becomes the first n entries of out
% - the second column of X become the next n entries of out
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

	[n, N] = size(X);
	out = reshape(X,n*N,1);
	%%below: the long version, works fine
	%for j=1:N
	%	startidx_nN = (j-1)*n+1;
	%	endidx_nN = j*n;
	%	out(startidx_nN:endidx_nN,1) = X(:,j);
	%end
end
