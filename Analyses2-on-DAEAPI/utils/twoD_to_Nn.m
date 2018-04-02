function out = twoD_to_Nn(X)
%function out = twoD_to_Nn(X)
% reshapes twoD format (ie, a matrix with n rows and N cols) to a vector in Nn format
% - the first row of X becomes the first n entries of out
% - the second row of X become the next n entries of out
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
	out = reshape(X.',n*N,1);
	%%below: the long version, works fine
	%for i=1:n
	%	startidx_Nn = (i-1)*N+1;
	%%	endidx_Nn = i*N;
	%	out(startidx_Nn:endidx_Nn,1) = X(i,:).';
	%end
end
