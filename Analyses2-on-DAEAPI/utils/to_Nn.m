function out = to_Nn(longvec_nN, n, N)
%function out = to_nN(longvec_Nn, n, N)
% rearranges a vector in nN format to Nn format.
% Recall:
% 1. nN format: N block-vectors (of F coeffs in standard order, or uniformly space time 
%    samples) stacked one after the other. Each block-vector is of size 
%    n (in the order of the DAE's unknowns)
% 2. Nn format: n block-vectors (in the order of the DAE's unknowns) stacked one after the other. 
%    Each block-vector is of size N (the F coeffs in standard order, or uniformly
%    spaced time samples) of the corresponding DAE unknown.

	% two reshapes are probably more efficient than a for loop
	X = reshape(longvec_nN, n, N);
	out = reshape(X.',n*N,1);
	%%below: the long version, works fine
	%startindexes_Nn = (0:(n-1))*N + 1; % of size size n
	%for j=1:N % jth harmonic or timepoints
	%	startidx_nN = (j-1)*n+1;
	%	endidx_nN = j*n;
	%	vecj = longvec_nN(startidx_nN:endidx_nN); % jth harmonics/timepoints of all DAE unknowns: size n vector
	%	indices_Nn = startindexes_Nn + (j-1);
	%	out(indices_Nn,1) = vecj;
	%end
end
