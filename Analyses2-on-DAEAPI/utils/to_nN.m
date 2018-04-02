function out = to_nN(longvec_Nn, n, N)
%function out = to_nN(longvec_Nn, n, N)
% rearranges a vector in Nn format to nN format.
% Recall:
% 1. nN format: N block-vectors (of F coeffs in standard order, or uniformly space time 
%    samples) stacked one after the other. Each block-vector is of size 
%    n (in the order of the DAE's unknowns)
% 2. Nn format: n block-vectors (in the order of the DAE's unknowns) stacked one after the other. 
%    Each block-vector is of size N (the F coeffs in standard order, or uniformly
%    spaced time samples) of the corresponding DAE unknown.

	% two reshapes are probably more efficient than a for loop
	X = reshape(longvec_Nn, N, n);
	out = reshape(X.',n*N,1);
	%%below: the long version, works fine
	%startindexes_nN = (0:(N-1))*n + 1; % of size N
	%for j=1:n % j = DAE unk index
	%	startidx_Nn = (j-1)*N+1;
	%	endidx_Nn = j*N;
	%	vecj = longvec_Nn(startidx_Nn:endidx_Nn); % all harmonics/timepoints of jth DAE unknowns: size N vector
	%	indices_nN = startindexes_nN + (j-1);
	%	out(indices_nN,1) = vecj;
	%end
end
