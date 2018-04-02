function out = make_periodic(ts, tddata)
	% expects:
	% 	tddata.T: period
	%	tddata.vals: nxN matrix of TD data (N is the # of tpts)
	%	tddata.tpts: length N vector of timepoints - in increasing
	%		order. last-first must be <= T
	% returns: a T-periodic version of the waveform, with linear
	% 	interpolation via interp1
	tpts = tddata.tpts;
	tpts = tpts - tpts(1); % start at 0
	N = size(tpts);
	NN = size(tddata.vals,2);
	if N ~= NN
		error('make_periodic: tpts and vals have different sizes');
	end
	if tpts(N) > tddata.T
		error('make_periodic: tpts(end) > T');
	end

	ts = mod(ts, tddata.T);
	outt = interp1(tpts, tddata.vals.', ts, 'linear', 'extrap');
	out = outt.';
end
% end of make_periodic
