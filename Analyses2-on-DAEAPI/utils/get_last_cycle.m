function [tddata, harmonics] = get_last_cycle(tran, T, cktfuncs, doplot)
	% now yank one cycle worth of waveforms out of tran
	% the periodic input to HB is over the range [0,T]
	if nargin < 4
		doplot = 0;
	end
	transol = getsolution(tran);
        % npts: 1001; tpts: [1x1001 double]; vals: [2x1001 double]
	lastcycle = floor(transol.tpts(end)/T);
	startpt = T*(lastcycle-1); endpt = T*lastcycle;
	indices = find((transol.tpts >= startpt) & (transol.tpts < endpt));

	ts = transol.tpts(indices);
	vals = transol.vals(:,indices);
	[n, Npts] = size(vals);
	%plot(ts,vals(1,:), '.-')

	% horrible hack: not clear what the following means if Npts is even
	harmonics = fft(vals.')/Npts; % fft is applied to each column of arg
	harmonics = harmonics'; % now rows are DAE vars, cols are harmonics
	tddata.tpts = ts;
	tddata.vals = vals;
	tddata.T = T;
	if 1 == doplot
		% plot it
		varnames = feval(cktfuncs.mpde_variable_names);
		for i = 1:size(harmonics,1) % one plot per DAE variable
			figure();
			xs = [0:(Npts-1)];
			ys = 20*log10([abs(harmonics(i,1)), ...
					2*abs(harmonics(i,2:Npts))]);
			stem(xs,ys);
			xlabel('harmonic number');
			ylabel('"voltage db" [20*log10(mag)]');
			ttl = sprintf('%s (fundamental freq is %0.6e)', ...
				varnames{i}, 1/T);
			title(ttl); grid on; axis tight;
		end
	end
end
% end get_last_cycle
