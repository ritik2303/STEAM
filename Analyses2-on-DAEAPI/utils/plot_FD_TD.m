function [figh, legends, colindex] = plot_FD_TD(titlestr, onames, frequency, fharms, out_FD_twoD, tpts, out_TD_twoD, ...
    time_units, figh, legends, colindex)
%function [figh, onames, colindex] = plot_FD_TD(titlestr, onames, frequency, fharms, out_FD_twoD, tpts, out_TD_twoD, ...
%    time_units, figh, legends, colindex)
%useful mostly for periodic waveforms - eg, in HB.plot
%
%Arguments:
%
%- titlestr: a name that goes into the title of the plots.
%- onames: cell array of output names, of length number of variables to be
%          plotted - onames can be a cell array of cell arrays, too,
%          corresponding to fharms/tpts/etc, below
%- frequency: fundamental frequency for the FD plots
%- fharms: array of harmonic numbers of the frequencies - of length N
%	       - fharms can be a cell array of such arrays, too. 
%- out_FD_twoD: the Fourier coeffs of the outputs in 2D matrix format (no x N)
%	            - out_FD_twoD can be a cell array of such matrices, too,
%	              corresponding to fharms.
%tpts: timepoints at which the time-domain data is specified - of length Nt
%	   - tpts can be a cell array of such arrays, too. 
%out_TD_twoD: time-domain data for the outputs in 2D matrix format (no x Nt)
%	          - out_TD_twoD can be a cell array of such matrices, too,
%	            corresponding to tpts.
%time_units: a string indicating the units of time. If not provided, assumed to
%            be 'sec'
%

%Author: Jaijeet Roychowdhury <jr@berkeley.edu> 2012/06/10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Roychowdhury.
% Copyright (C) 2008-2017 Jaijeet Roychowdhury <jr@berkeley.edu>. All rights
% reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if nargin < 8
		time_units = 'sec';
	end

	if 1 == isa(fharms, 'cell')
		n_data_sets = length(fharms);
	else
		n_data_sets = 1;
		fharms = {fharms};
		onames = {onames};
		out_FD_twoD = {out_FD_twoD};
		tpts = {tpts};
		out_TD_twoD = {out_TD_twoD};
	end

        if (nargin <  9) || (0 == sum(size(figh)))
            figh = cell(length(onames{1}),1);
        end

        if (nargin <  10) || (0 == sum(size(legends)))
            legends = cell(length(onames{1}),1);
        end

        if (nargin <  11) || (0 == sum(size(colindex)))
            colindex = cell(length(onames{1}),1);
        end

	if 0 == length(onames) || 0 == length(onames{1})
		n_o = size(out_FD_twoD{1},1);
	else
		n_o = length(onames{1});
	end
	for i=1:n_o
                if ( isa(figh{i}, 'matlab.ui.Figure') )
		    figure(figh{i}); hold on;
                else
                    figh{i} = figure;
                end

                if ( ~isa(legends{i}, 'cell') )
                    legends{i} = {}; 
                end

		subplot(2,1,1);
		hold on;
		for k=1:n_data_sets
			if 0 == length(onames) || 0 == length(onames{k})
				outputname = sprintf('output %d',i);
			else
				outnames = onames{k};
				outputname = outnames{i};
			end
			thiscol = getcolorfromindex(gca,length(legends{i})+k);
			FDdata = out_FD_twoD{k};
			y=20*log10([abs(FDdata(i,:))]);
			%stairs(x,y);
			%bar(x,y);
			stem(fharms{k},y,'color', thiscol, 'LineWidth', 2.0);
			legends{i} = {legends{i}{:}, escape_special_characters(outputname)};
		end % k
		xlabel('harmonic number');
		ylabel('"voltage db" [20*log10(mag)]');
		legend(legends{i});
                set(gca, 'FontSize', 28);
		if max(strcmp(time_units, {'s', 'sec', 'secs', 'second', 'seconds'})) > 0
			freq_units = 'Hz';
		else
			freq_units = sprintf('1/%s', time_units);
		end
		title(sprintf('%s\nFourier comp. magnitudes, freq=%0.6e (%s)', titlestr, frequency, freq_units));
		grid on; axis tight;

		subplot(2,1,2);
		hold on;
		for k=1:n_data_sets
			thiscol = getcolorfromindex(gca,length(legends{i})+k);
			TDdata = out_TD_twoD{k};
			plot(tpts{k}, TDdata(i,:),'.-', 'color', thiscol, 'LineWidth', 2.0);
		end % k=1:n_data_sets
		ylabel('value');
		xlabel(sprintf('time (%s)', time_units));
                set(gca, 'FontSize', 28);
		legend(legends{i});
		title(sprintf('%s\nTime-domain waveforms, period=%0.6e (%s)', titlestr, 1/frequency, time_units));
		grid on; axis tight;
	end
end
