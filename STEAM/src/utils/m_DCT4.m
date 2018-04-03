function [x_dct, T] = m_DCT4(x, mode)
% function [x_dct, T] = m_DCT4(x, mode)
% Author: Archit Gupta (September 06, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of DCT-IV. DCT transform of the fourth kind.
%   Refer WIKIPEDIA or <CS_DIR>/4_Notes/Discrete-cosine-transform.xoj for more
%   details. For a given sequence of values x[n], this computes the cosine
%   transform of an even extension about a point half a unit distance away from
%   the left end and odd extension about a point half a unit distance from the
%   right end.
%   
%   INPUTS:
%       x - A sequence of values for which DCT-4 has to be computed [Expect a
%       column vector. A row vector will produce an ERROR]
%       mode (optional) - This can be left empty if the regular DCT is to be
%       computed. However, if we want to compute the DCT matrix by using the DFT
%   OUTPUTS:
%       x_dct - A sequence of DCT coefficients (same length as x)
%       T (optional) - The tranform matrix associated with the size of input x
%
%   NOTE: This code is experimental. It does not attain the complexity of
%   O{nlog(n)}, where n is the length of the vector x. This has a complecity of
%   O(n^2), which is probably good enough for the small scale experiments that
%   we are doing. For better efficiency (atleast better asymptotic efficiency)
%   use DFT instead with an appropriately extended sequence.
%
%   TRANSFORM:
%       X(k) = sum_{n=0}^{N-1} x(n)*cos(pi*(n+0.5)*(k+0.5)/N)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin > 1)
    mode_extended = (mode == 'MODE_EXTENDED');
else
    mode_extended = 0;
end

[N, ~] = size(x);

if (N == 1)
    fprintf(2, ['WARNING: If row vector was given the output is wrong.\n', ...
        'Otherwise DCT4 is trivial, what is wrong with you!\n']);
    x_dct = x;
    T = 1;
    return;
end

% Let k denote the index of the DCT coefficient and let n denote the index of
% the time domain coefficient
k = [0:N-1]';
n = [0:N-1]';

% Baseline code: Corresponds to a DFT of a sequence of length 8N. TODO: Figure
% out the exact details. Anyways, the transform can be made orthogonal by
% multiplication by a constant factor sqrt(2/N)
T = cos(pi*(k+0.5)*(n+0.5)'/N);

if (~mode_extended)
    % Scaled implementation
    T = sqrt(2/N)*T;
end

x_dct = T*x;
