function [x_dct, T] = m_DCT3(x, mode)
% function [x_dct, T] = m_DCT3(x, mode)
% Author: Archit Gupta (September 06, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of DCT-III. DCT transform of the third kind.
%   Refer WIKIPEDIA or <CS_DIR>/4_Notes/Discrete-cosine-transform.xoj for more
%   details. For a given sequence of values x[n], this computes the cosine
%   transform of an even extension about the left end and odd extension about
%   the right end.
%   
%   INPUTS:
%       x - A sequence of values for which DCT-3 has to be computed [Expect a
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
%       X(k) = x0/2 + sum_{n=1}^{N-1} x(n)*cos(pi*n*(k+0.5)/N)
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
        'Otherwise DCT3 is trivial, what is wrong with you!\n']);
    x_dct = x;
    T = 1;
    return;
end

% Let k denote the index of the DCT coefficient and let n denote the index of
% the time domain coefficient
k = [0:N-1]';
n = [0:N-1]';

% Baseline code: MATLAB implementation (and several others) require scaling of
% the overall matrix and also the X0 term (so that the transform matrix is
% orthogonal). See scaled implemetation here to details

T = cos(pi*(k+0.5)*n'/N);
T(:,1) = T(:,1)/2;

if (~mode_extended)
    % Scaled implementation
    T = sqrt(2/N)*T;
    T(:,1) = sqrt(2)*T(:,1);
end

x_dct = T*x;
