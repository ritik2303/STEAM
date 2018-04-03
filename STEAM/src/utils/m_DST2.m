function [x_dst, T] = m_DST2(x)
% function [x_dst, T] = m_DST2(x)
% Author: Archit Gupta (September 13, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Implementation of DST-II. DST Transform of the second kind
%   Refer WIKIPEDIA or <CS_DIR>/4_Notes/Discrete-sine-transform.xoj for more
%   details. For a given sequence of values x[n], this computes the sine
%   transform of an odd extension about a point  which is half a unit distance
%   away from both the left and right ends.
%
%   INPUTS:
%       x - A sequence of values for which DST-2 has to be computed [Expect a
%       column vector or a matrix for which the DST will be computed columnwise.
%       A row vector will produce an error.]
%
%   OUTPUTS:
%       x_dst - A sequence of DST coefficients (same length/size as x)
%       T (optional) - The orthogonal transform matrix associated with the
%       length of the input
%
%   TRANSFORM
%       X(k) = sum_{n=0}^{N-1} x(n)*sin(pi*(n+0.5)*(k+1)/N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, ~] = size(x);

if (N == 1)
    fprintf(2, ['WARNING: If row vector was given, the output is wrong. \n', ...
        'Othewise DST1 is trivial, what''s wrong with you!\n']);
    x_dst = x;
    T = 1;
end

% Let k denote the index of the DST coefficient and let n denote the index of
% the time domain coefficient

k = [0:N-1]';
n = [0:N-1]';

% Baseline code, corresponding directly to a DFT of an odd extension of the
% series of points about the left and right end. The transform matrix is
% orthogonal. It has been made orthonormal in the code in use.
%T = sin(pi*(k+1)*(n+1)'/(N+1));

T = sqrt(2/N)*sin(pi*(k+1)*(n+0.5)'/N);
T(end,:) = T(end,:)/sqrt(2);
x_dst = T*x;
