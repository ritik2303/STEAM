function [x_dct, T] = m_DCT1(x, mode)
% function [x_dct, T] = m_DCT1(x, mode)
% Author: Archit Gupta (September 05, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implementation of DCT-I. DCT transform of the first kind.
%   Refer WIKIPEDIA or <CS_DIR>/4_Notes/Discrete-cosine-transform.xoj for more
%   details. For a given sequence of values x[n], this computes the cosine
%   transform of an even extension about both the left and right ends.
%   
%   INPUTS:
%       x - A sequence of values for which DCT-1 has to be computed [Expect a
%       column vector. A row vector will produce an ERROR. Any vector with
%       length < 2 will also produce an error]
%       
%       NOTE: For a matrix input, a 2D DCT is used. TODO: Doing a columnwise DCT
%       is easy but not supported at the moment. This can be done easily by
%       adding a mode
%
%       mode (optional) - This can be left empty if the regular DCT is to be
%       computed. However, if we want to compute the DCT matrix by using the DFT
%       formula (without normalization), pass 'MODE_EXTENDED'  in mode
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
%       X(k) = (x_{0} + ((-1)^k)*x_{N-1})/2 + ...
%           sum_{n=1}^{N-2} x(n)*cos(pi*n*k/(N-1))
%
%   MODIFIED: Feb 15, 2017
%   Extension to arbitrary number of dimensions (TODO:  Current implementation
%   takes O(n^2) time. We should get this down to O(n*log(n)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    s_x = size(x);
    if ( length(s_x) > 2 )
        fprintf(2, ['ERROR: DCT for inputs having more than 2 dimensions', ...
            ' not supported at the moment\n']);
        x_dct = 0;
        T = 0;
        return;
    end

    if (s_x(1) == 1)
        fprintf(2, ['ERROR: Either row vector was given or signal length < 2\n', ...
            'DCT1 not defined!\n']);
        x_dct = 0;
        T = 0;
        return;
    end

    % Baseline code: Corresponds directly to a DFT of the even extension of the
    % series of points about the left and right ends. However, the transform can be
    % made orthogonal by dividing the transform matrix T by sqrt(2/N-1)

    T = getDCTMat(s_x(1), mode);

    if ( s_x(2) > 2 )
        T_hat = getDCTMat(s_x(2), mode);
        x_dct =  T * x * T_hat';
    else
        x_dct = T * x;
    end
end

function T = getDCTMat(n_elems, extension)
    % Let k denote the index of the DCT coefficient and let n denote the index of
    % the time domain coefficient
    k = [0:n_elems-1]';
    n = [0:n_elems-1]';

    T = cos(pi*k*n'/(n_elems-1));
    T(:,1) = T(:,1)/2;
    T(:,end) = T(:,end)/2;

    if ( strcmp(extension, 'MODE_EXTENDED') )
        return;
    end

    % n_elemsormalization
    if ( strcmp(extension, 'CHEB_SERIES') )
        T(1,:) = T(1,:)/2;
        T(end,:) = T(end,:)/2;
        T = ( 2/(n_elems-1) ) * T;
        return;
    end

    T = sqrt(2/(n_elems-1))*T;
    return;
end
