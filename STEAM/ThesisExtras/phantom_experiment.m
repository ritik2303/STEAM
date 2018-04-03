% Assuming you have a recent version of MATLAB.
% This script generates the phantom, shows the image and its fourier transform.
% It then takes a sparse measurement of the Fourier samples along radial lines
% aand uses just these values to perform Total Variance minimization on the
% Fourier coefficients.
%
% Adopted from text in: @CITE-PAPER

% Create the Phantom
p_img = phantom;
figure, imshow(p_img)

% Get the Fourier coefficients (Take the logarithm to actually be able to see
% something), and plot them after normalizing the values.

fft_img = log(abs(fft2(p_img)));
fft_img = fft_img / max(max(fft_img));
figure, imshow(fftshift(fft_img));

% Generate RADIAL lines, along which the coefficients will be measured.
theta = [pi/12, pi/6, pi/4, pi/3, 5*pi/12, 7*pi/12, 2*pi/3, 3*pi/4, 5*pi/6, 11*pi/12];
mask  = zeros(size(fft_img));
orgn  = round([size(fft_img, 1)/2 size(fft_img, 2)/2]);

mask(orgn(1), :) = 1; % Horizontal Axis
mask(:, orgn(2)) = 1; % Vertical Axis

% Set the limits for x-value
x_pos_bound = size(fft_img, 1);
x_neg_bound = 1;

y_pos_bound = size(fft_img, 2);
y_neg_bound = 1;

for lin = 1:length(theta)
    p_val = orgn;
    pt    = orgn;  
    x_inc = cos(theta(lin));
    y_inc = sin(theta(lin));

    % Unmask values with positive x
    while (pt(1) < x_pos_bound) && (pt(1) > x_neg_bound)
        % Check that y is still in bounds
        if (pt(2) > y_pos_bound) || (pt(2) < y_neg_bound)
            break;
        end
        mask(pt(1), pt(2)) = 1;
        p_val(1) = p_val(1) + x_inc;
        p_val(2) = p_val(2) + y_inc; 
        pt = round(p_val);
    end

    p_val = orgn;
    pt    = orgn;
    % Unmask values with negative x
    while (pt(1) > x_neg_bound) && (pt(1) < x_pos_bound)
        % Check that y is still in bounds
        if (pt(2) > y_pos_bound) || (pt(2) < y_neg_bound)
            break;
        end

        mask(pt(1), pt(2)) = 1;
        p_val(1) = p_val(1) - x_inc;
        p_val(2) = p_val(2) - y_inc; 
        pt = round(p_val);
    end
end

figure, imshow(mask);

% Show the sampled Fourier Coefficients.
figure, surf(fftshift(fft_img) .* mask);
figure, imshow(fftshift(fft_img) .* mask);