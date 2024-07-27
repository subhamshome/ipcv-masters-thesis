function [wiener_sol, mse] = wiener_filter(Data, Mask, mu, TrueObject)   
    % WIENER_FILTER Function to perform Wiener filtering for image reconstruction.
    % Data: Fourier domain measurements.
    % Mask: Fourier mask.
    % mu: Regularization parameter.
    % TrueObject: Ground truth image for MSE calculation.

    % Define Difference Matrices for finite difference calculations
    d_col = [[0  0 0]
             [0 -1 1]
             [0  0 0]];  % Column-wise difference matrix
    d_row = d_col';  % Row-wise difference matrix (transpose of d_col)

    % Compute the FFT of the difference matrices
    d_col_fourier = MyFFT2RI(d_col, size(Mask, 1));  % FFT of column difference matrix
    d_row_fourier = MyFFT2RI(d_row, size(Mask, 2));  % FFT of row difference matrix

    % Compute conjugate of Mask in Fourier domain
    xNum = conj(Mask);

    % Compute denominator terms for Wiener filter
    HSquareFFT = abs(Mask) .^ 2;  % Square of Mask magnitude
    DSquareFFT = abs(d_col_fourier) .^ 2 + abs(d_row_fourier) .^ 2;  % Sum of squared difference matrices
    xDenom = HSquareFFT + mu * DSquareFFT;  % Combined denominator term

    % Calculate Wiener Gain
    wiener_gain = xNum ./ xDenom;

    % Apply Wiener Gain to the Data
    wiener_sol_fourier = wiener_gain .* Data;

    % Transform solution back to spatial domain
    wiener_sol_spatial = MyIFFT2(wiener_sol_fourier) ./ size(Mask, 1);

    % Ensure solution is real-valued
    wiener_sol = real(wiener_sol_spatial);

    % Output max imaginary part of Wiener solution (should be zero)
    fprintf("Wiener spatial:" + max(abs(imag(wiener_sol_spatial(:)))) + "\n");

    % Calculate Mean Squared Error (MSE) between reconstructed and true image
    mse = calculateMSE(TrueObject, wiener_sol);

    % Output MSE
    fprintf("Wiener MSE: " + mse + "\n");
end
