function [admm_sol, mse] = admm(Data, Mask, mu, num_iter, rho)
    % ADMM function for image reconstruction.
    % Data: Fourier domain measurements.
    % Mask: Fourier mask.
    % mu: Regularization parameter.
    % num_iter: Number of iterations.
    % rho: ADMM penalty parameter.

    % Define Difference Matrices for finite difference calculations
    d_col = [[0  0 0]
             [0 -1 1]
             [0  0 0]];  % Column-wise difference matrix
    d_row = d_col';  % Row-wise difference matrix (transpose of d_col)

    current_mse = zeros(num_iter, 1);  % Initialize MSE tracking array

    % Compute the FFT of the difference matrices
    d_col_fourier = MyFFT2RI(d_col, size(Mask, 1));  % FFT of column difference matrix
    d_row_fourier = MyFFT2RI(d_row, size(Mask, 2));  % FFT of row difference matrix

    % Compute conjugate of Mask in Fourier domain
    HConjFFT = conj(Mask);

    % Compute denominator terms for ADMM update equation
    HSquareFFT = abs(Mask) .^ 2;  % Square of Mask magnitude
    DSquareFFT = abs(d_col_fourier) .^ 2 + abs(d_row_fourier) .^ 2;  % Sum of squared difference matrices
    xDenom = HSquareFFT + mu * DSquareFFT + rho;  % Combined denominator term

    % Initialize lambda (Lagrange multipliers) and rec_image_max (max filtered image)
    lambda = zeros(size(Data));  % Initialize Lagrange multipliers
    lambda_fourier = MyFFT2(lambda);  % FFT of lambda
    rec_image_max = zeros(size(Data));  % Initialize max filtered image
    rec_image_max_fft = MyFFT2(rec_image_max);  % FFT of max filtered image

    % Initialize previous solution for MSE calculation
    previous_admm_sol = zeros(size(Data));

    % ADMM Iterative Algorithm
    for iteration = 1:num_iter 
        % Numerator term for ADMM update
        xNum = HConjFFT .* Data - lambda_fourier/2 + rho * rec_image_max_fft;
        
        % Compute updated solution in Fourier domain
        admm_sol_fourier = xNum ./ xDenom;
        
        % Transform solution back to spatial domain
        admm_sol_spatial = MyIFFT2(admm_sol_fourier);
        
        % Ensure solution is real-valued
        fprintf("Iteration " + iteration + " Positivity reconstruct: " + max(abs(imag(admm_sol_spatial(:)))) + "\n");
        admm_sol = real(admm_sol_spatial);
        
        % Update max filtered image
        rec_image_max = max(0, admm_sol + (lambda / rho));
        
        % Update Lagrange multipliers
        lambda = lambda + rho * (admm_sol - rec_image_max);
        lambda_fourier = MyFFT2(lambda);

        % Calculate Mean Squared Error with respect to previous iteration
        current_mse(iteration) = calculateMSE(previous_admm_sol, admm_sol);

        % Check for convergence
        if iteration > 1 && current_mse(iteration) < 1e-14
            fprintf('Converged at iteration %d\n', iteration);
            break;
        end

        % Update previous solution
        previous_admm_sol = admm_sol;
    end

    % Compute average MSE over all iterations
    mse = mean(current_mse);
end
