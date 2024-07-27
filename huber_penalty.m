function [huber_sol, mse] = huber_penalty(Data, Mask, T, alpha, mu, num_iter, showPlots)
    % HUBER_PENALTY Function to perform Huber penalty regularization for image reconstruction.
    % Data: Fourier domain measurements.
    % Mask: Fourier mask.
    % T: Threshold for Huber penalty.
    % alpha: Regularization parameter for Huber penalty.
    % mu: Regularization parameter.
    % num_iter: Number of iterations.
    % showPlots: Flag to display plots.
    
    % Define constants
    mu_p = mu / alpha;  % Adjusted regularization parameter
    mse = zeros(num_iter, 1);  % Preallocate MSE array
    
    % Define Difference Matrices for finite difference calculations
    d_col = [[0  0 0]
             [0 -1 1]
             [0  0 0]];  % Column-wise difference matrix
    d_row = d_col';  % Row-wise difference matrix (transpose of d_col)
    
    % Compute the FFT of the difference matrices
    d_col_fourier = MyFFT2RI(d_col, size(Mask, 1));  % FFT of column difference matrix
    d_row_fourier = MyFFT2RI(d_row, size(Mask, 2));  % FFT of row difference matrix
    
    % Compute denominator terms for Huber penalty
    HSquareFFT = abs(Mask) .^ 2;  % Square of Mask magnitude
    DSquareFFT = abs(d_col_fourier) .^ 2 + abs(d_row_fourier) .^ 2;  % Sum of squared difference matrices
    xDenom = HSquareFFT + mu_p * DSquareFFT;  % Combined denominator term

    % Compute initial interpixel differences
    delta_col = MyIFFT2(d_col_fourier .* Data);
    delta_row = MyIFFT2(d_row_fourier .* Data);
    
    % Update auxiliary variables based on Huber penalty
    aux_col = (1 - 2 * alpha * min(1, T ./ abs(delta_col))) .* delta_col;
    aux_row = (1 - 2 * alpha * min(1, T ./ abs(delta_row))) .* delta_row;
    
    % Compute conjugate of Mask and difference matrices in Fourier domain
    HConjFFT = conj(Mask);
    DColConjFFT = conj(d_col_fourier);
    DRowConjFFT = conj(d_row_fourier);
    
    % Compute FFT of auxiliary variables
    aux_col_fourier = MyFFT2(aux_col);
    aux_row_fourier = MyFFT2(aux_row);
    
    % Initialize solution
    huber_sol = real(MyIFFT2(Data));

    for iteration = 1:num_iter        
        % Compute numerator for Huber penalty minimization
        xNum = HConjFFT .* Data + mu_p * DColConjFFT .* aux_col_fourier + mu_p * DRowConjFFT .* aux_row_fourier;
        
        % Solve for Huber solution in Fourier domain
        xFFT = xNum ./ xDenom;
        
        % Update interpixel differences
        delta_col = MyIFFT2(d_col_fourier .* xFFT);
        delta_row = MyIFFT2(d_row_fourier .* xFFT);
        
        % Update auxiliary variables
        aux_col = (1 - 2 * alpha * min(1, T ./ abs(delta_col))) .* delta_col;
        aux_row = (1 - 2 * alpha * min(1, T ./ abs(delta_row))) .* delta_row;

        aux_col_fourier = MyFFT2(aux_col);
        aux_row_fourier = MyFFT2(aux_row);
        
        % Update solution
        previousData = huber_sol;
        huber_sol_spatial = MyIFFT2(xFFT);
        huber_sol = real(huber_sol_spatial);
        fprintf("Iteration " + iteration + " Huber spatial:" + max(abs(imag(huber_sol_spatial(:)))) + "\n");
        
        % Calculate Mean Squared Error (MSE)
        mse(iteration) = calculateMSE(previousData, huber_sol); 

        % Check for convergence
        if calculateMSE(huber_sol, previousData) < 1e-14
            break;
        end
    end
    
    % Display plots if required
    if showPlots        
        figure("Name", sprintf("alpha = %.3f", alpha));
        set(gcf, 'WindowState', 'maximized');
    
        subplot(132);
        imagesc(huber_sol);
        title('Reconstructed Image');
        colormap('gray'); colorbar; axis('square');
   
        subplot(133);
        plot(mse(mse > 0));
        title('Mean Squared Error');
        axis square;
        xlabel('Number of iterations');
        ylabel('MSE');
    
        sgtitle('Huber Penalty Convergence');
    end
    
    % Compute final mean squared error
    mse = mean(mse);
    
    fprintf("Huber spatial:" + max(abs(imag(huber_sol_spatial(:)))) + "\n");
    fprintf("Huber MSE: " + mse + "\n");
end
