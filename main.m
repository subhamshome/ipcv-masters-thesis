% Clean up workspace
clearvars; close all; clear; clc;
disp('Start');

% Define parameters
Size = 128;  % Size of the unknown image
theta_step = 2;  % Step size for theta
TheTheta = 1:theta_step:180;  % Theta values
NbreTheta = length(TheTheta);  % Number of theta values
Long = 256;  % Length of frequency grid

% Generate unknown image using Shepp-Logan phantom
Offset = 0;  % Offset value for padding
Coef = 1;  % Coefficient for phantom generation
TrueObject = phantom('Modified Shepp-Logan', Size);
PadSize = 0.5 * round(1 * Size / 2);  % Padding size
TrueObject = padarray(TrueObject, [PadSize PadSize], Offset, 'both');
fprintf("Max Imag Part - TrueObject: %f\n", max(abs(imag(TrueObject(:)))));

% Compute Fourier Transform of the true object
Size = size(TrueObject, 1);
TrueObjectTF = fftshift(fft2(TrueObject)); 
fprintf("Max Imag Part - TrueObjectTF: %f\n", max(abs(imag(TrueObjectTF(:)))));

% Display the true object and its Fourier Transform
figure(10); clf;
set(gcf, 'WindowState', 'maximized');
subplot(121);
    imagesc(TrueObject);
    axis square; colormap('gray'); colorbar; title('True');
subplot(122);
    TmpFreqGrid = linspace(-0.5, 0.5, size(TrueObjectTF, 1));
    imagesc(TmpFreqGrid, TmpFreqGrid, log(abs(TrueObjectTF))); 
    axis square; colormap('gray'); colorbar; title('LogAbs 2D-FFT of True');
sgtitle('True/Unknown Image');

% Perform Fourier mask (Radon transform like)
FreqGrid = 0.99 * linspace(-1/2, 1/2 - 1/Long, Long); 
TheKx = FreqGrid' * sind(TheTheta);
TheKy = FreqGrid' * cosd(TheTheta);

% Create mask using Radon transform
[Mask, Nbre] = ProcReGridPPV(TheKx, TheKy, ones(Long, NbreTheta), round(1 * Size));
Mask(find(Mask)) = 1;

% Display the mask and Nbre
figure(11); clf;
subplot(121);
    imagesc(Mask);
    axis square; colormap('gray'); colorbar; title('Mask');
subplot(122);
    imagesc(log(1 + Nbre));
    axis square; colormap('hot'); colorbar; title('Nbre');
impixelinfo;

% Compute Data (in the Fourier domain)
RePadSize = (size(Mask, 1) - size(TrueObjectTF, 1)) / 2;
PadTrueObjectTF = padarray(TrueObjectTF, [RePadSize RePadSize], 0, 'both');
fprintf("Max Imag Part - PadTrueObjectTF: %f\n", max(abs(imag(PadTrueObjectTF(:)))));
Data = PadTrueObjectTF .* Mask + 2 * MyFFT2(randn(size(Mask))); % Add noise to data
fprintf("Max Imag Part - Data after Mask Application: %f\n", max(abs(imag(Data(:)))));

% Reconstruct images
ReTrue = ifft2(fftshift(TrueObjectTF));
ReTruePad = real(ifft2(fftshift(PadTrueObjectTF))); 
ReMask = real(fftshift(ifft2(Mask)));
ReBuild = real(ifft2(fftshift(Data))); 

% Display reconstructed images
figure(12);
set(gcf, 'WindowState', 'maximized');
subplot(241);
    imagesc(ReTrue);
    axis square; colormap('gray'); colorbar; title('True');
subplot(242);
    imagesc(ReTruePad);
    axis square; colormap('gray'); colorbar; title('True Padded');
subplot(243);
    imagesc(log(abs(ReMask)));
    axis square; colormap('gray'); colorbar; title('Mask Re-Build');
subplot(244);
    imagesc(ReBuild);
    axis square; colormap('gray'); colorbar; title('Re-Build');
subplot(245);
    TmpFreqGrid = linspace(-0.5, 0.5, size(TrueObjectTF, 1));
    imagesc(TmpFreqGrid, TmpFreqGrid, log(abs(TrueObjectTF))); 
    axis square; colormap('gray'); colorbar; title('TF True');
subplot(246);
    TmpFreqGrid = linspace(-0.5, 0.5, size(PadTrueObjectTF, 1));
    imagesc(TmpFreqGrid, TmpFreqGrid, log(abs(PadTrueObjectTF))); 
    axis square; colormap('gray'); colorbar; title('IFFT Pad TF True');
subplot(247);
    TmpFreqGrid = linspace(-0.5, 0.5, size(TrueObjectTF, 1));
    imagesc(TmpFreqGrid, TmpFreqGrid, Mask); 
    axis square; colormap('gray'); colorbar; title('Mask');
subplot(248);
    TmpFreqGrid = linspace(-0.5, 0.5, size(PadTrueObjectTF, 1));
    imagesc(TmpFreqGrid, TmpFreqGrid, log(abs(Data))); 
    axis square; colormap('gray'); colorbar; title('Data');
sgtitle('Re-gridding Algorithm');
saveas(gcf, 'regridding.png');

% Wiener filter deconvolution
mu_wiener = 1e0; % Regularization parameter
[wiener_sol, wiener_mse] = wiener_filter(Data, Mask, mu_wiener, TrueObject);

% Display Wiener filter result
figure(13);
set(gcf, 'WindowState', 'maximized');
subplot(131);
    imagesc(ReTrue);
    axis square; colormap('gray'); colorbar; title('True');
subplot(132);
    imagesc(ReBuild);
    axis square; colormap('gray'); colorbar; title('Re-Build');
subplot(133);
    imagesc(wiener_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using Wiener Filter');
sgtitle('Wiener Filtering');

% Huber Penalty 
T = 0.1; % Threshold
alpha = 0.5; 
num_iter = 2000; 
mu_huber = 1e0;
[huber_sol, huber_mse] = huber_penalty(Data, Mask, T, alpha, mu_huber, num_iter, false);

% Display Huber penalty result
figure(15);
set(gcf, 'WindowState', 'maximized');
subplot(221);
    imagesc(ReTrue);
    axis square; colormap('gray'); colorbar; title('True');
subplot(222);
    imagesc(ReBuild);
    axis square; colormap('gray'); colorbar; title('Re-gridding output');
subplot(223);
    imagesc(wiener_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using Wiener Filter');
subplot(224);
    imagesc(huber_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using Huber Penalty');
sgtitle('Final Algorithm');
saveas(gcf, 'final_algorithm.png');

% ADMM Solution
mu_admm = 0.04;
rho = 0.7;
num_iter_admm = 2000;
[admm_sol, admm_mse] = admm(Data, Mask, mu_admm, num_iter_admm, rho);

% Display ADMM result
figure(19);
set(gcf, 'WindowState', 'maximized');
subplot(221);
    imagesc(ReTrue);
    axis square; colormap('gray'); colorbar; title('True');
subplot(222);
    imagesc(wiener_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using Wiener Filter');
subplot(223);
    imagesc(huber_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using Huber Penalty');
sgtitle('Huber Penalty');
subplot(224);
    imagesc(admm_sol);
    axis square; colormap('gray'); colorbar; title('Final Reconstructed Image using ADMM');
sgtitle('Wiener vs Huber vs ADMM');

%%
% Central Line Comparison
% Get the central horizontal line from each image
center_line_index = floor(size(ReTrue, 1) / 2) + 1;
line_ReTrue = ReTrue(center_line_index, :);
line_ReBuild = ReBuild(center_line_index, :);
line_wiener_sol = wiener_sol(center_line_index, :);
line_huber_sol = huber_sol(center_line_index, :);
line_admm_sol = admm_sol(center_line_index, :);

% Normalize Huber and ADMM center line to correct range
line_huber_sol = (line_huber_sol - min(line_huber_sol)) / (max(line_huber_sol) - min(line_huber_sol));
line_admm_sol = (line_admm_sol - min(line_admm_sol)) / (max(line_admm_sol) - min(line_admm_sol));

% Create a new figure for the central horizontal lines
figure(25);
set(gcf, 'WindowState', 'maximized');
hold on;
    plot(line_ReTrue, 'DisplayName', 'True', 'LineWidth', 2, 'Color', 'yellow'); 
    plot(line_ReBuild, 'DisplayName', 'ReBuild', 'LineWidth', 2, 'Color', 'red'); 
    plot(line_wiener_sol, 'DisplayName', 'Wiener Filter', 'LineWidth', 2, 'Color', 'green', 'LineStyle', '--'); 
    plot(line_huber_sol, 'DisplayName', 'Huber Penalty', 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '--'); 
    plot(line_admm_sol, 'DisplayName', 'ADMM', 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--'); 
hold off;

% Add labels, legend, and title
xlabel('Pixel Index');
ylabel('Intensity');
legend;
title('Central Horizontal Line Comparison');
grid on;
axis square;
ylim([-0.1, 1.1]);

