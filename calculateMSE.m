function normalizedSquaredError = calculateMSE(reconstructedImage, trueImage)
    % Compute the squared differences between reconstructed and true images
    squaredDifferences = (reconstructedImage - trueImage).^2;

    % Compute the normalized squared error
    normalizedSquaredError = sum(squaredDifferences, 'all') / sum(trueImage.^2, 'all');
end