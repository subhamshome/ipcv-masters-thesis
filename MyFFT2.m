function Frequentiel = MyFFT2(Spatial)	

% Calcul de la transformée et normalisation
	Frequentiel = fftshift( fft2(Spatial) ) / length(Spatial) ;

