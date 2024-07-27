function Spatial = MyIFFT2(Frequentiel)	

% Calcul de la transformée inverse et normalisation
	Spatial = ifft2( fftshift(Frequentiel) ) * length(Frequentiel);


