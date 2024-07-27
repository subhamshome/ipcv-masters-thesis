function Frequentiel = MyFFT2RI(Spatial,Long)	

% Complément - Bourrage
    Taille = length(Spatial);
    Ou = 1+Long/2-(Taille-1)/2 : 1+Long/2+(Taille-1)/2;
    SpatialComplet = zeros(Long,Long);
    SpatialComplet( Ou , Ou ) = Spatial;

% Calcul de la transformée
	Frequentiel = fftshift( fft2( fftshift(SpatialComplet) ) );
