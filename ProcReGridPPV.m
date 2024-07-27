% ProcReGridPPV.m 		Gio le 15-01-02
%
% Regridding sur grille cart�sienne au plus proche voisin
% (forme exacte, avec moyennage des donn�es "multiples"
% (forme rapide, sans boucle : round + find + diff (cumsum), ...
%
% Tailee -> size

function [DataReGri , MatNbre] = ProcReGridPPV(kx,ky,Zd,Taille)

% Re-gridde sur cartesien
	kx = round((kx+0.5)*Taille); 
	ky = round((ky+0.5)*Taille);
	kx = kx+1; 
	ky = ky+1;

%%%% Construction des nouvelles donn�es
	% Initialisation de matrices
		DataReGri = zeros(Taille,Taille);
		MatNbre = zeros(Taille,Taille);

	% Traitement des donn�es des (kx,ky) n'apparaissant qu'une fois
		DataReGri(kx+Taille*(ky-1)) = Zd;
		% MatNbre(kx+Taille*(ky-1)) = 1;

	% Colle tout en colonne	
		kx = kx(:); ky = ky(:); Zd = Zd(:);

	% Fabrique un indice complexe
		K = kx+sqrt(-1)*ky;
		clear kx ky 

	% Tri en fonction du module des indices complexes (et les donn�es sont tri�es avec)
 		Sort = sortrows([K Zd],1); 
		KSort = Sort(:,1);
		ZdSort = Sort(:,2);
		clear Sort %K Zd

	% D�tecte les d�buts (par 1) et les fins( par -1) de "s�quences" 
		DebutFin = zeros(size(KSort));
		DebutFin(diff(KSort)==0) = 1;

		DebutFin = diff([0;DebutFin]);
		DebutFin = [DebutFin ; Inf];

	% Extrait les K et les indices des d�buts et fins de s�quences
		KSortDebut = KSort(DebutFin==1);
		IndDebut = (find(DebutFin==1));
		IndFin =  (find(DebutFin==-1));
		clear DebutFin KSort

	% Calcul de la moyenne des donn�es corespondantes
		Moy = cumsum([0;ZdSort]);
		% Nbre = cumsum(ones(size(ZdSort)));

		Moy = Moy(IndFin+1) - Moy(IndDebut);
		% Nbre = Nbre(IndFin+1) - Nbre(IndDebut);
		% Moy = Moy ./ Nbre;
		clear IndFin IndDebut ZdSort

	% Mise en place 
		DataReGri( real(KSortDebut)+Taille*(imag(KSortDebut)-1) ) = Moy;
		DataReGri = reshape(DataReGri,Taille,Taille);
		% MatNbre( real(KSortDebut)+Taille*(imag(KSortDebut)-1) ) = Nbre;
		% MatNbre = reshape(MatNbre,Taille,Taille);
		clear Moy KSortDebut Nbre
