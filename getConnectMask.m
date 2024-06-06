function [mask] = getConnectMask(EEG, par)

n = EEG.nbchan;
mask = logical(eye(n));
radius = median([EEG.chanlocs.sph_radius]);
X = [EEG.chanlocs.X]./radius;
Y = [EEG.chanlocs.Y]./radius;
Z = [EEG.chanlocs.Z]./radius;
tri=delaunayn([X',Y',Z']);

if exist('par') && ~isnumeric(par) && strncmpi(par,'del',3)
    % Lambert projection: make the Z zero to be the bottom of a full sphere
    % (in case of a spherical headmodel.)
    rho = sqrt(2 * (1 - Z + max(Z)));
    X_proj = rho .* (X ./ sqrt(X.^2 + Y.^2));
    Y_proj = rho .* (Y ./ sqrt(X.^2 + Y.^2));

    % then Delauney 
    tri = delaunay(X_proj,Y_proj);
    for row=1:size(tri,1)
        mask(tri(row,1), tri(row,2)) = 1;
        mask(tri(row,2), tri(row,3)) = 1;
        mask(tri(row,1), tri(row,3)) = 1;
    end
else
    % do a simple distance threshold
    radius = median([EEG.chanlocs.sph_radius]);
    if ~exist('par')
        par = 4.5*radius/sqrt(n);
    end

    dX = bsxfun(@minus, X, X.');
    dY = bsxfun(@minus, Y, Y.');
    dZ = bsxfun(@minus, Z, Z.');
    distances = sqrt(dX.^2+dY.^2+dZ.^2);
    
    comb = nchoosek(1:size(tri,2),2);
    
    for p=1:size(tri,1)
        for i=1:size(comb,1)
            r = tri(p,comb(i,1));
            c = tri(p,comb(i,2));
            if distances(r,c)<par
                mask(r,c)=1;
            end
        end
    end
   
end

mask = mask|mask';