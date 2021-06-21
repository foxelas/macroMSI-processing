function [melMap, hbMap] = getMap2(melMap, hbMap, mapType, mask, id)
%     GETMAP returns an optical map showing chromophore components of the macropathology image
%
%     Usage:
%     [melMap, hbMap] = getMap(msi)
%     [melMap, hbMap] = getMap(msi, 'ding')
%
%     Available mapType values
%     ding: uses optical density
%     vasefi: uses absorption slope
%     ours: uses estimated reflectance and absorption slope

savedir = getSetting('savedir');
mapdir = getSetting('map');

if isstruct(id)
    id = id.Index;
end

switch lower(mapType)
    case 'ding'
        melBarTitle = 'Optical Density of Melanin (a.u.)';
        melSaveName = 'DingMel';
        hbBarTitle = 'Optical Density of Hemoglobin (a.u.)';
        hbMapSaveName = 'DingHb';

    case 'vasefi'
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        melBarTitle = 'Relative Melanin Concentration (a.u.)';
        melSaveName = 'VasefiMel';
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        hbBarTitle = 'Relative Total Hemoglobin Concentration (a.u.)';
        hbMapSaveName = 'VasefiHb';       

    case 'diebele'
        melBarTitle = 'Melanin Index (a.u.)';
        melSaveName = 'DiebeleMel';
        hbBarTitle = 'Erythema Index (a.u.)';
        hbMapSaveName = 'DiebeleHb';

    case 'kapsokalyvas'
        melBarTitle = 'Melanin Homogeneity (a.u.)';
        melSaveName = 'KapsokalyvasMel';
        setSetting('plotName', fullfile(savedir, mapdir, strcat('KapsokalyvasSupMel', '_', num2str(id), '_', 'ScaledMap.png')));
        hbBarTitle = 'Hemoglobin Homogeneity (a.u.)';
        hbMapSaveName = 'KapsokalyvasHb';
    otherwise
        disp('Unsupported type')
end


setSetting('hasLimit', true);
melMapLimits = [0,1];
hbMapLimits = [0,1];
setSetting('plotName', fullfile(savedir, mapdir, strcat(melSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plots(1, @plotMap, melMap, mask, [], false, melBarTitle, melMapLimits);
setSetting('plotName', fullfile(savedir, mapdir, strcat(hbMapSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plots(2, @plotMap, hbMap, mask, [], false, hbBarTitle, hbMapLimits);


end
