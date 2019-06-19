groups = findgroups([ID.MsiID]);
for g = 1:max(groups)
%% Read MSI
    gIdxs = find(groups == g);
    gMembers = ID(gIdxs);
    coordinates = [[gMembers.Originx]; [gMembers.Originy]]';
    idd = gMembers(1);
    files = {data([data.MsiID] == idd.MsiID).File};
    [~, whiteReference, ~] = readMSI(files); 
    labels = ~[gMembers.IsBenign];
    options.saveOptions.plotName = fullfile('..','..','..', 'output', 'general', strcat('labelled_', num2str(g), '.png'));
    options.saveOptions.saveImages = true;
    options.saveOptions.saveInHQ = false;
    plots('visual', 3, labels, '', 'Image', whiteReference, 'Overlay', zeros(size(whiteReference,1), size(whiteReference,2)), 'Cmap', 'jet', 'Alpha', 0.7, ...
    'Coordinates', coordinates,'saveOptions', options.saveOptions);

end