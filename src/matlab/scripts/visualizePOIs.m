setSetting('saveInHQ', false);

for g = 1:max([ID.Group])

    %% Read MSI
    infile = fullfile(getSetting('systemdir'), 'infiles', strcat('group_', num2str(g), '.mat'));
    load(infile, 'whiteReference', 'specimenMask');
    gIdxs = find([ID.Group] == g);
    gMembers = ID(gIdxs);
    coordinates = [[gMembers.Originx]; [gMembers.Originy]]';
    baseImage = whiteReference .* specimenMask;
    [m, n, ~] = size(baseImage);
    labels = {gMembers.Label};
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('labels'), strcat('labelled_bright_', num2str(g))));
    plotVisualResult(baseImage, zeros(m, n), '', labels, coordinates, 'jet', true, 2);
    %     for i = 1:length(coordinates)
    %         plotVisualResult(whiteReference, zeros(size(whiteReference,1),...
    %             size(whiteReference,2)), num2str(idd.Index + i), labels(i), coordinates(i,:), 'jet', true, 3);
    %     end
end
setSetting('saveInHQ');