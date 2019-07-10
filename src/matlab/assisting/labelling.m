load( fullfile('..','..','..', 'input', dataset, 'specimenMasks.mat'));

groups = findgroups([ID.MsiID]);
options.saveOptions.saveImages = true;
options.saveOptions.saveInHQ = false;

for g = 1:max(groups)
%% Read MSI
    gIdxs = find(groups == g);
    gMembers = ID(gIdxs);
    coordinates = [[gMembers.Originx]; [gMembers.Originy]]';
    idd = gMembers(1);
    files = {data([data.MsiID] == idd.MsiID).File};
    [~, whiteReference, ~] = readMSI(files); 
    baseImage = whiteReference .* specimenMasks{g};
    [m,n,~] = size(baseImage);
    labels = ~[gMembers.IsBenign];
    options.saveOptions.plotName = fullfile('..','..','..', 'output', ...
        dataset, '2-Labels', strcat('labelled_bright_', num2str(g)));
    plotVisualResult(baseImage, zeros(m,n), '', labels, coordinates, 'jet', true, 2, options.saveOptions);
%     for i = 1:length(coordinates)
%         plotVisualResult(whiteReference, zeros(size(whiteReference,1),...
%             size(whiteReference,2)), num2str(idd.Index + i), labels(i), coordinates(i,:), 'jet', true, 3, options.saveOptions);
%     end
end