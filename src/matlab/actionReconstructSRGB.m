%% Reconstruct sRGB from MSI
close all;
for k = 1:max([ID.Group])
    infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(k), '.mat'));
    load(infile, 'raw', 'whiteReference', 'specimenMask');
    figure(5); imshow(whiteReference);
    z = find([ID.Group] == k, 1);
    sRGB = createSRGB(raw, 'medium', ID(z), options, 'cmccat2000', specimenMask);
    outfile = fullfile(options.saveOptions.savedir, '10-sRGB', strcat('group_', num2str(k), '.mat'));
    save(outfile, 'sRGB');
end

% end of create srgb