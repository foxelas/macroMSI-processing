%% Reconstruct sRGB from MSI

for k = 20:max([ID.Group])
    infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(k), '.mat'));
    load(infile, 'raw', 'whiteReference');
    figure(5); imshow(whiteReference);
    z = find([ID.Group] == k, 1);
    sRGB = createSRGB(raw, 'original', ID(z), options);
end

% end of create srgb