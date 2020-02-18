%for k = 1:max([ID.Group])
k = 1;
infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(k), '.mat'));
load(infile, 'raw', 'whiteReference', 'specimenMask');
[channels, height, width, ~] = size(raw);
figure(5);
imshow(whiteReference);
z = find([ID.Group] == k, 1);

%sRGB = createSRGB(raw, 'medium', ID(z), options, 'cmccat2000', specimenMask);
%outfile = fullfile(options.saveOptions.savedir, outputFolderMap('sRGB'), strcat('group_', num2str(k), '.mat'));
%save(outfile, 'sRGB');
%end

reference  = raw2msi(white, 'max');

msi = raw2msi(raw, 'max');
foregroundMask = permute(repmat(double(specimenMask), 1, 1,  channels), [3 1 2]);
msi = bsxfun(@times, msi, foregroundMask);
cluster1 = whiteReference .* double(specimenMask);

[B, M, N] = size(msi);
opticalDensity = log(msi ./ reference);

plotMSI(msi, 1, options.saveOptions);
plotMSI(opticalDensity, 2, options.saveOptions);