%for k = 1:max([ID.Group])
k = 17;
msiType = 'max'; %'extended'; % 'max';
[msi, whiteReference, specimenMask, height, width, channels] = getImage(k, msiType, false);
systemdir = getSetting('systemDir');
reference = getReference(systemdir, height, width);
reference = raw2msi(reference, msiType);

%sRGB = createSRGB(raw, 'medium', ID(z), 'cmccat2000', specimenMask);
%outfile = fullfile(getSetting('savedir'), getSetting('sRGB'), strcat('group_', num2str(k), '.mat'));
%save(outfile, 'sRGB');
%end


%foregroundMask = permute(repmat(double(specimenMask), 1, 1,  channels), [3 1 2]);
%msi = bsxfun(@times, msi, foregroundMask);
opticalDensity = double(log10(msi./reference));

plots(1, @plotMSI, msi);
plots(2, @plotMSI, reference);
plots(3, @plotMSI, opticalDensity);

figure(4);
od630 = squeeze(opticalDensity(7, :, :));
subplot(3, 1, 1);
c = imagesc(od630);
c.Parent.Visible = 'off';
colorbar;
title('OD630nm - Melanin Map ')

odhg = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
subplot(3, 1, 2);
c = imagesc(odhg);
c.Parent.Visible = 'off';
colorbar;
title('OD575nm - 1.15 OD630nm - Hemoglobin Map ')

subplot(3, 1, 3);
od575 = squeeze(opticalDensity(5, :, :));
c = imagesc(od575);
c.Parent.Visible = 'off';
colorbar;
title('OD575nm')
