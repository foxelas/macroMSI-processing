for i = 113:8:length(data2)
    [msi, white] = readMSI({data2(i:(i+8-1)).File},[], [], [], [], [1, 450, 465, 505, 525, 575, 605, 630 ], true );
    imshow(white);%imshow(squeeze(msi(1,:,:,:)));
    title(i)
end