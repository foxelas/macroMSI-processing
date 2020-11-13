function [solaxSpec] = reconstructSolaxIoIlluminationSpectrum(savedir)
%savedir = "D:\temp\Google Drive\titech\research\experiments\output\5. Progress Reports\img";

    x = [350 360 370 380 400 410 425 450 460 480 520 525 530 575 620 625 630 640 660 690 695 700 725 750 800]';
    y = [0 0 0 0 20 63 40 82 80 58 80.5 79.5 80 68 73 72 73 70 60 39.5 40 30 15 9 0]';

    f1=fit(x,y,'linearinterp');
    options = fitoptions('Method','SmoothingSpline',...
                         'SmoothingParam',0.1);
    f2=fit(x,y,'smoothingspline', options);

    fig1 = figure(1);
    plot(x,y,'x','MarkerEdgeColor','black')
    hold on
    %plot(x, v)
    plot(f1,x,y)
    hold off
    grid on;
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel('Relative Spectrum (%)', 'FontSize', 15)
    title('Using [linearinterp] fitting', 'FontSize', 15)

    fig2 = figure(2);
    plot(x,y,'x','MarkerEdgeColor','black')
    hold on
    %plot(x, v)
    plot(f2,x,y)
    hold off
    grid on;
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel('Relative Spectrum (%)', 'FontSize', 15)
    title('Using [smoothingspline] fitting', 'FontSize', 15)

    setSetting('cropBorders', true);
    setSetting('plotName', fullfile(savedir, 'solaxSpectrum_lininterp.png'));
    savePlot(fig1);

    setSetting('plotName', fullfile(savedir,'solaxSpectrum_spline.png'));
    savePlot(fig2);

    wavelengths = [380:780]';
    solaxSpec = f2(wavelengths);
    solaxSpec = max(solaxSpec, 0);
    save('parameters/solax_reconstructed_spectrum.mat', 'solaxSpec');

end 
