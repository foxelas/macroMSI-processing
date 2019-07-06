function [] = plotReconstructionGIF(measured, reconstructions)

    w = 380:5:780;
    msiN = size(measured, 1);
    for i = 1:msiN
        figure(1);
        clf;
        hold on
        plot(w, measured(i,:), 'b', 'DisplayName', 'Measured', 'LineWidth', 2);
        plot(w, squeeze(reconstructions(i,114,:))', 'r', 'DisplayName', 'Simple', 'LineWidth', 2);
        plot(w, squeeze(reconstructions(i,3,:)), 'm', 'DisplayName', 'Spatiospectral', 'LineWidth', 2);
        plot(w, squeeze(reconstructions(i,792,:)), 'g', 'DisplayName', 'Spatial', 'LineWidth', 2);
        legend('Location', 'northwest', 'FontSize', 15)
        hold off
        xlabel('Wavelength')
        ylabel('Reflectance')
        pause(0.01)    
    end
end

