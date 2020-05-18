function [net, outputs] = getSom(inputs, dimension1,dimension2, saveOptions, inputName)

    outputFolderMap = getOutputDirectoryMap();

    net = selforgmap([dimension1 dimension2]);

    % Train the Network
    [net,tr] = train(net,inputs);

    % Test the Network
    outputs = net(inputs);

    % View the Network
    view(net)
    fig1 = figure(1);
    plotsomnd(net);
    saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('som'), ...
        strcat( inputName, '_', num2str(dimension1), 'x',num2str(dimension2), '_somd'));
    savePlot(fig1, saveOptions);  
    
    fig2 = figure(2);
    plotsomhits(net,inputs);
        saveOptions.plotName = fullfile(saveOptions.savedir, outputFolderMap('som'), ...
        strcat( inputName, '_', num2str(dimension1), 'x',num2str(dimension2), '_somh'));
    savePlot(fig2, saveOptions);  

end

