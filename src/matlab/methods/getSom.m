function [net, outputs] = getSom(inputs, dimension1, dimension2, inputName)
%     GETSOM returns som topology 
% 
%     Usage: 
%     [~, outputs] = getSom(inputs, dimension1, dimension2, 'spect');

net = selforgmap([dimension1, dimension2]);

% Train the Network
[net, ~] = train(net, inputs);

% Test the Network
outputs = net(inputs);

% View the Network
view(net)
fig1 = figure(1);
plotsomnd(net);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('som'), ...
    strcat(inputName, '_', num2str(dimension1), 'x', num2str(dimension2), '_somd')));
savePlot(fig1);

fig2 = figure(2);
plotsomhits(net, inputs);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('som'), ...
    strcat(inputName, '_', num2str(dimension1), 'x', num2str(dimension2), '_somh')));
savePlot(fig2);

end
