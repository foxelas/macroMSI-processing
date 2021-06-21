function [] = plotClassificationErrors(errors, fig)

n = min(20, length(errors));
accuracy = [errors(1:n).Accuracy];
typeI = [errors(1:n).TypeI];
typeII = [errors(1:n).TypeII];
hold on
scatter(1:n, accuracy, 'rx')
scatter(1:n, typeI, 'bo')
scatter(1:n, typeII, 'mo')
text(1:n, accuracy-2, arrayfun(@(x) sprintf('%.2f%%', x), accuracy, 'UniformOutput', false), 'FontSize', 8);
text(1:n, typeI+2, arrayfun(@(x) sprintf('%.2f%%', x), typeI, 'UniformOutput', false), 'FontSize', 8);
text(1:n, typeII-2, arrayfun(@(x) sprintf('%.2f%%', x), typeII, 'UniformOutput', false), 'FontSize', 8);

hold off
xlim([0, n + 1])
grid on
xticks(1:n)
ll = cell(n, 1);
for i = 1:n
    ll{i} = strcat(errors(i).Input, '|', errors(i).Projection, '|', num2str(errors(i).Neighbours), '-', errors(i).VoteRule, '|', errors(i).Distance);
end
xticklabels(ll)
xtickangle(45)
legend('Accuracy', 'False positive rate', 'False negative rate', 'Location', 'best')
ylabel('Classification Accuracy')
title(strcat('Classification results (', errors(i).Validation, ' validation)'));

savePlot(fig);

end
