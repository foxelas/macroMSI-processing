function [style] = getStyleName(labels)

if sum(contains(labels, 'Benign')) > 1
    style = 'class';
elseif sum(contains(labels, 'Unfixed')) > 1
    style = 'type';
else
    style = 'sample';
end
end
