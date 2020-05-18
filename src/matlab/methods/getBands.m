function bands = getBands(method)

switch method
    case 'extended'
        bands = {'450B', '465B', '465G', '505G', '525G', '575G', '605G', '605R', '630R'};
    case 'adjusted'
    case 'green'
    case 'rms'
    case 'max'
        bands = {'450', '465', '505', '525', '575', '605', '630'};
end

end