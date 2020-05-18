function outputFolderMap = getOutputDirectoryMap()

    key = {'backgroundRemoval', 'labels', 'segments', ...
        'segmentsForFeatureExtraction', 'reflectanceEstimation',...
        'reflectanceEstimationPerformance', 'lbpVisualization', ...
        'features', 'classifierPerformance', 'classifierPerformance-rgb',...
        'sRGB', 'visualTool', 'classifierPerformanceOnlySpect', 'som', 'dimred'};
    value = {'1-BackgroundRemoval', '2-Labels', '3-Segments', ...
        '4-SegmentsForFeatureExtraction', '5-ReflectanceEstimation',...
        '6-ReflectanceEstimationPerformance', '7-LBPVisualisation', ...
        '8-Features', '9-ClassifierPerformance', ...
        '9-ClassifierPerformance_rgb', '10-sRGB', '11-Visual Tool', ...
        '10-ClassifierPerformanceOnlySpect', '12-Som', '13-DimensionReduction'};
    
    outputFolderMap = containers.Map(key, value);   
end 
