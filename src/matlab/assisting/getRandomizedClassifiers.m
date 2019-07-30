fid = fopen('randomized.txt');
%count = 829
tline = fgetl(fid);
isNextLine = false;
count = 0;
while ischar(tline)
    if isNextLine 
        settings = strrep(tline, '{', "");
        settings = strrep(settings, '}', "");
        settings = strrep(settings, ", '", ", ");
        settings = strrep(settings, "':", ":");
        settings = strrep(settings, ":", "=");
        if contains(settings, "class_weight") && ~contains(settings, "balanced")
            settings = strrep(settings, "class_weight=", "class_weight={");
            settings = strrep(settings, "0= ", "0: ");
            settings = strrep(settings, "1= 1.0", "1: 1.0}");            
        end
        classifierSetting = strcat(classChar, "(", settings, additionalParams, ")");
        classifierSetting = strrep(classifierSetting, "('", "(");
        classifierSetting = strcat(classifierSetting, ",");
        disp(classifierSetting)
        count = count + 1;        
    end
    if contains(tline, "Best params for")
        isNextLine = true;
        if contains(tline, "Random Forest")
            classChar = "RandomForestClassifier";
            additionalParams = ", random_state=2, n_jobs = -1";
        elseif contains(tline, "SVM")
            classChar = "SVC";
            additionalParams = ", probability = True, random_state=2";
        elseif contains(tline, "KNN")
            classChar = "KNN";
            additionalParams = "";
        end
    else 
        isNextLine = false;
    end      
    tline = fgetl(fid);
end
fclose(fid);