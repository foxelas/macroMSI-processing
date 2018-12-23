function [ ] = showCroppedSection( image, imagePath, x, y, figtitle,  savedir, toClean, color, boxdim )

        if (nargin < 5) || isempty(figtitle)
            figtitle = '';
        end
        
        if (nargin < 6) || isempty(savedir)
            savedir = ''; 
            toSave = false;
        else 
            toSave = true;
        end
        
        if (nargin < 7)
            toClean = true;
        end
        
        if (nargin < 8)
            color = 'y';
        end
        
        if (nargin < 9)
            height = 5;
            width = 5;
        else
            height = boxdim(1);
            width = boxdim(2);
        end
        
        %to check that we read the correct region
        if toClean
            fig = figure(1);
            clf(fig);
        else 
            fig = gcf;
        end
        
        w = warning ('off','all');
        
        hold on
        
        if toClean
            if isempty(image)
                imshow(imread(imagePath), 'Border', 'tight');
            else
                imshow(image, 'Border', 'tight');
            end
        end
        
        rectangle('Position', [x, y, width, height], 'EdgeColor',color);
        title(figtitle)
        if toSave
            print(fig, savedir, '-djpeg');
        end
        
        hold off

        warning(w)

        pause(0.01);
   

end

