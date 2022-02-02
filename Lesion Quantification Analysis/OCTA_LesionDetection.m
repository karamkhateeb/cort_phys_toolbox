%% OCTA Lesion Detection
% This script is used to detect and measure OCTA lesion sizes
% Author: Karam Khateeb

%% 1.Store Post-Photothrombosis Stitched Lesion Images

imCell = cell(1,5);
imCell{1,1} = imread('PT2_postOCTAL.tif');
imCell{1,2} = imread('PT2_postOCTAR.tif');
imCell{1,3} = imread('PT3_postOCTAL_BGFixed.tif'); imCell3 = imread('PT3_postOCTAL.tif');
imCell{1,4} = imread('PT3_postOCTAR.tif');
imCell{1,5} = imread('PT4_postOCTAL.tiff'); imCell{1,5} = flip(imCell{1,5},2);
lesions = cell(1,5);

lesions{1}.im = imCell{1,1};
lesions{1}.ID = 'PT2_postOCTAL.tif';

lesions{2}.im = imCell{1,2};
lesions{2}.ID = 'PT2_postOCTAR.tif';

lesions{3}.im = imCell{1,3};
lesions{3}.ID = 'PT3_postOCTAL.tif and PT3_postOCTAL_BGFixed.tif';

lesions{4}.im = imCell{1,4};
lesions{4}.ID = 'PT3_postOCTAR.tif';

lesions{5}.im = imCell{1,5};
lesions{5}.ID = 'PT4_postOCTAL.tiff';
%% 2. Filter and Detect Lesion Locations and Sizes
for i = 1:length(imCell)
    currentIm = imCell{1,i};
    im = currentIm;

    if i == 1 || i == 2 || i == 4
        % increase contrast
        edgeThreshold = 1; amount = 1;
        im = localcontrast(im, edgeThreshold, amount);
        figure; imshow(im);

        % gaussian filter
        im = imgaussfilt(im, 3);
        imshow(im);
    elseif i == 3
        % gaussian filter
        im = imgaussfilt(im, 7);
        figure; imshow(im);
        
        % gaussian filter
        imCell3 = imgaussfilt(imCell3, 7);
        imshow(imCell3);
    elseif i == 5
        im = imgaussfilt(im,60);
        figure, imshow(im);
    end
    
    % binarize the images
    [counts, ~] = imhist(rgb2gray(imCell3), 256);
    if i == 5
        [counts,~] = imhist(rgb2gray(im), 256);
    end
    T = otsuthresh(counts);
    bw = imbinarize(rgb2gray(im), T);
    bw = imcomplement(bw);
    imshow(bw);

    % obtain stats about each lesion
    stats = regionprops('table',bw,'Centroid',...
        'MajorAxisLength','MinorAxisLength')

    centers = stats.Centroid(2:end,:);
    diameters = mean([stats.MajorAxisLength(2:end) stats.MinorAxisLength(2:end)],2);
    maxDs = stats.MajorAxisLength(2:end);
    minDs = stats.MinorAxisLength(2:end);
    radii = diameters/2;

    % visualize the lesion
    imshow(currentIm);hold on
    viscircles(centers,radii,'linewidth',1);
    hold off

    % Clean boundariesYL by removing registration artifacts, dust, noise,etc
    lesions{i}.Centers = []; lesions{i}.Radii = [];
    lesions{i}.majD = []; lesions{i}.minD = [];
    counter = 0;
    choosing = 1;
    disp('Choose a point along a boundary to keep it, press enter when done');
    [x, y] = ginput;
    for iSelxns = 1:length(x)
        located = 0;
        for iCircle = 1:numel(centers)/2
            circleX = centers(iCircle,1) + radii(iCircle)*cos(linspace(0,2*pi,500));
            circleY = centers(iCircle,2) + radii(iCircle)*sin(linspace(0,2*pi,500));
            for iPts = 1:numel(circleX)
                if abs(x(iSelxns)-circleX(iPts)) <= 15 && abs(y(iSelxns)-circleY(iPts)) <= 15 && ~located
                    located = 1;
                    counter = counter + 1;
                    lesions{i}.Centers = [lesions{i}.Centers; centers(iCircle,:)];
                    lesions{i}.Radii = [lesions{i}.Radii; radii(iCircle)];
                    lesions{i}.majD = [lesions{i}.majD; maxDs(iCircle)];
                    lesions{i}.minD = [lesions{i}.minD; minDs(iCircle)];
                end
            end
        end
    end
    
    imshow(currentIm);hold on
    viscircles(lesions{i}.Centers,lesions{i}.Radii,'linewidth',1);
    
    switch i
        case 1, ID = 'B_L';
        case 2, ID = 'B_R';
        case 3, ID = 'C_L';
        case 4, ID = 'C_R';
        case 5, ID = 'D_L';
    end
    
    imwrite(bw,['Binary' ID '.tif']);
    hold off
end
  

%% 3. Engrave lesions in to images to Convert to phyiscal measurements
% need to save the images with the circles on them as is to preserve
% resolution

% Draw each circle into the image itself in white
for i = 1:numel(lesions)
    lesions{i}.engraved = lesions{i}.im;
    for iCircle = 1:numel(lesions{i}.Radii)
        if i == 5
            circleRes = 1999;
        else
            circleRes = 1000;
        end
        circleX = lesions{i}.Centers(iCircle,1) + lesions{i}.Radii(iCircle)*cos(linspace(0,2*pi,circleRes));
        circleY = lesions{i}.Centers(iCircle,2) + lesions{i}.Radii(iCircle)*sin(linspace(0,2*pi,circleRes));
        
        indices = find(circleX > 0 );%& circleY > 0);
        circleX = circleX(indices);
        circleY = circleY(indices);
        
        for iPt = 1:numel(circleX)
            lesions{i}.engraved(round(circleY(iPt)),round(circleX(iPt)),:) = [255 255 255];
        end
    end
    
    figure, imshow(lesions{i}.engraved);
    
    filename = ['Engraved_' lesions{i}.ID];
    imwrite(lesions{i}.engraved,filename);
    
end

%% 4. Convert Everything after looking at the images in Adobe Illustrator to mm

resoln(1) = 21.5; % pixels/mm as measured in Illustrator
resoln(2) = 21.5;
resoln(3) = 21.5;
resoln(4) = 21.5;
resoln(5) = 1999/9; % FOV of single image is ~ 9 mm x 9 mm

for i = 1:numel(lesions)
    lesions{i}.mmRadii = lesions{i}.Radii/resoln(i);
    lesions{i}.mmMajD = lesions{i}.majD/resoln(i);
    lesions{i}.mmMinD = lesions{i}.minD/resoln(i);
end

%% 5. Save Lesion Indices
% here is a psuedo map of the expected lesions with their diameters in mm:
% - 2 -
% 2 - 1
% - 1 -
% 1 - .5
% - .5 -

% here is their corresponding indices for reference:
% - 3 -
% 1 - 6
% - 4 -
% 2 - 7
% - 5 -

lesions{1}.lesIDs = [1;4;5;6;7];
lesions{2}.lesIDs = [2;4;7];
lesions{3}.lesIDs = [1;2;3;4];
lesions{4}.lesIDs = [1;2;4];
lesions{5}.lesIDs = 8;

% identify which lesions that should have appeared, but did not show up
% these are different that the ones that were likely not in the FOV
lesions{1}.noIDs = 2;
lesions{2}.noIDs = [1;3;6];
lesions{3}.noIDs = [5;6;7];
lesions{4}.noIDs = [3;5;7];
lesions{5}.noIDs = [];
%% 6. Save Everything
 
save('OCTALesionInfo','lesions');
