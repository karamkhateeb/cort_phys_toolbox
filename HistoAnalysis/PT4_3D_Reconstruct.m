%% Histology Image Analysis for PT4 - 3D Lesion Reconstruction
% Author: Karam Khateeb and Mona Rahimi
% The purpose of this script is to reconstruct the lesions in 3D space
% based on the Nissl staining

% Images after rescaleing (0.075):
% 187pixels/cm, 53.5 microns per pixel, rounding to 50 microns. -determined from RULER062 
% Slices" 50 micron thick, 450 microns between slices
% for scaling, expand slices 10x

% Slide nomenclature: slide# (anterior -> posterior)
% boxes = 1:4;
% wells = [1 10 11 19 28 29];

slides = 1:24;
resizeFactor = 0.1;

filepath = [filesep 'Users' filesep 'Mona Rahimi' filesep 'Documents' filesep 'Yazdan Lab' filesep 'PT4'];
%% Register the Images - Crop and Rotate to Align Images

% First need a starter 'fixed' image -> start with 1.1.1
tic
% load the first 'fixed' image
fixedIm = imread(['PT4_Working_Slides' filesep '1.tif']);

% pick the area of interest and crop
cropRect = [2515 30 7096 6330]; % ROI determined manually [XYWH]

croppedFixedIm = imcrop(fixedIm, cropRect);
figure; imshow(croppedFixedIm);
croppedFixedIm = imresize(croppedFixedIm, resizeFactor);

% rotate for the left hemisphere
degreeRot = -5; % Determined manually
fixedIm = imrotate(croppedFixedIm, degreeRot); 
figure; imshow(fixedIm); title('Fixed Image');

%% Obtain the optimizer 
[optimizer, metric] = imregconfig('multimodal');

% Index to base registration off of every other image
register = 0;

for islide = 1:24
            tic
            % Load the 'moving' image
            % Checks if file exists, skips if DNE
            currentFile = ['PT4_Working_Slides' filesep num2str(islide) '.tif'];
            if exist(currentFile, 'file') == 0
                continue;
            end 
            currentIm = imread(['PT4_Working_Slides' filesep num2str(islide) '.tif']);
                          
            % pick the area of interest and crop
            cropRect = [2515 30 7096 6330]; % ROI determined manually [XYWH]
            croppedIm = imcrop(currentIm, cropRect);
            
            % Resize image for troubleshooting
            croppedIm = imresize(croppedIm, resizeFactor); 
                                    
            if (islide < 18 && islide > 11)
                [mp, fp] = cpselect(croppedIm, fixedIm, 'Wait', true);
                t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
                Rfixed = imref2d(size(fixedIm));
                registeredIm = imwarp(croppedIm, t, 'OutputView', Rfixed);
                % imshowpair(fixedIm, registeredIm, 'blend'); % for troubleshooting
            else  
                % figure; imshow(croppedIm); title('croppedIm') % for troubleshooting
                tform = imregtform(rgb2gray(croppedIm), rgb2gray(fixedIm), 'rigid', optimizer, metric);
                registeredIm = imwarp(croppedIm, tform, 'OutputView', imref2d(size(fixedIm)));
            end 
            
            
            % Grabs white patch
            patch = imcrop(croppedIm, [19 39 27 25]);  
            patchRed = patch(:,:,1);
            patchGreen = patch(:,:,2);
            patchBlue = patch(:,:,3);
            
            meanRed = mean(patchRed, 'all');
            meanGreen = mean(patchGreen, 'all');
            meanBlue = mean(patchBlue, 'all');
            
            % figure; imshow(registeredIm); title('registeredIm after patch');
            
            % Replace any black space with white space so that it doesn't
            % mess up feature identification and alignment
            redChannel = registeredIm(:,:,1);
            greenChannel = registeredIm(:,:,2);
            blueChannel = registeredIm(:,:,3);
            blackPixels = redChannel == 0 & greenChannel == 0 & blueChannel == 0;
            redChannel(blackPixels) = meanRed;
            greenChannel(blackPixels) = meanGreen;
            blueChannel(blackPixels) = meanBlue; 

            registeredIm = cat(3, redChannel, greenChannel, blueChannel);
            figure; imshow(registeredIm); title('registeredIm after replace black space');

            topCrop = imcrop(registeredIm, [1 1 760 421]);
            
            figure; imshow(topCrop); title('topCrop image');
            
            % resave the image in PT5_RegSlides folder
            imwrite(topCrop, [filepath filesep 'PT4_RegSlides' filesep num2str(islide) '.tif']);
   
            disp(['Saved registered ' num2str(islide) '.tif']);
           
            % Switch the fixed image to registered image every fourth image
            switch register 
                case 0
                    disp('Not switching fixedIm');
                    register = 1;
%                case 1
%                   disp('Not switching fixedIm');
%                   register = 2;
%                 case 2
%                     disp('Not switching fixedIm');
%                     register = 3;
                case 1
                    disp('SWITCH');
                    fixedIm = registeredIm; % Set the next fixedIm to be the previous registeredIm
                    register = 0;
              end
end
toc
disp('Done Saving'); 

% Pause here to enhance borders in Photoshop

 %% Isolate Brain Outline with and without Lesions (MATLAB Contrast Approach --> FAIL) 
% % edgeThreshold must be in [0,1]
% % amount must be in [-1,1]
% close all;
% currentIm = imread(['PT4_RegSlides' filesep '10.tif']);
% edgeThreshold = 1; amount = 1;
% contrIm = localcontrast(currentIm, edgeThreshold, amount);
% 
% edgeThreshold = 0.5; amount = -1;
% contrIm = localcontrast(contrIm, edgeThreshold, amount);
% 
% figure; imshow(contrIm);
% saveas(gcf, [filepath filesep 'PT4_FinalExclude' filesep '10.tif']); % Lesion outlined
% 
% edgeThreshold = 0.8; amount = 0.6;
% contrIm = localcontrast(currentIm, edgeThreshold, amount);
% 
% edgeThreshold = 0.3; amount = -0.8;
% contrIm = localcontrast(contrIm, edgeThreshold, amount);
% figure; imshow(contrIm);
% saveas(gcf, [filepath filesep 'PT4_FinalInclude' filesep '10.tif']); % Lesion blurred

%% Apply Filter to Isolate Brain Outline with and without Lesions
% Want to apply the same filtering schemes for both types of boundaries
% (with and without lesions) so that the volume comparisons can be more
% accurate and comparable

% Reset the for loops
% Slide nomenclature: slide# anterior -> posterior

slideIdx = 0;
fullBrain = {};
tic
for islide = 1:24
            % Load registered image
            currentIm = imread(['PT4_RegSlides' filesep num2str(islide) '.tif']);  
            noLesion = imread(['PT4_FinalInclude' filesep num2str(islide) '.tif']);
            yesLesion = imread(['PT4_FinalExclude' filesep num2str(islide) '.tif']);
            
            disp(['Reading ' num2str(islide) '.tif']);
            
%             % Increase contrast of the image
% %             edgeThreshold = 0.1; amount = 1;
% %             contrIm = localcontrast(currentIm, edgeThreshold, amount);
% %             figure; imshow(contrIm);
%             
% %             doubleIm = double(contrIm);
%             
% %             reshapedIm = reshape(doubleIm, size(doubleIm,1)*size(doubleIm,2),3);
%             
% %             coeff = pca(reshapedIm);
%             
% %             imTransformed = reshapedIm * coeff;
%             
% %             imPC1 = reshape(imTransformed(:,1), size(doubleIm,1), size(doubleIm, 2));
% %             imPC2 = reshape(imTransformed(:,2), size(doubleIm,1), size(doubleIm, 2));
% %             imPC3 = reshape(imTransformed(:,3), size(doubleIm,1), size(doubleIm, 2));
%             
% %             figure, imshow(imPC1,[]);
% %             figure, imshow(imPC2,[]);
% %             figure, imshow(imPC3,[]);
%             
%             
%             % Determine noise patch (non-brain tissue area of slide) to get variance
%             patch = imcrop(currentIm, [120 100 220 200]);
%             patchVar = std2(patch)^2;
%             
%             % Apply edge-preserving bilateral filter to smooth entire brain slice
%             DoS = 10*patchVar; % degree of smoothing greater than patch variance
%             spatialSigma = 50;
%             smoNL = imbilatfilt(noLesion, DoS, spatialSigma); % apply filter
%             smoYL = imbilatfilt(yesLesion, DoS, spatialSigma);
% %             figure; imshow(smoIm); title(['Smoothed ' num2str(ibox) '.' num2str(iwell) '.' num2str(ipass)]);
%             
% %             edgeThreshold = 0.4; amount = 0.5;
% %             smoIm = localcontrast(currentIm, edgeThreshold, amount);
% 
% %             % Apply Gaussian filter to image
% % %             smoIm = imgaussfilt(smoIm, 5);
%             
%             edgeThreshold = 0.9; amount = 0.9;
%             contrNL = localcontrast(smoNL, edgeThreshold, amount);
%             contrYL = localcontrast(smoYL, edgeThreshold, amount);
% 
%             % Apply Gaussian filter to image
%             gaussNL = imgaussfilt(contrNL, 3);
%             gaussYL = imgaussfilt(contrYL, 3);
% %             figure; imshow(smoIm);
%             
%             % do PCA on smoothed imaged
% %             doubleIm = double(smoIm);            
% %             reshapedIm = reshape(doubleIm, size(doubleIm,1)*size(doubleIm,2),3);
% %             coeff = pca(reshapedIm);
% %             imTransformed = reshapedIm * coeff;
% %             imPC1 = reshape(imTransformed(:,1), size(doubleIm,1), size(doubleIm, 2));
% %             figure; imshow(imPC1, []);
           
            % Use the same threshold for both NL and YL images to binarize
            [counts, x] = imhist(rgb2gray(noLesion), 256);
            T = otsuthresh(counts);
            
            % Binarize the smoothed image
            BW_NL = imbinarize(rgb2gray(noLesion), T);
            BW_YL = imbinarize(rgb2gray(yesLesion), T);
%             figure; imshow(BW);
            
            % Compute the boundaries
            boundariesNL = bwboundaries(BW_NL);
            boundariesYL = bwboundaries(BW_YL);
%             figure; imshow(BW, []); hold on; visboundaries(boundaries);
            
            figure; imshow(noLesion); hold on; visboundaries(boundariesNL);
            
            % Clean boundariesNL by removing registration artifacts, dust, noise,etc
            cleanBoundariesNL = {};
            counterNL = 0;                      
            choosing = 1;
            boundarySelected = 0;
            while choosing && ~boundarySelected
                 
                disp('Choose a point along a boundary to keep it');
                [x, y] = ginput(1);
                
                for iBound = 1:numel(boundariesNL) 
                    for iPts = 1:numel(boundariesNL{iBound})/2
                        if abs(x-boundariesNL{iBound}(iPts,2)) <= 15 && abs(y-boundariesNL{iBound}(iPts,1)) <= 15 && ~boundarySelected
                            boundarySelected = 1;
                            counterNL = counterNL + 1;
                            cleanBoundariesNL{counterNL} = boundariesNL{iBound,1}; 
                            continue
                        elseif boundarySelected
                            continue
                        end
                    end
                    if boundarySelected
                        continue
                    end
                end
                
                if ~boundarySelected
                    %                     prompt = ['Any more boundaries to select? '];
                    %                     choosing = input(prompt);
                    choosing = 0;
                    %                     boundarySelected = ~choosing;
                    %                 else
                    disp(['No boundary selected, try again']);
                end
            end
            close gcf;
            
            
            figure; imshow(yesLesion); hold on; visboundaries(boundariesYL);
            
            % Clean boundariesYL by removing registration artifacts, dust, noise,etc
            cleanBoundariesYL = {};
            counterYL = 0;                      
            choosing = 1;
            boundarySelected = 0;
            while choosing && ~boundarySelected
                disp('Choose a point along a boundary to keep it');
                [x, y] = ginput(1);
                for iBound = 1:numel(boundariesYL) 
                    for iPts = 1:numel(boundariesYL{iBound})/2
                        if abs(x-boundariesYL{iBound}(iPts,2)) <= 15 && abs(y-boundariesYL{iBound}(iPts,1)) <= 15 && ~boundarySelected
                            boundarySelected = 1;
                            counterYL = counterYL + 1;
                            cleanBoundariesYL{counterYL} = boundariesYL{iBound,1}; 
                            continue
                        elseif boundarySelected
                            continue
                        end
                    end
                    if boundarySelected
                        continue
                    end
                end
                
                if ~boundarySelected
                    disp(['No boundary selected, try again']);
                end
            end
            close gcf
            
            
            figure; imshow(currentIm); hold on; 
            visboundaries(cleanBoundariesNL); set(gca, 'ydir', 'reverse');
            visboundaries(cleanBoundariesYL); set(gca, 'ydir', 'reverse');
            
            saveas(gcf, [filepath filesep 'PT4_BrainBounds' filesep num2str(islide) '.tif']);
            slideIdx = slideIdx +1;
            fullBrain{slideIdx,1} = cleanBoundariesNL{1,1};
            fullBrain{slideIdx,2} = cleanBoundariesYL{1,1};
            toc
            close gcf
end

save('PT4_LesionBounds.mat', 'fullBrain', 'currentIm');

%% Prepare Slices into matrices for visualization
% Load in  fullBrain + image for width + height
% load(fullBrain.mat) or something
% currentIm = imread(['PT3_RegSlides' filesep num2str(ibox) '.'
% num2str(iwell) '.' num2str(ipass) '.tif']); This doesnt matter which
% image we just need to define the size of our matrix

dims = [size(currentIm,1), size(currentIm,2), numel(fullBrain)/2];
brain = zeros(dims);
lesions = zeros(dims);
lesionEdgeMatrix = zeros(dims);

% Fill in the outline for each slide
tic
for i = 1:numel(fullBrain)/2
    % Set up arrays for the slices to go into
    brainSlice = zeros([size(currentIm,1), size(currentIm,2)]);
    lesionSlice = zeros([size(currentIm,1), size(currentIm,2)]);
    % Trace outline to array for full brain
    for j = 1:length(fullBrain{i,1})
        %convert cells to pixels
        pixel = sub2ind(size(brainSlice),fullBrain{i,1}(j,1),fullBrain{i,1}(j,2));
        brainSlice(pixel) = 1;
    end
    
    %trace outline to arry for missing lesion
    for j = 1:length(fullBrain{i,2})
        pixel = sub2ind(size(brainSlice),fullBrain{i,2}(j,1),fullBrain{i,2}(j,2));
        lesionSlice(pixel) = 1;
    end
    
    % Fill in outline
    brainSlice = imfill(brainSlice,4,'holes');
%     imshow(brainSlice)
    % isolate lesions
    lesionSlice = imfill(lesionSlice,4,'holes');
    brain(:,:,i) = lesionSlice;
%     imshow(lesionSlice)
    % Find areas missing from full brain
    lesion = brainSlice - lesionSlice;
    lesion(lesion < 0) = 0;
    
    % remove noise
    lesionObj = bwconncomp(lesion, 4);
    lesionPixels = cellfun(@numel, lesionObj.PixelIdxList);
    lesionIdx = lesionPixels > 75;
    
    for j = 1:lesionObj.NumObjects
        if lesionIdx(j) == 0
            lesion(lesionObj.PixelIdxList{j}) = 0;
        end
    end
%     imshow(lesion)
%     pause(1)
    lesions(:,:,i) = lesion;
    lesionEdges{i} = bwboundaries(lesion);
    toc
end

% Get the outline of the lesions only for interpolation
% figure;
% for iSlice = 1:length(lesionEdges)
%     % Set up arrays for the slices to go into
%     lesionEdgeArray = zeros([size(currentIm,1), size(currentIm,2)]);
%     
%     % Trace outline to array for lesions
%     if ~isempty(lesionEdges{iSlice})
%         for j = 1:length(lesionEdges{iSlice,1}{1,1})
%             %convert cells to pixels if not empty
%             pixel = sub2ind(size(lesionEdgeArray),lesionEdges{iSlice,1}{1,1}(j,1),lesionEdges{iSlice,1}{1,1}(j,2));
%             lesionEdgeArray(pixel) = 1;
%         end
%     end
%     
%     imshow(lesionEdgeArray); pause(1);
%     
%     lesionEdgeMatrix(:,:,iSlice) = lesionEdgeArray;
% end

lesions = flip(lesions, 3);
save('PT4_LesionBounds.mat', 'lesions', '-append');
    
%% Visualize whole brain + lesions

% Smooth both surfaces
smoLesions = smooth3(lesions, 'gaussian', 5, 10);
smoBrain = smooth3(brain, 'gaussian', 5, 10);

% Now you need to scale it
% Set the scale in the AP axis (450um between each slice)

% The voxel volume
delta = 0.079; % in mm
dx = delta; dy = delta; dz = 0.5;
dV = dx*dy*dz; % Voxel volume in mm^3

save('PT4_LesionBounds.mat', 'dx', 'dy', 'dz','dV', '-append');

APScale = [0 0.45 0.90 1.35 2.25 2.70 3.15 3.6 4.05 4.50 4.95 5.85 6.30 6.75 7.70 8.1 8.60 9.45 9.90 10.35 11.25 11.7 12.15 12.2]; % in mm 

% Set the scale in the ML and DV axis (each pixel is 52.43 um per pixel)

MLScale = 0:delta:delta*(size(brain,1)-1);
DVScale = 0:delta:delta*(size(brain,2)-1); 
%%
figure; 
isosurface(DVScale, MLScale, APScale, smoLesions);
xlim([0 100]); ylim([0 100]); zlim([0 100]);
title('Scaled Lesions');
xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
saveas(gcf, [filepath filesep 'PT4_Figures' filesep 'PT4_ScaledLesionsFigure']);

figure; 
isosurface(DVScale, MLScale, APScale, smoBrain);
xlim([0 100]); ylim([0 100]); zlim([0 100]);
title('Scaled Brain');
xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
saveas(gcf, [filepath filesep 'PT4_Figures' filesep 'PT4_ScaledBrainFigure']);

% Interpolate between points
grid = {MLScale, DVScale, APScale};
LGI = griddedInterpolant(grid, smoLesions); % Default linear
BGI = griddedInterpolant(grid, smoBrain);

[xx, yy, zz] = ndgrid(MLScale, DVScale, APScale);

interpLesions = LGI(xx, yy, zz);
interpBrain = BGI(xx, yy, zz);

figure;
isosurface(DVScale, MLScale, APScale, interpLesions);
xlim([0 100]); ylim([0 100]); zlim([0 100]);
title('Scaled Interpolated Lesions');
xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
saveas(gcf, [filepath filesep 'PT4_Figures' filesep 'PT4_ScaledInterpLesionsFigure']);

figure;
isosurface(DVScale, MLScale, APScale, interpBrain);
xlim([0 100]); ylim([0 100]); zlim([0 100]);
title('Scaled Interpolated Brain');
xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
saveas(gcf, [filepath filesep 'PT4_Figures' filesep 'PT4_ScaledInterpBrainFigure']);

%% Overlay both figures and Capture a Rotating Video
figure;
lesionSurf= isosurface(DVScale, MLScale, APScale, interpLesions);
lesionPatch = patch(lesionSurf);
isonormals(DVScale,MLScale,APScale,interpLesions,lesionPatch);
set(lesionPatch,'FaceColor','red','EdgeColor','none','FaceAlpha',1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud

brainSurf = isosurface(DVScale, MLScale, APScale, interpBrain);
brainPatch = patch(brainSurf);
isonormals(DVScale,MLScale,APScale,interpBrain,brainPatch);
% set the color, mesh and transparency level of the surface
% set(brainPatch,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2); % Blue
set(brainPatch,'FaceColor',[255 231 214]./255,'EdgeColor','none','FaceAlpha',0.6);
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud

% xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
axis off
% title('Photothrombotic Lesions Reconstructed');
saveas(gcf, [filepath filesep 'PT4_Figures' filesep 'PT4_3DReconstructionOverlayTan']);
set(gcf,'color','w');
set(gca, 'fontsize', 18);
% Capture an MP4 movie file of the figure rotating
OptionZ.FrameRate=30; OptionZ.Duration=5; OptionZ.Periodic=true;
viewZ = [-50, 20;
    -37.5, 10;
    -25, 0;
    -12.5, -10;
    0, -20;
    12.5, -10;
    25, 0;
    37.5, 10
    50, 20;
    37.5, 30;
    25, 40;
    12.5, 50;
    0, 60;
    -12.5, 50;
    -25, 40;
    -37.5, 30;
    -50, 20]; % first column is pan, second is tilt
CaptureFigVid(viewZ, [filepath filesep 'PT4_Figures' filesep 'PT4_ReconstructTan'], OptionZ)
%% Overlay with Re-oriented Matrix
figure;
lesionSurf= isosurface(DVScale, APScale, MLScale, permute(flipud(interpLesions), [3 2 1]));
lesionPatch = patch(lesionSurf);
isonormals(DVScale,APScale,MLScale,permute(flipud(interpLesions), [3 2 1]),lesionPatch);
set(lesionPatch,'FaceColor','red','EdgeColor','none','FaceAlpha',1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud

brainSurf = isosurface(DVScale, APScale, MLScale, permute(flipud(interpBrain), [3 2 1]));
brainPatch = patch(brainSurf);
isonormals(DVScale,APScale,MLScale,permute(flipud(interpBrain), [3 2 1]),brainPatch);
% set the color, mesh and transparency level of the surface
% set(brainPatch,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2); % Blue
set(brainPatch,'FaceColor',[255 231 214]./255,'EdgeColor','none','FaceAlpha',0.6);
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud

% xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');
axis off
% title('Photothrombotic Lesions Reconstructed');
% saveas(gcf, 'PT33DReconstructionOverlayTanPerm');
set(gcf,'color','w');
set(gca, 'fontsize', 18);
% Capture an MP4 movie file of the figure rotating
OptionZ.FrameRate=30; OptionZ.Duration=20; OptionZ.Periodic=true;
viewZ = [-60, 30;
    -30, 30;
    0, 30;
    30, 30;
    60, 30;
    90, 30;
    120, 30;
    150, 30;
    180, 30;
    210, 30;
    240, 30;
    270, 30;
    300, 30;
    330, 30;
    360, 30;
    390, 30;
    420, 30]; % first column is pan, second is tilt
CaptureFigVid(viewZ, [filepath filesep 'PT4_Figures' filesep 'PT4_ReconstructTan'], OptionZ)

%% Calculate Volumes of the Lesions

% Need to get the isovalue that was used in isosurface
isovalLesions = isovalue(interpLesions); 
% isovalLesions = 0.2772;

% Get the combined volume of the lesions
totalVolume = sum((interpLesions >= isovalLesions),'all') * dV;


% Get the DV projection of all of the lesions
volume = sum((interpLesions >= isovalLesions), 1);
volVals = zeros(size(volume,2), size(volume,3));
for i = 1:size(volume,2)
    for j = 1:size(volume,3)
        volVals(i,j) = volume(1,i,j);
    end
end
figure; contourf(volVals); colorbar;

%% Nomenclature for beam diameters (mm) [2 1; 2 1 0.5; 1 .05] -> 
% PT5 had lesion in only RH

% Need to isolate the lesion (manually) create rectangle with x and y coordinates
% XVector; YVector for Right lesion1
L8X = 7:18; L8Y = 226:301;

L8vol = sum(volVals(L8Y, L8X), 'all') * dV;

beamDiameters = [1.5];

volArray = [L8vol];

% figure; scatter([beamDiameters beamDiameters], volArray);
save([filepath filesep 'PT4_Volumes'], 'beamDiameters', 'volArray');
