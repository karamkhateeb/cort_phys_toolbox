% Register OCTA images with cranial window images

filepath = 'C:\Users\kkhat\Documents\Everything\Documents\Yazdan Lab\PT Stroke Model\Figure Making\ManuscriptFigures\Figure2\Co-registration';

% [optimizer, metric] = imregconfig('multimodal');

for iPT = 2:3
    if iPT == 2 % PT2
        brainL = imread([filepath filesep 'PT2' filesep 'PT2_left.png']);
        brainR = imread([filepath filesep 'PT2' filesep 'PT2_right.png']);
        
        baseOCTAL = imread([filepath filesep 'PT2' filesep 'baseOCTA_left.tif']);
        baseOCTAL = flip(baseOCTAL,2);
        baseOCTAR = imread([filepath filesep 'PT2' filesep 'baseOCTA_right.tif']);
        baseOCTAR = flip(baseOCTAR,2);
        
        postOCTAL = imread([filepath filesep 'PT2' filesep 'postOCTA_left.tif']);
        postOCTAL = flip(postOCTAL,2);
        postOCTAR = imread([filepath filesep 'PT2' filesep 'postOCTA_right.tif']);
        postOCTAR = flip(postOCTAR,2);
    else
        brainL = imread([filepath filesep 'PT3' filesep 'PT3_left.png']);
        brainR = imread([filepath filesep 'PT3' filesep 'PT3_right.png']);
        
        baseOCTAL = imread([filepath filesep 'PT3' filesep 'baseOCTA_left.tif']);
        baseOCTAL = flip(baseOCTAL,2);
        baseOCTAR = imread([filepath filesep 'PT3' filesep 'baseOCTA_right.tif']);
        baseOCTAR = flip(baseOCTAR,2);
        
        postOCTAL = imread([filepath filesep 'PT3' filesep 'postOCTA_left.tif']);
        postOCTAL = flip(postOCTAL,2);
        postOCTAR = imread([filepath filesep 'PT3' filesep 'postOCTA_right.tif']);
        postOCTAR = flip(postOCTAR,2);
    end
    
    figure(1), imshow(brainL);
    figure(2), imshow(brainR);
    
    [mp, fp] = cpselect(baseOCTAL(:,:,1:3), brainL, 'Wait', true);
    t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
    Rfixed = imref2d(size(brainL));
    regBaseOCTAL = imwarp(baseOCTAL(:,:,1:3), t, 'OutputView', Rfixed);
    figure(1), imshowpair(brainL, regBaseOCTAL, 'blend');
    
    [mp, fp] = cpselect(postOCTAL(:,:,1:3), brainL, 'Wait', true);
    t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
    Rfixed = imref2d(size(brainL));
    regPostOCTAL = imwarp(postOCTAL(:,:,1:3), t, 'OutputView', Rfixed);
    figure(1), imshowpair(brainL, regPostOCTAL, 'blend');
    
    [mp, fp] = cpselect(baseOCTAR(:,:,1:3), brainR, 'Wait', true);
    t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
    Rfixed = imref2d(size(brainR));
    regBaseOCTAR = imwarp(baseOCTAR(:,:,1:3), t, 'OutputView', Rfixed);
    figure(2), imshowpair(brainR, regBaseOCTAR, 'blend');
    
    [mp, fp] = cpselect(postOCTAR(:,:,1:3), brainR, 'Wait', true);
    t = fitgeotrans(mp,fp,'nonreflectivesimilarity');
    Rfixed = imref2d(size(brainR));
    regPostOCTAR = imwarp(postOCTAR(:,:,1:3), t, 'OutputView', Rfixed);
    figure(2), imshowpair(brainR, regPostOCTAR, 'blend');
    
    switch iPT
        case 2
            baseNameL = 'PT2_baseOCTAL.tif';
            postNameL = 'PT2_postOCTAL.tif';
            
            baseNameR = 'PT2_baseOCTAR.tif';
            postNameR = 'PT2_postOCTAR.tif';
        case 3
            baseNameL = 'PT3_baseOCTAL.tif';
            postNameL = 'PT3_postOCTAL.tif';
            
            baseNameR = 'PT3_baseOCTAR.tif';
            postNameR = 'PT3_postOCTAR.tif';
    end
    
    imwrite(regBaseOCTAL,baseNameL);
    imwrite(regPostOCTAL,postNameL);
    imwrite(regBaseOCTAR,baseNameR);
    imwrite(regPostOCTAR,postNameR);
end

