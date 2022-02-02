% OCTA Analysis
% By: Karam Khateeb
% June 9, 2020
% October 5, 2020 - added depths and widths analysis
%
% Goal: Analyze and compare the OCTA-detected lesion sizes, lesion volumes
% as measured from histology, and the light intensities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Important Information
% Organization of Data:
% row 1 is monkey B left hemi - single FO (higher light intensity)
% row 2 is monkey B right hemi - triple FO (lower light intensity)
% row 3 is monkey C left hemi - single FO
% row 4 is monkey C right hemi - triple FO
% row 5 is monkey D left hemi - single FO no lens (even higher intensity?)
% row 6 is monkey E right hemi - single FO no lens

% here is a psuedo map of the lesions with their diameters in mm:
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

% Monkey D and E illumination diameter was 1.5, this will be given an index of 8

% make a key of matrix positions for each lesion detected with histology (to
% translate from Julien's loop scripts that get histological width and
% depth
histKey = ...
    [1 1; % PT2 L1
    1 2; % PT2 L2
    1 3; % PT2 L3
    1 4; % PT2 L4
    1 5; % PT2 L5
    1 6; % PT2 L6
    2 2; % PT2 R2
    2 4; % PT2 R4
    2 5; % PT2 R5
    2 6; % PT2 R6
    2 7; % PT2 R7
    3 6; % PT3 L6
    3 7; % PT3 L7
    3 3; % PT3 L3
    3 4; % PT3 L4
    3 1; % PT3 L1
    3 2; % PT3 L2
    4 1; % PT3 R1
    4 2; % PT3 R2
    4 4; % PT3 R4
    5 8; % PT4 lesion
    6 8; % PT5 lesion
    1 7; % PT2 L7
    2 1; % PT2 R1
    2 3; % PT2 R3
    3 5; % PT3 L5
    4 3; % PT3 R3
    4 5; % PT3 R5
    4 6; % PT3 R6
    4 7]; % PT R7
     
% make a key of matrix positions for each lesion undetected with histology (to
% translate from Julien's loop scripts that get histological width and
% depth
zeroHKey = ...
    [1 7;
    2 1;
    2 3;
    3 5;
    4 3;
    4 5;
    4 6;
    4 7];

% make a key of matrix positions for simulation data
simKey = ...
    [1 1;
    1 3;
    1 4;
    1 5;
    1 6;
    2 2;
    2 4;
    2 5;
    2 6;
    2 7;
    3 3;
    3 4;
    3 1;
    3 2;
    4 1;
    4 2;
    4 4;
    5 8;
    6 8;
    1 7;
    4 6;
    4 7];

% Make a key to denote which lesions are over a large blood vessel
vessel = zeros(6,8);
vesKey = ...
    [1 2;
    2 1;
    2 3;
    3 5;
    3 6;
    3 7;
    4 3;
    4 5];
for i = 1:size(vesKey,1)
    vessel(vesKey(i,1),vesKey(i,2)) = 1;
end

% Illumination diameters (note the diameters for monkey C left hemi are 
% going to be corrected always throughout this script such that they align 
% with every other hemisphere
apDs = [2 1 2 1 .5 1 .5 1.5]; apDs = repmat(apDs,6,1);

% load OCTA lesion information
load('OCTALesionInfo.mat');

% place OCTA-measured diameters in a matrix
octaDs = NaN(6,8);
for i = 1:numel(lesions)
    octaDs(i,lesions{i}.lesIDs) = lesions{i}.mmRadii*2;
    octaDs(i,lesions{i}.noIDs) = 0;
end

octaMaxDs = NaN(6,8);
for i = 1:numel(lesions)
    octaMaxDs(i,lesions{i}.lesIDs) = lesions{i}.mmMajD;
    octaMaxDs(i,lesions{i}.noIDs) = 0;
end

% OCTA Ds without zero values (cleaned)
octaMaxCleanDs = NaN(6,8);
for i = 1:numel(lesions)
    octaMaxCleanDs(i,lesions{i}.lesIDs) = lesions{i}.mmMajD;
    octaMaxCleanDs(i,lesions{i}.noIDs) = NaN;
end

octaMinDs = NaN(6,8);
for i = 1:numel(lesions)
    octaMinDs(i,lesions{i}.lesIDs) = lesions{i}.mmMinD;
    octaMinDs(i,lesions{i}.noIDs) = 0;
end

% OCTA Ds without large vessel values
octaMaxVDs = octaMaxDs;
octaMaxCleanVDs = octaMaxCleanDs;
for i = 1:size(vesKey,1)
    octaMaxVDs(vesKey(i,1),vesKey(i,2)) = NaN;
    octaMaxCleanVDs(vesKey(i,1),vesKey(i,2)) = NaN;
end

% Assess differences in OCTA area and aperture area
octaAreas = pi * octaMinDs .* octaMaxDs * .25;
apAreas = pi * (apDs.^2) * .25;
normAreas = octaAreas ./ apAreas;

% load lesion volumes from the histology 3D reconstruction
volumes = NaN(6,8); % volumes = zeros(6,8);
load('PT2_Volumes.mat','volArray'); 
volumes(1,1:7) = volArray(1:7); volumes(1,6) = 2.2189;
volumes(2,1:7) = volArray(8:14);

load('PT3_Volumes.mat','volArray'); 
volumes(3,1:7) = volArray(1:7);
volumes(4,1:7) = volArray(8:14);

load('PT4_Volumes.mat','volArray'); 
volumes(5,8) = volArray;

load('PT5_Volumes.mat','volArray');
volumes(6,8) = volArray;

volNaNs = volumes; volNaNs(find(volNaNs == 0)) = NaN;

% add depths - note: orientation of pt3_L lesions is flipped from JB
% orientation
% max depths
% depths = ...
%     [303 155 190 166 199 83 0 NaN;
%     0 126 0 228 171 157 200 NaN;
%     242 130 342 248 0 222 219 NaN;
%     227 136 0 279 0 0 0 NaN;
%     NaN NaN NaN NaN NaN NaN NaN 262;
%     NaN NaN NaN NaN NaN NaN NaN 233] * .01; % depths measured in mm (of representative Nissl slices)

% load metadata for lesions that includes complete histo lesion and (outdated) simulation data
hist_metadata = load(['data' filesep 'model_slices_w_nolesions.mat']);
hist_metadata = hist_metadata.full_metadata;
hist_metadata = hist_metadata([1:22 24:31]); % theres an extra PT2 L6 that shouldn't be marked as zero

% load metadata to be used for simulation analysis only (excludes zeros and
% vessels, but also includes corresponding histo data)
sim_metadata = load(['data' filesep 'final_full_metadata_nobigvessels.mat'],'scaled_metadata');
sim_metadata = sim_metadata.scaled_metadata;

d_w_a = zeros(numel(hist_metadata), 4);
for lesion_i = 1:numel(hist_metadata)
    lmd = hist_metadata(lesion_i);
    lesion = lmd.slice;
    if lesion == 0
        [depth, width, area, midwidth] = deal(0, 0, 0, 0);
    else
        [depth, width, area, midwidth] = get_lesion_maxdepths_and_avgwidths(lesion);
    end
    d_w_a(lesion_i, :) = [depth, width, area, midwidth];
end

depths = NaN(6,8);
for i = 1:size(d_w_a,1)
    depths(histKey(i,1),histKey(i,2)) = d_w_a(i,1) * .01; % in mm
end

depNaNs = depths;
for i = 1:size(zeroHKey,1)
    depNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

% Histo Ds without large vessel values
depVs = depths;
depNaNVs = depNaNs;
for i = 1:size(vesKey,1)
    depVs(vesKey(i,1),vesKey(i,2)) = NaN;
    depNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

% average depths
% depths = ...
%     [232.870198064079 118.328530092593 157.30553970464 146.2512211474 169.060426256552 72.1359347699476 0 NaN;
%     0 100.824811117107 0 192.628677055879 154.745786776962 130.66100433974 162.545536890646 NaN;
%     188.848801832291 102.930099408582 282.582883153261 194.413453284346 0 142.276459717297 177.736357475555 NaN;
%     182.345381452044 114.276017731506 0 216.528083159317 0 0 0 NaN;
%     NaN NaN NaN NaN NaN NaN NaN 220.754901761586;
%     NaN NaN NaN NaN NaN NaN NaN 201.924176660104] * .01; 
% 
% depNaNs = depths; depNaNs(find(depNaNs == 0)) = NaN;

widths = NaN(6,8);
for i = 1:size(d_w_a,1)
    widths(histKey(i,1),histKey(i,2)) = d_w_a(i,4) * .01; % in mm
end

widNaNs = widths;
for i = 1:size(zeroHKey,1)
    widNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

% Histo Ds without large vessel values
widVs = widths;
widNaNVs = widNaNs;
for i = 1:size(vesKey,1)
    widVs(vesKey(i,1),vesKey(i,2)) = NaN;
    widNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

% widths = ...
%     [227.927392739274 110.877419354839 236.573684210526 251.120481927711 188.64824120603 103.204819277108 0 NaN;
%     0 116.793650793651 0 247.535087719298 93.9824561403509 50.3694267515924 120.47 NaN;
%     241.095041322314 121.738461538462 233.067251461988 235.189516129032 0 180.630630630631 444.054794520548 NaN;
%     238.229074889868 180.014705882353 0 232.874551971326 0 0 0 NaN;
%     NaN NaN NaN NaN NaN NaN NaN 375.190839694656;
%     NaN NaN NaN NaN NaN NaN NaN 386.605150214592] * .01; % widths measured in mm (of representative Nissl slices)
% 
% widNaNs = widths; widNaNs(find(widNaNs == 0)) = NaN;

% get scaled simulated lesion dimensions
c_d_w_a = zeros(numel(sim_metadata), 4);
for lesion_i = 1:numel(sim_metadata)
    lmd = sim_metadata(lesion_i);
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 7.9e-05);
    cntr = lmd.scaled_simslice;
    [depth, width, area, midwidth] = get_lesion_maxdepths_and_avgwidths(cntr);
%     [depth, width, area, midwidth] = get_d_w_a_contour_slice(cntr);
    c_d_w_a(lesion_i, :) = [depth, width, area, midwidth];
end

simDeps = NaN(6,8);
for i = 1:size(c_d_w_a,1)
    simDeps(simKey(i,1),simKey(i,2)) = c_d_w_a(i,1) * .01; % in mm
end

simDepNaNs = simDeps;
for i = 1:size(zeroHKey,1)
    simDepNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

simDepVs = simDeps;
simDepNaNVs = simDepNaNs;
for i = 1:size(vesKey,1)
    simDepVs(vesKey(i,1),vesKey(i,2)) = NaN;
    simDepNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

simWids = NaN(6,8);
for i = 1:size(c_d_w_a,1)
    simWids(simKey(i,1),simKey(i,2)) = c_d_w_a(i,2) * .01; % in mm
end

simWidNaNs = simWids;
for i = 1:size(zeroHKey,1)
    simWidNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

simWidVs = simWids;
simWidNaNVs = simWidNaNs;
for i = 1:size(vesKey,1)
    simWidVs(vesKey(i,1),vesKey(i,2)) = NaN;
    simWidNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

% get NON-scaled simulated lesion dimensions
nsc_d_w_a = zeros(numel(sim_metadata), 4);
for lesion_i = 1:numel(sim_metadata)
    lmd = sim_metadata(lesion_i);
%     cntr = get_model_contour(lmd.aperture, lmd.intensity, 7.9e-05);
    cntr = lmd.simslice;
    [depth, width, area, midwidth] = get_lesion_maxdepths_and_avgwidths(cntr);
%     [depth, width, area, midwidth] = get_d_w_a_contour_slice(cntr);
    nsc_d_w_a(lesion_i, :) = [depth, width, area, midwidth];
end

nsSimDeps = NaN(6,8);
for i = 1:size(nsc_d_w_a,1)
    nsSimDeps(simKey(i,1),simKey(i,2)) = nsc_d_w_a(i,1) * .01; % in mm
end

nsSimDepNaNs = nsSimDeps;
for i = 1:size(zeroHKey,1)
    nsSimDepNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

nsSimDepVs = nsSimDeps;
nsSimDepNaNVs = nsSimDepNaNs;
for i = 1:size(vesKey,1)
    nsSimDepVs(vesKey(i,1),vesKey(i,2)) = NaN;
    nsSimDepNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

nsSimWids = NaN(6,8);
for i = 1:size(nsc_d_w_a,1)
    nsSimWids(simKey(i,1),simKey(i,2)) = nsc_d_w_a(i,2) * .01; % in mm
end

simWidNaNs = simWids;
for i = 1:size(zeroHKey,1)
    nsSimWidNaNs(zeroHKey(i,1),zeroHKey(i,2)) = NaN;
end

nsSimWidVs = nsSimWids;
nsSimWidNaNVs = nsSimWidNaNs;
for i = 1:size(vesKey,1)
    nsSimWidVs(vesKey(i,1),vesKey(i,2)) = NaN;
    nsSimWidNaNVs(vesKey(i,1),vesKey(i,2)) = NaN;
end

% % Simulation average depths and average widths
% simDeps = ...
%     [143.437486771687 95.9158916432751 140.785808740024 170.094874565163 79.877546525713 85.3326430235568 NaN NaN;
%     NaN 73.8088004190676 NaN 164.297755545317 70.0097060261155 58.583319032092 67.5943052561148 NaN;
%     143.437486771687 95.9158916432751 140.785808740024 170.094874565163 NaN 85.3326430235568 83.5975745325922 NaN;
%     129.094295463884 73.8088004190676 NaN 164.297755545317 NaN NaN NaN NaN;
%     NaN NaN NaN NaN NaN NaN NaN 241.099932987572;
%     NaN NaN NaN NaN NaN NaN NaN 241.099932987572] * .01;
% 
% simWids = ...
%     [374.829545454545 245.575 370.630057803468 395.591549295775 186.237113402062 219.742857142857 NaN NaN;
%     NaN 189.9 NaN 373.458128078818 167.482352941176 163.140845070423 163.036585365854 NaN;
%     374.829545454545 245.575 370.630057803468 395.591549295775 NaN 219.742857142857 194.940594059406 NaN;
%     348.83125 189.9 NaN 373.458128078818 NaN NaN NaN NaN;
%     NaN NaN NaN NaN NaN NaN NaN 606.655290102389;
%     NaN NaN NaN NaN NaN NaN NaN 606.655290102389] * .01;

% Load light intensities
% Nomenclature for beam diameters (mm) [2 1; 2 1 0.5; 1 .05; 1.5] -> 
% R/L[1 2; 3 4 5; 6 7; 8]
% L0pwrArr = [0 0 0]; R0pwrArr = [0 0 0]; % light intensity measurements with all apertures covered
% L1pwrArr = [2.06 1.27 2.13]; R1pwrArr = [1.13 1.21 1.65];
% L2pwrArr = [0.087 0.087 0.087]; R2pwrArr = [0.026 0.027 0.026];
% L3pwrArr = [1.83 1.87 1.87]; R3pwrArr = [1.39 1.45 1.47];
% L4pwrArr = [2.64 2.63 2.48]; R4pwrArr = [2.24 2.23 2.24];
% L5pwrArr = [0.116 0.141 0.171]; R5pwrArr = [0.063 0.063 0.62];
% L6pwrArr = [0.30 0.31 0.32]; R6pwrArr = [0.16 0.17 0.15];
% L7pwrArr = [0.216 0.212 0.217]; R7pwrArr = [0.095 0.104 0.103];
% L8pwrArr = [20.1 21.1 20.6]; 

load('LightIntensities.mat');

Lpwrs = [mean(L1pwrArr) mean([L2pwrArr L6pwrArr]) mean(L3pwrArr) ...
    mean(L4pwrArr) mean(L5pwrArr) mean([L2pwrArr L6pwrArr]) ...
    mean(L7pwrArr) mean(L8pwrArr)];
Rpwrs = [mean(R1pwrArr) mean([R2pwrArr R6pwrArr]) mean(R3pwrArr) ...
    mean(R4pwrArr) mean(R5pwrArr) mean([R2pwrArr R6pwrArr]) ...
    mean(R7pwrArr) mean(R8pwrArr)];

pwrs = [Lpwrs;Rpwrs;Lpwrs;Rpwrs;Lpwrs;Rpwrs];

% make everything as a single column for plotting
colVes =reshape(vessel,numel(vessel),1);
colApDs = reshape(apDs,numel(apDs),1);
colOCTADs = reshape(octaMaxDs,numel(octaMaxDs),1);
colOCTADsC = reshape(octaMaxCleanDs,numel(octaMaxCleanDs),1);
colOCTADsV = reshape(octaMaxVDs,numel(octaMaxVDs),1);
colOCTADsVC = reshape(octaMaxCleanVDs,numel(octaMaxCleanVDs),1);
colVols = reshape(volumes,numel(volumes),1);
colVolsC = reshape(volNaNs,numel(volNaNs),1);
colDeps = reshape(depths,numel(depths),1);
colDepsC = reshape(depNaNs,numel(depNaNs),1);
colDepsV = reshape(depVs,numel(depVs),1);
colDepsVC = reshape(depNaNVs,numel(depNaNVs),1);
colWids = reshape(widths,numel(widths),1);
colWidsC = reshape(widNaNs,numel(widNaNs),1);
colWidsV = reshape(widVs,numel(widVs),1);
colWidsVC = reshape(widNaNVs,numel(widNaNVs),1);
colPwrs = reshape(pwrs,numel(pwrs),1);
colDens = colPwrs ./ (pi * (colApDs/2).^2);
colSimDeps = reshape(simDeps,numel(simDeps),1);
colSimDepsC = reshape(simDepNaNs,numel(simDepNaNs),1);
colSimDepsV = reshape(simDepVs,numel(simDepVs),1);
colSimDepsVC = reshape(simDepNaNVs,numel(simDepNaNVs),1);
colSimWids = reshape(simWids,numel(simWids),1);
colSimWidsC = reshape(simWidNaNs,numel(simWidNaNs),1);
colSimWidsV = reshape(simWidVs,numel(simWidVs),1);
colSimWidsVC = reshape(simWidNaNVs,numel(simWidNaNVs),1);
colNSSimDeps = reshape(nsSimDeps,numel(nsSimDeps),1);
colNSSimDepsC = reshape(nsSimDepNaNs,numel(nsSimDepNaNs),1);
colNSSimDepsV = reshape(nsSimDepVs,numel(nsSimDepVs),1);
colNSSimDepsVC = reshape(nsSimDepNaNVs,numel(nsSimDepNaNVs),1);
colNSSimWids = reshape(nsSimWids,numel(nsSimWids),1);
colNSSimWidsC = reshape(nsSimWidNaNs,numel(nsSimWidNaNs),1);
colNSSimWidsV = reshape(nsSimWidVs,numel(nsSimWidVs),1);
colNSSimWidsVC = reshape(nsSimWidNaNVs,numel(nsSimWidNaNVs),1);


contApDs = linspace(min(colApDs),max(colApDs),100);
contPwrs = linspace(min(colPwrs),max(colPwrs),100);
contDens = linspace(min(colDens),max(colDens),100);
contOCTADs = linspace(min(colOCTADs),max(colOCTADs),100);
contOCTADsC = linspace(min(colOCTADsC),max(colOCTADsC),100);
contOCTADsV = linspace(min(colOCTADsV),max(colOCTADsV),100);
contOCTADsVC = linspace(min(colOCTADsVC),max(colOCTADsVC),100);
contWids = linspace(min(colWids),max(colWids),100);
contWidsC = linspace(min(colWidsC),max(colWidsC),100);
contWidsV = linspace(min(colWidsV),max(colWidsV),100);
contWidsVC = linspace(min(colWidsVC),max(colWidsVC),100);
contSimDeps = linspace(min(colSimDeps),max(colSimDeps),100);
contSimDepsC = linspace(min(colSimDepsC),max(colSimDepsC),100);
contSimDepsV = linspace(0,max(colSimDepsV),100);
contSimDepsVC = linspace(min(colSimDepsVC),max(colSimDepsVC),100);
contSimWids = linspace(min(colSimWids),max(colSimWids),100);
contSimWidsC = linspace(min(colSimWidsC),max(colSimWidsC),100);
contSimWidsV = linspace(0,max(colSimWidsV),100);
contSimWidsVC = linspace(min(colSimWidsVC),max(colSimWidsVC),100);
contNSSimDeps = linspace(min(colNSSimDeps),max(colNSSimDeps),100);
contNSSimDepsC = linspace(min(colNSSimDepsC),max(colNSSimDepsC),100);
contNSSimDepsV = linspace(0,max(colNSSimDepsV),100);
contNSSimDepsVC = linspace(min(colNSSimDepsVC),max(colNSSimDepsVC),100);
contNSSimWids = linspace(min(colNSSimWids),max(colNSSimWids),100);
contNSSimWidsC = linspace(min(colNSSimWidsC),max(colNSSimWidsC),100);
contNSSimWidsV = linspace(min(colNSSimWidsV),max(colNSSimWidsV),100);
contNSSimWidsVC = linspace(min(colNSSimWidsVC),max(colNSSimWidsVC),100);

wfun = 'ols';
cepts = true; % have intercepts? true or false
if ~cepts
    numCoFs = 1;
else
    numCoFs = 2;
end
%% 2. Model OCTA Lesion Diameters and Aperture Diameter

mdl_apDs_OCTADs = fitlm(colApDs,colOCTADs,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_OCTADs = mdl_apDs_OCTADs.Coefficients.Estimate; % store coefficients
if cepts
    reg_apDs_OCTADs = coF_apDs_OCTADs(1) + (coF_apDs_OCTADs(2) * contApDs); % regression line for the model
    est_apDs_OCTADs = coF_apDs_OCTADs(1) + (coF_apDs_OCTADs(2) * colApDs);
else
    reg_apDs_OCTADs = (coF_apDs_OCTADs * contApDs); 
    est_apDs_OCTADs = (coF_apDs_OCTADs * colApDs);
end
r2_apDs_OCTADs = computeR2(colOCTADs, est_apDs_OCTADs);

mdl_apDs_OCTADsC = fitlm(colApDs,colOCTADsC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_OCTADsC = mdl_apDs_OCTADsC.Coefficients.Estimate;
if cepts
    reg_apDs_OCTADsC = coF_apDs_OCTADsC(1) + (coF_apDs_OCTADsC(2) * contApDs);
    est_apDs_OCTADsC = coF_apDs_OCTADsC(1) + (coF_apDs_OCTADsC(2) * colApDs);
else
    reg_apDs_OCTADsC = (coF_apDs_OCTADsC * contApDs);
    est_apDs_OCTADsC = (coF_apDs_OCTADsC * colApDs);
end
r2_apDs_OCTADsC = computeR2(colOCTADsC, est_apDs_OCTADsC);

mdl_apDs_OCTADsV = fitlm(colApDs,colOCTADsV,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_OCTADsV = mdl_apDs_OCTADsV.Coefficients.Estimate;
if cepts
    reg_apDs_OCTADsV = coF_apDs_OCTADsV(1) + (coF_apDs_OCTADsV(2) * contApDs);
    est_apDs_OCTADsV = coF_apDs_OCTADsV(1) + (coF_apDs_OCTADsV(2) * colApDs);
else
    reg_apDs_OCTADsV = (coF_apDs_OCTADsV * contApDs);
    est_apDs_OCTADsV = (coF_apDs_OCTADsV * colApDs);
end
r2_apDs_OCTADsV = computeR2(colOCTADsV, est_apDs_OCTADsV);

mdl_apDs_OCTADsVC = fitlm(colApDs,colOCTADsVC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_OCTADsVC = mdl_apDs_OCTADsVC.Coefficients.Estimate;
if cepts
    reg_apDs_OCTADsVC = coF_apDs_OCTADsVC(1) + (coF_apDs_OCTADsVC(2) * contApDs);
    est_apDs_OCTADsVC = coF_apDs_OCTADsVC(1) + (coF_apDs_OCTADsVC(2) * colApDs);
else
    reg_apDs_OCTADsVC = (coF_apDs_OCTADsVC * contApDs);
    est_apDs_OCTADsVC = (coF_apDs_OCTADsVC * colApDs);
end
r2_apDs_OCTADsVC = computeR2(colOCTADsVC, est_apDs_OCTADsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_OCTADs)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Max. OCTA Diameter vs. Aperture Diameter')
% box off
% anova(mdl_apDs_OCTADsC,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_apDs_OCTADsC)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Max. OCTA Diameter (Detected) vs. Aperture Diameter')
% box off
% anova(mdl_apDs_OCTADsC,'Summary');

%% 3. Model Lesion Volumes and Aperture Diameter

mdl_apDs_vols = fitlm(colApDs,colVols,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_vols = mdl_apDs_vols.Coefficients.Estimate;
if cepts
    reg_apDs_vols = coF_apDs_vols(1) + (coF_apDs_vols(2) * contApDs);
    est_apDs_vols = coF_apDs_vols(1) + (coF_apDs_vols(2) * colApDs);
else
    reg_apDs_vols = (coF_apDs_vols * contApDs);
    est_apDs_vols = (coF_apDs_vols * colApDs);
end
r2_apDs_vols = computeR2(colVols, est_apDs_vols);

mdl_apDs_volsC = fitlm(colApDs,colVolsC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_volsC = mdl_apDs_volsC.Coefficients.Estimate;
if cepts
    reg_apDs_volsC = coF_apDs_volsC(1) + (coF_apDs_volsC(2) * contApDs);
    est_apDs_volsC = coF_apDs_volsC(1) + (coF_apDs_volsC(2) * colApDs);
else
    reg_apDs_volsC = (coF_apDs_volsC * contApDs);
    est_apDs_volsC = (coF_apDs_volsC * colApDs);
end
r2_apDs_volsC = computeR2(colVolsC, est_apDs_volsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_vols)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Volume vs. Aperture Diameter')
% box off
% anova(mdl_apDs_vols,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_apDs_volsC)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Volume (Detected) vs. Aperture Diameter')
% box off
% anova(mdl_apDs_vols,'Summary');

%% 4. Model OCTA Diameters and Light Intensity

mdl_pwrs_OCTADs = fitlm(colPwrs,colOCTADs,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_OCTADs = mdl_pwrs_OCTADs.Coefficients.Estimate;
if cepts
    reg_pwrs_OCTADs = coF_pwrs_OCTADs(1) + (coF_pwrs_OCTADs(2) * contPwrs);
    est_pwrs_OCTADs = coF_pwrs_OCTADs(1) + (coF_pwrs_OCTADs(2) * colPwrs);
else
    reg_pwrs_OCTADs = (coF_pwrs_OCTADs * contPwrs);
    est_pwrs_OCTADs = (coF_pwrs_OCTADs * colPwrs);
end
r2_pwrs_OCTADs = computeR2(colOCTADs, est_pwrs_OCTADs);

mdl_pwrs_OCTADsC = fitlm(colPwrs,colOCTADsC,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_OCTADsC = mdl_pwrs_OCTADsC.Coefficients.Estimate;
if cepts
    reg_pwrs_OCTADsC = coF_pwrs_OCTADsC(1) + (coF_pwrs_OCTADsC(2) * contPwrs);
    est_pwrs_OCTADsC = coF_pwrs_OCTADsC(1) + (coF_pwrs_OCTADsC(2) * colPwrs);
else
    reg_pwrs_OCTADsC = (coF_pwrs_OCTADsC * contPwrs);
    est_pwrs_OCTADsC = (coF_pwrs_OCTADsC * colPwrs);
end
r2_pwrs_OCTADsC = computeR2(colOCTADsC, est_pwrs_OCTADsC);

mdl_pwrs_OCTADsV = fitlm(colPwrs,colOCTADsV,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_OCTADsV = mdl_pwrs_OCTADsV.Coefficients.Estimate;
if cepts
    reg_pwrs_OCTADsV = coF_pwrs_OCTADsV(1) + (coF_pwrs_OCTADsV(2) * contPwrs);
    est_pwrs_OCTADsV = coF_pwrs_OCTADsV(1) + (coF_pwrs_OCTADsV(2) * colPwrs);
else
    reg_pwrs_OCTADsV = (coF_pwrs_OCTADsV * contPwrs);
    est_pwrs_OCTADsV = (coF_pwrs_OCTADsV * colPwrs);
end
r2_pwrs_OCTADsV = computeR2(colOCTADsV, est_pwrs_OCTADsV);

mdl_pwrs_OCTADsVC = fitlm(colPwrs,colOCTADsVC,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_OCTADsVC = mdl_pwrs_OCTADsVC.Coefficients.Estimate;
if cepts
    reg_pwrs_OCTADsVC = coF_pwrs_OCTADsVC(1) + (coF_pwrs_OCTADsVC(2) * contPwrs);
    est_pwrs_OCTADsVC = coF_pwrs_OCTADsVC(1) + (coF_pwrs_OCTADsVC(2) * colPwrs);
else
    reg_pwrs_OCTADsVC = (coF_pwrs_OCTADsVC * contPwrs);
    est_pwrs_OCTADsVC = (coF_pwrs_OCTADsVC * colPwrs);
end
r2_pwrs_OCTADsV = computeR2(colOCTADsVC, est_pwrs_OCTADsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_pwrs_OCTADs)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Light Intensity (mW)')
% title('Max. OCTA Diameter vs. Light Intensity')
% box off
% anova(mdl_pwrs_OCTADs,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_pwrs_OCTADsC)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Light Intensity (mW)')
% title('Max. OCTA Diameter (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_OCTADsC,'Summary');

%% 5. Model OCTA Diameters and ln(Light Intensity)

if cepts
    mdl_lnpwrs_OCTADs = fitlm(log(colPwrs),colOCTADs,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADs = mdl_lnpwrs_OCTADs.Coefficients.Estimate;
    reg_lnpwrs_OCTADs = coF_lnpwrs_OCTADs(1) + (coF_lnpwrs_OCTADs(2) * log(contPwrs));
    est_lnpwrs_OCTADs = coF_lnpwrs_OCTADs(1) + (coF_lnpwrs_OCTADs(2) * log(colPwrs));
else
    mdl_lnpwrs_OCTADs = fitlm(log(colPwrs+1),colOCTADs,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADs = mdl_lnpwrs_OCTADs.Coefficients.Estimate;
    reg_lnpwrs_OCTADs = (coF_lnpwrs_OCTADs * log(contPwrs+1));
    est_lnpwrs_OCTADs = (coF_lnpwrs_OCTADs * log(colPwrs+1));
end
r2_lnpwrs_OCTADs = computeR2(colOCTADs, est_lnpwrs_OCTADs);

if cepts
    mdl_lnpwrs_OCTADsC = fitlm(log(colPwrs),colOCTADsC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsC = mdl_lnpwrs_OCTADsC.Coefficients.Estimate;
    reg_lnpwrs_OCTADsC = coF_lnpwrs_OCTADsC(1) + (coF_lnpwrs_OCTADsC(2) * log(contPwrs));
    est_lnpwrs_OCTADsC = coF_lnpwrs_OCTADsC(1) + (coF_lnpwrs_OCTADsC(2) * log(colPwrs));
else
    mdl_lnpwrs_OCTADsC = fitlm(log(colPwrs+1),colOCTADsC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsC = mdl_lnpwrs_OCTADsC.Coefficients.Estimate;
    reg_lnpwrs_OCTADsC = (coF_lnpwrs_OCTADsC * log(contPwrs+1));
    est_lnpwrs_OCTADsC = (coF_lnpwrs_OCTADsC * log(colPwrs+1));
end
r2_lnpwrs_OCTADsC = computeR2(colOCTADsC, est_lnpwrs_OCTADsC);

if cepts
    mdl_lnpwrs_OCTADsV = fitlm(log(colPwrs),colOCTADsV,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsV = mdl_lnpwrs_OCTADsV.Coefficients.Estimate;
    reg_lnpwrs_OCTADsV = coF_lnpwrs_OCTADsV(1) + (coF_lnpwrs_OCTADsV(2) * log(contPwrs));
    est_lnpwrs_OCTADsV = coF_lnpwrs_OCTADsV(1) + (coF_lnpwrs_OCTADsV(2) * log(colPwrs));
else
    mdl_lnpwrs_OCTADsV = fitlm(log(colPwrs+1),colOCTADsV,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsV = mdl_lnpwrs_OCTADsV.Coefficients.Estimate;
    reg_lnpwrs_OCTADsV = (coF_lnpwrs_OCTADsV * log(contPwrs+1));
    est_lnpwrs_OCTADsV = (coF_lnpwrs_OCTADsV * log(colPwrs+1));
end
r2_lnpwrs_OCTADsV = computeR2(colOCTADsV, est_lnpwrs_OCTADsV);

if cepts
    mdl_lnpwrs_OCTADsVC = fitlm(log(colPwrs),colOCTADsVC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsVC = mdl_lnpwrs_OCTADsC.Coefficients.Estimate;
    reg_lnpwrs_OCTADsVC = coF_lnpwrs_OCTADsVC(1) + (coF_lnpwrs_OCTADsVC(2) * log(contPwrs));
    est_lnpwrs_OCTADsVC = coF_lnpwrs_OCTADsVC(1) + (coF_lnpwrs_OCTADsVC(2) * log(colPwrs));
else
    mdl_lnpwrs_OCTADsVC = fitlm(log(colPwrs+1),colOCTADsVC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_OCTADsVC = mdl_lnpwrs_OCTADsVC.Coefficients.Estimate;
    reg_lnpwrs_OCTADsVC = (coF_lnpwrs_OCTADsVC * log(contPwrs+1));
    est_lnpwrs_OCTADsVC = (coF_lnpwrs_OCTADsVC * log(colPwrs+1));
end
r2_lnpwrs_OCTADsVC = computeR2(colOCTADsVC, est_lnpwrs_OCTADsVC);

figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
subplot(1,2,1)
plot(mdl_lnpwrs_OCTADs)
ylabel('Max. OCTA Diameter (mm)')
xlabel('log(Light Intensity (mW))')
title('Max. OCTA Diameter vs. Light Intensity')
box off
anova(mdl_lnpwrs_OCTADs,'Summary');

subplot(1,2,2)
plot(mdl_lnpwrs_OCTADsC)
ylabel('Max. OCTA Diameter (mm)')
xlabel('log(Light Intensity (mW))')
title('Max. OCTA Diameter (Detected) vs. Light Intensity')
box off
anova(mdl_lnpwrs_OCTADsC,'Summary');

figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
subplot(1,2,1)
plot(mdl_lnpwrs_OCTADsV)
ylabel('Max. OCTA Diameter (mm)')
xlabel('log(Light Intensity (mW))')
title('Max. OCTA Diameter (Non-Vessel) vs. Light Intensity')
box off
anova(mdl_lnpwrs_OCTADsV,'Summary');

subplot(1,2,2)
plot(mdl_lnpwrs_OCTADsVC)
ylabel('Max. OCTA Diameter (mm)')
xlabel('log(Light Intensity (mW))')
title('Max. OCTA Diameter (Non-Vessel, Detected) vs. Light Intensity')
box off
anova(mdl_lnpwrs_OCTADsVC,'Summary');

%% 6. Model Lesion Volumes and Light Intensity

mdl_pwrs_vols = fitlm(colPwrs,colVols,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_vols = mdl_pwrs_vols.Coefficients.Estimate;
if cepts
    reg_pwrs_vols = coF_pwrs_vols(1) + (coF_pwrs_vols(2) * contPwrs);
    est_pwrs_vols = coF_pwrs_vols(1) + (coF_pwrs_vols(2) * colPwrs);
else
    reg_pwrs_vols = (coF_pwrs_vols * contPwrs);
    est_pwrs_vols = (coF_pwrs_vols * colPwrs);
end
r2_pwrs_vols = computeR2(colVols, est_pwrs_vols);

mdl_pwrs_volsC = fitlm(colPwrs,colVolsC,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_volsC = mdl_pwrs_volsC.Coefficients.Estimate;
if cepts
    reg_pwrs_volsC = coF_pwrs_volsC(1) + (coF_pwrs_volsC(2) * contPwrs);
    est_pwrs_volsC = coF_pwrs_volsC(1) + (coF_pwrs_volsC(2) * colPwrs);
else
    reg_pwrs_volsC = (coF_pwrs_volsC * contPwrs);
    est_pwrs_volsC = (coF_pwrs_volsC * colPwrs);
end
r2_pwrs_volsC = computeR2(colVolsC, est_pwrs_volsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_pwrs_vols)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Light Intensity (mW)')
% title('Lesion Volume vs. Light Intensity')
% box off
% anova(mdl_pwrs_vols,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_pwrs_volsC)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Light Intensity (mW)')
% title('Lesion Volume (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_volsC,'Summary');

%% 7. Model OCTA Diameters and Light Density

mdl_dens_OCTADs = fitlm(colDens,colOCTADs,'RobustOpts',wfun,'intercept',cepts)
coF_dens_OCTADs = mdl_dens_OCTADs.Coefficients.Estimate;
if cepts
    reg_dens_OCTADs = coF_dens_OCTADs(1) + (coF_dens_OCTADs(2) * contDens);
    est_dens_OCTADs = coF_dens_OCTADs(1) + (coF_dens_OCTADs(2) * colDens);
else
    reg_dens_OCTADs = (coF_dens_OCTADs * contDens);
    est_dens_OCTADs = (coF_dens_OCTADs * colDens);
end
r2_dens_OCTADs = computeR2(colOCTADs, est_dens_OCTADs);

mdl_dens_OCTADsC = fitlm(colDens,colOCTADsC,'RobustOpts',wfun,'intercept',cepts)
coF_dens_OCTADsC = mdl_dens_OCTADsC.Coefficients.Estimate;
if cepts
    reg_dens_OCTADsC = coF_dens_OCTADsC(1) + (coF_dens_OCTADsC(2) * contDens);
    est_dens_OCTADsC = coF_dens_OCTADsC(1) + (coF_dens_OCTADsC(2) * colDens);
else
    reg_dens_OCTADsC = (coF_dens_OCTADsC * contDens);
    est_dens_OCTADsC = (coF_dens_OCTADsC * colDens);
end
r2_dens_OCTADsC = computeR2(colOCTADsC, est_dens_OCTADsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_dens_OCTADs)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Light Density (mW/mm^2)')
% title('Max. OCTA Diameter vs. Light Density')
% box off
% anova(mdl_dens_OCTADs,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_dens_OCTADsC)
% ylabel('Max. OCTA Diameter (mm)')
% xlabel('Light Density (mW/mm^2)')
% title('Max. OCTA Diameter (Detected) vs. Light Density')
% box off
% anova(mdl_dens_OCTADsC,'Summary');

%% 8. Model Lesion Volumes and Light Density

mdl_dens_vols = fitlm(colDens,colVols,'RobustOpts',wfun,'intercept',cepts)
coF_dens_vols = mdl_dens_vols.Coefficients.Estimate;
if cepts
    reg_dens_vols = coF_dens_vols(1) + (coF_dens_vols(2) * contDens);
    est_dens_vols = coF_dens_vols(1) + (coF_dens_vols(2) * colDens);
else
    reg_dens_vols = (coF_dens_vols * contDens);
    est_dens_vols = (coF_dens_vols * colDens);
end
r2_dens_vols = computeR2(colVols, est_dens_vols);

mdl_dens_volsC = fitlm(colDens,colVolsC,'RobustOpts',wfun,'intercept',cepts)
coF_dens_volsC = mdl_dens_volsC.Coefficients.Estimate;
if cepts
    reg_dens_volsC = coF_dens_volsC(1) + (coF_dens_volsC(2) * contDens);
    est_dens_volsC = coF_dens_volsC(1) + (coF_dens_volsC(2) * colDens);
else
    reg_dens_volsC = (coF_dens_volsC * contDens);
    est_dens_volsC = (coF_dens_volsC * colDens);
end
r2_dens_volsC = computeR2(colVolsC, est_dens_volsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_dens_vols)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Light Density (mW/mm^2)')
% title('Lesion Volume vs. Light Density')
% box off
% anova(mdl_dens_vols,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_dens_volsC)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Light Density (mW/mm^2)')
% title('Lesion Volume (Detected) vs. Light Density')
% box off
% anova(mdl_dens_volsC,'Summary');

%% 9. See Correlation between OCTA Diameters and Lesion Volumes

mdl_OCTADs_vols = fitlm(colOCTADs,colVols,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_vols = mdl_OCTADs_vols.Coefficients.Estimate;
if cepts
    reg_OCTADs_vols = coF_OCTADs_vols(1) + (coF_OCTADs_vols(2) * contOCTADs);
    est_OCTADs_vols = coF_OCTADs_vols(1) + (coF_OCTADs_vols(2) * colOCTADs);
else
    reg_OCTADs_vols = (coF_OCTADs_vols * contOCTADs);
    est_OCTADs_vols = (coF_OCTADs_vols * colOCTADs);
end
r2_OCTADs_vols = computeR2(colVols, est_OCTADs_vols);

mdl_OCTADsC_vols = fitlm(colOCTADsC,colVols,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_vols = mdl_OCTADsC_vols.Coefficients.Estimate;
if cepts
    reg_OCTADsC_vols = coF_OCTADsC_vols(1) + (coF_OCTADsC_vols(2) * contOCTADsC);
    est_OCTADsC_vols = coF_OCTADsC_vols(1) + (coF_OCTADsC_vols(2) * colOCTADsC);
else
    reg_OCTADsC_vols = (coF_OCTADsC_vols * contOCTADsC);
    est_OCTADsC_vols = (coF_OCTADsC_vols * colOCTADsC);
end
r2_OCTADsC_vols = computeR2(colVols, est_OCTADsC_vols);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_OCTADs_vols)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Volume vs. Max. OCTA Diameter')
% box off
% anova(mdl_OCTADs_vols,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_OCTADsC_vols)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Volume vs. Max. OCTA Diameter (Detected)')
% box off
% anova(mdl_OCTADsC_vols,'Summary');

mdl_OCTADs_volsC = fitlm(colOCTADs,colVolsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_volsC = mdl_OCTADs_volsC.Coefficients.Estimate;
if cepts
    reg_OCTADs_volsC = coF_OCTADs_volsC(1) + (coF_OCTADs_volsC(2) * contOCTADs);
    est_OCTADs_volsC = coF_OCTADs_volsC(1) + (coF_OCTADs_volsC(2) * colOCTADs);
else
    reg_OCTADs_volsC = (coF_OCTADs_volsC * contOCTADs);
    est_OCTADs_volsC = (coF_OCTADs_volsC * colOCTADs);
end
r2_OCTADs_volsC = computeR2(colVolsC, est_OCTADs_volsC);

mdl_OCTADsC_volsC = fitlm(colOCTADsC,colVolsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_volsC = mdl_OCTADsC_volsC.Coefficients.Estimate;
if cepts
    reg_OCTADsC_volsC = coF_OCTADsC_volsC(1) + (coF_OCTADsC_volsC(2) * contOCTADsC);
    est_OCTADsC_volsC = coF_OCTADsC_volsC(1) + (coF_OCTADsC_volsC(2) * colOCTADsC);
else
    reg_OCTADsC_volsC = (coF_OCTADsC_volsC * contOCTADsC);
    est_OCTADsC_volsC = (coF_OCTADsC_volsC * colOCTADsC);
end
r2_OCTADsC_volsC = computeR2(colVolsC, est_OCTADsC_volsC);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(2,2,1)
% plot(mdl_OCTADs_vols), title('Lesion Volume vs Max. OCTA Diameter'), box off
% ylabel('Lesion Volume (mm^3)'), xlabel('')
% 
% subplot(2,2,2)
% plot(mdl_OCTADsC_vols), title('Lesion Volume vs Max. OCTA Diameter (Detected)'), box off
% xlabel(''), ylabel('');
% 
% subplot(2,2,3)
% plot(mdl_OCTADs_volsC), title('Lesion Volume (Detected) vs Max. OCTA Diameter'), box off
% ylabel('Lesion Volume (mm^3)'), xlabel('Max. OCTA Diameter (mm)')
% 
% subplot(2,2,4)
% plot(mdl_OCTADsC_volsC), title('Lesion Volume (Detected) vs Max. OCTA Diameter (Detected)'), box off
% xlabel('Max. OCTA Diameter (mm)'), ylabel('');

%% 10. Model Light Intensity and Aperture Diameters on OCTA Diameters 

tbl = [colPwrs colApDs];

mdl_apDs_pwrs_OCTADs = fitlm(tbl, colOCTADs, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (mm)'})
coF_apDs_pwrs_OCTADs = mdl_apDs_pwrs_OCTADs.Coefficients.Estimate;
reg_apDs_pwrs_OCTADs = coF_apDs_pwrs_OCTADs(1) + (coF_apDs_pwrs_OCTADs(2) * contPwrs) + (coF_apDs_pwrs_OCTADs(3) * contApDs);

mdl_apDs_pwrs_OCTADsC = fitlm(tbl, colOCTADsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Detected) (mm)'})
coF_apDs_pwrs_OCTADsC = mdl_apDs_pwrs_OCTADsC.Coefficients.Estimate;
reg_apDs_pwrs_OCTADsC = coF_apDs_pwrs_OCTADsC(1) + (coF_apDs_pwrs_OCTADsC(2) * contPwrs) + (coF_apDs_pwrs_OCTADsC(3) * contApDs);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_pwrs_OCTADs), box off
% title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter');
% 
% subplot(1,2,2)
% plot(mdl_apDs_pwrs_OCTADsC), box off
% title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter (Detected)');

[b,bint,r,rint,stats] = regress(reshape(octaMaxCleanDs,[],1),[tbl, ones(size(tbl,1),1)]);

%% 11. Model Light Intensity and Aperture Diameter on Lesion Volume

tbl = [colPwrs colApDs];

mdl_apDs_pwrs_vols = fitlm(tbl, colVols, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Lesion Volume (mm^3)'})
coF_apDs_pwrs_vols = mdl_apDs_pwrs_vols.Coefficients.Estimate;
reg_apDs_pwrs_vols = coF_apDs_pwrs_vols(1) + (coF_apDs_pwrs_vols(2) * contPwrs) + (coF_apDs_pwrs_vols(3) * contApDs);

mdl_apDs_pwrs_volsC = fitlm(tbl, colVolsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Detected) (mm)'})
coF_apDs_pwrs_volsC = mdl_apDs_pwrs_volsC.Coefficients.Estimate;
reg_apDs_pwrs_volsC = coF_apDs_pwrs_volsC(1) + (coF_apDs_pwrs_volsC(2) * contPwrs) + (coF_apDs_pwrs_volsC(3) * contApDs);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_pwrs_vols), box off
% title('Light Intensity and Aperture Diameter Effect on Lesion Volume');
% 
% subplot(1,2,2)
% plot(mdl_apDs_pwrs_volsC), box off
% title('Light Intensity and Aperture Diameter Effect on Lesion Volume (Detected)');

%% 12. Model Light Intensity and Aperture Diameter on Lesion Volume

tbl = [log(colPwrs) colApDs];

mdl_apDs_pwrs_wids = fitlm(tbl, colDepsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Lesion Depths (mm)'})
coF_apDs_pwrs_wids = mdl_apDs_pwrs_wids.Coefficients.Estimate;
reg_apDs_pwrs_vols = coF_apDs_pwrs_vols(1) + (coF_apDs_pwrs_vols(2) * contPwrs) + (coF_apDs_pwrs_vols(3) * contApDs);

mdl_apDs_pwrs_volsC = fitlm(tbl, colVolsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Detected) (mm)'})
coF_apDs_pwrs_volsC = mdl_apDs_pwrs_volsC.Coefficients.Estimate;
reg_apDs_pwrs_volsC = coF_apDs_pwrs_volsC(1) + (coF_apDs_pwrs_volsC(2) * contPwrs) + (coF_apDs_pwrs_volsC(3) * contApDs);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_pwrs_vols), box off
% title('Light Intensity and Aperture Diameter Effect on Lesion Volume');
% 
% subplot(1,2,2)
% plot(mdl_apDs_pwrs_volsC), box off
% title('Light Intensity and Aperture Diameter Effect on Lesion Volume (Detected)');


%% 13. Save r-squared and p vales for all models
r2_OCTADs = {r2_apDs_OCTADs;
    r2_lnpwrs_OCTADs;
    r2_dens_OCTADs;
    mdl_apDs_pwrs_OCTADs.Rsquared.Ordinary;
    []; []};

r2_OCTADsC = {r2_apDs_OCTADsC;
    r2_lnpwrs_OCTADsC;
    r2_dens_OCTADsC;
    mdl_apDs_pwrs_OCTADsC.Rsquared.Ordinary;
    []; []};

r2_vols = {r2_apDs_vols;
    r2_pwrs_vols;
    r2_dens_vols;
    mdl_apDs_pwrs_vols.Rsquared.Ordinary;
    r2_OCTADs_vols;
    r2_OCTADsC_vols};

r2_volsC = {r2_apDs_volsC;
    r2_pwrs_volsC;
    r2_dens_volsC;
    mdl_apDs_pwrs_volsC.Rsquared.Ordinary;
    r2_OCTADs_volsC;
    r2_OCTADsC_volsC};

r2valsTab = table(r2_OCTADs,r2_OCTADsC,r2_vols,r2_volsC,'VariableNames',...
    {'Max_OCTA_Diameter','Max_OCTA_Diameter_Detected','Volume','Volume_Detected'},...
    'RowNames',{'Aperture_Diameter','Light_Intensity','Light_Density',...
    'Aperture_Diameter_Light_Intensity','Max_OCTA_Diameter','Max_OCTA_Diameter_Detected'});

p_OCTADs = {mdl_apDs_OCTADs.Coefficients.pValue';
    mdl_lnpwrs_OCTADs.Coefficients.pValue';
    mdl_dens_OCTADs.Coefficients.pValue';
    mdl_apDs_pwrs_OCTADs.Coefficients.pValue';
    []; []};

p_OCTADsC = {mdl_apDs_OCTADsC.Coefficients.pValue';
    mdl_lnpwrs_OCTADsC.Coefficients.pValue';
    mdl_dens_OCTADsC.Coefficients.pValue';
    mdl_apDs_pwrs_OCTADsC.Coefficients.pValue';
    []; []};

p_vols = {mdl_apDs_vols.Coefficients.pValue';
    mdl_pwrs_vols.Coefficients.pValue';
    mdl_dens_vols.Coefficients.pValue';
    mdl_apDs_pwrs_vols.Coefficients.pValue';
    mdl_OCTADs_vols.Coefficients.pValue';
    mdl_OCTADsC_vols.Coefficients.pValue'};

p_volsC = {mdl_apDs_volsC.Coefficients.pValue';
    mdl_pwrs_volsC.Coefficients.pValue';
    mdl_dens_volsC.Coefficients.pValue';
    mdl_apDs_pwrs_volsC.Coefficients.pValue';
    mdl_OCTADs_volsC.Coefficients.pValue';
    mdl_OCTADsC_volsC.Coefficients.pValue'};

pvalsTab = table(p_OCTADs,p_OCTADsC,p_vols,p_volsC,'VariableNames',...
    {'Max_OCTA_Diameter','Max_OCTA_Diameter_Detected','Volume','Volume_Detected'},...
    'RowNames',{'Aperture_Diameter','Light_Intensity','Light_Density',...
    'Aperture_Diameter_Light_Intensity','Max_OCTA_Diameter','Max_OCTA_Diameter_Detected'});

save('ModelStats.mat','r2valsTab','pvalsTab');
%% 14. Overlay 3D Reconstruction Projection with OCTA Outline and Aperture Diameter

% load binarized image of OCTA image (Monkey C, left hemi, central lesion)
OCTA_bwIm = imread('BinaryC_L.tif');
OCTAbounds = bwboundaries(OCTA_bwIm);
% figure, imshow(OCTA_bwIm), hold on, visboundaries(OCTAbounds), daspect([1 1 1]);

% load lesion matrix
load('PT3_LesionBoundMatrix.mat','interpLesions');
% figure, isosurface(interpLesions), %daspect([1 1 1])

% The voxel volume  33.9 um/pixel in x and y
dx = 0.034; dy = 0.034; dz = 0.45;
dV = dx*dy*dz; % Voxel volume in mm^3

% Now you need to scale it
% Set the scale in the AP axis (450um between each slice)
APScale = 0:dz:dz*(size(interpLesions,3)-1); % in mm

% Set the scale in the ML and DV axis (each pixel is 52.5 um per pixel)
MLScale = 0:dy:dy*(size(interpLesions,2)-1); % in mm
DVScale = 0:dx:dx*(size(interpLesions,1)-1); 

% figure;
% isosurface(MLScale, DVScale, APScale, interpLesions); daspect([1 1 1]);
% xlim([0 100]); ylim([0 100]); zlim([0 100]);
% title('Scaled Interpolated Lesions');
% xlabel('Medial-lateral Axis (mm)'); ylabel('Dorsal-Ventral Axis (mm)'); zlabel('Anterior-Posterior Axis (mm)');

% Isolate central lesion (lesion 4) on left hemisphere
L4X = 24:35; L4Y = 401:532;

recon_CL4 = interpLesions(:,L4Y,L4X);
% figure, isosurface(MLScale(L4Y),DVScale,APScale(L4X),recon_CL4);
% daspect([1 1 1])

%get the correct orientation of the lesion
% first flip it
recon_CL4 = flip(recon_CL4,3);

% figure, isosurface(MLScale(L4Y),DVScale,APScale(L4X),recon_CL4);
% daspect([1 1 1])

flip_CL4 = permute(flipud(recon_CL4), [3 2 1]);
% figure, isosurface(MLScale(L4Y),APScale(L4X),DVScale,flip_CL4);
% daspect([1 1 1])

% [caz,cel] = view;

% obtain stats about each lesion
stats = regionprops('table',flip_CL4,'Centroid');
com = stats.Centroid;

% now rotate lesion so that we're looking at it straight on
AProt = -39;
rotAP_CL4 = imrotate3(flip_CL4,AProt,[0 1 0],'linear');
AP = 0:dz:dz*(size(rotAP_CL4,1)-1);
ML = 0:dy:dy*(size(rotAP_CL4,2)-1);
DV = 0:dx:dx*(size(rotAP_CL4,3)-1);
% figure, isosurface(ML,AP,DV,rotAP_CL4); %view(0,0); 
% ylabel('AP Axis')
% daspect([1 1 1]); view(0,0)

DVrot = -90;
rotDV_CL4 = imrotate3(rotAP_CL4,DVrot,[1 0 0],'linear');
AP = 0:dz:dz*(size(rotDV_CL4,3)-1);
ML = 0:dy:dy*(size(rotDV_CL4,1)-1);
DV = 0:dx:dx*(size(rotDV_CL4,2)-1);
% figure, isosurface(DV,ML,AP,rotDV_CL4); %view(0,0); 
% daspect([1 1 1]); view(0,0)
% zlabel('AP'), ylabel('DV'), xlabel('ML');

MLrot = 90;
rotML_CL4 = imrotate3(rotDV_CL4,MLrot,[0 1 0],'linear'); 
AP = 0:dz:dz*(size(rotML_CL4,2)-1);
ML = 0:dy:dy*(size(rotML_CL4,3)-1);
DV = 0:dx:dx*(size(rotML_CL4,1)-1);
% figure, isosurface(AP,DV,ML,rotML_CL4); %view(0,0);
% daspect([1 1 1]); view(0,0)
% xlabel('AP'), zlabel('ML')

isovalLesions = 0.3481;

% Get the combined volume of the lesions
totalVolume = sum((rotML_CL4 >= isovalLesions),'all') * dV;


% % Get the DV projection
% plane = zeros(size(rotML_CL4)); plane(:,find(DV == 9.928):find(DV == 10.68),:);
volume = sum(rotML_CL4 >= isovalLesions, 1);
volVals = zeros(size(volume,2), size(volume,3));
for i = 1:size(volume,2)
    for j = 1:size(volume,3)
        volVals(i,j) = volume(1,i,j);
    end
end
[x y] = meshgrid(AP,ML); 
% figure; contourf(x,y,volVals','lines','none'); daspect([1 1 1]); 
% colormap([0 0 0; 0 0 1]);

% next, draw aperture circle and add OCTA boundary to scale - piece of cake, right?

volBW = imbinarize(volVals');
% volBW = flip(volBW,1);
volBounds = bwboundaries(volBW,8);

% figure, imshow(volBW,'XData',AP,'YData',ML); daspect([1 1 1]), hold on
% visboundaries(OCTAbounds)
% 
% figure, imshow(volBW,'XData',AP,'YData',ML); %daspect([1 1 1]), hold on


% figure
% plot(volBounds{1}(:,2),volBounds{1}(:,1),'r','linewidth',2), hold on
% plot(OCTAbounds{7}(:,2),OCTAbounds{7}(:,1),'g','linewidth',2), daspect([1 1 1])
 
% scale volBounds to match OCTAbounds

OCTA_res = 21.5; % pixels/mm
ML_res = 1/(0.034); % pixels/mm
AP_res = 1/(0.45);

ML_factor = OCTA_res / ML_res;
AP_factor = OCTA_res / AP_res;

scaledVBounds = volBounds{1};
scaledVBounds(:,1) = volBounds{1}(:,1)*ML_factor;
scaledVBounds(:,2) = volBounds{1}(:,2)*AP_factor;

% figure
% plot(scaledVBounds(:,2),scaledVBounds(:,1),'r','linewidth',2), hold on
% plot(OCTAbounds{7}(:,2),OCTAbounds{7}(:,1),'g','linewidth',2), daspect([1 1 1])

% next, need to align the centers

% find center of mass of scaledVBounds
vcAP = mean(scaledVBounds(:,1));
vcML = mean(scaledVBounds(:,2));

% find center of mass of OCTAbounds
ocAP = mean(OCTAbounds{7}(:,1));
ocML = mean(OCTAbounds{7}(:,2));

% move everything to the origin
cVol = [(scaledVBounds(:,1) - vcAP), (scaledVBounds(:,2) - vcML)];
cOCT = [-(OCTAbounds{7}(:,1) - ocAP), (OCTAbounds{7}(:,2) - ocML)];

% also add aperture diameter (1 mm)
apX = (21.5/2) * cos(linspace(0,2*pi,1000));
apY = (21.5/2) * sin(linspace(0,2*pi,1000));

% figure('color','w')
% plot(cVol(:,2),cVol(:,1),'r','linewidth',2), hold on
% plot(cOCT(:,2),cOCT(:,1),'g','linewidth',2), daspect([1 1 1]),
% plot(apX,apY,'k','linewidth',2), box off, axis off

% maybe smooth the volume boundaries some more
smCVol = [smooth(cVol(:,1),5) smooth(cVol(:,2),20)];

% figure('color','w')
% plot(smCVol(:,2),smCVol(:,1),'r','linewidth',2), hold on
% plot(cOCT(:,2),cOCT(:,1),'g','linewidth',2), daspect([1 1 1]),
% plot(apX,apY,'k','linewidth',2), box off, axis off

% figure('color','w')
% patch(smCVol(:,2),smCVol(:,1),[248 117 117]./255,'FaceAlpha', 0.3,'EdgeColor','none'), hold on
% patch(cOCT(:,2),cOCT(:,1),[18 226 136]./255,'FaceAlpha',0.3,'EdgeColor','none'), daspect([1 1 1]),
% viscircles([0 0], 21.5/2, 'Color','w','LineStyle','--');
% axis off, box off

% separate cVol into just the columns to plot it as a series of lines
fcVol = [cVol(:,2) cVol(:,1)];
fcVol = sortrows(fcVol);

newSlice = false;
numSlices = 1;
slices = {[]};
for i = 1:numel(fcVol)/2
    if i ~=1, prevSlice = fcVol(i-1,1);
        currentSlice = fcVol(i,1);
        if currentSlice == prevSlice
            newSlice = false;
            slices{numSlices} = [slices{numSlices}; fcVol(i,:)];
        else
            newSlice = true;
            numSlices = numSlices + 1;
            slices{numSlices} = [fcVol(i,:)];
        end
    end
end

% figure
% for i = 1:length(slices)
%     plot(slices{1,i}(:,1), slices{1,i}(:,2), 'linewidth',4,'color',[248 117 117]./255), hold on
% end

% figure('color','w')
% for i = 1:length(slices)
%     plot(slices{1,i}(:,1), slices{1,i}(:,2), 'linewidth',4,'color',[[248 117 117]./255 .6]), hold on
% end
% patch(cOCT(:,2),cOCT(:,1),[18 226 136]./255,'FaceAlpha',0.3,'EdgeColor','none'), daspect([1 1 1]),
% viscircles([0 0], 21.5/2, 'Color','w','LineStyle','--');
% axis off, box off
%% 15. Make Composite Figure Drafts for Manuscript - Volume

mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

figure('color','w','units','normalize','outerposition',[0 0 .5 1])
xpos = [.067 .53 ]; ypos = [.6875 .375 .0625]+.04; wh = 0.4; h = 0.25;

% A. projection fig
subplot('position',[xpos(1) ypos(1) wh h])
patch(cOCT(:,2),cOCT(:,1),OCTColor,'FaceAlpha',alpha,'EdgeColor','none'), daspect([1 1 1]), hold on
for i = 1:length(slices)
    plot(slices{1,i}(:,1), slices{1,i}(:,2), 'linewidth',2.5,'color',[volColor alpha]), hold on
end
viscircles([0 0], 21.5/2, 'Color','w','LineStyle','--','linewidth',1);
axis off, box off

% B. vol vs OCTADs
subplot('position',[xpos(2) ypos(1) wh h])
plot(colOCTADs,colVols,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(contOCTADs,reg_OCTADs_vols,'color',mixColor,'linewidth',2),
plot(contOCTADsC,reg_OCTADsC_volsC, 'color',mixColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Volume (mm^3)'), xlabel('Max. OCTA Diameter (mm)')

Str_OCTADs_vols = ['r^2 = ' char(string(r2valsTab{5,3}))  newline 'p = ' char(string(pvalsTab.Volume{5}(numCoFs)))];
Str_OCTADsC_volsC = ['r^2 = ' char(string(r2valsTab{6,4})) newline 'p = ' char(string(pvalsTab.Volume_Detected{6}(numCoFs)))];
text(5, 20, Str_OCTADs_vols, 'fontsize', 8);% set(T, 'fontsize', 8);
text(5, 40, Str_OCTADsC_volsC, 'fontsize', 8);

% C. OCTADs vs apDs
subplot('position',[xpos(1) ypos(2) wh h])
plot(colApDs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',2),
plot(contApDs,reg_apDs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
text(1.75, 1.6, Str_apDs_OCTADs, 'fontsize', 8);
text(1.75, 4.6, Str_apDs_OCTADsC, 'fontsize', 8);

% D. vol vs apDs
subplot('position',[xpos(2) ypos(2) wh h])
plot(colApDs,colVols,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contApDs,reg_apDs_vols,'color',volColor,'linewidth',2),
plot(contApDs,reg_apDs_volsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Lesion Volume (mm^3)'), xlabel('Aperture Diameter (mm)')

Str_apDs_vols = ['r^2 = ' char(string(r2valsTab{1,3}))  newline 'p = ' char(string(pvalsTab.Volume{1}(numCoFs)))];
Str_apDs_volsC = ['r^2 = ' char(string(r2valsTab{1,4})) newline 'p = ' char(string(pvalsTab.Volume_Detected{1}(numCoFs)))];
text(2.1, 10, Str_apDs_vols, 'fontsize', 8);
text(2.1, 25, Str_apDs_volsC, 'fontsize', 8);

% E. OCTADs vs pwr
subplot('position',[xpos(1) ypos(3) wh h])
plot(log(colPwrs),colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(log(contPwrs),reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',2),
plot(log(contPwrs),reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
text(8, 2.5, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
text(5.5, 5, Str_pwrs_OCTADsC, 'fontsize', 8);

% F. vol vs pwr
subplot('position',[xpos(2) ypos(3) wh h])
plot(colPwrs,colVols,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contPwrs,reg_pwrs_vols,'color',volColor,'linewidth',2),
plot(contPwrs,reg_pwrs_volsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Volume (mm^3)'), xlabel('Light Intensity (mW)')

Str_pwrs_vols = ['r^2 = ' char(string(r2valsTab{2,3}))  newline 'p = ' char(string(pvalsTab.Volume{2}(numCoFs)))];
Str_pwrs_volsC = ['r^2 = ' char(string(r2valsTab{2,4})) newline 'p = ' char(string(pvalsTab.Volume_Detected{2}(numCoFs)))];
text(8, 10, Str_pwrs_vols, 'fontsize', 8);% set(T, 'fontsize', 8);
text(5, 30, Str_pwrs_volsC, 'fontsize', 8);

% OCTAD vs apDs, pwr
% subplot('position',[xPos(3) yPos(2) w h]);
% plot3(colPwrs,colApDs,colOCTADs,'.','markersize',10,'color',[14 170 102]./255), hold on
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADs,'color',[14 170 102]./255,'linewidth',2)
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADsC,'color',[14 170 102]./255,'linewidth',2,'linestyle','--')
% legend off, box off, grid on
% xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)'), zlabel('Max. OCTA Diameter (mm)')

% plot(mdl_apDs_pwrs_OCTADs)%,'.','markersize',10)%,'color',[14 170 102]./255)
% set('markersize',10,'color',[14 170 102]./255)
% vol vs apDs, pwr

% make a legend for dashed and solid lines
subplot('position',[.625 .06 0 0])
plot(1,'linestyle','-','color',legColor,'linewidth',1.5), hold on
plot(1,'linestyle','--','color',legColor,'linewidth',1.5)
legend('All Values', 'Detected Values','box','off');
box off


% print -dpdf -painters LesionQuantifications
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16. Model Lesion Depths and Aperture Diameter

mdl_apDs_deps = fitlm(colApDs,colDeps,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_deps = mdl_apDs_deps.Coefficients.Estimate; % store coefficients
if cepts
    reg_apDs_deps = coF_apDs_deps(1) + (coF_apDs_deps(2) * contApDs); % regression line for the model
    est_apDs_deps = coF_apDs_deps(1) + (coF_apDs_deps(2) * colApDs);
else
    reg_apDs_deps = (coF_apDs_deps * contApDs); % regression line for the model
    est_apDs_deps = (coF_apDs_deps * colApDs);
end
r2_apDs_deps = computeR2(colDeps, est_apDs_deps);

mdl_apDs_depsC = fitlm(colApDs,colDepsC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_depsC = mdl_apDs_depsC.Coefficients.Estimate;
if cepts
    reg_apDs_depsC = coF_apDs_depsC(1) + (coF_apDs_depsC(2) * contApDs);
    est_apDs_depsC = coF_apDs_depsC(1) + (coF_apDs_depsC(2) * colApDs);
else
    reg_apDs_depsC = (coF_apDs_depsC * contApDs);
    est_apDs_depsC = (coF_apDs_depsC * colApDs);
end
r2_apDs_depsC = computeR2(colDepsC, est_apDs_depsC);

mdl_apDs_depsV = fitlm(colApDs,colDepsV,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_depsV = mdl_apDs_depsV.Coefficients.Estimate;
if cepts
    reg_apDs_depsV = coF_apDs_depsV(1) + (coF_apDs_depsV(2) * contApDs);
    est_apDs_depsV = coF_apDs_depsV(1) + (coF_apDs_depsV(2) * colApDs);
else
    reg_apDs_depsV = (coF_apDs_depsV * contApDs);
    est_apDs_depsV = (coF_apDs_depsV * colApDs);
end
r2_apDs_depsV = computeR2(colDepsV, est_apDs_depsV);

mdl_apDs_depsVC = fitlm(colApDs,colDepsVC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_depsVC = mdl_apDs_depsVC.Coefficients.Estimate;
if cepts
    reg_apDs_depsVC = coF_apDs_depsVC(1) + (coF_apDs_depsVC(2) * contApDs);
    est_apDs_depsVC = coF_apDs_depsVC(1) + (coF_apDs_depsVC(2) * colApDs);
else
    reg_apDs_depsVC = (coF_apDs_depsVC * contApDs);
    est_apDs_depsVC = (coF_apDs_depsVC * colApDs);
end
r2_apDs_depsVC = computeR2(colDepsVC, est_apDs_depsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_apDs_deps)
% ylabel('Lesion Depth (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Depth vs. Aperture Diameter')
% box off
% anova(mdl_apDs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_apDs_depsC)
% ylabel('Lesion Depth (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Depth (Detected) vs. Aperture Diameter')
% box off
% anova(mdl_apDs_depsC,'Summary');

%% 17. Model Lesion Widths and Aperture Diameter

mdl_apDs_wids = fitlm(colApDs,colWids,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_wids = mdl_apDs_wids.Coefficients.Estimate; % store coefficients
if cepts
    reg_apDs_wids = coF_apDs_wids(1) + (coF_apDs_wids(2) * contApDs); % regression line for the model
    est_apDs_wids = coF_apDs_wids(1) + (coF_apDs_wids(2) * colApDs);
else
    reg_apDs_wids = (coF_apDs_wids * contApDs); % regression line for the model
    est_apDs_wids = (coF_apDs_wids * colApDs);
end
r2_apDs_wids = computeR2(colWids, est_apDs_wids);

mdl_apDs_widsC = fitlm(colApDs,colWidsC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_widsC = mdl_apDs_widsC.Coefficients.Estimate;
if cepts
    reg_apDs_widsC = coF_apDs_widsC(1) + (coF_apDs_widsC(2) * contApDs);
    est_apDs_widsC = coF_apDs_widsC(1) + (coF_apDs_widsC(2) * colApDs);
else
    reg_apDs_widsC = (coF_apDs_widsC * contApDs);
    est_apDs_widsC = (coF_apDs_widsC * colApDs);
end
r2_apDs_widsC = computeR2(colWidsC, est_apDs_widsC);

mdl_apDs_widsV = fitlm(colApDs,colWidsV,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_widsV = mdl_apDs_widsV.Coefficients.Estimate;
if cepts
    reg_apDs_widsV = coF_apDs_widsV(1) + (coF_apDs_widsV(2) * contApDs);
    est_apDs_widsV = coF_apDs_widsV(1) + (coF_apDs_widsV(2) * colApDs);
else
    reg_apDs_widsV = (coF_apDs_widsV * contApDs);
    est_apDs_widsV = (coF_apDs_widsV * colApDs);
end
r2_apDs_widsV = computeR2(colWidsV, est_apDs_widsV);

mdl_apDs_widsVC = fitlm(colApDs,colWidsVC,'RobustOpts',wfun,'intercept',cepts)
coF_apDs_widsVC = mdl_apDs_widsVC.Coefficients.Estimate;
if cepts
    reg_apDs_widsVC = coF_apDs_widsVC(1) + (coF_apDs_widsVC(2) * contApDs);
    est_apDs_widsVC = coF_apDs_widsVC(1) + (coF_apDs_widsVC(2) * colApDs);
else
    reg_apDs_widsVC = (coF_apDs_widsVC * contApDs);
    est_apDs_widsVC = (coF_apDs_widsVC * colApDs);
end
r2_apDs_widsVC = computeR2(colWidsVC, est_apDs_widsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
plot(mdl_apDs_wids)
ylabel('Lesion Width (mm)')
xlabel('Aperture Diameter (mm)')
title('Lesion Width vs. Aperture Diameter')
box off
anova(mdl_apDs_wids,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_apDs_widsC)
% ylabel('Lesion Width (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Width (Detected) vs. Aperture Diameter')
% box off
% anova(mdl_apDs_widsC,'Summary');

%% 18. Model Lesion Widths and Depths Together

mdl_wids_deps = fitlm(colWids,colDeps,'RobustOpts',wfun,'intercept',cepts)
coF_wids_deps = mdl_wids_deps.Coefficients.Estimate; % store coefficients
if cepts
    reg_wids_deps = coF_wids_deps(1) + (coF_wids_deps(2) * contWids); % regression line for the model
    est_wids_deps = coF_wids_deps(1) + (coF_wids_deps(2) * colWids);
else    
    reg_wids_deps = (coF_wids_deps * contWids); % regression line for the model
    est_wids_deps = (coF_wids_deps * colWids);
end
r2_wids_deps = computeR2(colDeps, est_wids_deps);

mdl_widsC_depsC = fitlm(colWidsC,colDepsC,'RobustOpts',wfun,'intercept',cepts)
coF_widsC_depsC = mdl_widsC_depsC.Coefficients.Estimate;
if cepts
    reg_widsC_depsC = coF_widsC_depsC(1) + (coF_widsC_depsC(2) * contWidsC);
    est_widsC_depsC = coF_widsC_depsC(1) + (coF_widsC_depsC(2) * colWidsC);
else
    reg_widsC_depsC = (coF_widsC_depsC * contWidsC);
    est_widsC_depsC = (coF_widsC_depsC * colWidsC);
end
r2_widsC_depsC = computeR2(colDepsC, est_widsC_depsC);

mdl_widsV_depsV = fitlm(colWidsV,colDepsV,'RobustOpts',wfun,'intercept',cepts)
coF_widsV_depsV = mdl_widsV_depsV.Coefficients.Estimate;
if cepts
    reg_widsV_depsV = coF_widsV_depsV(1) + (coF_widsV_depsV(2) * contWidsV);
    est_widsV_depsV = coF_widsV_depsV(1) + (coF_widsV_depsV(2) * colWidsV);
else
    reg_widsV_depsV = (coF_widsV_depsV * contWidsV);
    est_widsV_depsV = (coF_widsV_depsV * colWidsV);
end
r2_widsV_depsV = computeR2(colDepsV, est_widsV_depsV);

mdl_widsVC_depsVC = fitlm(colWidsVC,colDepsVC,'RobustOpts',wfun,'intercept',cepts)
coF_widsVC_depsVC = mdl_widsVC_depsVC.Coefficients.Estimate;
if cepts
    reg_widsVC_depsVC = coF_widsVC_depsVC(1) + (coF_widsVC_depsVC(2) * contWidsVC);
    est_widsVC_depsVC = coF_widsVC_depsVC(1) + (coF_widsVC_depsVC(2) * colWidsVC);
else
    reg_widsVC_depsVC = (coF_widsVC_depsVC * contWidsVC);
    est_widsVC_depsVC = (coF_widsVC_depsVC * colWidsVC);
end
r2_widsVC_depsVC = computeR2(colDepsVC, est_widsVC_depsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_wids_deps)
% ylabel('Lesion Depth (mm)')
% xlabel('Lesion Width (mm)')
% title('Lesion Width vs. Depth')
% box off
% anova(mdl_wids_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_widsC_depsC)
% ylabel('Lesion Depth (mm)')
% xlabel('Lesion Width (mm)')
% title('Lesion Width vs. Depth (Detected)')
% box off
% anova(mdl_widsC_depsC,'Summary');

%% 19. Model Lesion Depths and Light Intensity

mdl_pwrs_deps = fitlm(colPwrs,colDeps,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_deps = mdl_pwrs_deps.Coefficients.Estimate; % store coefficients
if cepts
    reg_pwrs_deps = coF_pwrs_deps(1) + (coF_pwrs_deps(2) * contPwrs); % regression line for the model
    est_pwrs_deps = coF_pwrs_deps(1) + (coF_pwrs_deps(2) * colPwrs);
else
    reg_pwrs_deps = (coF_pwrs_deps * contPwrs); % regression line for the model
    est_pwrs_deps = (coF_pwrs_deps * colPwrs);    
end
r2_pwrs_deps = computeR2(colDeps, est_pwrs_deps);

mdl_pwrs_depsC = fitlm(colPwrs,colDepsC,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_depsC = mdl_pwrs_depsC.Coefficients.Estimate;
if cepts
    reg_pwrs_depsC = coF_pwrs_depsC(1) + (coF_pwrs_depsC(2) * contPwrs);
    est_pwrs_depsC = coF_pwrs_depsC(1) + (coF_pwrs_depsC(2) * colPwrs);
else
    reg_pwrs_depsC = (coF_pwrs_depsC * contPwrs);
    est_pwrs_depsC = (coF_pwrs_depsC * colPwrs);
end
r2_pwrs_depsC = computeR2(colDepsC, est_pwrs_depsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_pwrs_deps)
% ylabel('Lesion Depth (mm)')
% xlabel('Ligth Intensity (mW)')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_pwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_pwrs_depsC)
% ylabel('Lesion Depth (mm)')
% xlabel('Aperture Diameter (mm)')
% title('Lesion Depth (Detected) vs. Aperture Diameter')
% box off
% anova(mdl_pwrs_depsC,'Summary');

%% 20. Model Lesion Depths (no PT4,5) and Light Intensity - Don't Use This

depsBC = depths; depsBC(5:6,8) = NaN;
colDepsBC = reshape(depsBC,numel(depsBC),1);

depsBCC = depNaNs; depsBCC(5:6,8) = NaN;
colDepsBCC = reshape(depsBCC,numel(depsBCC),1);

mdl_pwrs_depsBC = fitlm(colPwrs,colDepsBC)
coF_pwrs_depsBC = mdl_pwrs_depsBC.Coefficients.Estimate; % store coefficients
reg_pwrs_depsBC = coF_pwrs_depsBC(1) + (coF_pwrs_depsBC(2) * contPwrs); % regression line for the model

mdl_pwrs_depsBCC = fitlm(colPwrs,colDepsBCC)
coF_pwrs_depsBCC = mdl_pwrs_depsBCC.Coefficients.Estimate;
reg_pwrs_depsBCC = coF_pwrs_depsBCC(1) + (coF_pwrs_depsBCC(2) * contPwrs);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_pwrs_depsBC)
% ylabel('Lesion Depth (mm)')
% xlabel('Light Intensity (mW)')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_pwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_pwrs_depsBCC)
% ylabel('Lesion Depth (mm)')
% xlabel('Light Intensity (mW)')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_depsC,'Summary');

% subplot(2,2,3)
% h = plotResiduals(mdl_pwrs_depsBC)
% box off
% 
% subplot(2,2,4)
% h = plotResiduals(mdl_pwrs_depsBCC)
% box off

%% 21. Box-Cox Test with Depths

% figure('color','w')
% noNans = [colPwrs, colDepsBC];
% noNans(any(isnan(noNans),2),:) = [];
% newPwrs = noNans(:,1); newDeps = noNans(:,2);
% boxcoxlm(newPwrs,newDeps)
% box off

%% 22. Model Lesion Depths (no PT4,5) and (1/)sqrt(Light Intensity)

depsBC = depths; depsBC(5:6,8) = NaN;
colDepsBC = reshape(depsBC,numel(depsBC),1);

depsBCC = depNaNs; depsBCC(5:6,8) = NaN;
colDepsBCC = reshape(depsBCC,numel(depsBCC),1);

mdl_sqrpwrs_depsBC = fitlm(colPwrs.^(-.5),colDepsBC)
coF_sqrpwrs_depsBC = mdl_sqrpwrs_depsBC.Coefficients.Estimate; % store coefficients
reg_sqrpwrs_depsBC = coF_sqrpwrs_depsBC(1) + (coF_sqrpwrs_depsBC(2) * sqrt(contPwrs)); % regression line for the model

mdl_sqrpwrs_depsBCC = fitlm(colPwrs.^(-.5),colDepsBCC)
coF_sqrpwrs_depsBCC = mdl_sqrpwrs_depsBCC.Coefficients.Estimate;
reg_sqrpwrs_depsBCC = coF_sqrpwrs_depsBCC(1) + (coF_sqrpwrs_depsBCC(2) * sqrt(contPwrs));

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_sqrpwrs_depsBC)
% ylabel('Lesion Depth (mm)')
% xlabel('1/Sqrt(Light Intensity (mW))')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_pwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_sqrpwrs_depsBCC)
% ylabel('Lesion Depth (mm)')
% xlabel('1/Sqrt(Light Intensity (mW))')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_depsC,'Summary');

% subplot(2,2,3)
% plotResiduals(mdl_sqrpwrs_depsBC)
% box off
% 
% subplot(2,2,4)
% plotResiduals(mdl_sqrpwrs_depsBCC)
% box off

%% 23. Model Lesion Depths (no PT4,5) and ln(Light Intensity)

depsBC = depths; %depsBC(5:6,8) = NaN;
colDepsBC = reshape(depsBC,numel(depsBC),1);

depsBCC = depNaNs; %depsBCC(5:6,8) = NaN;
colDepsBCC = reshape(depsBCC,numel(depsBCC),1);

mdl_lnpwrs_depsBC = fitlm(log(colPwrs),colDepsBC,'RobustOpts','bisquare')
coF_lnpwrs_depsBC = mdl_lnpwrs_depsBC.Coefficients.Estimate; % store coefficients
reg_lnpwrs_depsBC = coF_lnpwrs_depsBC(1) + (coF_lnpwrs_depsBC(2) * log(contPwrs)); % regression line for the model

mdl_lnpwrs_depsBCC = fitlm(log(colPwrs),colDepsBCC,'RobustOpts','bisquare')
coF_lnpwrs_depsBCC = mdl_lnpwrs_depsBCC.Coefficients.Estimate;
reg_lnpwrs_depsBCC = coF_lnpwrs_depsBCC(1) + (coF_lnpwrs_depsBCC(2) * log(contPwrs));

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_lnpwrs_depsBC)
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_pwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_lnpwrs_depsBCC)
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_depsC,'Summary');

%% 24. Model Lesion Depths and ln(Light Intensity) - Robust Fit

% wfun = 'huber'; % weight function

% obtain coefficients
[wcoF_lnpwrs_deps s_lnpwrs_deps] = robustfit(log(colPwrs),colDeps,wfun,[],'on');
wreg_lnpwrs_deps = wcoF_lnpwrs_deps(1) + (wcoF_lnpwrs_deps(2) * log(contPwrs)); % regression line for the model

% calculate r squared
sse_lnpwrs_deps = s_lnpwrs_deps.dfe * s_lnpwrs_deps.robust_s^2;
ssr_lnpwrs_deps = norm(wreg_lnpwrs_deps - mean(wreg_lnpwrs_deps))^2;
r_lnpwrs_deps = 1 - (sse_lnpwrs_deps/(sse_lnpwrs_deps + ssr_lnpwrs_deps));
disp(['r^2 = ' num2str(r_lnpwrs_deps) ', p = ' num2str(s_lnpwrs_deps.p(2))])

% obtain coefficients
[wcoF_lnpwrs_depsC s_lnpwrs_depsC] = robustfit(log(colPwrs),colDepsC,wfun,[],'on');
wreg_lnpwrs_depsC = wcoF_lnpwrs_depsC(1) + (wcoF_lnpwrs_depsC(2) * log(contPwrs)); % regression line for the model

% calculate r squared
sse_lnpwrs_depsC = s_lnpwrs_depsC.dfe * s_lnpwrs_depsC.robust_s^2;
ssr_lnpwrs_depsC = norm(wreg_lnpwrs_depsC - mean(wreg_lnpwrs_depsC))^2;
r_lnpwrs_depsC = 1 - (sse_lnpwrs_depsC/(sse_lnpwrs_depsC + ssr_lnpwrs_depsC));
disp(['r^2 = ' num2str(r_lnpwrs_depsC) ', p = ' num2str(s_lnpwrs_depsC.p(2))])

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(contPwrs,wreg_lnpwrs_deps), hold on
% scatter(colPwrs,colDeps,'x')
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% 
% subplot(1,2,2)
% plot(contPwrs,wreg_lnpwrs_depsC), hold on
% scatter(colPwrs,colDepsC,'x')
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off

%% 25. Model Lesion Depths (no PT4,5) and ln(Light Intensity)

depsBC = depths; %depsBC(5:6,8) = NaN;
colDepsBC = reshape(depsBC,numel(depsBC),1);

depsBCC = depNaNs; %depsBCC(5:6,8) = NaN;
colDepsBCC = reshape(depsBCC,numel(depsBCC),1);

mdl_lnpwrs_depsBC = fitlm(log(colPwrs),colDepsBC,'RobustOpts','bisquare')
coF_lnpwrs_depsBC = mdl_lnpwrs_depsBC.Coefficients.Estimate; % store coefficients
reg_lnpwrs_depsBC = coF_lnpwrs_depsBC(1) + (coF_lnpwrs_depsBC(2) * log(contPwrs)); % regression line for the model

mdl_lnpwrs_depsBCC = fitlm(log(colPwrs),colDepsBCC,'RobustOpts','bisquare')
coF_lnpwrs_depsBCC = mdl_lnpwrs_depsBCC.Coefficients.Estimate;
reg_lnpwrs_depsBCC = coF_lnpwrs_depsBCC(1) + (coF_lnpwrs_depsBCC(2) * log(contPwrs));

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_lnpwrs_depsBC)
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_pwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_lnpwrs_depsBCC)
% ylabel('Lesion Depth (mm)')
% xlabel('ln(Light Intensity (mW))')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_depsC,'Summary');

%% 26. Model Lesion Depths and ln(Light Intensity)

% wfun = 'ols'; % weight function

if cepts
    mdl_lnpwrs_deps = fitlm(log(colPwrs),colDeps,'RobustOpts',wfun)
    coF_lnpwrs_deps = mdl_lnpwrs_deps.Coefficients.Estimate; % store coefficients
    reg_lnpwrs_deps = coF_lnpwrs_deps(1) + (coF_lnpwrs_deps(2) * log(contPwrs)); % regression line for the model
    est_lnpwrs_deps = coF_lnpwrs_deps(1) + (coF_lnpwrs_deps(2) * log(colPwrs));
else
    mdl_lnpwrs_deps = fitlm(log(colPwrs+1),colDeps,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_deps = mdl_lnpwrs_deps.Coefficients.Estimate; % store coefficients
    reg_lnpwrs_deps = (coF_lnpwrs_deps * log(contPwrs+1)); % regression line for the model
    est_lnpwrs_deps = (coF_lnpwrs_deps * log(colPwrs+1));
end
r2_lnpwrs_deps = computeR2(colDeps, est_lnpwrs_deps);

if cepts
    mdl_lnpwrs_depsC = fitlm(log(colPwrs),colDepsC,'RobustOpts',wfun)
    coF_lnpwrs_depsC = mdl_lnpwrs_depsC.Coefficients.Estimate;
    reg_lnpwrs_depsC = coF_lnpwrs_depsC(1) + (coF_lnpwrs_depsC(2) * log(contPwrs));
    est_lnpwrs_depsC = coF_lnpwrs_depsC(1) + (coF_lnpwrs_depsC(2) * log(colPwrs));
else
    mdl_lnpwrs_depsC = fitlm(log(colPwrs+1),colDepsC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_depsC = mdl_lnpwrs_depsC.Coefficients.Estimate;
    reg_lnpwrs_depsC = (coF_lnpwrs_depsC * log(contPwrs+1));
    est_lnpwrs_depsC = (coF_lnpwrs_depsC * log(colPwrs+1));
end
r2_lnpwrs_depsC = computeR2(colDepsC, est_lnpwrs_depsC);

if cepts
    mdl_lnpwrs_depsV = fitlm(log(colPwrs),colDepsV,'RobustOpts',wfun)
    coF_lnpwrs_depsV = mdl_lnpwrs_depsV.Coefficients.Estimate;
    reg_lnpwrs_depsV = coF_lnpwrs_depsV(1) + (coF_lnpwrs_depsV(2) * log(contPwrs));
    est_lnpwrs_depsV = coF_lnpwrs_depsV(1) + (coF_lnpwrs_depsV(2) * log(colPwrs));
else
    mdl_lnpwrs_depsV = fitlm(log(colPwrs+1),colDepsV,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_depsV = mdl_lnpwrs_depsV.Coefficients.Estimate;
    reg_lnpwrs_depsV = (coF_lnpwrs_depsV * log(contPwrs+1));
    est_lnpwrs_depsV = (coF_lnpwrs_depsV * log(colPwrs+1));
end
r2_lnpwrs_depsV = computeR2(colDepsV, est_lnpwrs_depsV);

if cepts
    mdl_lnpwrs_depsVC = fitlm(log(colPwrs),colDepsVC,'RobustOpts',wfun)
    coF_lnpwrs_depsVC = mdl_lnpwrs_depsVC.Coefficients.Estimate;
    reg_lnpwrs_depsVC = coF_lnpwrs_depsVC(1) + (coF_lnpwrs_depsVC(2) * log(contPwrs));
    est_lnpwrs_depsVC = coF_lnpwrs_depsVC(1) + (coF_lnpwrs_depsVC(2) * log(colPwrs));
else
    mdl_lnpwrs_depsVC = fitlm(log(colPwrs+1),colDepsVC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_depsVC = mdl_lnpwrs_depsVC.Coefficients.Estimate;
    reg_lnpwrs_depsVC = (coF_lnpwrs_depsVC * log(contPwrs+1));
    est_lnpwrs_depsVC = (coF_lnpwrs_depsVC * log(colPwrs+1));
end
r2_lnpwrs_depsVC = computeR2(colDepsVC, est_lnpwrs_depsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_lnpwrs_deps)
% ylabel('Lesion Depth (mm)')
% xlabel('log(Light Intensity (mW))')
% title('Lesion Depth vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_lnpwrs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_lnpwrs_depsC)
% ylabel('Lesion Depth (mm)')
% xlabel('log(Light Intensity (mW))')
% title('Lesion Depth (Detected) vs. Light Intensity')
% box off
% anova(mdl_lnpwrs_depsC,'Summary');

% [r2 p] = calcMdlStats(mdl_lnpwrs_deps)
% [r2 p] = calcMdlStats(mdl_lnpwrs_depsC)

%% 27. Model Lesion Widths and Light Intensity

mdl_pwrs_wids = fitlm(colPwrs,colWids,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_wids = mdl_pwrs_wids.Coefficients.Estimate; % store coefficients
if cepts
    reg_pwrs_wids = coF_pwrs_wids(1) + (coF_pwrs_wids(2) * contPwrs); % regression line for the model
    est_pwrs_wids = coF_pwrs_wids(1) + (coF_pwrs_wids(2) * colPwrs);
else
    reg_pwrs_wids = (coF_pwrs_wids * contPwrs); % regression line for the model
    est_pwrs_wids = (coF_pwrs_wids * colPwrs);
end
r2_pwrs_wids = computeR2(colWids, est_pwrs_wids);

mdl_pwrs_widsC = fitlm(colPwrs,colWidsC,'RobustOpts',wfun,'intercept',cepts)
coF_pwrs_widsC = mdl_pwrs_widsC.Coefficients.Estimate;
if cepts
    reg_pwrs_widsC = coF_pwrs_widsC(1) + (coF_pwrs_widsC(2) * contPwrs);
    est_pwrs_widsC = coF_pwrs_widsC(1) + (coF_pwrs_widsC(2) * colPwrs);
else
    reg_pwrs_widsC = (coF_pwrs_widsC * contPwrs);
    est_pwrs_widsC = (coF_pwrs_widsC * colPwrs);    
end
r2_pwrs_widsC = computeR2(colWidsC, est_pwrs_widsC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_pwrs_wids)
% ylabel('Lesion Width (mm)')
% xlabel('Ligth Intensity (mW)')
% title('Lesion Width vs. Light Intensity')
% box off
% anova(mdl_pwrs_wids,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_pwrs_widsC)
% ylabel('Lesion Width (mm)')
% xlabel('Ligth Intensity (mW)')
% title('Lesion Width (Detected) vs. Light Intensity')
% box off
% anova(mdl_pwrs_widsC,'Summary');

%% 28. Model Lesion Widths and Sqrt(Light Intensity)

% wfun = 'ols';

mdl_sqrpwrs_wids = fitlm(sqrt(colPwrs),colWids,'RobustOpts',wfun,'intercept',cepts)
coF_sqrpwrs_wids = mdl_sqrpwrs_wids.Coefficients.Estimate; % store coefficients
if cepts
    reg_sqrpwrs_wids = coF_sqrpwrs_wids(1) + (coF_sqrpwrs_wids(2) * sqrt(contPwrs)); % regression line for the model
    est_sqrpwrs_wids = coF_sqrpwrs_wids(1) + (coF_sqrpwrs_wids(2) * sqrt(colPwrs));
else
    reg_sqrpwrs_wids = (coF_sqrpwrs_wids * sqrt(contPwrs)); % regression line for the model
    est_sqrpwrs_wids = (coF_sqrpwrs_wids * sqrt(colPwrs));
end
r2_sqrpwrs_wids = computeR2(colWids, est_sqrpwrs_wids);

mdl_sqrpwrs_widsC = fitlm(sqrt(colPwrs),colWidsC,'RobustOpts',wfun,'intercept',cepts)
coF_sqrpwrs_widsC = mdl_sqrpwrs_widsC.Coefficients.Estimate;
if cepts
    reg_sqrpwrs_widsC = coF_sqrpwrs_widsC(1) + (coF_sqrpwrs_widsC(2) * sqrt(contPwrs));
    est_sqrpwrs_widsC = coF_sqrpwrs_widsC(1) + (coF_sqrpwrs_widsC(2) * sqrt(colPwrs));
else
    reg_sqrpwrs_widsC = (coF_sqrpwrs_widsC * sqrt(contPwrs));
    est_sqrpwrs_widsC = (coF_sqrpwrs_widsC * sqrt(colPwrs));
end
r2_sqrpwrs_widsC = computeR2(colWidsC, est_sqrpwrs_widsC);

figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
subplot(1,2,1)
plot(mdl_sqrpwrs_wids)
ylabel('Lesion Width (mm)')
xlabel('Sqrt(Ligth Intensity (mW)')
title('Lesion Width vs. Light Intensity')
box off
anova(mdl_pwrs_wids,'Summary');

subplot(1,2,2)
plot(mdl_sqrpwrs_widsC)
ylabel('Lesion Width (mm)')
xlabel('Sqrt(Ligth Intensity (mW))')
title('Lesion Width (Detected) vs. Light Intensity')
box off
anova(mdl_pwrs_widsC,'Summary');

%% 29. Model Lesion Widths and ln(Light Intensity)

% wfun = 'ols'; % weight function

if cepts
    mdl_lnpwrs_wids = fitlm(log(colPwrs),colWids,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_wids = mdl_lnpwrs_wids.Coefficients.Estimate; % store coefficients
    reg_lnpwrs_wids = coF_lnpwrs_wids(1) + (coF_lnpwrs_wids(2) * log(contPwrs)); % regression line for the model
    est_lnpwrs_wids = coF_lnpwrs_wids(1) + (coF_lnpwrs_wids(2) * log(colPwrs));
else
    mdl_lnpwrs_wids = fitlm(log(colPwrs+1),colWids,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_wids = mdl_lnpwrs_wids.Coefficients.Estimate; % store coefficients
    reg_lnpwrs_wids = (coF_lnpwrs_wids * log(contPwrs+1)); % regression line for the model
    est_lnpwrs_wids = (coF_lnpwrs_wids * log(colPwrs+1));
end
r2_lnpwrs_wids = computeR2(colWids, est_lnpwrs_wids);

if cepts
    mdl_lnpwrs_widsC = fitlm(log(colPwrs),colWidsC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsC = mdl_lnpwrs_widsC.Coefficients.Estimate;
    reg_lnpwrs_widsC = coF_lnpwrs_widsC(1) + (coF_lnpwrs_widsC(2) * log(contPwrs));
	est_lnpwrs_widsC = coF_lnpwrs_widsC(1) + (coF_lnpwrs_widsC(2) * log(colPwrs));
else
    mdl_lnpwrs_widsC = fitlm(log(colPwrs+1),colWidsC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsC = mdl_lnpwrs_widsC.Coefficients.Estimate;
    reg_lnpwrs_widsC = (coF_lnpwrs_widsC * log(contPwrs+1));
	est_lnpwrs_widsC = (coF_lnpwrs_widsC * log(colPwrs+1));
end
r2_lnpwrs_widsC = computeR2(colWidsC, est_lnpwrs_widsC);

if cepts
    mdl_lnpwrs_widsV = fitlm(log(colPwrs),colWidsV,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsV = mdl_lnpwrs_widsV.Coefficients.Estimate;
    reg_lnpwrs_widsV = coF_lnpwrs_widsV(1) + (coF_lnpwrs_widsV(2) * log(contPwrs));
	est_lnpwrs_widsV = coF_lnpwrs_widsV(1) + (coF_lnpwrs_widsV(2) * log(colPwrs));
else
    mdl_lnpwrs_widsV = fitlm(log(colPwrs+1),colWidsV,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsV = mdl_lnpwrs_widsV.Coefficients.Estimate;
    reg_lnpwrs_widsV = (coF_lnpwrs_widsV * log(contPwrs+1));
	est_lnpwrs_widsV = (coF_lnpwrs_widsV * log(colPwrs+1));
end
r2_lnpwrs_widsV = computeR2(colWidsV, est_lnpwrs_widsV);

if cepts
    mdl_lnpwrs_widsVC = fitlm(log(colPwrs),colWidsVC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsVC = mdl_lnpwrs_widsVC.Coefficients.Estimate;
    reg_lnpwrs_widsVC = coF_lnpwrs_widsVC(1) + (coF_lnpwrs_widsVC(2) * log(contPwrs));
	est_lnpwrs_widsVC = coF_lnpwrs_widsVC(1) + (coF_lnpwrs_widsVC(2) * log(colPwrs));
else
    mdl_lnpwrs_widsVC = fitlm(log(colPwrs+1),colWidsVC,'RobustOpts',wfun,'intercept',cepts)
    coF_lnpwrs_widsVC = mdl_lnpwrs_widsVC.Coefficients.Estimate;
    reg_lnpwrs_widsVC = (coF_lnpwrs_widsVC * log(contPwrs+1));
	est_lnpwrs_widsVC = (coF_lnpwrs_widsVC * log(colPwrs+1));
end
r2_lnpwrs_widsVC = computeR2(colWidsVC, est_lnpwrs_widsVC);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_lnpwrs_wids)
% ylabel('Lesion Width(mm)')
% xlabel('log(Light Intensity (mW))')
% title('Lesion Width vs. Light Intensity')
% % set(gca,'XScale','log')
% box off
% anova(mdl_lnpwrs_wids,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_lnpwrs_widsC)
% ylabel('Lesion Width (mm)')
% xlabel('log(Light Intensity (mW))')
% title('Lesion Width (Detected) vs. Light Intensity')
% box off
% anova(mdl_lnpwrs_widsC,'Summary');

[r2 p] = calcMdlStats(mdl_lnpwrs_wids)
[r2 p] = calcMdlStats(mdl_lnpwrs_widsC)

%% 30. See Correlation between OCTA Diameters and Lesion Widths

mdl_OCTADs_wids = fitlm(colOCTADs,colWids,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_wids = mdl_OCTADs_wids.Coefficients.Estimate;
if cepts
    reg_OCTADs_wids = coF_OCTADs_wids(1) + (coF_OCTADs_wids(2) * contOCTADs);
    est_OCTADs_wids = coF_OCTADs_wids(1) + (coF_OCTADs_wids(2) * colOCTADs);
else
    reg_OCTADs_wids = (coF_OCTADs_wids * contOCTADs);
    est_OCTADs_wids = (coF_OCTADs_wids * colOCTADs);
end
r2_OCTADs_wids = computeR2(colWids, est_OCTADs_wids);

mdl_OCTADsC_wids = fitlm(colOCTADsC,colWids,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_wids = mdl_OCTADsC_wids.Coefficients.Estimate;
if cepts
    reg_OCTADsC_wids = coF_OCTADsC_wids(1) + (coF_OCTADsC_wids(2) * contOCTADsC);
    est_OCTADsC_wids = coF_OCTADsC_wids(1) + (coF_OCTADsC_wids(2) * colOCTADsC);
else
    reg_OCTADsC_wids = (coF_OCTADsC_wids * contOCTADsC);
    est_OCTADsC_wids = (coF_OCTADsC_wids * colOCTADsC);
end
r2_OCTADsC_wids = computeR2(colWids, est_OCTADsC_wids);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_OCTADs_wids)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Volume vs. Max. OCTA Diameter')
% box off
% anova(mdl_OCTADs_wids,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_OCTADsC_wids)
% ylabel('Lesion Volume (mm^3)')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Volume vs. Max. OCTA Diameter (Detected)')
% box off
% anova(mdl_OCTADsC_wids,'Summary');

mdl_OCTADs_widsC = fitlm(colOCTADs,colWidsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_widsC = mdl_OCTADs_widsC.Coefficients.Estimate;
if cepts
    reg_OCTADs_widsC = coF_OCTADs_widsC(1) + (coF_OCTADs_widsC(2) * contOCTADs);
    est_OCTADs_widsC = coF_OCTADs_widsC(1) + (coF_OCTADs_widsC(2) * colOCTADs);
else
    reg_OCTADs_widsC = (coF_OCTADs_widsC * contOCTADs);
    est_OCTADs_widsC = (coF_OCTADs_widsC * colOCTADs);
end
r2_OCTADs_widsC = computeR2(colWidsC, est_OCTADs_widsC);

mdl_OCTADsC_widsC = fitlm(colOCTADsC,colWidsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_widsC = mdl_OCTADsC_widsC.Coefficients.Estimate;
if cepts
    reg_OCTADsC_widsC = coF_OCTADsC_widsC(1) + (coF_OCTADsC_widsC(2) * contOCTADsC);
    est_OCTADsC_widsC = coF_OCTADsC_widsC(1) + (coF_OCTADsC_widsC(2) * colOCTADsC);
else
    reg_OCTADsC_widsC = (coF_OCTADsC_widsC * contOCTADsC);
    est_OCTADsC_widsC = (coF_OCTADsC_widsC * colOCTADsC);
end
r2_OCTADsC_widsC = computeR2(colWidsC,est_OCTADsC_widsC);

mdl_OCTADsV_widsV = fitlm(colOCTADsV,colWidsV,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsV_widsV = mdl_OCTADsV_widsV.Coefficients.Estimate;
if cepts
    reg_OCTADsV_widsV = coF_OCTADsV_widsV(1) + (coF_OCTADsV_widsV(2) * contOCTADsV);
    est_OCTADsV_widsV = coF_OCTADsV_widsV(1) + (coF_OCTADsV_widsV(2) * colOCTADsV);
else
    reg_OCTADsV_widsV = (coF_OCTADsV_widsV * contOCTADsV);
    est_OCTADsV_widsV = (coF_OCTADsV_widsV * colOCTADsV);
end
r2_OCTADsV_widsV = computeR2(colWidsV,est_OCTADsV_widsV);

mdl_OCTADsVC_widsVC = fitlm(colOCTADsVC,colWidsVC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsVC_widsVC = mdl_OCTADsVC_widsVC.Coefficients.Estimate;
if cepts
    reg_OCTADsVC_widsVC = coF_OCTADsVC_widsVC(1) + (coF_OCTADsVC_widsVC(2) * contOCTADsVC);
    est_OCTADsVC_widsVC = coF_OCTADsVC_widsVC(1) + (coF_OCTADsVC_widsVC(2) * colOCTADsVC);
else
    reg_OCTADsVC_widsVC = (coF_OCTADsVC_widsVC * contOCTADsVC);
    est_OCTADsVC_widsVC = (coF_OCTADsVC_widsVC * colOCTADsVC);
end
r2_OCTADsVC_widsVC = computeR2(colWidsVC,est_OCTADsVC_widsVC);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(2,2,1)
% plot(mdl_OCTADs_wids), title('Lesion Width vs Max. OCTA Diameter'), box off
% ylabel('Lesion Width (mm)'), xlabel('')
% 
% subplot(2,2,2)
% plot(mdl_OCTADsC_wids), title('Lesion Width vs Max. OCTA Diameter (Detected)'), box off
% xlabel(''), ylabel('');
% 
% subplot(2,2,3)
% plot(mdl_OCTADs_widsC), title('Lesion Width (Detected) vs Max. OCTA Diameter'), box off
% ylabel('Lesion Width (mm)'), xlabel('Max. OCTA Diameter (mm)')
% 
% subplot(2,2,4)
% plot(mdl_OCTADsC_widsC), title('Lesion Width (Detected) vs Max. OCTA Diameter (Detected)'), box off
% xlabel('Max. OCTA Diameter (mm)'), ylabel('');

%% 31. See Correlation between OCTA Diameters and Lesion Depths

mdl_OCTADs_deps = fitlm(colOCTADs,colDeps,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_deps = mdl_OCTADs_deps.Coefficients.Estimate;
if cepts
    reg_OCTADs_deps = coF_OCTADs_deps(1) + (coF_OCTADs_deps(2) * contOCTADs);
    est_OCTADs_deps = coF_OCTADs_deps(1) + (coF_OCTADs_deps(2) * colOCTADs);
else
    reg_OCTADs_deps = (coF_OCTADs_deps * contOCTADs);
    est_OCTADs_deps = (coF_OCTADs_deps * colOCTADs);
end
r2_OCTADs_deps = computeR2(colDeps, est_OCTADs_deps);

mdl_OCTADsC_deps = fitlm(colOCTADsC,colDeps,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_deps = mdl_OCTADsC_deps.Coefficients.Estimate;
if cepts
    reg_OCTADsC_deps = coF_OCTADsC_deps(1) + (coF_OCTADsC_deps(2) * contOCTADsC);
    est_OCTADsC_deps = coF_OCTADsC_deps(1) + (coF_OCTADsC_deps(2) * colOCTADsC);
else
    reg_OCTADsC_deps = (coF_OCTADsC_deps * contOCTADsC);
    est_OCTADsC_deps = (coF_OCTADsC_deps * colOCTADsC);
end
r2_OCTADsC_deps = computeR2(colDeps, est_OCTADsC_deps);

% figure('color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
% subplot(1,2,1)
% plot(mdl_OCTADs_deps)
% ylabel('Lesion Depth (mm))')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Depth vs. Max. OCTA Diameter')
% box off
% anova(mdl_OCTADs_deps,'Summary');
% 
% subplot(1,2,2)
% plot(mdl_OCTADsC_deps)
% ylabel('Lesion Depth (mm))')
% xlabel('Max. OCTA Diameter (mm)')
% title('Lesion Depth vs. Max. OCTA Diameter (Detected)')
% box off
% anova(mdl_OCTADsC_deps,'Summary');

mdl_OCTADs_depsC = fitlm(colOCTADs,colDepsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADs_depsC = mdl_OCTADs_depsC.Coefficients.Estimate;
if cepts
    reg_OCTADs_depsC = coF_OCTADs_depsC(1) + (coF_OCTADs_depsC(2) * contOCTADs);
    est_OCTADs_depsC = coF_OCTADs_depsC(1) + (coF_OCTADs_depsC(2) * colOCTADs);
else
    reg_OCTADs_depsC = (coF_OCTADs_depsC * contOCTADs);
    est_OCTADs_depsC = (coF_OCTADs_depsC * colOCTADs);
end
r2_OCTADs_depsC = computeR2(colDepsC, est_OCTADs_depsC);

mdl_OCTADsC_depsC = fitlm(colOCTADsC,colDepsC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsC_depsC = mdl_OCTADsC_depsC.Coefficients.Estimate;
if cepts
    reg_OCTADsC_depsC = coF_OCTADsC_depsC(1) + (coF_OCTADsC_depsC(2) * contOCTADsC);
    est_OCTADsC_depsC = coF_OCTADsC_depsC(1) + (coF_OCTADsC_depsC(2) * colOCTADsC);
else
    reg_OCTADsC_depsC = (coF_OCTADsC_depsC * contOCTADsC);
    est_OCTADsC_depsC = (coF_OCTADsC_depsC * colOCTADsC);
end
r2_OCTADsC_depsC = computeR2(colDepsC, est_OCTADsC_depsC);

mdl_OCTADsV_depsV = fitlm(colOCTADsV,colDepsV,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsV_depsV = mdl_OCTADsV_depsV.Coefficients.Estimate;
if cepts
    reg_OCTADsV_depsV = coF_OCTADsV_depsV(1) + (coF_OCTADsV_depsV(2) * contOCTADsV);
    est_OCTADsV_depsV = coF_OCTADsV_depsV(1) + (coF_OCTADsV_depsV(2) * colOCTADsV);
else
    reg_OCTADsV_depsV = (coF_OCTADsV_depsV * contOCTADsV);
    est_OCTADsV_depsV = (coF_OCTADsV_depsV * colOCTADsV);
end
r2_OCTADsV_depsV = computeR2(colDepsV, est_OCTADsV_depsV);

mdl_OCTADsVC_depsVC = fitlm(colOCTADsVC,colDepsVC,'RobustOpts',wfun,'intercept',cepts)
coF_OCTADsVC_depsVC = mdl_OCTADsVC_depsVC.Coefficients.Estimate;
if cepts
    reg_OCTADsVC_depsVC = coF_OCTADsVC_depsVC(1) + (coF_OCTADsVC_depsVC(2) * contOCTADsVC);
    est_OCTADsVC_depsVC = coF_OCTADsVC_depsVC(1) + (coF_OCTADsVC_depsVC(2) * colOCTADsVC);
else
    reg_OCTADsVC_depsVC = (coF_OCTADsVC_depsVC * contOCTADsVC);
    est_OCTADsVC_depsVC = (coF_OCTADsVC_depsVC * colOCTADsVC);
end
r2_OCTADsVC_depsVC = computeR2(colDepsVC, est_OCTADsVC_depsVC);

% figure('color','w','units','normalize','outerposition',[0 0 1 1])
% subplot(2,2,1)
% plot(mdl_OCTADs_deps), title('Lesion Depth vs Max. OCTA Diameter'), box off
% ylabel('Lesion Depth (mm)'), xlabel('')
% 
% subplot(2,2,2)
% plot(mdl_OCTADsC_deps), title('Lesion Depth vs Max. OCTA Diameter (Detected)'), box off
% xlabel(''), ylabel('');
% 
% subplot(2,2,3)
% plot(mdl_OCTADs_depsC), title('Lesion Depth (Detected) vs Max. OCTA Diameter'), box off
% ylabel('Lesion Depth (mm)'), xlabel('Max. OCTA Diameter (mm)')
% 
% subplot(2,2,4)
% plot(mdl_OCTADsC_depsC), title('Lesion Depth (Detected) vs Max. OCTA Diameter (Detected)'), box off
% xlabel('Max. OCTA Diameter (mm)'), ylabel('');

%% 32. Save r-squared and p vales for all models
r2_deps = {r2_apDs_deps;
    r2_pwrs_deps;
    r2_wids_deps;
    r2_OCTADs_deps;
    r2_OCTADsC_deps;
    r2_lnpwrs_deps};

r2_depsC = {r2_apDs_depsC;
    r2_pwrs_depsC;
    r2_widsC_depsC;
    r2_OCTADs_depsC;
    r2_OCTADsC_depsC;
    r2_lnpwrs_depsC};

r2_wids = {r2_apDs_wids;
    r2_pwrs_wids;
    r2_wids_deps;
    r2_OCTADs_wids;
    r2_OCTADsC_wids;
    r2_lnpwrs_wids};

r2_widsC = {r2_apDs_widsC;
    r2_pwrs_widsC;
    r2_widsC_depsC;
    r2_OCTADs_widsC;
    r2_OCTADsC_widsC;
    r2_lnpwrs_widsC};

r2valsWDTab = table(r2_deps,r2_depsC,r2_wids,r2_widsC,'VariableNames',...
    {'Depth','Depth_Detected','Width','Width_Detected'},...
    'RowNames',{'Aperture_Diameter','Light_Intensity','Width'...
    'Max_OCTA_Diameter','Max_OCTA_Diameter_Detected','log(Power)'});

p_deps = {mdl_apDs_deps.Coefficients.pValue';
    mdl_pwrs_deps.Coefficients.pValue';
    mdl_wids_deps.Coefficients.pValue';
    mdl_OCTADs_deps.Coefficients.pValue';
    mdl_OCTADsC_deps.Coefficients.pValue';
    mdl_lnpwrs_deps.Coefficients.pValue'};

p_depsC = {mdl_apDs_depsC.Coefficients.pValue';
    mdl_pwrs_depsC.Coefficients.pValue';
    mdl_widsC_depsC.Coefficients.pValue';
    mdl_OCTADs_depsC.Coefficients.pValue';
    mdl_OCTADsC_depsC.Coefficients.pValue';
    mdl_lnpwrs_depsC.Coefficients.pValue'};

p_wids = {mdl_apDs_wids.Coefficients.pValue';
    mdl_pwrs_wids.Coefficients.pValue';
    mdl_wids_deps.Coefficients.pValue';
    mdl_OCTADs_wids.Coefficients.pValue';
    mdl_OCTADsC_wids.Coefficients.pValue';
    mdl_lnpwrs_wids.Coefficients.pValue'};

p_widsC = {mdl_apDs_widsC.Coefficients.pValue';
    mdl_pwrs_widsC.Coefficients.pValue';
    mdl_widsC_depsC.Coefficients.pValue';
    mdl_OCTADs_widsC.Coefficients.pValue';
    mdl_OCTADsC_widsC.Coefficients.pValue';
    mdl_lnpwrs_widsC.Coefficients.pValue'};

pvalsWDTab = table(p_deps,p_depsC,p_wids,p_widsC,'VariableNames',...
    {'Depth','Depth_Detected','Width','Width_Detected'},...
    'RowNames',{'Aperture_Diameter','Light_Intensity','Width'...
    'Max_OCTA_Diameter','Max_OCTA_Diameter_Detected','log(Power)'});

save('ModelStats.mat','r2valsWDTab','pvalsWDTab');


%% 33. Make Composite Figure Drafts for Manuscript - Depth

mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

figure('color','w','units','normalize','outerposition',[0 0 .5 1])
xpos = [.067 .53 ]; ypos = [.6875 .375 .0625]+.04; wh = 0.4; h = 0.25;

% A. projection fig
subplot('position',[xpos(1) ypos(1) wh h])
patch(cOCT(:,2),cOCT(:,1),OCTColor,'FaceAlpha',alpha,'EdgeColor','none'), daspect([1 1 1]), hold on
for i = 1:length(slices)
    plot(slices{1,i}(:,1), slices{1,i}(:,2), 'linewidth',2.5,'color',[volColor alpha]), hold on
end
viscircles([0 0], 21.5/2, 'Color','w','LineStyle','--','linewidth',1);
axis off, box off

% B. dep vs OCTADs
subplot('position',[xpos(2) ypos(1) wh h])
plot(colOCTADs,colDeps,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(contOCTADs,reg_OCTADs_deps,'color',mixColor,'linewidth',2),
plot(contOCTADsC,reg_OCTADsC_depsC, 'color',mixColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Depth (mm)'), xlabel('Max. OCTA Diameter (mm)')

Str_OCTADs_deps = ['r^2 = ' char(string(r2valsWDTab{4,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{4}(numCoFs)))];
Str_OCTADsC_depsC = ['r^2 = ' char(string(r2valsWDTab{5,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{5}(numCoFs)))];
text(4, 4, Str_OCTADs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
text(4.5, 2, Str_OCTADsC_depsC, 'fontsize', 8);

% C. OCTADs vs apDs
subplot('position',[xpos(1) ypos(2) wh h])
plot(colApDs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',2),
plot(contApDs,reg_apDs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
text(1.75, 1.6, Str_apDs_OCTADs, 'fontsize', 8);
text(1.75, 4.6, Str_apDs_OCTADsC, 'fontsize', 8);

% D. dep vs apDs
subplot('position',[xpos(2) ypos(2) wh h])
plot(colApDs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contApDs,reg_apDs_deps,'color',volColor,'linewidth',2),
plot(contApDs,reg_apDs_depsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Lesion Depth (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_deps = ['r^2 = ' char(string(r2valsWDTab{1,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{1}(numCoFs)))];
Str_apDs_depsC = ['r^2 = ' char(string(r2valsWDTab{1,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{1}(numCoFs)))];
text(1.25, 1, Str_apDs_deps, 'fontsize', 8);
text(2.1, 3.5, Str_apDs_depsC, 'fontsize', 8);

% E. OCTADs vs pwr
subplot('position',[xpos(1) ypos(3) wh h])
plot(colPwrs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(contPwrs,reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',2),
plot(contPwrs,reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
text(5, 3.2, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
text(4, 5.2, Str_pwrs_OCTADsC, 'fontsize', 8);

% F. dep vs pwr
subplot('position',[xpos(2) ypos(3) wh h])
plot(colPwrs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contPwrs,reg_lnpwrs_deps,'color',volColor,'linewidth',2),
plot(contPwrs,reg_lnpwrs_depsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Depth (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_deps = ['r^2 = ' char(string(r2valsWDTab{6,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{6}(numCoFs)))];
Str_pwrs_depsC = ['r^2 = ' char(string(r2valsWDTab{6,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{6}(numCoFs)))];
text(5, 2, Str_pwrs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
text(5, 3.2, Str_pwrs_depsC, 'fontsize', 8);

% OCTAD vs apDs, pwr
% subplot('position',[xPos(3) yPos(2) w h]);
% plot3(colPwrs,colApDs,colOCTADs,'.','markersize',10,'color',[14 170 102]./255), hold on
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADs,'color',[14 170 102]./255,'linewidth',2)
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADsC,'color',[14 170 102]./255,'linewidth',2,'linestyle','--')
% legend off, box off, grid on
% xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)'), zlabel('Max. OCTA Diameter (mm)')

% plot(mdl_apDs_pwrs_OCTADs)%,'.','markersize',10)%,'color',[14 170 102]./255)
% set('markersize',10,'color',[14 170 102]./255)
% vol vs apDs, pwr

% make a legend for dashed and solid lines
subplot('position',[.625 .06 0 0])
plot(1,'linestyle','-','color',legColor,'linewidth',1.5), hold on
plot(1,'linestyle','--','color',legColor,'linewidth',1.5)
legend('All Values', 'Detected Values','box','off');
box off


print -dpdf -painters LesionQuantificationsDepth

%% 34. Make Composite Figure Drafts for Manuscript - Width

mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

figure('color','w','units','normalize','outerposition',[0 0 .5 1])
xpos = [.067 .53 ]; ypos = [.6875 .375 .0625]+.04; wh = 0.4; h = 0.25;

% A. projection fig
subplot('position',[xpos(1) ypos(1) wh h])
patch(cOCT(:,2),cOCT(:,1),OCTColor,'FaceAlpha',alpha,'EdgeColor','none'), daspect([1 1 1]), hold on
for i = 1:length(slices)
    plot(slices{1,i}(:,1), slices{1,i}(:,2), 'linewidth',2.5,'color',[volColor alpha]), hold on
end
viscircles([0 0], 21.5/2, 'Color','w','LineStyle','--','linewidth',1);
axis off, box off

% B. wid vs OCTADs
subplot('position',[xpos(2) ypos(1) wh h])
plot(colOCTADs,colWids,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(contOCTADs,reg_OCTADs_wids,'color',mixColor,'linewidth',2),
plot(contOCTADsC,reg_OCTADsC_widsC, 'color',mixColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Max. OCTA Diameter (mm)')

Str_OCTADs_wids = ['r^2 = ' char(string(r2valsWDTab{4,3}))  newline 'p = ' char(string(pvalsWDTab.Width{4}(numCoFs)))];
Str_OCTADsC_widsC = ['r^2 = ' char(string(r2valsWDTab{5,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{5}(numCoFs)))];
text(3.4, 5.5, Str_OCTADs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
text(4.5, 3.45, Str_OCTADsC_widsC, 'fontsize', 8);

% C. OCTADs vs apDs
subplot('position',[xpos(1) ypos(2) wh h])
plot(colApDs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',2),
plot(contApDs,reg_apDs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
text(1.75, 1.6, Str_apDs_OCTADs, 'fontsize', 8);
text(1.75, 4.6, Str_apDs_OCTADsC, 'fontsize', 8);

% D. wid vs apDs
subplot('position',[xpos(2) ypos(2) wh h])
plot(colApDs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contApDs,reg_apDs_wids,'color',volColor,'linewidth',2),
plot(contApDs,reg_apDs_widsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_wids = ['r^2 = ' char(string(r2valsWDTab{1,3}))  newline 'p = ' char(string(pvalsWDTab.Width{1}(numCoFs)))];
Str_apDs_widsC = ['r^2 = ' char(string(r2valsWDTab{1,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{1}(numCoFs)))];
text(1.75, 2, Str_apDs_wids, 'fontsize', 8);
text(1.75, 5, Str_apDs_widsC, 'fontsize', 8);

% E. OCTADs vs pwr
subplot('position',[xpos(1) ypos(3) wh h])
plot(colPwrs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(contPwrs,reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',2),
plot(contPwrs,reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Max. OCTA Diameter (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
text(5, 3.2, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
text(4, 5.2, Str_pwrs_OCTADsC, 'fontsize', 8);

% F. wid vs pwr
subplot('position',[xpos(2) ypos(3) wh h])
plot(colPwrs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
plot(contPwrs,reg_lnpwrs_wids,'color',volColor,'linewidth',2),
plot(contPwrs,reg_lnpwrs_widsC, 'color',volColor,'linewidth',2,'linestyle','--')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_wids = ['r^2 = ' char(string(r2valsWDTab{6,3}))  newline 'p = ' char(string(pvalsWDTab.Width{6}(numCoFs)))];
Str_pwrs_widsC = ['r^2 = ' char(string(r2valsWDTab{6,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{6}(numCoFs)))];
text(5, 3.5, Str_pwrs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
text(5, 5.5, Str_pwrs_widsC, 'fontsize', 8);

% OCTAD vs apDs, pwr
% subplot('position',[xPos(3) yPos(2) w h]);
% plot3(colPwrs,colApDs,colOCTADs,'.','markersize',10,'color',[14 170 102]./255), hold on
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADs,'color',[14 170 102]./255,'linewidth',2)
% plot3(contPwrs,contApDs,reg_apDs_pwrs_OCTADsC,'color',[14 170 102]./255,'linewidth',2,'linestyle','--')
% legend off, box off, grid on
% xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)'), zlabel('Max. OCTA Diameter (mm)')

% plot(mdl_apDs_pwrs_OCTADs)%,'.','markersize',10)%,'color',[14 170 102]./255)
% set('markersize',10,'color',[14 170 102]./255)
% vol vs apDs, pwr

% make a legend for dashed and solid lines
subplot('position',[.625 .06 0 0])
plot(1,'linestyle','-','color',legColor,'linewidth',1.5), hold on
plot(1,'linestyle','--','color',legColor,'linewidth',1.5)
legend('All Values', 'Detected Values','box','off');
box off


print -dpdf -painters LesionQuantificationsWidths

%% 35. Create a Figure Showing All of Lesion Info for Each Hemisphere

% circles with orange edges show aperture diameter
% fill color of orange circle shows the light intensity on a log scale
% OCTA diameters represented by black circles
% fill color of black circles is redblue scale showing lesion depth

% downside of this type of figure is that for lesions that have a depth
% measurement but no OCTA diameter, we can't see it

% make map of holes in XY plane 
x = [6 6 12.5 12.5 12.5 19 19 12.5];
y = [16.625 8.375 20.25 12.5 4.75 16.625 8.375 12.5];
[xq, yq] = meshgrid(0:0.005:25,0:0.005:25);

nBin = 1000;

histColor = redblue(2*nBin+1); histColor = histColor(round((2*nBin+1)/2):end,:);
% divide depth values into bins
rbDeps = round((depths - min(depths,[],'all')) * nBin / (max(depths,[],'all') - min(depths,[],'all')))+1;

yGrad = [linspace(255,255,nBin+1);linspace(255,255,nBin+1);linspace(255,100,nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(pwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

xpos = [.25 .5 .25 .5 .25 .5];
ypos = [.8 .8 .6 .6 .4 .2];
wh = .2;
figure('color','w','units','normalize','outerposition',[.3 0 .3 1])
for iHemi = 1:size(apDs,1)
    
    subplot('position',[xpos(iHemi) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),octaMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor','k','facealpha',1);
        end
        
        % Show diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),apDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim([0 24])
        ylim([0 24])
        set(g,'edgecolor',[240 160 0]./255,'facealpha',1);
        axis off
    else
        for iAp = 1:size(apDs,2)

            % show max OCTA diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),octaMaxDs(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor','k','facealpha',1);
            end

            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),apDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim([0 25])
            ylim([0 24])
            set(g,'edgecolor',[240 160 0]./255,'facealpha',1);
            axis off
        end
    end
    switch iHemi
        case 1, title('Left','fontsize',11)
        case 2, title('Right','fontsize',11)
    end
end

% add scale bar
plot(linspace(10,20,100),linspace(2,2,100),'linewidth',1.5,'color','k')

[~,ht] = suplabel('Monkey','t',[-.25 .1323 .84 .84]);
set(ht,'fontsize',11);

[~,hb] = suplabel('B','t',[-.25 .06 .84 .84]);
set(hb,'fontsize',11);

[~,hc] = suplabel('C','t',[-.25 -.14 .84 .84]);
set(hc,'fontsize',11);

[~,hc] = suplabel('D','t',[-.25 -.34 .84 .84]);
set(hc,'fontsize',11);

[~,hc] = suplabel('E','t',[-.25 -.54 .84 .84]);
set(hc,'fontsize',11);

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1) ypos(6)-.03 2*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(2)+wh ypos(6)+.05 .033 3.5*wh];
caxis([min(depths,[],'all') max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 2.5;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

%% 36. Create Another Figure Showing All of Lesion Info for Each Hemisphere

% circles with orange edges show aperture diameter
% fill color of orange circle shows the light intensity on a log scale
% OCTA diameters represented by green circles
% black circles represent histology-measured width of lesions
% fill color of black circles is redblue scale showing lesion depth

% orientation of monkey C left hemisphere needs to be flipped
rePwrs = pwrs;
rePwrs(3,1:2) = pwrs(3,6:7); rePwrs(3,6:7) = pwrs(3,1:2);

reApDs = apDs;
reApDs(3,1:2) = apDs(3,6:7); reApDs(3,6:7) = apDs(3,1:2);

reOCTAMaxDs = octaMaxDs;
reOCTAMaxDs(3,1:2) = octaMaxDs(3,6:7); reOCTAMaxDs(3,6:7) = octaMaxDs(3,1:2);

reWidths = widths;
reWidths(3,1:2) = widths(3,6:7); reWidths(3,6:7) = widths(3,1:2);

reDepths = depths;
reDepths(3,1:2) = depths(3,6:7); reDepths(3,6:7) = depths(3,1:2);

% make map of holes in XY plane 
x = [6 6 12.5 12.5 12.5 19 19 12.5];
y = [16.625 8.375 20.25 12.5 4.75 16.625 8.375 12.5];
[xq, yq] = meshgrid(0:0.005:25,0:0.005:25);

nBin = 1000;

OCTColor = [14 170 102]./255; %[18 226 136]./255;

histColor = redblue(2*nBin+1); histColor = histColor(round((2*nBin+1)/2):end,:);
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;

yGrad = [linspace(255,255,nBin+1);linspace(255,255,nBin+1);linspace(255,100,nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(rePwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

xpos = [.25 .5 .25 .5 .25 .5];
ypos = [.8 .8 .6 .6 .4 .2];
wh = .2;
figure('color','w','units','normalize','outerposition',[.3 0 .3 1])
for iHemi = 1:size(reApDs,1)
    
    subplot('position',[xpos(iHemi) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        
        % show histology depth and width
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor','r','facealpha',1);
        end
        
        % show OCTA diameter
        [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
        h = fill(xunit,yunit,OCTColor); hold on
        set(h,'edgecolor',OCTColor,'facealpha',.3);
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim([0 24])
        ylim([0 24])
        set(g,'edgecolor',[240 160 0]./255,'facealpha',1);
        axis off
    else
        for iAp = 1:size(apDs,2)

            % show max histo lesion diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor','r','facealpha',1);
            end

            % show OCTA diameter
            [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,OCTColor); hold on
            set(h,'edgecolor',OCTColor,'facealpha',.3);
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim([0 25])
            ylim([0 24])
            set(g,'edgecolor',[240 160 0]./255,'facealpha',1);
            axis off
        end
    end
    switch iHemi
        case 1, title('Left','fontsize',11)
        case 2, title('Right','fontsize',11)
    end
end

% add scale bar
plot(linspace(10,20,100),linspace(2,2,100),'linewidth',1.5,'color','k')

[~,ht] = suplabel('Monkey','t',[-.25 .1323 .84 .84]);
set(ht,'fontsize',11);

[~,hb] = suplabel('B','t',[-.25 .06 .84 .84]);
set(hb,'fontsize',11);

[~,hc] = suplabel('C','t',[-.25 -.14 .84 .84]);
set(hc,'fontsize',11);

[~,hc] = suplabel('D','t',[-.25 -.34 .84 .84]);
set(hc,'fontsize',11);

[~,hc] = suplabel('E','t',[-.25 -.54 .84 .84]);
set(hc,'fontsize',11);

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1) ypos(6)-.03 2*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(2)+wh ypos(6)+.05 .033 3.5*wh];
caxis([min(depths,[],'all') max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 2.5;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

print -dpdf -painters TableFigure

%% 37. Create Another Figure Showing All of Lesion Info for Each Hemisphere

% circles with orange edges show aperture diameter
% fill color of orange circle shows the light intensity on a log scale
% OCTA diameters represented by green circles
% black circles represent histology-measured width of lesions
% fill color of black circles is redblue scale showing lesion depth

% OCTA and Histology circles will not overlap

% orientation of monkey C left hemisphere needs to be flipped
rePwrs = pwrs;
rePwrs(3,1:2) = pwrs(3,6:7); rePwrs(3,6:7) = pwrs(3,1:2);

reApDs = apDs;
reApDs(3,1:2) = apDs(3,6:7); reApDs(3,6:7) = apDs(3,1:2);

reOCTAMaxDs = octaMaxDs;
reOCTAMaxDs(3,1:2) = octaMaxDs(3,6:7); reOCTAMaxDs(3,6:7) = octaMaxDs(3,1:2);

reWidths = widths;
reWidths(3,1:2) = widths(3,6:7); reWidths(3,6:7) = widths(3,1:2);

reDepths = depths;
reDepths(3,1:2) = depths(3,6:7); reDepths(3,6:7) = depths(3,1:2);

% make map of holes in XY plane 
x = [6 6 12.5 12.5 12.5 19 19 12.5];
y = [16.625 8.375 20.25 12.5 4.75 16.625 8.375 12.5];
[xq, yq] = meshgrid(0:0.005:25,0:0.005:25);

nBin = 1000;

OCTColor = [14 170 102]./255; %[18 226 136]./255;
% volColor = [248 117 117]./255;
volColor = [194 10 10]./255;
yelColor = [255 183 21]./255;
% yelColor = [255 210 21]./255;

% rb = redblue(2*nBin+1); rb = rb(round((2*nBin+1)/2):end,:); % red histo lesion
histColor = [linspace(255,volColor(1)*255,nBin+1);linspace(255,volColor(2)*255,nBin+1);linspace(255,volColor(3)*255,nBin+1)]'./255;
reDepths(reDepths == 0) = min(colDeps(colDeps>0)); % remove zero values for color scale
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;


yGrad = [linspace(255,yelColor(1)*255,nBin+1);linspace(255,yelColor(2)*255,nBin+1);linspace(255,yelColor(3)*255,nBin+1)]'./255;
% yGrad = [linspace(255,yelColor(1),nBin+1);linspace(255,yelColor(2),nBin+1);linspace(255,yelColor(3),nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(rePwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

xpos = [.1 .3; .5 .7; .1 .3; .5 .7; .1 .3; .5 .7];
ypos = [.75 .75 .55 .55 .35 .15]-.05;
wh = .2;
figure('color','w','units','normalize','outerposition',[.2 0 .6 1])
for iHemi = 1:size(reApDs,1)
    
    % Show OCTA first
    subplot('position',[xpos(iHemi,1) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5
        iAp = 8;
        
        % show OCTA diameter
        [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
        h = fill(xunit,yunit,OCTColor); hold on
        set(h,'edgecolor',OCTColor,'facealpha',1);
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim([0 24])
        ylim([0 24])
        set(g,'edgecolor',yelColor,'facealpha',1);
        axis off
    elseif iHemi == 6
        axis off
    else
        for iAp = 1:size(apDs,2)-1

            % show OCTA diameter
            [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,OCTColor); hold on
            set(h,'edgecolor',OCTColor,'facealpha',1);
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim([0 24])
            ylim([0 24])
            set(g,'edgecolor',yelColor,'facealpha',1);
            axis off
        end
    end
    
    % show histology second
    subplot('position',[xpos(iHemi,2) ypos(iHemi) wh wh], 'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        
        % show histology depth and width
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor','none','facealpha',1);
        end
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim([0 24])
        ylim([0 24])
        set(g,'edgecolor',yelColor,'facealpha',1);
        axis off
        
    else
        for iAp = 1:size(apDs,2)-1

            % show max histo lesion diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor','none','facealpha',1);
            end
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim([0 24])
            ylim([0 24])
            set(g,'edgecolor',yelColor,'facealpha',1);
            axis off
        end
    end
end

% add scale bar (1 cm)
plot(linspace(9,19,100),linspace(2,2,100),'linewidth',1.5,'color','k')

[~,ht] = suplabel('Monkey','t',[-.35 .135 .84 .84]);
set(ht,'fontsize',11,'fontweight', 'normal');

[~,hb] = suplabel('B','t',[-.35 -.04 .84 .84]);
set(hb,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('C','t',[-.35 -.24 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('D','t',[-.35 -.44 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('E','t',[-.35 -.64 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hl] = suplabel('Left','t',[-.119 .135 .84 .84]);
set(hl,'fontsize',11,'fontweight','normal');

[~,hr] = suplabel('Right','t',[.285 .135 .84 .84]);
set(hr,'fontsize',11,'fontweight','normal');

[~,ho1] = suplabel('OCTA','t',[-.2215 .1 .84 .84]);
set(ho1,'fontsize',11,'fontweight','normal');

[~,ho2] = suplabel('OCTA','t',[.1823 .1 .84 .84]);
set(ho2,'fontsize',11,'fontweight','normal');

[~,hh1] = suplabel('Histology','t',[-.015 .1 .84 .84]);
set(hh1,'fontsize',11,'fontweight','normal');

[~,hh2] = suplabel('Histology','t',[.382 .1 .84 .84]);
set(hh2,'fontsize',11,'fontweight','normal');

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1)+.05 ypos(6)-.02 3.5*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(end)+wh ypos(6)+.05 .02 3.5*wh];
caxis([min(colDeps(colDeps>0)) max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 1.7;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

print -dpdf -painters TableFigurev3

%% 38. Show Example Lesion to Explain the Figure in the Previous Section

% show PT3 left lesion 4
iHemi = 3;
iAp = 4;

figure('color','white')

subplot('position',[0 .3 .5 .5],'color','none')
% show OCTA diameter
[xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
h = fill(xunit,yunit,OCTColor); hold on
set(h,'edgecolor',OCTColor,'facealpha',1);

% Show diameter and intensity of illumination
[xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
hold on
daspect([1 1 1])
xlim([10 16])
ylim([10 16])
set(g,'edgecolor',yelColor,'facealpha',1);
axis off
title('Max. OCTA Diameter','fontweight','normal')

subplot('position',[.5 .3 .5 .5],'color','none')
% show max histo lesion diameter with color showing depth
if ~isnan(depths(iHemi,iAp))
    [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
    h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
    set(h,'edgecolor','none','facealpha',1);
end

% Show diameter and intensity of illumination
[xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
hold on
daspect([1 1 1])
xlim([10 16])
ylim([10 16])
set(g,'edgecolor',yelColor,'facealpha',1);
axis off
title('Histological Diameter','fontweight','normal')

% add scale bar (5 mm)
plot(linspace(14,19,100),linspace(10,10,100),'linewidth',1.5,'color','k')

print -dpdf -painters SampleLesionTableFigv3

%% 39. Put Everything into one Figure (Composite of Composites)

figure('color','w','units','normalize','outerposition',[0 0 1 1])

OCTColor = [14 170 102]./255; %[18 226 136]./255;
% volColor = [248 117 117]./255;
volColor = [194 10 10]./255;
yelColor = [255 183 21]./255;
% yelColor = [255 210 21]./255;

% rb = redblue(2*nBin+1); rb = rb(round((2*nBin+1)/2):end,:); % red histo lesion
histColor = [linspace(255,volColor(1)*255,nBin+1);linspace(255,volColor(2)*255,nBin+1);linspace(255,volColor(3)*255,nBin+1)]'./255;
reDepths(reDepths == 0) = min(colDeps(colDeps>0)); % remove zero values for color scale
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;


yGrad = [linspace(255,yelColor(1)*255,nBin+1);linspace(255,yelColor(2)*255,nBin+1);linspace(255,yelColor(3)*255,nBin+1)]'./255;
% yGrad = [linspace(255,yelColor(1),nBin+1);linspace(255,yelColor(2),nBin+1);linspace(255,yelColor(3),nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(rePwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

% A. Circles figure covering the first half of figure
xpos = [.1 .3; .5 .7; .1 .3; .5 .7; .1 .3; .5 .7]./2;
ypos = [.75 .75 .55 .55 .35 .15];
wh = .2/2;

lims = [1 23];

cw = 0.25; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually

for iHemi = 1:size(reApDs,1)
    
    % Show OCTA first
    subplot('position',[xpos(iHemi,1) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5
        iAp = 8;
        
        % show OCTA diameter
        [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
        h = fill(xunit,yunit,OCTColor); hold on
        set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
    elseif iHemi == 6
        axis off
    else
        for iAp = 1:size(apDs,2)-1

            % show OCTA diameter
            [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,OCTColor); hold on
            set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
    
    % show histology second
    subplot('position',[xpos(iHemi,2) ypos(iHemi) wh wh], 'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        
        % show histology depth and width
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
        end
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
        
    else
        for iAp = 1:size(apDs,2)-1

            % show max histo lesion diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
            end
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
end

% add scale bar (1 cm)
plot(linspace(2.5,12.5,100),linspace(3,3,100),'linewidth',1.5,'color','k')

% now do the plots
mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

xpos1 = [.55 .7 .85]; ypos1 = [.6875 .375 .0625]+.06; wh1 = 0.1; h1 = 0.2;

lw = 1.5; % linewidth

% B. octa vs wid
subplot('position',[xpos1(1) ypos1(1) wh1 h1],'color','none')
% plot(colOCTADs,colWids,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(colOCTADsC,colWidsC,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(colOCTADsV,colWidsV,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(contOCTADs,reg_OCTADs_wids,'color',mixColor,'linewidth',lw,'linestyle','--'),
% plot(contOCTADsC,reg_OCTADsC_widsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
plot(contOCTADsV,reg_OCTADsV_widsV, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 6.2])
ylabel('Lesion Width (mm)'), xlabel('OCTA Diameter (mm)')

% Str_OCTADs_wids = ['r^2 = ' char(string(r2valsWDTab{4,3}))  newline 'p = ' char(string(pvalsWDTab.Width{4}(numCoFs)))];
% Str_OCTADsC_widsC = ['r^2 = ' char(string(r2valsWDTab{5,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{5}(numCoFs)))];
Str_OCTADsV_widsV = ['r^2 = ' num2str(r2_OCTADsV_widsV) newline 'p = ' num2str(mdl_OCTADsV_widsV.Coefficients.pValue(2))];
% text(2, 6, Str_OCTADs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2, 5, Str_OCTADsC_widsC, 'fontsize', 8);
text(2, 5, Str_OCTADsV_widsV, 'fontsize', 8);

% C. dep vs OCTADs
subplot('position',[xpos1(2) ypos1(1) wh1 h1],'color','none')
% plot(colOCTADs,colDeps,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(colOCTADsC,colDepsC,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(colOCTADsV,colDepsV,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(contOCTADs,reg_OCTADs_deps,'color',mixColor,'linewidth',lw,'linestyle','--'),
% plot(contOCTADsC,reg_OCTADsC_depsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
plot(contOCTADsV,reg_OCTADsV_depsV, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('OCTA Diameter (mm)')

% Str_OCTADs_deps = ['r^2 = ' char(string(r2valsWDTab{4,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{4}(numCoFs)))];
% Str_OCTADsC_depsC = ['r^2 = ' char(string(r2valsWDTab{5,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{5}(numCoFs)))];
Str_OCTADsV_depsV = ['r^2 = ' num2str(r2_OCTADsV_depsV) newline 'p = ' num2str(mdl_OCTADsV_depsV.Coefficients.pValue(2))];
% text(4, 4.2, Str_OCTADs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2.5, 0.75, Str_OCTADsC_depsC, 'fontsize', 8);
text(2.5, 0.75, Str_OCTADsV_depsV, 'fontsize', 8);

% D. wids vs deps
subplot('position',[xpos1(3) ypos1(1) wh1 h1],'color','none')
% plot(colWids,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colWidsC,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contWids,reg_wids_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
plot(contWidsC,reg_widsC_depsC,'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Lesion Width (mm)')

% Str_wids_deps = ['r^2 = ' char(string(r2valsWDTab{3,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{3}(numCoFs)))];
Str_widsC_depsC = ['r^2 = ' char(string(r2valsWDTab{3,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{3}(numCoFs)))];
% text(1.2, .5, Str_wids_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
text(2.5, 0.75, Str_widsC_depsC, 'fontsize', 8);

% E. OCTADs vs apDs
subplot('position',[xpos1(1) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(colApDs,colOCTADsC,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(colApDs,colOCTADsV,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
% Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
Str_apDs_OCTADsV = ['r^2 = ' num2str(r2_apDs_OCTADsV) newline 'p = ' num2str(mdl_apDs_OCTADsV.Coefficients.pValue(2))];
% text(1.5, 1.3, Str_apDs_OCTADs, 'fontsize', 8);
% text(1.3, 1, Str_apDs_OCTADsC, 'fontsize', 8);
text(1.3, 1, Str_apDs_OCTADsV, 'fontsize', 8);

% F. wids vs apDs
subplot('position',[xpos1(2) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colApDs,colWidsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colApDs,colWidsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contApDs,reg_apDs_wids,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_widsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_wids = ['r^2 = ' char(string(r2valsWDTab{1,3}))  newline 'p = ' char(string(pvalsWDTab.Width{1}(numCoFs)))];
% Str_apDs_widsC = ['r^2 = ' char(string(r2valsWDTab{1,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{1}(numCoFs)))];
Str_apDs_widsV = ['r^2 = ' num2str(r2_apDs_widsV) newline 'p = ' num2str(mdl_apDs_widsV.Coefficients.pValue(numCoFs))];
% text(1.5, 1.75, Str_apDs_wids, 'fontsize', 8);
% text(1.5, 1, Str_apDs_widsC, 'fontsize', 8);
text(1.5, 1, Str_apDs_widsV, 'fontsize', 8);

% G. deps vs apDs
subplot('position',[xpos1(3) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colApDs,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colApDs,colDepsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contApDs,reg_apDs_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_depsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_deps = ['r^2 = ' char(string(r2valsWDTab{1,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{1}(numCoFs)))];
% Str_apDs_depsC = ['r^2 = ' char(string(r2valsWDTab{1,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{1}(numCoFs)))];
Str_apDs_depsV = ['r^2 = ' num2str(r2_apDs_depsV) newline 'p = ' num2str(mdl_apDs_depsV.Coefficients.pValue(numCoFs))];
% text(1.25, 1, Str_apDs_deps, 'fontsize', 8);
% text(.25, 3.5, Str_apDs_depsC, 'fontsize', 8);
text(.25, 3.5, Str_apDs_depsV, 'fontsize', 8);

% H. OCTADs vs pwr
subplot('position',[xpos1(1) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(colPwrs,colOCTADsC,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(log(colPwrs),colOCTADsV,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(contPwrs,reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',lw,'linestyle','-')
plot(log(contPwrs),reg_lnpwrs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Log(Intensity (mW))')

% Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
% Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
Str_pwrs_OCTADsV = ['r^2 = ' num2str(r2_pwrs_OCTADsV) newline 'p = ' num2str(mdl_pwrs_OCTADsV.Coefficients.pValue(numCoFs))];
% text(5, 3.2, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_OCTADsC, 'fontsize', 8);
text(-1, 1, Str_pwrs_OCTADsV, 'fontsize', 8);

% I. wids vs pwr
subplot('position',[xpos1(2) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colPwrs,colWidsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(log(colPwrs),colWidsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contPwrs,reg_lnpwrs_wids,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(log(contPwrs),reg_lnpwrs_widsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Log(Intensity (mW))')

% Str_pwrs_wids = ['r^2 = ' char(string(r2valsWDTab{6,3}))  newline 'p = ' char(string(pvalsWDTab.Width{6}(numCoFs)))];
Str_pwrs_widsC = ['r^2 = ' char(string(r2valsWDTab{6,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{6}(numCoFs)))];
Str_pwrs_widsV = ['r^2 = ' num2str(r2_lnpwrs_widsV) newline 'p = ' num2str(mdl_lnpwrs_widsV.Coefficients.pValue(numCoFs))];
% text(5, 3.5, Str_pwrs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_widsC, 'fontsize', 8);
text(-1, 1, Str_pwrs_widsV, 'fontsize', 8);

% J. deps vs pwr
subplot('position',[xpos1(3) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colPwrs,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(log(colPwrs),colDepsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contPwrs,reg_lnpwrs_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(log(contPwrs),reg_lnpwrs_depsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Log(Intensity (mW))')

% Str_pwrs_deps = ['r^2 = ' char(string(r2valsWDTab{6,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{6}(numCoFs)))];
% Str_pwrs_depsC = ['r^2 = ' char(string(r2valsWDTab{6,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{6}(numCoFs)))];
Str_pwrs_depsV = ['r^2 = ' num2str(r2_lnpwrs_depsV) newline 'p = ' num2str(mdl_lnpwrs_depsV.Coefficients.pValue(numCoFs))];
% text(5, 1.8, Str_pwrs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 1, Str_pwrs_depsC, 'fontsize', 8);
text(-1, .75, Str_pwrs_depsV, 'fontsize', 8);

% make a legend for dashed and solid lines
subplot('position',[xpos1(3)-.04 .06 0 0])
plot(1,'linestyle','--','color',legColor,'linewidth',lw), hold on
plot(1,'linestyle','-','color',legColor,'linewidth',lw)
legend('All Values', 'Detected Values','box','off');
box off

[~,ht] = suplabel('Monkey','t',[-.37 .1 .84 .84]);
set(ht,'fontsize',11,'fontweight', 'normal');

[~,hb] = suplabel('B','t',[-.37 -.04 .84 .84]);
set(hb,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('C','t',[-.37 -.24 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('D','t',[-.37 -.44 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('E','t',[-.37 -.64 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hl] = suplabel('Left','t',[-.275 .1 .84 .84]);
set(hl,'fontsize',11,'fontweight','normal');

[~,hr] = suplabel('Right','t',[-.075 .1 .84 .84]);
set(hr,'fontsize',11,'fontweight','normal');

[~,ho1] = suplabel('OCTA','t',[-.3225 .05 .84 .84]);
set(ho1,'fontsize',11,'fontweight','normal');

[~,ho2] = suplabel('OCTA','t',[-.1225 .05 .84 .84]);
set(ho2,'fontsize',11,'fontweight','normal');

[~,hh1] = suplabel('Histology','t',[-.2225 .05 .84 .84]);
set(hh1,'fontsize',11,'fontweight','normal');

[~,hh2] = suplabel('Histology','t',[-.0225 .05 .84 .84]);
set(hh2,'fontsize',11,'fontweight','normal');

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1)+.05 ypos(6)-.02 3*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(end)+wh ypos(6)+(.5*wh) .011 3*(ypos(2)-ypos(3))];
caxis([min(colDeps(colDeps>0)) max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 4;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

% print -dpdf -painters Figure5V4

%% 40. Put Everything into one Figure (Composite of Composites) - Vessels Marked

figure('color','w','units','normalize','outerposition',[0 0 1 1])

OCTColor = [14 170 102]./255; %[18 226 136]./255;
% volColor = [248 117 117]./255;
volColor = [194 10 10]./255;
yelColor = [255 183 21]./255;
% yelColor = [255 210 21]./255;

% rb = redblue(2*nBin+1); rb = rb(round((2*nBin+1)/2):end,:); % red histo lesion
histColor = [linspace(255,volColor(1)*255,nBin+1);linspace(255,volColor(2)*255,nBin+1);linspace(255,volColor(3)*255,nBin+1)]'./255;
reDepths(reDepths == 0) = min(colDeps(colDeps>0)); % remove zero values for color scale
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;


yGrad = [linspace(255,yelColor(1)*255,nBin+1);linspace(255,yelColor(2)*255,nBin+1);linspace(255,yelColor(3)*255,nBin+1)]'./255;
% yGrad = [linspace(255,yelColor(1),nBin+1);linspace(255,yelColor(2),nBin+1);linspace(255,yelColor(3),nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(rePwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

% A. Circles figure covering the first half of figure
xpos = [.1 .3; .5 .7; .1 .3; .5 .7; .1 .3; .5 .7]./2;
ypos = [.75 .75 .55 .55 .35 .15];
wh = .2/2;

lims = [1 23];

cw = 0.25; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually

for iHemi = 1:size(reApDs,1)
    
    % Show OCTA first
    subplot('position',[xpos(iHemi,1) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5
        iAp = 8;
        
        % show OCTA diameter
        [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
        h = fill(xunit,yunit,OCTColor); hold on
        set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
    elseif iHemi == 6
        axis off
    else
        for iAp = 1:size(apDs,2)-1

            % show OCTA diameter
            [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,OCTColor); hold on
            set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
    
    % show histology second
    subplot('position',[xpos(iHemi,2) ypos(iHemi) wh wh], 'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        
        % show histology depth and width
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
        end
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
        
    else
        for iAp = 1:size(apDs,2)-1

            % show max histo lesion diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
            end
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
end

% add scale bar (1 cm)
plot(linspace(2.5,12.5,100),linspace(3,3,100),'linewidth',1.5,'color','k')

% now do the plots
mStyle = 'o';
mSize = 4;
xStyle = 'X';
xSize = 8;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

xpos1 = [.55 .7 .85]; ypos1 = [.6875 .375 .0625]+.06; wh1 = 0.1; h1 = 0.2;

lw = 1.5; % linewidth

% B. octa vs wid
subplot('position',[xpos1(1) ypos1(1) wh1 h1],'color','none')
for i = 1:length(colOCTADs)
    if colVes(i)
        plot(colOCTADs(i),colWids(i),xStyle,'markersize',xSize,'color',mixColor), hold on
%         plot(colOCTADsC(i),colWidsC(i),xStyle,'markersize',xSize,'color',mixColor), hold on
    else
        plot(colOCTADs(i),colWids(i),mStyle,'markersize',mSize,'color',mixColor), hold on
%         plot(colOCTADsC(i),colWidsC(i),mStyle,'markersize',mSize,'color',mixColor), hold on
    end
end    
plot(contOCTADs,reg_OCTADs_wids,'color',mixColor,'linewidth',lw,'linestyle','-'),
% plot(contOCTADsC,reg_OCTADsC_widsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 6.2])
ylabel('Lesion Width (mm)'), xlabel('OCTA Diameter (mm)')

Str_OCTADs_wids = ['r^2 = ' char(string(r2valsWDTab{4,3}))  newline 'p = ' char(string(pvalsWDTab.Width{4}(numCoFs)))];
% Str_OCTADsC_widsC = ['r^2 = ' char(string(r2valsWDTab{5,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{5}(numCoFs)))];
text(3, 1, Str_OCTADs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2, 5, Str_OCTADsC_widsC, 'fontsize', 8);

% C. dep vs OCTADs
subplot('position',[xpos1(2) ypos1(1) wh1 h1],'color','none')
for i = 1:length(colOCTADs)
    if colVes(i)
        plot(colOCTADs(i),colDeps(i),xStyle,'markersize',xSize,'color',mixColor), hold on
%         plot(colOCTADsC(i),colDepsC(i),xStyle,'markersize',xSize,'color',mixColor), hold on
    else
        plot(colOCTADs(i),colDeps(i),mStyle,'markersize',mSize,'color',mixColor), hold on
%         plot(colOCTADsC(i),colDepsC(i),mStyle,'markersize',mSize,'color',mixColor), hold on
    end
end
plot(contOCTADs,reg_OCTADs_deps,'color',mixColor,'linewidth',lw,'linestyle','-'),
% plot(contOCTADsC,reg_OCTADsC_depsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('OCTA Diameter (mm)')

Str_OCTADs_deps = ['r^2 = ' char(string(r2valsWDTab{4,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{4}(numCoFs)))];
% Str_OCTADsC_depsC = ['r^2 = ' char(string(r2valsWDTab{5,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{5}(numCoFs)))];
text(2.25, .75, Str_OCTADs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2.5, 0.75, Str_OCTADsC_depsC, 'fontsize', 8);

% D. wids vs deps
subplot('position',[xpos1(3) ypos1(1) wh1 h1],'color','none')
for i = 1:length(colWids)
    if colVes(i)
        plot(colWids(i),colDeps(i),xStyle,'markersize',xSize,'color',volColor), hold on
%         plot(colWidsC(i),colDepsC(i),xStyle,'markersize',xSize,'color',volColor), hold on
    else
        plot(colWids(i),colDeps(i),mStyle,'markersize',mSize,'color',volColor), hold on
%         plot(colWidsC(i),colDepsC(i),mStyle,'markersize',mSize,'color',volColor), hold on
    end
end
plot(contWids,reg_wids_deps,'color',volColor,'linewidth',lw,'linestyle','-'),
% plot(contWidsC,reg_widsC_depsC,'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Lesion Width (mm)')

Str_wids_deps = ['r^2 = ' char(string(r2valsWDTab{3,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{3}(numCoFs)))];
% Str_widsC_depsC = ['r^2 = ' char(string(r2valsWDTab{3,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{3}(numCoFs)))];
text(2, .75, Str_wids_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2.5, 0.75, Str_widsC_depsC, 'fontsize', 8);

% E. OCTADs vs apDs
subplot('position',[xpos1(1) ypos1(2) wh1 h1],'color','none')
for i = 1:length(colApDs)
    if colVes(i)
        plot(colApDs(i),colOCTADs(i),xStyle,'markersize',xSize,'color',OCTColor), hold on
%         plot(colApDs(i),colOCTADsC(i),xStyle,'markersize',xSize,'color',OCTColor), hold on
    else
        plot(colApDs(i),colOCTADs(i),mStyle,'markersize',mSize,'color',OCTColor), hold on
%         plot(colApDs(i),colOCTADsC(i),mStyle,'markersize',mSize,'color',OCTColor), hold on
    end
end
plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','-'),
% plot(contApDs,reg_apDs_OCTADsC, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
% Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
text(1.5, 1, Str_apDs_OCTADs, 'fontsize', 8);
% text(1.3, 1, Str_apDs_OCTADsC, 'fontsize', 8);

% F. wids vs apDs
subplot('position',[xpos1(2) ypos1(2) wh1 h1],'color','none')
for i = 1:length(colVes)
    if colVes(i)
        plot(colApDs(i),colWids(i),xStyle,'markersize',xSize,'color',volColor), hold on
%         plot(colApDs(i),colWidsC(i),xStyle,'markersize',xSize,'color',volColor), hold on
    else
        plot(colApDs(i),colWids(i),mStyle,'markersize',mSize,'color',volColor), hold on
%         plot(colApDs(i),colWidsC(i),mStyle,'markersize',mSize,'color',volColor), hold on
    end
end
plot(contApDs,reg_apDs_wids,'color',volColor,'linewidth',lw,'linestyle','-'),
% plot(contApDs,reg_apDs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_wids = ['r^2 = ' char(string(r2valsWDTab{1,3}))  newline 'p = ' char(string(pvalsWDTab.Width{1}(numCoFs)))];
% Str_apDs_widsC = ['r^2 = ' char(string(r2valsWDTab{1,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{1}(numCoFs)))];
text(1.5, 1, Str_apDs_wids, 'fontsize', 8);
% text(1.5, 1, Str_apDs_widsC, 'fontsize', 8);

% G. deps vs apDs
subplot('position',[xpos1(3) ypos1(2) wh1 h1],'color','none')
for i = 1:length(colVes)
    if colVes(i)
        plot(colApDs(i),colDeps(i),xStyle,'markersize',xSize,'color',volColor), hold on
%         plot(colApDs(i),colDepsC(i),xStyle,'markersize',xSize,'color',volColor), hold on
    else
        plot(colApDs(i),colDeps(i),mStyle,'markersize',mSize,'color',volColor), hold on
%         plot(colApDs(i),colDepsC(i),mStyle,'markersize',mSize,'color',volColor), hold on
    end
end    
plot(contApDs,reg_apDs_deps,'color',volColor,'linewidth',lw,'linestyle','-'),
% plot(contApDs,reg_apDs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Aperture Diameter (mm)')

Str_apDs_deps = ['r^2 = ' char(string(r2valsWDTab{1,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{1}(numCoFs)))];
% Str_apDs_depsC = ['r^2 = ' char(string(r2valsWDTab{1,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{1}(numCoFs)))];
text(1.25, .8, Str_apDs_deps, 'fontsize', 8);
% text(.25, 3.5, Str_apDs_depsC, 'fontsize', 8);

% H. OCTADs vs pwr
subplot('position',[xpos1(1) ypos1(3) wh1 h1],'color','none')
for i = 1:length(colVes)
    if colVes(i)
        semilogx(colPwrs(i),colOCTADs(i),xStyle,'markersize',xSize,'color',OCTColor), hold on        
%         semilogx(colPwrs(i),colOCTADsC(i),xStyle,'markersize',xSize,'color',OCTColor), hold on
    else
        semilogx(colPwrs(i),colOCTADs(i),mStyle,'markersize',mSize,'color',OCTColor), hold on        
%         semilogx(colPwrs(i),colOCTADsC(i),mStyle,'markersize',mSize,'color',OCTColor), hold on
    end
end
semilogx(contPwrs,reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','-'),
% semilogx(contPwrs,reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
% Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
text(.5, 1, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_OCTADsC, 'fontsize', 8);

% I. wids vs pwr
subplot('position',[xpos1(2) ypos1(3) wh1 h1],'color','none')
for i = 1:length(colVes)
    if colVes(i)
        semilogx(colPwrs(i),colWids(i),xStyle,'markersize',xSize,'color',volColor), hold on
%         semilogx(colPwrs(i),colWidsC(i),xStyle,'markersize',xSize,'color',volColor), hold on
    else
        semilogx(colPwrs(i),colWids(i),mStyle,'markersize',mSize,'color',volColor), hold on
%         semilogx(colPwrs(i),colWidsC(i),mStyle,'markersize',mSize,'color',volColor), hold on
    end
end
semilogx(contPwrs,reg_lnpwrs_wids,'color',volColor,'linewidth',lw,'linestyle','-'),
% semilogx(contPwrs,reg_lnpwrs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('Lesion Width (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_wids = ['r^2 = ' char(string(r2valsWDTab{6,3}))  newline 'p = ' char(string(pvalsWDTab.Width{6}(numCoFs)))];
% Str_pwrs_widsC = ['r^2 = ' char(string(r2valsWDTab{6,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{6}(numCoFs)))];
text(.5, 1, Str_pwrs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_widsC, 'fontsize', 8);

% J. deps vs pwr
subplot('position',[xpos1(3) ypos1(3) wh1 h1],'color','none')
for i = 1:length(colVes)
    if colVes(i)
        semilogx(colPwrs(i),colDeps(i),xStyle,'markersize',xSize,'color',volColor), hold on
%         semilogx(colPwrs(i),colDepsC(i),xStyle,'markersize',xSize,'color',volColor), hold on
    else
        semilogx(colPwrs(i),colDeps(i),mStyle,'markersize',mSize,'color',volColor), hold on
%         semilogx(colPwrs(i),colDepsC(i),mStyle,'markersize',mSize,'color',volColor), hold on
    end
end
semilogx(contPwrs,reg_lnpwrs_deps,'color',volColor,'linewidth',lw,'linestyle','-'),
% semilogx(contPwrs,reg_lnpwrs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylim([0 4])
ylabel('Lesion Depth (mm)'), xlabel('Light Intensity (mW)')

Str_pwrs_deps = ['r^2 = ' char(string(r2valsWDTab{6,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{6}(numCoFs)))];
% Str_pwrs_depsC = ['r^2 = ' char(string(r2valsWDTab{6,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{6}(numCoFs)))];
text(.5, .8, Str_pwrs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 1, Str_pwrs_depsC, 'fontsize', 8);

% make a legend for dashed and solid lines
subplot('position',[xpos1(3)-.04 .06 0 0])
plot(1,'x','color',legColor,'markersize',xSize), hold on
plot(1,'o','color',legColor,'markersize',mSize)
legend('Large Vessel Values', 'Non-large Vessel Values','box','off');
box off

[~,ht] = suplabel('Monkey','t',[-.37 .1 .84 .84]);
set(ht,'fontsize',11,'fontweight', 'normal');

[~,hb] = suplabel('B','t',[-.37 -.04 .84 .84]);
set(hb,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('C','t',[-.37 -.24 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('D','t',[-.37 -.44 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('E','t',[-.37 -.64 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hl] = suplabel('Left','t',[-.275 .1 .84 .84]);
set(hl,'fontsize',11,'fontweight','normal');

[~,hr] = suplabel('Right','t',[-.075 .1 .84 .84]);
set(hr,'fontsize',11,'fontweight','normal');

[~,ho1] = suplabel('OCTA','t',[-.3225 .05 .84 .84]);
set(ho1,'fontsize',11,'fontweight','normal');

[~,ho2] = suplabel('OCTA','t',[-.1225 .05 .84 .84]);
set(ho2,'fontsize',11,'fontweight','normal');

[~,hh1] = suplabel('Histology','t',[-.2225 .05 .84 .84]);
set(hh1,'fontsize',11,'fontweight','normal');

[~,hh2] = suplabel('Histology','t',[-.0225 .05 .84 .84]);
set(hh2,'fontsize',11,'fontweight','normal');

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1)+.05 ypos(6)-.02 3*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(end)+wh ypos(6)+(.5*wh) .011 3*(ypos(2)-ypos(3))];
caxis([min(colDeps(colDeps>0)) max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 4;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

% print -dpdf -painters FigureS1V1


%% 41. Put Everything into one Figure (Composite of Composites) - Non-vessels only

figure('color','w','units','normalize','outerposition',[0 0 1 1])

OCTColor = [14 170 102]./255; %[18 226 136]./255;
% volColor = [248 117 117]./255;
volColor = [194 10 10]./255;
yelColor = [255 183 21]./255;
% yelColor = [255 210 21]./255;

% rb = redblue(2*nBin+1); rb = rb(round((2*nBin+1)/2):end,:); % red histo lesion
histColor = [linspace(255,volColor(1)*255,nBin+1);linspace(255,volColor(2)*255,nBin+1);linspace(255,volColor(3)*255,nBin+1)]'./255;
reDepths(reDepths == 0) = min(colDeps(colDeps>0)); % remove zero values for color scale
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;


yGrad = [linspace(255,yelColor(1)*255,nBin+1);linspace(255,yelColor(2)*255,nBin+1);linspace(255,yelColor(3)*255,nBin+1)]'./255;
% yGrad = [linspace(255,yelColor(1),nBin+1);linspace(255,yelColor(2),nBin+1);linspace(255,yelColor(3),nBin+1)]'./255;
% divide power values into bins
logPwrs = log10(rePwrs);
yPwrs = round((logPwrs - min(logPwrs,[],'all')) * nBin / (max(logPwrs,[],'all') - min(logPwrs,[],'all')))+1;

% A. Circles figure covering the first half of figure
xpos = [.1 .3; .5 .7; .1 .3; .5 .7; .1 .3; .5 .7]./2;
ypos = [.75 .75 .55 .55 .35 .15];
wh = .2/2;

lims = [1 23];

cw = 0.25; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually

for iHemi = 1:size(reApDs,1)
    
    % Show OCTA first
    subplot('position',[xpos(iHemi,1) ypos(iHemi) wh wh],'color','none')
    
    if iHemi == 5
        iAp = 8;
        
        % show OCTA diameter
        [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
        h = fill(xunit,yunit,OCTColor); hold on
        set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
    elseif iHemi == 6
        axis off
    else
        for iAp = 1:size(apDs,2)-1

            % show OCTA diameter
            [xunit,yunit] = circle(x(iAp),y(iAp),reOCTAMaxDs(iHemi,iAp)/2);
            h = fill(xunit,yunit,OCTColor); hold on
            set(h,'edgecolor',OCTColor,'facealpha',1,'linewidth',cw);
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
    
    % show histology second
    subplot('position',[xpos(iHemi,2) ypos(iHemi) wh wh], 'color','none')
    
    if iHemi == 5 || iHemi == 6
        iAp = 8;
        
        % show histology depth and width
        if ~isnan(depths(iHemi,iAp))
            [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
            h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
            set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
        end
        
        % Show aperture diameter and intensity of illumination
        [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
        g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
        hold on
        daspect([1 1 1])
        xlim(lims)
        ylim(lims)
        set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
        axis off
        
    else
        for iAp = 1:size(apDs,2)-1

            % show max histo lesion diameter with color showing depth
            if ~isnan(depths(iHemi,iAp))
                [xunit,yunit] = circle(x(iAp),y(iAp),reWidths(iHemi,iAp)/2);
                h = fill(xunit,yunit,histColor(rbDeps(iHemi,iAp),:)); hold on
                set(h,'edgecolor',volColor,'facealpha',1,'linewidth',cw);
            end
            
            % Show diameter and intensity of illumination
            [xunit,yunit] = circle(x(iAp),y(iAp),reApDs(iHemi,iAp)/2);
            g = fill(xunit,yunit,yGrad(yPwrs(iHemi,iAp),:));
            hold on
            daspect([1 1 1])
            xlim(lims)
            ylim(lims)
            set(g,'edgecolor',yelColor,'facealpha',1,'linewidth',cw);
            axis off
        end
    end
end

% add scale bar (1 cm)
plot(linspace(2.5,12.5,100),linspace(3,3,100),'linewidth',1.5,'color','k')

% now do the plots
mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

xpos1 = [.55 .7 .85]; ypos1 = [.6875 .375 .0625]+.06; wh1 = 0.1; h1 = 0.2;

lw = 1.5; % linewidth

% B. octa vs wid
subplot('position',[xpos1(1) ypos1(1) wh1 h1],'color','none')
% plot(colOCTADs,colWids,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(colOCTADsC,colWidsC,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(colOCTADsV,colWidsV,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(contOCTADs,reg_OCTADs_wids,'color',mixColor,'linewidth',lw,'linestyle','--'),
% plot(contOCTADsC,reg_OCTADsC_widsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
plot(contOCTADsV,reg_OCTADsV_widsV, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 6.2])
ylabel('Histological Diameter (mm)'), xlabel('OCTA Diameter (mm)')

% Str_OCTADs_wids = ['r^2 = ' char(string(r2valsWDTab{4,3}))  newline 'p = ' char(string(pvalsWDTab.Width{4}(numCoFs)))];
% Str_OCTADsC_widsC = ['r^2 = ' char(string(r2valsWDTab{5,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{5}(numCoFs)))];
Str_OCTADsV_widsV = ['r^2 = ' num2str(r2_OCTADsV_widsV) newline 'p = ' num2str(mdl_OCTADsV_widsV.Coefficients.pValue(2))];
% text(2, 6, Str_OCTADs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2, 5, Str_OCTADsC_widsC, 'fontsize', 8);
text(2, 5, Str_OCTADsV_widsV, 'fontsize', 8);

% C. dep vs OCTADs
subplot('position',[xpos1(2) ypos1(1) wh1 h1],'color','none')
% plot(colOCTADs,colDeps,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(colOCTADsC,colDepsC,mStyle,'markersize',mSize,'color',mixColor), hold on
plot(colOCTADsV,colDepsV,mStyle,'markersize',mSize,'color',mixColor), hold on
% plot(contOCTADs,reg_OCTADs_deps,'color',mixColor,'linewidth',lw,'linestyle','--'),
% plot(contOCTADsC,reg_OCTADsC_depsC, 'color',mixColor,'linewidth',lw,'linestyle','-')
plot(contOCTADsV,reg_OCTADsV_depsV, 'color',mixColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Histological Depth (mm)'), xlabel('OCTA Diameter (mm)')

% Str_OCTADs_deps = ['r^2 = ' char(string(r2valsWDTab{4,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{4}(numCoFs)))];
% Str_OCTADsC_depsC = ['r^2 = ' char(string(r2valsWDTab{5,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{5}(numCoFs)))];
Str_OCTADsV_depsV = ['r^2 = ' num2str(r2_OCTADsV_depsV) newline 'p = ' num2str(mdl_OCTADsV_depsV.Coefficients.pValue(2))];
% text(4, 4.2, Str_OCTADs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2.5, 0.75, Str_OCTADsC_depsC, 'fontsize', 8);
text(2.5, 0.75, Str_OCTADsV_depsV, 'fontsize', 8);

% D. wids vs deps
subplot('position',[xpos1(3) ypos1(1) wh1 h1],'color','none')
% plot(colWids,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colWidsC,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colWidsV,colDepsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contWids,reg_wids_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contWidsC,reg_widsC_depsC,'color',volColor,'linewidth',lw,'linestyle','-')
plot(contWidsV,reg_widsV_depsV,'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
ylim([0 4])
ylabel('Histological Depth (mm)'), xlabel('Histological Diameter (mm)')

% Str_wids_deps = ['r^2 = ' char(string(r2valsWDTab{3,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{3}(numCoFs)))];
% Str_widsC_depsC = ['r^2 = ' char(string(r2valsWDTab{3,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{3}(numCoFs)))];
Str_widsC_depsV = ['r^2 = ' num2str(r2_widsV_depsV) newline 'p = ' num2str(mdl_widsV_depsV.Coefficients.pValue(2))];
% text(1.2, .5, Str_wids_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(2.5, 0.75, Str_widsC_depsC, 'fontsize', 8);
text(2.5, 0.75, Str_widsC_depsV, 'fontsize', 8);

% E. OCTADs vs apDs
subplot('position',[xpos1(1) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(colApDs,colOCTADsC,mStyle,'markersize',mSize,'color',OCTColor), hold on
plot(colApDs,colOCTADsV,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(contApDs,reg_apDs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_OCTADs = ['r^2 = ' char(string(r2valsTab{1,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{1}(numCoFs)))];
% Str_apDs_OCTADsC = ['r^2 = ' char(string(r2valsTab{1,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{1}(numCoFs)))];
Str_apDs_OCTADsV = ['r^2 = ' num2str(r2_apDs_OCTADsV) newline 'p = ' num2str(mdl_apDs_OCTADsV.Coefficients.pValue(2))];
% text(1.5, 1.3, Str_apDs_OCTADs, 'fontsize', 8);
% text(1.3, 1, Str_apDs_OCTADsC, 'fontsize', 8);
text(1.3, 1, Str_apDs_OCTADsV, 'fontsize', 8);

% F. wids vs apDs
subplot('position',[xpos1(2) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colApDs,colWidsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colApDs,colWidsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contApDs,reg_apDs_wids,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_widsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylabel('Histological Diameter (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_wids = ['r^2 = ' char(string(r2valsWDTab{1,3}))  newline 'p = ' char(string(pvalsWDTab.Width{1}(numCoFs)))];
% Str_apDs_widsC = ['r^2 = ' char(string(r2valsWDTab{1,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{1}(numCoFs)))];
Str_apDs_widsV = ['r^2 = ' num2str(r2_apDs_widsV) newline 'p = ' num2str(mdl_apDs_widsV.Coefficients.pValue(numCoFs))];
% text(1.5, 1.75, Str_apDs_wids, 'fontsize', 8);
% text(1.5, 1, Str_apDs_widsC, 'fontsize', 8);
text(1.3, 1, Str_apDs_widsV, 'fontsize', 8);

% G. deps vs apDs
subplot('position',[xpos1(3) ypos1(2) wh1 h1],'color','none')
% plot(colApDs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colApDs,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
plot(colApDs,colDepsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contApDs,reg_apDs_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contApDs,reg_apDs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
plot(contApDs,reg_apDs_depsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
xlim([0 2.5])
ylim([0 4])
ylabel('Histological Depth (mm)'), xlabel('Aperture Diameter (mm)')

% Str_apDs_deps = ['r^2 = ' char(string(r2valsWDTab{1,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{1}(numCoFs)))];
% Str_apDs_depsC = ['r^2 = ' char(string(r2valsWDTab{1,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{1}(numCoFs)))];
Str_apDs_depsV = ['r^2 = ' num2str(r2_apDs_depsV) newline 'p = ' num2str(mdl_apDs_depsV.Coefficients.pValue(numCoFs))];
% text(1.25, 1, Str_apDs_deps, 'fontsize', 8);
% text(.25, 3.5, Str_apDs_depsC, 'fontsize', 8);
text(1.3, 1, Str_apDs_depsV, 'fontsize', 8);

% H. OCTADs vs pwr
subplot('position',[xpos1(1) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colOCTADs,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(colPwrs,colOCTADsC,mStyle,'markersize',mSize,'color',OCTColor), hold on
semilogx(colPwrs,colOCTADsV,mStyle,'markersize',mSize,'color',OCTColor), hold on
% plot(contPwrs,reg_lnpwrs_OCTADs,'color',OCTColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_OCTADsC, 'color',OCTColor,'linewidth',lw,'linestyle','-')
semilogx(contPwrs,reg_lnpwrs_OCTADsV, 'color',OCTColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('OCTA Diameter (mm)'), xlabel('Light Intensity (mW)')

% Str_pwrs_OCTADs = ['r^2 = ' char(string(r2valsTab{2,1}))  newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter{2}(numCoFs)))];
% Str_pwrs_OCTADsC = ['r^2 = ' char(string(r2valsTab{2,2})) newline 'p = ' char(string(pvalsTab.Max_OCTA_Diameter_Detected{2}(numCoFs)))];
Str_pwrs_OCTADsV = ['r^2 = ' num2str(r2_pwrs_OCTADsV) newline 'p = ' num2str(mdl_pwrs_OCTADsV.Coefficients.pValue(numCoFs))];
% text(5, 3.2, Str_pwrs_OCTADs, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_OCTADsC, 'fontsize', 8);
text(1, 1, Str_pwrs_OCTADsV, 'fontsize', 8);

% I. wids vs pwr
subplot('position',[xpos1(2) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colWids,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colPwrs,colWidsC,mStyle,'markersize',mSize,'color',volColor), hold on
semilogx(colPwrs,colWidsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contPwrs,reg_lnpwrs_wids,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_widsC, 'color',volColor,'linewidth',lw,'linestyle','-')
semilogx(contPwrs,reg_lnpwrs_widsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylabel('Histological Diameter (mm)'), xlabel('Light Intensity (mW)')

% Str_pwrs_wids = ['r^2 = ' char(string(r2valsWDTab{6,3}))  newline 'p = ' char(string(pvalsWDTab.Width{6}(numCoFs)))];
Str_pwrs_widsC = ['r^2 = ' char(string(r2valsWDTab{6,4})) newline 'p = ' char(string(pvalsWDTab.Width_Detected{6}(numCoFs)))];
Str_pwrs_widsV = ['r^2 = ' num2str(r2_lnpwrs_widsV) newline 'p = ' num2str(mdl_lnpwrs_widsV.Coefficients.pValue(numCoFs))];
% text(5, 3.5, Str_pwrs_wids, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 2, Str_pwrs_widsC, 'fontsize', 8);
text(1, 1, Str_pwrs_widsV, 'fontsize', 8);

% J. deps vs pwr
subplot('position',[xpos1(3) ypos1(3) wh1 h1],'color','none')
% plot(colPwrs,colDeps,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(colPwrs,colDepsC,mStyle,'markersize',mSize,'color',volColor), hold on
semilogx(colPwrs,colDepsV,mStyle,'markersize',mSize,'color',volColor), hold on
% plot(contPwrs,reg_lnpwrs_deps,'color',volColor,'linewidth',lw,'linestyle','--'),
% plot(contPwrs,reg_lnpwrs_depsC, 'color',volColor,'linewidth',lw,'linestyle','-')
semilogx(contPwrs,reg_lnpwrs_depsV, 'color',volColor,'linewidth',lw,'linestyle','-')
legend off, box off
% xlim([0 2.5])
ylim([0 4])
ylabel('Histological Depth (mm)'), xlabel('Light Intensity (mW)')

% Str_pwrs_deps = ['r^2 = ' char(string(r2valsWDTab{6,1}))  newline 'p = ' char(string(pvalsWDTab.Depth{6}(numCoFs)))];
% Str_pwrs_depsC = ['r^2 = ' char(string(r2valsWDTab{6,2})) newline 'p = ' char(string(pvalsWDTab.Depth_Detected{6}(numCoFs)))];
Str_pwrs_depsV = ['r^2 = ' num2str(r2_lnpwrs_depsV) newline 'p = ' num2str(mdl_lnpwrs_depsV.Coefficients.pValue(numCoFs))];
% text(5, 1.8, Str_pwrs_deps, 'fontsize', 8);% set(T, 'fontsize', 8);
% text(5, 1, Str_pwrs_depsC, 'fontsize', 8);
text(1, .75, Str_pwrs_depsV, 'fontsize', 8);

% make a legend for dashed and solid lines
% subplot('position',[xpos1(3)-.04 .06 0 0])
% plot(1,'linestyle','--','color',legColor,'linewidth',lw), hold on
% plot(1,'linestyle','-','color',legColor,'linewidth',lw)
% legend('All Values', 'Detected Values','box','off');
% box off

[~,ht] = suplabel('Monkey','t',[-.37 .1 .84 .84]);
set(ht,'fontsize',11,'fontweight', 'normal');

[~,hb] = suplabel('B','t',[-.37 -.04 .84 .84]);
set(hb,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('C','t',[-.37 -.24 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('D','t',[-.37 -.44 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hc] = suplabel('E','t',[-.37 -.64 .84 .84]);
set(hc,'fontsize',11,'fontweight', 'normal');

[~,hl] = suplabel('Left','t',[-.275 .1 .84 .84]);
set(hl,'fontsize',11,'fontweight','normal');

[~,hr] = suplabel('Right','t',[-.075 .1 .84 .84]);
set(hr,'fontsize',11,'fontweight','normal');

[~,ho1] = suplabel('OCTA','t',[-.3225 .05 .84 .84]);
set(ho1,'fontsize',11,'fontweight','normal');

[~,ho2] = suplabel('OCTA','t',[-.1225 .05 .84 .84]);
set(ho2,'fontsize',11,'fontweight','normal');

[~,hh1] = suplabel('Histology','t',[-.2225 .05 .84 .84]);
set(hh1,'fontsize',11,'fontweight','normal');

[~,hh2] = suplabel('Histology','t',[-.0225 .05 .84 .84]);
set(hh2,'fontsize',11,'fontweight','normal');

% add color bar for light intensity on a log scale
yc = colorbar('horiz'); colormap(yGrad);
yc.TickDirection = 'in';
yc.Location = 'manual';
yc.Position = [xpos(1)+.05 ypos(6)-.02 3*wh .02];
set(yc, 'fontsize', 11);
caxis([min(logPwrs,[],'all') max(logPwrs,[],'all')]);
colorTitleHandle = get(yc,'xlabel'); titleString = 'Light Intensity (mW)';
set(colorTitleHandle,'String',titleString)
tickLabs = [1 10];
ticks = log10(tickLabs);
yc.Ticks = ticks; yc.TickLabels = tickLabs;

% add color bar for lesion depths
ax = axes; axis off
rbc = colorbar; colormap(ax,histColor);
rbc.TickDirection = 'in';
rbc.Location = 'manual';
rbc.Position = [xpos(end)+wh ypos(6)+(.5*wh) .011 3*(ypos(2)-ypos(3))];
caxis([min(colDeps(colDeps>0)) max(depths,[],'all')]);
colorTitleHandle = get(rbc,'ylabel'); titleString = 'Lesion Depth (mm)';
set(colorTitleHandle,'String',titleString,'rotation',-90);
pos = get(colorTitleHandle,'position'); pos(1) = pos(1) + 4;
set(colorTitleHandle,'position',pos);
set(rbc,'fontsize',11);

print -dpdf -painters Figure5V5

%% 42. Compare Simulation Depths with Histological Depth

mdl_simDeps_deps = fitlm(colSimDeps,colDeps,'intercept',false);
coF_simDeps_deps = mdl_simDeps_deps.Coefficients.Estimate;
reg_simDeps_deps = (coF_simDeps_deps(1) * contSimDeps);
est_simDeps_deps = (coF_simDeps_deps(1) * colSimDeps);
r2_simDeps_deps1 = calculateR2(colDeps, est_simDeps_deps);
r2_simDeps_deps2 = calculateR2(colDeps, colSimDeps);

mdl_simDeps_depsC = fitlm(colSimDepsC,colDepsC,'intercept',false);
coF_simDeps_depsC = mdl_simDeps_depsC.Coefficients.Estimate;
reg_simDeps_depsC = (coF_simDeps_depsC(1) * contSimDepsC);
est_simDeps_depsC = (coF_simDeps_depsC(1) * colSimDepsC);
r2_simDeps_depsC1 = calculateR2(colDepsC, est_simDeps_depsC);
r2_simDeps_depsC2 = calculateR2(colDepsC, colSimDepsC);

mdl_simDeps_depsV = fitlm(colSimDepsV,colDepsV,'intercept',false);
coF_simDeps_depsV = mdl_simDeps_depsV.Coefficients.Estimate;
reg_simDeps_depsV = (coF_simDeps_depsV(1) * contSimDepsV);
est_simDeps_depsV = (coF_simDeps_depsV(1) * colSimDepsV);
r2_simDeps_depsV1 = calculateR2(colDepsV, est_simDeps_depsV);
r2_simDeps_depsV2 = calculateR2(colDepsV, colSimDepsV);

mdl_simDeps_depsVC = fitlm(colSimDepsVC,colDepsVC,'intercept',false);
coF_simDeps_depsVC = mdl_simDeps_depsVC.Coefficients.Estimate;
reg_simDeps_depsVC = (coF_simDeps_depsVC(1) * contSimDepsVC);
est_simDeps_depsVC = (coF_simDeps_depsVC(1) * colSimDepsVC);
r2_simDeps_depsVC1 = calculateR2(colDepsVC, est_simDeps_depsVC);
r2_simDeps_depsVC2 = calculateR2(colDepsVC, colSimDepsVC);

figure('color','w')
subplot(1,2,1)
plot(colSimDepsVC,colDepsVC,'o'), hold on
plot(contSimDepsVC,contSimDepsVC,'--'),
plot(contSimDepsVC,reg_simDeps_depsVC)
box off
legend('Depths',['y = x; r^2 = ' num2str(r2_simDeps_depsVC2)],['y = ' num2str(coF_simDeps_depsVC) 'x; r^2 = ' num2str(r2_simDeps_depsC1)],'Location','southeast')
xlim([0 3]), ylim([0 3])
xlabel('Simulated Lesion Depth (mm)')
ylabel('Histological Lesion Depth (mm)')
title('Simulation vs Histological Lesion Depth')

subplot(1,2,2)
plot(colSimDepsV,colDepsV,'o'), hold on
plot(contSimDepsV,contSimDepsV,'--'),
plot(contSimDepsV,reg_simDeps_depsV)
box off
legend('Depths',['y = x; r^2 = ' num2str(r2_simDeps_depsV2)],['y = ' num2str(coF_simDeps_depsV) 'x; r^2 = ' num2str(r2_simDeps_depsV1)],'Location','southeast')
xlim([0 3]), ylim([0 3])
xlabel('Simulated Lesion Depth (mm)')
ylabel('Histological Lesion Depth (mm)')
title('Simulation vs Histological Lesion Depth (with zeros)')

%% 43. Compare Simulation Widths with Histological Width

mdl_simWids_wids = fitlm(colSimWids,colWids,'intercept',false);
coF_simWids_wids = mdl_simWids_wids.Coefficients.Estimate;
reg_simWids_wids = (coF_simWids_wids(1) * contSimWids);
est_simWids_wids = (coF_simWids_wids(1) * colSimWids);
r2_simWids_wids1 = calculateR2(colWids, est_simWids_wids);
r2_simWids_wids2 = calculateR2(colWids, colSimWids);

mdl_simWids_widsC = fitlm(colSimWidsC,colWidsC,'intercept',false);
coF_simWids_widsC = mdl_simWids_widsC.Coefficients.Estimate;
reg_simWids_widsC = (coF_simWids_widsC(1) * contSimWidsC);
est_simWids_widsC = (coF_simWids_widsC(1) * colSimWidsC);
r2_simWids_widsC1 = calculateR2(colWidsC, est_simWids_widsC);
r2_simWids_widsC2 = calculateR2(colWidsC, colSimWidsC);

mdl_simWids_widsV = fitlm(colSimWidsV,colWidsV,'intercept',false);
coF_simWids_widsV = mdl_simWids_widsV.Coefficients.Estimate;
reg_simWids_widsV = (coF_simWids_widsV(1) * contSimWidsV);
est_simWids_widsV = (coF_simWids_widsV(1) * colSimWidsV);
r2_simWids_widsV1 = calculateR2(colWidsV, est_simWids_widsV);
r2_simWids_widsV2 = calculateR2(colWidsV, colSimWidsV);

mdl_simWids_widsVC = fitlm(colSimWidsVC,colWidsVC,'intercept',false);
coF_simWids_widsVC = mdl_simWids_widsVC.Coefficients.Estimate;
reg_simWids_widsVC = (coF_simWids_widsVC(1) * contSimWidsVC);
est_simWids_widsVC = (coF_simWids_widsVC(1) * colSimWidsVC);
r2_simWids_widsVC1 = calculateR2(colWidsVC, est_simWids_widsVC);
r2_simWids_widsVC2 = calculateR2(colWidsVC, colSimWidsVC);

figure('color','w')
subplot(1,2,1)
plot(colSimWidsVC,colWidsVC,'o'), hold on
plot(contSimWidsVC,contSimWidsVC,'--'),
plot(contSimWidsVC,reg_simWids_widsVC)
box off
legend('Widths',['y = x; r^2 = ' num2str(r2_simWids_widsVC2)],['y = ' num2str(coF_simWids_widsVC) 'x; r^2 = ' num2str(r2_simWids_widsVC1)],'Location','southeast')
xlim([0 7]), ylim([0 7])
xlabel('Simulated Lesion Width (mm)')
ylabel('Histological Lesion Width (mm)')
title('Simulation vs Histological Lesion Width')

subplot(1,2,2)
plot(colSimWidsV,colWidsV,'o'), hold on
plot(contSimWidsV,contSimWidsV,'--'),
plot(contSimWidsV,reg_simWids_widsV)
box off
legend('Widths',['y = x; r^2 = ' num2str(r2_simWids_widsV2)],['y = ' num2str(coF_simWids_widsV) 'x; r^2 = ' num2str(r2_simWids_widsV1)],'Location','southeast')
xlim([0 7]), ylim([0 7])
xlabel('Simulated Lesion Width (mm)')
ylabel('Histological Lesion Width (mm)')
title('Simulation vs Histological Lesion Width (with zeros)')

%% 44. Compare Simulation Widths with OCTA Width

mdl_simWids_OCTADs = fitlm(colSimWids,colOCTADs,'intercept',false);
coF_simWids_OCTADs = mdl_simWids_OCTADs.Coefficients.Estimate;
reg_simWids_OCTADs = (coF_simWids_OCTADs(1) * contSimWids);
est_simWids_OCTADs = (coF_simWids_OCTADs(1) * colSimWids);
r2_simWids_OCTADs1 = calculateR2(colOCTADs, est_simWids_OCTADs);
r2_simWids_OCTADs2 = calculateR2(colOCTADs, colSimWids);

mdl_simWids_OCTADsC = fitlm(colSimWidsC,colOCTADsC,'intercept',false);
coF_simWids_OCTADsC = mdl_simWids_OCTADsC.Coefficients.Estimate;
reg_simWids_OCTADsC = (coF_simWids_OCTADsC(1) * contSimWidsC);
est_simWids_OCTADsC = (coF_simWids_OCTADsC(1) * colSimWidsC);
r2_simWids_OCTADsC1 = calculateR2(colOCTADsC, est_simWids_OCTADsC);
r2_simWids_OCTADsC2 = calculateR2(colOCTADsC, colSimWids);

mdl_simWids_OCTADsV = fitlm(colSimWidsV,colOCTADsV,'intercept',false);
coF_simWids_OCTADsV = mdl_simWids_OCTADsV.Coefficients.Estimate;
reg_simWids_OCTADsV = (coF_simWids_OCTADsV(1) * contSimWidsV);
est_simWids_OCTADsV = (coF_simWids_OCTADsV(1) * colSimWidsV);
r2_simWids_OCTADsV1 = calculateR2(colOCTADsV, est_simWids_OCTADsV);
r2_simWids_OCTADsV2 = calculateR2(colOCTADsV, colSimWidsV);

mdl_simWids_OCTADsVC = fitlm(colSimWidsVC,colOCTADsVC,'intercept',false);
coF_simWids_OCTADsVC = mdl_simWids_OCTADsVC.Coefficients.Estimate;
reg_simWids_OCTADsVC = (coF_simWids_OCTADsVC(1) * contSimWidsVC);
est_simWids_OCTADsVC = (coF_simWids_OCTADsVC(1) * colSimWidsVC);
r2_simWids_OCTADsVC1 = calculateR2(colOCTADsVC, est_simWids_OCTADsVC);
r2_simWids_OCTADsVC2 = calculateR2(colOCTADsVC, colSimWidsVC);

figure('color','w')
subplot(1,2,1)
plot(colSimWidsVC,colOCTADsVC,'o'), hold on
plot(contSimWidsVC,contSimWidsVC,'--'),
plot(contSimWidsVC,reg_simWids_OCTADsVC)
box off
legend('Widths',['y = x; r^2 = ' num2str(r2_simWids_OCTADsVC2)],['y = ' num2str(coF_simWids_OCTADsVC) 'x; r^2 = ' num2str(r2_simWids_OCTADsVC1)],'Location','southeast')
xlim([0 7]), ylim([0 7])
xlabel('Simulated Lesion Width (mm)')
ylabel('OCTA Lesion Width (mm)')
title('Simulation vs OCTA Lesion Width')

subplot(1,2,2)
plot(colSimWidsV,colOCTADsV,'o'), hold on
plot(contSimWidsV,contSimWidsV,'--'),
plot(contSimWidsV,reg_simWids_OCTADsV)
box off
legend('Widths',['y = x; r^2 = ' num2str(r2_simWids_OCTADsV2)],['y = ' num2str(coF_simWids_OCTADsV) 'x; r^2 = ' num2str(r2_simWids_OCTADsV1)],'Location','southeast')
xlim([0 7]), ylim([0 7])
xlabel('Simulated Lesion Width (mm)')
ylabel('OCTA Lesion Width (mm)')
title('Simulation vs OCTA Lesion Width (with zeros)')

%% 45. Plot the Residuals vs Fitted Values for Simulation vs Histo

mdl_nsSimDeps_depsV = fitlm(colNSSimDepsV,colDepsV,'intercept',false);
coF_nsSimDeps_depsV = mdl_nsSimDeps_depsV.Coefficients.Estimate;
reg_nsSimDeps_depsV = (coF_nsSimDeps_depsV(1) * contNSSimDepsV);
est_nsSimDeps_depsV = (coF_nsSimDeps_depsV(1) * colNSSimDepsV);
r2_nsSimDeps_depsV1 = calculateR2(colDepsV, est_nsSimDeps_depsV);
r2_nsSimDeps_depsV2 = calculateR2(colDepsV, colNSSimDepsV);

mdl_SQRTnsSimDeps_depsV = fitlm(sqrt(colNSSimDepsV),colDepsV,'intercept',false);
coF_SQRTnsSimDeps_depsV = mdl_SQRTnsSimDeps_depsV.Coefficients.Estimate;
reg_SQRTnsSimDeps_depsV = (coF_SQRTnsSimDeps_depsV(1) * sqrt(contNSSimDepsV));
est_SQRTnsSimDeps_depsV = (coF_SQRTnsSimDeps_depsV(1) * sqrt(colNSSimDepsV));
r2_SQRTnsSimDeps_depsV1 = calculateR2(colDepsV, est_SQRTnsSimDeps_depsV);
r2_SQRTnsSimDeps_depsV2 = calculateR2(colDepsV, sqrt(colNSSimDepsV));

mdl_nsSimWids_widsV = fitlm(colNSSimWidsV,colWidsV,'intercept',false);
coF_nsSimWids_widsV = mdl_nsSimWids_widsV.Coefficients.Estimate;
reg_nsSimWids_widsV = (coF_nsSimWids_widsV(1) * contNSSimWidsV);
est_nsSimWids_widsV = (coF_nsSimWids_widsV(1) * colNSSimWidsV);
r2_nsSimWids_widsV1 = calculateR2(colWidsV, est_nsSimWids_widsV);
r2_nsSimWids_widsV2 = calculateR2(colWidsV, colNSSimWidsV);

OCTColor = [14 170 102]./255; %[18 226 136]./255;
volColor = [248 117 117]./255;

cw = 1.; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually
mStyle = 'o';
mSize = 4;

figure('color','w','units','normalize','outerposition',[0 0 .65 1])
xPos = [0.05 0.37 0.69]; yPos = .3; w = .28; h = .28; %w = .2; h = .2;

subplot('position',[xPos(1) yPos w h])
plot(colNSSimDepsV,colDepsV,mStyle,'markersize',mSize,'linewidth',cw,'color',volColor), hold on
plot(contNSSimDepsV,reg_nsSimDeps_depsV,'linewidth',cw,'color',volColor)
box off
xlim([0 4]); ylim([0 5]);
xlabel('Simulated Depth (mm)'), ylabel('Histological Depth (mm)')
str_nsSimDeps_deps = ['r^2 = ' num2str(r2_nsSimDeps_depsV1)];
text(2.25, 1, str_nsSimDeps_deps, 'fontsize', 10);

subplot('position',[xPos(2) yPos w h])
plot(colNSSimDepsV,colDepsV,mStyle,'markersize',mSize,'linewidth',cw,'color',volColor), hold on
plot(contNSSimDepsV,reg_SQRTnsSimDeps_depsV,'linewidth',cw,'color',volColor)
box off
xlim([0 4]); ylim([0 5]);
xlabel('Simulated Depth (mm)'), ylabel('Histological Depth (mm)')
str_nsSimDeps_deps = ['r^2 = ' num2str(r2_SQRTnsSimDeps_depsV1)];
text(2.25, 1, str_nsSimDeps_deps, 'fontsize', 10);

subplot('position',[xPos(3) yPos w h])
plot(mdl_nsSimDeps_depsV.Fitted,mdl_nsSimDeps_depsV.Residuals.Raw,mStyle,'markersize',mSize,'linewidth',cw,'color',volColor), hold on
% plot(zeros(size(reg_nsSimDeps_depsV))','linewidth',cw,'--','color',volColor)
box off
% xlim([0 4]); ylim([0 5]);
yline(0,'--','linewidth',cw)
xlabel('Fitted Depth (mm)'), ylabel('Residuals (mm)')%, title('Residuals of Modeled vs Histological Depths')
% str_nsSimDeps_deps = ['r^2 = ' num2str(r2_SQRTnsSimDeps_depsV1)];
% text(2.25, 1, str_nsSimDeps_deps, 'fontsize', 10);

print -dpdf -painters FigureS2V2
%% 45. Make a Plot Comparing Simulation and Measured Lesion Sizes for Sim Figure

OCTColor = [14 170 102]./255; %[18 226 136]./255;
volColor = [248 117 117]./255;

cw = 1.; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually
mStyle = 'o';
mSize = 4;

figure('color','w','units','normalize','outerposition',[0 0 .65 1])
xPos = [0.05 0.37 0.69]; yPos = .3; w = .28; h = .28; %w = .2; h = .2;

subplot('position',[xPos(1) yPos w h])
plot(colSimWidsV,colWidsV,mStyle,'markersize',mSize,'linewidth',cw,'color',volColor), hold on
plot(contSimWidsV,contSimWidsV,'linewidth',cw,'color',volColor)
box off
xlim([0 4]); ylim([0 5]);
xlabel('Model Diameter (mm)'), ylabel('Histological Diameter (mm)')
str_simWids_wids = ['r^2 = ' num2str(r2_simWids_widsV2)];
text(2.25, 1, str_simWids_wids, 'fontsize', 10);

subplot('position',[xPos(2) yPos w h])
plot(colSimDepsV,colDepsV,mStyle,'markersize',mSize,'linewidth',cw,'color',volColor), hold on
plot(contSimDepsV,contSimDepsV,'linewidth',cw,'color',volColor)
box off
xlim([0 4]); ylim([0 4]);
xlabel('Model Depth (mm)'), ylabel('Histological Depth (mm)')
str_simDeps_deps = ['r^2 = ' num2str(r2_simDeps_depsV2)];
text(2.25, 1, str_simDeps_deps, 'fontsize', 10);

subplot('position',[xPos(3) yPos w h])
plot(colSimWidsV,colOCTADsV,mStyle,'markersize',mSize,'linewidth',cw,'color',OCTColor), hold on
plot(contSimWidsV,contSimWidsV,'linewidth',cw,'color',OCTColor)
box off
xlim([0 4]); ylim([0 4]);
xlabel('Model Diameter (mm)'), ylabel('OCTA Diameter (mm)')
str_simWids_OCTADs = ['r^2 = ' num2str(r2_simWids_OCTADsV2)];
text(2.25, 1, str_simWids_OCTADs, 'fontsize', 10);

print -dpdf -painters SimComparisonsFig
%% 46. Model Log(Intensity) and Aperture Diameters on OCTA Diameters 

cepts = true;
tbl = [log(colPwrs) colApDs];


mdl_apDs_lnpwrs_OCTADs = fitlm(tbl, colOCTADs, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (mm)'},'intercept',cepts)
coF_apDs_lnpwrs_OCTADs = mdl_apDs_lnpwrs_OCTADs.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_OCTADs = coF_apDs_lnpwrs_OCTADs(1) + (coF_apDs_lnpwrs_OCTADs(2) * log(contPwrs)) + (coF_apDs_lnpwrs_OCTADs(3) * contApDs);
    est_apDs_lnpwrs_OCTADs = coF_apDs_lnpwrs_OCTADs(1) + (coF_apDs_lnpwrs_OCTADs(2) * log(colPwrs)) + (coF_apDs_lnpwrs_OCTADs(3) * colApDs);
else
    reg_apDs_lnpwrs_OCTADs = (coF_apDs_lnpwrs_OCTADs(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_OCTADs(2) * contApDs);
    est_apDs_lnpwrs_OCTADs = (coF_apDs_lnpwrs_OCTADs(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_OCTADs(2) * colApDs);
end 
r2_apDs_lnpwrs_OCTADs = calculateR2(colOCTADs,est_apDs_lnpwrs_OCTADs);

mdl_apDs_lnpwrs_OCTADsC = fitlm(tbl, colOCTADsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_OCTADsC = mdl_apDs_lnpwrs_OCTADsC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_OCTADsC = coF_apDs_lnpwrs_OCTADsC(1) + (coF_apDs_lnpwrs_OCTADsC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_OCTADsC(3) * contApDs);
    est_apDs_lnpwrs_OCTADsC = coF_apDs_lnpwrs_OCTADsC(1) + (coF_apDs_lnpwrs_OCTADsC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_OCTADsC(3) * colApDs);
else
    reg_apDs_lnpwrs_OCTADsC = (coF_apDs_lnpwrs_OCTADsC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_OCTADsC(2) * contApDs);
    est_apDs_lnpwrs_OCTADsC = (coF_apDs_lnpwrs_OCTADsC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_OCTADsC(2) * colApDs);
end
r2_apDs_lnpwrs_OCTADsC = calculateR2(colOCTADsC,est_apDs_lnpwrs_OCTADsC);

mdl_apDs_lnpwrs_OCTADsV = fitlm(tbl, colOCTADsV, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Non-Vessels) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_OCTADsV = mdl_apDs_lnpwrs_OCTADsV.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_OCTADsV = coF_apDs_lnpwrs_OCTADsV(1) + (coF_apDs_lnpwrs_OCTADsV(2) * log(contPwrs)) + (coF_apDs_lnpwrs_OCTADsV(3) * contApDs);
    est_apDs_lnpwrs_OCTADsV = coF_apDs_lnpwrs_OCTADsV(1) + (coF_apDs_lnpwrs_OCTADsV(2) * log(colPwrs)) + (coF_apDs_lnpwrs_OCTADsV(3) * colApDs);
else
    reg_apDs_lnpwrs_OCTADsV = (coF_apDs_lnpwrs_OCTADsV(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_OCTADsV(2) * contApDs);
    est_apDs_lnpwrs_OCTADsV = (coF_apDs_lnpwrs_OCTADsV(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_OCTADsV(2) * colApDs);
end
r2_apDs_lnpwrs_OCTADsV = calculateR2(colOCTADsV,mdl_apDs_lnpwrs_OCTADsV.Fitted);

mdl_apDs_lnpwrs_OCTADsVC = fitlm(tbl, colOCTADsVC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Max. OCTA Diameter (Non-Vessels, Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_OCTADsVC = mdl_apDs_lnpwrs_OCTADsVC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_OCTADsVC = coF_apDs_lnpwrs_OCTADsVC(1) + (coF_apDs_lnpwrs_OCTADsVC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_OCTADsVC(3) * contApDs);
    est_apDs_lnpwrs_OCTADsVC = coF_apDs_lnpwrs_OCTADsVC(1) + (coF_apDs_lnpwrs_OCTADsVC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_OCTADsVC(3) * colApDs);
else
    reg_apDs_lnpwrs_OCTADsVC = (coF_apDs_lnpwrs_OCTADsVC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_OCTADsVC(2) * contApDs);
    est_apDs_lnpwrs_OCTADsVC = (coF_apDs_lnpwrs_OCTADsVC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_OCTADsVC(2) * colApDs);
end
r2_apDs_lnpwrs_OCTADsVC = calculateR2(colOCTADsVC,mdl_apDs_lnpwrs_OCTADsVC.Fitted);

figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(mdl_apDs_lnpwrs_OCTADs), box off
title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter');

subplot(2,2,2)
plot(mdl_apDs_lnpwrs_OCTADsC), box off
title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter (Detected)');

subplot(2,2,3)
plot(mdl_apDs_lnpwrs_OCTADsV), box off
title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter (Non-Vessel)');

subplot(2,2,4)
plot(mdl_apDs_lnpwrs_OCTADsVC), box off
title('Light Intensity and Aperture Diameter Effect on Max. OCTA Diameter (Non-Vessel, Detected)');

[b,bint,r,rint,stats] = regress(colOCTADs,[tbl, ones(size(tbl,1),1)]);
[b,bint,r,rint,stats] = regress(colOCTADsC,[tbl, ones(size(tbl,1),1)]);

%% 47. Model Log(Intensity) and Aperture Diameters on Histo Diameters 

cepts = true;
tbl = [log(colPwrs) colApDs];

mdl_apDs_lnpwrs_wids = fitlm(tbl, colWids, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Width (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_wids = mdl_apDs_lnpwrs_wids.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_wids = coF_apDs_lnpwrs_wids(1) + (coF_apDs_lnpwrs_wids(2) * log(contPwrs)) + (coF_apDs_lnpwrs_wids(3) * contApDs);
    est_apDs_lnpwrs_wids = coF_apDs_lnpwrs_wids(1) + (coF_apDs_lnpwrs_wids(2) * log(colPwrs)) + (coF_apDs_lnpwrs_wids(3) * colApDs);
else
    reg_apDs_lnpwrs_wids = (coF_apDs_lnpwrs_wids(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_wids(2) * contApDs);
    est_apDs_lnpwrs_wids = (coF_apDs_lnpwrs_wids(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_wids(2) * colApDs);
end
r2_apDs_lnpwrs_wids = calculateR2(colWids,mdl_apDs_lnpwrs_wids.Fitted);

mdl_apDs_lnpwrs_widsC = fitlm(tbl, colWidsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Width (Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_widsC = mdl_apDs_lnpwrs_widsC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_widsC = coF_apDs_lnpwrs_widsC(1) + (coF_apDs_lnpwrs_widsC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_widsC(3) * contApDs);
    est_apDs_lnpwrs_widsC = coF_apDs_lnpwrs_widsC(1) + (coF_apDs_lnpwrs_widsC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_widsC(3) * colApDs);
else
    reg_apDs_lnpwrs_widsC = (coF_apDs_lnpwrs_widsC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_widsC(2) * contApDs);
    est_apDs_lnpwrs_widsC = (coF_apDs_lnpwrs_widsC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_widsC(2) * colApDs);
end    
r2_apDs_lnpwrs_widsC = calculateR2(colWidsC,mdl_apDs_lnpwrs_widsC.Fitted);

mdl_apDs_lnpwrs_widsV = fitlm(tbl, colWidsV, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Width (Non-Vessel) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_widsV = mdl_apDs_lnpwrs_widsV.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_widsV = coF_apDs_lnpwrs_widsV(1) + (coF_apDs_lnpwrs_widsV(2) * log(contPwrs)) + (coF_apDs_lnpwrs_widsV(3) * contApDs);
    est_apDs_lnpwrs_widsV = coF_apDs_lnpwrs_widsV(1) + (coF_apDs_lnpwrs_widsV(2) * log(colPwrs)) + (coF_apDs_lnpwrs_widsV(3) * colApDs);
else
    reg_apDs_lnpwrs_widsV = (coF_apDs_lnpwrs_widsV(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_widsV(2) * contApDs);
    est_apDs_lnpwrs_widsV = (coF_apDs_lnpwrs_widsV(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_widsV(2) * colApDs);
end
r2_apDs_lnpwrs_widsV = calculateR2(colWidsV,mdl_apDs_lnpwrs_widsV.Fitted);

mdl_apDs_lnpwrs_widsVC = fitlm(tbl, colWidsVC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Width (Non-Vessel, Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_widsVC = mdl_apDs_lnpwrs_widsVC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_widsVC = coF_apDs_lnpwrs_widsVC(1) + (coF_apDs_lnpwrs_widsVC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_widsVC(3) * contApDs);
    est_apDs_lnpwrs_widsVC = coF_apDs_lnpwrs_widsVC(1) + (coF_apDs_lnpwrs_widsVC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_widsVC(3) * colApDs);
else
    reg_apDs_lnpwrs_widsVC = (coF_apDs_lnpwrs_widsVC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_widsVC(2) * contApDs);
    est_apDs_lnpwrs_widsVC = (coF_apDs_lnpwrs_widsVC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_widsVC(2) * colApDs);
end
r2_apDs_lnpwrs_widsVC = calculateR2(colWidsVC,mdl_apDs_lnpwrs_widsVC.Fitted);

figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(mdl_apDs_lnpwrs_wids), box off
title('Light Intensity and Aperture Diameter Effect on Histological Width');

subplot(2,2,2)
plot(mdl_apDs_lnpwrs_widsC), box off
title('Light Intensity and Aperture Diameter Effect on Histological Width (Detected)');

subplot(2,2,3)
plot(mdl_apDs_lnpwrs_widsV), box off
title('Light Intensity and Aperture Diameter Effect on Histological Width (Non-Vessel)');

subplot(2,2,4)
plot(mdl_apDs_lnpwrs_widsVC), box off
title('Light Intensity and Aperture Diameter Effect on Histological Width (Non-Vessel, Detected)');

[b,bint,r,rint,stats] = regress(colWids,[tbl, ones(size(tbl,1),1)]);
[b,bint,r,rint,stats] = regress(colWidsC,[tbl, ones(size(tbl,1),1)]);

%% 48. Model Log(Intensity) and Aperture Diameters on Histo Depth 

cepts = true;
tbl = [log(colPwrs) colApDs];

mdl_apDs_lnpwrs_deps = fitlm(tbl, colDeps, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Depth (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_deps = mdl_apDs_lnpwrs_deps.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_deps = coF_apDs_lnpwrs_deps(1) + (coF_apDs_lnpwrs_deps(2) * log(contPwrs)) + (coF_apDs_lnpwrs_deps(3) * contApDs);
    est_apDs_lnpwrs_deps = coF_apDs_lnpwrs_deps(1) + (coF_apDs_lnpwrs_deps(2) * log(colPwrs)) + (coF_apDs_lnpwrs_deps(3) * colApDs);
else
    reg_apDs_lnpwrs_deps = (coF_apDs_lnpwrs_deps(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_deps(2) * contApDs);
    est_apDs_lnpwrs_deps = (coF_apDs_lnpwrs_deps(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_deps(2) * colApDs);
end
r2_apDs_lnpwrs_deps = calculateR2(colDeps,mdl_apDs_lnpwrs_deps.Fitted);

mdl_apDs_lnpwrs_depsC = fitlm(tbl, colDepsC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Depth (Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_depsC = mdl_apDs_lnpwrs_depsC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_depsC = coF_apDs_lnpwrs_depsC(1) + (coF_apDs_lnpwrs_depsC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_depsC(3) * contApDs);
    est_apDs_lnpwrs_depsC = coF_apDs_lnpwrs_depsC(1) + (coF_apDs_lnpwrs_depsC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_depsC(3) * colApDs);
else
    reg_apDs_lnpwrs_depsC = (coF_apDs_lnpwrs_depsC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_depsC(2) * contApDs);
    est_apDs_lnpwrs_depsC = (coF_apDs_lnpwrs_depsC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_depsC(2) * colApDs);
end
r2_apDs_lnpwrs_depsC = calculateR2(colDepsC,mdl_apDs_lnpwrs_depsC.Fitted);

mdl_apDs_lnpwrs_depsV = fitlm(tbl, colDepsV, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Depth (Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_depsV = mdl_apDs_lnpwrs_depsV.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_depsV = coF_apDs_lnpwrs_depsV(1) + (coF_apDs_lnpwrs_depsV(2) * log(contPwrs)) + (coF_apDs_lnpwrs_depsV(3) * contApDs);
    est_apDs_lnpwrs_depsV = coF_apDs_lnpwrs_depsV(1) + (coF_apDs_lnpwrs_depsV(2) * log(colPwrs)) + (coF_apDs_lnpwrs_depsV(3) * colApDs);
else
    reg_apDs_lnpwrs_depsV = (coF_apDs_lnpwrs_depsV(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_depsV(2) * contApDs);
    est_apDs_lnpwrs_depsV = (coF_apDs_lnpwrs_depsV(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_depsV(2) * colApDs);
end
r2_apDs_lnpwrs_depsV = calculateR2(colDepsV,mdl_apDs_lnpwrs_depsV.Fitted);

mdl_apDs_lnpwrs_depsVC = fitlm(tbl, colDepsVC, 'VarNames', ...
    {'Light Intensity (mW)', 'Aperture Diameter (mm)', 'Histological Depth (Detected) (mm)'},'Intercept',cepts)
coF_apDs_lnpwrs_depsVC = mdl_apDs_lnpwrs_depsVC.Coefficients.Estimate;
if cepts
    reg_apDs_lnpwrs_depsVC = coF_apDs_lnpwrs_depsVC(1) + (coF_apDs_lnpwrs_depsVC(2) * log(contPwrs)) + (coF_apDs_lnpwrs_depsVC(3) * contApDs);
    est_apDs_lnpwrs_depsVC = coF_apDs_lnpwrs_depsVC(1) + (coF_apDs_lnpwrs_depsVC(2) * log(colPwrs)) + (coF_apDs_lnpwrs_depsVC(3) * colApDs);
else
    reg_apDs_lnpwrs_depsVC = (coF_apDs_lnpwrs_depsVC(1) * log(contPwrs+1)) + (coF_apDs_lnpwrs_depsVC(2) * contApDs);
    est_apDs_lnpwrs_depsVC = (coF_apDs_lnpwrs_depsVC(1) * log(colPwrs+1)) + (coF_apDs_lnpwrs_depsVC(2) * colApDs);
end
r2_apDs_lnpwrs_depsVC = calculateR2(colDepsVC,mdl_apDs_lnpwrs_depsVC.Fitted);

figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(mdl_apDs_lnpwrs_deps), box off
title('Light Intensity and Aperture Diameter Effect on Histological Depth');

subplot(2,2,2)
plot(mdl_apDs_lnpwrs_depsC), box off
title('Light Intensity and Aperture Diameter Effect on Histological Depth (Detected)');

subplot(2,2,3)
plot(mdl_apDs_lnpwrs_depsV), box off
title('Light Intensity and Aperture Diameter Effect on Histological Depth (Non-Vessel)');

subplot(2,2,4)
plot(mdl_apDs_lnpwrs_depsVC), box off
title('Light Intensity and Aperture Diameter Effect on Histological Depth (Non-Vessel, Detected)');

[b,bint,r,rint,stats] = regress(colDeps,[tbl, ones(size(tbl,1),1)]);
[b,bint,r,rint,stats] = regress(colDepsC,[tbl, ones(size(tbl,1),1)]);

%% 49. Show Partial Regression Plots for Each Lesion Dimension

% All data points included
figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,3,1), plotAdded(mdl_apDs_lnpwrs_OCTADs,3), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,4), plotAdded(mdl_apDs_lnpwrs_OCTADs,2), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,2), plotAdded(mdl_apDs_lnpwrs_wids,3), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,5), plotAdded(mdl_apDs_lnpwrs_wids,2), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,3), plotAdded(mdl_apDs_lnpwrs_deps,3), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
subplot(2,3,6), plotAdded(mdl_apDs_lnpwrs_deps,2), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
[~,ht] = suplabel('All Values','t');

% Only non-zero data points included
figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,3,1), plotAdded(mdl_apDs_lnpwrs_OCTADsC,3), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,4), plotAdded(mdl_apDs_lnpwrs_OCTADsC,2), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,2), plotAdded(mdl_apDs_lnpwrs_widsC,3), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,5), plotAdded(mdl_apDs_lnpwrs_widsC,2), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,3), plotAdded(mdl_apDs_lnpwrs_depsC,3), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
subplot(2,3,6), plotAdded(mdl_apDs_lnpwrs_depsC,2), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
[~,ht] = suplabel('Non-zero Values','t');

% Only non-large vessel data points included
figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,3,1), plotAdded(mdl_apDs_lnpwrs_OCTADsV,3), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,4), plotAdded(mdl_apDs_lnpwrs_OCTADsV,2), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,2), plotAdded(mdl_apDs_lnpwrs_widsV,3), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,5), plotAdded(mdl_apDs_lnpwrs_widsV,2), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,3), plotAdded(mdl_apDs_lnpwrs_depsV,3), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
subplot(2,3,6), plotAdded(mdl_apDs_lnpwrs_depsV,2), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
[~,ht] = suplabel('Non-large Vessel Values','t');

% Only non-large vessel and non-zero data points included
figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot(2,3,1), plotAdded(mdl_apDs_lnpwrs_OCTADsVC,3), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,4), plotAdded(mdl_apDs_lnpwrs_OCTADsVC,2), box off, ylabel('Adjusted OCTA Diameter (mm)'), title('')
subplot(2,3,2), plotAdded(mdl_apDs_lnpwrs_widsVC,3), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,5), plotAdded(mdl_apDs_lnpwrs_widsVC,2), box off, ylabel('Adjusted Histo Diameter (mm)'), title('')
subplot(2,3,3), plotAdded(mdl_apDs_lnpwrs_depsVC,3), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
subplot(2,3,6), plotAdded(mdl_apDs_lnpwrs_depsVC,2), box off, ylabel('Adjusted Histo Depth (mm)'), title('')
[~,ht] = suplabel('Non-zero, non-large Vessel Values','t');

%% 50. Make Pretty Partial Regression Plots for Paper

% now do the plots
mStyle = 'o';
mSize = 4;

volColor = [248 117 117]./255;
mixColor = [203 109 238]./255;
OCTColor = [14 170 102]./255; %[18 226 136]./255;
legColor = [17 17 17]./255;
alpha = 1;

xpos1 = [.55 .7 .85]; ypos1 = [.6875 .375 .0625]+.06; wh1 = 0.1; h1 = 0.2;

lw = 1.5; % linewidth
cw = 0.5;

figure('color','w','units','normalize','outerposition',[0 0 1 1])

% E. octa vs wid
subplot('position',[xpos1(1) ypos1(1) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_OCTADsV,3,'marker',mStyle,'color',OCTColor);
h(3).Visible = 'off';
set(h(2),'color',OCTColor,'linewidth',lw); set(h(3),'color',OCTColor);
legend off, box off
xlim([0 2]); ylim([0 4]);
ylabel('Adjusted OCTA Diameter (mm)'), xlabel('Adjusted Aperture Diameter (mm)')
title('')

% Str_OCTADsV_widsV = ['r^2 = ' num2str(r2_OCTADsV_widsV) newline 'p = ' num2str(mdl_OCTADsV_widsV.Coefficients.pValue(2))];
% text(2, 5, Str_OCTADsV_widsV, 'fontsize', 8);

% F. wids vs apDs
subplot('position',[xpos1(2) ypos1(1) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_widsV,3,'marker',mStyle,'color',volColor);
h(3).Visible = 'off';
set(h(2),'color',volColor,'linewidth',lw); set(h(3),'color',volColor);
legend off, box off
xlim([0 2])
ylabel('Adjusted Histo Width (mm)'), xlabel('Adjusted Aperture Diameter (mm)')
title('')

% Str_apDs_widsV = ['r^2 = ' num2str(r2_apDs_widsV) newline 'p = ' num2str(mdl_apDs_widsV.Coefficients.pValue(numCoFs))];
% text(1.3, 1, Str_apDs_widsV, 'fontsize', 8);

% G. deps vs apDs
subplot('position',[xpos1(3) ypos1(1) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_depsV,3,'marker',mStyle,'color',volColor);
h(3).Visible = 'off'; 
set(h(2),'color',volColor,'linewidth',lw); set(h(3),'color',volColor);
legend off, box off
xlim([0 2])
ylabel('Adjusted Histo Depth (mm)'), xlabel('Adjusted Aperture Diameter (mm)')
title('')

% Str_apDs_depsV = ['r^2 = ' num2str(r2_apDs_depsV) newline 'p = ' num2str(mdl_apDs_depsV.Coefficients.pValue(numCoFs))];
% text(1.3, 1, Str_apDs_depsV, 'fontsize', 8);

% H. OCTADs vs pwr
subplot('position',[xpos1(1) ypos1(2) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_OCTADsV,2); hold off
adjPwr1 = h(1).XData; adjOCT1 = h(1).YData;
adjPwr2 = h(2).XData; adjOCT2 = h(2).YData;
semilogx(exp(adjPwr1),adjOCT1,mStyle,'color',OCTColor); hold on
semilogx(exp(adjPwr2),adjOCT2,'color',OCTColor,'linewidth',lw)
% set(h(2),'color',OCTColor,'linewidth',lw); set(h(3),'color',OCTColor);
legend off, box off
% xlim([0 2.5])
ylabel('Adjusted OCTA Diameter (mm)'), xlabel('Adjusted Intensity (mW)')
title('')

% Str_pwrs_OCTADsV = ['r^2 = ' num2str(r2_pwrs_OCTADsV) newline 'p = ' num2str(mdl_pwrs_OCTADsV.Coefficients.pValue(numCoFs))];
% text(1, 1, Str_pwrs_OCTADsV, 'fontsize', 8);

% I. wids vs pwr
subplot('position',[xpos1(2) ypos1(2) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_widsV,2); hold off
adjPwr1 = h(1).XData; adjWid1 = h(1).YData;
adjPwr2 = h(2).XData; adjWid2 = h(2).YData;
semilogx(exp(adjPwr1),adjWid1,mStyle,'color',volColor); hold on
semilogx(exp(adjPwr2),adjWid2,'color',volColor,'linewidth',lw)
% set(h(2),'color',volColor,'linewidth',lw); set(h(3),'color',volColor);
legend off, box off
% xlim([0 2.5])
ylabel('Adjusted Histo Diameter (mm)'), xlabel('Adjusted Intensity (mW)')
title('')

% Str_pwrs_widsV = ['r^2 = ' num2str(r2_lnpwrs_widsV) newline 'p = ' num2str(mdl_lnpwrs_widsV.Coefficients.pValue(numCoFs))];
% text(1, 1, Str_pwrs_widsV, 'fontsize', 8);

% J. deps vs pwr
subplot('position',[xpos1(3) ypos1(2) wh1 h1],'color','none')
h = plotAdded(mdl_apDs_lnpwrs_depsV,2); hold off
adjPwr1 = h(1).XData; adjDep1 = h(1).YData;
adjPwr2 = h(2).XData; adjDep2 = h(2).YData;
semilogx(exp(adjPwr1),adjDep1,mStyle,'color',volColor); hold on
semilogx(exp(adjPwr2),adjDep2,'color',volColor,'linewidth',lw)
% set(h(2),'color',volColor,'linewidth',lw); set(h(3),'color',volColor);
legend off, box off
% xlim([0 2.5])
ylim([0 4])
ylabel('Adjusted Histo Depth (mm)'), xlabel('Adjusted Intensity (mW)')
title('')

% Str_pwrs_depsV = ['r^2 = ' num2str(r2_lnpwrs_depsV) newline 'p = ' num2str(mdl_lnpwrs_depsV.Coefficients.pValue(numCoFs))];
% text(1, .75, Str_pwrs_depsV, 'fontsize', 8);

wf = 2; % scaling factor for plotting measured widths/diameters of markers
df = 4; % scaling factor for plotting measured depth

% show OCTA diameters
subplot('position',[xpos1(1) ypos1(3) wh1 h1])
[eApD,eOCT] = equalIs(colApDs,colOCTADsC);
[ePwr,eOCT] = equalIs(colPwrs,colOCTADsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'o','markersize',eOCT(iLes)*wf,...
        'color',OCTColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
ylabel('Aperture Diameter (mm)')
title('OCTA Diameter (mm)')

% show Histo diameters
subplot('position',[xpos1(2) ypos1(3) wh1 h1])
[eApD,eWid] = equalIs(colApDs,colWidsC);
[ePwr,eWid] = equalIs(colPwrs,colWidsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'o','markersize',eWid(iLes)*wf,...
        'color',volColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
% ylabel('Aperture Diameter (mm)')
title('Histological Diameter (mm)')

% show Histo depths
subplot('position',[xpos1(3) ypos1(3) wh1 h1])
[eApD,eDep] = equalIs(colApDs,colDepsC);
[ePwr,eDep] = equalIs(colPwrs,colDepsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'s','markersize',eDep(iLes)*df,...
        'color',volColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
% ylabel('Aperture Diameter (mm)')
title('Histological Depth (mm)')

% make a scale bar for 5 mm
subplot('position',[xpos1(2)+(.45*wh1) ypos1(3)-.1 .1*wh1 .1*h1])
% subplot('position',[xPos(1) yPos+(1.3*h) .1*w .1*h])
plot(0,0,'o','markersize',5*wf,'color','k','linewidth',cw), hold on
% plot(0,0,'s','markersize',5*2,'color','k','linewidth',cw)
box on, axis off

print -dpdf -painters PartialRegressionsV2

%% 51. 3D Surface Plots of Multivariate Models

% Create a meshgrid of log(power) and aperture diameter
[X,Y] = meshgrid(colPwrs,colApDs);

% Plot for all values
Z_OCTA = coF_apDs_lnpwrs_OCTADs(1) + (coF_apDs_lnpwrs_OCTADs(3)*Y) + (coF_apDs_lnpwrs_OCTADs(2)*log(X));
Z_Wids = coF_apDs_lnpwrs_wids(1) + (coF_apDs_lnpwrs_wids(3)*Y) + (coF_apDs_lnpwrs_wids(2)*log(X));
Z_Deps = coF_apDs_lnpwrs_deps(1) + (coF_apDs_lnpwrs_deps(3)*Y) + (coF_apDs_lnpwrs_deps(2)*log(X));
figure('color','w','units','normalize','outerposition',[0 0 1 1])
subplot('position',[.02 .3 .3 .3]), surf(log(X),Y,Z_OCTA), title('OCTA'), xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)')
subplot('position',[.35 .3 .3 .3]), surf(log(X),Y,Z_Wids), title('Histo Width'), xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)')
subplot('position',[.68 .3 .3 .3]), surf(log(X),Y,Z_Deps), title('Histo Depth'), xlabel('Light Intensity (mW)'), ylabel('Aperture Diameter (mm)')


%% 52. Plot Aperture Diameter and Light Intensity with Marker Size as Lesion Size

OCTColor = [14 170 102]./255; %[18 226 136]./255;
volColor = [248 117 117]./255;
% volColor = [194 10 10]./255;
yelColor = [255 183 21]./255;
% yelColor = [255 210 21]./255;

% rb = redblue(2*nBin+1); rb = rb(round((2*nBin+1)/2):end,:); % red histo lesion
histColor = [linspace(255,volColor(1)*255,nBin+1);linspace(255,volColor(2)*255,nBin+1);linspace(255,volColor(3)*255,nBin+1)]'./255;
reDepths(reDepths == 0) = min(colDeps(colDeps>0)); % remove zero values for color scale
% divide depth values into bins
rbDeps = round((reDepths - min(reDepths,[],'all')) * nBin / (max(reDepths,[],'all') - min(reDepths,[],'all')))+1;

cw = 1; % circle marker width - I don't think it makes a difference maybe? I might be seeing a difference in illustrator actually
figure('color','w','units','normalize','outerposition',[0 0 .8 1])
xPos = [0.05 0.37 0.69]; yPos = .3; w = .28; h = .28; w = .2; h = .2;

wf = 2; % scaling factor for plotting measured widths/diameters of markers
df = 4; % scaling factor for plotting measured depth

% show OCTA diameters
subplot('position',[xPos(1) yPos w h])
[eApD,eOCT] = equalIs(colApDs,colOCTADsC);
[ePwr,eOCT] = equalIs(colPwrs,colOCTADsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'o','markersize',eOCT(iLes)*wf,...
        'color',OCTColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
ylabel('Aperture Diameter (mm)')
title('OCTA Diameter (mm)')

% show Histo diameters
subplot('position',[xPos(2) yPos w h])
[eApD,eWid] = equalIs(colApDs,colWidsC);
[ePwr,eWid] = equalIs(colPwrs,colWidsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'o','markersize',eWid(iLes)*wf,...
        'color',volColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
% ylabel('Aperture Diameter (mm)')
title('Histological Diameter (mm)')

% show Histo depths
subplot('position',[xPos(3) yPos w h])
[eApD,eDep] = equalIs(colApDs,colDepsC);
[ePwr,eDep] = equalIs(colPwrs,colDepsC);
for iLes = 1:length(eApD)    
    semilogx(ePwr(iLes),eApD(iLes),'s','markersize',eDep(iLes)*df,...
        'color',volColor,'linewidth',cw), hold on
    ylim([0 2.5])
    box off
end
xlabel('Light Intensity (mW)')
% ylabel('Aperture Diameter (mm)')
title('Histological Depth (mm)')

% make a scale bar for 5 mm
subplot('position',[xPos(2)+(.45*w) yPos-.1 .1*w .1*h])
% subplot('position',[xPos(1) yPos+(1.3*h) .1*w .1*h])
plot(0,0,'o','markersize',5*wf,'color','k','linewidth',cw), hold on
% plot(0,0,'s','markersize',5*2,'color','k','linewidth',cw)
box on, axis off

print -dpdf -painters MultivariateV2

%% 53. Correlation Coefficient Analyses -> Correlation Heat Maps

labels = {'OCTA Width','Histo Width','Histo Depth','Aperture','log(Power)'};

[R,P] = corrcoef([colOCTADs,colWids,colDeps,colApDs,log(colPwrs)],'Rows','Pairwise');
figure('color','w'), imagesc(R), box off
xticks(1:5); yticks(1:5);
set(gca,'XTickLabel',labels,'YTickLabel',labels);
caxis([0 1]), h = colorbar; xtickangle(60), colormap(redblue);
title('Correlations')
ylabel(h, 'Correlation Coefficient');
daspect([1 1 1]);

[RC,PC] = corrcoef([colOCTADsC,colWidsC,colDepsC,colApDs,log(colPwrs)],'Rows','Pairwise');
figure('color','w'), imagesc(RC), box off
xticks(1:5); yticks(1:5);
set(gca,'XTickLabel',labels,'YTickLabel',labels);
caxis([0 1]), h = colorbar; xtickangle(60), colormap(redblue);
title('Correlations with Non-zeros')
ylabel(h, 'Correlation Coefficient');
daspect([1 1 1]);

[RV,PV] = corrcoef([colOCTADsV,colWidsV,colDepsV,colApDs,log(colPwrs)],'Rows','Pairwise');
figure('color','w'), imagesc(RV), box off
xticks(1:5); yticks(1:5);
set(gca,'XTickLabel',labels,'YTickLabel',labels);
caxis([0 1]), h = colorbar; xtickangle(60), colormap(redblue);
title('Correlations with Non-Vessels')
ylabel(h, 'Correlation Coefficient');
daspect([1 1 1]);

[RVC,PVC] = corrcoef([colOCTADsVC,colWidsVC,colDepsVC,colApDs,log(colPwrs)],'Rows','Pairwise');
figure('color','w'), imagesc(RVC), box off
xticks(1:5); yticks(1:5);
set(gca,'XTickLabel',labels,'YTickLabel',labels);
caxis([0 1]), h = colorbar; xtickangle(60), colormap(redblue);
title('Correlations with Non-zeros and Non-Vessels')
ylabel(h, 'Correlation Coefficient');
daspect([1 1 1]);

%% Funciones Importantes
function [xunit,yunit] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
% hold off
end

function [r2,p] = calcMdlStats(mdl)
r2 = mdl.SSR/mdl.SST;
p = mdl.Coefficients.pValue;
end

% function r2 = computeR2(y,y_estimate,w)
% if nargin==2, w = ones(size(y)); end 
% sse = sum( w .* (y - y_estimate).^2 ,'omitnan'); % Sum of Squares due to Error //// Sum of Squares of residuals
% sst = sum( w .* (y - mean(y,'omitnan')).^2 ,'omitnan'); % total sum of squares
% ssr = sum( w .* (y_estimate-mean(y,'omitnan')).^2 ,'omitnan'); % sum of squares of regression 
% r2 = 1-sse/sst;
% end

function r2 = computeR2(y,y_estimate)

% if no weigths specified, then w = ones
% if nargin == 3, w = ones(size(y)); end

% Only use relevant data to calculate R2
yy = NaN(size(y));
y_est = NaN(size(y));
for i = 1:length(y)
    if isnan(y_estimate(i))
        yy(i) = NaN;
    else
        yy(i) = y(i);
    end
    
    if isnan(y(i))
        y_est(i) = NaN;
    else
        y_est(i) = y_estimate(i);
    end
    
end
yy = removeNaNs(yy); y_est = removeNaNs(y_est);
sse = sum((yy - y_est).^2 ,'omitnan'); % Sum of Squares due to Error //// Sum of Squares of residuals
sst1 = sum((yy - mean(yy,'omitnan')).^2 ,'omitnan'); % total sum of squares
ssr = sum((y_est-mean(yy,'omitnan')).^2 ,'omitnan'); % sum of squares of regression
sst2 = ssr + sse;
r2 = 1-sse/sst1;
end

function r2 = computeR2V2(y,y_estimate)

% if no weigths specified, then w = ones
% if nargin == 3, w = ones(size(y)); end

% Only use relevant data to calculate R2
yy = NaN(size(y));
for i = 1:length(y)
    if isnan(y_estimate(i))
        yy(i) = NaN;
    else
        yy(i) = y(i);
    end    
end

sse = sum((yy - y_estimate).^2 ,'omitnan'); % Sum of Squares due to Error //// Sum of Squares of residuals
sst1 = sum((yy - mean(yy,'omitnan')).^2 ,'omitnan'); % total sum of squares
ssr = sum((y_estimate-mean(yy,'omitnan')).^2 ,'omitnan'); % sum of squares of regression
sst2 = ssr + sse;
r2 = 1-sse/sst2;
end

function R2 = calculateR2(z,z_est)
% calcuateR2 Cacluate R-squared
% R2 = calcuateR2(z,z_est) takes two inputs - The observed data z and its
% estimation z_est (may be from a regression or other model), and then
% compute the R-squared value a measure of goodness of fit. R2 = 0
% corresponds to the worst fit, whereas R2 = 1 corresponds to the best fit.
% 
% Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)

yy = NaN(size(z));
for i = 1:length(z)
    if isnan(z_est(i))
        yy(i) = NaN;
    else
        yy(i) = z(i);
    end    
end

z = yy;

xx = NaN(size(z));
for i = 1:length(z)
    if isnan(z(i))
        xx(i) = NaN;
    else
        xx(i) = z_est(i);
    end    
end

z_est = xx;

z = removeNaNs(z);
z_est = removeNaNs(z_est);

r = z-z_est;
normr = norm(r);
SSE = normr.^2;
SST = norm(z-mean(z))^2;
R2 = 1 - SSE/SST;

end

% this function makes x and y vectors so that same elements are NaNs, then
% removes NaNs
function [xx, yy] = equalIs(x,y)

xx = NaN(size(x));
for i = 1:length(x)
    if isnan(y(i))
        xx(i) = NaN;
    else
        xx(i) = x(i);
    end    
end

yy = NaN(size(y));
for i = 1:length(y)
    if isnan(xx(i))
        yy(i) = NaN;
    else
        yy(i) = y(i);
    end    
end

xx = removeNaNs(xx);
yy = removeNaNs(yy);
end

function colVec = removeNaNs(colVec)

% remove NaN elements
rows = any(isnan(colVec),2);
colVec(rows) = [];

end
