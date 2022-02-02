load('D:\data\pt\PT2_LesionBounds.mat')

pt2_L = 256*lesions(200:800, 600:1400, :);

%%
pt2_L1 = pt2_L(100:400, 300:550, 36);
pt2_L1 = imrotate(pt2_L1, -21);
pt2_L1 = pt2_L1(123:end, :);
pt2_L1 = imresize(pt2_L1, 2.5, 'nearest');
pt2_L1 = pt2_L1(5:504, 35:834);
image(pt2_L1)

%%
pt2_L2 = pt2_L(100:550, 1:350, 36);
pt2_L2 = imrotate(pt2_L2, -47);
pt2_L2 = pt2_L2(300:end, :);
pt2_L2 = imresize(pt2_L2, 2.5, 'nearest');
pt2_L2 = pt2_L2(46:545, 270:1069);
image(pt2_L2)

%%
pt2_L3 = pt2_L(11:400, 450:end, 20);
pt2_L3 = imrotate(pt2_L3, -22);
pt2_L3 = pt2_L3(177:end, :);
pt2_L3 = imresize(pt2_L3, 2.5, 'nearest');
pt2_L3 = pt2_L3(1:500, 201:1000);
image(pt2_L3)

%%
pt2_L4 = 256*lesions(403:710, 600:1200, 21);
pt2_L4 = imrotate(pt2_L4, -33);
pt2_L4 = pt2_L4(206:end-150, :);
pt2_L4 = imresize(pt2_L4, 2.5, 'nearest');
pt2_L4 = pt2_L4(2:501, 451:1250);
image(pt2_L4)
% intensity: 2.64

%%
pt2_L5 = pt2_L(350:600,1:300,20);
pt2_L5 = imrotate(pt2_L5, -38);
pt2_L5 = pt2_L5(133:end, :);
pt2_L5 = imresize(pt2_L5, 2.5, 'nearest');
pt2_L5 = pt2_L5(2:501, 110:909);
image(pt2_L5)


%%
% pt2_L7 = 0;

%% PT2 R
pt2_R = 256*lesions(1:800, 1400:end, :);
% vis_bounds_pt3(pt2_R)

%%
pt2_R2 = pt2_R(500:end, 300:end, 1);
pt2_R2 = imrotate(pt2_R2, 31);
pt2_R2 = pt2_R2(150:end, :);
pt2_R2 = imresize(pt2_R2, 2.5, 'nearest');
pt2_R2 = pt2_R2(630:1129, 810:1609);
image(pt2_R2)
%%

pt2_R4 = 256*lesions(400:800, 1900:2200, 20);
pt2_R4 = imrotate(pt2_R4, 34);
pt2_R4 = imresize(pt2_R4, 2.5, 'nearest');
pt2_R4 = pt2_R4(400:899, 120:919);
image(pt2_R4)
%%
pt2_R5 = pt2_R(600:end, 600:end, 19);
pt2_R5 = imrotate(pt2_R5, 30);
pt2_R5 = pt2_R5(295:end, :);
pt2_R5 = imresize(pt2_R5, 2.5, 'nearest');
pt2_R5 = pt2_R5(7:506, 271:1070);
image(pt2_R5)

%%
pt2_R6 = pt2_R(1:410, 100:550, 35);

pt2_R6 = imrotate(pt2_R6, 29);
pt2_R6 = pt2_R6(378:end, :);
pt2_R6 = imresize(pt2_R6, 2.5, 'nearest');
pt2_R6 = pt2_R6(1:end-3, 425:1224);
image(pt2_R6)

%%
pt2_R7 = pt2_R(290:800, 550:800, 31);
pt2_R7 = imrotate(pt2_R7, 40);
pt2_R7 = pt2_R7(233:end, :);
pt2_R7 = imresize(pt2_R7, 2.5, 'nearest');
pt2_R7 = pt2_R7(21:520, 230:1029);
image(pt2_R7);

% pixel size was .025
% intensity: 2.24

%%
load("D:\data\pt\PT2L6_LesionBounds.mat")
pt2_L6 = 256*lesions(100:300, 200:600);
pt2_L6 = imrotate(pt2_L6,  -34);
pt2_L6 = imresize(pt2_L6, 3.96, 'nearest');
pt2_L6 = pt2_L6(567:1066, 372:1171);
image(pt2_L6)
daspect([1 1 1])

%%
load('D:\data\pt\PT3_LesionBounds.mat')

pt3_L = 256*lesions(:, 1:800, :);
pt3_R = 256*lesions(150:400, 900:end, :);

%%
pt3_L4 = imrotate(pt3_L(300:430, 310:700, size(pt3_L, 3)+1-30), -36);
pt3_L4 = imresize(pt3_L4, 3.39, 'nearest');
pt3_L4 = pt3_L4(400:899, 230:1029);
image(pt3_L4)

%%
pt3_L2 = pt3_L(200:400, 380:700, end);
pt3_L2 = imresize(pt3_L2, 3.39, 'nearest');;
pt3_L2 = imrotate(pt3_L2, -41);
pt3_L2 = pt3_L2(565:1064, 251:1050);
image(pt3_L2)

%%
pt3_L6 = pt3_L(200:400, 500:700,13);
pt3_L6 = imresize(pt3_L6, 3.39, 'nearest');;
pt3_L6 = imrotate(pt3_L6, -28);
pt3_L6 = pt3_L6(420:919, 1:800);
image(pt3_L6)

%%
pt3_L7 = pt3_L(350:500, 200:600, 10);
pt3_L7 = imresize(pt3_L7, 3.39, 'nearest');;
pt3_L7 = imrotate(pt3_L7, -28);
pt3_L7 = pt3_L7(411:910, 355:1154);
image(pt3_L7)


%%
pt3_L1 = pt3_L(1:350, 620:end, 44);
pt3_L1 = imresize(pt3_L1, 3.39, 'nearest');;
pt3_L1 = imrotate(pt3_L1, -22);
pt3_L1 = pt3_L1(741:1240, 79:878);
image(pt3_L1)

%%
pt3_R1 = pt3_R(1:200,1:500,19);
pt3_R1 = imrotate(pt3_R1, 25);
pt3_R1 = pt3_R1(183:end, :);
pt3_R1 = imresize(pt3_R1, 3.39, 'nearest');
pt3_R1 = pt3_R1(1:500, 751:1550);
image(pt3_R1)

%%
pt3_R2 = pt3_R(150:end,250:end,17);
pt3_R2 = imrotate(pt3_R2, 23);
pt3_R2 = pt3_R2(246:end, :);
pt3_R2 = imresize(pt3_R2, 3.39, 'nearest');
pt3_R2 = pt3_R2(7:506, 583:1382);
pt3_R2(:, 1:100) = 0;
image(pt3_R2)

%%
pt3_R4 = 256*imrotate(lesions(250:411, 1200:1400, size(lesions, 3)+1-24), 21);
pt3_R4 = pt3_R4(50+15:end-12, 27-16:26+220, :);
pt3_R4 = imresize(pt3_R4, 3.39, 'nearest');
pt3_R4 = pt3_R4(3:end, 1:end-1);
image(pt3_R4)

%%
load('D:\data\pt\1.10.2_PT3L3_Edge.mat')
pt3_L3 = 256*lesions(1:600, 550:800);
pt3_L3 = imresize(pt3_L3, 3.39, 'nearest');
pt3_L3 = imrotate(pt3_L3, -20);
pt3_L3 = pt3_L3(883:1382, 316:1115);
image(pt3_L3)

%%
load('D:\data\pt\PT4_LesionBounds.mat')

pt4_L4 = 256*lesions(:,:,11);
pt4_L4 = imrotate(pt4_L4, -20);
pt4_L4 = imresize(pt4_L4, 7.9, 'nearest');
pt4_L4 = pt4_L4(1761:2260, 2154:2953);
image(pt4_L4)

%%
load('D:\data\pt\PT5_LesionBounds.mat')

pt5_R4 = imrotate(256*lesions(212:400, 1100:1319,  size(lesions, 3)+1-17), 14);
pt5_R4 = imresize(pt5_R4, 3.9, 'nearest');
pt5_R4 = pt5_R4(240:739, 140:939);
image(pt5_R4)