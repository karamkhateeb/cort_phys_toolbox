get_all_slices;
load_lesion_metadata;
%%

global model_C_data_p5mm
global model_C_data_1mm
global model_C_data_1p5mm
global model_C_data_2mm
% 
% model_C_data_p5mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_p5mm\gridsearchbestvalid_p5mm.mc2", 2);
% model_C_data_1mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_1mm\gridsearchbestvalid_1mm.mc2", 2);
% model_C_data_1p5mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_1p5mm\gridsearchbestvalid_1p5mm.mc2", 2);
% model_C_data_2mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_2mm\gridsearchbestvalid_2mm.mc2", 2);

grid_i = 228;
model_C_data_p5mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\p5mm\grid_"+ ...
    grid_i+"\p5mmgrid_"+grid_i+".mc2", 1);
model_C_data_1mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\1mm\grid_"+ ...
    grid_i+"\1mmgrid_"+grid_i+".mc2", 1);
model_C_data_1p5mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\1p5mm\grid_"+ ...
    grid_i+"\1p5mmgrid_"+grid_i+".mc2", 1);
model_C_data_2mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\2mm\grid_"+ ...
    grid_i+"\2mmgrid_"+grid_i+".mc2", 1);

%%
loop_c_depths_and_widths;
% c_d_w_a = readmatrix("Users/julienbloch/lab/pt_3d_reconstruction/modeling/c_depths_and_widths_and_areas.csv");
loop_depths_and_widths;
% d_w_a = readmatrix("Users/julienbloch/lab/pt_3d_reconstruction/modeling/depths_and_widths_and_areas.csv");
%%

c_depths = c_d_w_a(:, 1);
h_depths = d_w_a(:, 1);
sqrt_c_depths = c_depths.^.5;
sqp_c_depths = c_depths.^1.5;
depth_mdl = fitlm(sqrt_c_depths, h_depths, 'intercept', false);
depth_scaling = depth_mdl.Fitted./c_depths;

%%

c_widths = c_d_w_a(:, 2);
h_widths = d_w_a(:, 2);
sqp_c_widths = c_widths.^1.5;

width_mdl = fitlm(c_widths, h_widths, 'intercept', false);
width_scaling = width_mdl.Fitted./c_widths;

%%
for i=1:1
    h_vals = d_w_a(:, i);
    bads = ~isnan(h_vals);
    h_vals = h_vals(bads);
    c_vals = c_d_w_a(:, i);
    c_vals = c_vals(bads);
    log_c_vals = log(c_vals);
    sqrt_c_vals = c_vals.^.5;

%     boxcoxlm(c_vals, h_vals);
    
    disp(i)
    disp(calculateR2(h_vals, c_vals));
    
%     mdl = fitlm(log_c_vals, h_vals, 'intercept', false);    
%     disp(calculateR2(h_vals, mdl.Fitted));
    
    mdl = fitlm(sqrt_c_vals, h_vals, 'intercept', false);
    disp(calculateR2(h_vals, mdl.Fitted));
%     

%     mdl = fitlm(sqrt_c_vals, h_vals, 'purequadratic', 'intercept', true);
%     disp(calculateR2(h_vals, mdl.Fitted));
    
    mdl = fitlm(c_vals, h_vals, 'intercept', false);
    disp(calculateR2(h_vals, mdl.Fitted));
%     
%     mdl = fitlm(c_vals, h_vals, 'intercept', true);
%     disp(calculateR2(h_vals, mdl.Fitted));
    
%     plot(mdl)
%     hold on
    
    %     mdl = fitlm(c_vals, h_vals, 'intercept', true);
%     disp(mdl.Rsquared.Ordinary)
%     disp('')

    scatter(c_vals, h_vals)
    hold on
    plot(1:350, predict(mdl, [(1:350)].'));
%     plot((1:350), predict(mdl, [sqrt(1:350)].'));
%     plot(mdl)
end