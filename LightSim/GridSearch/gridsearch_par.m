% % Find best threshold to explain all lesions
% % generate an average IOU over slices for each threshold
get_all_slices;
load_lesion_metadata;

%%
grid_IOUs = zeros(2, 400);
%%
apertures = ["1mm", "1p5mm", "2mm", "p5mm"];
fun = @get_best_threshold;
parfor (grid_i = 1:400, 2)

    model_C_data_p5mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\p5mm\grid_"+ ...
        grid_i+"\p5mmgrid_"+grid_i+".mc2", 1);
    model_C_data_1mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\1mm\grid_"+ ...
        grid_i+"\1mmgrid_"+grid_i+".mc2", 1);
    model_C_data_1p5mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\1p5mm\grid_"+ ...
        grid_i+"\1p5mmgrid_"+grid_i+".mc2", 1);
    model_C_data_2mm = load_and_clean_contours("D:\lab\mcxoutput\gridsearch\2mm\grid_"+ ...
        grid_i+"\2mmgrid_"+grid_i+".mc2", 1);
    
    close all
    fun = @(threshold)get_threshold_acc(threshold, lesion_metadata, ...
        model_C_data_p5mm, model_C_data_1mm, model_C_data_1p5mm, model_C_data_2mm);

    % call bayesopt on the function
    threshs = optimizableVariable('threshold', [1e-7 1e-2],'Transform','log');
    
    
    results = bayesopt(fun, [threshs], ...
        'IsObjectiveDeterministic',true, ...
        'MaxObjectiveEvaluations', 20, ...
        'ExplorationRatio', 0.4, ...
        'AcquisitionFunctionName', 'expected-improvement');
    
    newvals = [results.MinObjective, results.XAtMinObjective{1, 1}];
    grid_IOUs(:, grid_i) = newvals;
end

%%
function obj = get_threshold_acc(threshold, lesion_metadata, ...
        model_C_data_p5mm, model_C_data_1mm, model_C_data_1p5mm, model_C_data_2mm)
%     IOUs = [];
    c_maxdepths = [];
    maxdepths = [];
    c_avgwidths = [];
    avgwidths = [];
    
    for les_i = 1:22
        lmd = lesion_metadata(les_i);
        
        switch lmd.aperture
            case 0.5
                model_C_data = model_C_data_p5mm;
            case 1
                model_C_data = model_C_data_1mm;
            case 1.5
                model_C_data = model_C_data_1p5mm;
            case 2
                model_C_data = model_C_data_2mm;
        end

        cntr = get_model_contour(lmd.aperture, lmd.intensity, threshold.(1), model_C_data);
        
        [c_maxdepth, c_avgwidth, ~, ~] = get_d_w_a_contour_slice(cntr);
        [maxdepth, avgwidth, ~, ~] = get_lesion_maxdepths_and_avgwidths(lmd.slice);
%         show_IOU(cntr, lmd.slice);

        c_maxdepths = [c_maxdepths c_maxdepth];
        maxdepths = [maxdepths maxdepth];
        c_avgwidths = [c_avgwidths c_avgwidth];
        avgwidths = [avgwidths avgwidth];
        
        
    end
    depth_r2 = calculateR2(maxdepths, c_maxdepths);
    width_r2 = calculateR2(avgwidths, c_avgwidths);
    obj = ((depth_r2-1)^2 + (width_r2-1)^2);
end