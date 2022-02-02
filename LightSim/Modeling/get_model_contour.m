function best_contour = get_model_contour(aperture, intensity, threshold, model_C_data)
    %% Get contour lines from model to compare to histology
    % assumes contours have been loaded into memory for speedup
%     switch aperture
%         case 0.5
%             global model_C_data_p5mm
%             model_C_data = model_C_data_p5mm;
%         case 1
%             global model_C_data_1mm
%             model_C_data = model_C_data_1mm;
%         case 1.5
%             global model_C_data_1p5mm
%             model_C_data = model_C_data_1p5mm;
%         case 2
%             global model_C_data_2mm
%             model_C_data = model_C_data_2mm;
%     end

    %% get modified threshold
    mod_threshold = threshold / intensity;

    
    %% find contour closest to modified threshold
    closest_level_error = inf;
    for model_c_i = 1:size(model_C_data, 2)
        model_c = model_C_data(model_c_i);
        if abs(model_c.level - mod_threshold) < closest_level_error
            best_contour = model_c;
            closest_level_error = abs(model_c.level - mod_threshold);
        end
    end

%     plot(best_contour.xdata, best_contour.ydata)
%     xlim([1, 800])
%     ylim([1, 500])
end