function model_good_contours = load_and_clean_contours(filename, numtimesteps)
    model_data = loadmc2(filename, [800, 800, 563, numtimesteps]);
    %%
    model_data = squeeze(model_data(:, 400, :, :)) + ...
        squeeze(model_data(400, :, :, :)) + ...
        squeeze(model_data(401, :, :, :)) + squeeze(model_data(:, 401, :, :));
    model_data = sum(model_data, 3);
    model_data = model_data / sum(model_data(376:425, 61)); % the total light intensity entering is normalized to 1
    
    blurred_model_data = imgaussfilt(model_data, 4) + 1e-10;
    
    log_model_data = log10(blurred_model_data(:, 62:end-2));
    flat_log_model_data = imrotate(log_model_data, 90);
    conts = contourc(flat_log_model_data, 500);
    model_C_data = contourdata(conts);

    %% clean up contours
    model_good_contours = struct();
    i = 1;
    for model_c_i = 100:size(model_C_data, 2)
        model_c = model_C_data(model_c_i);

        % skip if closed
        if ~model_c.isopen
            continue
        end

        % skip if not passing through center
        if min(model_c.xdata) > 400 || max(model_c.xdata) < 400
            continue
        end

        % can't reach the full depth
        if min(model_c.ydata) == 1
            continue
        end
        
        % make sure it connects to origin
        model_c.xdata = [min(model_c.xdata); model_c.xdata; max(model_c.xdata)];
        model_c.ydata = [500; model_c.ydata; 500];

        % assign values
        model_good_contours(i).xdata = model_c.xdata;
        model_good_contours(i).ydata = model_c.ydata;
        model_good_contours(i).level = 10^(model_c.level);

        i = i + 1;
    end
end