%
% c_d_w_a = readmatrix("C:\Users\Julien\lab\pt_3d_reconstruction\modeling\c_depths_and_widths_and_areas.csv");
% d_w_a = readmatrix("C:\Users\Julien\lab\pt_3d_reconstruction\modeling\depths_and_widths_and_areas.csv");

for i=2:5
    h_vals = d_w_a(:, i);
    bads = ~isnan(h_vals);
    h_vals = h_vals(bads);
    c_vals = c_d_w_a(:, i);
    c_vals = c_vals(bads);
    log_c_vals = log(c_vals);
    
    
    disp(i)
    disp(calculateR2(h_vals, c_vals));
    mdl = fitlm(log_c_vals, h_vals, 'intercept', false);
    disp(mdl.Rsquared.Ordinary)
    
%     plot(mdl)
%     hold on
    
    mdl = fitlm(c_vals, h_vals, 'intercept', false);
    disp(mdl.Rsquared.Ordinary)
    disp('')
    
    
    mdl = fitlm(c_vals, h_vals, 'intercept', true);
    disp(mdl.Rsquared.Ordinary)
    disp('')
%     plot(mdl)
    
end