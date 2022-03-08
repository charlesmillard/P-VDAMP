function [mask, prob_map] = load_mask(nx, ny, delta, sample_type)
    % load mask and probability map
    
    mask_loc = ['masks/', num2str(nx), 'x', num2str(ny), '/', num2str(round(1/delta))];
    try
        load([mask_loc, '/prob_map.mat']) % load in variable density prob
    catch
        error('For this demo there is only masks at R=5 and R=10. For others you will need to make your own.')
    end

    if strcmp(sample_type, 'bern')
        mask = binornd(1, prob_map, nx,ny);
    elseif strcmp(sample_type, 'pd')
        load([mask_loc, '/mask.mat'])
    else
        error('Invalid sample_type - please pick bern or pd')
    end

end
