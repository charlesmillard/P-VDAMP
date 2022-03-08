function tau = tau_update(xi, z, prob_map, mask, psdH, psdV, psdD, psdA, sigma_meas)
    % update wavelet-domain covariance model
    scales = size(xi, 1);
    nc = size(z, 3);
    
    nx_pad = 2*size(xi{1}.H, 1);
    ny_pad = 2*size(xi{1}.H, 2);
    
    tau = cell(scales, 1); 
    for s =1:scales
        tau{s} = struct('H', 0, 'V', 0, 'D', 0);
    end
    tau{scales}.A = 0;
    
    for s = 1:scales  
        flat_csm = reshape(xi{s}.H, [nx_pad*ny_pad/2^(2*s), nc]);
        y_cov = multicoil_z_cov(z, mask, prob_map, psdH(:,:,s), sigma_meas);
        tau{s}.H = real(reshape(sum(conj(flat_csm) .* (flat_csm *y_cov), 2), [nx_pad/2^s ny_pad/2^s]));           

        flat_csm = reshape(xi{s}.V, [nx_pad*ny_pad/2^(2*s), nc]);
        y_cov = multicoil_z_cov(z, mask, prob_map, psdV(:,:,s), sigma_meas);
        tau{s}.V = real(reshape(sum(conj(flat_csm) .* (flat_csm *y_cov), 2), [nx_pad/2^s ny_pad/2^s]));     
        
        flat_csm = reshape(xi{s}.D, [nx_pad*ny_pad/2^(2*s), nc]);
        y_cov = multicoil_z_cov(z, mask, prob_map, psdD(:,:,s), sigma_meas);
        tau{s}.D = real(reshape(sum(conj(flat_csm) .* (flat_csm *y_cov), 2), [nx_pad/2^s ny_pad/2^s]));
    end
    
    flat_csm = reshape(xi{scales}.A, [nx_pad*ny_pad/2^(2*scales), nc]);
    y_cov = multicoil_z_cov(z, mask, prob_map, psdA, sigma_meas);
    tau{scales}.A = real(reshape(sum(conj(flat_csm) .* (flat_csm *y_cov), 2), [nx_pad/2^s ny_pad/2^s]));
end