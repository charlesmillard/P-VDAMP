function z_cov = multicoil_z_cov(dcoil, mask, prob_map, psd, sigma_meas)
    % compute covariance term which is a part of tau update
    weight = (1./prob_map - 1)./prob_map;

    nc = size(dcoil, 3);
    ns = nnz(mask);
    
    sigma_meas = permute(repmat(sigma_meas, [1, 1, ns]), [3, 1, 2]);  
    sigma_meas = sigma_meas./prob_map(logical(mask)) ;

    mask = logical(mask.*ones(size(dcoil)));
    weight = weight.*ones(size(dcoil));
    
    y = reshape(dcoil(mask), [ns, nc]);
    w2 = reshape(weight(mask), [ns, nc]);
    p2 = repmat(psd(mask(:,:,1)), [1, nc]);

    z_cov = (p2 .* w2 .*  y )' * y + squeeze(sum(p2.*sigma_meas, 1));
end

