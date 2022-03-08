function [x_hat, r_hat, hist, hist_r] = P_VDAMP(dcoil, mask, prob_map, x0, csm, opts)
%     The Parallel Variable Density Approximate Message Passing (P-VDAMP) algorithm for
%     reconstruction of multi-coil MRI data
%     IN:
%         dcoil: (nx*ny*nc) measurement data i.e. noisy undersampled Fourier coefficients...
%             unsampled coefficients have entry zero
%         mask: (nx*ny) sampling mask of zeros and ones
%         prob_map: (nx*ny) expectation of random mask
%         x0: (nx*ny) reference image for tracking reconstruction progress
%         csm: (nx*ny*nc) coil sensitivities
%         opts: options object with attributes
%             maxIter: maximum number of allowed iterations; default 20
%             maxTime: maximum allowed time in second; default 50
%             verbose: if 1, print progress metrics; default 0
%             scales: number of wavelet decomposition scales; default 4 
%             wavType: choose mother wavelet; default 'db4'
%             rho: damping factor 0<rho<1, where 1 is no damping; default 1
%             saveHist: if 1, save detailed history of reconstruction (see hist object below); default 0
%             tauStop: if 1, stop iteration when tau increases or changes by less than 0.01%; default 1
%             sigmaMeas: (nc*nc) measurement noise covariance; default 0
%             
%      OUT:
%         x_hat: P-VDAMP's estimate of x0 
%         hist: 
%           tau_mean: (maxIter) mean of tau, MSE estimate of r
%           x_mse: (maxIter) MSE of biased estimate
%           timer: (maxIter) timer of biased estimate
%           if saveHist == 1, also saves
%               r: (nx*ny*maxIter) image representation of vector subject to thresholding
%               tauPy: (nx*ny*maxIter) pyramid representation of wavelet-domain covariance              
%         hist_r:
%           timer: (maxIter) timer of unbiased estimate
%           x_mse: (maxIter) MSE of unbiased estimate
%            
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2022  Charles Millard
% Copyright (C) 2022  Siemens Healthineers           
                
    [nx, ny, ~] = size(dcoil);
     
    opts = setDefaultOpts(opts);

    % prepare data types
    tau = cell(opts.scales, 1); 
    err = cell(opts.scales, 1);
    C_tilde= cell(opts.scales, 1);
    alpha = cell(opts.scales, 1);
    for s =1:opts.scales
        tau{s} = struct('H', 0, 'V', 0, 'D', 0);
        err{s} = struct('H', 0, 'V', 0, 'D', 0);
        C_tilde{s} = struct('H', [], 'V', [], 'D', []);
        alpha{s} = struct('H', 0, 'V', 0, 'D', 0);        
    end  
    tau{opts.scales}.A = 0;
    err{opts.scales}.A = 0;
    C_tilde{opts.scales}.A = [];
    alpha{opts.scales}.A = 0;   
    
    % compute amount of padding to cleanly divide wavelets
    coarse_wavsz = 2^opts.scales;
    nx_pad = nx + (coarse_wavsz - mod(nx, coarse_wavsz))*(mod(nx, coarse_wavsz) > 0);
    ny_pad = ny + (coarse_wavsz - mod(ny, coarse_wavsz))*(mod(ny, coarse_wavsz) > 0);
    pad = @(A) padarray(A, [(nx_pad-nx)/2, (ny_pad-ny)/2]);
    unpad = @(A) A((nx_pad-nx)/2 + 1:end - (nx_pad-nx)/2, (ny_pad-ny)/2 + 1:end - (ny_pad-ny)/2, :);

    csm = pad(csm);

    % allocate memory for hist attributes
    if opts.saveHist
        hist.r = zeros(nx_pad, ny_pad, opts.maxIter);
        hist.tauPy = zeros(nx_pad, ny_pad, opts.maxIter);
    end
    C_thr_tracker = zeros(nx_pad, ny_pad, opts.maxIter);

    % for masking MSE computations
    im_mask = (abs(x0) > max(abs(x0(:)))/20); 

    % compute wavelet spectra
    [psdH, psdV, psdD, psdA] = WavSpec2d(opts.scales, opts.wavType, nx, ny);
    
    % compute xi for error propagation function tau_update
    xi = WavXi(csm, opts.wavType, opts.scales);

    %% **** start P-VDAMP **** 
    time_init = tic; 
    
    % unbiased initialisation
    inv_p = prob_map.^(-1);
    r = pad(ifftnc(inv_p.*mask.*dcoil));
    r = sum(conj(csm) .*r, 3);
    
    % error propagation
    tau = tau_update(xi, dcoil, prob_map, mask, psdH, psdV, psdD, psdA, opts.sigmaMeas);
    
    for iter = 1:opts.maxIter

        % save hist attributes
        tau_py = pyramid(tau);
        hist.tau_mean(iter) = mean(tau_py(:));
        hist_r.x_mse(iter) = immse(im_mask.*unpad(r),im_mask.*x0);
        hist_r.timer(iter) = toc(time_init);
        if opts.saveHist
            hist.r(:,:,iter) = r; 
            hist.tauPy(:,:,iter) = tau_py;
        end

        if opts.verbose
            disp(['tau mean at iteration ', num2str(iter), ' is ', num2str(hist.tau_mean(iter))])
        end

        % *** SURE-optimal soft thresholding ***
        C = multiscaleDecomp(r, opts.scales, opts.wavType);
        [C_thr, ~, df] = multiscaleSUREsoft(C, tau);
        
        % damping
        if iter > 1 
            C_damp = opts.rho*pyramid(C_thr) + (1-opts.rho)*C_damp;
        else
            C_damp = pyramid(C_thr);
        end
        C_thr = pyramidInv(C_damp, opts.scales);            

        % *** Wavelet domain per-subband Onsager correction ***
        for s=1:opts.scales
            subbands = fieldnames(C_thr{s});           
            for i = 1:numel(subbands)
                b_name = subbands{i}; % subband name
                if iter > 1 % damped
                    alpha{s}.(b_name) = opts.rho*mean(df{s}.(b_name))/2;
                else
                    alpha{s}.(b_name) = mean(df{s}.(b_name))/2;
                end
                C_tilde{s}.(b_name) =  C_thr{s}.(b_name) - alpha{s}.(b_name).*C{s}.(b_name);
                C_tilde{s}.(b_name) = C_tilde{s}.(b_name)./(1 - alpha{s}.(b_name));
            end
        end
 
        % *** Density compensated gradient descent ***
        r_tilde = multiscaleRecon(C_tilde, opts.wavType);
        z = mask.*(dcoil - fftnc(unpad(csm .* r_tilde)));
        r = r_tilde + sum(conj(csm) .*pad(ifftnc(inv_p .* z)),3); 
        
        % save denoiser wavelets and timer
        C_thr_tracker(:,:,iter) = pyramid(C_thr);
        hist.timer(iter) = toc(time_init);   

        % check time stoppping criterion
        if hist.timer(iter) > opts.maxTime
            break
        end  
        
        % check tau stopping criterion
        if opts.tauStop && iter>1
            if hist.tau_mean(iter) > hist.tau_mean(iter-1)
                disp(['Stopping at iteration ', num2str(iter), ' due to tau increase'])
                break
            elseif abs(hist.tau_mean(iter) - hist.tau_mean(iter-1))/hist.tau_mean(iter-1) < 0.001
                disp(['Stopping at iteration ', num2str(iter), ' as tau changed by less than 0.01%'])
                break
            end
        end 
        
        % coloured noise power spectrum re-estimation
        tau = tau_update(xi, z, prob_map, mask, psdH, psdV, psdD, psdA, opts.sigmaMeas);
    end
    
    % retrospectively calculate MSE of biased estimate from unweighted gradient descent
    % (done retrospectively so timer is fairer reflection of reality)
    for it = 1:iter
        C_thr = pyramidInv(C_thr_tracker(:,:,it), opts.scales);
        xk_tilde = multiscaleRecon(C_thr, opts.wavType);
        z = mask.*(dcoil - fftnc(unpad(csm .* xk_tilde)));
        xk = xk_tilde + sum(conj(csm) .*pad(ifftnc(z)),3); 
        xk = unpad(xk);
        hist.x_mse(it) = immse(im_mask.*xk, im_mask.*x0);
    end

    % trim hist attributes in case a stopping criterion reached
    if opts.saveHist
        hist.r = hist.r(:, :, hist.x_mse>0);
        hist.tauPy = hist.tauPy(:, :, hist.x_mse>0);
    end
    
    % unbiased estimate
    r_hat = unpad(r);

    % biased estimate
    x_hat = xk;
end


