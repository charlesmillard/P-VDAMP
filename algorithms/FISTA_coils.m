function [x_hat, hist] = FISTA_coils(dcoil, csm, mask, x0, opts)
%    The FISTA or SURE-IT algorithm for
%    image reconstruction from multi-coil k-space measurements, where
%    a sparse model on wavelet coefficients is assumed
%    IN:
%         dcoil: (nx*ny*nc) measurement data i.e. noisy undersampled Fourier coefficients...
%             unsampled coefficients have entry zero
%         csm: (nx*ny*nc) coil sensitivity maps
%         mask: (nx*ny) sampling mask of zeros and ones
%         x0: (nx*ny) reference image for tracking reconstruction progress
%         opts: options object with attributes
%             maxIter: maximum number of allowed iterations; default 50
%             maxTime: maximum allowed time in second; default 50
%             tauStop: if 1, stop iteration when tau increases or changes by less than 0.01%; default 1
%             verbose: if 1, print progress metrics; default 0
%             scales: number of wavelet decomposition scales; default 4    
%             wavType: choose mother wavelet; default 'db4'
%             saveHist: if 1, save detailed history of reconstruction; default 0
%             SURE: if 1, use SURE to automatically tune thresholds; default 1
%             lambda: if not using SURE, specify (global) sparse weighting
%     OUT:
%         x_hat: FISTA's estimate x0 
%         var_final_est: final message precision in k-space
%         hist: history objects with attributes
%           x_mse: (maxIter) ground truth MSE of x   
%           tau_mean: (maxIter) MSE of pre-thresholded image estimate
%           timer: (maxIter) timer of estimate
%               if saveHist == 1, also
%                r: (nx*ny*maxIter) image representation of vector subject to thresholding
%                x: (nx*ny*maxIter) image estimate
%                        
    
    [nx, ny, ~] = size(dcoil);

    opts = setDefaultOptsFISTA(opts);

    % compute amount of padding to cleanly divide wavelets
    coarse_wavsz = 2^opts.scales;
    nx_pad = nx + (coarse_wavsz - mod(nx, coarse_wavsz))*(mod(nx, coarse_wavsz) > 0);
    ny_pad = ny + (coarse_wavsz - mod(ny, coarse_wavsz))*(mod(ny, coarse_wavsz) > 0);
    pad = @(A) padarray(A, [(nx_pad-nx)/2, (ny_pad-ny)/2]);
    unpad = @(A) A((nx_pad-nx)/2 + 1:end - (nx_pad-nx)/2, (ny_pad-ny)/2 + 1:end - (ny_pad-ny)/2, :);
    
    csm = pad(csm);
    
    % allocate memory for hist attributes
    if opts.saveHist
        hist.r = zeros(nx_pad,ny_pad,opts.maxIter);
        hist.x = zeros(nx_pad,ny_pad,opts.maxIter);
    end

    im_mask = (abs(x0) > max(abs(x0(:)))/20);
    
    
    time_init = tic; 

    r = sum(conj(csm).*pad(ifftnc(mask.*dcoil)), 3);

    var_est = immse(x0, unpad(r));
    tau = cell(opts.scales, 1); 
    lambda = cell(opts.scales, 1);
    for s =1:opts.scales
        tau{s} = struct('H', var_est, 'V', var_est, 'D', var_est);
        if opts.SURE == 0 
            lambda{s} = struct('H', opts.lambda, 'V', opts.lambda, 'D', opts.lambda);
        end 
    end
    tau{opts.scales}.A = var_est;
    lambda{opts.scales}.A = opts.lambda;
    
    t_new = 1; % fista weighting initialisation
    x = 0;
    %% *** start FISTA *** 
    for iter = 1:opts.maxIter
        if opts.saveHist
            hist.r(:,:,iter) = r;
        end
        
        hist.tau_mean(iter) = var_est;
        if opts.verbose
            disp(['tau mean at iteration ', num2str(iter), ' is ', num2str(hist.tau_mean(iter))])
        end
        
        % thresholding
        C = multiscaleDecomp(r, opts.scales, opts.wavType);
        if opts.SURE
            C_thr = multiscaleSUREsoft(C, tau);
        else      
            C_thr = multiscaleComplexSoft(C, tau, lambda);
        end
        
        x_new = multiscaleRecon(C_thr, opts.wavType);

        % FISTA weighting
        x_old = x;
        t_old = t_new;
        t_new = (1+sqrt(1+4*t_old^2))/2;
        x = x_new + (t_old -1)/t_new*(x_new -x_old);

        % gradient descent
        z = dcoil - mask.*fftnc(unpad(csm.*x));       
        r = x+sum(conj(csm).*pad(ifftnc(mask.*z)), 3);     
        
        hist.x_mse(iter) = immse(unpad(x_new).*im_mask, x0.*im_mask);
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
        
        var_est = immse(x0, unpad(r));
        
        if opts.saveHist
            hist.x(:,:,iter) = x_new;
        end
          
        for s =1:opts.scales
            tau{s} = struct('H', var_est, 'V', var_est, 'D', var_est);
        end
        tau{opts.scales}.A = var_est;
    end
    
    % trim hist attributes in case a stopping criterion reached
    if opts.saveHist
        hist.r = hist.r(:,:,hist.x_mse>0); 
        hist.x = hist.x(:,:,hist.x_mse>0);
    end

    x_hat = unpad(x_new);
end

