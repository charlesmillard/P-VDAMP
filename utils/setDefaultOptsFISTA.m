function opts = setDefaultOptsFISTA(opts)

        % maximum allowed number of iterations
    if ~(isfield(opts,'maxIter') && ~isempty(opts.maxIter))
        opts.maxIter = 20;
    end   
    
    % maximum allowed time
    if ~(isfield(opts,'maxTime') && ~isempty(opts.maxTime))
        opts.maxTime = 50;
    end
    
    % stop when aliasing model tau increases?
    if ~(isfield(opts,'tauStop') && ~isempty(opts.tauStop)) 
        opts.tauStop = 1;
    end
    
    % is verbose
    if ~(isfield(opts,'verbose') && ~isempty(opts.verbose))
        opts.verbose = 0;
    end
    
    % number of wavelet scales
    if ~(isfield(opts,'scales') && ~isempty(opts.scales))
        opts.scales = 4;
    end

    % save history
    if ~(isfield(opts,'saveHist') && ~isempty(opts.saveHist))
        opts.saveHist = 0;
    end
    
    % wavelet family
    if ~(isfield(opts,'wavType') && ~isempty(opts.wavType)) 
        opts.wavType = 'db4';
    end
    
    % use SURE to choose thresholds?
    if ~(isfield(opts,'SURE') && ~isempty(opts.SURE)) 
        opts.SURE = 1;
    end

    % set lambda
    if ~(isfield(opts,'lambda') && ~isempty(opts.lambda)) 
        opts.lambda = 1;
    end
    
end