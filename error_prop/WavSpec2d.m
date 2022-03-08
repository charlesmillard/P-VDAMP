function [psdH, psdV, psdD, psdA] = WavSpec2d(scales, wavType, nx, ny)
    % compute 1D wavelet spectra
    specX = WavSpec(scales, nx, wavType);
    specY = WavSpec(scales, ny, wavType);
    
    % 2D spectra from 1D
    psdH = zeros(nx, ny, scales);
    psdV = zeros(nx, ny, scales);
    psdD = zeros(nx, ny, scales);
    for s = 1:scales
        psdH(:,:,s) = specX(:,s,2) * specY(:,s,1)';         
        psdV(:,:,s) = specX(:,s,1) * specY(:,s,2)'; 
        psdD(:,:,s) = specX(:,s,2) * specY(:,s,2)'; 
    end
    psdA = specX(:,scales,1) * specY(:,scales,1)'; 
end

function spec = WavSpec(scales, sample, wavType)
    % Compute power spectrum of wavelet filters at multiple scales
    
    % coarse_wavsz = 2^scales;
    % sample_pad = (coarse_wavsz - mod(sample, coarse_wavsz))*(mod(sample, coarse_wavsz) > 0);
    
    [filt(:,1) , filt(:,2)] = wfilters(wavType);
 
    spec = zeros(sample, scales, 2);
    
    L = zeros(sample, 1);
    H = zeros(sample, 1);
    
    L(1:size(filt, 1)) = filt(:, 1);
    H(1:size(filt, 1)) = filt(:, 2);
    
    spec(:,1,1) = abs(fft(L)).^2;
    spec(:,1,2) = abs(fft(H)).^2;
    
    sigLen = sample;
    numBlock = 1;
    
    for s=2:scales
        sigLen = sigLen / 2;
        numBlock = numBlock *2;
        
        % Multiply with the previous low-pass spectrum with aliasing
        spec_block = reshape(spec(1:numBlock:end,1,1)*ones(1,numBlock), [],1);
        if numel(spec_block) ~= sample
            % trim to approximate spectrum if wavelet does not exactly fit
            % - could be improved
            p = (numel(spec_block) - sample)/2;
            spec_block = spec_block(p+1:end-p);
        end
        spec(:, s, 1) = spec(:,s-1,1) .* spec_block;
        
        spec_block = reshape(spec(1:numBlock:end,1,2)*ones(1,numBlock), [],1);
        if numel(spec_block) ~= sample
            p = (numel(spec_block) - sample)/2;
            spec_block = spec_block(p+1:end-p);
        end     
        spec(:, s, 2) = spec(:,s-1,1) .*spec_block;
    end
    
    spec = fftshift(spec, 1) / sample;
end


    