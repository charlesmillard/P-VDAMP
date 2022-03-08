function xi = WavXi(smaps, wavType, scales)
    % compute scalar used to approximate wavelet spectra
    nc = size(smaps, 3);
    xi = cell(scales, 1); 
    
    [filt1 , filt2] = wfilters(wavType);

    dwtmode('per','nodisp')
    for c = 1:nc
        smap_coarse = conj(smaps(:, :, c));
        for s = 1:scales
            [c_wav,s_wav] = wavedec2(smap_coarse, 1, abs(filt1).^2, abs(filt2).^2);
            [xi{s}.H(:,:,c),xi{s}.V(:,:,c),xi{s}.D(:,:,c)] = detcoef2('all',c_wav,s_wav,1);
            smap_coarse = dwt2(smap_coarse./2, wavType, 'mode', 'per');
        end
        xi{scales}.A(:,:,c) = appcoef2(c_wav,s_wav,abs(filt1).^2,abs(filt2).^2);
    end
end

    
    

