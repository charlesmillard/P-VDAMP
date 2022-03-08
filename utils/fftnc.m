function res = fftnc(x)

% res = fftnc(x)
% 
% orthonormal forward N-dimensional FFT.
% Does *not* do 3D fft. For dimension >2, does series of 2D ffts 
%

ndim = length(size(x));

if ndim < 3
    res = fftshift(fftn(fftshift(x)));
    res = res/sqrt(length(x(:))); % normalization
elseif ndim == 3
    res = complex(zeros(size(x)));
    for c=1:size(x, 3)
        res(:,:,c) = fftnc(x(:,:,c));
    end
elseif ndim == 4 
    res = complex(zeros(size(x)));
    for t = 1:size(x,4)
        res(:,:,:,t) = fftnc(x(:,:,:,t));
    end
else
    error('maximum 4D (3D + coils) allowed for ifftnc')
end

