function res = ifftnc(x)
% orthonormal centered N-dimensional ifft

ndim = length(size(x));

if ndim < 3
    res = ifftshift(ifftn(ifftshift(x)));
    res = sqrt(length(x(:)))*res;
elseif ndim == 3
    res = complex(zeros(size(x)));
    for c=1:size(x, 3)
        res(:,:,c) = ifftnc(x(:,:,c));
    end
elseif ndim == 4
    res = complex(zeros(size(x)));
    for c=1:size(x, 4)
        res(:,:,:,c) = ifftnc(x(:,:,:,c));
    end
else 
    error('maximum 4D (3D + coils) allowed for ifftnc')
end

