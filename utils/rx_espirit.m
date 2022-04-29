function sens = rx_espirit(calib, imsize, kernel, eig_thresh, mask_thresh)
%
% sens = tx_espirit(calib, imsize, kernel, eig_thresh, mask_thresh)
%
% Inputs:
%           calib:      [NKx, NKy, z, NRx] complex Rx k-space
%           imsize:     [Nx, Ny] image dimensions for output
%           kernel:     [kx, ky] kernel size
%                       (default [5,5])
%           eig_thresh: scalar threshold for eigenvalues 
%                       if < 1, interpreted as s ≤ s(1)*eig_thresh
%                       if > 1, interpreted as s_r, r ≤ eig_thresh
%                       (default 0.02)
%           mask_thresh:scalar threshold for sensitivity masking
%                       (default is 0.99)
%
% Outputs:
%           sens:       [Nx, Ny, NRx] Rx sensitivity maps
%
% Implementation of receive sensitivity mapping based on 
% ESPIRiT (MRM 2014)
%
% MChiew (mark.chiew@ndcn.ox.ac.uk)

% Default params
if nargin < 5
    mask_thresh =   0.99;
end
if nargin < 4
    eig_thresh  =   0.02;
end
if nargin < 3
    kernel      =   [5,5];
end

if length(size(calib)) == 3
    calib_temp(:,:,1, :) = calib;
    calib = calib_temp;
    clear calib_temp
end

% Get dimensions
dims    =   size(calib(:,:,1,:));
Nz      =   size(calib,3);

% Initialise outputs
sens    =   zeros([imsize Nz dims(4)]);

for z = 1:Nz
    % Generate calibration matrix
    H       =   fold_rx(Hankel_fwd(calib(:,:,z,:), kernel, dims));
    
    % Get left singular (kernel x coil) vectors, above singular value threshold
    [U,S,~] =   svd(H, 'econ');
    if eig_thresh < 1
        U       =   U(:,diag(S) > S(1)*eig_thresh);
    else
        U       =   U(:,1:round(eig_thresh));
    end

    % Get zero-padded space x coil x component kernel images
    U       =   reshape(U,kernel(1),kernel(2),dims(4),[]);
    U       =   padarray(U,[ceil((imsize-kernel)/2) 0 0],'pre');
    U       =   padarray(U,[floor((imsize-kernel)/2) 0 0],'post');
    tmp     =   fftshift(fftshift(ifft2(fftshift(fftshift(U,1),2)),1),2);

    % Perform SVD on coil x component matrices, voxelwise
    % Keep the first component
    m       =   zeros(imsize);
    for i = 1:imsize(1)
        for j = 1:imsize(2)
            [U,S,~]     =   svd(squeeze(tmp(i,j,:,:)), 'econ');        
            if S(1)*prod(imsize)/sqrt(prod(kernel)) > mask_thresh
                sens(i,j,z,:) =   U(:,1);
            end
        end
    end
end
% Rotate phase relative to first coil
sens    =   sens.*exp(-1j*angle(sens(:,:,:,1)));
sens = squeeze(sens);
end

function x = fold_rx(x, dims)
    if nargin < 2
        dims(1) = size(x,1)*size(x,4);
        dims(2) = size(x,2)*size(x,3);
    end
    x = reshape(permute(x, [1, 4, 2, 3]), dims);
end

function h = Hankel_fwd(x, kernel, dims)
    if numel(dims) == 3
        dims(4) = 1;
    end
    Nx  = dims(1);
    Ny  = dims(2);
    N1  = dims(3);   
    N2  = dims(4);
        
    h   = zeros(prod(kernel), prod([Nx,Ny]-kernel+1), N1, N2);

    idx = 0;
    for kx = 1:Nx - kernel(1) + 1
        for ky = 1:Ny - kernel(2) + 1
            idx = idx + 1;
            h(:, idx, :, :) = reshape(x(kx:kx+kernel(1)-1, ky:ky+kernel(2)-1,:,:),prod(kernel),1,N1,N2);
        end
    end
end
