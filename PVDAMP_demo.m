%% Runs the P-VDAMP algorithm only

clear all;
addpath(genpath('.'))
rng(777)

%% load demo data
load('demo_brain_data.mat')
[nx, ny, nc] = size(smaps);

%% generate sampling mask
sample_type = 'bern'; % 'bern' for bernoulli or 'pd' for poisson disc
R = 5; % acceleration factor. For this demo you can choose from R=5, 10

[mask, prob_map] = load_mask(nx, ny, 1/R, sample_type);

%% generate data
dcoil = mask.*dcoil_full;

% zero-filled
x_zf = ifftnc(dcoil./prob_map);
x_zf = sum(conj(smaps) .*x_zf, 3);

%% run P-VDAMP
opts.maxIter = 30; 
opts.maxTime = 1000;
opts.verbose = 1;
opts.scales = 4;
opts.saveHist = 1;
opts.wavType = 'db4'; 
opts.tauStop = 1;
opts.rho = 0.75;
opts.sigmaMeas = 0;

disp('Running P-VDAMP...')
[x_hat, r_hat, hist, hist_r]  = P_VDAMP(dcoil, mask, prob_map, x0, smaps, opts);

%%
x_mnsq = mean(abs(x0(:)).^2);
NMSE = @(x_mse) 10*log10(x_mse./x_mnsq);

%% plot NMSE vs time
figure('Name','NMSE vs time');
hold on;
plot(hist.timer, NMSE(hist.x_mse), '-', 'LineWidth', 2, 'Color', [1, 0, 0] , 'MarkerSize', 10, 'Marker', '.')
plot(hist_r.timer, NMSE(hist_r.x_mse), '-', 'LineWidth', 2, 'Color', [0, 0.5, 1], 'MarkerSize', 10, 'Marker', '.')
xlabel('Time (s)');
ylabel('NMSE (dB)');
legend('Biased estiamte', 'Unbiased Estimate')
set(gca, 'FontName', 'Times' );
grid on
hold off

%% view all reconstructions

figure;
subplot(2,4,1);
imshow(abs(x0),[]);
title('reference')

subplot(2,4,2);
imshow(abs(x_zf),[]);
title('Zero-filled, density compensated')

subplot(2,4,3);
imshow(abs(x_hat),[]);
title('P-VDAMP recon')

subplot(2,4,4);
imshow(abs(r_hat),[]);
title('Unbiased P-VDAMP recon')

subplot(2,4,5);
imshow(mask,[]);
title('mask')

mx_error = max(abs(x0(:) - x_zf(:)))/2;
subplot(2,4,6);
imshow(abs(x0 - x_zf),[0, mx_error]);
title('error')
set(gca, 'colormap', parula)

subplot(2,4,7);
imshow(abs(x0 - x_hat),[0, mx_error]);
title('error')
set(gca, 'colormap', parula)

subplot(2,4,8);
imshow(abs(x0 - r_hat),[0, mx_error]);
title('error')
set(gca, 'colormap', parula)

%% evidence of state evolution

it_choice = [1,2,3,4];

figure('Name', 'Evidence of State Evolution');
[ha, pos] = tight_subplot(4, numel(it_choice), 0.04, 0.11,0.05);

coarse_wavsz = 2^opts.scales; % pad to cleanly divide wavelets if needed
nx_pad = nx + (coarse_wavsz - mod(nx, coarse_wavsz))*(mod(nx, coarse_wavsz) > 0);
ny_pad = ny + (coarse_wavsz - mod(ny, coarse_wavsz))*(mod(ny, coarse_wavsz) > 0); 
pad = @(A) padarray(A, [(nx_pad-nx)/2, (ny_pad-ny)/2]);
unpad = @(A) A((nx_pad-nx)/2 + 1:end - (nx_pad-nx)/2, (ny_pad-ny)/2 + 1:end - (ny_pad-ny)/2, :);

w0 = multiscaleDecomp(pad(x0), opts.scales, opts.wavType);
I0 = pyramid(w0);

sig = 1/sqrt(2);
ex = @(x) exp(-x.^2/(2*sig^2))/(sig*sqrt(2*pi));

for iter = 1:numel(it_choice)  
    ii = it_choice(iter);
       
    axes(ha(iter));
    C = multiscaleDecomp(hist.r(:,:,ii),opts.scales, opts.wavType);
    IC = pyramid(C);
    diff = abs(IC-I0);
    
    if iter == 1
        thr = 0.3*max(diff(:));
    end
  
    imagesc(diff, [0 thr]); colormap jet;
    title(['|r_{', num2str(ii-1), '} - w_0|']);
    colorbar off;
    axis image off;
    set(gca,'FontName','times','FontSize', 12)
    
    axes(ha(iter+numel(it_choice)));
    modl = sqrt(hist.tauPy(:,:,ii));
    imagesc(modl, [0 thr]); colormap parula;
    title(['\tau_{', num2str(ii-1), '}^{1/2}']);
    colorbar off;
    axis image off;
    set(gca,'FontName','times','FontSize', 12)
    
    axes(ha(iter+2*numel(it_choice)));
    imagesc(abs(diff./modl),[0,1.5]); axis image off;
    title(['(r_{', num2str(ii-1), '} - w_0)/\tau_{', num2str(ii-1), '}^{1/2}'])

    d = IC-I0;
    normd = real(d(:))./modl(:);
    axes(ha(iter+3*numel(it_choice))); hold on;

    histogram(normd,  'Normalization', 'pdf', BinWidth=0.1); hold on;
    xlim([-3,3])
    plot(-3:0.01:3, ex(-3:0.01:3), 'r', 'LineWidth',3)
    hold off;   
end



