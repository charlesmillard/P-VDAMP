%% Runs P-VDAMP, Optimal FISTA, Shared FISTA, SURE-IT and A-FISTA

clear all;
addpath(genpath('.'))

% if you want to run A-FISTA you'll need to download their code from
% https://gitlab.com/cfmm/datasets/auto-lambda-for-cs-recons
% and add it to the path:
addpath(genpath('your_own_path/auto-lambda-for-cs-recons')); 

do_afista = exist('auto_l1_recon', 'dir'); 
%% load data
disp('loading data...')
load('demo_brain_data.mat') % load 2D data in format (nx, ny, ncoils). Please feel free to try your own!

[nx, ny, nc] = size(dcoil_full);
smaps = rx_espirit(dcoil_full, [nx, ny], [5,5], 0.02, 0);
x0 = sum(conj(smaps).*ifftnc(dcoil_full), 3);
%% generate sampling mask
sample_type = 'bern'; % 'bern' for bernoulli or 'pd' for poisson disc
R = 5; % acceleration factor. For this demo choose from R=5, 10

[mask, prob_map] = load_mask(nx, ny, 1/R, sample_type);

% generate data
dcoil = mask.*dcoil_full;

%% run P-VDAMP
opts.maxIter = 200; 
opts.maxTime = 1000;
opts.verbose = 1;
opts.scales = 4;
opts.saveHist = 0;
opts.wavType = 'db4'; 
opts.tauStop = 1;
opts.rho = 0.75;
opts.sigmaMeas = 0;

disp('Running P-VDAMP...')
[x_hat{1}, x_hat{2}, hist{1}, hist{2}]  = P_VDAMP(dcoil, mask, prob_map, x0, smaps, opts);
alg_names{1} = 'P-VDAMP';
alg_names{2} = 'Unbiased P-VDAMP';

%% Optimal FISTA
opts.SURE = 0;
load(['lambda_us', num2str(R), sample_type]);
opts.lambda = lambda_optimal;
disp('Running optimal FISTA...')
[x_hat{3}, hist{3}] = FISTA_coils(dcoil, smaps, mask, x0, opts); 
alg_names{3} = 'Optimal FISTA';

%% Shared FISTA
opts.lambda = 5.9e4;
disp('Running shared FISTA...')
[x_hat{4}, hist{4}] = FISTA_coils(dcoil, smaps, mask, x0, opts); 
alg_names{4} = 'Shared FISTA';

%% SURE-IT
opts.SURE = 1;
disp('Running SURE-IT...')
[x_hat{5}, hist{5}]  = FISTA_coils(dcoil, smaps, mask, x0, opts); 
alg_names{5} = 'SURE-IT';

%% A-FISTA
disp('Running auto...')
if do_afista
    [x_hat{6}, hist{6}] = complex_mFISTA2_adapted(dcoil, mask, smaps, x0, opts, 'Lambda','multiple','Levels',4);
    alg_names{6} = 'A-FISTA';
else
    warning('You do not have A-FISTA code on path, so not running A-FISTA')
end

%%
x_zf = ifftnc(dcoil./prob_map);
x_zf = sum(conj(smaps) .*x_zf, 3);

x_mnsq = mean(abs(x0(:)).^2);
NMSE = @(x_mse) 10*log10(x_mse./x_mnsq);

%% show NMSE vs iteration
line_colour =  [0.1, 0.6, 0; 1, 0.8, 0; 1, 0, 0; 0, 0.5, 1; 0,0,0; 0.8,0,0.8];
line_type = {'-', '-', '--','--', '--', '--'};

figure;
hold on;
for ii = 1:numel(x_hat)
    plot(hist{ii}.timer, NMSE(hist{ii}.x_mse), '-', 'LineWidth', 0.1, 'Color', line_colour(ii,:), 'MarkerSize', 4, 'Marker', '.')
end
xlabel('Time (s)');
ylabel('NMSE (dB)');

legend(alg_names, 'NumColumns', 3, 'FontSize', 20);
set(gca, 'FontName', 'Times' );
grid on
hold off

%% show all reconstructions
figure;
im_tile = abs(x0);
error_tile = mask;
t = 'reference';
for ii = 1:numel(x_hat)
    im_tile = [im_tile, abs(x_hat{ii})];
    error_tile = [error_tile, abs(x0 - x_hat{ii})];
    t = string(t) +  ' || ' + string(alg_names{ii});
end

subplot(2,1,1)
imshow(im_tile, []);
title(t)
mx_error = max(abs(x0(:) - x_zf(:)))/2;
subplot(2,1,2)
imshow(error_tile, [0, mx_error]);
set(gca, 'colormap', parula)

%% display quantitative quality scores
im_mask = (abs(x0) > max(abs(x0(:)))/20);

disp('*** NMSE ***')
disp(['zero-filled: ' + string(nmse(x_zf.*im_mask, x0.*im_mask))])
for ii = 1:numel(x_hat)
    disp([string(alg_names{ii}) +  ': ' + string(nmse(x_hat{ii}.*im_mask, x0.*im_mask))])
end

disp('*** SSIM ***')
mx = max(abs(x0(:)));
disp(['zero-filled: ' + string(ssim(abs(x_zf).*im_mask/mx, abs(x0).*im_mask/mx))])
for ii = 1:numel(x_hat)
    disp([string(alg_names{ii}) +  ': ' + string(ssim(abs(x_hat{ii}).*im_mask/mx, abs(x0).*im_mask/mx))])
end

disp('*** HFEN ***')
disp(['zero-filled: ' + string(hfen(x_zf.*im_mask, x0.*im_mask))])
for ii = 1:numel(x_hat)
    disp([string(alg_names{ii}) +  ': ' + string(hfen(x_hat{ii}.*im_mask, x0.*im_mask))])
end

disp('*** Time ***')
for ii = 1:numel(x_hat)
    disp([string(alg_names{ii}) +  ': ' + string(hist{ii}.timer(end))])
end
