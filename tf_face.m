loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\ByGrp\face\';
loadfex = '.mat'; id = dir([loadpath,'*',loadfex]);
savepath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\face\';
savefex = '.mat';
list = cellfun(@(x) x(1:end-length(loadfex)), {id.name}, 'UniformOutput', false);

strs = squeeze(split(list, '_'));
cond = unique(strs(:, 1));
resp = unique(strs(:, 2));

fs = 120;
wavename = 'cgau4';
Fc = centfrq(wavename);
totalscal = fs/2;
c = 2*Fc*totalscal;
scals = c./(1:totalscal);
freq = scal2frq(scals, wavename, 1/fs);
nfreq = length(freq);

count = 0;
for m = 1:length(cond)
    for n = 1:length(resp)
        trfs = importdata([loadpath, cond{m},'_',resp{n}, loadfex]);
        
        tf = [];
        nsubj = length(trfs.subj); nelec = length(trfs.elec); ntime = length(trfs.time);
        data = nan(nsubj, nelec, nfreq, ntime);
        for j = 1:nsubj
            for k = 1:nelec
                data(j,k,:,:) = abs(cwt(trfs.data{j}(:,k), scals, wavename));
            end
        end
        tf.data = data;
        tf.subj = trfs.subj;
        tf.elec = trfs.elec;
        tf.freq = freq;
        tf.time = trfs.time;
        
        save([savepath,cond{m},'_',resp{n},'_tf',savefex], 'tf')
        count = count + 1;
        disp(count/length(list));
    end
end

%%
clear;
close all;
resp = 'R';

cond = 'link'; % {}
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\face\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
link = avgs - mean(avgs(:,:,ibase), 3);

cond = 'back'; % {}
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\face\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
back = avgs - mean(avgs(:,:,ibase), 3);

data = [];
data.powspctrm = link-back;
data.time = tf.time;
data.freq = tf.freq;
data.dimord = 'chan_freq_time';
data.label = tf.elec;

cfg = [];
cfg.layout   = 'easycapM1.mat';
cfg.fontsize = 15;
% cfg.style    = 'straight';
cfg.marker   = 'on';
cfg.comment  = 'no';

figure
tmpcfg = cfg;
ft_multiplotTFR(tmpcfg, data);
title([cond,' ',resp]);
colorbar('FontSize',5);
axis tight;

%% Topo
close all
xlimit = [0.4, 0.5];
ylimit = [8,12];
zlimit = [-1.5,1.5];

% xlimit = [0.1, 0.15];
% ylimit = [13, 22];
% zlimit = [-1.5,1.5];

% xlimit = [0.2, 0.6];
% ylimit = [8, 12];
% zlimit = [-1.5,1.5];

mode = {'link', 'back'};
for i = 1:2
    
    cond = mode{i}; % {}
    loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\face\';
    load([loadpath, cond, '_', resp, '_tf.mat']);
    avgs = squeeze(mean(tf.data, 1));
    ibase = tf.time < 0;
    link = avgs - mean(avgs(:,:,ibase), 3);
    
    data = [];
    data.powspctrm = link;
    data.time = tf.time;
    data.freq = tf.freq;
    data.dimord = 'chan_freq_time';
    data.label = tf.elec;
    
    figure
    tmpcfg = cfg;
    tmpcfg.xlim = xlimit; % time
    tmpcfg.ylim = ylimit; % frequency
    tmpcfg.zlim = zlimit; % colorbar
    tmpcfg.comment = 'auto';
    
    ft_topoplotTFR(tmpcfg, data);
    title([cond, ' ', resp]);
    colorbar('FontSize',8);
    axis tight;
end
%%
cond = 'link';
resp = 'LR';
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
avgs_bc = avgs - mean(avgs(:,:,ibase), 3);

subplot(211)
hold on
ichan = strcmp(tf.elec, 'PO3');
PO3 = squeeze((mean(avgs_bc(ichan, 8:12, :), 2)));
plot(tf.time, PO3);
ichan = strcmp(tf.elec, 'PO4');
PO4 = squeeze((mean(avgs_bc(ichan, 8:12, :), 2)));
plot(tf.time, PO4);
plot(tf.time, PO4-PO3);
hold off
title([cond,' ',resp]);
legend({'PO3', 'PO4', 'PO4-PO3'});

cond = 'back';
resp = 'LR';
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
avgs_bc = avgs - mean(avgs(:,:,ibase), 3);

subplot(212)
hold on
ichan = strcmp(tf.elec, 'PO3');
PO3 = squeeze((mean(avgs_bc(ichan, 8:12, :), 2)));
plot(tf.time, PO3);
ichan = strcmp(tf.elec, 'PO4');
PO4 = squeeze((mean(avgs_bc(ichan, 8:12, :), 2)));
plot(tf.time, PO4);
plot(tf.time, PO4-PO3);
hold off
title([cond,' ',resp]);
legend({'PO3', 'PO4', 'PO4-PO3'});

% ft_singleplotER
% ft_multiplotER
%% diff
cond = 'link';

resp = 'L';
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
avgs_bc = avgs - mean(avgs(:,:,ibase), 3);
ltrf = avgs_bc;

resp = 'R';
loadpath = 'D:\MATLAB\SocialFace\DatAna\TRF\TF\';
load([loadpath, cond, '_', resp, '_tf.mat']);
avgs = squeeze(mean(tf.data, 1));
ibase = tf.time < 0;
avgs_bc = avgs - mean(avgs(:,:,ibase), 3);
rtrf = avgs_bc;

% % data structure
% Li = [];
% Li.powspctrm = ltrf;
% Li.time = trfs.time;
% Li.freq = trfs.freq;
% Li.dimord = 'chan_freq_time';
% Li.label = trfs.elec;
%
% Ri = []; % data structure
% Ri.powspctrm = rtrf;
% Ri.time = trfs.time;
% Ri.freq = trfs.freq;
% Ri.dimord = 'chan_freq_time';
% Ri.label = trfs.elec;

Li_Ri = []; % data structure
Li_Ri.powspctrm = rtrf - ltrf;
Li_Ri.time = tf.time;
Li_Ri.freq = tf.freq;
Li_Ri.dimord = 'chan_freq_time';
Li_Ri.label = tf.elec;

cfg = [];
cfg.layout   = 'easycapM1.mat';
cfg.fontsize = 18;
% cfg.style    = 'straight';
cfg.marker   = 'on';
cfg.comment  = 'no';

% % left
% figure
% tmpcfg = cfg;
% ft_multiplotTFR(tmpcfg, Li);
% title('left');
% colorbar('FontSize',18);
% axis tight;
%
% % right
% figure
% tmpcfg = cfg;
% ft_multiplotTFR(tmpcfg, Ri);
% title('right');
% colorbar('FontSize',18);
% axis tight;

% diff
figure
tmpcfg = cfg;
ft_multiplotTFR(tmpcfg, Li_Ri);
title(cond);
colorbar('FontSize',18);
axis tight

%% plot
xlimit = [0.25, 0.35];
ylimit = [8,12];
zlimit = [-1,1];

figure
tmpcfg = cfg;
tmpcfg.xlim     = xlimit;%time
tmpcfg.ylim     = ylimit; % frequency
tmpcfg.zlim     = zlimit; % colorbar
ft_topoplotTFR(tmpcfg, Li);
colorbar('FontSize',18);
axis tight;

figure
tmpcfg = cfg;
tmpcfg.xlim     = xlimit;%time
tmpcfg.ylim     = ylimit; % frequency
tmpcfg.zlim     = zlimit; % colorbar
ft_topoplotTFR(tmpcfg, Ri);
title('right');
colorbar('FontSize',18);
axis tight;

figure
tmpcfg = cfg;
tmpcfg.xlim     = xlimit;%time
tmpcfg.ylim     = ylimit; % frequency
tmpcfg.zlim     = zlimit; % colorbar
ft_topoplotTFR(tmpcfg, Li_Ri);
if strcmp(cond, 'lt')
    title('left-right');
elseif strcmp(cond, 'rt')
    title('right-left');
else
    title('left-right');
end
colorbar('FontSize',18);
axis tight

%%
% % plot
% figure;
% subplot(3,1,1);imagesc(time,f,squeeze(tf_lt_tlb_bs(59,:,:)));
% colorbar;
% xlabel('time(s)');
% ylabel('frequency(Hz)');
% title('lt_L');
% % caxis([-0.3,0.8]);
% axis xy;
% % set(gca,'ylim',[0 30]);
% subplot(3,1,2);imagesc(time,f,squeeze(tf_lt_trb_bs(59,:,:)));
% colorbar;
% xlabel('time(s)');
% ylabel('frequency(Hz)');
% title('lt_R');
% % caxis([-0.3,0.8]);
% % set(gca,'ylim',[0 30]);
% axis xy;
% subplot(3,1,3);imagesc(time,f,squeeze(tf_lt_tlb_bs(59,:,:)-tf_lt_trb_bs(59,:,:)));
% colorbar;
% xlabel('time(s)');
% ylabel('frequency(Hz)');
% title('Diff');
% % caxis([-0.5,0.5]);
% % set(gca,'ylim',[0 30]);
% axis xy;