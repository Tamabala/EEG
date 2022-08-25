% loadpath = 'D:\MATLAB\NewPacman\EEGDat\RawLocal\';
% savepath = 'D:\MATLAB\NewPacman\EEGDat\ClearLocal\';

loadpath = 'D:\MATLAB\NewPacman\EEGDat\RawGlobal\';
savepath = 'D:\MATLAB\NewPacman\EEGDat\ClearGlobal\';

loadfex = '.vhdr'; id = dir([loadpath,'*',loadfex]);
loadlist = cellfun(@(x) x(1:end-length(loadfex)), {id.name}, 'UniformOutput', false);
savefex = '_neat.mat'; id = dir([savepath,'*',savefex]);
savelist = cellfun(@(x) x(1:end-length(savefex)), {id.name}, 'UniformOutput', false);

id_list = setdiff(loadlist, savelist);
nsub = length(id_list);
if isempty(id_list), error('NO new data'); end


subj = id_list{1};

%% ft_definetrial : Define trials
cfg = [];% reset
cfg.dataset = [loadpath,subj,'.vhdr'];
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = {'S  1', 'S  2', 'S  3'};
cfg.trialdef.prestim    = 0;
cfg.trialdef.poststim   = 5;
cfg = ft_definetrial(cfg);

% ft_databrowser(cfg);

% ft_preprocessing : Filter
cfg.bpfilter =  'yes';% lp/hp/bs/dft/median,'yes' or 'no'
cfg.bpfreq = [0.1,40];
cfg.bpfiltord = 2;
cfg.padding = 15;
cfg.padtype = 'data';
eeg = ft_preprocessing(cfg);% or (cfg, data)

% ft_channelrepair : Channel repair
% cfg = [];
% cfg.method = 'spline' ;%  'weighted'(defauit), 'average', 'spline', 'slap' or 'nan'.
% cfg.missingchannel = {'Oz'};
% cfg.layout = 'easycapM1.mat';% ft_prepare_layout : output cfg.elec
% eeg = ft_channelrepair(cfg, eeg);

% ft_preprocessing : Rereference
cfg = [];
cfg.channel     = {'all', '-IO'}; % chantype only for NeuroOmega, default=[1:31,33:62].
cfg.reref       = 'yes'; % Rereference
cfg.refchannel  = {'all', '-HEOG', '-VEOG', '-IO'};
cfg.implicitref = 'FCz'; % the implicit (non-recorded) reference channel is added to the data representation,default = [].
cfg.demean      = 'no';% whether to apply baseline correction (default = 'no')
eeg = ft_preprocessing(cfg, eeg);

save([savepath, subj,'.mat'], 'eeg');

% ica
cfg = [];
cfg.method = 'fastica';
eeg_ica = ft_componentanalysis(cfg, eeg);
save([savepath, subj, '_ica.mat'], 'eeg_ica');

ncomp = length(eeg_ica.label);
for i = 1:ceil(ncomp/20)
    figure;
    cfg = [];
    cfg.component = 1+20*(i-1):min(20*i,ncomp);% specify the component(s) that should be plotted
    cfg.layout    = 'easycapM1.mat'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, eeg_ica);
end

rm_com = input('input the component number to remove:');
cfg = [];
cfg.component = rm_com; % to be removed component(s)
cfg.demean    = 'no';
eeg_neat      = ft_rejectcomponent(cfg, eeg_ica, eeg);%

ft_databrowser([], eeg_neat);

% cfg = [];
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.2, 0];
% eeg_neat = ft_preprocessing(cfg, eeg_neat);
save([savepath, subj, '_neat.mat'], 'eeg_neat');
