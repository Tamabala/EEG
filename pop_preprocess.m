%% load and eeglab
loadpath = 'D:\MATLAB\a_experiment\Pacman\EEGDat\new_g\';
savepath = 'D:\MATLAB\a_experiment\Pacman\EEGDat\neat\global\';
id = dir(strcat(loadpath,'*.vhdr'));
id_list = cell(1,length(id));
for i = 1:length(id)
    id_list{i} = id(i).name(1:(length(id(i).name)-5));
end

%% prepro
subj = id_list{2};

EEG = pop_loadbv(loadpath, [subj,'.vhdr']);
EEG = pop_chanedit(EEG, 'lookup','D:\MATLAB\R2019a\toolbox\eeglab\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 30 );

layout = readlocs('D:\MATLAB\R2019a\toolbox\eeglab\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
Oz = layout(strcmp({layout(:).labels}, 'Oz'));
Oz.type = [];
EEG = pop_interp(EEG, Oz, 'spherical');

FCz = layout(strcmp({layout(:).labels}, 'FCz'));
arg = {EEG.nbchan+1, '', FCz.labels, FCz.sph_radius, FCz.sph_theta, FCz.sph_phi,...
    FCz.theta, FCz.radius, FCz.X, FCz.Y, FCz.Z, 'FCz', [], 0, FCz.sph_theta_besa, FCz.sph_phi_besa};
EEG = pop_chanedit(EEG, 'add', arg);

EEG = pop_select(EEG, 'nochannel',{'IO'});
EEG = pop_reref(EEG, [], 'refloc', EEG.chaninfo.nodatchans);
% pop_saveset(EEG, 'filename', [subj,'.set'], 'filepath', savepath);

% EEG_tmp = pop_epoch(EEG ,{'S  8'}, [-0.2,1] );
% figure; pop_spectopo(EEG_tmp, 1, [-200,998], 'EEG', 'percent', 50, 'freq', [6 10 22], 'freqrange',[1,30])
% repair = input('input the bad channel to repair:');
% EEG = pop_interp(EEG, repair, 'spherical');

%% ica
EEG_ica = pop_runica(EEG, 'icatype', 'fastica', 'stabilization','on','approach', 'symm');% ,,'approach', 'symm'£¬ 'numOfIC', 50,
% pop_saveset(EEG_ica, 'filename', [subj,'_ica.set'], 'filepath', savepath);

% inspect by eye
pop_selectcomps(EEG_ica,1:size(EEG_ica.icaweights,1));

% % iclabel
% EEG_label = pop_iclabel( EEG_ica, 'default');
% % pop_selectcomps
% pop_viewprops(EEG_label,0,1:size(EEG_label.icaact,1),{'freqrange', [1 30]});
% eyeic = input('input the bad ic to remove:');
% 
% EEG_iflag = pop_icflag( EEG_label, [NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
% a_label = find(EEG_iflag.reject.gcompreject == 1);
% 
% % sasic
% cfg = [];
% % cfg.ADJUST.enable = 1;
% cfg.SNR.enable = 1;
% cfg.focalcomp.enable = 1;
% % cfg.trialfoc.enable = 1;
% cfg.autocorr.enable = 1;
% cfg.opts.noplot = 1;
% EEG_sasic = eeg_SASICA(EEG_getic, cfg);
% a_sasic = find(EEG_sasic.reject.gcompreject == 1);
% 
% % mara
% [a_mara, info]= MARA(EEG_getic);
% 
% autoic = unique([a_label', a_sasic, a_mara]);
% 
% check = setxor(autoic,eyeic);
% EEG_epoch = pop_epoch( EEG_ica ,{'S  8'}, [-0.2,1] );
% for i = 1:length(check)
%     pop_prop(EEG_epoch,0,check(i),[],{'freqrange', [1 30]});% second input 0 for ica
% end
% exic = input('input the bad ic to exclude:');
% badic = setxor(unique([autoic,eyeic]),exic);
% if length(badic) > 10
%     fprintf(['warning: ',num2str(length(badic)),' ic will be removed from data\n'])
% end

badic = [1,2,3,5,32,50];
EEG_neat = pop_subcomp(EEG_ica, badic);
eegplot( EEG_ica.data, 'srate', EEG_ica.srate, 'title', 'Pre: K / Post: R',...
    'limits', [EEG_ica.xmin, EEG_ica.xmax], 'data2', EEG_neat.data);
EEG_neat = pop_epoch(EEG_neat, {'S  8'}, [-0.3, 5]);
EEG_neat = pop_rmbase(EEG_neat, [-300, 0]);
pop_saveset(EEG_neat, 'filename', [subj,'_neat.set'], 'filepath', savepath);
