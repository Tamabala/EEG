task = 'global'; fixl = 0;

switch task
    case 'global'
        loadpath = 'D:\MATLAB\NewPacman\EEGDat\ClearGlobal\';
        if fixl == 0
            savepath = 'D:\MATLAB\NewPacman\TRFDat\global\';
        else
            savepath = 'D:\MATLAB\NewPacman\TRFFix\global\';
        end
        capture = 'syncro';
    case 'local'
        loadpath = 'D:\MATLAB\NewPacman\EEGDat\ClearLocal\';
        if fixl == 0
            savepath = 'D:\MATLAB\NewPacman\TRFDat\local\';
        else
            savepath = 'D:\MATLAB\NewPacman\TRFFix\local\';
        end
        capture = 'existx';
end

loadfex = '_neat.mat'; id = dir([loadpath,'*',loadfex]);
loadlist = cellfun(@(x) x(1:end-length(loadfex)), {id.name}, 'UniformOutput', false);
savefex = '.mat'; id = dir([savepath,'*',savefex]);
savelist = cellfun(@(x) x(1:end-length(savefex)), {id.name}, 'UniformOutput', false);

id_list = setdiff(loadlist, savelist);
nsub = length(id_list);
if isempty(id_list), error('NO new data'); end

cfgs = importdata('D:\MATLAB\NewPacman\GenCfg\cfgs.mat');
cfgs = cfgs.(task);

fs = 120; time = 5;
npnts = fs*time; ncut = 0.5*fs;
want = ncut+1:npnts-ncut;
nwant = length(want);

need = cellfun(@(x) x == '0', cfgs.(capture));
stim = cfgs.lumset(need);
stim = cellfun(@(x) x(want), stim, 'UniformOutput', false);

T = cellfun(@(x) [x.top], stim, 'UniformOutput', 0); T = [T{:}];
L = cellfun(@(x) [x.left], stim, 'UniformOutput', 0); L = [L{:}];
R = cellfun(@(x) [x.right], stim, 'UniformOutput', 0); R = [R{:}];
B = cellfun(@(x) [x.bottom], stim, 'UniformOutput', 0); B = [B{:}];
TL = T.*L; LB = L.*B; BR = B.*R; RT = R.*T;
TLB = T.*L.*B; TRB = T.*R.*B;
stim = [T; L; R; B; TL; LB; BR; RT; TLB; TRB];

stimName = {'T', 'L', 'R', 'B', 'TL', 'LB', 'BR', 'RT', 'TLB', 'TRB'};
nstim = length(stimName);
stim = mat2cell(stim, nstim, nwant*ones(1, sum(need)));

% stim = [];
% for i = 1:nstim
%     var = stimName{i};
%     if length(var) ~= 1 
%         eval([var,' = 1;'])
%         for j = 1:length(var)
%             eval([var, '=', var, '.*', var(j),';'])
%         end
%     end
%     stim = cellfun(@(x) [x, var], stim, 'UniformOutput', 0); %B = [B{:}]';
% %     eval(['stim = [stim, ', var, '];']) 
% end
% stim = cellfun(@(x) cell2mat(struct2cell(x)), stim, 'UniformOutput', false);

mode = {'lt','rt','nt'};
for i = 1:nsub
    subj = id_list{i};
    EEG = importdata([loadpath, subj, loadfex]);
    
    cfg = [];
    cfg.resamplefs = fs;
    EEG = ft_resampledata(cfg, EEG);
    
    cfg = [];
    cfg.polyremoval = 'yes';
    cfg.polyorder = 1;
    EEG = ft_preprocessing(cfg, EEG);
    
    nlabel = length(EEG.label);
    ntrial = length(EEG.trialinfo);
    resp = EEG.trial;
    resp = cellfun(@(x) x(:,want), resp, 'UniformOutput', false);% same
    
    if fixl == 0
        stii = stim; stii = [stii{:}];
        stii = zscore(stii, 1, 2);
        stii = mat2cell(stii', nwant*ones(1, ntrial), nstim);
        
        resi = resp; resi = [resi{:}];
        resi = zscore(resi, 1, 2);
        resi = mat2cell(resi', nwant*ones(1, ntrial), nlabel);
        
        tmin = -200; tmax = 1000; dirTRF = 1; Lambdas = 10.^(-5:5);
        cv = mTRFcrossval(stii,resi,fs,dirTRF,tmin,tmax,Lambdas,'verbose',0,'zeropad',1,'fast',1);
        [Rmax,LambdaIndx] = max(mean(cv.r,[1,3],'omitnan')); lambda = Lambdas(LambdaIndx);
    elseif fixl == 1
        lambda = 1;
    end
    
    data = cell(length(mode),1);
    for j = 1:length(mode)
        if strcmp(mode{j}, 'lt')
            marker = 1;
        elseif strcmp(mode{j}, 'rt')
            marker = 2;
        elseif strcmp(mode{j}, 'nt')
            marker = 3;
        end
        idx = find(EEG.trialinfo == marker);
        
        stii = stim(idx); stii = [stii{:}];
        stii = zscore(stii, 1, 2);
        stii = mat2cell(stii', nwant*ones(1, length(idx)), nstim);
        
        resi = resp(idx); resi = [resi{:}];
        resi = zscore(resi, 1, 2);
        resi = mat2cell(resi', nwant*ones(1, length(idx)), nlabel);
       
        model = mTRFtrain(stii,resi,fs,dirTRF,tmin,tmax,lambda,'zeropad',1,'method','ridge','type','multi','verbose',0);
        data{j} = model.w;
    end
    
    trf = [];
    trf.mode = mode;
    trf.data = data;
    trf.stim = stimName;
    trf.time = -0.2:1/120:1;
    trf.elec = EEG.label;
    save([savepath, subj, savefex], 'trf');
end
