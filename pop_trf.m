loadpath = 'D:\MATLAB\a_experiment\Pacman\EEGDat\neat\global\';
savepath = 'D:\MATLAB\a_experiment\Pacman\TRFDat\global\';
fex = '.set'; id = dir([loadpath,'*',fex]);
id_list = cellfun(@(x) x(1:end-length(fex)), {id.name}, 'UniformOutput', false);

load('D:\MATLAB\a_experiment\Pacman\GenCfg\cfgs.mat')
varName = fieldnames(cfgs.global.lumset{1});
nvar = length(varName); npnt = 600;
idx  = cellfun(@(x) x == '0', cfgs.global.capture);
stim = cfgs.global.lumset(idx);
stim = cellfun(@(x) cell2mat(struct2cell(x)), stim, 'UniformOutput', false);
stim = cellfun(@(x) x(:, 61:540), stim, 'UniformOutput', false); stim = [stim{:}]; 
stim = zscore(stim, 1, 2);
stim = mat2cell(stim', 480*ones(1, sum(idx)), nvar);

subj = id_list{2};
EEG = pop_loadset([subj, fex], loadpath);
EEG = pop_resample(EEG, 120);
nlabel = EEG.nbchan;
resp = EEG.data(:,end-(npnt-1):end,:);
resp = resp(:,61:540,idx);
if size(resp,3) ~= 120
    disp('Trial-dim mismatch')
    return;
end
resp = zscore(reshape(resp, 64, []), 1, 2);
resp = mat2cell(resp', 480*ones(1, sum(idx)), nlabel);

%
fs = EEG.srate; tmin = -200; tmax = 1000; dirTRF = 1; Lambdas = 10.^(-4:4);
cv = mTRFcrossval(stim,resp,fs,dirTRF,tmin,tmax,Lambdas,'verbose',0,'zeropad',1,'fast',1);
[Rmax,LambdaIndx] = max(mean(cv.r,[1,3],'omitnan')); lambda = Lambdas(LambdaIndx);
model = mTRFtrain(stim,resp,fs,dirTRF,tmin,tmax,lambda,'zeropad',1,'method','ridge','type','multi','verbose',0);

numFig = ceil(nlabel/25);
time = tmin:1000/fs:tmax;
for i = 1:numFig
    pool = 1+25*(i-1):min(25*i, nlabel);
    f1 = figure;
    f1.Name = [num2str(pool(1)),' - ', num2str(pool(end))];
    f1.Position = [40,100,800,600];
    for j = pool
        subplot(5,5,j-25*(i-1));
        hold on;
        for k = 3%1:nvar
            plot(time, squeeze(model.w(k,:,j)));
        end
%         ylim([-0.1, 0.1]);
        title([EEG.chanlocs(j).labels])
    end
end
hold off;
