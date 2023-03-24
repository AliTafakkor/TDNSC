function varargout = TDNSC_auditorymodels(q, varargin)
%% ---------------- Add Paths ----------------
addpath('../utils');
addpath('./nsltools');
%% --------------- Parameters ----------------
p.saveflag = false;
p.savepath = './figs';
p.datapath = '/Windows/EEGDataset/data';
p.stimulipath = '~/Workspace/Thesis/analysis/TDNSC/stimuli';

p.n_stim = 80;
p.n_cat = 4;
p.cat_names = {'Animals', 'People', 'Objects', 'Scenes'};
p.n_withincat = 20;

%% -------------------------------------------
switch q
    case 'load:stimuli'
        vararginparse(varargin, {'spath', 'fpath'});
        
        if exist('spath')
            p.stimulipath = spath;
        end
        if exist('fpath')
            p.savepath = fpath;
        end

        categories = ["animals", "people", "objects", "scenes"];
        ite = 1;
        for c = 1:length(categories)
            category = categories(c);
            files = dir(fullfile(p.stimulipath, category, '*.wav'));
            for i=1:size(files,1)
                stim(ite).name = files(i).name(1:end-4);
                stim(ite).category = category;
                stim(ite).ID = ite;
                stim(ite).path = fullfile(files(i).folder, files(i).name);
                % Read file
                [y, fs] = audioread(stim(ite).path);
                stim(ite).sig = y;
                stim(ite).fs = fs;
                stim(ite).length = size(y,1)/fs;
                stim(ite).energy = sum(y(:,1).^2);
                stim(ite).maxamp = max(y(:,1));
    
                ite = ite + 1; 
            end
        end
    
        % Load oddball
        stim(ite).name = 'noise';
        stim(ite).category = 'oddball';
        stim(ite).ID = ite;
        stim(ite).path = fullfile(p.stimulipath,'noise.wav');
        [y, fs] = audioread(stim(ite).path);
        stim(ite).sig = y;
        stim(ite).fs = fs;
        stim(ite).length = size(y,1)/fs;
        stim(ite).energy = sum(y(:,1).^2);
        stim(ite).maxamp = max(y(:,1));
    
        varargout{1} = stim; 
    case ''
end
