p.n_stim = 80;
p.n_cat = 4;
p.categories = {'animals', 'people', 'objects', 'scenes'};
p.cat_names = {'Animals', 'People', 'Objects', 'Scenes'};
p.n_withincat = 20;
p.fs = 512;
p.first_timepoint = -200;
p.last_timepoint = 2000;
p.ttime = 1300; % truncate time
p.ttime_ind = floor((p.ttime-p.first_timepoint)/1000*p.fs)+1;
p.time = p.first_timepoint:(1/p.fs)*1000:p.last_timepoint;
p.time = p.time(1:p.ttime_ind);
p.n_timepoints = length(p.time);
p.alpha = 0.01;
p.alpha_corrected = p.alpha/p.n_timepoints;%0.00005; % Bonferroni corrected

% Paths
p.figpath = './fig';
p.matpath = './mat';
p.data_path = '/Windows/EEGDataset/data';
p.stimulipath = '~/Workspace/Thesis/analysis/TDNSC/stimuli';

p.subj_ids = [1, 2, 3, 4, 7, 11, 14, 18, 19, 23, 24, 26:37];
p.n_subjects = length(p.subj_ids);
p.dvw = false;
p.pca = false;

p.ow = true; % over-write results?

p.sound_ch = 1;
p.trig_ch = 2;

p.colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
