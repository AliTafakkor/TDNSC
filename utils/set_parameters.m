%% --------------- Paths ----------------
p.path.root = '~/Workspace/TDNSC';									        % path to the root of repository
p.path.data = '~/Workspace/TDNSC/data';							            % path to the dataset
p.path.stim = fullfile(p.path.root, 'stimuli');						        % 
p.path.fig = './fig';												        %
p.path.mat = './mat';												        %

%% -------------- Stimuli ---------------
p.stim.num = 80;													        % number of stimuli
p.stim.num_cat = 4;													        % number of categories
p.stim.num_withincat = 20;											        % number of stimuli within category
p.stim.categories = {'animals', 'people', 'objects', 'scenes'};		        % category names
p.stim.categories_eeg = {'animals', 'objects', 'scenes', 'people'};	        % category order in eeg experiment (note the difference)

%% ---------------- Subjects -------------
p.subj.all_dir = dir(fullfile(p.path.data,'S*'));                           % direcotry structure of all subjects
p.subj.all_id = cellfun(@(x) str2num(x(2:end)), ...                         % id list of all subjects 
                        {p.subj.all_dir.name});

ind = cellfun(@(x) isfolder(fullfile(p.path.data,x,'eeg','epoched')), ...
    {p.subj.all_dir.name});
p.subj.epoched_dir = p.subj.all_dir(ind);
p.subj.epoched_id = p.subj.all_id(ind);

%% ---------------- EEG -----------------
p.eeg.fs = 512;														        % sampling frequency (Hz)
p.eeg.channels = {'A','B'};                                                 % channels for analysis 
                                                                            % ('A*' and 'B*' are the eeg channels in BioSemi 64-channels)

p.eeg.trial_start = -200; 											        % start of the trial relative to onset (ms)
p.eeg.trial_end = 2000;												        % end of the trial relative to onset (ms)
p.eeg.time = p.eeg.trial_start:(1000/p.eeg.fs):p.eeg.trial_end;             % time vector: corresponding time to each sample relative to onset (ms)

p.eeg.trial_trancate = 1300; 									            % truncate time (ms)
% p.eeg.ttime_ind = floor((p.eeg.trial_trancate-p.eeg.trial_start)/1000*p.eeg.fs)+1;       % 
% 
% p.eeg.time = p.eeg.time(1:p.ttime_ind);
% p.eeg.n_timepoints = length(p.time);
% 
% 
% p.eeg.alpha = 0.01;
% p.eeg.alpha_corrected = p.alpha/p.n_timepoints;                     % Bonferroni corrected
% p.eeg.ids_subj = [1, 2, 3, 4, 7, 11, 14, 18, 19, 23, 24, 26:37];
% 
% 
% 
% p.eeg.num_subj = length(p.eeg.ids_subj);
% p.eeg.dvw = false;
% p.eeg.pca = false;
% 
% p.ow = true; % over-write results?
% 
% p.sound_ch = 1;
% p.trig_ch = 2;

p.colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
