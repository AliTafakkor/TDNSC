function [trials, time] = brainstorm_load_trials(trialsdir, condition, channelname, verbose)
% BRAINSTORM_LOAD_TRIALS Loads all trials for a specified subject and condition, 
% excludes bad trials, and loads data for the specified channels only.
%
% INPUTS:
%   trialsdir   (required) - String specifying the path to the directory containing the trial .mat files.
%   condition   (required) - Integer specifying the condition for which to load the trials.
%   channelname (required) - Cell array of strings specifying the names of channels to load,
%                            e.g., {'MEG', 'EEG', 'A', 'B'}.
%   verbose     (optional) - Integer specifying verbose level, default 
%
% OUTPUTS:
%   trials      - 3D array (channels x time x trials) containing the loaded trials data after excluding bad trials.
%   time        - Vector specifying the time points associated with the trials.
%
% EXAMPLE USAGE:
%   [trials, time] = brainstorm_load_trials('path/to/trials', 1, {'MEG', 'EEG'});
%
% NOTES:
%   - This function requires a brainstormstudy.mat file in the trials directory 
%     and a channel.mat file either in the same folder or under @default_study in the subject folder.
%   - It automatically excludes bad trials from the loaded data.
%   - Only the specified channels are imported for analysis.
%
%   Ali Tafakkor (atafakko@uwo.ca), University of Western Ontario

% Set default value of verbose to zero (no messages)
if nargin==3, verbose = 0; end

% Get subject directory
[subjdir, ~, ~] = fileparts(trialsdir);

% Find channel file
try 
    channelfile = dir(fullfile(trialsdir,'channel*.mat'));
    Channel = load(fullfile(channelfile.folder,channelfile.name));
    % Message
    if verbose, disp('subject-specific channel found!'); end
catch
    channelfile = dir(fullfile(subjdir,'@default_study','channel*.mat'));
    Channel = load(fullfile(subjdir,'@default_study',channelfile.name));
end 
Channel = Channel.Channel;

% Find channelnames
for ch = 1:length(channelname)
    pos = strfind({Channel.Name},channelname{ch});
    ndx{ch} = find(~cellfun(@isempty,pos));
end
channelndx = unique([ndx{:}],'stable');

% Find trial files
ftrials = dir(fullfile(trialsdir,sprintf('*data_%d_trial*',condition)));

% Remove bad trials from trials found
bststudy = load(fullfile(trialsdir,'brainstormstudy.mat'));
ndx = [];
for i = 1:length(bststudy.BadTrials)
    pos = strfind({ftrials.name},bststudy.BadTrials{i});
    bad_ind = find(~cellfun(@isempty,pos));
    if ~isempty(bad_ind)
        ndx = [ndx, bad_ind];
    end
end
ftrials(ndx) = [];

% Read trials
for i = 1:length(ftrials)
    data = load(fullfile(trialsdir,ftrials(i).name));
    trials{i} = data.F(channelndx,:);
    time = data.Time;

    % get import time from history
    [nrow,ncol] = size(data.History);
    for j = 1:nrow
        if strcmp(data.History(j,2),'import_time')
            t = sscanf(data.History{j,3}(6:end),['%f']);
            trial_import_time(i) = t;
        end
    end
    
end

% Reorder trials based on import time
try
[~,I] = sort(trial_import_time);
trial_import_time = trial_import_time(I);
trials = trials(I);
end

% Convert to matrix
trials = single(cat(3,trials{:}));
