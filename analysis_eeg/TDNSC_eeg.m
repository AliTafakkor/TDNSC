function varargout = TDNSC_eeg(q, varargin)
% TDNSC_EEG Summarizes all analyses of EEG data, including preprocessing,
% decoding, permutation tests, building RDMs, and more.
%
% INPUTS:
%   q          - String specifying the requested analysis.
%   varargin   - Variable input arguments (name-value pairs) required or
%                optional depending on the specific analysis requested.
%
% OUTPUTS:
%   varargout  - Variable output arguments corresponding to the results of the analysis.
%
% EXAMPLE USAGE:
%   args = {'subj_id', 1, 'tempgen', true, 'dsfactor', 4, 'kfold', 10, 'nperm', 50};
%   TDNSC_eeg('analysis:decode_single', args{:})
%
% NOTES:
%   - This function performs various EEG analyses, each with specific preprocessing
%     steps, decoding procedures, and permutation tests.
%   - The analyses generate RDMs (Representational Dissimilarity Matrices) as part of the workflow.
%   - Ensure that all necessary preprocessing steps are completed before running the decoding analyses. 
%   - The function supports flexible input arguments to customize analysis settings.
%
% Ali Tafakkor (atafakko@uwo.ca), University of Western Ontario

%% ------------------------------- Add Paths ------------------------------
addpath(genpath('../utils'));

%% ------------------------------ Parameters ------------------------------
set_parameters; 

%% ------------------------------- Analysis -------------------------------
switch q
    case 'analysis:decode_single'
    % Decodes epoched data for a given subject (pairwise and temporal generalization analysis)
    %
    % INPUTS:
    %   subj_id   (required) - Integer specifying subject id
    %
    %   epochsdir (optional) - String specifying the data folder name, default is 'epoched'
    %   outputdir (optional) - String specifying the output folder name, default is 'decodings'
    %   level     (optional) - Decoding 'stimuli' or 'category' of the sounds, default is 'stimuli'
    %   channels  (optional) - String or cell array of channels to be used for decoding, default is all BioSemi channels
    %   tempgen   (optional) - Logical (true/false) specifying temporal generalization analysis, default is false
    %   pca       (optional) - Logical (true/false) specifying PCA before decoding, default is false
    %   dvw       (optional) - Logical (true/false) specifying decision value weighted accuracy, default is false
    %   nperm     (optional) - Integer specifying the number of permutations for pseudo trials, default is 100
    %   kfold     (optional) - Integer specifing number of pseudo trials, default is 5
    %   dsfactor  (optional) - Integer specifying the downsampling factor, default is 1 (no down sampling)

    % Parse arguments
    vararginparse(varargin, {'subj_id'}, {'epochsdir', 'outputdir', 'level', 'channels', 'tempgen', 'pca', 'dvw', 'nperm', 'kfold', 'dsfactor'}, ...
                                         {'epoched', 'decodings', 'stimuli', p.eeg.channels, false, false, false, 100, 5, 1})

    % Temporal generalization argument
    if tempgen, method = 'temporalgen';
    else,       method = 'pairwise';
    end

    % Path to epoched data of the subject
    trialsdir = fullfile(p.path.data, sprintf('S%02d', subj_id), 'eeg', epochsdir);
    % Check if the directory exists
    if ~isfolder(trialsdir), error('Epoched data for subject S%02d could not be found under the path: %s', subj_id, trialsdir); end

    % Message
    fprintf("\n%s analysis for subject S%02d\n", upper(method), subj_id);
    fprintf("Loading data from '%s'\n", trialsdir)
    fprintf("Started at %s\n", datetime)

    % Load available trials
    trials = cell(1,p.stim.num);
    for c = 1:p.stim.num, [trials{c},time] = brainstorm_load_trials(trialsdir, c, channels); end

    % Number of trials for each condition/stimuli
    num_trials = cellfun(@(x) size(x,3), trials);
    
    % Determine labels based on the decoding level
    if     strcmp(level,'stimuli'),   labels = sprintfc('%02g',1:p.stim.num);
    elseif strcmp(level, 'category'), labels = repelem(p.stim.categories_eeg, p.stim.num_withincat);
    else, error("'%s' is not a valid option for 'level', choose either 'stimuli' or 'category'.", level)
    end

    % Concatenate and repeat labels
    labels = repelem(labels, num_trials);
    data = cat(3,trials{:});

    % Message
    fprintf("Finished at %s\n", datetime)

    % Compress data using pca
    if pca, [data,k,explained] = fl_pca(data,99.99); end
    
    % Change kfold parameter in case there are not enough trials in one condition
    kfold = min(min(num_trials),kfold);

    % Message
    if dvw, fprintf("Decoding %s (decision value weighted) started at %s\n", level, datetime);
    else fprintf("Decoding %s started at %s\n", level, datetime); end

    % SVM decoding
    D = fl_decodesvm(data(:,1:dsfactor:end,:),labels,'method',method,'numpermutation',nperm,'kfold',kfold,'dvweighted',dvw,'verbose',2);
    
    % Append parameters 
    D.level = level;
    D.time = time(1:dsfactor:end);
    D.pca = pca;
    D.dsfactor = dsfactor;

    % Message
    fprintf("Saving... \nStarted at %s\n", datetime)

    % Output directory
    output_path = fullfile(p.path.data, sprintf('S%02d', subj_id), 'eeg', outputdir);
    if ~isfolder(output_path), mkdir(output_path); end

    % Results file
    resfile = fullfile(output_path, sprintf('decodings_level-%s_tempgen-%d_pca-%d_dvw-%d_nperm-%d_folds-%d_dsfactor-%d.mat', ...
                                             level, tempgen, pca, dvw, nperm, kfold, dsfactor));
    % Save
    save(resfile, 'D');
    
    % Message
    fprintf("Results saved in '%s'. \n", resfile)
    fprintf('Finished at %s\n', datetime)
    
    case 'analysis:decode_multiple'
    % Decodes epoched data for multiple subjects given a subject list, or all subjects with available data by default
    % (pairwise and temporal generalization analysis)
    %
    % INPUTS:
    %   subjects  (optional) - List of subjects' id, default is all subjects with epoched eeg data
    %   epochsdir (optional) - String specifying the data folder name, default is 'epoched'
    %   level     (optional) - Decoding 'stimuli' or 'category' of the sounds, default is 'stimuli'
    %   channels  (optional) - String or cell array of channels to be used for decoding, default is all BioSemi channels
    %   tempgen   (optional) - Logical (true/false) specifying temporal generalization analysis, default is false
    %   pca       (optional) - Logical (true/false) specifying PCA before decoding, default is false
    %   dvw       (optional) - Logical (true/false) specifying decision value weighted accuracy, default is false
    %   nperm     (optional) - Integer specifying the number of permutations for pseudo trials, default is 100
    %   kfold     (optional) - Integer specifing number of pseudo trials, default is 5
    %   dsfactor  (optional) - Integer specifying the downsampling factor, default is 1 (no down sampling)
    
    % Parse arguments
    vararginparse(varargin, {}, {'subjects', 'epochsdir', 'outputdir', 'level', 'channels', 'tempgen', 'pca', 'dvw', 'nperm', 'kfold', 'dsfactor'}, ...
                                {p.subj.epoched_id, 'epoched', 'decodings', 'stimuli', p.eeg.channels, false, false, false, 100, 5, 1})
    
    args = {'epochsdir', epochsdir, ...
            'outputdir', outputdir, ...
            'level',     level, ...
            'channels',  channels, ...
            'tempgen',   tempgen, ...
            'pca',       pca, ...
            'dvw',       dvw, ...
            'nperm',     nperm, ...
            'kfold',     kfold, ...
            'dsfactor',  dsfactor ...
            };
    
    % Run decode_single on each subject
    for s = subjects, TDNSC_eeg('analysis:decode_single', 'subj_id', s, args{:}); end

    case 'analysis:decoding2RDM_single'
    % Converts pairwise decodings to RDMs for a single subject
    %
    % INPUTS:
    %   subj_id     (required) - Integer specifying subject id
    %
    %   decodingdir (optional) - String specifying the input folder name, default is 'decodings'
    %   rdmdir      (optional) - String specifying the output folder name, default is 'rdms'
    %   level       (optional) - String specifying decoding level 'stimuli' or 'category' of the sounds
    %   dopca       (optional) - Logical (true/false) specifying PCA before decoding
    %   dvw         (optional) - Logical (true/false) specifying decision value weighted accuracy
    %   nperm       (optional) - Integer specifying the number of permutations for pseudo trials
    %   kfold       (optional) - Integer specifing number of pseudo trials
    %   dsfactor    (optional) - Integer specifying the downsampling factor
    %   
    %   NOTES: 
    %       - Where not specified, default value for optional input arguments is '*' meaning any value will be matched.
    
    % Parse arguments
    vararginparse(varargin, {'subj_id'}, {'decodingdir', 'rdmdir', 'level', 'dopca', 'dvw', 'nperm', 'kfold', 'dsfactor'}, ...
                                         {'decodings', 'rdms', '*', '*', '*', '*', '*', '*', '*', '*'})

    % Path to epoched data of the subject
    decodingdir = fullfile(p.path.data, sprintf('S%02d', subj_id), 'eeg', decodingdir);
    % Check if the directory exists
    if ~isfolder(decodingdir), error('Decoding data for subject S%02d could not be found under the path: %s', subj_id, decodingdir); end
    
    % Path to output folder
    rdmdir = fullfile(p.path.data, sprintf('S%02d', subj_id), 'eeg', rdmdir);
    % Create output directory if it doesn't exits
    if ~isfolder(rdmdir), mkdir(rdmdir); end

    % Message
    fprintf("Generating RDMs form decodings for subject S%02d.\n", subj_id);

    % Search pattern
    filepattern = fullfile(decodingdir, sprintf('decodings_level-%s_tempgen-0_pca-%s_dvw-%s_nperm-%s_folds-%s_dsfactor-%s.mat', ...
                                         level, dopca, dvw, nperm, kfold, dsfactor));
    files = dir(filepattern);
    
    % Throw error if no file were found
    if isempty(files), error("No decoding files found for subject 'S%02d' with specified parameters!", subj_id); end

    % Convert and save files
    for i = 1:length(files)
        % Full file path
        ffpath = fullfile(files(i).folder, files(i).name);
        
        % Message
        fprintf("Decoding file found: '%s'.\n", ffpath)

        % Load decoding struct
        D = load(ffpath, 'D').('D');
        % Get labels, decodings, and time
        labels = D.condlabel;
        d = D.d;
        time = D.time;
        % Get number of timepoints and conditions
        ntp = size(d,1);
        ncond = length(labels);

        % Initialize
        RDMs = zeros(ncond, ncond, ntp);

        % Convert RDM for each time point
        for t = 1:ntp-1
            RDM = squareform(d(i,:));
            % Reorder conditions in the RDM
            RDMs(:,:,t) = TDNSC_reorder_RDM(p, RDM);
        end
        [RDMs(:,:,ntp), labels] = TDNSC_reorder_RDM(p, squareform(d(ntp,:)), labels);

        % Generate filename
        rdmfile = strsplit(files(i).name, '_');
        rdmfile{1} = 'rdms';
        rdmfile = strjoin(rdmfile, '_');
        % Save output
        save(fullfile(rdmdir, rdmfile), 'RDMs', 'labels', 'time');

        % Message
        fprintf("RDMs file saved: '%s'.\n", fullfile(rdmdir, rdmfile))
    end
    
    case 'analysis:decoding2RDM_multiple'
    % Converts pairwise decodings to RDMs for multiple subjects given a subject list,
    % or all subjects with available data by default
    %
    % INPUTS:
    %   subjects    (optional) - List of subjects' id, default is all subjects with epoched eeg data
    %   decodingdir (optional) - String specifying the input folder name, default is 'decodings'
    %   rdmdir      (optional) - String specifying the output folder name, default is 'rdms'
    %   level       (optional) - String specifying decoding level 'stimuli' or 'category' of the sounds
    %   dopca       (optional) - Logical (true/false) specifying PCA before decoding
    %   dvw         (optional) - Logical (true/false) specifying decision value weighted accuracy
    %   nperm       (optional) - Integer specifying the number of permutations for pseudo trials
    %   kfold       (optional) - Integer specifing number of pseudo trials
    %   dsfactor    (optional) - Integer specifying the downsampling factor
    %   
    %   NOTES: 
    %       - Where not specified, default value for optional input arguments is '*' meaning any value will be matched.
    
    % Parse arguments
    vararginparse(varargin, {}, {'subjects', 'decodingdir', 'rdmdir', 'level', 'dopca', 'dvw', 'nperm', 'kfold', 'dsfactor'}, ...
                                {p.subj.epoched_id, 'decodings', 'rdms', '*', '*', '*', '*', '*', '*', '*', '*'})
    
    args = {'decodingdir', decodingdir, ...
            'rdmdir',      rdmdir, ...
            'level',       level, ...
            'dopca',       dopca, ...
            'dvw',         dvw, ...
            'nperm',       nperm, ...
            'kfold',       kfold, ...
            'dsfactor',    dsfactor ...
            };
    
    % Run decoding2RDM_single on each subject
    for s = subjects, TDNSC_eeg('analysis:decoding2RDM_single', 'subj_id', s, args{:}); end

    case 'dir2targz'
    % Turns an 'eeg' sub-directory into tar file and compresses, 
    % applies on specified subjects or all subjects who have the specified directory
    %
    % INPUTS:
    %   subdirname  (required) - String specifying the name of the sub-directory
    %
    %   subjects    (optional) - List of integers specifying the ids of subjects to process,
    %                            default is all subjects with the specified sub-directory
    %
    % OUTPUTS:
    %   This case doesn't return any outputs, but generates <folder>.tar.gz files in the EEG directory of each subject
    
    % Parse arguments
    vararginparse(varargin, {'subdirname'}, {'subjects'}, {[]})

    % Automatically resolving subjects
    if isempty(subjects) 
        pattern = fullfile(p.path.data, 'S*', 'eeg', subdirname);
        dirs = unique({dir(pattern).folder});
    % Processing user provided subjects list
    else
        dirs = cellfun(@(s) fullfile(p.path.data,sprintf('S%02d',s),'eeg',subdirname), num2cell(subjects), 'UniformOutput', false);
        for d = dirs, if ~isfolder(d{1}), error("Directory '%s' does not exist!", d{1}); end, end
    end

    for d = dirs, tar(strcat(d{1}, '.tar.gz'), d{1}); end

end