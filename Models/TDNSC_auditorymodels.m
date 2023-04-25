function varargout = TDNSC_auditorymodels(q, varargin)
%% --------------------- Add Paths ---------------------
addpath('../utils');
addpath('./utils')
addpath('./nsltools');

%% -------------------- Parameters ---------------------
set_parameters;
p.savefig = false;

%% -----------------------------------------------------
switch q
    case 'load:stimuli'
        vararginparse(varargin, {}, {'spath', 'fpath'});
        
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

    case 'load:spec_cochleograms'
        load(fullfile(p.matpath, 'spectrotemporal-features.mat'), 'cochs', 'f', 't');
        varargout = {cochs, t, f};
    case 'load:NSLfeatures'
        load(fullfile(p.matpath, 'NSL-features.mat'));
        varargout = {cochs, time, freq, stm, rate, scale};

    case 'align_stimuli'
        vararginparse(varargin, {'stim'}, {});

        for i = 1:p.n_stim+1
            fs = stim(i).fs;
            trig_ch = stim(i).sig(:,p.trig_ch);
            onset_ind = find(trig_ch,1);
            aud = stim(i).sig(onset_ind:end,p.sound_ch);
            % pad to 1 sec
            aud = [aud; zeros(onset_ind-1,1)];
            stim(i).sig = aud;
        end

        varargout{1} = stim;

    
    case 'spectrotemporal:coch'
        vararginparse(varargin, {'stim'}, {});

        addpath('./spectrotemporal-synthesis-v2')
        addpath('./spectrotemporal-synthesis-v2/Sound_Texture_Synthesis_Toolbox');
        load('./spectrotemporal-synthesis-v2/parameters_PLoSBio2018.mat', 'P');

        ME = nan(p.n_stim, length(P.spec_mod_rates), length(P.temp_mod_rates)-1);
        cochs = nan(p.n_stim, 400, 217);

        duration_sec = 1;
        parfor i = 1:p.n_stim-1
            wav = stim(i).sig(:,p.sound_ch);
            fs = stim(i).fs;
            if fs ~= P.audio_sr
                wav = resample(wav, P.audio_sr, fs);
            end

            [coch, f, t, R_orig] = cochleagram_wrapper(wav, duration_sec, P);
            ME(i,:,:) = modulation_energy(coch, P, 1:size(coch,1));
            cochs(i,:,:) = coch;
            fprintf('%d / %d compeleted \n',i,p.n_stim);
        end

        i = p.n_stim;
        wav = stim(i).sig(:,p.sound_ch);
        fs = stim(i).fs;
        if fs ~= P.audio_sr
            wav = resample(wav, P.audio_sr, fs);
        end

        [coch, f, t, R_orig] = cochleagram_wrapper(wav, duration_sec, P);
        ME(i,:,:) = modulation_energy(coch, P, 1:size(coch,1));
        cochs(i,:,:) = coch;
        fprintf('%d / %d compeleted \n',i,p.n_stim);

        save(fullfile(p.matpath, 'spectrotemporal-features.mat'), ...
            'cochs', 'f', 't', 'ME');
        varargout = {cochs, ME};
    
    case 'nsl:coch+stm'
        vararginparse(varargin, {'stim'}, {'octshft'})
        
        if ~exist('octshft', 'var'), octshft = 0; end;
    
        addpath('nsltools/')
        loadload;
        paras(4) = octshft;							% octave shift
        paras(1:2) = paras(1:2) / 2 ^(octshft+1);	% frame size and step
        
        rv=2.^(2:1:8);% rate vector (Hz)
        sv=2.^(-2:.5:3); %scale vector (cycle/octave)

        path = fullfile(p.matpath, 'nsl-corfiles');
        if ~isfolder(path), mkdir(path); end;

        cochs = nan(p.n_stim, 250, 128);
        stm = nan(p.n_stim, 250, 128, length(sv), 2*length(rv));

        for i = 1:p.n_stim
            x = stim(i).sig(:,p.sound_ch);
            fs = stim(i).fs;
            if fs ~= 16000
                x = resample(x, 16000, fs);
            end
            
            x=unitseq(x); % UNITSEQ.m sequence normalization to N(0, 1)
            
            aud = wav2aud(x, paras);
            fname=fullfile(path, sprintf('stim%d.cor',i));
            cor=aud2cor(aud, paras, rv, sv, fname);

            cochs(i,:,:) = aud;
            stm(i,:,:,:,:) = permute(cor, [3 4 1 2]);
        end

        N = 250; %time X freq 
        time=(1:N)' * paras(1);
        freq = 440 * 2.^ ([(0:128)-31]/24);	% for 16 kHz 
        %last channel is lost, 
        % because of linear differentiation between adjacent channels
        % see COCHFILT.txt
        scale=sv;
        rate=rv;
               
%         figure;
%         imagesc(flipud(aud'), [0 0.5]);
%         colormap("parula");
        save(fullfile(p.matpath, 'NSL-features.mat'), ...
            'cochs', 'freq', 'time', 'stm', 'scale', 'rate', ...
            '-v7.3');
        varargout{1} = cochs;
        varargout{2} = stm;

    case 'rdm:features2RDM'
        vararginparse(varargin, {'features', 'featname'}, {'method'})
        % features are expected to be a matrix of size: n_stim x any size
        % and dimension for feature
        % 

        % Default for optional argins
        if ~exist('method')
            method = 'euclidean';
        end

        features = reshape(features, p.n_stim, []);
        rdm = zeros(p.n_stim);
        
        switch method
            case 'pearson'
                rdm = 1-corrcoef(features');
            case 'euclidean'
                for i = 1:p.n_stim
                    for j= i+1:p.n_stim
                        rdm(i,j) = sum((features(i,:,:)-features(j,:,:)).^2, 'all');
                    end
                end
                rdm = sqrt(rdm);
                rdm = rdm + rdm';
        end
        
        fname = sprintf('rdm_%s_%c.mat', featname, method(1));
        save(fullfile(p.matpath, fname), 'rdm');
        varargout{1} = rdm;

    case 'rdm:timewindowed-coch2RDM' 
        vararginparse(varargin, {'coch', 't'}, {'method'})
        % features are expected to be a matrix of size: n_stim x any size
        % and dimension for feature
        % 

        % Default for optional argins
        if ~exist('method')
            method = 'euclidean';
        end

        RDMs = zeros(p.n_timepoints,p.n_stim,p.n_stim);
        %t = t *1000;

        parfor tp = 1:p.n_timepoints
            if p.time(tp)<=0
                win_ind = 1:length(t); % use whole cochleagram
            elseif p.time(tp)<1000
                [~,ind] = min(abs(t-p.time(tp)));
                win_ind = 1:ind;
            else
                win_ind = 1:length(t); % use whole cochleagram
            end
            
            coch_win = coch(:,win_ind,:);
            features = reshape(coch_win, p.n_stim, []);
                
            switch method
                case 'pearson'
                    rdm = 1-corrcoef(features');
                case 'euclidean'
                    rdm = zeros(p.n_stim);
                    for i = 1:p.n_stim
                        for j= i+1:p.n_stim
                            rdm(i,j) = sum((features(i,:,:)-features(j,:,:)).^2, 'all');
                        end
                    end
                    rdm = sqrt(rdm);
                    rdm = rdm + rdm';
            end

            RDMs(tp,:,:) = rdm;
        end
        fname = sprintf('rdm_timewindowed_cochs_%c.mat', method(1));
        save(fullfile(p.matpath, fname), 'RDMs');
        varargout{1} = RDMs;

    case 'rdm:timewindowed-stm2RDM' % TODO
        vararginparse(varargin, {'stm', 't'}, {'method'})
        % features are expected to be a matrix of size: n_stim x any size
        % and dimension for feature
        % 

        % Default for optional argins
        if ~exist('method')
            method = 'euclidean';
        end

        RDMs = zeros(p.n_timepoints,p.n_stim,p.n_stim);

        for tp = 1:p.n_timepoints
            if p.time(tp)<=0
                win_ind = 1:length(t); % use whole cochleagram
            elseif p.time(tp)<1000
                [~,ind] = min(abs(t-p.time(tp)));
                win_ind = 1:ind;
            else
                win_ind = 1:length(t); % use whole cochleagram
            end
            
            stm_win = stm(:,win_ind,:,:,:);
            features = reshape(stm_win, p.n_stim, []);
                
            switch method
                case 'pearson'
                    rdm = 1-corrcoef(features');
                case 'euclidean'
                    rdm = zeros(p.n_stim);
                    for i = 1:p.n_stim
                        for j= i+1:p.n_stim
                            rdm(i,j) = sum((features(i,:,:)-features(j,:,:)).^2, 'all');
                        end
                    end
                    rdm = sqrt(rdm);
                    rdm = rdm + rdm';
                otherwise
                   rdm = zeros(p.n_stim);
            end

            RDMs(tp,:,:) = rdm;
        end
        fname = sprintf('rdm_timewindowed_stm_%c.mat', method(1));
        save(fullfile(p.matpath, fname), 'RDMs');
        varargout{1} = RDMs;

    case 'plot:stimuli_overview'
        vararginparse(varargin, {'stim'}, {});

        f = figure;
        f.WindowState = "maximized";
        
        fs = stim(1).fs;
        t = 0:1/fs:1-1/fs;

        for i=1:80
            y = stim(i).sig(:,p.sound_ch);
            h = subplot(8,10,i);
            hold on
            plot(t*1000,y,'k');
            
            if(size(stim(i).sig,2)>1)
                plot(0,0,'.', 'MarkerSize', 20, 'Color','r');
                plot(t(find(stim(i).sig(:,p.trig_ch),1))*1000,0,'.', 'MarkerSize', 20, 'Color', 'g');
            end

            axis([0 1000 -1 1])
            if i == 80
                xlabel('Time (ms)');
                h.YAxis.Visible = 'off';
            else
                axis off;
            end
        end
        
        %pause(2)
        %saveas(f,'./overview.svg')
        %close(f);

    case 'plot:cochleograms' % TODO
        vararginparse(varargin, {'cochs'}, {})
        
        f=figure;
        f.WindowState = 'maximized';

        for i = 1:p.n_stim
            subplot(8,10,i)
            C = squeeze(cochs(i,:,:));
            imagesc(flipud(C'));
            colormap("parula");
            if(i~=p.n_stim)
                axis off;
            else
                xlabel('time')
                ylabel('freq (Hz)');
            end
        end
    otherwise
        error('Undefined command!')
end
