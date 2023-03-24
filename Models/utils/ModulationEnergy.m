clear; clc;
addpath('./spectrotemporal-synthesis-v2')
addpath('./spectrotemporal-synthesis-v2/Sound_Texture_Synthesis_Toolbox');
load('parameters_PLoSBio2018.mat', 'P');

stimset_path = './stimset';
directory = dir(stimset_path);
k = 1;
stimuli = [];
stimuli_name = {};
duration_sec = 1;
for i = 3:length(directory)
    cat_path = fullfile(directory(i).folder, directory(i).name);
    stim_path = dir(cat_path);
    for j = 3:length(stim_path)
        stim = fullfile(cat_path, stim_path(j).name);

        [wav, sr] = audioread(stim);
        wav = wav(1:duration_sec*sr,1);
        if sr ~= P.audio_sr
            wav = resample(wav, P.audio_sr, sr);
        end

        stimuli(k, :) = wav;
        stimuli_name{k} = stim_path(j).name(1:end-4);
        k = k+1;
    end
end

%% Modulation energy of cochleograms (stage 2)
ME = nan(size(stimuli,1), length(P.spec_mod_rates), length(P.temp_mod_rates)-1);
Cochs = nan(size(stimuli,1), 400, 217);
parfor i = 1:size(stimuli,1)
    [coch, f, t, R_orig] = cochleagram_wrapper(squeeze(stimuli(i,:)), duration_sec, P);
    ME(i,:,:) = modulation_energy(coch, P, 1:size(coch,1));
    Cochs(i,:,:) = coch;
    fprintf('%d / %d compeleted \n',i,size(stimuli,1));
end
%%
save('ME.mat', 'ME')
save('Cochleagrams.mat', 'Cochs')
%%
load('ME.mat', 'ME'); load("Cochleagrams.mat","Cochs");
%%
ME_norm = zeros(size(ME));
parfor i = 1:size(ME,1)
    ME_norm(i,:,:) = ME(i,:,:)/max(ME(i,:,:),[],'all');
end
%% plot modulation energy examples
figure;
selected = [1, 21, 41, 61];
for i = 1:4
    subplot(2,4,i)
    s = selected(i);
    
    f = P.f; t = P.t;
    C = squeeze(Cochs(s,:,:));
    dims = size(C);
    n_t = dims(1);
    n_f = dims(2);
    imagesc(flipud(C'));
    % y-axis ticks
    yi = round( linspace( 1, n_f, 5 ) );
    f_flip = flip(f);
    set(gca, 'YTick', yi, 'YTickLabel', num2cellstr( f_flip(yi), '%.0f' ) );
    clear yi f_flip;
    
    % x-axis ticks
    xi = round( linspace(1, n_t, 5) );
    set(gca, 'XTick', xi, 'XTickLabel',  num2cellstr( t(xi) , '%.1f') );
    clear xi
    
    % x-label and y-label
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    box off;
    title(stimuli_name{s});

    subplot(2,4,i+4)
    pos_temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);
    contour(pos_temp_mod_rates, P.spec_mod_rates, squeeze(ME_norm(s,:,:)),...
        'LineWidth', 2)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'YTick', [0.25, 1, 4], 'YTickLabel', [0.25, 1, 4]);
    set(gca, 'XTick', [1, 4, 16, 64], 'XTickLabel', [1, 4, 16, 64])
    %set(gca, 'FontSize', 12);
    xlabel({'Temporal modulation'; 'rate (Hz)'});
    ylabel({'Spectral modulation'; 'scale (cycles/octave)'});
    title(stimuli_name{s});
end

%% mean modulation energy for categories
categories = {'animals', 'objects', 'people', 'scenes'};
categories_ind = reshape(1:80, 20, 4)';

for i = 1:4
    subplot(2,2,i)
    pos_temp_mod_rates = P.temp_mod_rates(P.temp_mod_rates>0);
    contour(pos_temp_mod_rates, P.spec_mod_rates,...
            squeeze(mean(ME_norm(categories_ind(i,:),:,:), 1)),...
            'LineWidth', 2)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'YTick', [0.25, 1, 4], 'YTickLabel', [0.25, 1, 4]);
    set(gca, 'XTick', [1, 4, 16, 64], 'XTickLabel', [1, 4, 16, 64])
    %set(gca, 'FontSize', 12);
    xlabel({'Temporal modulation'; 'rate (Hz)'});
    ylabel({'Spectral modulation'; 'scale (cycles/octave)'});
    title(sprintf('Average modulation energy of %s', categories{i}));
end
