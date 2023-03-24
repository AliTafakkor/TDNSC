clear; clc;
% Use Right channel as trigger wavefrom
% 1: Left
sound_ch = 1;
% 2: Right
trig_ch = 2;


%% parameters
stimulipath = '.';

%% Load mp3 files
categories = ["animals", "people", "objects", "scenes"];
ite = 1;
for c = 1:length(categories)
    category = categories(c);
    files = dir(fullfile(stimulipath, category, '*.mp3'));
    for i=1:size(files,1)
        stim(ite).name = files(i).name(1:end-4);
        stim(ite).category = category;
        stim(ite).ID = ite;
        stim(ite).path = fullfile(files(i).folder, files(i).name);
        % Read file
        [y, fs] = audioread(stim(ite).path);
        stim(ite).fs = fs;
        
        % Clip to 1s
        stim(ite).sig = y(1:fs,:);

        stim(ite).length = size(stim(ite).sig,1)/fs;
        stim(ite).energy = sum(y(:,sound_ch).^2);
        stim(ite).maxamp = max(y(:,sound_ch));

        ite = ite + 1; 
    end
end

% Load oddball
stim(ite).name = 'noise';
stim(ite).category = 'oddball';
stim(ite).ID = ite;
stim(ite).path = fullfile(stimulipath,'noise.wav');
[y, fs] = audioread(stim(ite).path);
stim(ite).sig = y;
stim(ite).fs = fs;
stim(ite).length = size(y,1)/fs;
stim(ite).energy = sum(y(:,1).^2);
stim(ite).maxamp = max(y(:,1));

%%
target_energy = 100;
for i = 1:81
    comp_coeff = sqrt(target_energy/stim(i).energy);
    stim(i).sig(:,sound_ch) = stim(i).sig(:,sound_ch)*comp_coeff;
    stim(i).energy = sum(stim(i).sig(:,sound_ch).^2);
    stim(i).maxamp = max(stim(i).sig(:,sound_ch));
    %audiowrite([stim(i).path(1:end-3), 'wav'], stim(i).sig, fs);
end

%%
absolute_threshold = 0.005;
t = 0:1/fs:1-1/fs;
for i=1:81
    y = stim(i).sig(:,sound_ch);

    % find onset based on amplitude thresholding
    amp_threshold = absolute_threshold; %0.03*stim(i).maxamp;
    onset_ind = find(amp_threshold<abs(y),1);
    amp_onset = t(onset_ind);
    stim(i).sig(:,trig_ch) = ones(fs,1);
    stim(i).sig(1:onset_ind,trig_ch) = 0;

    f = figure;
    hold on
    plot(t*1000,y,'k');
    plot(t*1000,stim(i).sig(:,trig_ch), 'r')
    plot([0],[0],'.', 'MarkerSize', 20, 'Color','r');
    plot([amp_onset]*1000,[0],'.', 'MarkerSize', 20, 'Color', 'g');

    axis([0 1000 -1 1])

    legend({'Sound Waveform (Left channel)',...
            'Trigger waveform (Right channel)',...
            'Play onset', 'Amplitude onset'},...
           'Location', 'northeastoutside')
    xlabel('Time (ms)');
    h = gca;
    h.YAxis.Visible = 'off';
    title(sprintf('stimulus #%d',i))
    saveas(f,sprintf('./figs/%d.jpg',i))
    close(f);
end

%%
f = figure;
f.WindowState = "maximized";

for i=1:80
    y = stim(i).sig(:,sound_ch);
    h = subplot(8,10,i);
    hold on
    plot(t*1000,y,'k');
    
    plot(0,0,'.', 'MarkerSize', 20, 'Color','r');
    plot(t(find(stim(i).sig(:,trig_ch),1))*1000,0,'.', 'MarkerSize', 20, 'Color', 'g');

    axis([0 1000 -1 1])
    if i == 80
        xlabel('Time (ms)');
        h.YAxis.Visible = 'off';
    else
        axis off;
    end
end

pause(2)
saveas(f,'./overview.svg')
close(f);

%% write to .wav files
for i = 1:length(stim)
    audiowrite([stim(i).path(1:end-3), 'wav'], stim(i).sig, fs);
end