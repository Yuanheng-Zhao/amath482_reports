%{ 
AMATH 482
Professor: Craig Gin
HW#4: Music Genre Classification  
Jonathan Zhao
%}

clear; close all; clc
%% Song/Music by artist Ketsa (Instrumental)
urls = ["https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_12_-_Green_Man.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_11_-_Slow_Vibing.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_07_-_The_Stork.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_02_-_Seeing_You_Again.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequency/Ketsa_-_08_-_Multiverse.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_13_-_Mission_Ready.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_10_-_Memories_Renewed.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_09_-_Life_Illusion.mp3"];

Ketsa_clips = [];  % 40 clips from the songs/music by Ketsa
dur = 5;  % 5 seconds clip
for i = 1:length(urls)
    [y, Fs] = webread(urls(i));
    for j = 1:5
        [y1, y2] = extract_clip2(y, Fs, 30 + j*20, dur);
        y_ave = (y1 + y2) / 2;
        Ketsa_clips = [Ketsa_clips y_ave];
    end
end
%%
Ketsa_dat = wavelet_spectrogram(Ketsa_clips, Fs);

%% Song/Music by artist Lately Kind of Yeah (Rock)
urls = ["https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_15_-_Heart_Feel.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_16_-_Johnny_Mathis.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_10_-_Tubescreamer.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_12_-_Kingdom_Come.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_02_-_Geist_IX.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_14_-_Goner.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_13_-_Two-Face.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_11_-_Long_Live.mp3"];
    
Lky_clips = [];  % 40 clips from the songs/music by Ketsa
dur = 5;  % 5 seconds clip
for i = 1:length(urls)
    [y, Fs] = webread(urls(i));
    for j = 1:5
        [y1, y2] = extract_clip2(y, Fs, 30 + j*20, dur);
        y_ave = (y1 + y2) / 2;
        Lky_clips = [Lky_clips y_ave];
    end
end

%%
Lky_dat = wavelet_spectrogram(Lky_clips, Fs);

%% Song/Music by artist Dee Yan-Key (Latin)
urls = ["https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_06_-_Beguine_Sailing_Trip.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_08_-_Montuno_Evening_Mood.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_07_-_Tango_Good_Air.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_10_-_Fast_Bossa_Nova_Falling_Stars.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_11_-_Rumba_Windblown.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_03_-_Paso_Doble_Boat_Ahoy.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_05_-_Habanera_By_the_Sea.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_01_-_Salsa_After-Work_Party.mp3"];

Dyk_clips = [];  % 40 clips from the songs/music by Ketsa
dur = 5;  % 5 seconds clip
for i = 1:length(urls)
    [y, Fs] = webread(urls(i));
    for j = 1:5
        [y1, y2] = extract_clip2(y, Fs, 30 + j*30, dur);
        y_ave = (y1 + y2) / 2;
        Dyk_clips = [Dyk_clips y_ave];
    end
end

%%
Dyk_dat = wavelet_spectrogram(Dyk_clips, Fs);

%%
[U,S,V,w,vband1,vband2,vband3] = genre_trainer(Ketsa_dat,Lky_dat,Dyk_dat,10);

%%
plot(vband1,1,'ro');hold on
plot(vband2,2,'go');
plot(vband3,3,'bo');

sort1 = sort(vband1);
sort2 = sort(vband2);
sort3 = sort(vband3);

threshold1 = get_threshold(sort2,sort3);  % 7.0377e+03
threshold2 = get_threshold(sort3,sort1);  % 1.4806e+04
%%
figure
subplot(2,2,1)
histogram(sort1,20); hold on, plot([threshold2 threshold2], [0 10], 'r')
title('band1 - Instrumental')
subplot(2,2,2)
histogram(sort3,20); hold on, plot([threshold2 threshold2], [0 10], 'r')
title('band3 - Rock')
subplot(2,2,3)
histogram(sort3,20); hold on, plot([threshold1 threshold1], [0 10], 'r')
title('band3 - Rock')
subplot(2,2,4)
histogram(sort2,20); hold on, plot([threshold1 threshold1], [0 10], 'r')
title('band2 - Latin')


%% Test the classifier
urls = ["https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_03_-_Dusty_Hills.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_01_-_Dreaming_Days.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_05_-_Crescents.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_14_-_Goner.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_01_-_Specific_Ocean.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Decoder_Magazine/Lately_Kind_of_Yeah/Poindexter/Lately_Kind_of_Yeah_-_03_-_The_Little_Red_School_House.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_12_-_Samba_Pop_Pancake_Tuesday.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_04_-_Samba_Rose_Monday.mp3";
    "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Dee_Yan-Key/latin_summer/Dee_Yan-Key_-_09_-_Cha-Cha-Cha_Beach_Life.mp3"];

test_clips = [];  % 40 clips from the songs/music by Ketsa
dur = 5;  % 5 seconds clip
slice_time = [30 30 30 20 20 10 35 40 40];
for i = 1:length(urls)
    [y, Fs] = webread(urls(i));
    for j = 1:5
        [y1, y2] = extract_clip2(y, Fs, 30 + j*slice_time(i), dur);
        y_ave = (y1 + y2) / 2;
        test_clips = [test_clips y_ave];
    end
end

%%
TestNum = size(test_clips,2)
Test_wave = wavelet_spectrogram(test_clips, Fs); % wavelet transformation 
TestMat = U'*Test_wave;  % PCA projection
pval = w'*TestMat;  % LDA projection
%%
hidden_labels = repelem([1 2 3], 15); % Answer to our test data

ResVec = zeros(1,TestNum); % Answer by the classifier model
for i = 1:TestNum
    if pval(i) > threshold2
        ResVec(i) = 1; % band1
    elseif pval(i) < threshold1
        ResVec(i) = 2; % band2
    else
        ResVec(i) = 3; % band3
    end
end

err_num = 0;
for i = 1:TestNum
    if (ResVec(i) ~= hidden_labels(i))
        err_num = err_num + 1;
    end
end
% Report the testing result
disp('Number of mistakes')
err_num
disp('Rate of success');
sucRate = 1-err_num/TestNum

%%
function musData = wavelet_spectrogram(clips_data, Fs)
    [n,nclips] = size(clips_data);
    L = n / Fs;  % 5 secs
    t2 = linspace(0,L,n+1); 
    t = t2(1:n);
    tslide = 0:0.1:L; 
    musData = zeros(length(tslide) * n, nclips); 
    for i = 1:nclips  % for each clip, we create a corresponding spectrogram
        y = clips_data(:,i);
        v = y';
    %     k = ((2*pi)/L)*[0:(n-1)/2 -(n-1)/2:-1];  % length of n is odd
    %     ks = fftshift(k);
        vgt_spec = zeros(length(tslide), n);
        for j = 1:length(tslide)
            g = exp(-100*(t-tslide(j)).^2); 
            vg = g.* v; 
            vgt = fft(vg);
            vgt_spec(j,:) = fftshift(abs(vgt));
        end
        musData(:, i) = reshape(vgt_spec,length(tslide)*n,1);
    end
end


% Load a clip of a song with a specifed url, start time and duration in
% seconds. Assume the input time in the range of the duration of the song.
function [y1, y2, Fs] = extract_clip(url, start_time, duration)  

    % y is an N x 2 array: the first column represents the left channel 
    % and the second column represents the right channel.
    [y, Fs] = webread(url);
    mus_length = length(y) / Fs;  % duration of the song in seconds
    if (start_time > mus_length) || (start_time + duration > mus_length)
        disp(url)
        error("ERROR.\nThe specified range of the clip exceeds the length of the song")
    end
    start_point = start_time * Fs;
    end_point = start_point + duration * Fs;
    y = y(start_point:end_point, :);
    y1 = y(:,1); y2 = y(:,2);
    
end


% 2nd ver: Load a clip based on input signal.
function [y1, y2] = extract_clip2(y, Fs, start_time, duration)  

    % y is an N x 2 array: the first column represents the left channel 
    % and the second column represents the right channel.
    mus_length = length(y) / Fs;  % duration of the song in seconds
    if (start_time > mus_length) || (start_time + duration > mus_length)
        error("ERROR: The specified range of the clip exceeds the length of the song")
    end
    start_point = start_time * Fs;
    end_point = start_point + duration * Fs;
    y = y(start_point:end_point, :);
    y1 = y(:,1); y2 = y(:,2);
    
end


function [U,S,V,w,vband1,vband2,vband3] = genre_trainer(band1,band2,band3,feature)
    n1 = size(band1,2);
    n2 = size(band2,2);
    n3 = size(band3,2);
    [U,S,V] = svd([band1 band2 band3],'econ');
    genres = S*V'; % projection onto principal components
    U = U(:,1:feature);
    band11 = genres(1:feature,1:n1);
    band22 = genres(1:feature,n1+1:n1+n2);
    band33 = genres(1:feature,n1+n2+1:n1+n2+n3);
    
    mband1 = mean(band11,2);
    mband2 = mean(band22,2);
    mband3 = mean(band33,2);
    mtotal = mean([mband1 mband2 mband3],2);
    
    Sw = 0; % within class variances
    for k = 1:n1
        Sw = Sw + (band11(:,k)-mband1)*(band11(:,k)-mband1)';
    end
    for k = 1:n2
        Sw = Sw + (band22(:,k)-mband2)*(band22(:,k)-mband2)';
    end
    for k = 1:n3
        Sw = Sw + (band33(:,k)-mband3)*(band33(:,k)-mband3)';
    end
    Sb = 40*(mband1-mtotal)*(mband1-mtotal)' + 40*(mband2-mtotal)*(mband2-mtotal)'...
        + 40*(mband3-mtotal)*(mband3-mtotal)';

    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    vband1 = w'*band11; 
    vband2 = w'*band22;
    vband3 = w'*band33;
    
end


% method helps to determine a threshold between two adjacent categories
% Assume values of s1 are just lower than those of s2
function threshold = get_threshold(s1, s2)
    t1 = length(s2);
    t2 = 1;
    while s2(t1) > s1(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    threshold = (s2(t1)+s1(t2))/2;
end


