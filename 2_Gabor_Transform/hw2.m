%{ 
AMATH 482
Professor: Craig Gin
HW#1: Gabor Transform
Jonathan Zhao
%}
%% Part I - Starter Code
clear; close all; clc
load handel  % Fs = the sampling frequency in Hz
             % y = the audio signal amplitude as a single column vector
v = y';
plot((1:length(v))/Fs, v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

% p8 = audioplayer(v,Fs);
% playblocking(p8);   

%% Part I - Exploring different window widths
L = length(v)/Fs;  % 8.925 seconds
n = length(v);
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k=((2*pi)/L)*[0:(n-1)/2 -(n-1)/2:-1];  % length of n is odd
ks=fftshift(k);

a = [1 100 1000];  % parameters for Gaussian Filter
tslide = 0:0.1:L;  
vgt_spec = zeros(length(tslide), n);
vgst_spec = zeros(length(tslide), n);  % spectrogram w/ small window width
vglt_spec = zeros(length(tslide), n);  % spectrogram w/ large window width
for j = 1:length(tslide)
    g = exp(-a(2)*(t-tslide(j)).^2); 
    g_narroww = exp(-a(3)*(t-tslide(j)).^2);  % narrower window
    g_widew = exp(-a(1)*(t-tslide(j)).^2);  % wider window
    
    vg = g.* v; 
    vgt = fft(vg); 
    vgs = g_narroww .* v;
    vgst = fft(vgs);
    vgl = g_widew .* v;
    vglt = fft(vgl);
    
    vgt_spec(j,:) = fftshift(abs(vgt));
    vgst_spec(j,:) = fftshift(abs(vgst));
    vglt_spec(j,:) = fftshift(abs(vglt));
end

% Spectrograms with different window widths
figure(1)
subplot(3,1,1)
pcolor(tslide,ks/(2*pi),vgt_spec.'), 
shading interp
title('Normal Window: a=100')  % The 'a' here is not the width, but the parameter in 
                % the Gaussian filter to change the width.
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,2)
pcolor(tslide,ks/(2*pi),vgst_spec.'), 
shading interp 
title('Narrow Window: a=1000')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,3)
pcolor(tslide,ks/(2*pi),vglt_spec.'), 
shading interp 
title('Wide Window: a=1')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

%% Part I - Oversampling & Undersampling
tslide_over = 0:0.01:L;  % using very small translations of the Gabor window
vgtspec_over = zeros(length(tslide_over), n);

for j = 1:length(tslide_over)
    g = exp(-100*(t-tslide_over(j)).^2); 
    vg = g.* v;
    vgt = fft(vg);
    vgtspec_over(j,:) = fftshift(abs(vgt));
end

tslide_under = 0:1:L;  %  using very course/large translations
vgtspec_under = zeros(length(tslide_under), n);

for j = 1:length(tslide_under)
    g = exp(-100*(t-tslide_under(j)).^2); 
    vg = g.* v;
    vgt = fft(vg);
    vgtspec_under(j,:) = fftshift(abs(vgt));
end

% Spectrograms with different time translations
figure(2)
subplot(3,1,1)
pcolor(tslide,ks/(2*pi),vgt_spec.'), 
shading interp 
title('length(tslide)=91')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,2)
pcolor(tslide_over,ks/(2*pi),vgtspec_over.'), 
shading interp 
title('Oversampling: length(tslide)=901')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,3)
pcolor(tslide_under,ks/(2*pi),vgtspec_under.'), 
shading interp 
title('Undersampling: length(tslide)=10')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)


%% Part I - Using Different Gabor windows

vmht_spec = zeros(length(tslide), n); 
vsht_spec = zeros(length(tslide), n);
for j = 1:length(tslide)
    %  Mexican hat wavelet
    sigma = 0.05;
    mex_hat = (2/(sqrt(3*sigma)*(pi^0.25))).*(1-((t-tslide(j))/sigma).^2)...
        .* exp(-((t-tslide(j)).^2)/(2*sigma^2));
    % mex_hat = (1-((t-tslide(j)).^2)) .* exp(-((t-tslide(j)).^2)/2);
    
    % Step-function (Shannon) window
    width = 0.05;
    shn = abs(t - tslide(j)) <= width /2;
    
    vmh = mex_hat .* v;
    vmht = fft(vmh);
    
    vsh = shn .* v;
    vsht = fft(vsh);

    vmht_spec(j,:) = fftshift(abs(vmht));
    vsht_spec(j,:) = fftshift(abs(vsht));
end

figure(3)
subplot(3,1,1)
pcolor(tslide,ks/(2*pi),vgt_spec.'), 
shading interp 
title('Gaussian')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,2)
pcolor(tslide,ks/(2*pi),vmht_spec.'), 
shading interp 
title('Mexican hat')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)

subplot(3,1,3)
pcolor(tslide,ks/(2*pi),vsht_spec.'), 
shading interp 
title('Shannon (step-function)')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(hot)


%% Part II - Starter Code
clear; close all; clc

[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (piano)');
% p8 = audioplayer(y,Fs); playblocking(p8);


%% Part II - 1

% We first focus on the spectrogram of the song on the piano 
% to reproduce the music notes.
v = y';
n = length(v);
t2 = linspace(0,tr_piano,n+1); 
t_p = t2(1:n); 
k_p = (2*pi/tr_piano)*[0:n/2-1 -n/2:-1]; 
ks_p = fftshift(k_p);

tslide_p = 0:0.2:tr_piano;
spec_p = zeros(length(tslide_p), n);
notes_piano = zeros(1, length(tslide_p));
for j = 1:length(tslide_p)
    g = exp(-100*(t_p-tslide_p(j)).^2);  % Use Gaussian Filter
    
    vgp = g .* v;
    vgpt = fft(vgp);
    [M, I] = max(vgpt);  % Get the index with strongest frequency (music score)
    
    notes_piano(1,j) = abs(k_p(I))/(2*pi);
    spec_p(j,:) = fftshift(abs(vgpt));
end

% The length and the range of frequencies are different between the 
% music score on the piano and recorder.
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (recorder)');
% p8 = audioplayer(y,Fs); playblocking(p8);

v = y';
n = length(v);
t2 = linspace(0,tr_rec,n+1); 
t_rec = t2(1:n); 
k_rec = (2*pi/tr_rec)*[0:n/2-1 -n/2:-1]; 
ks_rec = fftshift(k_rec);

tslide_rec = 0:0.2:tr_rec;
spec_rec = zeros(length(tslide_rec), n);
notes_rec = zeros(1, length(tslide_rec));
for j = 1:length(tslide_rec)
    g = exp(-100*(t_rec-tslide_rec(j)).^2);  % Use Gaussian Filter
    
    vgr = g .* v;
    vgrt = fft(vgr);
    [M, I] = max(vgrt);  % Get the index with strongest signal (music score)
    
    notes_rec(1,j) = abs(k_rec(I))/(2*pi);
    spec_rec(j,:) = fftshift(abs(vgrt));
end


figure(4)
subplot(1,2,1)
pcolor(tslide_p,(ks_p/(2*pi)), spec_p .'), 
shading interp
title ("Piano")
xlabel('Time (s)'), ylabel('Frequency (Hz)')
ylim ([100 400])
colormap(hot)

subplot(1,2,2)
pcolor(tslide_rec,(ks_rec/(2*pi)), spec_rec .'), 
shading interp
title ("Recorder")
xlabel('Time (s)'), ylabel('Frequency (Hz)')
ylim ([600 1200])
colormap(hot)

%% Reproduce the music score/note on the piano
figure(5)
subplot(1,2,1)
plot(tslide_p, notes_piano,'o');
title ("Piano music score (200~350 Hz)")
xlabel('Time (s)'), ylabel('Music Note')
ylim ([200 350])
yticks([220.00 246.94 261.63 293.66 329.63 349.23])
yticklabels({'A4=220','B4=247','C4=262','D4=294','E4=330','F4=349'})

subplot(1,2,2)
plot(tslide_rec, notes_rec,'o');
title ("Recorder music score (600~1200 Hz)")
xlabel('Time (s)'), ylabel('Music Note')
ylim ([600 1200])
yticks([659.26 698.46 783.99 880.00 987.77 1046.5 1174.7])
yticklabels({'E5=659','F5=698','G5=784','A6=880','B6=988','C6=1046','D6=1175'})



