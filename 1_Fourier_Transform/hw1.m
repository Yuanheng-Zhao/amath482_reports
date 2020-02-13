%{ 
AMATH 482
Professor: Craig Gin
HW#1: An ultrasound problem 
Jonathan Zhao
%}
%% Part1: Initialization 
clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); 
x=x2(1:n); 
y=x; 
z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%% Part2-1: Adding up the signals of all realizations
ave = zeros(n, n, n);
for realize=1:20
    Un(:,:,:)=reshape(Undata(realize,:),n,n,n);
    Utn = fftn(Un);
    ave = ave + Utn;
    
%     close all, isosurface(X,Y,Z,abs(Un),0.4)
%     axis([-20 20 -20 20 -20 20]), grid on, drawnow
%     xlabel('x'); ylabel('y'); zlabel('z');
%     pause
end

%% Part2-2: Averaging & Finding the center frequency
ave = abs(fftshift(ave));
flatAve = reshape(ave, n^3, 1);
strongest = max(flatAve);
ave = ave/strongest;
close all, isosurface(Kx,Ky,Kz,ave,0.8)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]), grid on, drawnow
xlabel('Kx'); ylabel('Ky'); zlabel('Kz');

ave = ifftshift(ave); % shift back frequency domain
% look for the strongest signal
for m = 1:n
    [j, i] = find(ave(:,:,m) == 1); 
    if isempty(i)~=1
        center_frequency = [i, j, m];
        break
    end
end
kx=k(center_frequency(1)); 
ky=k(center_frequency(2)); 
kz=k(center_frequency(3)); 

%% Part3-1: Create a Filter
tau = 1.0;
filter = exp(-tau*(fftshift(Kx)-kx).^2).*...
    exp(-tau*(fftshift(Ky)-ky).^2).*...
    exp(-tau*(fftshift(Kz)-kz).^2);

%% Part3-2: Filtering & locating the marble in each realization
position = zeros(20, 3);
for realize = 1:20
    % apply filter on each realization
    Un(:,:,:)=reshape(Undata(realize,:),n,n,n);
    Utn = fftn(Un);
    unft = filter.* Utn;
    unf = ifftn(unft); % inverse multi-dimensional FFT
    unf = abs(unf);
    flat_unf = reshape(unf, n^3, 1);
    strongest = max(flat_unf);
    % look for the strongest signal
    for m = 1:n
        [j, i] = find(unf(:,:,m) == strongest);
        if isempty(i)~=1
            position(realize,:) = [x(i), y(j), z(m)];
            break
        end
    end

    isosurface(X,Y,Z,unf/strongest,0.8)
    axis([-L L -L L -L L]), grid on, drawnow
    xlabel('x'); ylabel('y');zlabel('z');
    pause(1)
end

%% Part4: Solution
% plot the path of the marble 
plot3(position(:,1), position(:,2), position(:,3), '.-', 'MarkerSize',...
    10, 'Linewidth', 2)
axis([-L L -L L -L L]), grid on, drawnow
xlabel('x'); ylabel('y'); zlabel('z');

% Center frequency
[kx, ky, kz]
% The position of the marble at the 20th data measurement.
solution = [position(20,1), position(20,2),position(20,3)];
