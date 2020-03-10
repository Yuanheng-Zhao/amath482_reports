%% Load the data (videos) for Test 1/Case 1
clear all; close all; clc

load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')
% implay(vidFrames1_1)
%%
filter = zeros(480 ,640);  
filter(200:430, 300:400) = 1;  % Filter for Case1 - Camera1
data1 = load_cropped_data(vidFrames1_1, filter, 240);
%%
filter = zeros(480 ,640);  
filter(100:390, 240:350) = 1;  % Filter for Case1 - Camera2
data2 = load_cropped_data(vidFrames2_1, filter, 240);
%%
filter = zeros(480 ,640);  
filter(230:330, 260:480) = 1;  % Filter for Case1 - Camera3
data3 = load_cropped_data(vidFrames3_1, filter, 240);
%%
collected_data = collect(data1, data2, data3);
%% Results for Case 1
[m,n]= size(collected_data);  % compute data size

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  % subtract mean
[U,S,V]= svd(collected_data'/sqrt(n-1));  % perform the SVD
lambda = diag(S).^2;  % produce diagonal variances
Y = collected_data' * V;  % produce the principal components projection

figure(1)
plot(1:6, lambda/sum(lambda), '-*', 'Linewidth', 2);
title("Case 1 - Dynamics of Interests");
xlabel("Principal Components"); ylabel("Diagonal Variances");
y = lambda/sum(lambda);
for i=1:2   
    text(i,y(i),num2str(y(i)));
end

figure(2)
subplot(2,1,1)
plot(1:216, collected_data(2,:), 1:216, collected_data(1,:), 'Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
title("Case 1 - Original displacement");
legend("Z", "XY")
subplot(2,1,2)
plot(1:216, Y(:,1),'r','Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
title("Case 1 - Principal Component Projection(s)");
legend("Principal Component1")



%% Load the data for Test 2/Case 2
clear all; close all; clc;

load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

%%
filter = zeros(480 ,640);  
filter(200:400, 300:430) = 1;  % Filter for Case2 - Camera1
data1 = load_cropped_data(vidFrames1_2, filter, 240);
%%
filter = zeros(480 ,640);  
filter(50:420, 180:440) = 1;  % Filter for Case2 - Camera2
data2 = load_cropped_data(vidFrames2_2, filter, 240);
%%
filter = zeros(480 ,640);  
filter(180:330, 280:470) = 1;  % Filter for Case2 - Camera3
data3 = load_cropped_data(vidFrames3_2, filter, 240);
%%
collected_data = collect(data1, data2, data3);
%% Results for Test 2
[m,n]= size(collected_data);  % compute data size

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  % subtract mean
[U,S,V]= svd(collected_data'/sqrt(n-1));  % perform the SVD
lambda = diag(S).^2;  % produce diagonal variances
Y = collected_data' * V;  % produce the principal components projection

figure(3)
plot(1:6, lambda/sum(lambda), '-*', 'Linewidth', 2);
title("Case 2 - Dynamics of Interests") ;
xlabel("Principal Components"); ylabel("Diagonal Variances");
y = lambda/sum(lambda);
for i=1:4   
    text(i,y(i),num2str(y(i)));
end

figure(4)
subplot(2,1,1)
plot(1:295, collected_data(2 ,:), 1:295, collected_data(1 ,:) ,'Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
legend("Z", "XY")
title("Case 2 - Original displacement");
subplot(2,1,2)
plot(1:295, Y(:,1), 1:295, Y(:,2), 1:295, Y(:,3), 'r', 'Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
title("Case 2 - Principal Component Projection(s)");
legend("PC1" , " PC2")



%% Load the data for Test 3/Case 3
clear all; close all; clc;

load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')
%%
filter = zeros(480 ,640);  
filter(240:420, 280:380) = 1;  % Filter for Case3 - Camera1
data1 = load_cropped_data(vidFrames1_3, filter, 240);
%%
filter = zeros(480 ,640);  
filter(180:380, 240:400) = 1;  % Filter for Case3 - Camera2
data2 = load_cropped_data(vidFrames2_3, filter, 240);
%%
filter = zeros(480 ,640);  
filter(180:330, 240:480) = 1;  % Filter for Case3 - Camera3
data3 = load_cropped_data(vidFrames3_3, filter, 240);
%%
collected_data = collect(data1, data2, data3);
%% Results for Test 3
[m,n]= size(collected_data);  % compute data size

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  % subtract mean
[U,S,V]= svd(collected_data'/sqrt(n-1));  % perform the SVD
lambda = diag(S).^2;  % produce diagonal variances
Y = collected_data' * V;  % produce the principal components projection

figure(5)
plot(1:6, lambda/sum(lambda), '-*', 'Linewidth', 2);
title("Case 3 - Dynamics of Interests");
xlabel("Principal Components"); ylabel("Diagonal Variances");
y = lambda/sum(lambda);
for i=1:4   
    text(i,y(i),num2str(y(i)));
end

figure(6)
subplot(2,1,1)
plot(1:235, collected_data(2 ,:), 1:235, collected_data(1,:),'Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
title("Case 3 - Original displacement");
legend("Z", "XY")
subplot(2,1,2)
plot(1:235,Y(:,1),1:235,Y(:,2),1:235,Y(:,3),1:235,Y(:,4),'r','Linewidth',2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)");
title("Case 3 - Principal Component Projection(s)");
legend("PC1", "PC2", "PC3", "PC4")



%% Load the data for Test 4/Case 4
clear all; close all; clc;

load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')
%%
filter = zeros(480 ,640);  
filter(230:440, 330:460) = 1;  % Filter for Case4 - Camera1
data1 = load_cropped_data(vidFrames1_4, filter, 240);
%%
filter = zeros(480 ,640);  
filter(100:360, 230:410) = 1;  % Filter for Case4 - Camera2
data2 = load_cropped_data(vidFrames2_4, filter, 240);
%%
filter = zeros(480 ,640);  
filter(140:280, 320:500) = 1;  % Filter for Case4 - Camera3
data3 = load_cropped_data(vidFrames3_4, filter, 230);
%%
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);

% Trim the data to make them a consistent length.
% For Test4, the video recorded by camera 3 is the shortest. 
% Thus trim the other two as the length of video 3.
data1 = data1(1:length(data3), :);
data2 = data2(1:length(data3), :);

collected_data = [data1'; data2'; data3'];

%% Results for Test 4
[m,n]= size(collected_data);  % compute data size

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  % subtract mean
[U,S,V]= svd(collected_data'/sqrt(n-1));  % perform the SVD
lambda = diag(S).^2;  % produce diagonal variances
Y = collected_data' * V;  % produce the principal components projection

figure(7)
plot(1:6, lambda/sum(lambda), '-*', 'Linewidth', 2);
title("Case 4 - Dynamics of Interests");
xlabel("Principal Components"); ylabel("Diagonal Variances");
y = lambda/sum(lambda);
for i=1:4
    text(i,y(i),num2str(y(i)));
end

figure(8)
subplot(2,1,1)
plot(1:376 , collected_data(2,:), 1:376, collected_data(1,:),'Linewidth', 2)
ylabel("Displacement(pixels)") ; xlabel("Time(frames)") ;
title("Case 4 - Original displacement") ;
legend("Z", "XY")
subplot(2,1,2)
plot(1:376, Y(:,1), 1:376, Y(:,2), 1:376, Y(:,3), 'Linewidth', 2)
ylabel("Displacement(pixels)"); xlabel("Time(frames)") ;
title("Case 4 - Principal Component Projection(s)") ;
legend("PC1", "PC2", "PC3")

