close all;
clc;
clear;

load lighthouse;    %loads the lightouse image into variable xx

figure(1);  %Figure 1 shows the orignal image right next to the down-sampled image
hold on;
subplot(1, 2, 1);
%title('Original Image');
imshow(xx); %Displays original image

xp = xx(1:2:end,1:2:end);   %Down-sample the original image by a factor of 2
size(xp);   %Compensate for the size of the down-sampled image
subplot(1, 2, 2);
%title('Down-Sample Image');
imshow(xp);  %Display down-sampled image
hold off;

figure(2);  %Display original image by itself
%title('Original Image');
imshow(xx);

figure(3);  %Display down-sampled image by itself
%title('Down-Sample Image');
imshow(xp);

%Plots for the 2D-FFT of the images for visual representation of the
%frequency domain

% xf = fft2(xx);
% figure(3);
% hold on;
% imshow(abs(fftshift(xf)));
% hold off;
% 
% xfp = fft2(xp);
% figure(4);
% hold on;
% imshow(abs(fftshift(xfp)));
% hold off;
% 
% figure(5);
% hold on;
% imagesc(abs(fftshift(xf)))
% hold off;
% 
% figure(6);
% hold on;
% imagesc(abs(fftshift(xfp)))
% hold off;