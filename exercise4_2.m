close all
clear 
clc

im_ams = imread('amsterdam.bmp');

figure (1)
imagesc(im_ams)
axis image
axis on
grid on
grid minor

hold on
plot([383 365], [333 74], 'r','linewidth',3);
plot(377, 491, 'm*', 'markersize', 10);

max_red = max(max(im_ams(:,:,1)));
max_green = max(max(im_ams(:,:,1)));
max_blue = max(max(im_ams(:,:,1)));

[row_red, col_red,val_r] = find (im_ams(:,:,1) == 254);
[row_green, col_green,val_g] = find(im_ams(:,:,2) == 254);
[row_blue, col_blue,val_b] = find(im_ams(:,:,3) == 254);

max_ind_red = randi (length(row_red),1,1);
max_ind_green = randi (length(row_red),1,1);
max_ind_blue = randi (length(row_red),1,1);

plot(row_red(max_ind_red),col_red(max_ind_red),'ro');

plot(row_green(max_ind_green),col_green(max_ind_green),'go');

plot(row_blue(max_ind_blue),col_blue(max_ind_blue),'bo');

