data = csvread('CON_s9305_Da_mean.csv');
data2 = csvread('CON_s9305_Dr_mean.csv');
figure;
mask = data > 0;
imagesc(data2);
axis square
colorbar
figure;
plot(data(mask), data2(mask), 'o')
