clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

nx = 150;
nt = 1000;

u00t = load("u00.txt");
u01t = load("u01.txt");
u10t = load("u10.txt");
uat = load("ua.txt");

for j = 1 : nt
    u00(j, :) = u00t((j-1)*(nx+1)+1 : j*(nx+1));
    u01(j, :) = u01t((j-1)*(nx+1)+1 : j*(nx+1));
    u10(j, :) = u10t((j-1)*(nx+1)+1 : j*(nx+1));
    ua(j, :) = uat((j-1)*(nx+1)+1 : j*(nx+1));
end

delta = 0.1;
dt = 0.05;
ox = dt*(1:1:nt);
oy = delta*(0:1:nx);

fig = figure;
surf(ox, oy, u00');
shading("flat");
view(2);
xlabel("x", "FontSize", 14);
ylabel("t", "FontSize", 14);
title("\alpha = 0, \beta = 0", "FontSize", 14);
saveas(fig, "u00.png");

fig = figure;
surf(ox, oy, u01');
shading("flat");
view(2);
xlabel("x", "FontSize", 14);
ylabel("t", "FontSize", 14);
title("\alpha = 0, \beta = 0.1", "FontSize", 14);
saveas(fig, "u01.png");

fig = figure;
surf(ox, oy, u10');
shading("flat");
view(2);
xlabel("x", "FontSize", 14);
ylabel("t", "FontSize", 14);
title("\alpha = 0, \beta = 1", "FontSize", 14);
saveas(fig, "u10.png");

e00t = load("e00.txt");
e01t = load("e01.txt");
e10t = load("e10.txt");
eat = load("ea.txt");

fig = figure;
plot(ox, eat);
xlabel("t", "FontSize", 14);
ylabel("E", "FontSize", 14);
title("\alpha = 1, \beta = 1", "FontSize", 14);
saveas(fig, "ea.png");

fig = figure;
surf(ox, oy, ua');
shading("flat");
view(2);
xlabel("x", "FontSize", 14);
ylabel("t", "FontSize", 14);
title("\alpha = 1, \beta = 1", "FontSize", 14);
saveas(fig, "ua.png");
