clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

local1 = load("src/loc_s_1.0.txt");
local2 = load("src/loc_s_1.4.txt");
local3 = load("src/loc_s_1.8.txt");
local4 = load("src/loc_s_1.9.txt");

glob_v1 = load("src/glob_v_0.6.txt");
glob_v2 = load("src/glob_v_1.0.txt");
glob_err1 = load("src/glob_err_0.6.txt");
glob_err2 = load("src/glob_err_1.0.txt");

glob_s1 = load("src/glob_s_0.6.txt");
glob_s2 = load("src/glob_s_1.0.txt");

liter1 = local1(:, 1);
liter2 = local2(:, 1);
liter3 = local3(:, 1);
liter4 = local4(:, 1);
ls1 = local1(:, 2);
ls2 = local2(:, 2);
ls3 = local3(:, 2);
ls4 = local4(:, 2);

gx1 = glob_v1(:, 1);
gx2 = glob_v2(:, 1);
gy1 = glob_v1(:, 2);
gy2 = glob_v2(:, 2);
gerr1t = glob_err1(:, 3);
gerr2t = glob_err2(:, 3);
gerr1 = zeros(151, 101);
gerr2 = zeros(151, 101);
gv1t = glob_v1(:, 3);
gv2t = glob_v2(:, 3);
gv1 = zeros(151, 101);
gv2 = zeros(151, 101);
for i = 1 : 151
    gv1(i, :) = gv1t((i-1)*101+1 : i*101);
    gv2(i, :) = gv2t((i-1)*101+1 : i*101);
    gerr1(i, :) = gerr1t((i-1)*101+1 : i*101);
    gerr2(i, :) = gerr2t((i-1)*101+1 : i*101);
end

git1 = glob_s1(:, 1);
git2 = glob_s2(:, 1);
gs1 = glob_s1(:, 2);
gs2 = glob_s2(:, 2);

fig = figure;
semilogx(liter1, ls1, liter2, ls2, liter3, ls3, liter4, ls4);
xlabel("it", "FontSize", 14);
ylabel("S(it)", "FontSize", 14);
title("Relaksacja lokalna", "FontSize", 14);
legend("\omega=1.0, it="+size(liter1, 1), ...
    "\omega=1.4, it="+size(liter2, 1), ...
    "\omega=1.8, it="+size(liter3, 1), ...
    "\omega=1.9, it="+size(liter4, 1), ...
    "FontSize", 14);
saveas(fig, "img/rel_loc.png");

fig = figure;
semilogx(git1, gs1, git2, gs2);
xlabel("it", "FontSize", 14);
ylabel("S(it)", "FontSize", 14);
title("Relaksacja globalna", "FontSize", 14);
legend("\omega=0.6, it="+size(git1, 1), ...
    "\omega=1.0, it="+size(git2, 1), ...
    "FontSize", 14);
saveas(fig, "img/rel_glob.png");

x = 0:.1:15;
y = 0:.1:10;
figure; contourf(x, y, gv1');
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("Relaksacja globalna - zrelaksowany potencjal", "FontSize", 14);
saveas(fig, "img/map_0_6.png");
figure; contourf(x, y, gv2');
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("Relaksacja globalna - zrelaksowany potencjal", "FontSize", 14);
saveas(fig, "img/map_1_0.png");
figure; contourf(x, y, gerr1');
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("Relaksacja globalna - blad rozwiazania", "FontSize", 14);
saveas(fig, "img/map_err_0_6.png");
figure; contourf(x, y, gerr2');
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("Relaksacja globalna - blad rozwiazania", "FontSize", 14);
saveas(fig, "img/map_err_1_0.png");