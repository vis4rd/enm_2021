clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

grid_s1 = load("bin/grid_k1_s.txt");
grid_s2 = load("bin/grid_k2_s.txt");
grid_s4 = load("bin/grid_k4_s.txt");
grid_s8 = load("bin/grid_k8_s.txt");
grid_s16 = load("bin/grid_k16_s.txt");

gs1_it = grid_s1(:, 1);
gs2_it = grid_s2(:, 1);
gs4_it = grid_s4(:, 1);
gs8_it = grid_s8(:, 1);
gs16_it = grid_s16(:, 1);

gs1_s = grid_s1(:, 2);
gs2_s = grid_s2(:, 2);
gs4_s = grid_s4(:, 2);
gs8_s = grid_s8(:, 2);
gs16_s = grid_s16(:, 2);

fig = figure;
plot(gs1_it, gs1_s, gs2_it, gs2_s, gs4_it, gs4_s, gs8_it, gs8_s, gs16_it, gs16_s);
xlabel("it", "FontSize", 14);
ylabel("S(it)", "FontSize", 14);
title("Relaksacja wielosiatkowa", "FontSize", 14);
legend("k = 1", "k = 2", "k = 4", "k = 8", "k = 16", "FontSize", 14);
saveas(fig, "img/gridrel_s_it.png");

print_fn(1);
print_fn(2);
print_fn(4);
print_fn(8);
print_fn(16);

function [] = print_fn(k)
    grid_k = load("bin/grid_k"+k+"_v.txt");
    gkv = zeros(129, 129);
    for i = 1 : 129
        gkv(i, :) = grid_k((i-1)*129+1 : i*129);
    end
    fig = figure;
    surf(0:.2*k:25.6, 0:.2*k:25.6, gkv(1:k:end, 1:k:end)', 'FaceColor','texturemap')
    shading("flat");
    view(2);
    xlabel("x", "FontSize", 14)
    ylabel("y", "FontSize", 14)
    xlim([0, 25.6]);
    ylim([0, 25.6]);
    title("Mapa zrelaksowanego potencjalu V(x, y) dla k = "+k, "FontSize", 14)
    saveas(fig, "img/map_k_"+k+".png")
end
