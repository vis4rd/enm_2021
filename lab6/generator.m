clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

print_fn(50, 10, 1, 1);
print_fn(100, 10, 1, 1);
print_fn(200, 10, 1, 1);
print_fn(100, 0, 1, 1);
print_fn(100, 0, 1, 2);
print_fn(100, 0, 1, 10);

function [s_v] = print_fn(nx, V, eps1, eps2)
    sparse = load("bin/V_"+nx+"_eps1_"+eps1+"_eps2_"+eps2+"_V_"+V+"_data.txt");
    nx = nx+1;
    s_v = zeros(nx, nx);
    sparse_v = sparse(:, 4);
    for i = 1 : nx
        s_v(i, :) = sparse_v((i-1)*nx+1 : i*nx);
    end
    fig = figure;
    surf(0:0.1:(nx-1)/10, 0:0.1:(nx-1)/10, s_v, 'FaceColor','texturemap');
    shading("flat");
    view(2);
    xlabel("x", "FontSize", 14)
    ylabel("y", "FontSize", 14)
    xlim([0, (nx-1)/10]);
    ylim([0, (nx-1)/10]);
    title("Mapa potencjalu V(x, y) dla nx,ny = "+(nx-1)+", eps1 = "+eps1+", eps2 = "+eps2, "FontSize", 14)
    saveas(fig, "map_nx_"+(nx-1)+"eps1_"+eps1+"eps2_"+eps2+"V_"+V+".png")
end
