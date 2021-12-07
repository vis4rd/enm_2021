clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

setup(201, 91, -1000);
setup(201, 91, -4000);
setup(201, 91, 4000);

function[] = print(nx, ny, target, ntitle, xmin, xmax, ymin, ymax, q, drawlines)
    fig = figure;
    if drawlines == false
        surf(0:.01:(nx-1)/100, 0:.01:(ny-1)/100, target', 'FaceColor', 'texturemap');
        shading("flat");
    else
        contourf(0:.01:(nx-1)/100, 0:.01:(ny-1)/100, target', 34);
    end
    view(2);
    title(ntitle+"(x, y)", "FontSize", 14);
    xlabel("x", "FontSize", 14)
    ylabel("y", "FontSize", 14)
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    saveas(fig, "map_"+ntitle+"_q_"+q+".png");
end

function[] = setup(nx, ny, Q)
    res = load("build/results_"+Q+".txt");
    psi_t = res(:, 1);
    zeta_t = res(:, 2);
    v_x_t = res(:, 3);
    v_y_t = res(:, 4);
    
    v_x = zeros(nx, ny);
    v_y = zeros(nx, ny);
    psi = zeros(nx, ny);
    zeta = zeros(nx, ny);
    
    for i = 1 : nx
        v_x(i, :) = v_x_t((i-1)*ny+1 : i*ny);
        v_y(i, :) = v_y_t((i-1)*ny+1 : i*ny);
        psi(i, :) = psi_t((i-1)*ny+1 : i*ny);
        zeta(i, :) = zeta_t((i-1)*ny+1 : i*ny);
    end
    print(201, 91, psi, "psi", 0.01, 2, 0, .90, Q, true);
    print(201, 91, zeta, "zeta", 0.01, 2, 0, .90, Q, true);
    print(201, 91, v_x, "u", 0.01, 2, 0, .90, Q, false);
    print(201, 91, v_y, "v", 0.01, 2, 0, .90, Q, false);
end
