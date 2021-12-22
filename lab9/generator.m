clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

NX = 39;
NY = 39;

for it = [100 200 500 1000 2000]
    file = fopen(sprintf('%s%d%s', 'results/T_', it, '.txt'), 'r');
    temp = (fscanf(file, '%f'))';
    fclose(file);
    T = zeros(NX);
    for i = 1 : NX
        T(i, :) = temp((i-1) * NX+1 : i * NX);
    end

    fig = figure;
    surf(0:1:NX-1, 0:1:NY-1, T', 'FaceColor', 'TextureMap', 'EdgeColor', 'None');
    xlim([0, NX-1]);
    ylim([0, NY-1]);
    view(2);
    title(sprintf('%s_{it=%d}', 'T', it), "FontSize", 14);
    xlabel('x', "FontSize", 14);
    ylabel('y', "FontSize", 14);
    saveas(fig, sprintf('%s%d%s', 'img/T_', it, '.png'));
    
    file = fopen(sprintf('%s%d%s', 'results/grad2T_', it, '.txt'), 'r');
    temp = (fscanf(file, '%f'))';
    fclose(file);
    T = zeros(NX);
    for i = 1 : NX
        T(i, :) = temp((i-1) * NX+1 : i * NX);
    end

    fig = figure;
    surf(0:1:NX-1, 0:1:NY-1, T', 'FaceColor', 'TextureMap', 'EdgeColor', 'None');
    xlim([0, NX-1]);
    ylim([0, NY-1]);
    view(2);
    title(sprintf('%s_{it=%d}', '\DeltaT', it), "FontSize", 14);
    xlabel('x', "FontSize", 14);
    ylabel('y', "FontSize", 14);
    saveas(fig, sprintf('%s%d%s', 'img/grad2T_', it, '.png'));
end
