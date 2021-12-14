clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

drawU(401, 91, 0);
drawU(401, 91, 1);

cx0 = load("build/c_avgx_D_0.0.txt");
c0 = cx0(:, 1);
avgx0 = cx0(:, 2);

cx1 = load("build/c_avgx_D_0.1.txt");
c1 = cx1(:, 1);
avgx1 = cx1(:, 2);

t = 0:0.6/10000:0.6-0.00001;

figure;
plot(t, c0, t, avgx0);
ylabel("x_{st} (c_{sr})", "FontSize", 14);
xlabel("t", "FontSize", 14);
title("D = 0", "FontSize", 14);
legend("c_{sr}", "x_{st}", "FontSize", 14);
ylim([0, 4]);

figure;
plot(t, c1, t, avgx1);
ylabel("x_{st} (c_{sr})", "FontSize", 14);
xlabel("t", "FontSize", 14);
title("D = 0.1", "FontSize", 14);
legend("c_{sr}", "x_{st}", "FontSize", 14);
ylim([0, 4]);

vx_t = load("build/vx_D_0.0.txt");
vy_t = load("build/vy_D_0.0.txt");
vx = zeros(401, 91);
vy = zeros(401, 91);
for j = 1 : 401
    vx(j, :) = vx_t((j-1)*91+1 : j*91);
    vy(j, :) = vy_t((j-1)*91+1 : j*91);
end
fig = figure;
surf(0:.01:(401-1)/100, 0:.01:(91-1)/100, vx');
shading("flat");
view(2);
saveas(fig, "build/vx.png");
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("V_x(x, y)", "FontSize", 14);
fig = figure;
surf(0:.01:(401-1)/100, 0:.01:(91-1)/100, vy');
shading("flat");
view(2);
saveas(fig, "build/vy.png");
xlabel("x", "FontSize", 14);
ylabel("y", "FontSize", 14);
title("V_y(x, y)", "FontSize", 14);

function[] = print(nx, ny, target, iter, D)
    fig = figure;
    surf(0:.01:(nx-1)/100, 0:.01:(ny-1)/100, target');
    shading("flat");
    view(2);
    title("u(x, y), t = "+(iter+1)/5+"t_{MAX}", "FontSize", 14);
    xlabel("x", "FontSize", 14);
    ylabel("y", "FontSize", 14);
    saveas(fig, "build/AD_t_"+(iter+1)/5+"_D_"+D+".png");
end

function[] = drawU(nx, ny, D)
    for i = 0 : 4
        u_t = load("build/u_t_"+i+"_D_0."+D+".txt");
        u = zeros(nx, ny);
        for j = 1 : nx
            u(j, :) = u_t((j-1)*ny+1 : j*ny);
        end
        print(nx, ny, u, i, D);
    end
end
