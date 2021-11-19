clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

trapez1 = load("trapez.txt");
t = trapez1(:, 1);
dt = trapez1(:, 2);
x = trapez1(:, 3);
v = trapez1(:, 4);

trapez2 = load("trapez2.txt");
t2 = trapez2(:, 1);
dt2 = trapez2(:, 2);
x2 = trapez2(:, 3);
v2 = trapez2(:, 4);

figure;
plot(t, v, t2, v2);
xlabel("t", "FontSize", 14);
ylabel("v(t)", "FontSize", 14);
title("Metoda trapezow", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(t, x, t2, x2);
xlabel("t", "FontSize", 14);
ylabel("x(t)", "FontSize", 14);
title("Metoda trapezow", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(t, dt, t2, dt2);
xlabel("t", "FontSize", 14);
ylabel("\Deltat(t)", "FontSize", 14);
title("Metoda trapezow", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(x, v, x2, v2);
xlabel("x", "FontSize", 14);
ylabel("v(x)", "FontSize", 14);
title("Metoda trapezow", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

trapez1 = load("RK2.txt");
t = trapez1(:, 1);
dt = trapez1(:, 2);
x = trapez1(:, 3);
v = trapez1(:, 4);

trapez2 = load("RK22.txt");
t2 = trapez2(:, 1);
dt2 = trapez2(:, 2);
x2 = trapez2(:, 3);
v2 = trapez2(:, 4);

figure;
plot(t, v, t2, v2);
xlabel("t", "FontSize", 14);
ylabel("v(t)", "FontSize", 14);
title("Metoda RK2", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(t, x, t2, x2);
xlabel("t", "FontSize", 14);
ylabel("x(t)", "FontSize", 14);
title("Metoda RK2", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(t, dt, t2, dt2);
xlabel("t", "FontSize", 14);
ylabel("\Deltat(t)", "FontSize", 14);
title("Metoda RK2", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);

figure;
plot(x, v, x2, v2);
xlabel("x", "FontSize", 14);
ylabel("v(x)", "FontSize", 14);
title("Metoda RK2", "FontSize", 14);
legend("TOL = 10^{-2}", "TOL = 10^{-5}", "FontSize", 14);