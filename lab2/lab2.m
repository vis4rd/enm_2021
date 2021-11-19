clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

N = 500; % populacja (N = u + z)
beta = .001; % czestosc kontaktow osob zarazonych ze zdrowymi
gamma = .1; % sredni czas trwania choroby
tmax = 100;
dt = .1;
t = linspace(0, tmax, tmax/dt);
% z = 0; % osoby zdrowe

alpha = beta*N - gamma;
micron_max = 20;
TOL = 10e-6;
u = ones(1, tmax/dt); % osoby zarazone
u_prev = ones(1, tmax/dt); % poprzednia iteracja

isTOL = @(u_t, u_prev, i) (abs(u_t(i) - u_prev(i)) < TOL);
f = @(U, mi) ( alpha*U(mi) - beta*U(mi)^2 );

% Metoda trapezow

% Metoda Picarda
for iter = 2 : micron_max
    for i = 2 : tmax/dt
        u(i) = u(i-1) + dt/2 * ( f(u, i-1) + f(u_prev, i-1) );
        if isTOL(u, u_prev, i) == true
            break;
        end
        u_prev = u;
    end
end

figure; plot(t, u, t, N-u);
title("Metoda niejawna trapezow - Metoda Picarda", 'FontSize', 14);
xlabel("t", 'FontSize', 14);
ylabel("u(t), z(t)", 'FontSize', 14);
legend("u - zarazeni", "z - zdrowi", 'FontSize', 12);

% Iteracja Newtona
u = ones(micron_max, tmax/dt); % osoby zarazone

for iter = 2 : micron_max
    for i = 2 : tmax/dt
        u(iter, i) = u(iter-1, i) - ( (u(iter-1, i) - u(iter, i-1) - dt/2 * (f(u(iter, :), i-1) + f(u(iter-1, :), i))) /...
            (1 - dt/2 * (alpha - 2*beta*u(iter-1, i))) );
        u_prev = u;
    end
end

figure; plot(t, u(micron_max, :), t, N-u(micron_max, :));
title("Metoda niejawna trapezow - Iteracja Netwona", 'FontSize', 14);
xlabel("t", 'FontSize', 14);
ylabel("u(t), z(t)", 'FontSize', 14);
legend("u - zarazeni", "z - zdrowi", 'FontSize', 12);

% Metoda RK2
c1 = .5 - sqrt(3)/6;
c2 = .5 + sqrt(3)/6;
a11 = .25;
a12 = .25 - sqrt(3)/6;
a21 = .25 + sqrt(3)/6;
a22 = .25;
b1 = .5;
b2 = .5;

m11 = @(U1, mi) ( 1 - dt*a11*(alpha-2*beta*U1(mi)) );
m12 = @(U2, mi) ( -dt*a12*(alpha-2*beta*U2(mi)) );
m21 = @(U1, mi) ( -dt*a21*(alpha-2*beta*U1(mi)) );
m22 = @(U2, mi) ( 1 - dt*a22*(alpha-2*beta*U2(mi)) );

F1 = @(U1, U2, u, mi) ( U1(mi) - u - dt*(a11*f(U1, mi) + a12*f(U2, mi)) );
F2 = @(U1, U2, u, mi) ( U2(mi) - u - dt*(a21*f(U1, mi) + a22*f(U2, mi)) );

div = @(U1, U2, mi) ( m11(U1, mi)*m22(U2, mi) - m12(U2, mi)*m21(U1, mi) );
DU1 = @(U1, U2, u, mi) ( (F2(U1, U2, u, mi)*m12(U2, mi) - F1(U1, U2, u, mi)*m22(U2, mi)) / div(U1, U2, mi) );
DU2 = @(U1, U2, u, mi) ( (F1(U1, U2, u, mi)*m21(U1, mi) - F2(U1, U2, u, mi)*m11(U1, mi)) / div(U1, U2, mi) );

next_u = @(u, i, U1, U2, mi) ( u(i) + dt*(b1*f(U1, round(mi + c1*dt)+1) + b2*f(U2, round(mi + c2*dt)+1)) );
next_U1 = @(U1, U2, u, mi) ( U1(mi) + DU1(U1, U2, u, mi) );
next_U2 = @(U1, U2, u, mi) ( U2(mi) + DU2(U1, U2, u, mi) );

U1 = zeros(1, micron_max);
U2 = zeros(1, micron_max);
u = zeros(1, tmax/dt); u(1) = 1;
up = u;
for mi = 2 : micron_max
    for i = 2 : tmax/dt
        U1(mi) = next_U1(U1, U2, u(i-1), mi-1);
        U2(mi) = next_U2(U1, U2, u(i-1), mi-1);
        u(i) = next_u(u, i-1, U1, U2, mi-1);
        if isTOL(u, up, i) == true
            break;
        end
        up = u;
    end
end

figure; plot(t, u, t, N-u);
title("Metoda niejawna - RK2", 'FontSize', 14);
xlabel("t", 'FontSize', 14);
ylabel("u(t), z(t)", 'FontSize', 14);
legend("u - zarazeni", "z - zdrowi", 'FontSize', 12);
