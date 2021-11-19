clear; %czysci wszystkie zmienne w srodowisku itp
clc; %czysci konsole
clf; %clear figure
close all; %zamyka wszystkie podokna

tmin = 0;
tmax = 5;
dt1 = 0.01;
dt2 = 0.1;
dt3 = 1;
n1 = ((tmax - tmin) / dt1) + 1;
n2 = ((tmax - tmin) / dt2) + 1;
n3 = ((tmax - tmin) / dt3) + 1;

t = linspace(0, 5, n1);
t2 = linspace(0, 5, n2);
t3 = linspace(0, 5, n3);

y = zeros(1, n1); y(1) = 1;
y2 = zeros(1, n2); y2(1) = 1;
y3 = zeros(1, n3); y3(1) = 1;

lambda = -1;
f = @(y, i)(lambda * y(i));

for i = 2 : n1
    y(i) = y(i-1) + dt1 * f(y, i-1);
end
for i = 2 : n2
    y2(i) = y2(i-1) + dt2 * f(y2, i-1);
end
for i = 2 : n3
    y3(i) = y3(i-1) + dt3 * f(y3, i-1);
end

figure;
title('Metoda jawna Eulera', 'FontSize', 16);
hold on;
plot(t, y, '-', 'LineWidth', 2);
plot(t2, y2, '-', 'LineWidth', 2);
plot(t3, y3, '-', 'LineWidth', 2);
plot(t, exp(lambda * t), '--', 'LineWidth', 4);
xlabel('t');
ylabel('y(t)');
legend('\Deltat=0.01', '\Deltat=0.1', '\Deltat=1', 'e^\lambda^t', 'FontSize',16);
set(gca,'FontSize',16);
hold off;

figure;
subplot(1, 3, 1);
plot(t, y - exp(lambda * t), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
subplot(1, 3, 2);
plot(t2, y2 - exp(lambda * t2), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
sgtitle('Blad globalny \delta(t) = y_{num}(t) - y_{dok}(t)', 'FontSize', 16);
subplot(1, 3, 3);
plot(t3, y3 - exp(lambda * t3), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');

T1 = table(t', y', (y - exp(lambda*t))');
T2 = table(t2', y2', (y2 - exp(lambda*t2))');
T3 = table(t3', y3', (y3 - exp(lambda*t3))');
writetable(T1, 'eulerd1.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T1, 'eulerd2.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T1, 'eulerd3.txt', 'Delimiter',' ', 'WriteVariableNames', false);

% RK2

k1 = @(y, i)(f(y, i));
k2 = @(y, i, dt)(lambda*(y(i) + dt*k1(y, i)));

for i = 2 : n1
    y(i) = y(i-1) + dt1/2 * (k1(y, i-1) + k2(y, i-1, dt1));
end
for i = 2 : n2
    y2(i) = y2(i-1) + dt2/2 * (k1(y2, i-1) + k2(y2, i-1, dt2));
end
for i = 2 : n3
    y3(i) = y3(i-1) + dt3/2 * (k1(y3, i-1) + k2(y3, i-1, dt3));
end

figure;
title('Metoda RK2', 'FontSize', 16);
hold on;
plot(t, y, '-', 'LineWidth', 2);
plot(t2, y2, '-', 'LineWidth', 2);
plot(t3, y3, '-', 'LineWidth', 2);
plot(t, exp(lambda * t), '--', 'LineWidth', 4);
xlabel('t');
ylabel('y(t)');
legend('\Deltat=0.01', '\Deltat=0.1', '\Deltat=1', 'e^\lambda^t', 'FontSize',16);
set(gca,'FontSize',16);
hold off;

figure;
sgtitle('RK2 - Blad globalny \delta(t) = y_{num}(t) - y_{dok}(t)', 'FontSize', 16);
subplot(1, 2, 1);
hold on;
plot(t, y - exp(lambda * t), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
plot(t2, y2 - exp(lambda * t2), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
hold off;
subplot(1, 2, 2);
plot(t3, y3 - exp(lambda * t3), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');

T1 = table(t', y', (y - exp(lambda*t))');
T2 = table(t2', y2', (y2 - exp(lambda*t2))');
T3 = table(t3', y3', (y3 - exp(lambda*t3))');
writetable(T1, 'rk2d1.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T2, 'rk2d2.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T3, 'rk2d3.txt', 'Delimiter',' ', 'WriteVariableNames', false);

% RK4

k1 = @(y, i)(f(y, i));
k2 = @(y, i, dt)(lambda*(y(i) + dt/2*k1(y, i)));
k3 = @(y, i, dt)(lambda*(y(i) + dt/2*k2(y, i, dt)));
k4 = @(y, i, dt)(lambda*(y(i) + dt*k3(y, i, dt)));

for i = 2 : n1
    y(i) = y(i-1) + dt1/6 * (k1(y, i-1) + 2*k2(y, i-1, dt1) + 2*k3(y, i-1, dt1) + k4(y, i-1, dt1));
end
for i = 2 : n2
    y2(i) = y2(i-1) + dt2/6 * (k1(y2, i-1) + 2*k2(y2, i-1, dt2) + 2*k3(y2, i-1, dt2) + k4(y2, i-1, dt2));
end
for i = 2 : n3
    y3(i) = y3(i-1) + dt3/6 * (k1(y3, i-1) + 2*k2(y3, i-1, dt3) + 2*k3(y3, i-1, dt3) + k4(y3, i-1, dt3));
end

figure;
title('Metoda RK4', 'FontSize', 16);
hold on;
plot(t, y, '-', 'LineWidth', 2);
plot(t2, y2, '-', 'LineWidth', 2);
plot(t3, y3, '-', 'LineWidth', 2);
plot(t, exp(lambda * t), '--', 'LineWidth', 4);
xlabel('t');
ylabel('y(t)');
legend('\Deltat=0.01', '\Deltat=0.1', '\Deltat=1', 'e^\lambda^t', 'FontSize',16);
set(gca,'FontSize',16);
hold off;

figure;
sgtitle('RK4 - Blad globalny \delta(t) = y_{num}(t) - y_{dok}(t)', 'FontSize', 16);
subplot(1, 2, 1);
hold on;
plot(t, y - exp(lambda * t), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
plot(t2, y2 - exp(lambda * t2), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');
hold off;
subplot(1, 2, 2);
plot(t3, y3 - exp(lambda * t3), '-', 'LineWidth', 2);
xlabel('t');
ylabel('\delta(t)');

T1 = table(t', y', (y - exp(lambda*t))');
T2 = table(t2', y2', (y2 - exp(lambda*t2))');
T3 = table(t3', y3', (y3 - exp(lambda*t3))');
writetable(T1, 'rk4d1.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T2, 'rk4d2.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T3, 'rk4d3.txt', 'Delimiter',' ', 'WriteVariableNames', false);
