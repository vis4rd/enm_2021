clear;
clc;
clf;
close all;

dt = 10^-4;
R = 100.0;
L = .1;
C = 10^-3;
omega_0 = 1.0 / sqrt(L*C);
T_0 = 2 * pi / omega_0;

tmin = 0;
tmax = 4 * T_0;
n = ceil((tmax - tmin)/dt);

omega_v1 = .5 * omega_0;
omega_v2 = .8 * omega_0;
omega_v3 = omega_0;
omega_v4 = 1.2 * omega_0;

f = @(I, i, k)(I(i) + k);
g = @(V, I, Q, i, kq, ki)( (V / L) - ((Q(i) + kq) / (L * C)) - (R * (I(i) + ki) / L));
V = @(om, i, shift)(10 * sin(om * i * dt));

I = zeros(1, n);
Q = zeros(1, n);
I2 = zeros(1, n);
Q2 = zeros(1, n);
I3 = zeros(1, n);
Q3 = zeros(1, n);
I4 = zeros(1, n);
Q4 = zeros(1, n);

kq1 = @(i, I) (f(I, i, 0.0));
ki1 = @(i, om, Q, I) (g(V(om, i), I, Q, i, 0.0, 0.0));
kq2 = @(i, om, Q, I) (f(I, i, dt/2 * ki1(i, om, Q, I)));
ki2 = @(i, om, Q, I) (g(V(om, i+.5), I, Q, i, dt/2 * kq1(i, I), dt/2 * ki1(i, om, Q, I)));
kq3 = @(i, om, Q, I) (f(I, i, dt/2 * ki2(i, om, Q, I)));
ki3 = @(i, om, Q, I) (g(V(om, i+.5), I, Q, i, dt/2 * kq2(i, om, Q, I), dt/2 * ki2(i, om, Q, I)));
kq4 = @(i, om, Q, I) (f(I, i, dt * ki3(i, om, Q, I)));
ki4 = @(i, om, Q, I) (g(V(om, i+1), I, Q, i, dt * kq3(i, om, Q, I), dt * ki3(i, om, Q, I)));

for i = 2 : n
    Q(i) = Q(i-1) + dt/6 * (kq1(i-1, I) + 2*kq2(i-1, omega_v1, Q, I) + 2*kq3(i-1, omega_v1, Q, I) + kq4(i-1, omega_v1, Q, I));
    Q2(i) = Q2(i-1) + dt/6 * (kq1(i-1, I2) + 2*kq2(i-1, omega_v2, Q2, I2) + 2*kq3(i-1, omega_v2, Q2, I2) + kq4(i-1, omega_v2, Q2, I2));
    Q3(i) = Q3(i-1) + dt/6 * (kq1(i-1, I3) + 2*kq2(i-1, omega_v3, Q3, I3) + 2*kq3(i-1, omega_v3, Q3, I3) + kq4(i-1, omega_v3, Q3, I3));
    Q4(i) = Q4(i-1) + dt/6 * (kq1(i-1, I4) + 2*kq2(i-1, omega_v4, Q4, I4) + 2*kq3(i-1, omega_v4, Q4, I4) + kq4(i-1, omega_v4, Q4, I4));
    
    I(i) = I(i-1) + dt/6 * (ki1(i-1, omega_v1, Q, I) + 2*ki2(i-1, omega_v1, Q, I) + 2*ki3(i-1, omega_v1, Q, I) + ki4(i-1, omega_v1, Q, I));
    I2(i) = I2(i-1) + dt/6 * (ki1(i-1, omega_v2, Q2, I2) + 2*ki2(i-1, omega_v2, Q2, I2) + 2*ki3(i-1, omega_v2, Q2, I2) + ki4(i-1, omega_v2, Q2, I2));
    I3(i) = I3(i-1) + dt/6 * (ki1(i-1, omega_v3, Q3, I3) + 2*ki2(i-1, omega_v3, Q3, I3) + 2*ki3(i-1, omega_v3, Q3, I3) + ki4(i-1, omega_v3, Q3, I3));
    I4(i) = I4(i-1) + dt/6 * (ki1(i-1, omega_v4, Q4, I4) + 2*ki2(i-1, omega_v4, Q4, I4) + 2*ki3(i-1, omega_v4, Q4, I4) + ki4(i-1, omega_v4, Q4, I4));
end

t = linspace(tmin, tmax, n);
figure;
hold on;
plot(t, Q); title('RRZ 2. rzedu dla obwodu RLC (metoda RK4)', 'FontSize', 14);
plot(t, Q2);
plot(t, Q3);
plot(t, Q4);
hold off;
axis([0 0.2516 -2.0e-3 3.3e-3])
xlabel('t', 'FontSize', 14);
ylabel('Q(t)', 'FontSize', 14);
legend('\omega_v=0.5\omega_0', '\omega_v=0.8\omega_0', '\omega_v=\omega_0', '\omega_v=1.2\omega_0', 'FontSize', 12);

figure;
hold on;
plot(t, I); title('RRZ 2. rzedu dla obwodu RLC (metoda RK4)', 'FontSize', 14);
plot(t, I2);
plot(t, I3);
plot(t, I4);
hold off;
axis([0 0.2516 -0.11 0.12])
xlabel('t', 'FontSize', 14);
ylabel('I(t)', 'FontSize', 14);
legend('\omega_v=0.5\omega_0', '\omega_v=0.8\omega_0', '\omega_v=\omega_0', '\omega_v=1.2\omega_0', 'FontSize', 12);

T1 = table(t', I', Q');
T2 = table(t', I2', Q2');
T3 = table(t', I3', Q3');
T4 = table(t', I4', Q4');
writetable(T1, 'rlcd1.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T2, 'rlcd2.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T3, 'rlcd3.txt', 'Delimiter',' ', 'WriteVariableNames', false);
writetable(T4, 'rlcd4.txt', 'Delimiter',' ', 'WriteVariableNames', false);
