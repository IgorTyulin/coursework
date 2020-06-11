clear all;
clc;

nano_multiplier = 10^9; %to transform plot x axis into acceptable form

a = 10; % aperture radius
c = 3e+8; % speed of light, m/sec
one_t_h = 9*pi/180;
    
lambda = a/10;
N_FFT = 8192*2;
T = 100/c*lambda;
d_t = T/(N_FFT - 1);

t_h = (0:0.2:10)*pi/180; % Theta angle
t = -T/2:d_t:T/2;
t_demonstrate = -5.5*(10^-9):0.01*(10^-9):5.5*(10^-9); % For reasonable plot display

%Ищем значения для формулы с подставленной единицей и рисуем графики
%calculation
% normal
t_actual = t_demonstrate;
ro = one_t_h;
f_norm = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_norm(i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2);

% falling n = 1
t_actual = t_demonstrate;
ro = one_t_h;
n = 1;
f_fal = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_fal(i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;

% falling n = 2
t_actual = t_demonstrate;
ro = one_t_h;
n = 2;
f_fal_2 = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_fal_2(i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;

% 1st plot (norm, falling n=1, falling n=2)
% Transient far fields of circular flat aperture for different field amplitude distributions over the amplitude
figure
plot(t_demonstrate*nano_multiplier, f_norm, 'k')
hold on
plot(t_demonstrate*nano_multiplier, f_fal, '--k')
plot(t_demonstrate*nano_multiplier, f_fal_2, '-.k')
grid on
hold off
title('Transient far fields of circular flat aperture. a = 10 m. \theta = 9^o')
xlabel('Time, ns')
ylabel('Amplitude')
legend('Uniform distribution', 'Falling distribution. n=1', 'Falling distribution. n=2')

% 2nd plot (norm, falling n=1, falling n=2)
% Derivatives of transient far fields or three electric field distributions
figure
n = size(t_demonstrate,2);

% norm
f_der = diff(f_norm) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier, f_der, 'k')

hold on

% falling 1
f_der = diff(f_fal) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier,  f_der, '--k')

%falling 2
f_der = diff(f_fal_2) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier,  f_der, '-.k')

grid on
hold off    
title('Derivatives of transient far fields. a = 10 m. \theta = 9^o')
xlabel('Time, ns')
ylabel('Normilized Amplitude')
legend('Uniform distribution', 'Falling distribution. n=1', 'Falling distribution. n=2')


% 3rd plot (norm, falling n=1, falling n=2)
% Far field antenna patterns for monochromatic signal
for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;  
    
    % norm
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2);
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);
    
    % falling 1
    n = 1;
    E_f_fal_1(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_fal_1(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_fal_1(j, :) = fft(E_f_fal_1(j, :), N_FFT);
    
    % falling 2
    n = 2;
    E_f_fal_2(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_fal_2(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_fal_2(j, :) = fft(E_f_fal_2(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_fal_1(:, p))/max(abs(E_f_fal_1(:, p)))), '--k')
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_fal_2(:, p))/max(abs(E_f_fal_2(:, p)))), '-.k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('Uniform distribution', 'Falling distribution. n=1', 'Falling distribution. n=2')


%Подставляем в формулу cos (вместо cos^2) и рисуем графики    
%calculation
t_actual = t_demonstrate;
ro = one_t_h;

% normal
f_fal_3 = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_fal_3(i) = (cos(ro)/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2);

% falling n = 1
n = 1;
f_fal_4 = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_fal_4(i) = (cos(ro)/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;

% falling n = 2
n = 2;
f_fal_5 = zeros(1, length(t_actual));
i = find(abs(c*t_actual) < a * sin(ro));
f_fal_5(i) = (cos(ro)/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;

% 1st plot (falling n=0, falling n=1, falling n=2)
% Transient far fields of circular flat aperture for different field amplitude distributions over the amplitude
figure
plot(t_demonstrate*nano_multiplier, f_fal_3, 'k')
hold on
plot(t_demonstrate*nano_multiplier, f_fal_4, '--k')
plot(t_demonstrate*nano_multiplier, f_fal_5, '-.k')
grid on
hold off
title('Transient far fields of circular flat aperture. a = 10 m. \theta = 9^o')
xlabel('Time, ns')
ylabel('Amplitude')
legend('Falling distribution.', 'Falling distribution. n=1', 'Falling distribution. n=2')

% 2nd plot (falling n=0, falling n=4, falling n=5)
% Derivatives of transient far fields or three electric field distributions
figure
n = size(t_demonstrate, 2);

% falling n = 0
f_der = diff(f_fal_3) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier,  f_der, 'k')

hold on

% falling n = 1
f_der = diff(f_fal_4) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier,  f_der, '--k')

% falling n = 2
f_der = diff(f_fal_5) ./ diff(t_demonstrate);
f_der = f_der ./ norm(f_der, n);
plot(t_demonstrate(1: end - 1)*nano_multiplier,  f_der, '-.k')

grid on
hold off    
title('Derivatives of transient far fields. a = 10 m. \theta = 9^o')
xlabel('Time, ns')
ylabel('Normilized Amplitude')
legend('Falling distribution. ', 'Falling distribution. n=1', 'Falling distribution. n=2')

% 3rd plot (norm, falling n=1, falling n=2)
% Far field antenna patterns for monochromatic signal
for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % norm
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2);
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
    
    % falling n=1
    n = 1;
    E_cos_fal_1(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_fal_1(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_fal_1(j, :) = fft(E_cos_fal_1(j, :), N_FFT);
    
    % falling n=2
    n = 2;
    E_cos_fal_2(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_fal_2(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_fal_2(j, :) = fft(E_cos_fal_2(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_fal_1(:, p))/max(abs(E_cos_fal_1(:, p)))), '--k')
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_fal_2(:, p))/max(abs(E_cos_fal_2(:, p)))), '-.k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('Falling distribution.', 'Falling distribution. n=1', 'Falling distribution. n=2')    


%Comparison of graphs between 1 and cos
for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % norm with 1
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_f_norm(j, i) = (1/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2);
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    %norm with cos
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_cos_norm(j, i) = (cos(ro)/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2);
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 0')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')

for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % falling 1 with 1
    n = 1;
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    % falling 1 with cos
    n = 1;
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 1')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')

for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % falling 2 with 1
    n = 2;
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    % falling 2 with cos
    n = 2;
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_fal_2(:, p))/max(abs(E_f_fal_2(:, p)))), '-.k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 2')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')

for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % falling 3 with 1
    n = 3;
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    % falling 3 with cos
    n = 3;
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 3')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')

for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % falling 4 with 1
    n = 4;
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    % falling 4 with cos
    n = 4;
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 4')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')

for j = 1 : size(t_h, 2)
    i_t_h = t_h(j);
    t_actual = t/2;
    ro = i_t_h;
    
    % falling 5 with 1
    n = 5;
    E_f_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro));
    E_f_norm(j, i) = (1/(pi*sin(ro)^2))*sqrt((a*sin(ro))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro))^2))).^n;
    E_f_norm(j, :) = fft(E_f_norm(j, :), N_FFT);

    % falling 5 with cos
    n = 5;
    E_cos_norm(j, :) = zeros(1, length(t_actual));
    i = find(abs(c*t_actual) < a * sin(ro+0.01));
    E_cos_norm(j, i) = (cos(ro+0.01)/(pi*sin(ro+0.01)^2))*sqrt((a*sin(ro+0.01))^2-(c*t_actual(i)).^2) .* (1-(((c*t_actual(i)).^2)/((a*sin(ro+0.01))^2))).^n;
    E_cos_norm(j, :) = fft(E_cos_norm(j, :), N_FFT);
end

d_f = 1/N_FFT/d_t;
k = 1:1:N_FFT/2;
[F, p] = min(abs(d_f*k - 3e+8/lambda));
figure
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_f_norm(:, p))/max(abs(E_f_norm(:, p)))), 'k')
hold on
plot(4*a/lambda*sin(t_h), 20*log10(abs(E_cos_norm(:, p))/max(abs(E_cos_norm(:, p)))), '--k')
grid on
hold off
title('Far field antenna patterns for monochromatic signal. a/\lambda = 10. n = 5')
xlabel('2*a/\lambda*sin(\Theta)')
ylabel('Amplitude, dB')
legend('With 1', 'With cos')