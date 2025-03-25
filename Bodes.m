% BODE CON DISTINTOS H Y U

close all
clear all
clc

% Configuración:
s = tf('s');

optionss = bodeoptions;
optionss.MagVisible = 'on';
optionss.PhaseMatching = 'on';
optionss.PhaseMatchingValue = -180;
optionss.PhaseMatchingFreq = 1;
optionss.Grid = 'on';

% Constantes (en metros):
Qi = 8 * 0.001 / 60;  % Caudal constante de entrada (m³/s)
diam = 10.65 * 0.001; % Diámetro de la cañería de salida (m)
l_chico = 0.1; % Lado chico del tanque 
l_grande = 0.4; % Lado grande del tanque
h_tanque = 0.9; % Altura del tanque
a_salida = pi * (diam / 2)^2; % Área de salida
g = 9.81;  % Gravedad

% Vector de alturas
h_vec = [0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80];

% Calcular u para cada h en h_vec
u_vec = Qi ./ (a_salida .* sqrt(2 * g * h_vec));


figure(); hold on
title('Diagrama de Bode para distintos pares de h y u')

% Iterar sobre cada par (h, u)
for i = 1:length(h_vec)
    h0 = h_vec(i);
    u0 = u_vec(i);

    orden = 2;
    x = sym('x', [orden 1], 'real');
    u = sym('u', 'real');

    % Punto de equilibrio (x' = 0)
    u_e = u0;
    x_e = [h0; 0];

    % Definir ecuaciones
    f1 = x(2);
    f2 = ((Qi - (u * a_salida * sqrt(2 * g * x(1)))) / ...
        ((l_chico)^2 + (((2 * l_chico * ((l_grande) - (l_chico))) * x(1)) / h_tanque) + ...
        ((((l_grande) - (l_chico)) / h_tanque) * x(1))^2));
    

    f = [f1; f2];

    % Salida (altura del agua)
     y = x(1);

    % Linealización
    A = jacobian(f, x);
    A = double(subs(A, {x(1), x(2), u}, {x_e(1), x_e(2), u_e}));

    B = jacobian(f, u);
    B = double(subs(B, {x(1), x(2), u}, {x_e(1), x_e(2), u_e}));

    C = jacobian( y, x);
    C = double(subs(C, {x(1), x(2), u}, {x_e(1), x_e(2), u_e}));

    D = jacobian( y, u);
    D = double(subs(D, {x(1), x(2), u}, {x_e(1), x_e(2), u_e}));

    % Función de transferencia de la planta linealizada
    P = tf(ss(A, B, C, D));

    % Graficar Bode
    bode(P, optionss);
end

legend(arrayfun(@(h) sprintf('h = %.2f m', h), h_vec, 'UniformOutput', false))
hold off
