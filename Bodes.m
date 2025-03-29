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

% Iterar sobre cada par (h, u)
for i = 1:length(h_vec)
    h0 = h_vec(i);
    u0 = u_vec(i);

    orden = 1;
    x=sym('x',[orden 1],'real');
    u=sym('u','real');
    
    % Punto de equlibrio (x'=0)
    u_e = u0;
    x_e = h0;
    
    %x punto
    f = ((Qi - (u * a_salida * sqrt(2 * g * x))) / ...
        ((l_chico)^2 + (((2 * l_chico * ((l_grande) - (l_chico))) * x) / h_tanque) + ...
        ((((l_grande) - (l_chico)) / h_tanque) * x)^2));
    
    
    
    %salida (Altura del agua)
    y = x;
    
    A = jacobian(f,x);
    %la funcion subs cambia las ocurrencias de {x,u} por {x_e,u_e}
    A = double(subs(A,{x,u},{x_e,u_e}));
    
    B = jacobian(f,u);
    B = double(subs(B,{x,u},{x_e,u_e}));
    
    C = jacobian(y,x);
    C = double(subs(C,{x,u},{x_e,u_e}));
    
    D = jacobian(y,u);
    D = double(subs(D,{x,u},{x_e,u_e}));
    
    % Trasnferencia de la Planta Linealizada
    P = tf(ss(A,B,C,D))
    
    Avals=eig(A)
    
    bode(P,optionss);

end

legend(arrayfun(@(h) sprintf('h = %.2f m', h), h_vec, 'UniformOutput', false))
hold off
