         close all
        clear all
        clc
        
        
        % Config:
        s = tf('s');
        
        optionss=bodeoptions;
        optionss.MagVisible='on';
        optionss.PhaseMatching='on';
        optionss.PhaseMatchingValue=-180;
        optionss.PhaseMatchingFreq=1;
        optionss.Grid='on';
        
        %Constantes (en metros):
        
        Qi = 8 * 0.001 / 60;  %Caudal cte de entrada (en m3/s)
        diam = 10.65 * 0.001; %diametro de la cañeria de salida
        l_chico = 0.1; %lado chico del tanque 
        l_grande = 0.4; % " grande del tanque
        h_tanque = 0.9; % altura del tanque
        a_salida= pi*(diam/2)^2; %area de salida
        g = 9.81;  % Gravedad
        h0=0.45; %eq

        u0=0.504 %CALCULADO EN EL ONENOTE
        
    
    %-----------------------------------
    
    orden = 2;
    x=sym('x',[orden 1],'real');
    u=sym('u','real');
    
    % Punto de equlibrio (x'=0)
    u_e = u0;
    x_e = [h0 ; v0];
    
    %vector de x punto
    f1 = x(2);
    f2 = ((Qi-(u*a_salida*(sqrt(g*2)*x(1)))) / ((l_chico)^2 + (((2*l_chico*((l_grande)-(l_chico)))*x(1))/h_tanque) + ((((l_grande)-(l_chico))/h_tanque)*x(1))^2));

    
    
    f = [f1;f2];
    
    %salida (Altura del agua)
    g = x(1);
    
    A = jacobian(f,x);
    %la funcion subs cambia las ocurrencias de {x(1),x(2),u} por {x_e(1),x_e(2),u_e}
    A = double(subs(A,{x(1),x(2),u},{x_e(1),x_e(2),u_e}));
    
    B = jacobian(f,u);
    B = double(subs(B,{x(1),x(2),u},{x_e(1),x_e(2),u_e}));
    
    C = jacobian(g,x);
    C = double(subs(C,{x(1),x(2),u},{x_e(1),x_e(2),u_e}));
    
    D = jacobian(g,u);
    D = double(subs(D,{x(1),x(2),u},{x_e(1),x_e(2),u_e}));
    
    % Trasnferencia de la Planta Linealizada
    P = zpk(ss(A,B,C,D));
    
    Avals=eig(A);
    
    %figure(); hold on
    %bode(P,optionss);




