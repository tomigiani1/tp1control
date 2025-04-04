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

        u0= Qi/(a_salida*sqrt(2*g*h0)) 
        
    
    
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
    
    %figure(); hold on
    %bode(P,optionss);

    %Diseño el controlador. Lo quiero lo mas simple posible, por lo que diseñaré un PI

    k=db2mag(14.9);
    C= (-1*k)*(s+0.02)/s;

    
    
    L= minreal(P*C);

    %figure(); hold on
    %bode(L,optionss);