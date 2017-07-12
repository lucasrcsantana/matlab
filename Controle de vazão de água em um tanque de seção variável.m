% DEFINI��ES
esp = 9.5;
R1= 0.32;
R2=0.32;
R3= 0.32;
Rac= 0.5;
R = 9.5;

for H = 1:21
fprintf('LOOP: %d \n',H)

 syms h3 q3 x3 S A3;   
    % TANQUE 3
    
    if H < 3    
        A3 = 6*esp;
    else
        x3 = sqrt(R^2-(R-(h3-2.1))^2);
        A3=2*x3*esp;
    end
    
    F3 = q3/A3 - sqrt(h3)/(R3*A3); 

    %Matrizes
    a11 = diff(F3,h3);
    b11 = diff(F3,q3);  
    
    h3 = H;

    %Pontos de opera��o
    q3 = sqrt(h3)/R3;

    %SUBS da Matriz em x3
    
    x3 = subs(x3,sym('h3'),h3);

    %SUBS da Matriz A2 a A3 
    A3 = subs(A3,[sym('x3'),sym('h3'),sym('h3')],[x3,h3,q3]);

    %SUBS da Matriz A
    a11 = subs(a11,[sym('x3'),sym('h3'),sym('q3')],[x3,h3,q3]);
    if H < 3
        a11= eval(a11);
    end

    %SUBS da Matriz B
    b11 = subs(b11,[sym('x3'),sym('h3'),sym('q3')],[x3,h3,q3]);

    %�rea do tanque 3
    A3 ;

    %Matriz A, B, C e D
    A = a11;
    B = b11;
    C = 1;
    D = 0;

    %Calculo de G(s) e T(s)
    G=C*(inv(((S*eye(1))-A))*B)+D; 
    G = vpa(G,3);
    %G = eval(G,'S')
    eig(G);
  

    [NG1,DG1] = ss2tf(A,B,C,D,1);
    G11 = tf(NG1,DG1)

    [z1,p1,kp1] = tf2zp(NG1,DG1); % encontrar polos e zeros

    Z1 = -1.05*p1; %determina��o de zeros do controlador
    
    kpc = (1/50)/(kp1); %retirado do root locus
    
    fprintf('Controlador e processo sem ganho Kp');
    G11C=G11*tf([1 Z1],[1 0])
    
   
%     figure('Name','Lugar das Ra�zes')
%     subplot(2,1,1)
%     rlocus(G11)   % LR sem controlador
%     title('Lugar das ra�zes sem controlador')
%     subplot(2,1,2)
%     rlocus(G11C)  % LR com controlador
%     title('Lugar das ra�zes com controlador')
%    
%     figure('Name','Resposta ao degrau unit�rio')
%     subplot(3,1,1)   
%     step(feedback(G11,1))   
%     title('Resposta ao degrau sem controlador')
%     xlabel('t [s]') % x-axis label
%     ylabel('Ganho') % y-axis label
%     axis([0,900,0,1])
%     subplot(3,1,2) 
%     step(feedback(G11C,1))  
%     title('Resposta ao degrau com controlador sem ganho KP')
%     xlabel('t [s]') % x-axis label
%     ylabel('Ganho') % y-axis label
%     axis([0,900,0,1])
%     subplot(3,1,3) 
%     step(feedback(kpc*G11C,1))  
%     title('Resposta ao degrau com controlador com ganho KP')
%     xlabel('t [s]') % x-axis label
%     ylabel('Ganho') % y-axis label
%     axis([0,900,0,1])
    
    fprintf('Controlador e processo com ganho Kp')
    G11C= kpc*G11*tf([1 Z1],[1 0])
    
    vazao(H)=q3;
    kp(H) = kp1;
    z1c(H) = Z1;
    kp(H) = kpc;
    ki(H)= kpc*Z1;
    %ki(H)=(h3-b11)/a11;
    areas(H) = A3;
    polos(H) = DG1(2);
    zeros(H) = NG1(2);
    altura(H) = H;
end

% GR�FICOS

figure('Name','Vaz�o')
plot(altura)
title('Vaz�o [cm�/s] x Altura')
axis([-inf,inf,-inf,inf])

figure('Name','Vaz�o')
subplot(2,1,1)
plot(vazao)
title('Vaz�o [cm�/s] x Altura')
xlabel('H [cm]') % x-axis label
ylabel('Vaz�o [cm�/s]') % y-axis label
axis([-inf,inf,-inf,inf])
subplot(2,1,2)
plot(areas,vazao)
title('Vaz�o [cm�/s] x �rea')
xlabel('�reas [cm�]') % x-axis label
ylabel('Vaz�o [cm�/s]') % y-axis label
axis([-inf,inf,-inf,inf])

figure('Name','�reas')
plot(areas)
title('�reas do Tanque 3 x Altura')
xlabel('H [cm]') % x-axis label
ylabel('�reas [cm�]') % y-axis label
axis([-inf,inf,-inf,inf])

figure('Name','Ganhos de Ki')
subplot(2,1,1)
plot(ki)
title('Ganhos de Ki x Altura')
xlabel('H [cm]') % x-axis label
ylabel('Ganhos') % y-axis label
axis([-inf,inf,-inf,inf])
subplot(2,1,2)
plot(areas,ki)
title('Ganhos de Ki x �rea')
xlabel('�reas [cm�]') % x-axis label
ylabel('Ganhos') % y-axis label
axis([-inf,inf,-inf,inf])

figure('Name','Ganhos de KP')
subplot(2,1,1)
plot(kp)
title('Ganhos de Kp x Altura')
xlabel('H [cm]') % x-axis label
ylabel('Ganhos') % y-axis label
axis([-inf,inf,-inf,inf])
subplot(2,1,2)
plot(areas,kp)
title('Ganhos de Kp x �rea')
xlabel('�reas [cm�]') % x-axis label
ylabel('Ganhos') % y-axis label
axis([-inf,inf,-inf,inf])

figure('Name','P�los')
subplot(2,1,1)
plot(polos)
title('P�los da planta x Altura')
xlabel('H [cm]') % x-axis label
ylabel('P�los') % y-axis label
axis([-inf,inf,-inf,inf])
subplot(2,1,2)
plot(areas,polos)
title('P�los da planta x �rea')
xlabel('�reas [cm�]') % x-axis label
ylabel('P�los') % y-axis label
axis([-inf,inf,-inf,inf])

figure('Name','Zeros')
subplot(2,1,1)
plot(zeros)
title('Zeros do Controlador x Altura')
xlabel('H [cm]') % x-axis label
ylabel('Zeros') % y-axis label
axis([-inf,inf,-inf,inf])
subplot(2,1,2)
plot(areas,zeros)
title('Zeros do Controlador x �rea')
xlabel('�reas [cm�]') % x-axis label
ylabel('Zeros') % y-axis label
axis([-inf,inf,-inf,inf])

fprintf('FIM DO PROGRAMA %d \n\n')