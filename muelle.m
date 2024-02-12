%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empleo del método de Newton para resolver el problema no lineal de la
% ecuación del péndulo -u''= sin(x) en el intervalo (0,2*pi)

clear

eI = 0.0; %extremo izquierdo del intervalo
eD = 2*pi; %extremo derecho del intervalo
N  = 2^8; % N+1 subintervalos
h  = (eD - eI)/(N+1);
x  = (eI:h:eD)';

xInt = x(2:N+1);

% La matriz
E = sparse(2:N, 1:N-1, ones(N-1,1), N, N);
A = 2*speye(N,N) - E - E';
A = A/h^2; %Matriz

% uCero = 0.7*cos(xInt); %Dato inicial (a)
uCero = 0.7 + sin(xInt/2); %Dato inicial (b)

delta = uCero;
iter = 0;

while max(abs(delta))/max(abs(uCero)) > 10e-6
    
    iter = iter +1 % contador de iteraciones

    DG = A - sparse(1:N, 1:N, cos(uCero), N, N); % matriz del sistema de Newton
    
    % término independiente del sistema de Newton
    b = A*uCero - sin(uCero);
    b(1) = b(1) - 0.7/h^2;
    b(N) = b(N) - 0.7/h^2;
    
    
    delta = -DG\b; % resolución del sistema de Newton

    uOut = uCero + delta; % cálculo del nuevo iterando 

    plot(x, [0.7; uOut; 0.7])
    %axis([0 2*pi -1 1]) %Ejes para el caso (a)
    axis([0 2*pi -1 12]) %Ejes para el caso (b)
    hold on
    pause(2)
    uCero = uOut;
    if iter > 15; break 
    end
end