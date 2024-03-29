%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolución de la ecuación de ondas unidimensional por el metodo de diferencias finitas en una cuerda no homogénea
clc
clear
rho1 = 1; % densidad en medio 1
mu1 = 1;  % tensión en medio 1

rho2 = 10; % densidad en medio 2
mu2 = 10;  % tensión en medio 2

N = 299;  % número de pasos en espacio
h = 4/(N+1);
L = 500;  % número de pasos en tiempo
dt = 2/L; % resolvemos en (0, T] en tiempo
x = (-2+h:h:2-h)'; % resolvemos en (-2, 2) en espacio

r1 = (x <= 0).*rho1 + (x > 0).*rho2;
r2 = (x < 0).*rho1 + (x >= 0).*rho2;
r = (r1 + r2)./2;

m1 = (x <= 0).*mu1 + (x > 0).*mu2;
m2 = (x < 0).*mu1 + (x >= 0).*mu2;
m = (m1 + m2)./2;


u0 = 0.6*exp(-80*(x + 0.675).^2); % posición inicial (velocidad inicial = 0)
plot([-2; x; 2], [0; u0; 0])
axis([-2 2 -1 1])
pause

e = (dt/h)^2*(m./r);
e1 = (dt/h)^2*(m1./r); 
e2 = (dt/h)^2*(m2./r);

A = sparse(1:N, 1:N, 2*(1 - e), N, N)+ ...
      sparse(1:N-1, 2:N, e2(1:N-1), N, N) + ...
        sparse(2:N,1:N-1, e1(2:N),N,N);

u1 = (A*u0)/2; % posición en t=dt

for j = 1:L-1
   u=A*u1 - u0;
plot([-2; x; 2], [0; u; 0])
axis([-2 2 -1 1])
pause(0.05)
u0 = u1;
u1 = u;
end