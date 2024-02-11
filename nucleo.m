function n = nucleo(p,q,inte,alfa,beta)
    %-----
    % Calcula el núcleo del operador integral asociado a una ecuación de
    % Sturm-Liouville con coeficientes constantes
    %-----
    % Entradas:
    %   p         = Función p de la ecuación de Sturm-Liouville (constante)
    %   q         = Función q de la ecuación de Sturm-Liouville (constante)
    %   inte      = Intervalo de definición de la EDO
    %   alfa      = vector con dos coordenadas con las condiciones en el
    %   extremo izquierdo.
    %   beta      = vector con dos coordenadas con las condiciones en el
    %   extremo derecho.
    %-----
    % Salidas
    %   n = la función inline con el núcleo
    %-----

    % defino variables simbólicas
    syms x t y(t)
    % y la expresión de la ecuación diferencial
    Dy(t) = diff(y(t),t);
    D2y(t) = diff(y(t),t,2);
    Expr(t) = p*D2y(t)+q*y(t);

    %   cálculo de u
    %   hay que resolver la EDO p*y''+q*y=0, con y(a)=-alfa(2); y'(a) = alfa(1)
    u(t) = dsolve(Expr(t)==0,[y(inte(1))==-alfa(2),Dy(inte(1))==alfa(1)]);
    u(t) = simplify(u(t));

    %-----
    %   cálculo de v
    %   hay que resolver la EDO p*y''+q*y=0, con y(b)=-beta(2); y'(b) = beta(1)
    v(t) = dsolve(Expr(t)==0,[y(inte(2))==-beta(2),Dy(inte(2))==beta(1)]);
    v(t) = simplify(v(t));

    %-----
    %   cálculo de J
    J(t) = det([u(t) v(t); diff(u(t),t) diff(v(t),t)]);

    % finalmente se monta el núcleo
    n1 = matlabFunction(u(t)*v(x)/J(t)/p);
    n2 = matlabFunction(v(t)*u(x)/J(t)/p);

    n = @(t,x) n1(t,x).*(inte(1)<=t & t<=x)+n2(t,x).*(x<t & t<=inte(2));
end