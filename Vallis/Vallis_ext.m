function f = Vallis_ext(t,X)
global A B C U T_deep T dx omega 
f = zeros(12,1);
x = X(1);
y = X(2);
z = X(3);
Y= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];
f(1) = B*(z-y)/2/dx-C*(x-U);
f(2) = x*(T_deep-z)/2/dx-A*(y-T);
f(3) = x*(y-T_deep)/2/dx-A*(z-T);
Jac = [-C, -B/2/dx, B/2/dx;...
       (T_deep-z)/2/dx, -A, -x/2/dx;...
       (y-T_deep)/2/dx, x/2/dx, -A];
f(4:12)=Jac*Y;
