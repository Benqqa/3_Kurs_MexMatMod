clear all, clc, close all
E = 2e+11;
h = 0.140;
b = 0.058;
t1 = 0.0049;
t2 = 0.0049;
J = (b.*(h.^3) - (b - t2).*((h - 2.*t1).^3))./12;
%J=0.000000372506
N_nodes = 11;
N_elements = 10;
q = 2;
l = 0.1;
Elements = zeros(N_elements, 4);
for i = 1:N_elements
for k = 1:4
Elements(i,k) = 2.*i + k - 2;
end
end
ke = ((E.*J)./(l.^3))*[12, 6.*l, -12, 6.*l;
6.*l, 4.*(l.^2), -6.*l, 2.*(l.^2);
-12, -6.*l, 12, -6.*l;
6.*l, 2.*(l.^2), -6.*l, 4.*(l.^2)];
F = zeros(q.*N_nodes, 1);
F(22) = 10000;
K = zeros(q.*N_nodes, q.*N_nodes);
for i = 1:N_elements
K(Elements(i,1), Elements(i,1)) =K(Elements(i,1), Elements(i,1)) + ke(1, 1);
K(Elements(i,2), Elements(i,1)) =K(Elements(i,2), Elements(i,1)) + ke(2, 1);
K(Elements(i,2), Elements(i,2)) =K(Elements(i,2), Elements(i,2)) + ke(2, 2);
K(Elements(i,3), Elements(i,1)) =K(Elements(i,3), Elements(i,1)) + ke(3, 1);
K(Elements(i,3), Elements(i,2)) =K(Elements(i,3), Elements(i,2)) + ke(3, 2);
K(Elements(i,3), Elements(i,3)) =K(Elements(i,3), Elements(i,3)) + ke(3, 3);
K(Elements(i,4), Elements(i,1)) =K(Elements(i,4), Elements(i,1)) + ke(4, 1);
K(Elements(i,4), Elements(i,2)) =K(Elements(i,4), Elements(i,2)) + ke(4, 2);
K(Elements(i,4), Elements(i,3)) =K(Elements(i,4), Elements(i,3)) + ke(4, 3);
K(Elements(i,4), Elements(i,4)) =K(Elements(i,4), Elements(i,4)) + ke(4, 4);
end
K = K + K' - diag(diag(K));
K(1, :) = 0; K(:, 1) = 0; K(1, 1) = 1;
K(2, :) = 0; K(:, 2) = 0; K(2, 2) = 1;
%K(9, :) = 0; K(:, 9) = 0; K(9, 9) = 1;
U = linsolve(K,F);
Deformation = zeros(N_nodes, 1);
for i = 1:N_nodes
Deformation(i) = U(2.*i - 1);
end
Deformation
dx = 0.015625;
x1 = [0 : dx : l];
x2 = [l : dx : 2.*l];
x3 = [2.*l : dx : 3.*l];
x4 = [3.*l : dx : 4.*l];
N1 =@(p) 0.25.*((1 - p).^2).*(2 + p);
N2 =@(p) (l./8).*((1 - p).^2).*(1 + p);
N3 =@(p) 0.25.*((1 + p).^2).*(2 - p);
N4 =@(p) -(l./8).*((1 + p).^2).*(1 - p);
u1 = N1((2./l).*x1 - 1).*U(1) + N2((2./l).*x1 - 1).*U(2)+ N3((2./l).*x1 - 1).*U(3) + N4((2./l).*x1 - 1).*U(4);
u2 = N1((2./l).*x1 - 1).*U(3) + N2((2./l).*x1 - 1).*U(4)+ N3((2./l).*x1 - 1).*U(5) + N4((2./l).*x1 - 1).*U(6);
u3 = N1((2./l).*x1 - 1).*U(5) + N2((2./l).*x1 - 1).*U(6)+ N3((2./l).*x1 - 1).*U(7) + N4((2./l).*x1 - 1).*U(8);
u4 = N1((2./l).*x1 - 1).*U(7) + N2((2./l).*x1 - 1).*U(8)+ N3((2./l).*x1 - 1).*U(9) + N4((2./l).*x1 - 1).*U(10);
figure
plot([0 : l : 1],Deformation)
figure
hold on
plot(x1,u1)
plot(x1,u2)
plot(x1,u3)
plot(x1,u4)