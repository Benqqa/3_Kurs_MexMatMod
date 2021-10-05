%suka
%углы
syms alfa fi psi;
%угловая скорость
syms w_alfa w_fi w_psi;
%параметры
syms R m N v_c dv_c g mu f t;
%единичные ветора
% j=[1,0,0];
% i=[0,1,0];
i=[1,0,0];
j=[0,1,0];
k=[0,0,1];
%тензора поаворота по осям
P_i=create_P(fi, i);
P_j=create_P(alfa, j);
P_k=create_P(psi, k);
%Итоговый тензор поворота
P=P_k*P_j*P_i
%Радиус вектор
R_alfa_psi=(P_k*P_j*k')'
%угловая скорость
w=w_psi*k+(P_k*w_alfa*j')'+(P_k*P_j*w_fi*i')'
%Тензор симметрии
Teta=(m*R^2)/2*kron(i,i')+(m*R^2)/4*(kron(j,j')+kron(k,k'))
%скорость точки
v_b=v_c+cross(w,R_alfa_psi)
%силат рения от скорости
if(v_b ~= 0)
    F_tr=-mu*N*v_b/norm(v_b)
else
    F_tr=f
end
%2й закн ньютона
dv_c=(F_tr+N*R_alfa_psi-m*g*k)/m
%Баланс (B)'=A
B=(P*Teta*P'*w')'
dB=int(B,t)
A=cross(R_alfa_psi,F_tr)
fi=A\dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Найдем v_c
v_c=int(dv_c,t)

function [P]=create_P(fi, v)
    P=zeros(3,3);
    E=eye(3);
    P=kron(v,v')+cos(fi)*(E-kron(v,v'))+sin(fi)*v_x_M(v,E);
end
function [T]=v_x_M(v,M)
    %векторное происзведение вектора на матрицу
    T=[cross(v,M(1,:));cross(v,M(2,:));cross(v,M(3,:))];
end


