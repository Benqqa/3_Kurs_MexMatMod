%suka
%����
syms alfa fi psi;
%������� ��������
syms w_alfa w_fi w_psi;
%���������
syms R m N v_c dv_c g mu f t;
%��������� ������
% j=[1,0,0];
% i=[0,1,0];
i=[1,0,0];
j=[0,1,0];
k=[0,0,1];
%������� ��������� �� ����
P_i=create_P(fi, i);
P_j=create_P(alfa, j);
P_k=create_P(psi, k);
%�������� ������ ��������
P=P_k*P_j*P_i
%������ ������
R_alfa_psi=(P_k*P_j*k')'
%������� ��������
w=w_psi*k+(P_k*w_alfa*j')'+(P_k*P_j*w_fi*i')'
%������ ���������
Teta=(m*R^2)/2*kron(i,i')+(m*R^2)/4*(kron(j,j')+kron(k,k'))
%�������� �����
v_b=v_c+cross(w,R_alfa_psi)
%����� ����� �� ��������
if(v_b ~= 0)
    F_tr=-mu*N*v_b/norm(v_b)
else
    F_tr=f
end
%2� ���� �������
dv_c=(F_tr+N*R_alfa_psi-m*g*k)/m
%������ (B)'=A
B=(P*Teta*P'*w')'
dB=int(B,t)
A=cross(R_alfa_psi,F_tr)
fi=A\dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������ v_c
v_c=int(dv_c,t)

function [P]=create_P(fi, v)
    P=zeros(3,3);
    E=eye(3);
    P=kron(v,v')+cos(fi)*(E-kron(v,v'))+sin(fi)*v_x_M(v,E);
end
function [T]=v_x_M(v,M)
    %��������� ������������� ������� �� �������
    T=[cross(v,M(1,:));cross(v,M(2,:));cross(v,M(3,:))];
end


