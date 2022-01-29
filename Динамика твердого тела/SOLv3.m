clc
%углы
syms t alfa(t) fii(t) psii(t) vx(t) vy(t) vz(t) w1(t) w2(t) w3(t);
%углова€ скорость
w_alfa=diff(alfa,t)
w_fii =diff(fii,t)
w_psii =diff(psii,t)
%единичные ветора
i=[1,0,0];
j=[0,1,0];
k=[0,0,1];
%пам-параметры
R=10;
m=1;
g=10;
mu=0.5;
N=m*g;
f=1
w0=[w1 w2 w3]
v_c=vx*i+vy*j+vz*k
dv_c=diff(v_c,t)
%тензора поаворота по ос€м
P_i=create_P(fii, i);
P_j=create_P(alfa, j);
P_k=create_P(psii, k);
%»тоговый тензор поворота
P=P_k*P_j*P_i;
%–адиус вектор
R_alfa_psii=(P_k*P_j*(-R*k.')).';
%углова€ скорость
w=w_psii*k+(P_k*w_alfa*j.').'+(P_k*P_j*w_fii*i.').';
%инерци€
Teta=(m*R^2)/2*kron(i,i.')+(m*R^2)/4*(kron(j,j.')+kron(k,k.'));
%скорость точки б
v_b=v_c+cross(w0,R_alfa_psii);
%силат рени€ от скорости
F_tr=-mu*N*v_b/norm(v_b)
%доп переменнные
B=(P*Teta*P.'*w0.')
A=cross(R_alfa_psii,F_tr).';
%%%%%%%%%%%%%%%%%%%%%%%%
%формируем основные уравнени€
mass_dif1 = diff(B,t) == -A
mass_dif2 = m*diff(v_c,t) == F_tr-m*g*k+N*R_alfa_psii
mass_dif3 = w0 == w
syms d_vx d_vy d_vz d_w1 d_w2 d_w3 d_alfa d_psii d_fii dd_alfa dd_psii dd_fii;
%приводим перевенные в пор€док
a1=toSym(mass_dif3.',t, alfa, fii, psii, vx, vy, vz, w1, w2, w3);
a2=toSym(mass_dif2.',t, alfa, fii, psii, vx, vy, vz, w1, w2, w3);
a3=toSym(mass_dif1,t, alfa, fii, psii, vx, vy, vz, w1, w2, w3);
%выражаем старшие производные
[res_d_alfa, res_d_psii, res_d_fii]=solve(a1,[d_alfa d_psii d_fii]);
[res_d_vx, res_d_vy, res_d_vz]=solve(a2,[d_vx d_vy d_vz]);
[res_d_w1, res_d_w2, res_d_w3]=solve(a3,[d_w1 d_w2 d_w3]);
%решаем исстему
t_interval=[0 1];
Start_cond=[15;10;20; 0;0;0; 0;0;30; 10;10;10;];
[t,X]=ode45(@(t,X)new_solver(t,X,res_d_alfa,res_d_psii,res_d_fii,res_d_vx,res_d_vy,res_d_vz,res_d_w1,res_d_w2,res_d_w3,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii),t_interval,Start_cond)
%вырезка
res_x=X(:,1:3)
res_fi=X(:,4:6)
res_v=X(:,7:9)
res_w=X(:,10:12)
%посмотрим что там получилось?0_0
figure
plot(t,X(:,1:3))
legend('x(t)','y(t)','z(t)')
figure
plot(t,X(:,4:6))
legend('alfa(t)','fii(t)','psii(t)')
figure
plot(t,X(:,7:9))
legend('vx(t)','vy(t)','vz(t)')
figure
plot(t,X(:,10:12))
legend('w1(t)','w2(t)','w3(t)')
function dXdt= new_solver(t,X,res_d_alfa,res_d_psii,res_d_fii,res_d_vx,res_d_vy,res_d_vz,res_d_w1,res_d_w2,res_d_w3,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii)
    t
    %что-то здесь не так... пока наблюдаем....
    dx1=double(X(7))
    dx2=double(X(8))
    dx3=double(X(9))
    dx4=double(SetX(res_d_alfa,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,1,1,1));
    dx5=double(SetX(res_d_fii,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,1,1,1));
    dx6=double(SetX(res_d_psii,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,1,1,1));
    dx7=double(SetX(res_d_vx,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    dx8=double(SetX(res_d_vy,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    dx9= double(SetX(res_d_vz,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    dx10=double(SetX(res_d_w1,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    dx11=double(SetX(res_d_w2,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    dx12=double(SetX(res_d_w3,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6));
    
dXdt=[double(dx1);double(dx2);double(dx3);double(dx4);double(dx5);double(dx6);double(dx7);double(dx8);double(dx9);double(dx10);double(dx11);double(dx12)]

end
function e1 = SetX(res,X,alfa, fii, psii, vx, vy, vz, w1, w2, w3, d_vx, d_vy, d_vz, d_w1, d_w2 ,d_w3 ,d_alfa ,d_psii ,d_fii ,dd_alfa ,dd_psii ,dd_fii,dx1,dx2,dx3,dx4,dx5,dx6)
    %из символа в число
    e1 = subs(res,alfa,X(4));
    e1 = subs(e1,fii,X(5));
    e1 = subs(e1,psii,X(6));
    e1 = subs(e1,vx,X(7));
    e1 = subs(e1,vy,X(8));
    e1 = subs(e1,vz,X(9));
    e1 = subs(e1,w1,X(10));
    e1 = subs(e1,w2,X(11));
    e1 = subs(e1,w3,X(12));
    e1 = subs(e1,d_alfa,dx4);
    e1 = subs(e1,d_fii,dx5);
    e1 = subs(e1,d_psii,dx6);
end
function a = toSym(mass,t,alfa, fii, psii, vx, vy, vz, w1, w2, w3)
    a=mass;
    syms d_vx d_vy d_vz d_w1 d_w2 d_w3 d_alfa d_psii d_fii dd_alfa dd_psii dd_fii;
    a=subs(a,diff(alfa(t), t),d_alfa);
    a=subs(a,diff(psii(t), t),d_psii);
    a=subs(a,diff(fii(t), t),d_fii);
    a=subs(a,diff(vx(t), t),d_vx);
    a=subs(a,diff(vy(t), t),d_vy);
    a=subs(a,diff(vz(t), t),d_vz);
    a=subs(a,diff(w1(t), t),d_w1);
    a=subs(a,diff(w2(t), t),d_w2);
    a=subs(a,diff(w3(t), t),d_w3);
    a=subs(a,diff(alfa(t), t,t),dd_alfa);
    a=subs(a,diff(psii(t), t,t),dd_psii);
    a=subs(a,diff(fii(t), t,t),dd_fii);
end
function [P]=create_P(fii, v)
    %объ€влем поворот
    P=zeros(3,3);
    E=eye(3);
    P=kron(v,v.')+cos(fii)*(E-kron(v,v.'))+sin(fii)*v_x_M(v,E);
end
function [T]=v_x_M(v,M)
    %векторное происзведение вектора на матрицу
    T=[cross(v,M(1,:));cross(v,M(2,:));cross(v,M(3,:))];
end