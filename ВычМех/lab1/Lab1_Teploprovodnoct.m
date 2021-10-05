T_x_0 = @(x) cos(2*x)
T_0_t = @(t) 1-6*t
T_06_t= @(t) 0.3624

dt=0.001
h=0.1
lambda = 1
p=1
C_v=1

x_0=0
x_end=0.6
t_0=0
t_end=0.1
dx=t_end/6

X=[x_0:h:x_end]
Time=[t_0:dt:t_end]
T=zeros(length(X),length(Time))
%нач заполнение
T(1,:)=T_0_t(Time)
T(length(X),:)=T_06_t(Time)
T(:,1)=T_x_0(X)
T1=ne_yav_metod(T,X,Time,h,dt,lambda,p,C_v)
T2=yav_metod(T,X,Time,h,dt,lambda,p,C_v)
figure
hold on
xlabel('t') 
ylabel('x') 
zlabel('T') 
surf(T1)
surf(T2)
colorbar

figure
hold on
xlabel('t') 
ylabel('x') 
zlabel('T') 
surf(T2)
colorbar

figure
hold on
xlabel('t') 
ylabel('x') 
zlabel('T') 
surf(T1)
colorbar

function [T]=ne_yav_metod(T,X,Time,h,dt,lambda,p,C_v)
    for k=1:length(Time)-1
        for i=2:length(X)-1
            T(i,k+1)=((T(i+1,k)+T(i-1,k))*lambda/h^2+(p*C_v/dt-2*lambda/h^2)*T(i,k))*dt/(p*C_v);
        end
    end
end
function [T]=yav_metod(T,X,Time,h,dt,lambda,p,C_v)
    X_l=length(X);
    T_l=length(Time);
    A=lambda/h^2;
    C=lambda/h^2;
    B=(p*C_v*h^2+2*lambda*dt)/(dt*h^2);
    for k=1:T_l-1
        for i=1:X_l
            F(i,k)=p*C_v*T(i,k)/dt;
        end
        P(1)=C/B;
        Q(1)=F(1,k)/B;
        %прямой ход
        for i=2:X_l
           P(i)=C/(B-A*P(i-1));
           Q(i)=(F(i,k)+A*Q(i-1))/(B-A*P(i-1));
        end
        %обратный ход
        for i=X_l-1:-1:2
            T(i,k+1)=P(i)*T(i+1,k+1)+Q(i);
        end
    end
end

