syms P(x,y)
P(x,y)= [1,x,y]

x1=0;
y1=0;
x2=0;
y2=1;
x3=1;
y3=0;

X=ones(3,3);
X(1,2)=x1;
X(1,3)=y1;
X(2,2)=x2;
X(2,3)=y2;
X(3,2)=x3;
X(3,3)=y3;

A=inv(X)

N=P*A
L=formula(N);

N1=subs(L(1,1))
N2=subs(L(1,2))
N3=subs(L(1,3))

% Z=N1
% ezmesh(N1)
% fsurf(x,y,N1) 
% figure
hold on

% ezplot(N2)
% plot()
% ezplot(N3)

% 
% figure
% O=[0.5:0.1:1]
% N2
% N3
% N1
% double(subs(N2,y,{O}))
% [X,Y]=meshgrid(double(subs(N3,[x],{[0:0.1:1]})),double(subs(N2,y,{[0:0.1:1]})));
% X
% Y
% % % 
% Z=double(subs(N1,[x,y],{X,Y}))
% % Z(Z>0.5)=0
% % % ezmesh(Z)
% % % % ezplot(Z)
% % % Z=1-X-Y
% 
% % surfc(X,Y,Z)   
% surfc(Z) 
% figure
% fill3(double(subs(N3,[x],{[0:0.1:1]})),double(subs(N2,[y],{[0:0.1:1]})),double(subs(N1,[x,y],{[0:0.1:1],[0:0.1:1]})),'r')
% [X,Y]=meshgrid(2:8,2:8,1:6);
% Z=X+Y;
% ezsurf(X,Y,Z);
% interval=[0 1 0 1 0 1];
% Z=double(N1,x,1);  
% fimplicit3(Z,interval)
% hold on
% N1+' ='+N2
% solve(eq(N1,N2))
fill3([x1,x2,x3],[y1,y2,y3],[double(subs(N1,[x,y],[x1,y1])),double(subs(N1,[x,y],[x2,y2])),double(subs(N1,[x,y],[x3,y3]))], 'r')
fill3([x1,x2,x3],[y1,y2,y3],[double(subs(N2,[x,y],[x1,y1])),double(subs(N2,[x,y],[x2,y2])),double(subs(N2,[x,y],[x3,y3]))], 'b')
fill3([x1,x2,x3],[y1,y2,y3],[double(subs(N3,[x,y],[x1,y1])),double(subs(N3,[x,y],[x2,y2])),double(subs(N3,[x,y],[x3,y3]))], 'g')