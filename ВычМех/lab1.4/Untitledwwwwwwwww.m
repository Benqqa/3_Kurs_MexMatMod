
syms P(x,y,z)
P(x,y,z)= [1,x,y,z,x*z,y*z]

x1=0;
y1=0;
x2=0;
y2=1;
x3=1;
y3=0;
prisma=ones(6,3)
t1=[0,0,0];
t2=[0,1,0];
t3=[1,0,0];
t4=[0,0,1];
t5=[1,0,1];
t6=[0,1,1];
prisma(1,:)=t1;
prisma(2,:)=t2;
prisma(3,:)=t3;
prisma(4,:)=t4;
prisma(5,:)=t5;
prisma(6,:)=t6
X=ones(6,6);
k=1
for i =[1:6]
    X(k,1)=1;
    X(k,2)=prisma(i,1);
    X(k,3)=prisma(i,2);
    X(k,4)=prisma(i,3);
    X(k,5)=prisma(i,1)*prisma(i,3);
    X(k,6)=prisma(i,2)*prisma(i,3);
    
    k=k+1;
end
X

A=inv(X)

N=P*A
L=formula(N);
for i=[1:6]
N1=subs(L(1,i))
figure
hold on
C1 = [double(subs(N1,[x,y,z],[t1(1),t1(2),t1(3)])) 0 ;
     double(subs(N1,[x,y,z],[t2(1),t2(2),t2(3)])) 0;
     double(subs(N1,[x,y,z],[t6(1),t6(2),t6(3)])) 0;
     double(subs(N1,[x,y,z],[t4(1),t4(2),t4(3)])) 0];
C2 = [double(subs(N1,[x,y,z],[t1(1),t1(2),t1(3)])) 0 ;
     double(subs(N1,[x,y,z],[t3(1),t3(2),t3(3)])) 0;
     double(subs(N1,[x,y,z],[t5(1),t5(2),t5(3)])) 0;
     double(subs(N1,[x,y,z],[t4(1),t4(2),t4(3)])) 0];
C3 = [double(subs(N1,[x,y,z],[t2(1),t2(2),t2(3)])) 0 ;
     double(subs(N1,[x,y,z],[t3(1),t3(2),t3(3)])) 0;
     double(subs(N1,[x,y,z],[t5(1),t5(2),t5(3)])) 0;
     double(subs(N1,[x,y,z],[t6(1),t6(2),t6(3)])) 0];
C4 = [double(subs(N1,[x,y,z],[t1(1),t1(2),t1(3)])) 0 ;
     double(subs(N1,[x,y,z],[t3(1),t3(2),t3(3)])) 0;
     double(subs(N1,[x,y,z],[t2(1),t2(2),t2(3)])) 0;];
C5 = [double(subs(N1,[x,y,z],[t4(1),t4(2),t4(3)])) 0 ;
     double(subs(N1,[x,y,z],[t5(1),t5(2),t5(3)])) 0;
     double(subs(N1,[x,y,z],[t6(1),t6(2),t6(3)])) 0;];
fill3([t1(1), 1; t2(1),1;t6(1),1;t4(1),1],[t1(2),1;t2(2),1;t6(2),1;t4(2),1],[t1(3),1;t2(3),1;t6(3),1;t4(3),1;],C1)

fill3([t1(1),1;t3(1),1;t5(1),1;t4(1),1],[t1(2),1;t3(2),1;t5(2),1;t4(2),1],[t1(3),1;t3(3),1;t5(3),1;t4(3),1], C2)

fill3([t2(1),1;t3(1),1;t5(1),1;t6(1),1],[t2(2),1;t3(2),1;t5(2),1;t6(2),1],[t2(3),1;t3(3),1;t5(3),1;t6(3),1], C3)

fill3([t1(1),1;t3(1),1;t2(1),1],[t1(2),1;t3(2),1;t2(2),1],[t1(3),1;t3(3),1;t2(3),1],C4)

fill3([t4(1),1;t5(1),1;t6(1),1],[t4(2),1;t5(2),1;t6(2),1],[t4(3),1;t5(3),1;t6(3),1], C5)

% fill3([t1(1),t2(1),t6(1),t4(1); 1 1 1 1],[t1(2),t2(2),t6(2),t4(2); 1 1 1 1],[t1(3),t2(3),t6(3),t4(3); 1 1 1 1],C')
%orig
% fill3([t1(1),t2(1),t6(1),t4(1)],[t1(2),t2(2),t6(2),t4(2)],[t1(3),t2(3),t6(3),t4(3)],C)
% fill3([t1(1),t3(1),t5(1),t4(1)],[t1(2),t3(2),t5(2),t4(2)],[t1(3),t3(3),t5(3),t4(3)], 'r')
% fill3([t2(1),t3(1),t5(1),t6(1)],[t2(2),t3(2),t5(2),t6(2)],[t2(3),t3(3),t5(3),t6(3)], 'r')
% fill3([t1(1),t3(1),t2(1)],[t1(2),t3(2),t2(2)],[t1(3),t3(3),t2(3)], C)
% fill3([t4(1),t5(1),t6(1)],[t4(2),t5(2),t6(2)],[t4(3),t5(3),t6(3)], 'r')
%  colormap(cmap)
 
% figure
% fill3([x1,x2,x3],[y1,y2,y3],[double(subs(N2,[x,y],[x1,y1])),double(subs(N2,[x,y],[x2,y2])),double(subs(N2,[x,y],[x3,y3]))], 'b')
% figure
% fill3([x1,x2,x3],[y1,y2,y3],[double(subs(N3,[x,y],[x1,y1])),double(subs(N3,[x,y],[x2,y2])),double(subs(N3,[x,y],[x3,y3]))], 'g')
% 
end







