clear clc;
x1 = 0;y1 = 0;z1 = 0; %точка 1
x2 = 8;y2 = 8;z2 = 6;  %точка 2
x3 = x1;y3 = y1;z3=z2; 
figure
hold on
aaaaa(0,1,1,0,0,0,0,0,1)
aaaaa(0,0,1,0,0,1,1,0,0)
aaaaa(0,0,0,0,1,1,1,0,1)

function []=aaaaa(x1,x2,x3,y1,y2,y3,z1,z2,z3)

x_a = [x1 x2 x3];
y_a = [x1 x2 x3];
z_a = [z1 z2 z3];
A1 = [x2-x1 y2-y1 z2-z1];
A2 = [x3-x2 y3-y2 z3-z2];
C=cross(A1,A2);
l=10;
h_x=(max(y_a)-min(y_a))/l;
h_y=(max(y_a)-min(y_a))/l;
[X,Y]=meshgrid(min(x_a):h_x:max(x_a), min(y_a):h_y:max(y_a));
if C(3)==0
    Y=-C(1)*X/C(2);
    hz=(max(z_a)-min(z_a))/(length(X)-1);
    Z1=min(z_a):hz:max(z_a);
    for j=1:length(Y)
        Z(:,j)=Z1;
    end
 else
    Z=(C(1)*X+C(2)*Y)/C(3);
end

surf(X,Y,Z);
xlabel('x');ylabel('y');zlabel('z')
end