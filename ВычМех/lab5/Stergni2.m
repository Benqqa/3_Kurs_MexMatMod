Mass_node=[6.78988647,           0.;
3.75697184,  0.741771281;
6.78988647,          1.5;
9.38103199,   2.14778638;
12.7898865,           0.;
12.7898865,           3.;
15.0334682,   3.56089544;
18.7898865,           0.;
18.7898865,          4.5;
21.7898865,         5.25;
24.7898865,           0.;
24.7898865,           6.;
0.789886594,           0.;
 ]
Mass_Element=[
1, 2;
3, 2;
1, 3;
4, 1;
5, 4;
5, 6;
7, 6;
7, 5;
5, 8;
8, 7;
9, 7;
8, 9;
10,  8;
11, 10;
11, 12;
12, 10;
10,  9;
2, 13;
4, 3;
1, 5;
13,  1;
6, 4;
8, 11;
    ]
GU_el=[13,11] %закреп
N_no=length(Mass_node);
N_el=length(Mass_Element);
E=2e11; %юнг
A=0.0001; % площадь
l_e=getLenEls(Mass_Element,Mass_node) %длина каждого элемента
massK_loc=getK_loc(E,A,l_e) %*[1,-1;-1,1] не забвть
%сосредоточеная сила
% F=zeros(length(Mass_node),2);
% F(2,2)=1000 ; F(3,2)= 1000; F(5,2)=1000 ; F(9,2) = 1000; F(11,2) = 1000;
% F1_gl = getF_gl(F,Mass_Element,Mass_node, l_e )
F=zeros(2*length(Mass_node),1);
% F(2*2)=1000 ;
% F(3*2)= -1000; F(5*2)=-1000 ; F(9*2) = -1000; %F(11*2) = 1000;
F(2*2)= -1000; F(3*2)=-1000 ; F(4*2)=-1000 ; F(6*2)=-1000 ; F(7*2)=-1000 ; F(9*2)=-1000 ; F(10*2)=-1000 ; F(12*2)=-1000 ;
F1_gl=F

%
%матрица жесткости
K_gl=getK_gl(Mass_Element,massK_loc,Mass_node, l_e )
%вектор сил
% F1_gl= getF_gl(F,Mass_Element,Mass_node, l_e )
length(F1_gl)
length(K_gl(1,:))
%задние ГУ
[guF_gl,guK_gl]=setGuToFK(F1_gl,K_gl , GU_el)
det(guK_gl)
% guF_gl=guF_gl(1:length(guF_gl)-1)
% guK_gl=guK_gl(1:length(guK_gl)-1,:)
% guF_gl=[guF_gl(1:11),guF_gl(13:length(guK_gl)-1)]
% guK_gl=[guK_gl(1:11,:),guK_gl(13:length(guK_gl)-1,:)]
%перемещение
U=guK_gl\guF_gl
%деформации
new_len=zeros(length(Mass_Element),1);
deffs=zeros(length(Mass_Element),1);
for i=1:length(Mass_Element)
    nodes=Mass_Element(i,:);
    X1=Mass_node(nodes(1),1)+U(2*nodes(1)-1);
    X2=Mass_node(nodes(2),1)+U(2*nodes(2)-1);
    Y1=Mass_node(nodes(1),2)+U(2*nodes(1));
    Y2=Mass_node(nodes(2),2)+U(2*nodes(2));
    new_len(i)=sqrt((X2-X1)^2 +(Y2-Y1)^2);
    deffs(i)=(new_len(i)-l_e(i))/l_e(i);
    
end
stesses=deffs*E
Forces=stesses*A

function [F1_gl,K_gl]=setGuToFK(F1_gl,K_gl , GU_el)
    for i=1:length(F1_gl)
        for j=1:length(GU_el)
            if(i == 2*(GU_el(j)-1)+1)
                F1_gl(i)  = 0;
                
                F1_gl(i+1)  = 0;
                
                K_gl(i,:) = 0;
                K_gl(:,i) = 0;
                
                
                K_gl(i+1,:) = 0;
                K_gl(:,i+1) = 0;
      
                K_gl(i,i) = 1;
                K_gl(i+1,i+1) = 1;
            end
        end
    end
end
function mass_k_loc=getK_loc(E,A,l_e)
    mass_k_loc=zeros(length(l_e),1);
    for i=1: length(l_e)
        mass_k_loc(i)=E*A/l_e(i);
    end
end
function mass_l=getLenEls(Mass_Element,Mass_node)
    mass_l=zeros(length(Mass_Element),1);
    for i=1:length(Mass_Element)
        nodes=Mass_Element(i,:)
        mass_l(i)=sqrt((Mass_node(nodes(2),1)-Mass_node(nodes(1),1))^2 +(Mass_node(nodes(2),2)-Mass_node(nodes(1),2))^2);
    end
    mass_l
end
function T=getT(x_loc,y_loc,l)
    l_loc=(x_loc(2)-x_loc(1))/l;
    m_loc=(y_loc(2)-y_loc(1))/l;
    T=[l_loc, m_loc, 0    0;
        0,       0, l_loc, m_loc]
end
function K_gl = getK_gl(Mass_Element,massK_loc,Mass_node, l_e )
K_gl=zeros(2*length(Mass_node));
    for i =1:length(Mass_Element)
        nodes=Mass_Element(i,:)
        T=getT([Mass_node(nodes(1),1),Mass_node(nodes(2),1)],[Mass_node(nodes(1),2),Mass_node(nodes(2),2)], l_e(i));
        k_loc =[1,-1;-1,1].*massK_loc(i)
        T
        k_gl=T'*k_loc*T
        for k=1:2
            for j=1:2
                K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+0)     =K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+0)                  +k_gl(2*(k-1)+1+0,2*(j-1)+1+0);
                K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+1)     =K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+1)                  +k_gl(2*(k-1)+1+0,2*(j-1)+1+1);               
                K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+0)     =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+0)                  +k_gl(2*(k-1)+1+1,2*(j-1)+1+0);
                K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1)     =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1)                  +k_gl(2*(k-1)+1+1,2*(j-1)+1+1);
            end
        end
    end
    K_gl
end
