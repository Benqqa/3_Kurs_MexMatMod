Mass_node=[ 21,-6;
            24,0;
            18,0;
            15,-6;
            12,0;
            18,-6;
            12,-6;
            9,-6;
            6,0;
            3,-6;
            0,0;
            6,-6 ]
Mass_Element=[
    1, 2;
    3, 2;
    4, 3;
    5, 4;
    5, 3;
    3, 1;
    6, 1;
    4, 6;
    7, 4;
    8, 7;
    8, 5;
    9, 5;
    10,  9;
    11, 10;
    11,  9;
    9, 8;
    12,  8;
    10, 12;
    ]

GU_el=[5] %закреп
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
F(11*2)= 1000;
F(9*2)= 1000;
F(5*2)= 1000;
F(3*2)= 1000; F(2*2)= 1000; 
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
%gtht
U=guK_gl\guF_gl


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
% Mass_el=[];
% for i=1:length(Mass_Element)
%     for j=1:length(Mass_Element(1,:))
%         Mass_el=[Mass_el;Mass_Element(i,j)];
%     end
% end
% Mass_el

    for i =1:length(Mass_Element)
        nodes=Mass_Element(i,:)
        T=getT([Mass_node(nodes(1),1),Mass_node(nodes(2),1)],[Mass_node(nodes(1),2),Mass_node(nodes(2),2)], l_e(i));
        k_loc =[1,-1;-1,1].*massK_loc(i)
        T
        k_gl=T'*k_loc*T
        for k=1:2
            for j=1:2
                k
                
%                 K_gl(2*nodes(k)-1,2*nodes(k)-1)=K_gl(2*nodes(k)-1,2*nodes(k)-1)+k_gl(1,1);
%                 K_gl(2*nodes(k),2*nodes(k))=K_gl(2*nodes(k),2*nodes(k))+k_gl(2,2);
%                 K_gl(2*nodes(j)-1,2*nodes(j)-1)=K_gl(2*nodes(j)-1,2*nodes(j)-1)+k_gl(3,3);
%                 K_gl(2*nodes(j),2*nodes(j))=K_gl(2*nodes(j),2*nodes(j))+k_gl(4,4);
                K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+0)     =K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+0)                  +k_gl(2*(k-1)+1+0,2*(j-1)+1+0);
                
                K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+1)     =K_gl(2*(nodes(k)-1)+1+0,2*(nodes(j)-1)+1+1)                  +k_gl(2*(k-1)+1+0,2*(j-1)+1+1);               
                K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+0)     =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+0)                  +k_gl(2*(k-1)+1+1,2*(j-1)+1+0);
                K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1)     =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1)                  +k_gl(2*(k-1)+1+1,2*(j-1)+1+1);
                if(2*(nodes(k)-1)+1+1 == 24)
                    wwdawdwda1=k_gl(2*(k-1)+1+1,2*(j-1)+1+1)
                    k_gl
                    2*(k-1)+1+1
                    2*(j-1)+1+1
                    k_loc
                    T
                end
                if(2*(nodes(k)-1)+1+1 == 24)
                    wwdawdwda2=k_gl(2*(k-1)+1+1,2*(j-1)+1+1)
                    k_gl
                    2*(k-1)+1+1
                    2*(j-1)+1+1
                    k_loc
                    T
                end
%                  K_gl(2*(nodes(k)-1)+1,2*(nodes(j)-1)+1)     =K_gl(2*(nodes(k)-1)+1,2*(nodes(j)-1)+1)                  +i;
%                 K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1) =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1+1)              +i;
%                 K_gl(2*(nodes(k)-1)+1,2*(nodes(j)-1)+1+1)   =K_gl(2*(nodes(k)-1)+1,2*(nodes(j)-1)+1+1)                +i;
%                 K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1)   =K_gl(2*(nodes(k)-1)+1+1,2*(nodes(j)-1)+1)                +i;
            end
        end
    end
    K_gl
end
