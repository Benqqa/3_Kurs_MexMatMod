Mass_node=[21.,          -6.;24.,           0.;18.,           0.;15.,          -6.;12.,           0.;18.,          -6.;12.,          -6.;9.,          -6.;6.,           0.;3.,          -6.;0.,           0.;6.,          -6. ;   ]
Mass_Element=[1, 2;3, 2;4, 3;5, 4;5, 3;3, 1;6, 1;4, 6;7, 4;8, 7;8, 5;9, 5;10,  9;11, 10;11,  9;9, 8;12,  8;10, 12;]
GU_el=[1,10] %закреп
N_no=length(Mass_node);
N_el=length(Mass_Element);
E=2e11; %юнг
A=0.0001; % площадь
l_e=getLenEls(Mass_Element,Mass_node) %длина каждого элемента

F=zeros(2*length(Mass_Element),1);
F(2*2+1)=1000 ; F(3*2+1)= 1000; F(5*2+1)=1000 ; F(9*2+1) = 1000; F(11*2+1) = 1000;
F1_gl=F
% getK(Mass_node,Mass_Element,l_e)


% [guF_gl,guK_gl]=setGuToFK(F1_gl,K_gl , GU_el)
% function K=getK(Mass_node,Mass_Element,l_e)
%    for i=1:length(Mass_Element)
%     C=[Mass_node(Mass_Element(i,1),:); Mass_node(Mass_Element(i,2),:)]
%     IC=inv(C)
%     B=[1,-1;-1,1]./l_e(i)
%     D=
%     
% %     K=B'*D*B*E*A
%    end
% end
function mass_l=getLenEls(Mass_Element,Mass_node)
    mass_l=zeros(length(Mass_Element),1);
    for i=1:length(Mass_Element)
        mass_l(i)=3
        %sqrt((Mass_node(Mass_Element(i,2),1)-Mass_node(Mass_Element(i,1),1))^2 +(Mass_node(Mass_Element(i,2),2)-Mass_node(Mass_Element(i,1),2))^2);
    end
end
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