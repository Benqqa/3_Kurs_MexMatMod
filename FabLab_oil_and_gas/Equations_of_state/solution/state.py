# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 23:18:48 2021

@author: maxbo
"""
import math
import sympy as sp
import pandas as pd
class Solution():
    def __init__(self,components_table,mix,p,T):
        self.components_table=components_table
        self.R=8.314
        self.eps=10**(-5)
        self.mix=mix
        self.p=p #давление
        self.T=T #температура
        for component in self.components_table:
            component.f_K(self.T, self.p)
        self.V=self.F_V()
        self.TableOfPairCollisionRate=self.getTableOfPairCollisionRate()
        self.z_y=self.SupercompressibilityCoefficient(self.a_m(),self.b_m(),self.c_m(),self.d_m())
        self.z_x=self.SupercompressibilityCoefficient(self.x_a_m(),self.x_b_m(),self.x_c_m(),self.x_d_m())
        
    def getTableOfPairCollisionRate(self):#коэф парного столкновения
        table = pd.read_excel('./коэф парного взаимодействия.xlsx',index_col='c(ij)')
        return table
       
    def F_V(self):
        
        V=sp.symbols('V')
        _sum=0
        for component in self.components_table:
            _sum=_sum+((component.z*(component.K-1))/(V*(component.K-1)+1))
        res=list(sp.solveset(sp.Eq((_sum),0),V,sp.Reals))
        print("f_v",res)
        return res
    
    def F_v(self,):
        for value in self.V:
            print("value== ",value)
            if(value<0):
                print("V ненасыщенном жидком состоянии")
            elif(value==0):
                print("V насыщенное жидкое состояние (точка кипения).")
            elif(value>0 and value<1):
                print("V двухфазное парожидкостное состояние.")
            elif(value==1):
                print("V однофазное насыщенное паровое (газовое) состояние (точка росы).")
            elif(value>1):
                print("V однофазное ненасыщенное газовое состояние.")
        _sum1=0
        for component in self.components_table:
            _sum1+=component.z*component.K
        _sum2=0
        for component in self.components_table:
            _sum2+=+component.z/component.K
        print(_sum1)
        print(_sum2)
        if(_sum1<1):
            return "ненасыщенном жидком состоянии"
        elif(_sum1==1):
            return "насыщенное жидкое состояние (точка кипения)." 
        elif(_sum2>1 and _sum1>1 ):
            return "двухфазное парожидкостное состояние." 
        elif(_sum2==1):
            return "однофазное насыщенное паровое (газовое) состояние (точка росы)." 
        elif(_sum2<1):
            return "однофазное ненасыщенное газовое состояние."
    def C(self, row_name,column_name): # коэффициентов парного взаимодействия
        table = self.TableOfPairCollisionRate
        i=table.columns.values.tolist().index(row_name)
        return table[column_name].iloc[i]
    #какие то коэфициенты
    def a_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                for component_j in self.components_table:
                    _sum+=component_i.y(self.T, self.p, value)*component_j.y(self.T, self.p, value)*(1-self.C(component_i.name,component_j.name))*math.sqrt(component_i.a(self.T)*component_j.a(self.T))
            res.append(_sum)
        return res
    def b_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.y(self.T, self.p, value)*component_i.b()
            res.append(_sum)
        return res
    def c_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.y(self.T, self.p, value)*component_i.c()
            res.append(_sum)
        return res
    def d_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.y(self.T, self.p, value)*component_i.d()
            res.append(_sum)
        return res
    #
    #какие то коэфициенты2
    def x_a_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                for component_j in self.components_table:
                    _sum+=component_i.x(self.T, self.p, value)*component_j.y(self.T, self.p, value)*(1-self.C(component_i.name,component_j.name))*math.sqrt(component_i.a(self.T)*component_j.a(self.T))
            res.append(_sum)
        return res
    def x_b_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.x(self.T, self.p, value)*component_i.b()
            res.append(_sum)
        return res
    def x_c_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.x(self.T, self.p, value)*component_i.c()
            res.append(_sum)
        return res
    def x_d_m(self):
        res=[]
        for value in self.V:    
            _sum=0
            for component_i in self.components_table:
                    _sum+=component_i.y(self.T, self.p, value)*component_i.d()
            res.append(_sum)
        return res
    #
   # поиск -Коэффициенту сверхсжимаемости
    
    #что-то вспомогательное
    def A_m(self,a_m):
        res=[]
        for value in a_m:  
            res.append(value*self.p/((self.R**2)*(self.T)**2))
        return res
    def B_m(self, b_m):
        res=[]
        for value in b_m:  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    def C_m(self,c_m):
        res=[]
        for value in c_m:  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    def D_m(self,d_m):
        res=[]
        for value in d_m:  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    #конец чего-то вспомогательного
    def SupercompressibilityCoefficient(self,a_m,b_m,c_m,d_m):
        mass=[]
        A1=self.A_m(a_m)
        B1=self.B_m(b_m)
        C1=self.C_m(c_m)
        D1=self.D_m(d_m)
        for i in range(len(A1)):  
            z=sp.symbols('z')
            A=A1[i]
            B=B1[i]
            C=C1[i]
            D=D1[i]
            _sum=z**3+(C+D-B-1)*(z**2)+(A-B*C+C*D-B*D-D-C)*z-(B*C*D+C*D+A*B)
            Equ=sp.Eq((_sum),0)
            res=list(sp.solveset(Equ,z,sp.Reals))
            mass.append(max(res))
        print("max",max(mass))
        return max(mass)
   #    
   #логарифм летучести компонентов в паровой фазе
    def ln_f(self,component,a_m,b_m,c_m,d_m,z):
        A1=self.A_m(a_m)
        B1=self.B_m(b_m)
        C1=self.C_m(c_m)
        D1=self.D_m(d_m)
        a1=a_m
        b1=b_m
        c1=c_m
        d1=d_m
        res=[]
        for i in range(len(A1)):  
            value=self.V[i]
            A=A1[i]
            B=B1[i]
            C=C1[i]
            D=D1[i]
            a=a1[i]
            b=b1[i]
            c=c1[i]
            d=d1[i]
            _sum=0
            for component_j in self.components_table:
                _sum+=component_j.y(self.T,self.p,value)*(1-self.C(component.name,component_j.name))*math.sqrt(component.a(self.T)*component_j.a(self.T))
            if(component.y(self.T,self.p,value)>=0):
                res.append( math.log1p(component.y(self.T,self.p,value)*self.p) - math.log1p(z-B)-(A)*(((2*_sum)/a))-((component.c()-component.d())/(c-d))*math.log1p((z+C)/(z+D))/(C-D) +B/(z-B)-A*((C/(z+C))-(D/(z+D)))/(C-D) )# кто такой Bi

        return res
   #
   #Корректируют значения коэффициентов распределения
    def K_2(self):
        flag=True
        while(flag):
            for component in self.components_table:
                ln_f_L=max(self.ln_f(component, self.a_m(), self.b_m(), self.c_m(), self.d_m(), self.z_y))
                ln_f_V=max(self.ln_f(component, self.x_a_m(), self.x_b_m(), self.x_c_m(), self.x_d_m(), self.z_x))
                print(ln_f_L)
                print('ln_f_V',ln_f_V)
                f_L=math.exp(ln_f_L)
                f_V=math.exp(ln_f_V)
                print(f_L)
                print(f_V)
                if(f_V != 0):
                    if(math.fabs((f_L/f_V)-1)>self.eps):
                        print(self.F_v())
                        flag =False
                        break
                    
                    component.K=component.K*f_L/f_V
class Component():
    def __init__(self,name,T_c,p_c,w,R,z):
        self.T_c=T_c #критические температуру 
        self.p_c=p_c #давление
        self.w=w #ацентрический фактор 
        self.R=8.314 # а что это?
        self.z=z #мольная доля компонента в смеси
        self.name=name
        self.K=1
        
        #self.T=T # а это откуда?
    def Z_c(self):
        return 0.3357-0.0294*self.w
    def Omega_c(self):
        return 0.75001
    def Psy(self):
        if(self.w < 0.4489):
            return 1.050+0.105*self.w+0.482*(self.w**2)
        elif(self.w>0.4489):
            return 0.429+1.004*self.w+1.56*(self.w**2)
        else:
            return 1.194
   #coefficients of the equation of state
    def alfa(self):
        return self.Omega_c()**3
    def betta(self):
        return self.Z_c()+self.Omega_c()-1
    def sigma(self):
        return -1*self.Z_c()+self.Omega_c()*(0.5+(self.Omega_c()-0.75)**(0.5))
    def delta(self):
        return -1*self.Z_c()+self.Omega_c()*(0.5-(self.Omega_c()-0.75)**(0.5))
   #
   #values
    def a(self,T):
        a_c=(self.alfa()*((self.R)**2)*((self.T_c)**2))/self.p_c
        calc=(1+self.Psy()*(1-((T/self.T_c)**(0.5))))**2
        return (a_c*calc)
    def b(self):
        return self.betta()*self.R*self.T_c/self.p_c
    def c(self):
        return self.sigma()*self.R*self.T_c/self.p_c
    def d(self):
        return self.delta()*self.R*self.T_c/self.p_c
   #
   #distribution coefficients of components
    def p_Si(self,T):
        return math.exp(5.373*(1+self.w)*(1-self.T_c/T))*self.p_c
    def f_K(self,T,p):
        self.K=self.p_Si(T)/p
        #return self.p_Si(T)/p
   #
   #уравнением фазовой концентрации компонентов смеси.
    def y(self,T,p,V):
        return (self.z*self.K)/(V*(self.K-1)+1)
   
   #
   #мольные доли компонентов смеси в жидкой фазе. 
    def x(self,T,p,V):
        return self.z/(V*(self.K-1)+1)
        
   #
    
class Mix():
     def __init__(self,mixes,consts):
         self.mixes=mixes
         self.consts=consts
     def CreateMix(self,mix_number):
        table_mix = pd.read_excel('./'+self.mixes)
        headers_mix=table_mix.columns.values.tolist()
        values_mix=table_mix.iloc[mix_number]
        mass=[]
        #табличка с означениями температуры и давления + вырезать из класса части R
        for i in range(len(values_mix)):
            mass.append(Component(headers_mix[i],T,p,w,R,values_mix[i]))
        return table
    
mass=[Component("CH4",0.0001,1,0.0001,0.0001,0.05),Component("CH4",1,0.0001,1,2,0.03),Component("CH4",1,0.0001,1,2,0.03),Component("CH4",1,0.0001,1,2,0.43),Component("CH4",1,1,0.0001,1,0.15)]
# Test=Solution(mass,1,5,5)
# TEst_C=Component("CH4",0.0001,1,0.0001,0.0001,0.05)
# Test.K_2()
Mix1=Mix("mix.xlsx")
a=Mix1.CreateMix(1)

# print("rfif",TEst_C.y(1,1,Test.V[1]))
#Test.F_V()
# print(Test.ln_f(TEst_C,Test.a_m(),Test.b_m(),Test.c_m(),Test.d_m(),Test.z_y))
# print(Test.ln_f(TEst_C,Test.x_a_m(),Test.x_b_m(),Test.x_c_m(),Test.x_d_m(),Test.z_x))