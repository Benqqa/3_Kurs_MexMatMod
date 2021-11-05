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
        self.mix=mix
        self.p=p #давление
        self.T=T #температура
        self.V=self.F_V()
        self.TableOfPairCollisionRate=self.getTableOfPairCollisionRate()
        self.z=self.SupercompressibilityCoefficient()
    def getTableOfPairCollisionRate(self):#коэф парного столкновения
        table = pd.read_excel('./коэф парного взаимодействия.xlsx',index_col='c(ij)')
        return table
       
    def F_V(self):
        V=sp.symbols('V')
        _sum=0
        for component in self.components_table:
           # _sum=
            _sum=_sum+((component.z*(component.K(self.T,self.p)-1))/(V*(component.K(self.T,self.p)-1)+1))
        res=list(sp.solveset(sp.Eq((_sum),0),V,sp.Reals))
        print("f_v",res)
        return res
    def F_v(self,V):
        for value in V:
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
            _sum1+=component.z*component.K(self.T,self.p)
        _sum2=0
        for component in self.components_table:
            _sum2+=+component.z/component.K(self.T,self.p)
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
   # поиск -Коэффициенту сверхсжимаемости
    
    #что-то вспомогательное
    def A_m(self):
        res=[]
        for value in self.a_m():  
            res.append(value*self.p/((self.R**2)*(self.T)**2))
        return res
    def B_m(self):
        res=[]
        for value in self.b_m():  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    def C_m(self):
        res=[]
        for value in self.c_m():  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    def D_m(self):
        res=[]
        for value in self.d_m():  
            res.append(value*self.p/((self.R)*(self.T)))
        return res
    #конец чего-то вспомогательного
    def SupercompressibilityCoefficient(self):
        mass=[]
        A1=self.A_m()
        B1=self.B_m()
        C1=self.C_m()
        D1=self.D_m()
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
    def ln_f(self,component):
        A=self.A_m()
        B=self.B_m()
        C=self.C_m()
        D=self.D_m()
        res=[]
        for value in self.V:
            _sum=0
            for component_j in self.components_table:
                _sum+=component_j.y(self.T,self.p,value)*(1-self.C(component.name,component_j.name))*math.sqrt(component.a(self.T)*component_j.a(self.T))
            res.append( math.log1p(component.y(self.T,self.p,value)*self.p) - math.log1p(self.z-B)-(A)*(((2*_sum)/(self.a_m()))-((component.c()-component.d())/(self.c_m()-self.d_m())))*math.log1p((self.z+C)/(self.z+D))/(C-D) +B/(self.z-B)-A*((C/(self.z+C))-(D/(self.z+D)))/(C-D) )# кто такой Bi
        return res
   #
class Component():
    def __init__(self,name,T_c,p_c,w,R,z):
        self.T_c=T_c #критические температуру 
        self.p_c=p_c #давление
        self.w=w #ацентрический фактор 
        self.R=8.314 # а что это?
        self.z=z #мольная доля компонента в смеси
        self.name=name
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
    def K(self,T,p):
        return self.p_Si(T)/p
   #
   #уравнением фазовой концентрации компонентов смеси.
    def y(self,T,p,V):
        return (self.z*self.K(T, p))/(V*(self.K(T, p)-1)+1)
        


   #
    
class Mix():
     def __init__(self,mix_table):
         self.mix_table=mix_table

Test=Solution([Component("CH4",0.0001,1,0.0001,0.0001,0.05),Component("CH4",1,0.0001,1,2,0.03),Component("CH4",1,0.0001,1,2,0.03),Component("CH4",1,0.0001,1,2,0.43),Component("CH4",1,1,0.0001,1,0.15)],1,5,10)
TEst_C=Component("CH4",0.0001,1,0.0001,0.0001,0.05)
print("rfif",TEst_C.y(1,1,Test.V[1]))
#Test.F_V()
print(Test.F_v(Test.V))