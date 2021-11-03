# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 23:18:48 2021

@author: maxbo
"""
import math
import sympy as sp

class Solution():
    def __init__(self,components_table,mix,p,T):
        self.components_table=components_table
        self.mix=mix
        self.p=p #давление
        self.T=T #температура
    def F_V(self):
        V=sp.symbols('V')
        _sum=0
        for component in self.components_table:
           # _sum=
            _sum=_sum+((component.z*(component.K(self.T,self.p)-1))/(V*(component.K(self.T,self.p)-1)+1))
        res=list(sp.solveset(sp.Eq((_sum),0),V,sp.Reals))
        print(res[0])
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
            _sum1=_sum1+component.z*component.K(self.T,self.p)
        _sum2=0
        for component in self.components_table:
            _sum2=_sum2+component.z/component.K(self.T,self.p)
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
    
        

class Component():
    def __init__(self,T_c,p_c,w,R,z):
        self.T_c=T_c #критические температуру 
        self.p_c=p_c #давление
        self.w=w #ацентрический фактор 
        self.R=R # а что это?
        self.z=z #мольная доля компонента в смеси
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
        res=[]
        for value in V:
            res.append((self.z*self.K(T, p))/(value*(self.K(T, p)-1)+1))
        return res
   #
    
class Mix():
     def __init__(self,mix_table):
         self.mix_table=mix_table

Test=Solution([Component(0.0001,1,0.0001,0.0001,0.05),Component(1,0.0001,1,2,0.03),Component(1,0.0001,1,2,0.03),Component(1,0.0001,1,2,0.43),Component(1,1,0.0001,1,0.15)],1,5,10)
TEst_C=Component(0.0001,1,0.0001,0.0001,0.05)
print("rfif",TEst_C.y(1,1,Test.F_V()))
#Test.F_V()
print(Test.F_v(Test.F_V()))