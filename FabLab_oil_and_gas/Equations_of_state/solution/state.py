# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 23:18:48 2021

@author: maxbo
"""
import math
import symPy

class Solution():
    def __init__(self,components_table,mix,p,T):
        self.components_table=components_table
        self.mix=mix
        self.p=p #давление
        self.T=T #температура
    def F_V(self):
        V=symPy.symbols('V')
        _sum=0
        for component in self.components_table:
            _sum=
           # _sum=_sum+((component.z*(component.K(self.T,self.p)-1))/(V*(component.K(self.T,self.p)-1)+1))
        

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
    
class Mix():
     def __init__(self,mix_table):
         self.mix_table=mix_table

        