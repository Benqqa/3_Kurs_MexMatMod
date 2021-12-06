from scipy import integrate
import numpy as np
# на будущее
import math
# ...
import matplotlib.pyplot as plt


# константы интегрирования
a = 0
b = 0.4
step = 0.01
# начальные условия
alpha=1
beta=0
#######
k1=0.015
k2=1720
k=2e11
nu=0.001002
M=0.0048
d=0.0001
S=0.077*0.0195
ro=1000
m=S*ro*d
#######
fi=(k1+k2*nu)*(M+m)
print(fi)
fi=20.68
omega=k/(M+m)
omega=27.707*2

def se_solve(Y,r):
    """
    Функция, задающая систему дифференциальных уравнений
    """
    return [Y[1], -fi*2*Y[1]-(omega**2)*Y[0]]
def main():
    # массивы для решения задачи
    a_x=np.arange(a,b,step)
    y=np.arange(a,b,step)
    # интегрируем систему уравнений
    asol=integrate.odeint(se_solve,[alpha,beta],a_x)
    print(asol)
    # на графике будем изображать квадрат модуля волновой функции
    # несущий физический смысл плотности вероятности
    for i in range(0,len(asol)):
        y[i]=asol[i][0]
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(a_x,y)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.set_title("Численное решение")

    plt.show()
    #  2d
    # theta = m *y
    # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # ax.plot(theta, y)
    # ax.set_rmax(2)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    # ax.grid(True)
    # ax.set_title("A line plot on a polar axis", va='bottom')
    # plt.show()
    
    # #3d
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    
    # # Create the mesh in polar coordinates and compute corresponding Z.
    # r = np.linspace(a,b,len(asol))
    # p = np.linspace(0,2*np.pi, len(asol))
    # R, P = np.meshgrid(r, p)
    # Z = y*np.cos(m*P)
    
    # # Express the mesh in the cartesian system.
    # X, Y = R*np.cos(P), R*np.sin(P)
    
    # # Plot the surface.
    # ax.plot_surface(X, Y, Z, cmap=plt.cm.winter)
    
    # # Tweak the limits and add latex math labels.
    # ax.set_zlim(0, 100)
    
    # ax.set_xlabel("x")
    # ax.set_ylabel("y")
    # ax.set_zlabel("w(r,phi)")
    
    # plt.show()

if __name__ == '__main__':
    main()
    x = np.arange(0, 0.4, 0.01)
    # fi=2
    fig = plt.figure()
    ax = plt.axes()
    ax.set_xlabel("t")
    ax.set_ylabel("x")
   
    
    plt.plot(x,np.sqrt(alpha**2+((beta+fi*alpha)/(omega))**2)*np.exp(-fi*x)*np.cos(np.sqrt(omega**2-fi**2)*x-np.arctan(beta/(omega*alpha)+fi/omega)))
    plt.plot(x,-np.sqrt(alpha**2+((beta+fi*alpha)/(omega))**2)*np.exp(-fi*x),linestyle = '--',
        linewidth = 1,
        color = 'crimson')
    plt.plot(x,np.sqrt(alpha**2+((beta+fi*alpha)/(omega))**2)*np.exp(-fi*x),linestyle = '--',
        linewidth = 1,
        color = 'crimson')
    plt.plot(x,x*0,linestyle = '-',linewidth = 0.5,color='black')
    plt.show()
    
import sympy as sym
sym.init_printing()
  # Определить символические константы x и f (x)
x = sym.Symbol('x', real=True)
f = sym.symbols('f', cls=sym.Function, real=True)
  # Используйте команду diffeq для представления дифференциального уравнения: f '' (x) + 3f '(x) + 2f (x) = exp (-x)
diffeq = sym.Eq(f(x).diff(x, x) +fi* f(x).diff(x) + omega*f(x), 0)
  # Вызов функции dsolve, возврат объекта Eq, точность управления подсказкой
s=sym.dsolve(diffeq, f(x))

print(s)

# sym.plot(s,show=False)