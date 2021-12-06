from scipy import integrate
import numpy as np
# на будущее
import math
# ...
import matplotlib.pyplot as plt


# константы интегрирования
a = 0.0001
b = 3
step = 0.01
# начальные условия
alpha=1
beta=0
# неиспользуемые дополнительные переменные
F=100
S=0.002827
ro=7850
alfa=F/(S*ro)
m=2

def se_solve(Y,r):
    """
    Функция, задающая систему дифференциальных уравнений
    """
    return [Y[1], -1/r*Y[1]-(alfa*alfa*r*r-m*m)*Y[0]/(r*r)]
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
    plt.plot(a_x,y)
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
    
    #3d
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    # Create the mesh in polar coordinates and compute corresponding Z.
    r = np.linspace(a,b,len(asol))
    p = np.linspace(0,2*np.pi, len(asol))
    R, P = np.meshgrid(r, p)
    Z = y*np.cos(m*P)
    
    # Express the mesh in the cartesian system.
    X, Y = R*np.cos(P), R*np.sin(P)
    
    # Plot the surface.
    ax.plot_surface(X, Y, Z, cmap=plt.cm.winter)
    
    # Tweak the limits and add latex math labels.
    ax.set_zlim(0, 10**8)
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("w(r,phi)")
    
    plt.show()

if __name__ == '__main__':
    main()
    
# import sympy as sym
# sym.init_printing()
#   # Определить символические константы x и f (x)
# x = sym.Symbol('x', real=True)
# f = sym.symbols('f', cls=sym.Function, real=True)
#   # Используйте команду diffeq для представления дифференциального уравнения: f '' (x) + 3f '(x) + 2f (x) = exp (-x)
# diffeq = sym.Eq(x*x*f(x).diff(x, x) + x*f(x).diff(x) + (x*x-100)*f(x), 0)
#   # Вызов функции dsolve, возврат объекта Eq, точность управления подсказкой
# s=sym.dsolve(diffeq, f(x))

# print(s)

# sym.plot(s,show=False)