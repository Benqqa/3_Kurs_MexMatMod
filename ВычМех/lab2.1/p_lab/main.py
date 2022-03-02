import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.linalg import solve


class Shveller:
    pass


class Balcka:
    def __init__(self, x, l, E, rho, h, w, Mass_node, Mass_Element, P):
        print("__init__ Balcka")
        self.x = x
        self.l = l
        self.E = E
        self.rho = rho
        self.g = 9.81
        self.P = P
        self.Mass_node = Mass_node
        self.Mass_Element = Mass_Element
        self.type = self.Shveller(0.140, 0.058, 0.0049, 0.0049, 2)

        self.S = self.type.S
        # self.S = self.S(h, w)
        self.k_e = self.CreateMatrix_k_e()
        self.getK_global()
        self.getF_global()
        # self.Solve()

    class Shveller:
        def __init__(self, h, b, s, t, variant):
            self.i_y = None
            self.i_x = None
            self.I_y = None
            self.I_x = None
            self.W_y = None
            self.W_X = None
            self.S = None
            print("__init__ Shveller")
            self.h = h
            self.b = b
            self.s = s
            self.t = t
            self.variant = variant
            self.Calc()

        def Calc(self):
            if self.variant == 1:
                self.S = 0.001206  # Площадь швеллера
                self.W_X = 0.000049943  # Момент сопротивления Wx
                self.W_y = 0.000008655  # Момент сопротивления Wy
                self.I_x = 0.000003496003  # Момент инерции Ix
                self.I_y = 0.000000372506  # Момент инерции Iy
                self.i_x = 0.054  # Радиус инерции ix
                self.i_y = 0.018  # Радиус инерции iy
            elif self.variant == 2:
                self.I_x = 0.000003496003101266667  # Момент инерции Ix
                # self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (self.b - 0) ** 3) / 3
                self.I_y = 2.422756617262131 * (10 ** (-7))  # Момент инерции Iy
                self.I_x = (self.b * self.h ** 3 - (self.b - self.t) * (self.h - 2 * self.t) ** 3) / 12
                self.S = 0.00092218  # Площадь швеллера
                # self.S = 0.001206  # Площадь швеллера
            else:
                self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (
                            self.b - 0) ** 3) / 3

    def S(self, h, w):
        return h * w

    def J(self):
        return self.type.I_x

    def N_i(self, eta):
        return 1 / 4 * ((1 - eta) ** 2) * (2 + eta)

    def N_i_theta(self, eta):
        return 1 / 8 * self.l*(1 + eta) * ((1 - eta) ** 2)

    def N_j(self, eta):
        return 1 / 4 * (2 - eta) * ((1 + eta) ** 2)

    def N_j_theta(self, eta):
        return -1 / 8 * self.l*(1 - eta) * ((1 + eta) ** 2)

    def CreateMatrix_B(self, eta):
        B = np.zeros((1, 4))
        B[0, 0] = 6 * eta / self.l
        B[0, 1] = 3 * eta - 1
        B[0, 2] = -6 * eta / self.l
        B[0, 3] = 3 * eta + 1
        B = B / self.l
        # print(B)
        return B

    def CreateMatrix_k_e(self):
        k_e = np.zeros((4, 4))
        k_e[0, :] = [12, 6 * self.l, -12, 6 * self.l]
        k_e[1, :] = [6 * self.l, 4 * self.l ** 2, -6 * self.l, 2 * self.l ** 2]
        k_e[2, :] = [-12, -6 * self.l, 12, - 6 * self.l]
        k_e[3, :] = [6 * self.l, 2 * self.l ** 2, -6 * self.l, 4 * self.l ** 2]
        k_e = k_e * self.E * self.J() / (self.l ** 3)
        print("k_e")
        # print(k_e)
        return k_e

    def CreateRow_f_e(self, eta):
        f_e = np.zeros((1, 4))
        # print(f_e)
        N = self.CreateColumn_N(eta)
        # P = self.P.transpose()
        # N.transpose() * P +self.l / 2 * P * np.matrix([[1], [self.l / 6], [1], [-self.l / 6]]).transpose()+
        f_e = N.transpose() * P + self.l / 2 * P * np.matrix(
            [[1], [self.l / 6], [1], [-self.l / 6]]).transpose() + self.l / 2 * self.rho * self.S * self.g * np.matrix(
            [[1], [self.l / 6], [1], [-self.l / 6]]).transpose()
        print(f_e)
        return f_e

    def CreateColumn_N(self, eta):
        N = np.zeros((4, 1))
        N[:, 0] = [self.N_i(eta), self.N_i_theta(eta), self.N_j(eta), self.N_j_theta(eta)]
        # print(N)
        return N

    def Eta(self, x):
        return 2 * x / self.l - 1

    def getK_global(self):
        K = np.zeros((2 * len(self.Mass_node), 2 * len(self.Mass_node)))
        for i in range(len(self.Mass_Element)):
            nodes = self.Mass_Element[i, :]
            for k in range(2):
                for j in range(2):
                    K[2 * nodes[0, k], 2 * nodes[0, j]] += self.k_e[2 * k, 2 * j]
                    K[(2 * nodes[0, k]) + 1, (2 * nodes[0, j]) + 1] += self.k_e[(2 * k) + 1, (2 * j) + 1]
                    K[2 * nodes[0, k], (2 * nodes[0, j]) + 1] += self.k_e[2 * k, (2 * j) + 1]
                    K[(2 * nodes[0, k]) + 1, 2 * nodes[0, j]] += self.k_e[(2 * k) + 1, 2 * j]
        # print(K)
        # print('\n'.join(' '.join(str(col) for col in row) for row in K))
        return K

    def getF_global(self):
        X = [0, 1 / 10, 2 / 10, 3 / 10, 4 / 10, 5 / 10, 6 / 10, 7 / 10, 8 / 10, 9 / 10, 10 / 10]
        F = 0
        # for i in range(len(Mass_Element)):
        #     x = X[i]
        #     F += self.CreateRow_f_e(self.Eta(x))
        F2 = self.P
        # for i in range(len(self.M))
        F = self.P  # +self.l / 2 * self.rho * self.S * self.g * np.matrix([[1], [self.l / 6], [1], [-self.l / 6]]).transpose()
        return F

    def setGuToK(self, GU):
        K_new = self.getK_global()
        for i in range(len(K_new)):
            for j in range(len(GU)):
                if (i == 2 * (GU[j] - 1)):

                    K_new[i, :] = 0
                    K_new[:, i] = 0
                    K_new[i + 1, :] = 0
                    K_new[:, i + 1] = 0
                    K_new[i, i] = 1
                    K_new[i + 1, i + 1] = 1

        return K_new

    def Solve(self, GU):

        K_new = self.getK_global()
        K_gu = self.setGuToK(GU)
        print("Fssssssss")
        # print(self.getF_global())
        # print(K_new)
        print(self.getF_global())
        F = self.getF_global().transpose()
        print("pjzw")
        print('\n'.join(' '.join(str(col) for col in row) for row in K_gu))
        print("det", np.linalg.det(K_gu))
        # print(K_gu)
        U = solve(K_gu, F)
        # print(U)
        [X,Y]=self.getDeff(U)
        [X1, U_new] = self.getUfromN(U,0)


        
        stresses=X*self.E
        forces=stresses*self.type.S
        plt.subplot(2,2,1)
        self.Graph(X1, U_new)
        plt.subplot(2,2,2)
        plt.plot(np.linspace(0, self.x, 11), X)
        plt.scatter(np.linspace(0, self.x, 11), X, color='orange', s=40, marker='o')
        plt.subplot(2,2,3)
        plt.plot(np.linspace(0, self.x, 11), stresses)
        plt.scatter(np.linspace(0, self.x, 11), stresses, color='orange', s=40, marker='o')
        plt.subplot(2, 2, 4)
        plt.plot(np.linspace(0, self.x, 11), forces)
        plt.scatter(np.linspace(0, self.x, 11), forces, color='orange', s=40, marker='o')
        plt.show()

        return U
    def getDeff(self,U):
        X = np.zeros((int(len(U) / 2)))
        Y = np.zeros((int(len(U) / 2)))
        k = 0
        j = 0
        for i in range(len(U)):
            # print(U[i])
            if (i % 2) != 0:
                Y[k] = U[i]
                k += 1
            else:
                X[j] = U[i]
                j += 1
        self.X = X
        self.Y=Y
        print(Y, "Y")
        print(X, "X")
        return [X,Y]
    def getUfromN(self, U,number):
        X = np.zeros((int(len(U) / 2 - 1), 20))
        U_new = np.zeros((int(len(U) / 2 - 1), 20))
        for i in range(int(len(U) / 2 - 1)):
            X[i] = np.linspace(i * self.l, (i + 1) * self.l, 20)
            # print(np.linspace(i*self.l, (i+1)*self.l, 20))
        for i in range(int(len(U) / 2 - 1)):
            U_new[i] = self.N_i(self.Eta(X[number])) * U[2 * i] + self.N_i_theta(self.Eta(X[number])) * U[2 * i + 1] + self.N_j(self.Eta(X[number])) * U[2 * i + 2] + self.N_j_theta(self.Eta(X[number])) * U[2 * i + 3]

            print(U_new[i], "U")
        return [X, U_new]

    def Graph(self, X, U):
        for i in range(len(X)):
            plt.plot(X[i], U[i])
            # plt.scatter(X[i], U[i], color='orange', s=40, marker='o')
        plt.scatter(np.linspace(0, self.x, 11), self.X, color='orange', s=40, marker='o')



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    E = 2 * (10 ** 11)  # сталь
    rho = 7700
    M = 10000
    l = 0.1
    x = 1
    w = 0.05
    h = 0.1
    Mass_node = np.matrix(
        [[0, 0], [0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 0], [0.6, 0], [0.7, 0], [0.8, 0],
         [0.9, 0], [1, 0]])
    Mass_Element = np.matrix(
        [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10]])
    P = np.zeros((1, len(Mass_Element) * 2 + 2)
                 )  # кто такое Р для элеменитов, длиной 4 - це вопрос важный, и почему Р разные, хмхмххмхмхмм
    P[0, len(Mass_Element) * 2 + 1] = M
    # P=M
    print(np.matrix([[1], [1], [1], [1]]))

    aaa = Balcka(x, l, E, rho, h, w, Mass_node, Mass_Element, P)
    aaa.CreateMatrix_B(1)
    U = aaa.Solve([1])

