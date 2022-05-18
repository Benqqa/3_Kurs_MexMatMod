import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.linalg import solve


class dynamic_beam:
    def __init__(self, x, l, E, rho, Mass_node, Mass_Element, P):
        print("__init__ Balcka")
        self.x = x
        self.l = l
        self.E = E
        self.rho = rho
        self.g = 9.8
        self.P = P
        self.type = self.Shveller(0.140, 0.058, 0.0049, 0.0049, 1)
        self.Mass_node = Mass_node
        self.Mass_Element = Mass_Element
        self.sym_eta = sym.Symbol('x')
        self.sym_N_i = (1 / 4) * ((1 - self.sym_eta) ** 2) * (2 + self.sym_eta)
        self.sym_N_i_theta = (1 / 8) * self.l * ((1 - self.sym_eta) ** 2) * (1 + self.sym_eta)
        self.sym_N_j = (1 / 4) * ((1 + self.sym_eta) ** 2) * (2 - self.sym_eta)
        self.sym_N_j_theta = -(1 / 8) * self.l * ((1 + self.sym_eta) ** 2) * (1 - self.sym_eta)

        self.k_e = self.CreateMatrix_k_e()
        self.M_e = self.CreateMatrix_M_e()

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
                self.S = 0.001200  # Площадь швеллера
                self.W_X = 0.000049943  # Момент сопротивления Wx
                self.W_y = 0.000008655  # Момент сопротивления Wy
                self.I_x = 0.00000349  # Момент инерции Ix
                self.I_y = 0.0000003725  # Момент инерции Iy
                self.i_x = 0.054  # Радиус инерции ix
                self.i_y = 0.018  # Радиус инерции iy
            elif self.variant == 2:
                self.I_x = 0.000003496003101266667  # Момент инерции Ix
                # self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (self.b - 0) ** 3) / 3
                self.I_y = 2.422756617262131 * (10 ** (-7))  # Момент инерции Iy
                self.I_x = (self.b * (self.h ** 3) - (self.b - self.t) * (self.h - 2 * self.t) ** 3) / 12
                self.S = 0.00092218  # Площадь швеллера
                # self.S = 0.001206  # Площадь швеллера
            elif self.variant == 3:
                self.I_x = 0.000003496  # Момент инерции Ix
                # self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (self.b - 0) ** 3) / 3
                self.I_y = 3.72506 * (10 ** -7)
                # self.I_x = (self.b * self.h ** 3 - (self.b - self.t) * (self.h - 2 * self.t) ** 3) / 12
                self.S = 0.00092218  # Площадь швеллера
                # self.S = 0.001206  # Площадь швеллера
            elif self.variant == 4:
                # self.I_x = 0.0000039051  # Момент инерции Ix
                self.I_x = 0.00000082952  # Момент инерции Ix
                # self.I_x=0.046666689718487794
                # self.I_x = 0.0000041  # Момент инерции Ix
                # self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (self.b - 0) ** 3) / 3
                self.I_y = 1.60618 * (10 ** -7)
                # self.I_x = (self.b * self.h ** 3 - (self.b - self.t) * (self.h - 2 * self.t) ** 3) / 12
                self.S = 0.00092218  # Площадь швеллера
                # self.S = 0.001206  # Площадь швеллера
            else:
                self.I_y = (self.h * 0 ** 3 - (self.h - 2 * self.t) * (0 - self.t) ** 3 + self.t * (
                        self.b - 0) ** 3) / 3

    def J(self):
        return self.type.I_x

    def N_i(self, eta):
        return 1 / 4 * ((1 - eta) ** 2) * (2 + eta)

    def d_N_i(self, eta):
        return (3 * (-1 + eta ** 2)) / 4

    def dd_N_i(self, eta):
        return (3 * eta) / 2

    def ddd_N_i(self, eta):
        return 3 / 2

    def N_i_theta(self, eta):
        return 1 / 8 * self.l * (1 + eta) * ((1 - eta) ** 2)

    def d_N_i_theta(self, eta):
        return (self.l * (-1 - 2 * eta + 3 * eta ** 2)) / 8

    def dd_N_i_theta(self, eta):
        return (self.l * (-1 + 3 * eta)) / 4

    def ddd_N_i_theta(self, eta):
        return (self.l * 3) / 4

    def N_j(self, eta):
        return 1 / 4 * (2 - eta) * ((1 + eta) ** 2)

    def d_N_j(self, eta):
        return (-3 * (-1 + eta ** 2)) / 4

    def dd_N_j(self, eta):
        return (-3 * eta) / 2

    def ddd_N_j(self, eta):
        return (-3) / 2

    def N_j_theta(self, eta):
        return -1 / 8 * self.l * (1 - eta) * ((1 + eta) ** 2)

    def d_N_j_theta(self, eta):
        return (self.l * (-1 + 2 * eta + 3 * eta ** 2)) / 8

    def dd_N_j_theta(self, eta):
        return (self.l * (1 + 3 * eta)) / 4

    def ddd_N_j_theta(self, eta):
        return (self.l * 3) / 4

    def CreateMatrix_k_e(self):  # локальная матрица жесткости
        k_e = np.zeros((4, 4))
        k_e[0, :] = [12, 6 * self.l, -12, 6 * self.l].copy()
        k_e[1, :] = [6 * self.l, 4 * (self.l ** 2), -6 * self.l, 2 * (self.l ** 2)].copy()
        k_e[2, :] = [-12, -6 * self.l, 12, - 6 * self.l].copy()
        k_e[3, :] = [6 * self.l, 2 * (self.l ** 2), -6 * self.l, 4 * self.l ** 2].copy()
        k_e *= self.E * self.J() / (self.l ** 3)
        print("k_e")
        # print(k_e)
        return k_e

    def CreateMatrix_m_e(self):  # локальная матрица жесткости
        k_e = np.zeros((4, 4))
        k_e[0, :] = [12, 0, 0, 0].copy()
        k_e[1, :] = [0, (self.l ** 2), 0, 0].copy()
        k_e[2, :] = [0, 0, 12, 0].copy()
        k_e[3, :] = [0, 0, 0, self.l ** 2].copy()
        k_e *= self.rho * self.type.S * self.l / 24
        print("k_e")
        # print(k_e)
        return k_e

    # def N_N_t(self, eta):
    #     # return self.N_i(eta) * self.N_i(eta) + self.N_i_theta(eta) * self.N_i_theta(eta) + self.N_j(eta) * self.N_j(
    #     #     eta) + self.N_j_theta(eta) * self.N_j_theta(eta)
    #     N = np.matrix([self.N_i(eta), self.N_i_theta(eta), self.N_j(eta), self.N_j_theta(eta)])
    #     return N * N.transpose()

    # def CreateMatrix_M_e(self):  # локальная Матрица массс
    #     print("CreateMatrix_M_e")
    #     M_e = np.zeros((1, len(self.Mass_node)))
    #     for index, node in enumerate(self.Mass_node):
    #         print(node, index)
    #         print(node[0, 0])
    #         print(self.Eta(node[0, 0]))
    #         v, err = integrate.quad(self.N_N_t, -1, 1)
    #         print(v, err)
    #         M_in_node = self.rho * self.l / 2 * v
    #         print(M_in_node, "M_in_node")
    #         M_e[0, index] = 2 * M_in_node
    #     M_e[0, 0] = M_e[0, 0] / 2
    #     M_e[0, len(self.Mass_node)-1] = M_e[0, len(self.Mass_node)-1] / 2
    #     print(M_e)
    #     return M_e
    def CreateMatrix_M_e(self):
        print("CreateMatrix_M_e2")
        N = np.matrix([self.sym_N_i, self.sym_N_i_theta, self.sym_N_j, self.sym_N_j_theta]).copy()
        MMM = sym.simplify(np.dot(N.T, N)).copy()
        M_e = np.zeros((4, 4)).copy()
        for i, row in enumerate(MMM):
            for j, element in enumerate(row):
                M_e[i, j] = sym.integrate(element, (self.sym_eta, -1, 1))
        M_e *= self.rho * self.l / 2   * self.type.S
        # print(M_e, "M_e")
        return M_e

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

    def get_global(self, Matrix_e):
        K = np.zeros((2 * len(self.Mass_node), 2 * len(self.Mass_node)))
        for i in range(len(self.Mass_Element)):
            nodes = self.Mass_Element[i, :]
            for k in range(2):
                for j in range(2):
                    K[2 * nodes[0, k], 2 * nodes[0, j]] += Matrix_e[2 * k, 2 * j].copy()
                    K[(2 * nodes[0, k]) + 1, (2 * nodes[0, j]) + 1] += Matrix_e[(2 * k) + 1, (2 * j) + 1].copy()
                    K[2 * nodes[0, k], (2 * nodes[0, j]) + 1] += Matrix_e[2 * k, (2 * j) + 1].copy()
                    K[(2 * nodes[0, k]) + 1, 2 * nodes[0, j]] += Matrix_e[(2 * k) + 1, 2 * j].copy()
        # print(K)
        # print('\n'.join(' '.join(str(col) for col in row) for row in K))
        return K

    def Eta(self, x):
        return 2 * x / self.l - 1

    def getDeff(self, U):
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
        self.Y = Y
        # print(Y, "Y")
        # print(X, "X")
        return [X, Y]

    def Solve(self, GU, time, dt):
        K_gl = self.get_global(self.k_e.copy()).copy()
        M_gl = self.get_global(self.M_e.copy()).copy()
        K_gu = self.setGu(GU, K_gl.copy()).copy()
        M_gu = self.setGu(GU, M_gl.copy()).copy()

        #
        # c=1
        # dt=(self.l**2)/c*(self.type.S/(48*self.J()))**(1/2)
        #
        F = self.getF_global().transpose().copy()
        print("K_gl")
        # print('\n'.join(' '.join(str(col) for col in row) for row in K_gl))
        print("M_gl")
        # print('\n'.join(' '.join(str(col) for col in row) for row in M_gl))

        [U, dU,ddU] = self.explicit_integration_scheme(time, dt, K_gu.copy(), M_gu.copy(), np.zeros((len(F), 1)), F)
        # [X, Y] = self.getDeff(U)
        # print(U)
        # plt.subplot(2, 2, 2)

        step_frame = int(time / dt / 10)
        step_frame2 = int(time / dt / 20)
        step_table = int(0.1 / dt)
        t1 = int(time / dt / 2)
        step_t1 = int(t1 / 10)
        print(step_table)
        # for index, el in enumerate(U):
        #     # if index<3:
        #     plt.plot(self.getDeff(el)[0])
        # plt.plot(self.getDeff(U[0])[0])
        # перемещениея
        plt.subplot(3, 2, 1)
        Legend_names = []
        for i in range(6):
            plt.plot(self.getDeff(U[i * step_table])[0])
            Legend_names.append('t=0.' + str(i))
            print("time-", i * step_table * dt, self.getDeff(U[i * step_table])[0], "U_step05")
            # print(i*step_table*dt)
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('U,[м]')
        plt.title('Перемещения вдоль оси OY t=0-0.5')
        plt.subplot(3, 2, 2)
        Legend_names = []
        for i in range(4):
            plt.plot(self.getDeff(U[(6 + i) * step_table])[0])
            Legend_names.append('t=0.' + str(6 + i))
            print("time-", (6 + i) * step_table * dt, self.getDeff(U[(6 + i) * step_table])[0], "U_step51")
        plt.plot(self.getDeff(U[len(U) - 1])[0])
        Legend_names.append('t=1')
        print("time-", (len(U) - 1) * dt, self.getDeff(U[len(U) - 1])[0], "U_step51")
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('U,[м]')
        plt.title('Перемещения вдоль оси OY t=0.5-1')

        # скорости
        plt.subplot(3, 2, 3)
        Legend_names = []
        for i in range(6):
            plt.plot(self.getDeff(dU[i * step_table])[0])
            Legend_names.append('t=0.' + str(i))
            print("time-", i * step_table * dt, self.getDeff(dU[i * step_table])[0], "dU_step05")
            # print(i*step_table*dt)
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('dU,[м/c]')
        plt.title('Скорости вдоль оси OY t=0-0.5')
        plt.subplot(3, 2, 4)
        Legend_names = []
        for i in range(4):
            plt.plot(self.getDeff(dU[(6 + i) * step_table])[0])
            Legend_names.append('t=0.' + str(6 + i))
            print("time-", (6 + i) * step_table * dt, self.getDeff(dU[(6 + i) * step_table])[0], "dU_step51")
        plt.plot(self.getDeff(dU[len(dU) - 1])[0])
        Legend_names.append('t=1.')
        print("time-", (len(dU) - 1) * dt, self.getDeff(dU[len(dU) - 1])[0], "dU_step51")
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('dU,[м/c]')
        plt.title('Скорости вдоль оси OY t=0.5-1')

        # ускорения
        plt.subplot(3, 2, 5)
        Legend_names = []
        for i in range(6):
            plt.plot(self.getDeff(ddU[i * step_table])[0])
            Legend_names.append('t=0.' + str(i))
            print("time-", i * step_table * dt, self.getDeff(ddU[i * step_table])[0], "ddU_step05")
            # print(i*step_table*dt)
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('ddU,[м/c^2]')
        plt.title('Ускорения вдоль оси OY t=0-0.5')
        plt.subplot(3, 2, 6)
        Legend_names = []
        for i in range(4):
            plt.plot(self.getDeff(ddU[(6 + i) * step_table])[0])
            Legend_names.append('t=0.' + str(6 + i))
            print("time-", (6 + i) * step_table * dt, self.getDeff(ddU[(6 + i) * step_table])[0], "ddU_step51")
        plt.plot(self.getDeff(ddU[len(ddU) - 1])[0])
        Legend_names.append('t=1')
        print("time-", (len(ddU) - 1) * dt, self.getDeff(ddU[len(ddU) - 1])[0], "ddU_step51")
        plt.legend(Legend_names, loc="lower left")
        plt.xlabel('x,[м]')
        plt.ylabel('ddU,[м/c^2]')
        plt.title('Ускорения вдоль оси OY t=0.5-1')

        plt.show()

    def explicit_integration_scheme2(self, time, dt, K_gl, M_gl, U_0, F_0):
        # print(F_0, "F_0")

        arr = np.arange(0, time, dt)
        dF = F_0.copy() / (time / dt / 2)
        # print(dF, 'dF')
        F_0 = np.zeros((len(F_0), 1))
        U_new = np.zeros((len(arr), len(U_0)))
        DU_new = np.zeros((len(arr), len(U_0)))
        # print(U_new, "U_new")
        # print(U_0)
        U_new[0] = U_0.T.copy()
        ddU = np.zeros((len(U_0), 1))
        dU_new = np.zeros((len(U_0), 1))
        ddU_new = np.zeros((len(U_0), 1))
        dU = np.zeros((len(U_0), 1))
        a0 = 1 / (dt ** 2)
        a1 = 1 / (2 * dt)
        a2 = 2 * a0
        a3 = 1 / a2
        U_old = np.matrix([U_new[0].copy()]).copy() - dt * np.matrix([U_new[0].copy()]).copy() + a3 * np.matrix(
            [U_new[0].copy()]).copy()
        M_e2 = M_gl * a0
        M_tr = self.Gauss(M_e2)

        for index, el in enumerate(arr):
            U = np.matrix([U_new[index - 1].copy()]).copy()
            if index == 0:
                R = -1 * (K_gl - a2 * M_gl) * U_old.T - (a0 * M_gl) * U.T
                U_new[index, :] = np.squeeze(np.asarray(np.linalg.inv(M_e2.copy()) * R))
            else:
                R = F_0 - 1 * (K_gl - a2 * M_gl) * np.matrix([U_new[index - 2].copy()]).copy().T - (a0 * M_gl) * U.T
                U_new[index, :] = np.squeeze(np.asarray(np.linalg.inv(M_e2.copy()) * R))
        for el in U_new:
            print(self.getDeff(el)[0], 'U2')
        return [U_new, DU_new]

    def explicit_integration_scheme3(self, time, dt, K_gl, M_gl, U_0, F_0):
        # print(F_0, "F_0")

        arr = np.arange(0, time, dt)
        dF = F_0.copy() / (time / dt / 2)
        # print(dF, 'dF')
        F_0 = np.zeros((len(F_0), 1))
        U_new = np.zeros((len(arr), len(U_0)))
        DU_new = np.zeros((len(arr), len(U_0)))
        # print(U_new, "U_new")
        # print(U_0)
        U_new[0] = U_0.T.copy()
        ddU = np.zeros((len(U_0), 1))
        dU_new = np.zeros((len(U_0), 1))
        ddU_new = np.zeros((len(U_0), 1))
        dU = np.zeros((len(U_0), 1))

        for index, el in enumerate(arr):
            U = np.matrix([U_new[index - 1].copy()]).copy()
            # U_new[index, :] = np.squeeze(np.asarray(np.linalg.inv(M_gl.copy())*(dt**2)/2*(F_0-(K_gl-2/(dt**2) * M_gl)*U.T)))
            U_new[index, :] = np.squeeze(
                np.asarray(solve(2 / (dt ** 2) * M_gl, (F_0 - (K_gl - 2 / (dt ** 2) * M_gl) * U.T))))
            # print(F_0, "F_03")
            if index * dt >= time / 2:
                F_0 = np.zeros((len(F_0), 1))
            else:
                F_0 += dF.copy()
        for el in U_new:
            print(self.getDeff(el)[0], 'U3')
        return [U_new, DU_new]

    def Gauss(self, A):
        for k in range(len(A) - 1):
            At = A.copy()
            for i in range(len(A)):
                for j in range(len(A)):
                    if (i <= k):
                        A[i][j] = At[i][j]
                    elif (i > k and j > k):
                        A[i][j] = round(At[i][j] - (At[i][k] / At[k][k]) * At[k][j], 4)
                    elif (i > k and j <= k):
                        A[i][j] = 0
        return A

    def explicit_integration_scheme(self, time, dt, K_gl, M_gl, U_0, F_0):
        # print(F_0, "F_0")

        arr = np.arange(0, time, dt)
        dF = F_0.copy() / (time / dt / 2)
        # print(dF, 'dF')
        F_0 = np.zeros((len(F_0), 1))
        U_new = np.zeros((len(arr), len(U_0)))
        DU_new = np.zeros((len(arr), len(U_0)))
        DDU_new = np.zeros((len(arr), len(U_0)))
        # print(U_new, "U_new")
        # print(U_0)
        U_new[0] = U_0.T.copy()
        ddU = np.zeros((len(U_0), 1))
        dU_new = np.zeros((len(U_0), 1))
        ddU_new = np.zeros((len(U_0), 1))
        dU = np.zeros((len(U_0), 1))
        # print(ddU, "ddU")
        # print(dU, "dU")
        m = np.linalg.inv(M_gl.copy())
        for index, el in enumerate(arr):
            U = np.matrix([U_new[index - 1].copy()]).copy()
            ddU_new = m * (F_0.copy() - (K_gl.copy() * U.T.copy()))
            dU_new = dU.copy() + (dt * ddU_new.copy())
            DU_new[index, :] = np.squeeze(np.asarray(dU_new.T.copy()))
            DDU_new[index, :] = np.squeeze(np.asarray(ddU_new.T.copy()))
            U_new[index, :] = np.squeeze(np.asarray(U.T.copy() + (dt * dU_new.copy())))
            dU = dU_new.copy()
            if index * dt >= time / 2:
                F_0 = np.zeros((len(F_0), 1))
            else:
                F_0 += dF.copy()
        # for el in U_new:
        #     print(self.getDeff(el)[0], 'U')
        # print(U_new, "U_new")
        return [U_new, DU_new,DDU_new]

    def explicit_integration_scheme4(self, time, dt, K_gl, M_gl, U_0, F_0):
        # print(F_0, "F_0")

        arr = np.arange(0, time, dt)
        dF = F_0.copy() / (time / dt / 2)
        # print(dF, 'dF')
        F_0 = np.zeros((len(F_0), 1))
        U_new = np.zeros((len(arr), len(U_0)))
        DU_new = np.zeros((len(arr), len(U_0)))
        # print(U_new, "U_new")
        # print(U_0)
        U_new[0] = U_0.T.copy()
        ddU = np.zeros((len(U_0), 1))
        dU_new = np.zeros((len(U_0), 1))
        ddU_new = np.zeros((len(U_0), 1))
        dU = np.zeros((len(U_0), 1))
        # print(ddU, "ddU")
        # print(dU, "dU")
        m = np.linalg.inv(M_gl.copy())
        for index, el in enumerate(arr):
            U = np.matrix([U_new[index - 1].copy()]).copy()
            ddU_new = m * (F_0.copy() - (K_gl.copy() * U.T.copy()))
            # ddU_new = np.linalg.inv(M_gl.copy()) * (F_0.copy() - (K_gl.copy() * U.T.copy()))
            # print(ddU_new, "ddU_new")
            dU_new = dU.copy() + (dt * ddU.copy())
            # if(index+1 <= len(arr)-1):
            DU_new[index, :] = np.squeeze(np.asarray(dU_new.T.copy()))
            U_new[index, :] = np.squeeze(np.asarray(U.T.copy() + (dt * dU_new.copy())))
            dU = dU_new.copy()
            ddU = ddU_new.copy()
            # print(F_0, "F_0")
            if index * dt >= time / 2:
                F_0 = np.zeros((len(F_0), 1))
            else:
                F_0 += dF.copy()
        for el in U_new:
            print(self.getDeff(el)[0], 'U')
        # print(U_new, "U_new")
        return [U_new, DU_new]

    def setGu(self, GU, Matrix):
        for i in range(len(Matrix)):
            for j in range(len(GU)):
                if i == 2 * (GU[j] - 1):
                    Matrix[i, :] = 0
                    Matrix[:, i] = 0
                    Matrix[i + 1, :] = 0
                    Matrix[:, i + 1] = 0
                    Matrix[i, i] = 1
                    Matrix[i + 1, i + 1] = 1

        return Matrix

    def getFfromNandU(self, U, number):
        X = np.zeros((int(len(U) / 2 - 1), 22))
        U_new = np.zeros((int(len(U) / 2 - 1), 22))
        for i in range(int(len(U) / 2 - 1)):
            X[i] = np.linspace(i * self.l, (i + 1) * self.l, 22)
        for i in range(int(len(U) / 2 - 1)):
            U_new[i] = self.ddd_N_i(self.Eta(X[number])) * U[2 * i] + self.ddd_N_i_theta(self.Eta(X[number])) * U[
                2 * i + 1] + self.ddd_N_j(self.Eta(X[number])) * U[2 * i + 2] + self.ddd_N_j_theta(
                self.Eta(X[number])) * U[
                           2 * i + 3]
        U_res = self.E * self.J() * 8 * U_new / (self.l ** 3)
        # print(U_res, "U_F2")
        return [X, U_res]

    def getMfromNandU(self, U, number):
        X = np.zeros((int(len(U) / 2 - 1), 22))
        U_new = np.zeros((int(len(U) / 2 - 1), 22))
        for i in range(int(len(U) / 2 - 1)):
            X[i] = np.linspace(i * self.l, (i + 1) * self.l, 22)
            # print(np.linspace(i*self.l, (i+1)*self.l, 20))
        for i in range(int(len(U) / 2 - 1)):
            U_new[i] = self.dd_N_i(self.Eta(X[number])) * U[2 * i] + self.dd_N_i_theta(self.Eta(X[number])) * U[
                2 * i + 1] + self.dd_N_j(self.Eta(X[number])) * U[2 * i + 2] + self.dd_N_j_theta(self.Eta(X[number])) * \
                       U[
                           2 * i + 3]

        print(self.E * self.J() * 4 * U_new / (self.l ** 2), "U_M")
        return [X, self.E * self.J() * 4 * U_new / (self.l ** 2)]

    def getUfromN(self, U, number):
        X = np.zeros((int(len(U) / 2 - 1), 20))
        U_new = np.zeros((int(len(U) / 2 - 1), 20))
        for i in range(int(len(U) / 2 - 1)):
            X[i] = np.linspace(i * self.l, (i + 1) * self.l, 20)
            # print(np.linspace(i*self.l, (i+1)*self.l, 20))
        for i in range(int(len(U) / 2 - 1)):
            U_new[i] = self.N_i(self.Eta(X[number])) * U[2 * i] + self.N_i_theta(self.Eta(X[number])) * U[
                2 * i + 1] + self.N_j(self.Eta(X[number])) * U[2 * i + 2] + self.N_j_theta(self.Eta(X[number])) * U[
                           2 * i + 3]

            # print(U_new[i], "U")
        return [X, U_new]


if __name__ == '__main__':
    E = 2 * (10 ** 11)  # сталь
    rho = 7800
    M = 10000
    l = 0.1
    x = 1
    Mass_node = np.matrix(
        [[0, 0], [0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 0], [0.6, 0], [0.7, 0], [0.8, 0],
         [0.9, 0], [1, 0]])
    Mass_Element = np.matrix(
        [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10]])
    P = np.zeros((1, len(Mass_Element) * 2 + 2))  # нагрузка
    P[0, 21] = M
    new_Beam = dynamic_beam(x, l, E, rho, Mass_node, Mass_Element, P)
    GU = [1]
    new_Beam.Solve(GU, 1, 0.000001)
