import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.linalg import solve


class thermal_conductivity:
    def __init__(self, x, l, Mass_node, Mass_Element, Mass_Element_Priming, Mass_node_T_air, Mass_node_T_water,
                 Lambda_Concrete, Lambda_Priming, T_air, T_water):
        print("__init__ termal")
        self.x = x
        self.l = l
        self.E = E
        # self.rho = rho
        self.g = 9.8
        # self.P = P
        self.T_air = T_air
        self.T_water = T_water
        self.Lambda_Concrete = Lambda_Concrete
        self.Lambda_Priming = Lambda_Priming
        self.Mass_node = Mass_node
        self.Mass_Element = Mass_Element
        self.Mass_Element_Priming = Mass_Element_Priming
        self.Mass_node_T_air = Mass_node_T_air
        self.Mass_node_T_water = Mass_node_T_water
        self.sym_eta = sym.Symbol('x')
        self.sym_ksi = sym.Symbol('y')
        self.sym_N_i = 1 - self.sym_ksi - self.sym_eta
        self.sym_N_j = self.sym_eta
        self.sym_N_k = self.sym_ksi
        self.B = np.matrix([[-1, 0, 1], [-1, 1, 0]])

        # self.k_e = self.CreateMatrix_k_e()
        # self.M_e = self.CreateMatrix_m_e()

    def CreateMatrix_С_e(self):
        print("CreateMatrix_M_e2")
        N = np.matrix([self.sym_N_i, self.sym_N_j, self.sym_N_k]).copy()
        CCC = sym.simplify(np.dot(N.T, N)).copy()
        C_e = np.zeros((4, 4)).copy()
        for i, row in enumerate(CCC):
            for j, element in enumerate(row):
                C_e[i, j] = sym.integrate(element, (self.sym_eta, -1, 1))
        C_e *= self.rho * self.l / 2
        # print(M_e, "M_e")
        return C_e

    def get_global(self, ind, Matrix_e):
        K = np.zeros((len(self.Mass_node), len(self.Mass_node)))
        K[ind[0], ind[0]] = Matrix_e[0, 0]
        K[ind[0], ind[1]] = Matrix_e[0, 1]
        K[ind[1], ind[0]] = Matrix_e[1, 0]
        K[ind[1], ind[2]] = Matrix_e[1, 2]
        K[ind[2], ind[1]] = Matrix_e[2, 1]
        K[ind[0], ind[2]] = Matrix_e[0, 2]
        K[ind[2], ind[0]] = Matrix_e[2, 0]
        K[ind[1], ind[1]] = Matrix_e[1, 1]
        K[ind[2], ind[2]] = Matrix_e[2, 2]
        return K

    def get_stiffness_matrix(self, coords, Lambda):
        Coord = np.matrix([[coords[0, 0], coords[0, 1]], [coords[1, 0], coords[1, 1]], [coords[2, 0], coords[2, 1]]])
        J = self.B * Coord
        # print(J)
        # print(self.B)
        Res_B = np.zeros((3, 2))
        # new_B = np.matrix([[-1, -1], [0, 1], [1, 0]])

        for i in range(3):
            vec_B = np.array([self.B[0, i], self.B[1, i]])
            # print('vec_B', vec_B)
            B = solve(J, vec_B)
            Res_B[i, 0] = B[0]
            Res_B[i, 1] = B[1]
            # print('B', B)
        print("Res_B", Res_B)
        print(np.dot(Res_B, Res_B.T))
        Ki = Lambda * np.dot(Res_B, Res_B.T) * np.linalg.det(J) / 2

        # print('Ki', Res_Ki)
        print('Ki', Ki)
        return Ki

    def Solve(self):  # , GU, time, dt
        print('solve')
        N = len(self.Mass_node)
        # print(N)
        K = np.zeros((N, N))

        for i in range(len(self.Mass_Element)):
            el = self.Mass_Element[i]
            # print('rrrrrr', i, el)
            Enodes = np.matrix([self.Mass_node[el[0]], self.Mass_node[el[1]], self.Mass_node[el[2]]])
            print(Enodes)
            if i in Mass_Element_Priming:
                # print("Lambda_Priming")
                _lambda = self.Lambda_Priming
            else:
                _lambda = self.Lambda_Concrete
                # print('Lambda_Concrete')
            Ki = self.get_stiffness_matrix(Enodes, _lambda)
            Ki = self.get_global(self.Mass_Element[i], Ki)
            K = K + Ki
        F = np.zeros((N, 1))

        for i in range(len(self.Mass_node_T_water)):
            K[self.Mass_node_T_water[i], :] = 0  # np.zeros((1, N))
            # K[:, self.Mass_node_T_water[i]] = 0
            K[self.Mass_node_T_water[i], self.Mass_node_T_water[i]] = 1
            F[self.Mass_node_T_water[i]] = self.T_water
        for i in range(len(self.Mass_node_T_air)):
            K[self.Mass_node_T_air[i], :] = 0
            # K[:,self.Mass_node_T_air[i]] = 0
            K[self.Mass_node_T_air[i], self.Mass_node_T_air[i]] = 1
            F[self.Mass_node_T_air[i]] = self.T_air


        print(K)
        # for i in range(len(K)):
        #     for j in range(len(K[i])):
        #         print(K[i][j], end=' ')
        T = solve(K, F)
        #np.savetxt('test1.txt', T, fmt='%.7f')
        np.savetxt('test1.txt', T, fmt='%.8e')
        print(T)
        print(len(T))


if __name__ == '__main__':
    node = open('nodes.txt', 'r')
    Mass_node = [[float(i) for i in (line.replace(" ", '').split(",")[1:])] for line in node.read().splitlines()]
    elem_all = open('elem_nodes.txt', 'r')
    Mass_Element = [[int(i) - 1 for i in line.replace(" ", '').split(",")[1:]] for line in elem_all.read().splitlines()]
    np.savetxt('MMMMM.txt', Mass_node, fmt='%d')
    # elem_Concrete = open('elem_Concrete.txt', 'r')
    # Mass_Element_Concrete=[int(i) for i in elem_Concrete.read().replace("\n", ',').replace(" ", '').split(",")]
    elem_Priming = open('elem_Concrete.txt', 'r')
    Mass_Element_Priming = [int(i) - 1 for i in elem_Priming.read().replace("\n", ',').replace(" ", '').split(",")]
    print('aaa', Mass_Element_Priming)
    node_T_air = open('node_T_air.txt', 'r')
    Mass_node_T_air = [int(i) - 1 for i in node_T_air.read().replace("\n", ',').replace(" ", '').split(",")]
    node_T_water = open('node_T_water.txt', 'r')
    Mass_node_T_water = [int(i) - 1 for i in node_T_water.read().replace("\n", ',').replace(" ", '').split(",")]

    print(Mass_node_T_water)

    print('start')

    Lambda_Concrete = 1.5
    Lambda_Priming = 1.75
    T_air = 25
    T_water = 5
    E = 2 * (10 ** 11)  # сталь
    # rho = 7800
    # M = 10000
    l = 0.1
    x = 1

    test = thermal_conductivity(x, l, Mass_node, Mass_Element, Mass_Element_Priming, Mass_node_T_air, Mass_node_T_water,
                                Lambda_Concrete, Lambda_Priming, T_air, T_water)
    test.Solve()
