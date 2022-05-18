import math

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.linalg import solve


class thermal_conductivity:
    def __init__(self, x, l, Mass_node, Mass_Element, Mass_Element_Priming, Mass_Element_C, Mass_Element_C2, E_Priming,
                 nu_Priming, E_C,
                 nu_C, E_C2, nu_C2, Mass_node_B_X, Mass_node_B_Y):
        print("__init__ termal")
        self.x = x
        self.l = l
        self.E_Priming = E_Priming
        self.E_C = E_C
        self.E_C2 = E_C2
        self.nu_Priming = nu_Priming
        self.nu_C = nu_C
        self.nu_C2 = nu_C2
        # self.rho = rho
        self.g = 9.8
        self.B_X = -1 * self.g
        # self.P = P
        # self.T_air = T_air
        # self.T_water = T_water
        # self.Lambda_Concrete = Lambda_Concrete
        # self.Lambda_Priming = Lambda_Priming

        self.Mass_node_B_X = Mass_node_B_X
        self.Mass_node_B_Y = Mass_node_B_Y

        self.Mass_node = Mass_node
        self.Mass_Element = Mass_Element
        self.Mass_Element_Priming = Mass_Element_Priming
        self.Mass_Element_C = Mass_Element_C
        self.Mass_Element_C2 = Mass_Element_C2
        # self.Mass_node_T_air = Mass_node_T_air
        # self.Mass_node_T_water = Mass_node_T_water
        self.sym_eta = sym.Symbol('x')
        self.sym_ksi = sym.Symbol('y')
        self.sym_N_i = 1 - self.sym_ksi - self.sym_eta
        self.sym_N_j = self.sym_eta
        self.sym_N_k = self.sym_ksi
        self.B = np.matrix([[-1, 0, 0, 0, 1, 0],
                            [-1, 0, 1, 0, 0, 0],
                            [-1, -1, 1, 0, 0, 1]])

        self.B_t = np.matrix([[-1, 0, 1],
                              [-1, 1, 0]])

        # self.D = (self.E * (1 - self.nu) / ((1 - 2 * self.nu) * (1 + self.nu))) * np.matrix(
        #     [[1, self.nu / (1 - self.nu), 0], [self.nu / (1 - self.nu), 1, 0],
        #      [0, 0, (1 - 2 * self.nu) / (2 * (1 - self.nu))]])
        # self.k_e = self.CreateMatrix_k_e()
        # self.M_e = self.CreateMatrix_m_e()

    def get_global(self, ind, Matrix_e):
        print('Matrix_e', Matrix_e)  # проверить ее размер
        print(len(Matrix_e))
        K = np.zeros((len(self.Mass_node) * 2, len(self.Mass_node) * 2))
        print(len(K))

        # for i in range(len(Matrix_e)):
        #     for j in range(len(Matrix_e) - 1):
        #         if (i % 2 == 0):
        #             firtst = ind[math.floor(i / 2)] * 2
        #         else:
        #             firtst = (ind[math.floor(i / 2)] * 2) + 1
        #         if (j % 2 == 0):
        #             twice = ind[math.floor(j / 2)] * 2
        #         else:
        #             twice = (ind[math.floor(j / 2)] * 2) + 1
        #         K[firtst, twice] = Matrix_e[i, j]

        K[ind[0] * 2, ind[0] * 2] = Matrix_e[0, 0]
        K[(ind[0] * 2) + 1, (ind[0] * 2) + 1] = Matrix_e[0 + 1, 0 + 1]
        K[(ind[0] * 2), (ind[0] * 2) + 1] = Matrix_e[0, 0 + 1]
        K[(ind[0] * 2) + 1, (ind[0] * 2)] = Matrix_e[0 + 1, 0]

        K[ind[0] * 2, ind[1] * 2] = Matrix_e[0, 1 * 2]
        K[(ind[0] * 2) + 1, (ind[1] * 2) + 1] = Matrix_e[0 + 1, 1 * 2 + 1]
        K[(ind[0] * 2), (ind[1] * 2) + 1] = Matrix_e[0, 1 * 2 + 1]
        K[(ind[0] * 2) + 1, (ind[1] * 2)] = Matrix_e[0 + 1, 1 * 2]

        K[ind[1] * 2, ind[0] * 2] = Matrix_e[1 * 2, 0]
        K[(ind[1] * 2) + 1, (ind[0] * 2) + 1] = Matrix_e[1 * 2 + 1, 0 + 1]
        K[(ind[1] * 2), (ind[0] * 2) + 1] = Matrix_e[1 * 2, 0 + 1]
        K[(ind[1] * 2) + 1, (ind[0] * 2)] = Matrix_e[1 * 2 + 1, 0]

        K[ind[1] * 2, ind[2] * 2] = Matrix_e[1 * 2, 2 * 2]
        K[(ind[1] * 2) + 1, (ind[2] * 2) + 1] = Matrix_e[1 * 2 + 1, 2 * 2 + 1]
        K[(ind[1] * 2), (ind[2] * 2) + 1] = Matrix_e[1 * 2, 2 * 2 + 1]
        K[(ind[1] * 2) + 1, (ind[2] * 2)] = Matrix_e[1 * 2 + 1, 2 * 2]

        K[ind[2] * 2, ind[1] * 2] = Matrix_e[2 * 2, 1 * 2]
        K[(ind[2] * 2) + 1, (ind[1] * 2) + 1] = Matrix_e[2 * 2 + 1, 1 * 2 + 1]
        K[(ind[2] * 2), (ind[1] * 2) + 1] = Matrix_e[2 * 2, 1 * 2 + 1]
        K[(ind[2] * 2) + 1, (ind[1] * 2)] = Matrix_e[2 * 2 + 1, 1 * 2]

        K[ind[0] * 2, ind[2] * 2] = Matrix_e[0, 2 * 2]
        K[(ind[0] * 2) + 1, (ind[2] * 2) + 1] = Matrix_e[0 + 1, 2 * 2 + 1]
        K[(ind[0] * 2), (ind[2] * 2) + 1] = Matrix_e[0, 2 * 2 + 1]
        K[(ind[0] * 2) + 1, (ind[2] * 2)] = Matrix_e[0 + 1, 2 * 2]

        K[ind[2] * 2, ind[0] * 2] = Matrix_e[2 * 2, 0]
        K[(ind[2] * 2) + 1, (ind[0] * 2) + 1] = Matrix_e[2 * 2 + 1, 0 + 1]
        K[(ind[2] * 2), (ind[0] * 2) + 1] = Matrix_e[2 * 2, 0 + 1]
        K[(ind[2] * 2) + 1, (ind[0] * 2)] = Matrix_e[2 * 2 + 1, 0]

        K[ind[1] * 2, ind[1] * 2] = Matrix_e[1 * 2, 1 * 2]
        K[(ind[1] * 2) + 1, (ind[1] * 2) + 1] = Matrix_e[1 * 2 + 1, 1 * 2 + 1]
        K[(ind[1] * 2), (ind[1] * 2) + 1] = Matrix_e[1 * 2, 1 * 2 + 1]
        K[(ind[1] * 2) + 1, (ind[1] * 2)] = Matrix_e[1 * 2 + 1, 1 * 2]

        K[ind[2] * 2, ind[2] * 2] = Matrix_e[2 * 2, 2 * 2]
        K[(ind[2] * 2) + 1, (ind[2] * 2) + 1] = Matrix_e[2 * 2 + 1, 2 * 2 + 1]
        K[(ind[2] * 2), (ind[2] * 2) + 1] = Matrix_e[2 * 2, 2 * 2 + 1]
        K[(ind[2] * 2) + 1, (ind[2] * 2)] = Matrix_e[2 * 2 + 1, 2 * 2]

        # K[ind[0], ind[1]] = Matrix_e[0, 1]
        # K[ind[1], ind[0]] = Matrix_e[1, 0]
        # K[ind[1], ind[2]] = Matrix_e[1, 2]
        # K[ind[2], ind[1]] = Matrix_e[2, 1]
        # K[ind[0], ind[2]] = Matrix_e[0, 2]
        # K[ind[2], ind[0]] = Matrix_e[2, 0]
        # K[ind[1], ind[1]] = Matrix_e[1, 1]
        # K[ind[2], ind[2]] = Matrix_e[2, 2]

        return K
    def get_J(self, coords):
        Coord = np.matrix([[coords[0, 0], coords[0, 1]], [coords[1, 0], coords[1, 1]], [coords[2, 0], coords[2, 1]]])
        J = self.B_t * Coord
        return J
    def get_stiffness_matrix(self, coords, D):
        Coord = np.matrix([[coords[0, 0], coords[0, 1]], [coords[1, 0], coords[1, 1]], [coords[2, 0], coords[2, 1]]])
        print('Coord', Coord)
        print('self.B', self.B)
        J = self.B_t * Coord  # разобраться какой B использовать
        print('J', J)

        # print(self.B)
        Res_B = np.zeros((3, 6))
        # new_B = np.matrix([[-1, -1], [0, 1], [1, 0]])

        for i in range(3):
            print('i', i)
            vec_B = np.array([self.B_t[0, i], self.B_t[1, i]])
            # print('vec_B', vec_B)
            B = solve(J, vec_B)

            Res_B[0, i * 2] = B[0]
            Res_B[1, (i * 2) + 1] = B[1]
            Res_B[2, i * 2] = B[1]
            Res_B[2, (i * 2) + 1] = B[0]
            # print('B', B)
        print("Res_B", Res_B)
        print(np.dot(Res_B, Res_B.T))
        Ki = np.dot(np.dot(Res_B.T, D), Res_B) * np.linalg.det(J) / 2

        # print('Ki', Res_Ki)
        print('Ki', Ki)
        return Ki

    def Solve(self):  # , GU, time, dt
        print('solve')
        N = len(self.Mass_node)
        # print(N)
        K = np.zeros((N * 2, N * 2))
        D_Prinmig = (self.E_Priming * (1 - self.nu_Priming) / (
                (1 - 2 * self.nu_Priming) * (1 + self.nu_Priming))) * np.matrix(
            [[1, self.nu_Priming / (1 - self.nu_Priming), 0],
             [self.nu_Priming / (1 - self.nu_Priming), 1, 0],
             [0, 0, (1 - 2 * self.nu_Priming) / (2 * (1 - self.nu_Priming))]
             ])
        D_C = (self.E_C * (1 - self.nu_C) / (
                (1 - 2 * self.nu_C) * (1 + self.nu_C))) * np.matrix(
            [[1, self.nu_C / (1 - self.nu_C), 0], [self.nu_C / (1 - self.nu_C), 1, 0],
             [0, 0, (1 - 2 * self.nu_C) / (2 * (1 - self.nu_C))]])
        D_C2 = (self.E_C2 * (1 - self.nu_C2) / (
                (1 - 2 * self.nu_C2) * (1 + self.nu_C2))) * np.matrix(
            [[1, self.nu_C2 / (1 - self.nu_C2), 0], [self.nu_C2 / (1 - self.nu_C2), 1, 0],
             [0, 0, (1 - 2 * self.nu_C2) / (2 * (1 - self.nu_C2))]])

        for i in range(len(self.Mass_Element)):
            el = self.Mass_Element[i]
            # print('rrrrrr', i, el)
            Enodes = np.matrix([self.Mass_node[el[0]], self.Mass_node[el[1]], self.Mass_node[el[2]]])
            print(Enodes)

            if i in self.Mass_Element_Priming:
                # print("Lambda_Priming")
                D = D_Prinmig
            elif i in self.Mass_Element_C:
                D = D_C
            else:
                D = D_C2
                # print('Lambda_Concrete')

            Ki = self.get_stiffness_matrix(Enodes, D)
            Ki = self.get_global(self.Mass_Element[i], Ki)
            K = K + Ki
        F = np.zeros((N * 2, 1))

        def p_hydro(y):
            return 9800 * (115.171875 - y)

        def S_trapeze(a, b, h):
            return 1/2*(a + b)*h


        def distance(coord_1, coord_2):
            return float(math.sqrt((coord_1[0] - coord_2[0]) ** 2 + (coord_1[1] - coord_2[1]) ** 2))



        # Гу x
        for i in range(len(self.Mass_node_B_X)):
            print('ii', i)
            K[self.Mass_node_B_X[i] * 2, :] = 0
            # K[(self.Mass_node_B_X[i] * 2) + 1, :] = 0
            # K[:,self.Mass_node_T_air[i]] = 0
            K[self.Mass_node_B_X[i] * 2, self.Mass_node_B_X[i] * 2] = 1
            # K[(self.Mass_node_B_X[i] * 2) + 1, (self.Mass_node_B_X[i] * 2) + 1] = 1
            F[self.Mass_node_B_X[i] * 2] = 0
        # ГУ y
        for i, index in enumerate(range(len(self.Mass_node_B_Y))):
            print('ii', i)
            print('index', index)
            # K[self.Mass_node_B_Y[i] * 2, :] = 0
            K[(self.Mass_node_B_Y[i] * 2) + 1, :] = 0
            # K[:,self.Mass_node_T_air[i]] = 0
            # K[self.Mass_node_B_Y[i] * 2, self.Mass_node_B_Y[i] * 2] = 1
            K[(self.Mass_node_B_Y[i] * 2) + 1, (self.Mass_node_B_Y[i] * 2) + 1] = 1
            F[(self.Mass_node_B_Y[i] * 2) + 1] = 0
        nodes = self.Mass_node
        # # давление на горизонтальное верхнее
        # F[2 * 12] += S_trapeze(p_hydro(nodes[12][1]), p_hydro(1 / 2 * (nodes[12][1] + nodes[9][1])),
        #                        1 / 2 * (nodes[12][1] - nodes[9][1]))
        # F[2 * 9] += S_trapeze(p_hydro(nodes[9][1]), p_hydro(1 / 2 * (nodes[9][1] + nodes[12][1])),
        #                        1 / 2 * (nodes[12][1] - nodes[9][1]))
        # давление на наклонное чудо
        F[2 * 13] += S_trapeze(p_hydro(nodes[13][1]), p_hydro(1 / 2 * (nodes[13][1] + nodes[135][1])),
                                1 / 2 * distance(nodes[13], nodes[135])) * 0.995
        F[2 * 13+1] += -S_trapeze(p_hydro(nodes[13][1]), p_hydro(1 / 2 * (nodes[13][1] + nodes[135][1])),
                                 1 / 2 * distance(nodes[13], nodes[135])) * 1.395
        F[2*135] += S_trapeze(p_hydro(nodes[135][1]), p_hydro(1/2*(nodes[135][1] + nodes[13][1])), 1/2*distance(nodes[13], nodes[135]))*0.995
        F[2 * 135+1] += -S_trapeze(p_hydro(nodes[135][1]), p_hydro(1 / 2 * (nodes[135][1] + nodes[13][1])),
                               1 / 2 * distance(nodes[135], nodes[13])) * 1.395
        F[2 * 135] += S_trapeze(p_hydro(nodes[135][1]), p_hydro(1 / 2 * (nodes[135][1] + nodes[136][1])),
                                1 / 2 * distance(nodes[136], nodes[135])) * 0.995
        F[2 * 135 + 1] += -S_trapeze(p_hydro(nodes[135][1]), p_hydro(1 / 2 * (nodes[135][1] + nodes[136][1])),
                                     1 / 2 * distance(nodes[135], nodes[136])) * 1.395

        F[2 * 136] += S_trapeze(p_hydro(nodes[136][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[135][1])),
                                1 / 2 * distance(nodes[136], nodes[135])) * 0.995
        F[2 * 136+1] += -S_trapeze(p_hydro(nodes[136][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[135][1])),
                                 1 / 2 * distance(nodes[136], nodes[135])) * 1.395
        F[2 * 136] += S_trapeze(p_hydro(nodes[136][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[137][1])),
                                1 / 2 * distance(nodes[136], nodes[137])) * 0.995
        F[2 * 136+1] += -S_trapeze(p_hydro(nodes[136][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[137][1])),
                                 1 / 2 * distance(nodes[136], nodes[137])) * 1.395

        F[2 * 137] += S_trapeze(p_hydro(nodes[137][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[137][1])),
                                1 / 2 * distance(nodes[136], nodes[137])) * 0.995
        F[2 * 137+1] += -S_trapeze(p_hydro(nodes[137][1]), p_hydro(1 / 2 * (nodes[136][1] + nodes[137][1])),
                                 1 / 2 * distance(nodes[136], nodes[137])) * 1.395
        F[2 * 137] += S_trapeze(p_hydro(nodes[137][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[137][1])),
                                1 / 2 * distance(nodes[138], nodes[137])) * 0.995
        F[2 * 137+1] += -S_trapeze(p_hydro(nodes[137][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[137][1])),
                                 1 / 2 * distance(nodes[138], nodes[137])) * 1.395

        F[2 * 138] += S_trapeze(p_hydro(nodes[138][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[137][1])),
                                1 / 2 * distance(nodes[138], nodes[137])) * 0.995
        F[2 * 138 + 1] += -S_trapeze(p_hydro(nodes[138][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[137][1])),
                                     1 / 2 * distance(nodes[138], nodes[137])) * 1.395
        F[2 * 138] += S_trapeze(p_hydro(nodes[138][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[139][1])),
                                1 / 2 * distance(nodes[138], nodes[139])) * 0.995
        F[2 * 138 + 1] += -S_trapeze(p_hydro(nodes[138][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[139][1])),
                                     1 / 2 * distance(nodes[138], nodes[139])) * 1.395

        F[2 * 139] += S_trapeze(p_hydro(nodes[139][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[139][1])),
                                1 / 2 * distance(nodes[138], nodes[139])) * 0.995
        F[2 * 139 + 1] += -S_trapeze(p_hydro(nodes[139][1]), p_hydro(1 / 2 * (nodes[138][1] + nodes[139][1])),
                                     1 / 2 * distance(nodes[138], nodes[139])) * 1.395
        F[2 * 139] += S_trapeze(p_hydro(nodes[139][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[139][1])),
                                1 / 2 * distance(nodes[140], nodes[139])) * 0.995
        F[2 * 139 + 1] += -S_trapeze(p_hydro(nodes[139][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[139][1])),
                                     1 / 2 * distance(nodes[140], nodes[139])) * 1.395

        F[2 * 140] += S_trapeze(p_hydro(nodes[140][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[139][1])),
                                1 / 2 * distance(nodes[140], nodes[139])) * 0.995
        F[2 * 140 + 1] += -S_trapeze(p_hydro(nodes[140][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[139][1])),
                                     1 / 2 * distance(nodes[140], nodes[139])) * 1.395
        F[2 * 140] += S_trapeze(p_hydro(nodes[140][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[141][1])),
                                1 / 2 * distance(nodes[140], nodes[141])) * 0.995
        F[2 * 140 + 1] += -S_trapeze(p_hydro(nodes[140][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[141][1])),
                                     1 / 2 * distance(nodes[140], nodes[141])) * 1.395

        F[2 * 141] += S_trapeze(p_hydro(nodes[141][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[141][1])),
                                1 / 2 * distance(nodes[140], nodes[141])) * 0.995
        F[2 * 141 + 1] += -S_trapeze(p_hydro(nodes[141][1]), p_hydro(1 / 2 * (nodes[140][1] + nodes[141][1])),
                                     1 / 2 * distance(nodes[140], nodes[141])) * 1.395
        F[2 * 141] += S_trapeze(p_hydro(nodes[141][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[141][1])),
                                1 / 2 * distance(nodes[142], nodes[141])) * 0.995
        F[2 * 141 + 1] += -S_trapeze(p_hydro(nodes[141][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[141][1])),
                                     1 / 2 * distance(nodes[142], nodes[141])) * 1.395

        F[2 * 142] += S_trapeze(p_hydro(nodes[142][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[141][1])),
                                1 / 2 * distance(nodes[142], nodes[141])) * 0.995
        F[2 * 142 + 1] += -S_trapeze(p_hydro(nodes[142][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[141][1])),
                                     1 / 2 * distance(nodes[142], nodes[141])) * 1.395
        F[2 * 142] += S_trapeze(p_hydro(nodes[142][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[143][1])),
                                1 / 2 * distance(nodes[142], nodes[143])) * 0.995
        F[2 * 142 + 1] += -S_trapeze(p_hydro(nodes[142][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[143][1])),
                                     1 / 2 * distance(nodes[142], nodes[143])) * 1.395

        F[2 * 143] += S_trapeze(p_hydro(nodes[143][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[143][1])),
                                1 / 2 * distance(nodes[142], nodes[143])) * 0.995
        F[2 * 143 + 1] += -S_trapeze(p_hydro(nodes[143][1]), p_hydro(1 / 2 * (nodes[142][1] + nodes[143][1])),
                                     1 / 2 * distance(nodes[142], nodes[143])) * 1.395
        F[2 * 143] += S_trapeze(p_hydro(nodes[143][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[143][1])),
                                1 / 2 * distance(nodes[144], nodes[143])) * 0.995
        F[2 * 143 + 1] += -S_trapeze(p_hydro(nodes[143][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[143][1])),
                                     1 / 2 * distance(nodes[144], nodes[143])) * 1.395

        F[2 * 144] += S_trapeze(p_hydro(nodes[144][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[143][1])),
                                1 / 2 * distance(nodes[144], nodes[143])) * 0.995
        F[2 * 144 + 1] += -S_trapeze(p_hydro(nodes[144][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[143][1])),
                                     1 / 2 * distance(nodes[144], nodes[143])) * 1.395
        F[2 * 144] += S_trapeze(p_hydro(nodes[144][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[5][1])),
                                1 / 2 * distance(nodes[144], nodes[5])) * 0.995
        F[2 * 144 + 1] += -S_trapeze(p_hydro(nodes[144][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[5][1])),
                                     1 / 2 * distance(nodes[144], nodes[5])) * 1.395

        F[2 * 5] += S_trapeze(p_hydro(nodes[5][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[5][1])),
                                1 / 2 * distance(nodes[144], nodes[5])) * 0.995
        F[2 * 5 + 1] += -S_trapeze(p_hydro(nodes[5][1]), p_hydro(1 / 2 * (nodes[144][1] + nodes[5][1])),
                                     1 / 2 * distance(nodes[144], nodes[5])) * 1.395
        # давление на горизонтальное ниженне
        p_vert = p_hydro(0)
        our_vert_value = p_vert * (
            abs(nodes[42][0] - nodes[5][0]))
        F[5 * 2 + 1] += -our_vert_value / 2
        F[42 * 2 + 1] += -our_vert_value
        F[43 * 2 + 1] += -our_vert_value
        F[44 * 2 + 1] += -our_vert_value
        F[45 * 2 + 1] += -our_vert_value
        F[46 * 2 + 1] += -our_vert_value
        F[47 * 2 + 1] += -our_vert_value
        F[48 * 2 + 1] += -our_vert_value
        F[49 * 2 + 1] += -our_vert_value
        F[50 * 2 + 1] += -our_vert_value
        F[51 * 2 + 1] += -our_vert_value
        F[6 * 2 + 1] += -our_vert_value/2

        # for i in range(len(self.Mass_node_B_X)):
        #     print('ii', i)
        #     K[self.Mass_node_B_X[i] * 2, :] = 0
        #     K[(self.Mass_node_B_X[i] * 2) + 1, :] = 0
        #     # K[:,self.Mass_node_T_air[i]] = 0
        #     K[self.Mass_node_B_X[i] * 2, self.Mass_node_B_X[i] * 2] = 1
        #     K[(self.Mass_node_B_X[i] * 2) + 1, (self.Mass_node_B_X[i] * 2) + 1] = 1
        #     F[self.Mass_node_B_X[i] * 2] = 0
        #     if i == 0 or i == len(self.Mass_node_B_X) - 1:
        #         F[(self.Mass_node_B_X[i] * 2) + 1] = -our_vert_value / 2
        #     else:
        #         F[(self.Mass_node_B_X[i] * 2) + 1] = -our_vert_value

        rho_1 = 2500
        rho_2 = 2200
        g = 9.8
        #массовые силы
        for number_el in self.Mass_Element_C:
            el = self.Mass_Element[number_el]
            # print('rrrrrr', i, el)
            Enodes = np.matrix([self.Mass_node[el[0]], self.Mass_node[el[1]], self.Mass_node[el[2]]])
            nodes_current = self.Mass_Element[i]
            F_current = abs(np.linalg.det(self.get_J(Enodes))) / 2 * rho_1 * g / 3
            for j in range(3):
                F[2 * (nodes_current[j] - 1) + 1] += -F_current

        for number_el in self.Mass_Element_C2:
            el = self.Mass_Element[number_el]
            # print('rrrrrr', i, el)
            Enodes = np.matrix([self.Mass_node[el[0]], self.Mass_node[el[1]], self.Mass_node[el[2]]])
            nodes_current = self.Mass_Element[i]
            F_current = abs(np.linalg.det(self.get_J(Enodes))) / 2 * rho_2 * g / 3
            for j in range(3):
                F[2 * (nodes_current[j] - 1) + 1] += -F_current


        print(K)
        # print('\n'.join('\t'.join(map(str, row)) for row in matrix))
        # for i in range(len(K)):
        #     for j in range(len(K[i])):
        #         print(K[i][j], end=' ')
        # print(K)
        T = solve(K, F)
        # np.savetxt('test1.txt', T, fmt='%.7f')
        # np.savetxt('test1.txt', [T,T], fmt='%.8e')
        my_file = open("uuuu12.txt", 'w')
        X = []
        Y = []
        for i, index in enumerate(T):
            if i % 2 == 0:
                X.append((index[0]))
            else:
                Y.append(((index[0])))
        XY = np.zeros((len(X), 3))
        for i in range(len(X)):
            XY[i, 0] = X[i]
            XY[i, 1] = Y[i]
            XY[i, 2] = 0
        print('qwer', X, Y, XY)
        my_file.write('\n'.join([str(i[0]) + ' ' + str(i[1]) + ' 0' for i in XY]))
        # '\n'.join([i[1:-1] for i in ','.join(map(str,T[0])).split(",")])
        my_file.close()
        print(T)
        print(len(T))


if __name__ == '__main__':
    node = open('nodes.txt', 'r')
    Mass_node = [[float(i) for i in (line.replace(" ", '').split(",")[1:])] for line in node.read().splitlines()]
    elem_all = open('elem_nodes.txt', 'r')
    Mass_Element = [[int(i) - 1 for i in line.replace(" ", '').split(",")[1:]] for line in elem_all.read().splitlines()]
    np.savetxt('MMMMM.txt', Mass_node, fmt='%d')

    el_Priming = [int(i) - 1 for i in
                  open('Priming_el_range.txt', 'r').read().replace("\n", ',').replace(" ", '').split(",")]
    Rage_el_Priming = [i for i in range(el_Priming[0], el_Priming[1] + 1)]
    str_P = ", ".join(map(str, Rage_el_Priming))
    # print(str_P)
    np.savetxt('elem_Priming.txt', Rage_el_Priming, fmt='%i')
    print(Rage_el_Priming)
    elem_Priming = open('elem_Priming.txt', 'r')
    # int(i) - 1
    Mass_Element_Priming = [int(i) - 1 for i in elem_Priming.read()[:-1].replace("\n", ',').replace(" ", '').split(",")]
    print('aaa', Mass_Element_Priming)

    elem_Concrete = open('elem_Concrete.txt', 'r')
    Mass_Element_C = [int(i) - 1 for i in elem_Concrete.read().replace("\n", ',').replace(" ", '').split(",")]

    el_C2 = [int(i) - 1 for i in open('C2_el_range.txt', 'r').read().replace("\n", ',').replace(" ", '').split(",")]
    Rage_el_C2 = [i for i in range(el_C2[0], el_C2[1] + 1)]
    np.savetxt('elem_Concrete2.txt', Rage_el_C2, fmt='%i')
    print(Rage_el_C2)
    elem_Concrete2 = open('elem_Concrete2.txt', 'r')
    Mass_Element_C2 = [int(i) - 1 for i in elem_Concrete2.read()[:-1].replace("\n", ',').replace(" ", '').split(",")]

    Mass_node_B_X = [int(i) - 1 for i in
                     open('nodes_BC_X.txt', 'r').read().replace("\n", ',').replace(" ", '').split(",")]

    Mass_node_B_Y = [int(i) - 1 for i in
                     open('nodes_BC_Y.txt', 'r').read().replace("\n", ',').replace(" ", '').split(",")]
    E_Priming = 1.7e+10
    nu_Priming = 0.2
    E_C = 2.5e+10
    nu_C = 0.2
    E_C2 = 2.2e+10
    nu_C2 = 0.2
    # rho = 7800
    # M = 10000
    l = 0.1
    x = 1
    test = thermal_conductivity(x, l, Mass_node, Mass_Element, Mass_Element_Priming, Mass_Element_C, Mass_Element_C,
                                E_Priming,
                                nu_Priming, E_C, nu_C, E_C2, nu_C2, Mass_node_B_X, Mass_node_B_Y)
    test.Solve()
