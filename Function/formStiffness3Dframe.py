import numpy as np
import math as m


def formStiffness3Dframe(GDof, numberElements, elementNodes, numberNodes, nodeCoordinates, E, A, Iz, Iy, G, J):
    stiffness = np.zeros((GDof, GDof))
    for i in range(0, numberElements):
        indice = elementNodes[i]
        elementDof = [6*indice[0, 0]-6, 6*indice[0, 0]-5, 6*indice[0, 0]-4,
                      6*indice[0, 0]-3, 6*indice[0, 0]-2, 6*indice[0, 0]-1, 6*indice[0, 1]-6, 6*indice[0, 1]-5, 6*indice[0, 1]-4, 6*indice[0, 1]-3, 6*indice[0, 1]-2, 6*indice[0, 1]-1]
        x1 = nodeCoordinates[indice[0, 0]-1, 0]
        y1 = nodeCoordinates[indice[0, 0]-1, 1]
        z1 = nodeCoordinates[indice[0, 0]-1, 2]
        x2 = nodeCoordinates[indice[0, 1]-1, 0]
        y2 = nodeCoordinates[indice[0, 1]-1, 1]
        z2 = nodeCoordinates[indice[0, 1]-1, 2]
        L = m.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1))
        k1 = E*A/L
        k2 = 12*E*Iz/(L*L*L)
        k3 = 6*E*Iz/(L*L)
        k4 = 4*E*Iz/L
        k5 = 2*E*Iz/L
        k6 = 12*E*Iy/(L*L*L)
        k7 = 6*E*Iy/(L*L)
        k8 = 4*E*Iy/L
        k9 = 2*E*Iy/L
        k10 = G*J/L

        a = np.matrix([[k1, 0, 0], [0, k2, 0], [0, 0, k6]])
        b = np.matrix([[0, 0, 0], [0, 0, k3], [0, -k7, 0]])
        c = np.matrix([[k10, 0, 0], [0, k8, 0], [0, 0, k4]])
        d = np.matrix([[-k10, 0, 0], [0, k9, 0], [0, 0, k5]])

        k = np.vstack((np.hstack((a, b, -a, b)), np.hstack((b.T, c, b, d)),
                      np.hstack(((-a).T, b.T, a, -b)), np.hstack((b.T, d.T, (-b).T, c))))
        if x1 == x2 & y1 == y2:
            if z2 > z1:
                Lambda = np.matrix([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
            else:
                Lambda = np.matrix([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        else:
            CXx = (x2-x1)/L
            CYx = (y2-y1)/L
            CZx = (z2-z1)/L
            D = m.sqrt(CXx*CXx + CYx*CYx)
            CXy = -CYx/D
            CYy = CXx/D
            CZy = 0
            CXz = -CXx*CZx/D
            CYz = -CYx*CZx/D
            CZz = D
            Lambda = np.matrix(
                [[CXx, CYx, CZx], [CXy, CYy, CZy], [CXz, CYz, CZz]])

        R = np.vstack((np.hstack((Lambda, np.zeros((3, 9)))), np.hstack((np.zeros((3, 3)), Lambda, np.zeros(
            (3, 6)))), np.hstack((np.zeros((3, 6)), Lambda, np.zeros((3, 3)))), np.hstack((np.zeros((3, 9)), Lambda))))
        idx = np.ix_(elementDof, elementDof)
        stiffness[idx] = stiffness[idx]+R.T*k*R
    return stiffness, R,k,elementDof