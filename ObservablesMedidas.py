import CNYT2 as complejos
import numpy as np


def probabilidadposicion(ket, position):
    n = complejos.Modulo(ket[position]) ** 2
    d = complejos.normaVector(ket) ** 2
    return round(n/d, 3)


def bra(ket):
    for num in ket:
        if isinstance(num, list):
            num[1] *= -1
        else:
            num *= -1
    return ket


def transicion(ket1, ket2):
    braket2 = bra(complejos.traspuesta(ket2)[0])
    norma1 = complejos.normaVector(complejos.traspuesta(ket1)[0])
    norma2 = complejos.normaVector(complejos.traspuesta(ket2)[0])
    norma = norma1 * norma2
    aux = complejos.traspuesta(ket1)
    prob = complejos.productoMatrices(braket2, complejos.traspuesta(ket1)[0])
    ans = main.Producto([1/norm, 0], prob)
    return ans


def media(observable, ket):
    bra_ket = bra(ket)
    r1 = complejos.accion(observable, ket)
    r2 = complejos.productoMatrices(res1, bra_ket)
    return r2


def varianza(observable, ket):
    braket = bra(ket)
    miu = media(observable, ket)
    ident_miu = [[(0, 0) for j in range(len(observable[0]))] for i in range(len(observable))]
    for i in range(len(observable)):
        for j in range(len(observable[i])):
            if i == j:
                ident_miu[i][j] = complejos.inversaMatriz(miu)
    ident_miu = complejos.AdicionMatricesComplejas(ident_miu, observable)
    s = complejos.productoMatrices(ident_miu, ident_miu)
    a1 = complejos.accion(s, ket)
    a2 = complejos.productoMatrices(a1, braket)
    return a2


def eigenValuVect(matrix):
    evalues, evectors = np.linalg.eig(matrix)
    values = []
    vectors = []
    for i in range(len(evalues)):
        values += [(evalues[i].real, evalues[i].imag)]
    for i in range(len(evectors)):
        vector = []
        for j in range(len(evectors[0])):
            vector += [(evectors[i][j].real, evectors[i][j].imag)]
        vectors += [vector]
    return values, vectors
def unitaria(matrix):
    if len(matrix) == len(matrix[0]):
        id = [[[0, 0] for j in range(len(matrix[0]))] for i in range(len(matrix))]
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i == j:
                    id[i][j] = [1, 0]
        aux = complejos.matrizAdjunta(matrix)
        p = complejos.productoMatrices(aux, matrix)
        ans = True
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if p[i][j] != id[i][j]:
                    ans = False
        if ans:
            return True
        else:
            return False
    else:
        return False


# Ejercicios
# 4.3.1
def ex1():
    vector = [[[1, 0]], [[0, 0]]]
    observable = [[[0, 0], [1, 0]], [[1, 0], [0, 0]]]
    observacion = complejos.accion(observable, vector)
    values, vectors = eigenValuVect(observable)
    print('Resultado de la observacion', observacion)
    print('EigenValues: ', values)
    print('EigenVectors: ', vectors)


# 4.3.2
def ex2():
    vector = [[[1, 0]], [[0, 0]]]
    observable = [[[0, 0], [1, 0]], [[1, 0], [0, 0]]]
    values, vectors = eigenValuVect(observable)
    for i in range(len(vectors)):
        print(transicion(vector, vectors[i]))


# 4.4.1
def ex3():
    matrix1 = [[[0, 0], [1, 0]], [[1, 0], [0, 0]]]
    matrix2 = [[[(2**(1/2))/2, 0], [(2**(1/2))/2, 0]], [[(2**(1/2))/2, 0], [-(2**(1/2))/2, 0]]]
    if unitaria(matrix1) and unitaria(matrix2):
        print(complejos.unitaria(complejos.productoMatrices(matrix1, matrix2)))

