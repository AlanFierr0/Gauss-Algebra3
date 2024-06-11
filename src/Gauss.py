import random
from datetime import time, datetime
import numpy as np

class Gauss:

    def _parteDecendente(self, a, b):
        n = len(b)
        for k in range(n):
            aux = a[k][k]
            if aux == 0:
                continue  # Skip division by zero
            # Normalize the pivot row
            for j in range(k, n):  # Note: should start from k
                a[k][j] = a[k][j] / aux
            b[k] = b[k] / aux
            # Eliminate the current column
            for i in range(k + 1, n):
                factor = a[i][k]
                for j in range(k, n):
                    a[i][j] = a[i][j] - factor * a[k][j]
                b[i] = b[i] - factor * b[k]

    def _parteDescendenteOpt(self, a, b):
        n = len(b)
        for k in range(n):
            aux = a[k][k]
            if aux == 0:
                continue  # Skip division by zero
            if k > 0:
                a[k][k - 1] = a[k][k - 1] / aux
            a[k][k] /= aux
            if k < n - 1:
                a[k][k + 1] = a[k][k + 1] / aux
            b[k] = b[k] / aux

            if k < n - 1:
                factor = a[k + 1][k]
                a[k + 1][k] = 0.0
                a[k + 1][k + 1] = a[k + 1][k + 1] - factor * a[k][k + 1]
                b[k + 1] = b[k + 1] - factor * b[k]

    def _parteAcendente(self, a, b):
        n = len(b)
        x = [0] * n  # Initialize x with zeros
        x[n - 1] = b[n - 1]
        for i in range(n - 2, -1, -1):
            suma = 0
            for j in range(i + 1, n):
                suma = suma + a[i][j] * x[j]
            x[i] = b[i] - suma
        return x

    def gauss(self, a, b, tridiagonal=False):
        a = np.array(a, dtype=float)  # Convert to numpy array for easier handling
        b = np.array(b, dtype=float)
        if tridiagonal:
            self._parteDescendenteOpt(a, b)
        else:
            self._parteDecendente(a, b)
        return self._parteAcendente(a, b)

    def gaussConPivote(self, a, b):
        maximo = a[0][0]
        index = 0
        primeraFila = a[0]
        primerResultado = b[0]
        for i in range(1, len(a)):
            if a[i][0] > maximo:
                maximo = a[i][0]
                index = i
        a[0] = a[index]
        a[index] = primeraFila
        b[0] = b[index]
        b[index] = primerResultado
        return self.gauss(a, b)

    def gaussOptTrdiagonal(self, a, b):
        self._parteDescendenteOpt(a, b)
        return self._parteAcendente(a, b)

class TestGauss:

    def __init__(self):
        self.gauss = Gauss()

    def random_tridiagonal_matrices(self, size):
        a = [[random.randint(1, 10) if abs(i-j) <= 1 else 0 for j in range(size)] for i in range(size)]
        b = [random.randint(1, 10) for _ in range(size)]
        return a, b

    def random_full_matrices(self, size):
        a = np.random.randint(10, size=(size, size))
        b = np.random.randint(10, size=(size, 1))
        return a, b

    def test_tridiagonal_matrices(self):
        for size in [100, 1000]:
            a, b = self.random_tridiagonal_matrices(size)
            result = self.gauss.gaussOptTrdiagonal(a, b)
            A = np.array(a)
            R = np.array(result)
            # Multiplicar las matrices
            resultado = A @ R

            assert np.allclose(resultado, np.array(b)), "La matriz no es igual"
            print(f"Funciona con {size}")

    def test_full_matrices(self):
        for size in [30, 50]:
            a, b = self.random_full_matrices(size)
            result = self.gauss.gauss(a, b)
            A = np.array(a)
            R = np.array(result)
            # Multiplicar las matrices
            resultado = A @ R
            assert np.allclose(resultado, np.array(b)), "La matriz no es igual"
            print(f"Funciona con {size}")

    def test_time_matrx(self):
        size = 800
        a, b = self.random_tridiagonal_matrices(size)
        initial_nonOpt_time = datetime.now()
        self.gauss.gauss(a, b)
        final_nonOpt_time = datetime.now()
        delta_non_opt = final_nonOpt_time - initial_nonOpt_time
        initial_opt_time = datetime.now()
        self.gauss.gaussOptTrdiagonal(a, b)
        final_opt_time = datetime.now()
        delta_opt_time = final_opt_time - initial_opt_time
        print(f"Funcion no optimizada tardo {delta_non_opt.total_seconds()} segundos")
        print(f"Funcion optimizada tardo {delta_opt_time.total_seconds()} segundos")

    def test_full_matrices_precision(self):
        for size in [30, 50]:
            a, b = self.random_full_matrices(size)
            result_no_pivot = self.gauss.gauss(a, b)
            result_pivot = self.gauss.gaussConPivote(a, b)
            A = np.array(a)
            b = np.array(b)
            residual_no_pivot = np.linalg.norm(A @ result_no_pivot - b)
            residual_pivot = np.linalg.norm(A @ result_pivot - b)
            print(f"Size: {size} - Residual without pivot: {residual_no_pivot}, with pivot: {residual_pivot}")

tester = TestGauss()
tester.test_tridiagonal_matrices()
tester.test_time_matrx()
tester.test_full_matrices()
tester.test_full_matrices_precision()
