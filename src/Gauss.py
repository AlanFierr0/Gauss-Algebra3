import random

class Gauss:

    def _parteDecendente(self, a, b):
        n = len(b)
        for k in range(n):
            aux = a[k][k]
            if aux == 0:
                continue  # Skip division by zero
            for j in range(k):
                a[k][j] = a[k][j] / aux
            b[k] = b[k] / aux
            for i in range(k + 1, n):
                aux = a[i][k]
                for j in range(k, n):
                    a[i][j] = a[i][j] - aux * a[k][j]
                b[i] = b[i] - aux * b[k]

    def _parteDescendenteOpt(self, a, b):
        n = len(b)
        for k in range(n):
            aux = a[k][k]
            if aux == 0:
                continue  # Skip division by zero
            if k > 0:
                a[k][k - 1] = a[k][k - 1] / aux
            a[k][k] /= a[k][k]
            if k < n - 1:
                a[k][k + 1] = a[k][k + 1] / aux
            b[k] = b[k] / aux

            if k < n - 1:
                aux = a[k + 1][k]
                a[k + 1][k] = 0.0
                a[k + 1][k + 1] = a[k + 1][k + 1] - aux * a[k][k + 1]
                b[k + 1] = b[k + 1] - aux * b[k]

    def _parteAcendente(self, a, b):
        n = len(b)
        x = list(range(n))
        x[n - 1] = b[n - 1]
        for i in range(n - 2, -1, -1):
            suma = 0
            for j in range(i + 1, n):
                suma = suma + a[i][j] * x[j]
            x[i] = b[i] - suma
        return x

    def gauss(self, a, b, opt=False):
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

        scale_factor = 10
        a = [[random.randint(1, 10) // scale_factor for _ in range(size)] for _ in range(size)]
        b = [random.randint(1, 10) // scale_factor for _ in range(size)]
        return a, b

    def test_tridiagonal_matrices(self):

        for size in [100, 1000]:
            a, b = self.random_tridiagonal_matrices(size)
            print("Matriz original (a):")
            for row in a:
                print(row)
            print("Vector b original:")
            print(b)

            result = self.gauss.gaussOptTrdiagonal(a, b)
            print("Resultado después de aplicar la eliminación gaussiana:")
            print(result)


            assert len(result) == len(b), "La longitud de la solución no coincide con el vector b"

    def test_full_matrices(self):

        for size in [30, 50]:
            a, b = self.random_full_matrices(size)
            print("Matriz original (a):")
            for row in a:
                print(row)
            print("Vector b original:")
            print(b)

            result = self.gauss.gauss(a, b)
            print("Resultado después de aplicar la eliminación gaussiana:")
            print(result)
            assert len(result) == len(b), "La longitud de la solución no coincide con el vector b"



tester = TestGauss()
tester.test_tridiagonal_matrices()
tester.test_full_matrices()
