class Gauss:

    def _parteDecendente(self, a, b):
        n = len(b)
        for k in range(n):
            aux = a[k][k]
            for j in range(k):
                a[k][j] = a[k][j]/aux
            b[k] = b[k]/aux
            for i in range(k+1, n):
                aux = a[i][k]
                for j in range(k, n):
                    a[i][j] = a[i][j] - aux * a[k][j]
                b[i] = b[i] - aux * b[k]

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



matriz = [
    [1, 2, 2],
    [2, 5, 8],
    [3, 6, 9]
]

resultado = [1, 2, 3]

gauss = Gauss()

# print(gauss.gauss(matriz, resultado))

# print(gauss.gaussConPivote(matriz, resultado))