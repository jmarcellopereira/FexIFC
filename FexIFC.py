from scipy.misc import derivative
from uncertainties import *
from uncertainties.umath import *
from uncertainties import unumpy
from numpy import *
##from math import *

def XraizN(f, x0 = 1, max_passos = 100, tolerancia = 1e-6):
    
    xn = x0
    passo = 1
    
    while(passo <= max_passos and abs(f(xn)) > tolerancia and abs(derivada(f, xn)) > tolerancia):
        xn -= f(xn)/derivada(f, xn)
        passo += 1
        
    if(abs(derivada(f, xn)) <= tolerancia):
        print("Erro! O valor da derivada da função está próximo de zero!", xn)
        
    elif (passo >= max_passos):
        print("Erro! Limite de passos estourados sem encontrar o resultado")
        
    else:
        return xn
    
def Xmod_simpson(f, a, b, passos=1000):
    h = (b - a)/float(passos)
    a1 = a + h/2.0
    
    s1 = sum(f(a1 + i*h) for i in range(0, passos))
    s2 = sum(f(a + i*h)  for i in range(1, passos))
    
    return (h/6.0)*(f(a) + f(b) + 4.0*s1 + 2.0*s2)


def XinterpL(x, x_dados, y_dados):
    pontos = tuple(zip(x_dados, y_dados))
    n = len(pontos) - 1
    
    def L(i, x):
        Li = 1
        for j in range(n + 1) :
            if(i != j):
                Li = Li*(x - pontos[j][0])/(pontos[i][0] - pontos[j][0])
        return Li

    def Lagrange(x):
        p = 0
        for i in range(n + 1):       
            p = p + pontos[i][1]*L(i, x)
        return p
            
    return Lagrange(x)

raiz         = wrap(XraizN)
derivada     = wrap(derivative)
integral     = wrap(Xmod_simpson)
interpolacao = wrap(XinterpL) 