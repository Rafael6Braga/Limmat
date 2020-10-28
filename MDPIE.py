# Importar as seguintes bibliotecas
import sympy as sp
import numpy as np
# Definir a letra x como parametro
x = sp.Symbol('x')
##-------- Início dos métodos auxiliares --------##
# Método para construção do numerador
def padeExponencialNumerador( grauDoNumerador, grauDoDenominador ):    
    # definir parâmetros e objetos auxiliares
    Numerador = sp.Function('Numerador')
    N = sp.Matrix(np.zeros( ( 1, grauDoNumerador + 1 ) ))    
    # calcular os coeficientes    
    for i in range( 0, grauDoNumerador + 1 ):
        numerador = sp.factorial( grauDoNumerador + grauDoDenominador - i ) * sp.factorial( grauDoNumerador )
        denominador = sp.factorial( grauDoNumerador + grauDoDenominador ) * sp.factorial( i ) * sp.factorial( grauDoNumerador - i )        
        N[i] = (numerador / denominador ) * x**i    
    # somar os monómios -> polinómio      
    Numerador = sum( N )    
    return ( Numerador )

# Método para construção do denominador
def padeExponencialDenominador( grauDoNumerador, grauDoDenominador ):
    # definir parâmetros e objectos auxiliares    
    x = sp.Symbol('x')
    Denominador = sp.Function('Denominador')
    D = sp.Matrix(np.zeros( ( 1, grauDoDenominador + 1 ) ))    
    # calcular os coeficientes    
    for i in range( 0, grauDoDenominador + 1 ):
        numerador = sp.factorial( grauDoNumerador + grauDoDenominador - i ) * sp.factorial( grauDoDenominador )
        denominador = sp.factorial( grauDoNumerador + grauDoDenominador ) * sp.factorial(i) * sp.factorial( grauDoDenominador - i )        
        D[i] = (numerador / denominador ) *  ( (-1)**(i) ) * x**i    
    # somar os monómios -> polinómio     
    Denominador = sum(D)    
    return ( Denominador )
##-------- Fim dos métodos auxiliares --------##
# Método direto para construção de um aproximante de Padé da função exponencial
def padeExponencial( grauDoNumerador, grauDoDenominador):    
    #Aproximante de Padé
    Numerador = padeExponencialNumerador( grauDoNumerador, grauDoDenominador)
    Denominador = padeExponencialDenominador( grauDoNumerador, grauDoDenominador)
    Pade = Numerador / Denominador
    return( Pade )
