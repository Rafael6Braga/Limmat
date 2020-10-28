# Importar as seguintes bibliotecas
import sympy as sp
import numpy as np
# Importar o método 'inv'
from numpy.linalg import inv
# Definir a letra x como parametro
x = sp.Symbol('x')
##-------- Início dos métodos auxiliares --------##
# Método para a construção da matriz dos coeficientes do sistema Ab = a
def matrizdosCoeficientes( objeto, grauDoNumerador, grauDoDenominador):
    # Matriz dos coeficientes    
    A = sp.Matrix( np.zeros( ( grauDoDenominador, grauDoDenominador ) ) ) 
    # Se o obeto é uma lista
    if (type( objeto ) == list or type( objeto ) == tuple ):
        # Se a lista contém os primeiros n + 1 coeficientes da série
        if ( grauDoNumerador + grauDoDenominador + 1 <= len(objeto) ):
            for linha in range( 0, grauDoDenominador ):  
                for coluna in range(0, grauDoDenominador ):
                    # Para todo o n menor que zero n = 0
                    if ( grauDoNumerador - grauDoDenominador + coluna + 1 + linha < 0 ):
                        A[ linha, coluna ] = 0
                    # Coeficiente a_0
                    elif ( grauDoNumerador - grauDoDenominador + coluna + 1 + linha == 0 ):
                        A[ linha, coluna ] = objeto[0]
                    # Coeficientes a_n para n maior ou igual a 1
                    else:
                        an = objeto[ grauDoNumerador - grauDoDenominador + coluna + 1 + linha ]
                        A[ linha, coluna ] =  an
            return(A)
        # erro
        else: return('O número de coeficientes não é suficiente.')
    # Se o objeto é uma função        
    else:
        for linha in range( 0, grauDoDenominador ):  
            for coluna in range( 0, grauDoDenominador ):               
                # Para todo o n menor que zero n = 0
                if ( grauDoNumerador - grauDoDenominador + coluna + linha + 1 < 0 ):
                    A[ linha, coluna ] = 0               
                # Coeficiente a_0 da série de Taylor da função
                elif ( grauDoNumerador - grauDoDenominador + coluna + 1 + linha == 0 ):
                    derivada = sp.diff( objeto, x, 0 )
                    derivadaEmZero = derivada.subs( x, 0 )                    
                    A[ linha, coluna ] = derivadaEmZero / sp.Rational(sp.factorial( 0 ) )                  
                # Coeficiente a_n da série de Taylor da função, para n maior ou igual a 1    
                else:
                    derivada = sp.diff( objeto, x, grauDoNumerador - grauDoDenominador + coluna + linha + 1  )
                    derivadaEmZero =  derivada.subs( x, 0 )                   
                    A[ linha, coluna ] = derivadaEmZero / sp.Rational(sp.factorial( grauDoNumerador - grauDoDenominador + coluna + linha + 1 ) )    
        return ( A )
    return
    
# Método para a construção da matriz dos termos independentes do sistema Ab = a
def matrizDosTermosIndependentes( objeto, grauDoNumerador, grauDoDenominador ): 
    # Matriz dos termos independentes
    a = sp.Matrix( np.zeros(( 1, grauDoDenominador )) ) 
    # Se o objeto é uma lista 
    if (type( objeto ) == list or type( objeto ) == tuple ):       
        # Se a lista contém os primeiros n + 1 coeficientes da série
        if grauDoNumerador + grauDoDenominador + 1 <= len(objeto):            
            for coluna in range( 0, grauDoDenominador ):
                an = objeto[ grauDoNumerador + coluna +1]                  
                # Coeficientes da matriz dos termos independentes
                a[ 0, coluna ] = ( -1 ) * an                
            return ( a )        
        # erro
        else: return('O número de coeficientes não é suficiente.')
    # Se o objeto é uma função 
    else:
        for coluna in range( 0, grauDoDenominador ):            
            # Cálculo dos coeficientes da série de Taylor da função, em x = 0
            derivada = sp.diff( objeto, x, grauDoNumerador + coluna + 1 )
            derivadaEmZero = derivada.subs( x, 0 )            
            # Coeficientes da matriz dos termos independentes
            a[ 0, coluna ] = (-1) *( derivadaEmZero / sp.Rational(sp.factorial( grauDoNumerador + coluna + 1 ) ) )           
        return ( a )
    return
    
# Método auxiliar para calculo das primeiras p + 1 equações
def coeficientesDaFuncao( objeto, grauDoNumerador, grauDoDenominador ):  
    # Matriz dos termos independentes
    cA = sp.Matrix(np.zeros(( 1, grauDoNumerador + 1 )))  
    # Se o objeto é uma lista
    if (type( objeto) == list or type( objeto) == tuple):      
        # Se a lista contém os primeiros n + 1 coeficientes da série
        if ( grauDoNumerador + grauDoDenominador + 1 <= len( objeto ) ):          
            for coluna in range( 0, grauDoNumerador + 1 ):
                if ( grauDoNumerador - coluna < 0 ):
                    cA[ coluna ] = 0             
                # Coeficientes da matriz dos termos independentes 
                else:
                    an = objeto[ grauDoNumerador - coluna ] 
                    cA[ coluna ] = an      
            return( cA )  
        # erro
        else: return('O número de coeficientes não é suficiente.')  
    # Se o objeto é uma função        
    else: 
        for coluna in range( 0, grauDoNumerador + 1 ):
            if ( grauDoNumerador - coluna < 0 ):
                cA[ coluna ] = 0               
            else:               
                # Cálculo dos coeficientes da série de Taylor de função, em x = 0 
                fn = sp.diff( objeto, x, grauDoNumerador - coluna )
                cA[ coluna ] = fn.subs( x, 0 )/sp.Rational( sp.factorial( grauDoNumerador - coluna ) )                
        return ( cA )
    return
##-------- Fim dos métodos auxiliares --------##

# Método direto para construção de um aproximante de Padé
def padeDiretoPrecisaoInfinita( objeto, grauDoNumerador, grauDoDenominador ):  
    # Matriz dos coeficientes do sistema Ab = a
    A = matrizdosCoeficientes( objeto, grauDoNumerador, grauDoDenominador )
    # erro
    if(A is None): return
    if (type(A) is str) : return(A)    
    # Se a matriz é nula
    if(A == 0): return("A matriz dos coeficientes não é invertível. ")        
    # Se a matriz A é invertível    
    if ( A.det() != 0 ):    
     ##-- Construção do denominador --##   
        # Cálculo da matiz inversa de A
        A1 = A.inv("LU")       
        # Cálculo dos coeficientes do denominador do aproximante de Padé
        Bn = matrizDosTermosIndependentes( objeto, grauDoNumerador, grauDoDenominador ) * A.inv("LU")        
        # Potências de x do denominador 
        Dx = sp.Matrix( np.zeros( ( grauDoDenominador + 1, 1 ) ) )
        for linha in range( 0, grauDoDenominador + 1 ):
                Dx[ linha ] = x**(linha)        
        # Vetor para os coeficientes do denominador
        bn = sp.Matrix( np.zeros(( 1, grauDoDenominador + 1 )) )        
        #  Definir b_0 = 1 para a função racional em x = 0 tomar o valor N(0) = c_0
        bn[0] = 1       
        # Matriz transposta dos coeficientes do denominador
        for coluna in range( 1, grauDoDenominador + 1 ):
            bn[ coluna ] = Bn[ -coluna ]            
        # Denominador do aproximante de Padé
        Denominador = sp.Function( 'Denominador' )
        Denominador = bn * Dx        
     ##-- Construção do numerador --##    
        # Vetor para os coeficientes do numerador
        cn = sp.Matrix( np.zeros( ( 1, grauDoNumerador + 1 ) ) )        
        # Potências do numerador
        Nx = sp.Matrix( np.zeros( ( grauDoNumerador + 1 ) ) )
        for coluna in range( 0, grauDoNumerador + 1 ):
                Nx[ coluna ] = x**( coluna )   
        # Se o grau do numerador é menor do que o grau do denominador
        if ( grauDoNumerador < grauDoDenominador ):            
            coluna = 1
            while( coluna <= grauDoNumerador + 1 ):                
                An = sp.Matrix( coeficientesDaFuncao( objeto, grauDoNumerador, grauDoDenominador )[ -coluna: ] )
                Bn = sp.Matrix( bn[ 0:coluna ] )
                cn[ coluna-1 ] = sp.Transpose( An ) * Bn                
                coluna += 1               
            # Numerador do aproximante de Padé    
            Numerador = sp.Function( 'Numerador' )
            Numerador = cn * Nx           
            # Função racional
            R = sp.Function( 'R' )
            R = Numerador / Denominador           
            # Aproximante de Padé
            Pade = sp.Function( 'Pade')
            Pade = R[0]           
        # Se o grau do mumerador é maior do que o grau do denominador
        else:           
            # Calculo dos coeficientes c_n até n = grau do denominador
            coluna = 1
            while( coluna <= grauDoDenominador + 1  ):               
                    An = sp.Matrix( coeficientesDaFuncao( objeto, grauDoNumerador, grauDoDenominador )[ -coluna: ] )
                    Bn = sp.Matrix( bn[ 0:coluna ] )
                    cn[ coluna - 1 ]=sp.Transpose( An ) * Bn                    
                    coluna += 1                  
            # Cálculo dos coeficiente c_n para n > grau do denominador        
            j = -1
            while ( coluna <= grauDoNumerador + 1):                
                    An = sp.Matrix( coeficientesDaFuncao( objeto, grauDoNumerador, grauDoDenominador )[ -coluna:j ] )
                    Bn = sp.Matrix( bn[ 0: ] )
                    cn[ coluna - 1 ] = sp.Transpose( An ) * Bn
                    coluna += 1
                    j -= 1          
            # Numerador do aproximante de Padé         
            Numerador = sp.Function( 'Numerador' )
            Numerador = cn * Nx            
            # Função racional
            R = sp.Function( 'R' )
            R = Numerador/Denominador            
            # Aproximante de Padé
            Pade = sp.Function( 'Pade')
            Pade = R[0]            
        return ( Pade )    
    # Se a matriz não é invertível
    else: return( 'A matrix não é invertível.' )   
    return
