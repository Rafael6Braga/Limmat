# Importar as seguintes bibliotecas
import sympy as sp
import numpy as np
from mpmath import *
# Definir a letra x como parametro
x = sp.Symbol('x')
##-------- Início do método auxiliar --------##
# Método para a construção da série f(x) truncada em n
def serieTruncada( objeto, n, prec ):    
    # Fixar precisão aritmética de 'prec' algarismos significativos
    mp.dps = prec    
    # Vetor dos coeficientes
    an = sp.Matrix( np.zeros( (n + 1) ) )    
    # Vetor do parâmetro x
    fx = sp.Matrix(np.zeros( (1, n + 1) ) )    
    # Se o objeto é uma lista
    if ( type( objeto ) == list or type( objeto ) == tuple ):        
        # Se a lista contém os primeiros n + 1 coeficientes da série
        if ( len( objeto ) >= n + 1):            
            # Construir a série f(x) truncada em n
            coluna = 0
            iteracoes = 0
            while ( coluna <= n ):                   
            # valor o coeficiente com 'prec' algarismos significativos
                aN =  objeto[ coluna] 
                if(aN.is_integer == False):
                        st_aN = str(aN)
                        fraq = 0
                        while( st_aN[fraq] != '/'):
                                fraq += 1
                        an[coluna] = an_float = mpf(st_aN[0:fraq])/ mpf(st_aN[fraq+1:])
                else: 
                    an[coluna] = mpf(aN)                    
                fx[coluna] = x**coluna
                coluna += 1                
            Sn = fx * an           
            return( Sn )      
        # erro
        else: return('Número de coeficientes insuficiente.')        
    # Se o objeto é uma função         
    else:        
        # Construir a série de Taylor truncada em n
        coluna = 0
        while(coluna <= n ):              
            derivada = sp.diff( objeto, x, coluna )
            derivadaEmZero = derivada.subs( x, mpf(0) ) 
            an[coluna] = mpf( sp.N( derivadaEmZero, prec )) / mpf(sp.factorial(coluna) ) 
            fx[coluna] = x**coluna
            coluna += 1           
        Sn = fx * an
        return( Sn )    
    return
##-------- Fim do método auxiliar --------##
# Algoritmo de Baker para a construção de um ou vários aproximantes de Padé
def padeRecursivoPrecisaoFinita( objeto, n, iterações, percurso, prec ):    
    # Fixar precisão aritmética de 'prec' algarismos significativos
    mp.dps = prec    
    # Validar número de iterações
    max_iteraçoes = 2 * n -1
    if( iterações > max_iteraçoes ): return( 'Número de iterações inválido.' )    
   # Vetores para os numeradores e denominadores dos aproximantes de Padé 
    numerador = sp.Matrix( np.zeros( (2 + iterações) ) )
    denominador = sp.Matrix(np.zeros( (2 + iterações) ) )    
    # Inputs para inicializar o algoritmo    
    # ( f_(n) , 1 )    
    f_n = serieTruncada( objeto, n, prec )
    numerador[0] = f_n[0]
    denominador[0] = 1    
    # ( f_(n-1) , 1 )
    f_n_1 = serieTruncada( objeto, n - 1, prec )
    numerador[1] = f_n_1[0]
    denominador[1] = 1   
    # Vetor para guardar os aproximantes de Padé construídos
    pade = sp.Matrix( np.zeros( (2 + iterações) ) )    
    #Padé [n/0]
    pade[0] = numerador[0] / denominador[0] 
    #Padé [(n-1)/0]
    pade[1] = numerador[1] / denominador[1]        
    # Iniciar o algoritmo de Baker
    i = 2
    j = 2
    while ( i < j + iterações ):        
        # Coeficiente da potência de maior grau do numerador da penúltima iteração
        cc0 = sp.Poly( numerador[i-2], x )
        c0 = cc0.coeffs()[0]        
         # Coeficiente da potência de maior grau do numerador da última iteração     
        cc1 = sp.Poly( numerador[i-1], x )
        c1 = cc1.coeffs()[0]        
        # Expressões recursivas de Baker para obtenção do aproximante [p-i/i]
        numerador[i] = (1/c1) * sp.expand( sp.simplify( ((  (c1) * numerador[i-2] ) - (c0) * x  * numerador[i-1] )))
        denominador[i] = (1/c1) * sp.expand( sp.simplify( (  (c1) * denominador[i-2] - (c0)* x  * denominador[i-1] )))        
        pade[i] = numerador[i] / denominador[i]                    
        i += 1        
        if ( i < j + iterações ):            
            # Coeficiente da potência de maior grau do numerador da penúltima iteração
            cc0 = sp.Poly( numerador[i-2], x )
            c0 = cc0.coeffs()[0]
            # Coeficiente da potência de maior grau do numerador da última iteração  
            cc1 = sp.Poly( numerador[i-1], x )
            c1 = cc1.coeffs()[0]
            c2 = c1 - c0
            # Expressões recursivas de Baker para obtenção do aproximante [p-i-1/i]
            numerador[i] = (1/c2)*sp.expand( sp.simplify( ( ((c1) * numerador[i-2]) - (c0) * numerador[i-1] )  ) )
            denominador[i] = (1/c2)*sp.expand( sp.simplify( ( (c1) * denominador[i-2] - (c0) * denominador[i-1] )  ) )             
            pade[i] = numerador[i] / denominador[i]                        
        i += 1        
    # Todos os aproximantes de Padé construídos   
    if ( percurso == 1 ): return ( pade )    
    # Último aproximante de Padé construído
    elif ( percurso == 0 ): return(pade[-1])
    #erro
    else: return('O parâmetro percurso toma apenas os valores: 0 ou 1.')
    return
