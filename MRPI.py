#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Importar as seguintes librarias
import sympy as sp
import numpy as np
# Definir a letra x como parametro
x = sp.Symbol('x')
# Método para a construção da série f(x) truncada em n
def serieTruncada( objeto, n ):    
    # Vetor dos coeficientes
    an = sp.Matrix( np.zeros( (n + 1) ) )    
    # Vetor do parâmetro x
    fx = sp.Matrix(np.zeros( (1, n + 1) ) )    
    # Se o objeto é uma lista
    if ( type( objeto ) == list or type( objeto ) == tuple ):        
        # Se a lista contém os primeiros n + 1 coeficientes da série
        if ( len(objeto) >= n + 1 ):            
            # Construir a série f(x) truncada em n
            coluna = 0
            while( coluna <= n ):       
                an[coluna] = objeto[ coluna] 
                fx[coluna] = x**coluna
                coluna += 1             
            Sn = fx * an            
            return( Sn )        
        # erro
        else: return('O número de coeficientes não é suficiente.')            
    # Se o objeto é uma função         
    else:        
        # Construir a série de Taylor truncada em n
        coluna = 0
        while(coluna <= n ):           
            derivada = sp.diff( objeto, x, coluna )
            derivadaEmZero = derivada.subs( x, 0 ) 
            an[coluna] = derivadaEmZero / sp.Rational(sp.factorial(coluna)) 
            fx[coluna] = x**coluna
            coluna += 1         
        Sn = fx * an            
        return( Sn )   
    return

# Algoritmo de Baker para a construção de um ou vários aproximantes de Padé
def padeRecursivoPrecisaoInfinita( objeto, n, iterações, percurso ):    
    # Validar número de iterações
    max_iteraçoes = 2 * n -1
    if( iterações > max_iteraçoes ): return( 'Número de iterações inválido.' )    
    # Vetores para os numeradores e denominadores dos aproximantes de Padé 
    numerador = sp.Matrix( np.zeros( (2 + iterações) ) )
    denominador = sp.Matrix( np.zeros( (2 + iterações) ) )    
    # Definir o nome dos numerados e denominadores como funções 
    for i in range( 0, 2+iterações ):        
            Ni = sp.Function('N'+str(i))
            numerador[i] = Ni            
            Di = sp.Function('D'+str(i)) 
            denominador[i] = Di
    # Inputs para inicializar o algoritmo    
    # ( f_(n) , 1 )
    p = serieTruncada( objeto, n)
    numerador[0] = p[0]
    denominador[0] = 1    
    # ( f_(n-1) , 1 )
    q = serieTruncada( objeto, n-1 )
    numerador[1] = q[0]
    denominador[1] = 1   
    # Vetor para guardar os aproximantes de Padé construídos
    pade = sp.Matrix( np.zeros( (2 + iterações) ) )    
    #Padé [n/0]
    pade[0] = numerador[0] / denominador[0] 
    #Padé [(n-1)/0]
    pade[1] = numerador[1] / denominador[1]    
    #Validar nomalidade
    if ( str( pade[1] ) == str( pade[0] ) ): return('A tabela de Padé é não normal.')   
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
        numerador[i] = sp.expand( sp.simplify( ( c1 * numerador[i-2] - x * c0 * numerador[i-1] ) / c1 ) )
        denominador[i] = sp.expand(sp.simplify( ( c1 * denominador[i-2] - x * c0 * denominador[i-1] ) / c1 ) )        
        pade[i] = numerador[i] / denominador[i]        
        # Validar nomalidade
        if (str(pade[i]) == str(pade[i-1])): return('A tabela de Padé é não normal.')            
        i += 1        
        if ( i < j + iterações ):            
            # Coeficiente da potência de maior grau do numerador da penúltima iteração
            cc0 = sp.Poly( numerador[i-2], x )
            c0 = cc0.coeffs()[0]            
            # Coeficiente da potência de maior grau do numerador da última iteração
            cc1 = sp.Poly( numerador[i-1], x )
            c1 = cc1.coeffs()[0]            
            # Expressões recursivas de Baker para obtenção do aproximante [p-i-1/i]
            numerador[i] = sp.expand( sp.simplify( ( c1 * numerador[i-2] - c0 * numerador[i-1] ) / (c1 - c0) ) )
            denominador[i] = sp.expand( sp.simplify( ( c1 * denominador[i-2] - c0 * denominador[i-1] ) / (c1 - c0) ) )            
            pade[i] = numerador[i] / denominador[i]           
            # Validar nomalidade
            if ( str( pade[i] ) == str( pade[i-1] ) ): return( 'A tabela de Padé é não normal.' )            
        i += 1       
    # Todos os aproximantes de Padé construídos    
    if ( percurso == 1 ): return ( pade )    
    # Último aproximante de Padé construído
    elif ( percurso == 0 ): return( pade[-1] )
    #erro
    else: return( 'O parâmetro percurso toma apenas os valores: 0 ou 1.' )    
    return

