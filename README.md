## Portuguese version

### Rotinas para a construção de um ou vários aproximantes de Padé.

##### Estas rotinas foram desenvolvidas para a elaboração da tese para obtenção do grau de Mestre em Engenharia Matemática na Faculdade de Ciências da Universidade do Porto. Em breve será colocado um link que remete para o documento em formato pdf.

#### A1. MDPI - Método Direto Precisão Infinita

##### Para a construção de um aproximante de Padé, usando o método direto com precisão aritmética infinita, basta executar todos os métodos definidos em MDPI.py, e por último executar o comando:
                          padeDiretoPrecisaoInfinita(objeto, grauDoNumerador, grauDoDenominador) 
##### onde:
    objeto = uma função ou uma lista de coeficientes,
    grauDoNumerador = grau do numerador do aproximante de Padé [p/q](x),
    grauDoDenominador = grau do denominador do aproximante de Padé [p/q](x).

##### No caso de escolher uma função, terá de utilizar o prefixo sp. (sympy) seguido do nome da função disponível na biblioteca sympy, por exemplo,

                                      objeto = sp.exp(x), 
				      
##### ou criar uma função, por exemplo,

                                      objeto = (x+1)/sp.sqrt(x**2+4). 
                                                       
##### O parâmetro por defeito é a letra x. Se pretender mudar para outra letra deve alterar todas as ocorrências da variável x em todos os métodos. A rotina irá construir a série de Taylor centrada em 0 dessa função e utilizar os primeiros p+q+1 coeficientes para a construção do aproximante de Padé \[p/q\](x). Alternativamente, pode definir uma lista com os primeiros p+q+1 coeficientes da série formal de potências, para a construção do aproximante de Padé \[p/q\]f(x). Nessa lista deverá colocar os valores, por exemplo, da seguinte forma:

                    objeto = [sp.Rational(1)/sp.Rational(234), sp.Rational(8)/sp.Rational(3), ... ]
                                      
##### garantido a precisão infinita dos coeficientes.


#### A2. MDPIE - Método Direto Precisão Infinita Exponencial
##### Para a construção de um aproximante de Padé, usando o método direto com precisão aritmética finita, basta executar todos os métodos definidos em MDPIE.py, e por último executar o comando:

                       padeDiretoPrecisaoFinita(objeto, grauDoNumerador, grauDoDenominador, prec)
		       
##### onde:

    objeto = uma função ou uma lista de coeficientes (ver comentários em A.1),
    grauDoNumerador = grau do numerador do aproximante de Padé [p/q](x),
    grauDoDenominador = grau do denominador do aproximante de Padé [p/q](x),
    prec = número de algarismos significativos.
    

#### A3. MDPI - Método Recursivo Precisão Infinita

##### Para a construção de um ou vários aproximante de Padé usando o algoritmo de Baker com precisão aritmética infinita, basta executar todos o métodos definidos em MDPI.py, e por último executar o comando:

                       padeRecursivoPrecisaoInfinita(objeto, n, iterações, percurso)
		       
##### onde: 
     objeto = uma função ou uma lista de coeficientes (ver comentários em A.1),
     n = série formal de potências f(x) truncada em x^n, 
     iterações = no máximo 2*n-1, ou seja, se por exemplo n=3, então o algoritmo começa em [3/0] e [2/0], ao fim de 
    2*-1 iterações o último aproximante de Padé construído será o [0/3],
     percurso = 1 para visualizar todos os aproximante de Padé construídos até à iteração k; 0 para visualizar o aproximante de 
     Padé construído na iteração k.
     

#### A4. MDPI - Método Direto Precisão Finita

##### Para a construção de um aproximante de Padé, usando o método direto com precisão aritmética finita, basta executar todos os métodos definidos em MDPI.py, e por último executar o comando:
                        padeDiretoPrecisaoFinita(objeto, grauDoNumerador, grauDoDenominador, prec) 
##### onde:

      objeto = uma função ou uma lista de coeficientes (ver comentários em A.1),
      grauDoNumerador = grau do numerador do aproximante de Padé [p/q](x),
      grauDoDenominador = grau do denominador do aproximante de Padé [p/q](x),
      prec = número de algarismos significativos.
      

#### A5. MDPI - Método Recursivo Precisão Finita

##### Para a construção de um ou vários aproximante de Padé usando o algoritmo de Baker com precisão aritmética finita, basta executar todos os métodos definidos em MDPI.py, e por último executar o comando:
                         padeRecursivoPrecisaoInfinita(objeto, n, iterações, percurso) 
##### onde:
     objeto = uma função ou uma lista de coeficientes (ver comentários em A.1),
     n = série formal de potências f(x) truncada em x^n,
     iterações = no máximo 2n-1, ou seja, se por exemplo n=3, então o algoritmo começa em [3/0] e [2/0], ao fim de 
    2*3-1 iterações o último aproximante de Padé construído será o [0/3],
     percurso = 1 para visualizar todos os aproximante de Padé construídos até à iteração k; 0 para visualizar o aproximante de
    Padé construído na iteração k,
     prec = número de algarismos significativos.

