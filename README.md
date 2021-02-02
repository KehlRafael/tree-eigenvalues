# tree-eigenvalues
 Implementação do algoritmo para localização de autovalores em grafos árvore publicado na revista Science Direct, disponível [aqui](https://doi.org/10.1016/j.laa.2010.08.006). O programa solicita ao usuário que ele digite o caminho para um arquivo de texto formatado e possibilita ao usuário localizar os autovalores do grafo árvore dado, além de visualizar o passo a passo da execução do algoritmo.
 
## Execução
 Basta compilar `main.c` e dar como input um arquivo de text no formato:
 ```
 n
 i j 1
 i j 1
 i j 1
 ```
 Onde `n` é o número de vértices, `i` e `j` são, respectivamente, a linha e a coluna da matriz onde o valor `1` será inserido. Como a matriz de adjacências para os grafos contemplados pelo algoritmo são sempre simétricas, o programa preenche as entradas `(i,j)` e `(j,i)` da matriz.
 
## Apresentação
 Esse trabalho foi apresentado no Salão de Iniciação Científica da UFRGS em 2017. A publicação do trabalho no repositório digital da UFRGS pode ser acessada [aqui](https://lume.ufrgs.br/handle/10183/176156). O [resumo](SIC/Resumo.pdf), [apresentação](SIC/Apresentacao.pdf) e [poster](SIC/Poster.pdf) também estão disponíveis no repositório junto com os [exemplos](Examples) utilizados na apresentação.
 
