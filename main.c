#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <dos.h>
#include <dir.h>
#include <math.h>

#define EPS 0.00000001

void Diagonalize(double **adj, double alpha, int dim);
double **OperaNotNull(double **adj, int i, int j, int dim);
double **OperaNull(double **adj, int i, int j, int k, int dim, double alpha);
void MostraMatriz(double **adj, int dim, int k, int l);
void CompactaMatriz(double **adj, int compactMatrix[][2], int dim);
double *DiagonalizeCompact(int compactMatrix[][2], double *D, int nnull);
void NoIntervalo(double **adj, int compactMatrix[][2], double *D, int nnull, int dim, double a, double b);
void SetColor(int ForgC);

int main(int argc, char* argv[]){
    FILE *file=NULL;
    int i, j=1, dim, val, nnull=0, k=0, maior=0, menor=0, igual=0;
    double **adj;
    double alpha, a=0, b=0;
    char arq[16]="";

    SetColor(10);
    while(file==NULL){
        printf ("Digite o nome completo do arquivo com os dados da matriz: \n");
        SetColor(11);
        scanf ("%s", arq);
        SetColor(10);
        printf ("\n");
        file = fopen(arq, "r");

        if(file==NULL){
            SetColor(12);
            printf("Cannot open file.\n");
            SetColor(10);
        }
    }

    fscanf(file,"%i \n",&dim);

    double D[dim];
    adj = malloc(dim*sizeof(double));

    for (i=0; i<dim; i++)
        adj[i] = malloc(dim*sizeof(double));

    //printf("%i \n",dim);

    for(int i=0; i<dim; i++){
        for(int j=i; j<dim; j++){
            adj[i][j]=0;
            adj[j][i]=0;
        }
    }

    while ((fscanf(file,"%i %i %i\n",&i, &j, &val))!=EOF){
        //printf("%i %i %i\n",i, j, val);
        adj[i-1][j-1]=val;
        adj[j-1][i-1]=val;
        if(val==1){
            nnull++;
        }
    }

    fclose(file);

    nnull=2*nnull;

    int compactMatrix[nnull][2];

    CompactaMatriz(adj, compactMatrix, dim);

    printf ("Executar passo a passo na matriz? (0-Nao, 1-Sim) \n");
    SetColor(11);
    scanf ("%i", &k);
    printf("\n");
    SetColor(10);
    if(k==1){
        printf ("Digite o valor de alpha: \n");
        SetColor(11);
        scanf ("%lf", &alpha);
        printf("\n");
        SetColor(10);

        Diagonalize(adj, alpha, dim);

    } else {
        printf ("Deseja saber quantos autovalores ha em um intervalo? (0-Nao, 1-Sim) \n");
        SetColor(11);
        scanf ("%i", &k);
        printf("\n");
        SetColor(10);

        if(k==1){
            printf ("Digite o intervalo, separando o inicio e fim por um espaco: \n");
            SetColor(11);
            scanf ("%lf %lf", &a, &b);
            printf("\n");
            SetColor(10);

            SetColor(14);
            NoIntervalo(adj, compactMatrix, D, nnull, dim, a, b);
            SetColor(10);
        } else {
            printf ("Digite o valor de alpha: \n");
            SetColor(11);
            scanf ("%lf", &alpha);
            printf("\n");
            SetColor(10);

            for(i = 0; i < dim; i++){
                D[i]=-alpha;
            }

            DiagonalizeCompact(compactMatrix, D, nnull);

            for(i = 0; i < dim; i++){
                if(fabs(D[i]) < EPS){
                    igual++;
                } else {
                    if(D[i] > 0){
                        maior++;
                    } else {
                        menor++;
                    }
                }
            }

            printf ("A arvore dada possui: \n");
            SetColor(14);
            printf ("%i autovalor(es) maior(es) que %f \n%i autovalor(es) menor(es) que %f \n",
                    maior, alpha, menor, alpha);
            if(igual != 0){
                printf ("%f e autovalor com multiplicidade %i \n", alpha, igual);
            }
            SetColor(10);
        }
    }

    for (i=0; i<dim; i++){
        free(adj[i]);
    }
    free(adj);
    for (i=0; i<2; i++){
        free(compactMatrix[i]);
    }
    free(compactMatrix);
    //free(D);

    printf("\n");
    system("pause");

    //return 0;
}

void Diagonalize(double **adj, double alpha, int dim){
    int i=0, j=0, k=0, l=0, maior=0, menor=0, igual=0;

    for(i = 0; i < dim; i++){
        adj[i][i] = -alpha;
    }

    MostraMatriz(adj, dim, -1, -1);
    printf("Matriz inicial, com -alpha na diagonal. \n");
    system("pause");
    printf("\n");

    for(i = dim-1; i >= 0; i--){
        for(j = dim-1; j > i; j--){
            if(adj[i][j] == 1){
                if(fabs(adj[j][j]) < EPS){
                    if(i != 0){
                        for(l = 0; l < i; l++){
                            if(adj[i][l]==1){
                                k = l;
                                l = i;
                            }
                        }
                    } else {
                        k = -1;
                    }

                    OperaNull(adj, i, j, k, dim, alpha);
                    k = 0;

                } else {
                    OperaNotNull(adj, i, j, dim);
                }
            }
        }
    }

    MostraMatriz(adj, dim, -1, -1);
    printf("A matriz diagonalizada. \n");

    for(i = 0; i < dim; i++){
        if(fabs(adj[i][i]) < EPS){
            igual++;
        } else {
            if(adj[i][i] > 0){
                maior++;
            } else {
                menor++;
            }
        }
    }

    printf ("A arvore dada possui: \n");
    SetColor(14);
    printf ("%i autovalor(es) maior(es) que %f \n%i autovalor(es) menor(es) que %f \n",
                    maior, alpha, menor, alpha);
    if(igual != 0){
        printf ("%f e autovalor com multiplicidade %i \n", alpha, igual);
    }
    SetColor(10);

}

double **OperaNotNull(double **adj, int i, int j, int dim){
    int l;

    for(l = 0; l < dim; l++){
        adj[i][l] = adj[i][l]-(1/adj[j][j])*adj[j][l];
        if(fabs(adj[i][l]) < EPS){
            adj[i][l] = fabs(adj[i][l]);
        }
    }
    for(l = 0; l < dim; l++){
        adj[l][i] = adj[l][i]-(1/adj[j][j])*adj[l][j];
        if(fabs(adj[l][i]) < EPS){
            adj[l][i] = fabs(adj[l][i]);
        }
    }
    MostraMatriz(adj, dim, i, j);
    printf("A linha  %i recebeu a l(%i)-(1/%f)*l(%i) \n", i, i, adj[j][j], j);
    printf("A coluna %i recebeu a c(%i)-(1/%f)*c(%i) \n", i, i, adj[j][j], j);
    system("pause");
    printf("\n");

    return adj;
}

double **OperaNull(double **adj, int i, int j, int k, int dim, double alpha){
    int l, m;

    for(l = j-1; l > i; l--){ // elimina as entradas nas linhas e colunas dos outros filhos
        if(adj[i][l] == 1){
            for(m = 0; m < dim; m++){ // linha l - linha j, coluna l - coluna j
                adj[l][m] = adj[l][m]-adj[j][m];
                adj[m][l] = adj[m][l]-adj[m][j];

                if(fabs(adj[l][m]) < EPS){
                    adj[l][m] = fabs(adj[l][m]);
                }
                if(fabs(adj[m][l]) < EPS){
                    adj[m][l] = fabs(adj[m][l]);
                }
            }
            MostraMatriz(adj, dim, l, j);
            printf("A linha  %i recebeu a l(%i)-l(%i) \n", l, l, j);
            printf("A coluna %i recebeu a c(%i)-c(%i) \n", l, l, j);
            system("pause");
            printf("\n");
        }
    }

    /****************************************************/
    if(k >= 0){
        for(l = 0; l < dim; l++){ // linha k - linha j, coluna k - coluna j
            adj[k][l] = adj[k][l]-adj[j][l];
            adj[l][k] = adj[l][k]-adj[l][j];

            if(fabs(adj[k][l]) < EPS){
                adj[k][l] = fabs(adj[k][l]);
            }
            if(fabs(adj[l][k]) < EPS){
                adj[l][k] = fabs(adj[l][k]);
            }
        }
        MostraMatriz(adj, dim, k, j);
        printf("A linha  %i recebeu a l(%i)-l(%i) \n", k, k, j);
        printf("A coluna %i recebeu a c(%i)-c(%i) \n", k, k, j);
        system("pause");
        printf("\n");
    }

    /****************************************************/

    if(fabs(alpha) > EPS){
        for(l = 0; l < dim; l++){ // linha i + alpha/2 linha j, coluna i + alpha/2 coluna j
            adj[i][l] = adj[i][l]+(alpha/2)*adj[j][l];
            adj[l][i] = adj[l][i]+(alpha/2)*adj[l][j];

            if(fabs(adj[i][l]) < EPS){
                adj[i][l] = fabs(adj[i][l]);
            }
            if(fabs(adj[l][i]) < EPS){
                adj[l][i] = fabs(adj[l][i]);
            }
        }
        MostraMatriz(adj, dim, k, j);
        printf("A linha  %i recebeu a l(%i)+(%f/2)*l(%i) \n", i, i, alpha, j);
        printf("A coluna %i recebeu a c(%i)+(%f/2)*c(%i) \n", k, k, alpha, j);
        system("pause");
        printf("\n");
    }

    /****************************************************/

    for(l = 0; l < dim; l++){ // linha j + linha i, coluna j + coluna i
        adj[j][l] = adj[j][l]+adj[i][l];
        adj[l][j] = adj[l][j]+adj[l][i];

        if(fabs(adj[j][l]) < EPS){
            adj[j][l] = fabs(adj[j][l]);
        }
        if(fabs(adj[l][j]) < EPS){
            adj[l][j] = fabs(adj[l][j]);
        }
    }
    MostraMatriz(adj, dim, j, i);
    printf("A linha  %i recebeu a l(%i)+l(%i) \n", j, j, i);
    printf("A coluna %i recebeu a c(%i)+c(%i) \n", j, j, i);
    system("pause");
    printf("\n");

    /****************************************************/

    for(l = 0; l < dim; l++){ // linha i - linha j/2
        adj[i][l] = adj[i][l]-(adj[j][l]/2);

        if(fabs(adj[i][l]) < EPS){
            adj[i][l] = fabs(adj[i][l]);
        }
    }
    for(l = 0; l < dim; l++){ // coluna i - coluna j/2
        adj[l][i] = adj[l][i]-(adj[l][j]/2);

        if(fabs(adj[l][i]) < EPS){
            adj[l][i] = fabs(adj[l][i]);
        }
    }
    MostraMatriz(adj, dim, i, j);
    printf("A linha  %i recebeu a l(%i)-(1/2)*l(%i) \n", i, i, j);
    printf("A coluna %i recebeu a c(%i)-(1/2)*c(%i) \n", i, i, j);
    system("pause");
    printf("\n");

    return adj;
}

void MostraMatriz(double **adj, int dim, int k, int l){
    SetColor(15);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if((j == l || i == k) ||
               (j == k || i == l)){
                SetColor(8);
            } else {
                SetColor(15);
            }
            if (adj[i][j] < 0){
                printf("%f\t", adj[i][j]);
            } else {
                printf("%f\t", adj[i][j]);
            }
        }
        printf("\n");
    }
    SetColor(10);
}

void CompactaMatriz(double **adj, int compactMatrix[][2], int dim){
    int i,j,k=0;
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++){
            if (adj[i][j] != 0){
                compactMatrix[k][0] = i;
                compactMatrix[k][1] = j;
                k++;
            }
        }
    }
}

double *DiagonalizeCompact(int compactMatrix[][2], double *D, int nnull){
    int i=0, k=0, aux=0;
    double sum=0;

    for(i=nnull-1; i>=0; i--){
        if(compactMatrix[i][0] < compactMatrix[i][1] && i > 0){
            if(fabs(D[compactMatrix[i][1]]) < EPS){
                D[compactMatrix[i][1]]=2;
                D[compactMatrix[i][0]]=-0.5;

                aux=compactMatrix[i][0];
                while(compactMatrix[i][0]==aux){
                    i--;
                    if(compactMatrix[i][0] > compactMatrix[i][1]){
                        for(k=0; k < i; k++){
                            if(compactMatrix[i][0] == compactMatrix[k][1] &&
                               compactMatrix[i][1] == compactMatrix[k][0]){

                                compactMatrix[k][1] = -1;
                                compactMatrix[i][1] = -1;
                                sum = 0;
                                k = i;
                            }
                        }
                    }
                }
                i++;
            } else {
                sum = sum + 1/D[compactMatrix[i][1]];
            }
        } else {
            if(i==0){
                if(compactMatrix[0][1]>=0){
                    if(fabs(D[compactMatrix[0][1]]) < EPS){
                        D[compactMatrix[0][1]]=2;
                        D[compactMatrix[0][0]]=-0.5;
                    } else {
                        sum = sum + 1/D[compactMatrix[0][1]];
                        D[compactMatrix[0][0]] = D[compactMatrix[0][0]] - sum;
                    }
                } else {
                    D[compactMatrix[0][0]] = D[compactMatrix[0][0]] - sum;
                }
                sum = 0;
            } else {
                if(compactMatrix[i][1] >= 0 ||
                   compactMatrix[i][0] != compactMatrix[i-1][0]){
                    D[compactMatrix[i][0]] = D[compactMatrix[i][0]] - sum;
                    sum = 0;
                }
            }
        }
    }

    return D;
}

void NoIntervalo(double **adj, int compactMatrix[][2], double *D, int nnull, int dim, double a, double b){
    int maior1=0, menor1=0, igual1=0, maior2=0, menor2=0, igual2=0, i;

    /*printf("%f\t%f\n", a,b);*/

    CompactaMatriz(adj, compactMatrix, dim);
    for(i = 0; i < dim; i++){
        D[i]=-a;
    }
    DiagonalizeCompact(compactMatrix, D, nnull);

    for(i = 0; i < dim; i++){
        if(fabs(D[i]) < EPS){
            igual1++;
        } else {
            if(D[i] > 0){
                maior1++;
            } else {
                menor1++;
            }
        }
    }

    CompactaMatriz(adj, compactMatrix, dim);
    for(i = 0; i < dim; i++){
        D[i]=-b;
    }
    DiagonalizeCompact(compactMatrix, D, nnull);

    for(i = 0; i < dim; i++){
        if(fabs(D[i]) < EPS){
            igual2++;
        } else {
            if(D[i] > 0){
                maior2++;
            } else {
                menor2++;
            }
        }
    }

    printf("No intervalo [%f,%f] existe(m) %i autovalor(es) \n", a,b,maior1-maior2+igual1);
}

void SetColor(int ForgC){
    WORD wColor;

    /**
    Name         | Value
                 |
    Black        |   0
    Blue         |   1
    Green        |   2
    Cyan         |   3
    Red          |   4
    Magenta      |   5
    Brown        |   6
    Light Gray   |   7
    Dark Gray    |   8
    Light Blue   |   9
    Light Green  |   10
    Light Cyan   |   11
    Light Red    |   12
    Light Magenta|   13
    Yellow       |   14
    White        |   15
    **/

    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO csbi;

    //We use csbi for the wAttributes word.
    if(GetConsoleScreenBufferInfo(hStdOut, &csbi)){
        //Mask out all but the background attribute, and add in the forgournd color
        wColor = (csbi.wAttributes & 0xF0) + (ForgC & 0x0F);
        SetConsoleTextAttribute(hStdOut, wColor);
    }
    return;
}
