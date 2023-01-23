// klou Anas and Cachin--Bernard Gwendal
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct{
    double *values;
    uint8_t dim;
}matrix;

/*
 * This function return the value of matrix mat at position row 'i' and column 'j'
 */
double getValue(matrix mat,  uint8_t i, uint8_t j){
    // Check that the position to access is allowed or no
    if(i > mat.dim || j > mat.dim){
        printf("Tried to access unreachable value in the matrix : (%u, %u) with size of matrix : (%u, %u)", i, j, mat.dim, mat.dim);
        exit(1);
    }
    return mat.values[i*mat.dim+j];
}

/*
 * This function set a value in matrix mat at position row 'i' and column 'j'
 */
double setValue(matrix mat,  uint8_t i, uint8_t j, double value){

    // Check that the position to access is allowed or no
    if(i > mat.dim || j > mat.dim){
        printf("Tried to modify unreachable value in the matrix : (%u, %u) with size of matrix : (%u, %u)", i, j, mat.dim, mat.dim);
        exit(1);
    }
    mat.values[i*mat.dim+j] = value;
}

void printMat(matrix mat){
    for(uint8_t i = 0; i < mat.dim; i ++){
        for(uint8_t j = 0; j < mat.dim; j ++){
            printf(" (%u, %u) : %f\t", i, j, getValue(mat, i, j));
        }
        printf("\n");
    }
}

/*
 * This function computes the maximum among maxi = {mat_i,j} where i>= indice and return i st maxi == mat_i,j
 */
int getLineMax(matrix mat, uint8_t indice, uint8_t j){
    uint8_t lineMax = indice;
    double maxValue = getValue(mat, indice, j);
    for(uint8_t i = indice + 1; i < mat.dim; i++){
        if(getValue(mat, i, j) > maxValue){
            maxValue = getValue(mat, i, j);
            lineMax = i;
        }
    }
    return lineMax;
}

/*
 * This function devide the row i by a denominator
 */
void divideLine(matrix mat, uint8_t i, double denominator){
    for(uint8_t j = 0; j < mat.dim; j ++){
        mat.values[mat.dim*i + j] /= denominator;
    }
}

/*
 * swapLines swaps row1 and row2 of matrix mat
 */
void swapLines(matrix mat, uint8_t row1, uint8_t row2){
    double tmp = 0;
    for(uint8_t j = 0; j < mat.dim; j ++){
        tmp = mat.values[mat.dim*row1 + j];
        mat.values[mat.dim*row1 + j] = mat.values[mat.dim*row2 + j];
        mat.values[mat.dim*row2 + j] = tmp;
    }
}


/*
 *  SubstractLines does substraction between two row row1 and row2  starting from index j=currentColumn
 */
void substractLines(matrix mat, uint8_t row1, uint8_t row2, uint8_t currentColumn){
    double multiplier = mat.values[mat.dim*row1 + currentColumn];
    for(uint8_t j = currentColumn; j < mat.dim; j++){
        mat.values[mat.dim*row1 + j] -= mat.values[mat.dim*row2 + j] * multiplier;
    }
}


/*
 *  This function computes the scalar product between column col1 of matrix mat 1 and column col 2 of matrix mat2
 *  Here we assume that the two columns have the same size.
 */
double scalarProduct(matrix mat1, matrix mat2, uint8_t col1, uint8_t col2){
    double ret = 0.0;
    for(uint8_t i = 0; i < mat1.dim; i++){
        ret += getValue(mat1, i, col1) * getValue(mat2, i, col2);
    }
    return ret;
}

/*
 *  Compute the determinant of matrix mat using Gaus elimination
 */
double determinant(matrix* mat){
    // Here we create a copy for our matrix, because we don't want to modify its content
    matrix *copy = malloc(sizeof(matrix));
    memcpy(copy, mat, sizeof(matrix));
    copy->values = malloc(sizeof(double) *copy->dim * copy->dim);
    // Copy the content of mat into matrix named "copy"
    memcpy(copy->values, mat->values, sizeof(double) *copy->dim * copy->dim);
    uint8_t r = 0; // r is index of  line of the last pivot found
    
    double ret = 1.0; // the value of determinant
    for(uint8_t j = 0; j < copy->dim; j ++){
        uint8_t k = getLineMax(*copy, r, j); // the index of the row with the largest value in column j
        double pivot = getValue(*copy, k, j);
        if(abs(pivot) != 0.0){
            divideLine(*copy, k, pivot);// We normalize the rox of pivot in order pivot take 1
            ret *= pivot; // Since we divide a line by a contant so normaly we multiply  determinant by the same contant

            if(k != r){
                swapLines(*copy, k, r);
                ret *= -1.0; // because we did a swap
            }
            for(uint8_t i = 0; i < copy->dim; i++){
                if(i!=r){
                    substractLines(*copy, i, r, j); // Substract line i of line r in order to annul amt[i,j]
                }
            }
            r = r +1;
        }
    }
    for(uint8_t i = 0; i < copy->dim; i++){
        ret *= getValue(*copy, i , i); // we compute the determinant of the triangular matrix
    }
    return ret;
}


/*
 *  orthogonolizing the base mat using Gram Schmidt process and put this orthogonal set in a matrix gmbasis
 *  @mat     : the set of vectors we want to orthogonolize
 *  @gmbasis : the generated orthogonal set 
 *  @mu      : the intermidiates coeffictients during Gram Schmidt process
 *  @G       : the set of norms of vectors in gmbasis basis
 */
void gramSchmidt(matrix* mat, matrix* gmbasis, matrix* mu, double *G){ 

    // Create a copy named gmbasisof matrix mat
    memcpy(gmbasis, mat, sizeof(matrix));
    gmbasis->values = malloc(sizeof(double) *gmbasis->dim * gmbasis->dim);
    memcpy(gmbasis->values, mat->values, sizeof(double) *gmbasis->dim * gmbasis->dim);

    // The matrix that will contain all intermidiates coeffictients
    memcpy(mu, mat, sizeof(matrix));
    mu->values = malloc(sizeof(double) *mu->dim * mu->dim);

    G[0] = 1.0;
    for(uint8_t i = 0; i < mat->dim; i++){

        for(uint8_t j = 0; j < i; j++){
            setValue(*mu, i, j, scalarProduct(*mat, *gmbasis, i, j)/G[j]);
            for(uint8_t k = 0; k < gmbasis->dim; k++){
                setValue(*gmbasis, k, i, getValue(*gmbasis, k, i) - getValue(*mu, i, j) * getValue(*gmbasis, k, j));
            }
        }
        G[i] = scalarProduct(*gmbasis, *gmbasis, i, i);
        if(G[i] == 0){
            printf("Non invertible matrix");
            exit(1);
        }
    }
}

/*
 *  This the reduction function in LLL algorithm (replace fk by fk-qfl)
 *  @mat     : the set of vectors we want to reduce (f1,f2,..,fn) mat's columns
 *  @mu      : the intermidiates coeffictients computed during Gram Schmidt process
 *  k        : the index of the vector we want to reduce
 *  l        : the index of the vector that will be used to reduce vector k
 */

void red(matrix* mat, matrix* mu, uint8_t k, uint8_t l){
    // Check if the vector k need reduction
    if(getValue(*mu, k, l) <= 0.5){
        return;
    }else{
        int64_t q = floor(getValue(*mu, k, l)+0.5);
        for(uint8_t i = 0; i < mat->dim; i++){
            setValue(*mat, i, k, getValue(*mat, i, k) - q* getValue(*mat, i, l));
        }
        // Update Gram-Schmidt coefficients
        setValue(*mu, k, l, getValue(*mu, k, l)-q);
        for(uint8_t i = 0; i < l; i++){
            setValue(*mu, k, i, getValue(*mu, k, i) - q* getValue(*mu, l, i));
        }
    }
}

/*
 * Exchange  column indexed k with the one indexed k-1
 * And make the corresponding updates to Gram-Schmidt coefficients
 */
void swap(matrix* mat, matrix* mu, double* G, uint8_t k){
    for(uint8_t i = 0; i < mat->dim; i ++){
        double tmp = mat->values[mat->dim*i + k];
        mat->values[mat->dim*i + k] = mat->values[mat->dim*i + k -1];
        mat->values[mat->dim*i + k -1] = tmp;
    }
    // Update Gram-Schmidt coefficients
    if(k > 1){
        for(uint8_t i = 0; i < k -1; i++){
            double tmp = getValue(*mu, k-1, i);
            setValue(*mu, k-1, i, getValue(*mu, k, i));
            setValue(*mu, k, i, tmp);
        }
    }
    double muLocal = getValue(*mu, k, k-1);
    double GLocal = G[k] + pow(muLocal, 2)*G[k-1];
    setValue(*mu, k, k-1, muLocal*G[k-1]/GLocal);
    G[k] = G[k-1]*G[k]/GLocal;
    G[k-1] = GLocal;
    for(uint8_t i = k+1; i < mat->dim; i++){
        double t = getValue(*mu, i, k);
        setValue(*mu, i, k, getValue(*mu, i, k-1)- muLocal*t);
        setValue(*mu, i, k-1, t + getValue(*mu, k, k-1) * getValue(*mu, i, k));
    }

}

/*
 * Verification property: It checks if the basis is lll reduced or no
 */
int verifyLLLBasis(matrix *gmbasis_org, matrix *mu, double c){
    matrix *gmbasis = malloc(sizeof(matrix));
    gmbasis->dim = gmbasis_org->dim;
    memcpy(gmbasis, gmbasis_org, sizeof(matrix));
    gmbasis->values = malloc(sizeof(double) *gmbasis->dim * gmbasis->dim);
    memcpy(gmbasis->values, gmbasis_org->values, sizeof(double) *gmbasis->dim * gmbasis->dim);
    for(uint8_t i = 1; i < mu->dim; i++){
        for(uint8_t j = 0; j < mu->dim; j++){
            setValue(*gmbasis,j,i, getValue(*gmbasis_org, j, i) +  getValue(*mu, i, i-1)*getValue(*gmbasis_org, j, i-1));
        }
    }

    for(uint8_t i = 0; i < mu->dim; i++){
        for(uint8_t j = i+1; j < mu->dim; j++){

            if(abs(getValue(*mu, i, j)) > 0.5){
                printf("1. The basis is not LLL reduced! \n");
                free(gmbasis->values);
                free(gmbasis);
                gmbasis = NULL;
                return 0;
            }
        }
        if(i > 1 && c*scalarProduct(*gmbasis_org, *gmbasis_org, i-1, i-1) > scalarProduct(*gmbasis, *gmbasis, i, i)){
            free(gmbasis->values);
            free(gmbasis);
            gmbasis = NULL;
            printf("2. The basis is not LLL reduced! \n");
            return 0;
        }
    }
    free(gmbasis->values);
    free(gmbasis);
    gmbasis = NULL;
    printf("Verified : The basis is  LLL reduced! \n");
    return 1;
}

/*
 * It reduces a given basis `mat` using LLL algorithm with a parameter c
 * @mat : (f1,f2,..,fn) are mat's column 
 */
void LLL(matrix* mat, double c){
    matrix *gmbasis = malloc(sizeof(matrix));
    matrix *mu = malloc(sizeof(matrix));
    double G[mat->dim];
    // Compute the corresponding basis using Gram-Schmidt process
    gramSchmidt(mat, gmbasis, mu, G);

    
    //verifyLLLBasis(gmbasis, mu, c); // Décommentez pour voir  si la base est initialement réduite ou non

    uint8_t k = 1;

    while(k < mat->dim){
        // replace fk by fk-m[k,k-1]*f(k-1)
        red(mat, mu, k, k-1);

        // Check if we had to exchange fk and f(k-1)
        if(G[k] < (c - pow(getValue(*mu, k, k-1), 2)) * G[k-1]){
            swap(mat, mu, G, k);
            k = 1 > k -1 ? 1 : k-1;

            // goto : red(mat, mu, k, k-1);
            continue;
        }

        if(k > 1){
            // reduce completly column from k
            for(uint8_t l = k - 2; l >= 0 && l <= k - 2; l--){ // On a ajouté la 2eme condition (l <= k - 2) à cause de Overflow si l==0 => l-1==255
                red(mat, mu, k, l);
            }
        }

        k = k+1;
    }

    // Pour tester si la base est réduite décommentez les deux lignes suivantes

    // gramSchmidt(mat, gmbasis, mu, G);
    // verifyLLLBasis(gmbasis, mu, c);
}

/*
 * It generates a random matrix of size rank x rank
 */
void generateRandMatrix(matrix *mat, uint8_t rank){
    mat->dim = rank;
    mat->values = malloc(sizeof(double) *mat->dim * mat->dim);
    for(uint8_t i = 0; i < rank; i++){
        for(uint8_t j = 0; j < rank; j++){
            setValue(*mat, i, j, (float)rand()/(float)(RAND_MAX/200)-100.0);
        }
    }
}


/*
 * Test determinant function on a small matrix of rank 3
 */
void test_dtr(){
    printf("==> Running test on determinant function \n");
    matrix* mat = malloc(sizeof(matrix));
    mat->dim = 3;
    mat->values = malloc(sizeof(double) *mat->dim * mat->dim);
    mat->values[0] = 1;
    mat->values[1] = 2;
    mat->values[2] = 0;
    mat->values[3] = 1;
    mat->values[4] = 0;
    mat->values[5] = 3;
    mat->values[6] = 2;
    mat->values[7] = 1;
    mat->values[8] = 3;
    double det = determinant(mat);
    double toFind = 3.0;
    printf("Found : %f\n", det);
    printf("To find : %f\n", toFind);

} 

void small_test_lll(){
    printf("==> Running a small test of lll \n");
    matrix* mat = malloc(sizeof(matrix));
    mat->dim = 2;
    mat->values = malloc(sizeof(double) *mat->dim * mat->dim);
    mat->values[0] = 201;
    mat->values[1] = 1648;
    mat->values[2] = 37;
    mat->values[3] = 297;
    printMat(*mat);
    LLL(mat,3.0/4);
    printf("After LLL transformation: \n");
    printMat(*mat);

} 

/*
 * Test to verify if our implémentation work or no
 * To see the result please decomment the two last lines LLL() function
 */
void test_LLL(){
    matrix* mat = malloc(sizeof(matrix));
 
    uint8_t rank = 100;
    double c = 3.0/4;
    generateRandMatrix(mat, rank);

    LLL(mat,c);


}


/*
 * Test of performance
*/
void test_performance(){
    printf("==> Running performance test of lll \n");
    matrix* mat = malloc(sizeof(matrix));
 
    uint8_t rank = 100;
    double c = 3.0/4;
    generateRandMatrix(mat, rank);

    matrix* copy = malloc(sizeof(matrix));
    for(c = 0.25; c < 1.0; c+=0.01){
        memcpy(copy, mat, sizeof(matrix));
        clock_t time_beg = clock();
        LLL(copy, c);
        double time_computation = (double)(clock()- time_beg)/(CLOCKS_PER_SEC);
        matrix *gmbasis = malloc(sizeof(matrix));
        matrix *mu = malloc(sizeof(matrix));
        double G[copy->dim];
        gramSchmidt(copy, gmbasis, mu, G);
        printf("For c = %f\nThe time needed is : %fsec\nThe shortest lenghth is : %lf\n", c, time_computation, scalarProduct(*mat, *mat, 0, 0));
    }
    

}

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    test_dtr();
    small_test_lll();
    test_LLL();
    //test_performance();
    return 0;
}

