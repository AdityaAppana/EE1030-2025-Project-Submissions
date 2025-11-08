#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include<math.h>

__attribute__((visibility("default")))

void free_matrix(double** mat, int rows) {

    for (int i = 0; i < rows; i++) {
        free(mat[i]);
    }
    free(mat);
}

double** matmul(double**A, double**B, int Arows, int Acols, int Brows, int Bcols){
    
    if(Acols != Brows){
        printf("ERROR");
        return NULL;
    }
    
    else{
        
        double **C = (double**)malloc(Arows*sizeof(double*));
        for(int i =0; i<Arows; i++){
            C[i] = (double*)calloc(Bcols,sizeof(double));
        }
        
        
        for(int i = 0; i<Arows; i++){
            
            
            for( int k =0; k<Bcols; k++){
                for(int j =0; j<Acols; j++){
                    C[i][k] += A[i][j] * B[j][k];
                }
                
            }
        }
        
        return C;
        
    }
}

double norm (double *vec1, int dim){
    
    double squares = 0;
    
    for(int i =0; i<dim; i++){
        squares += vec1[i]*vec1[i];
    }
    double total = sqrt(squares);
    return total;
}



double* scalarmul( double*A, int dim, double cons){
    
    double*Anew = (double*)malloc(dim*sizeof(double));
    for(int i =0; i<dim; i++){
        Anew[i] = A[i]*cons;
    }
    
    return Anew;
}

double* qscalarmul( double*A, int dim, double cons){
    
    for(int i =0; i<dim; i++){
        A[i] = A[i]*cons;
    }
    
    return A;
}


double innerprod(double *v1, double*v2, int dim){
    double dot = 0;
    for(int i =0; i<dim; i++){
        dot+= v1[i]*v2[i];
    }
    
    return dot;
}

double* vecsub(double * v1, double*v2, int dim){
    
    double*Anew = (double*)malloc(dim*sizeof(double));
    for(int i =0; i<dim; i++){
        Anew[i] = v1[i]-v2[i];
    }
    return Anew;
}

double* qvecsub(double * v1, double*v2, int dim){
    
    for(int i =0; i<dim; i++){
        v1[i] = v1[i]-v2[i];
    }
    return v1;
}

double* vecadd(double * v1, double*v2, int dim){
    double*Anew = (double*)malloc(dim*sizeof(double));
    for(int i =0; i<dim; i++){
        Anew[i] = v1[i]+v2[i];
    }
    return Anew;
}

double** transpose( double** A, int Arows, int Acols){
    
    double **AT = (double**)malloc(Acols*sizeof(double*));
    for(int i =0; i<Acols; i++){
        AT[i] = (double*)malloc(Arows*sizeof(double));
    }
    
    for(int i =0; i<Arows; i++){
        for(int j = 0; j<Acols; j++){
            AT[j][i] = A[i][j];
        }
    }
    
    return AT;
}

double** qrdecompQ (double** A, int Arows, int Acols){
    
    double **Q = (double**)malloc(Acols*sizeof(double*));
    for(int i =0; i<Acols; i++){
        Q[i] = (double*)malloc(Arows*sizeof(double));
    }

    double**AT = transpose(A, Arows, Acols); //This will have dimensions Acols * Arows
    
    int n = Arows;
    
    double normval = norm(AT[0], n);
        if (normval == 0.0) normval = 1.0; 
        for (int k = 0; k < n; k++) {
            Q[0][k] = AT[0][k] / normval;
        }
    
    for( int i = 1; i<Acols; i++){
        for(int k =0; k<Arows; k++){
            Q[i][k] = AT[i][k];
        }
        
        for(int j = i; j>0; j--){

            double* proj = scalarmul( Q[j-1], n, innerprod(Q[j-1],AT[i],Arows));
            
        qvecsub(Q[i], proj , Arows);
            free(proj);
        }

        double normqi = norm(Q[i],n);
        if (normqi == 0.0){
            normqi = 1.0;
        }

        qscalarmul(Q[i], n, 1.0 / normqi);

    }
    
    double** QT = transpose(Q, Acols, Arows);
    free_matrix(Q, Acols);
    free_matrix(AT, Acols);
    return QT;
}

double** qrdecompR (double** A, double**Q, int Arows, int Acols){
    
    double**QT = transpose(Q, Arows, Acols);
    double **R = matmul(QT, A, Acols, Arows, Arows, Acols);
    
    free_matrix(QT, Acols);
return R;
    
}

double** matrixadd( double**A, double**B, int rows, int columns){
    
    for(int i =0; i<rows; i++){
        for(int j = 0; j<columns; j++){
            
            A[i][j] = A[i][j] + B[i][j];
    
        }
    }

    return A;
}

//So we can reconstruct the matrix more easily
double** outerproduct(double* u, double* v, int udim, int vdim) {
    
    double** C = (double**)malloc(udim * sizeof(double*));

    for (int i = 0; i < udim; i++) {
        C[i] = (double*)malloc(vdim * sizeof(double));

        for (int j = 0; j < vdim; j++) {
            C[i][j] = u[i] * v[j];
            
}
}

return C;

}

// So we can make the other eigen-vector matrix without running QR again
double*matrixvectormult (double** X, int Xrows, int Xcols, double * vector, double scalar){
    
    double * final = (double*)calloc( Xrows ,sizeof(double));
    for(int i =0; i<Xrows; i++){
        for(int j =0; j<Xcols; j++){
            final[i] += X[i][j]*vector[j];
        }
    }
    
    qscalarmul(final,Xrows,scalar);
    return final;
}


void copy_matrix_to_flat(double** matrix, double* flatarray, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            flatarray[i * cols + j] = matrix[i][j];
        }
    }
}

double** create_matrix_from_flat(double* flatarray, int rows, int cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = flatarray[i * cols + j];
        }
    }
    return matrix;
}
//start of main function


void svd(double* ogflat, double* sum_flat,int rows,int columns, int rank){


    double** og = create_matrix_from_flat(ogflat, rows, columns);

        double**buffer = transpose(og,rows,columns);
        
        double** XTX = matmul(buffer, og, columns,rows, rows,columns);

        
        free_matrix(buffer, columns);
        
        

double **Ixtx = (double**)malloc(columns*sizeof(double*));
            
        for(int i =0; i<columns; i++){
            Ixtx[i] = (double*)calloc(columns, sizeof(double));
                Ixtx[i][i] = 1;

        }
            //QR Iteration on XTX which gives eigenvectors of XTX
        //Obtained V
            
            for(int i =0; i<15; i++){
                
                double** Q = qrdecompQ(XTX, columns, columns);
                double** R = qrdecompR(XTX, Q, columns, columns);
                
                double** old_Ixtx = Ixtx;
                double** old_XTX = XTX;
                
                //Ixtx becomes the eigen-vector matrix
                Ixtx = matmul(old_Ixtx, Q, columns, columns, columns, columns);
                
                //XTX converges to the diagonal matrix whose diagonal elements are the squares of eigen-values
                XTX = matmul(R, Q, columns, columns, columns, columns);
        
                free_matrix(old_Ixtx, columns);
                free_matrix(old_XTX, columns);
                free_matrix(Q, columns);
                free_matrix(R, columns);
                
            }
        
        
        double*singvalues = (double*)malloc(sizeof(double)*columns);
        int*order= (int*)malloc(sizeof(int)*columns);
        
    
            
            for(int i =0; i<columns; i++){
                singvalues[i] = sqrt(fmax(0.0,XTX[i][i]));
                order[i] = i;
            
        }
        
        
        
        //Sorting the eigen-values
        
        for(int i =0; i<columns; i++){
            
            for(int j =0; j<columns - 1 -i; j++){
                
                if (singvalues[j] < singvalues[j + 1]){
                        double hold = singvalues[j];
                        singvalues[j] = singvalues[j + 1];
                        singvalues[j + 1] = hold;
                    
                    int spotholder = order[j];
                    order[j] = order[j+1];
                    order[j+1] = spotholder;
                }
            }
        }
        
        double ** sum = (double**)malloc(rows*sizeof(double));
    for(int l = 0; l<rows; l++){
        sum[l] = (double*)calloc(columns, sizeof(double));
    }

        for(int i =0; i<rank; i++){
            
        int p = order[i];
        double sigma = singvalues[i];

        if (sigma < 1e-15) continue;
        

        double* v_p = (double*)malloc(columns * sizeof(double));
    if (v_p == NULL) {
        continue; // or exit
    }

        for (int k = 0; k < columns; k++) {
            v_p[k] = Ixtx[k][p];
        }

        double* u_p = matrixvectormult(og, rows, columns, v_p, 1.0/sigma);
        double** rankapprox = outerproduct(u_p, v_p, rows, columns);
        
        for (int t = 0; t < rows; t++) {
            qscalarmul(rankapprox[t], columns, sigma);
        }

        matrixadd(sum, rankapprox, rows, columns);
        
        free(u_p);
        free_matrix(rankapprox, rows);
        free(v_p);
        

        }

        copy_matrix_to_flat(sum, sum_flat, rows, columns);

        free(singvalues);
    free(order);
    free_matrix(sum, rows);
    free_matrix(og, rows);
    free_matrix(XTX, columns);
    
    free_matrix(Ixtx, columns);
    }

