// Gruppo: 
// Dellepiane Emanuele - 4876072
// Manini Filippo - 4798004
// Miggiano Davide - 4840761

#include <iostream>
#include <cmath>

using namespace std;

double** createMatrix(int row, int column) {
	if(row < 1)
		return NULL;

	double** matrice = new double*[row];

	for(int i = 0; i < row; ++i)
		matrice[i] = new double[column];

	return matrice;
}

double** fillMatrix(double **matrix, int row, int column, double *vector){
    if(row < 1 || column < 0 )
        return NULL;

    int counter = 0;

    for(int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            matrix[i][j] = vector[counter];
            counter++;
        }
    }
    return matrix;
}

void copymatrix(double **matrix, double **matrixCopy, int row, int column){
    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < column; ++j)
        {
            matrixCopy[i][j] = matrix[i][j];
        }
    }
}

void printMatrix(double **matrix, int row, int column, int type){
	for(int i = 0; i < row; ++i)
	{
		for(int j = 0; j < column; ++j)
		{
            if(type == 0)
			    cout << matrix[i][j] << "\t";
            if(type == 1)
                cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

double fattoriale(int x){
	if(x<0)
		return -1;

	if(x == 0)
		return 1;
	else
		return x*fattoriale(x-1);
}

double findMax(double* b, int elem){
    double max = b[0];
    for(int i = 1; i <= elem; ++i)
    {
        if(b[i] > max)
            max = b[i];
    }
    return max;
}

double* pascalVector(int row, int column){
    double *vector = new double[row*column];
	int counter = 0;

	for(int i = 1; i <= row; i++)
    {
        for (int j = 1; j <= column; j++)
        {
            vector[counter] = (double)(fattoriale(i+j-2) / (fattoriale(i-1) * fattoriale(j-1)));
            counter++;
        }
    }
	return vector;
}

double* tridiagonalVector(int row, int column){
    double *vector = new double[row*column];
    int counter = 0;

    for(int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            if(i == j)
                vector[counter] = 2.0;

            else if(fabs(i-j) == 1)
                vector[counter] = -1.0;

            else
                vector[counter] = 0.0;
            
            counter++;
        }
    }
    return vector;
}

double* systemSolution(int num, int dim){
    double *x = new double[dim];
    for(int i = 0; i< dim; ++i)
        x[i] = num;
    return x;
}

double* termineNoto(double **matrix, double* x, int row, int column){
    double *b = new double[row];
    int counter = 0;

    for (int i = 0; i < row; ++i)
    {
        b[i] = 0;
        for (int j= 0; j < column; ++j)
        {
            b[counter] += ( x[j] * matrix[i][j]);
        }
        ++counter;
    }
    return b;
}

double* termineNotoPerturbato(double* x, int row){
    double normaX = findMax(x,row);
    double *tNoto = new double[row];
    
    for(int i = 0; i < row; ++i)
    {
        if(i%2 == 0)
            tNoto[i] = normaX * (-0.01) + x[i];
        else
            tNoto[i] = normaX * (0.01) + x[i] ;
        
    }
    return tNoto;
}

void stampaTermineNoto(double *b, int row){
    for(int i = 0; i < row; ++i)
        cout<<b[i]<<" ";
    cout<<"\n";
}

double** mergeMatrixVector(double** matrix, double* vector, int row, int column, int dimVector ){
    if(row != dimVector)
        return NULL;
    
    double** mergeMatrix = createMatrix(row,column+2);

    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < column; ++j)
            mergeMatrix[i][j] = matrix[i][j];
        
        mergeMatrix[i][column] = vector[i];
    }
    return mergeMatrix;
}

void switchRow(double** matrix, int row, int column, int rowPivot, int columnPivot){
    if(columnPivot > rowPivot)
        return;
    
    for(int i = rowPivot+1; i <row; ++i)
    {
        if(matrix[i][columnPivot] != 0)
        {
            double *temp = new double[column];
            for(int j = 0; j < column; ++j)
            {
                temp[j] = matrix[rowPivot][j];
                matrix[rowPivot][j] = matrix[i][j];
                matrix[i][j] = temp[j];
            }
            break;
        }
    }
}

double** gauss(double** matrix, double* vector, int row, int column){
    double** matrixAB = mergeMatrixVector(matrix,vector,row,column,row);
    float pivot,m;

    for(int k=0; k<row-1; ++k)
    {
        pivot = matrixAB[k][k];
        if(pivot != 0)
        {
            for(int i = k+1; i < row; ++i)
            {
                m = matrixAB[i][k] / pivot;
                //Per evitare cicli inutili (caso elemento colonna già ridotta)
                if(m != 0)
                {
                    for(int j = k; j < column+1; ++j)
                    {
                        matrixAB[i][j] -=  (m * matrixAB[k][j]);
                    }
                }
            }
        }
        else
        {
            cout<<"Invalid Pivot...switching rows\n";
            // PIVOT NULLO, DEVO FARE SOSTITUZIONE RIGA!!
            //switchRow(matrix, row, column, k,k);
        }
        
    }
    return matrixAB;
}

float* backSubstitution(double** matrix, int row, int column){
    // Prende una matrice ridotta precedentemente da Gauss e ne calcola le soluzioni

    float* solution = new float[row];
    float adder;

    // Estraggo il vettore che avevo unito alla matrice
	for(int i = 0; i < row; ++i)
		solution[i] = matrix[i][column];

    //sostituzione all'indietro
    for(int i = row-1; i >= 0; --i)
    {
        adder = 0.0;
        for(int j = i+1; j < row; ++j)
            adder += matrix[i][j] * solution[j];
        
        solution[i] = (solution[i] - adder) / matrix[i][i];
    }

    return solution;
}

void printSolution(float *vector, int dim, int label){
	switch(label)
	{
		case 11:
			cout << "\n- Soluzioni Gauss Matrice A1 -\n";
			break;
        case 12:
            cout << "\n- Soluzioni Gauss Matrice A1 Perturbata -\n";
			break;
		case 21:
			cout << "\n- Soluzioni Gauss Matrice A2 -\n";
			break;
        case 22:
			cout << "\n- Soluzioni Gauss Matrice A2 Perturbata -\n";
			break;
		case 31:
			cout << "\n- Soluzioni Gauss Matrice Pascal -\n";
			break;
        case 32:
			cout << "\n- Soluzioni Gauss Matrice Pascal Perturbata -\n";
			break;
		case 41:
			cout << "\n- Soluzioni Gauss Matrice Tridiagonale -\n";
			break;
		case 42:
			cout << "\n- Soluzioni Gauss Matrice Tridiagonale Perturbata -\n";
			break;
		default:
			cout << "\n- Soluzioni Gauss Matrice -\n";
			break;
	}

	for(int i = 0; i < dim; ++i)
	{
		cout << "x" << i+1 << " = " << vector[i] << "\n";
	}
	cout << endl;
}

void clean(double** matrix,double* vector){
    delete [] matrix;
    delete [] vector;
}



// ----------------------------------------------------------------------------
int main()
{
    int d0 = 1, d1 = 6;
    const int dimMatrixA = 4;               //dimensione matrici A
    const int dimMatrixP = 10;              //dimensione matrice di Pascal
    const int dimMatrixT = 10*(d1+1)+d0;    //dimensione matrice Tridiagonale

    double **A1Copy = createMatrix(dimMatrixA,dimMatrixA);
    double **A2Copy = createMatrix(dimMatrixA,dimMatrixA);
    double **PCopy = createMatrix(dimMatrixP,dimMatrixP);
    double **TCopy = createMatrix(dimMatrixT,dimMatrixT);

    //vettori contenenti i valori per le matrici A
    double a1[16] = {3, 1, -1, 0, 0, 7, -3, 0, 0, -3, 9, -2, 0, 0, 4, -10};
    double a2[16] = {2, 4, -2, 0, 1, 3, 0, 1, 3, -1, 1, 2, 0, -1, 2, 1};

    //Matrici A
    double **matrixA1 = fillMatrix(createMatrix(dimMatrixA,dimMatrixA), dimMatrixA, dimMatrixA,a1);
    double **matrixA2 = fillMatrix(createMatrix(dimMatrixA,dimMatrixA), dimMatrixA, dimMatrixA,a2);
    //Matrice di Pascal
    double **matrixP = fillMatrix(createMatrix(dimMatrixP,dimMatrixP), dimMatrixP, dimMatrixP,pascalVector(dimMatrixP,dimMatrixP));
    //Matrice Tridiagonale
    double **matrixT = fillMatrix(createMatrix(dimMatrixT,dimMatrixT), dimMatrixT, dimMatrixT,tridiagonalVector(dimMatrixT,dimMatrixT));


    //Calcolo Termine Noto per ogni matrice
    double *termineNotoA1 = termineNoto(matrixA1, systemSolution(1,dimMatrixA),dimMatrixA,dimMatrixA);
    double *termineNotoPert_A1 = termineNotoPerturbato(termineNotoA1,dimMatrixA);  
    
    double *termineNotoA2 = termineNoto(matrixA2, systemSolution(1,dimMatrixA),dimMatrixA,dimMatrixA);
    double *termineNotoPert_A2 = termineNotoPerturbato(termineNotoA2,dimMatrixA);
    
    double *termineNotoP = termineNoto(matrixP, systemSolution(1,dimMatrixP),dimMatrixP,dimMatrixP);
    double *termineNotoPert_P = termineNotoPerturbato(termineNotoP,dimMatrixP);

    double *termineNotoT = termineNoto(matrixT, systemSolution(1,dimMatrixT),dimMatrixT,dimMatrixT);
    double *termineNotoPert_T = termineNotoPerturbato(termineNotoT,dimMatrixT);

    //Copia matrici per successiva sostituzione all'indietro
    copymatrix(matrixA1,A1Copy,dimMatrixA,dimMatrixA);
    copymatrix(matrixA2,A2Copy,dimMatrixA,dimMatrixA);
    copymatrix(matrixP,PCopy,dimMatrixP,dimMatrixP);
    copymatrix(matrixT,TCopy,dimMatrixT,dimMatrixT);


    //Risoluzione sistema Ax = b tramite Gauss in precisione singola (float)
    matrixA1 = gauss(matrixA1, termineNotoA1,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(matrixA1,dimMatrixA,dimMatrixA),dimMatrixA,11);
    A1Copy = gauss(A1Copy, termineNotoPert_A1,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(A1Copy,dimMatrixA,dimMatrixA),dimMatrixA,12);


    matrixA2 = gauss(matrixA2, termineNotoA2,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(matrixA2,dimMatrixA,dimMatrixA),dimMatrixA,21);
    A2Copy = gauss(A2Copy, termineNotoPert_A2,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(A2Copy,dimMatrixA,dimMatrixA),dimMatrixA,22);


    matrixP = gauss(matrixP, termineNotoP,dimMatrixP,dimMatrixP);
    printSolution(backSubstitution(matrixP,dimMatrixP,dimMatrixP),dimMatrixP,31);
    PCopy = gauss(PCopy,termineNotoPert_P,dimMatrixP,dimMatrixP);
    printSolution(backSubstitution(PCopy,dimMatrixP,dimMatrixP),dimMatrixP,32);

    matrixT = gauss(matrixT, termineNotoT,dimMatrixT,dimMatrixT);
    printSolution(backSubstitution(matrixT,dimMatrixT,dimMatrixT),dimMatrixT,41);
    TCopy = gauss(TCopy, termineNotoPert_T,dimMatrixT,dimMatrixT);
    printSolution(backSubstitution(TCopy,dimMatrixT,dimMatrixT),dimMatrixT,42);

    clean(matrixA1,termineNotoA1);
    clean(matrixA2,termineNotoA2);
    clean(matrixP,termineNotoP);
    clean(matrixT,termineNotoT);
    return 0;

}
