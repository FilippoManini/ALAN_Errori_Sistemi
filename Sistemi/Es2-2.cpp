// Gruppo: 
// Dellepiane Emanuele - 4876072
// Manini Filippo - 4798004
// Miggiano Davide - 4840761

#include <iostream>
#include <cmath>

using namespace std;

double** createMatrix(int row, int column) 
{
	if(row < 1)
		return NULL;

	double** matrice = new double*[row];

	for(int i = 0; i < row; ++i)
		matrice[i] = new double[column];

	return matrice;
}

double** fillMatrix(double **matrix, int row, int column, double *vector)
{
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

void printMatrix(double **matrix, int row, int column, int type)
{
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

double fattoriale(int x)
{
	if(x<0)
		return -1;

	if(x == 0)
		return 1;
	else
		return x*fattoriale(x-1);
}

double* pascalVector(int row, int column)
{
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

double* tridiagonalVector(int row, int column)
{
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

double* systemSolution(int num, int dim)
{
    double *x = new double[dim];
    for(int i = 0; i< dim; ++i)
        x[i] = num;
    return x;
}

double* termineNoto(double **matrix, double* x, int row, int column)
{
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

void stampaTermineNoto(double *b, int row)
{
    for(int i = 0; i < row; ++i)
        cout<<b[i]<<" ";
    cout<<"\n";
}

double** mergeMatrixVector(double** matrix, double* vector, int row, int column, int dimVector)
{
    if(row != dimVector)
        return NULL;
    
    double** mergeMatrix = createMatrix(row,column+1);

    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < column; ++j)
            mergeMatrix[i][j] = matrix[i][j];
        
        mergeMatrix[i][column] = vector[i];
    }
    return mergeMatrix;
}

void switchRow(double** matrix, int row, int column, int rowPivot, int columnPivot)
{
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

double** gauss(double** matrix, double* vector, int row, int column)
{
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

float* backSubstitution(double** matrix, int row, int column)
{
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

void printSolution(float *vector, int dim, int label)
{
	switch(label)
	{
		case 1:
			cout << "\n- Soluzioni Gauss Matrice A1 -\n";
			break;
		case 2:
			cout << "\n- Soluzioni Gauss Matrice A2 -\n";
			break;
		case 3:
			cout << "\n- Soluzioni Gauss Matrice B -\n";
			break;
		case 4:
			cout << "\n- Soluzioni Gauss Matrice C -\n";
			break;
		default:
			cout << "\n- Soluzioni Gauss Matrice -\n";
			break;
	}

	for(int i = 0; i < dim; ++i)
	{
		cout << "x" << i+1 << " = " << vector[i] << "\t";
	}
	cout << endl;
}

void clean(double** matrix,double* vector)
{
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

    //Stampa delle diverse matrici
    cout<< "\nMatrice A1:\n";
    printMatrix(matrixA1,dimMatrixA,dimMatrixA,0);

    cout<< "\nMatrice A2:\n";
    printMatrix(matrixA2,dimMatrixA,dimMatrixA,0);

    cout<< "\nMatrice Pascal:\n";
    printMatrix(matrixP,dimMatrixP,dimMatrixP,0);

    cout<< "\nMatrice Tridiagonale:\n";
    printMatrix(matrixT,dimMatrixT,dimMatrixT,1);

    //Calcolo Termine Noto per ogni matrice
    double *termineNotoA1 = termineNoto(matrixA1, systemSolution(1,dimMatrixA),dimMatrixA,dimMatrixA);
    cout << "\nCalcolo termine noto matrice A1 = ";
    stampaTermineNoto(termineNotoA1,dimMatrixA);
    
    double *termineNotoA2 = termineNoto(matrixA2, systemSolution(1,dimMatrixA),dimMatrixA,dimMatrixA);
    cout << "\nCalcolo termine noto matrice A2 = ";
    stampaTermineNoto(termineNotoA2,dimMatrixA);

    double *termineNotoP = termineNoto(matrixP, systemSolution(1,dimMatrixP),dimMatrixP,dimMatrixP);
    cout << "\nCalcolo termine noto matrice Pascal = ";
    stampaTermineNoto(termineNotoP,dimMatrixP);

    double *termineNotoT = termineNoto(matrixT, systemSolution(1,dimMatrixT),dimMatrixT,dimMatrixT);
    cout << "\nCalcolo termine noto matrice Tridiagonale = ";
    stampaTermineNoto(termineNotoT,dimMatrixT);

    //Risoluzione sistema Ax = b tramite Gauss in precisione singola (float)
    matrixA1 = gauss(matrixA1, termineNotoA1,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(matrixA1,dimMatrixA,dimMatrixA),dimMatrixA,1);

    matrixA2 = gauss(matrixA2, termineNotoA2,dimMatrixA,dimMatrixA);
    printSolution(backSubstitution(matrixA2,dimMatrixA,dimMatrixA),dimMatrixA,2);

    matrixP = gauss(matrixP, termineNotoP,dimMatrixP,dimMatrixP);
    printSolution(backSubstitution(matrixP,dimMatrixP,dimMatrixP),dimMatrixP,3);

    matrixT = gauss(matrixT, termineNotoT,dimMatrixT,dimMatrixT);
    printSolution(backSubstitution(matrixT,dimMatrixT,dimMatrixT),dimMatrixT,4);

    clean(matrixA1,termineNotoA1);
    clean(matrixA2,termineNotoA2);
    clean(matrixP,termineNotoP);
    clean(matrixT,termineNotoT);
    return 0;

}
