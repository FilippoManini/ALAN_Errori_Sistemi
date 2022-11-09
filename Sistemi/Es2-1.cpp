// Gruppo: 
// Dellepiane Emanuele, 4876072
// Manini Filippo, 4798004
// Miggiano Davide, 4840761

/*
Norma infinito di una matrice:
max(somma elementi riga della matrice)
*/

#include <iostream>
#include <math.h>
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

void algoritmo1(double **matrix, int row, int column)
{
    double sum = 0.0, sumMax = 0.0;

    cout<< "Enter the values in the matrix:"<<endl;

    for(int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            cout<< "Enter A["<<i+1<<"]["<<j+1<<"]: ";
            cin >> matrix[i][j];
            //fabs: modulo
            sum += fabs(matrix[i][j]);
        }
        if(sum > sumMax)
            sumMax = sum;
        
        sum = 0.0;
    }
    cout<<"La norma infinito della matice A e': "<< sumMax <<endl;
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

void pascalMatrix(double **matrix,int row, int column)
{
    double sum, sumMax;

    for(int i = 1; i <= row; i++)
    {
        for (int j = 1; j <= column; j++)
        {
            matrix[i-1][j-1] = (double)(fattoriale(i+j-2) / (fattoriale(i-1) * fattoriale(j-1)));
            sum += fabs(matrix[i-1][j-1]);
        }
        if(sum > sumMax)
            sumMax = sum;
        
        sum = 0.0;
    }
    cout<<"La norma infinito della matice P e': "<< sumMax <<endl;
}

void tridiagonalMatrix(double **matrix, int row, int column)
{
    double sum, sumMax;

    for(int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            if(i == j)
                matrix[i][j] = 2.0;

            else if(fabs(i-j) == 1)
                matrix[i][j] = -1.0;

            else
                matrix[i][j] = 0.0;
            
            sum += fabs(matrix[i][j]);
        }
        if(sum > sumMax)
            sumMax = sum;
        
        sum = 0.0;
    }

    cout<<"La norma infinito della matice T e': "<< sumMax <<endl;
}

int main()
{
    double **matrix;
    int choice = 0, m, n, d0 = 1, d1 = 6;

    while(1)
    {
        cout<<"Scegliere di quali matrici calcolare la norma infinito:\n0-Termina\n1-Matrici inserite dall'utente\n2-Matrice di Pascal 10x10\n3-Matrice Tridiagonale nxn\n";
        cin >> choice;
        cout << "\n\n";

        if (choice == 0) break;

        switch(choice)
        {
            case 1:
                cout << "-------------------------------------------"<<endl;
                cout << "Enter no. of rows: ";
                cin >> m;
                cout << "Enter no. of columns: ";
                cin >> n;

                matrix = createMatrix(m,n);
                if(matrix != NULL)
                    algoritmo1(matrix,m,n);
                else
                    cout<<"Errore dimensione matrice\n";
                break;

            case 2:
                cout << "-------------------------------------------"<<endl;
                matrix = createMatrix(10,10);
                pascalMatrix(matrix,10,10);
                printMatrix(matrix,10,10,0);
                break;

            case 3:
                cout << "-------------------------------------------"<<endl;
                n = 10 * (d1+1) + d0;
                matrix = createMatrix(n,n);
                tridiagonalMatrix(matrix,n,n);
                printMatrix(matrix,n,n,1);
                break;

            default:
                break;

        }
        cout << "\n---------------------------------------"<<endl;
        delete [] matrix;
    }
    return 0;
}

