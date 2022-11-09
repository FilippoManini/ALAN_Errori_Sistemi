// Gruppo: 
// Dellepiane Emanuele, 4876072
// Manini Filippo, 4798004
// Miggiano Davide, 4840761

#include <iostream>
#include <cmath>
using namespace std;

//f = funzione
//x = valori della funzione
//N = grado [max] del polinomio 

double fattoriale(int x)
{
	if(x<0)
		return -1;

	if(x == 0)
		return 1;
	else
		return x*fattoriale(x-1);
}

double polinomioTaylor(int N, double x)
{
    double ris = 0.0;
    for(int n=0; n <= N; n++)
        ris+= (double)(pow(x,n))/fattoriale(n);

    return ris;
}

double erroreAssoluto(double valorePerturbato, double valore)
{
	return fabs(valorePerturbato-valore);
}

double erroreRelativo(double erroreAssoluto, double valore)
{
	if(valore == 0)
		return(-1);	// Divisione per 0!

	return erroreAssoluto/valore;
}

void algoritmo1(int N, double x)
{
    double ris1 = polinomioTaylor(N,x);
    double ris2 = exp(x);
    //cout << "Risultato Polinomio di Taylor: " << ris1 << endl;
    cout << "Risultato Funzione exp: " << ris2 << endl; 

    cout << "Errore assoluto: " << erroreAssoluto(ris1, ris2) << endl;
    cout << "Errore relativo: " << erroreRelativo(erroreAssoluto(ris1,ris2), ris2) << endl;
}

void algoritmo2(int N, double x)
{
    double ris1 = polinomioTaylor(N,-x);
    double ris2 = 1/ris1;
    //cout << "Risultato Polinomio di Taylor: " << ris1 << endl;
    cout << "Risulatato Reciproco del Polinomio di Taylor: "<< ris2 <<endl;

    cout << "Errore assoluto: " << erroreAssoluto(ris2, exp(x)) << endl;
    cout << "Errore relativo: " << erroreRelativo(erroreAssoluto(ris2, exp(x)), ris1) << endl;
}

int main()
{
    int N = 0, choice = 0;
    double x = 0.0;

    while(1)
    {
        cout << "Scegli quale algoritmo eseguire:\n1-Algoritmo 1\n2-Algoritmo 2\n0-Termina"<<endl;
        cin >> choice;

        if (choice == 0) break;

        switch (choice)
        {
            case 1:
                cout << "Inserisci il grado: ";
                cin >> N;
                cout << "Inserisci il punto: ";
                cin >> x;
                algoritmo1(N,x);
                cout << "\n-------------------------------------"<< endl;
                break;
            
            case 2:
                cout << "Inserisci il grado: ";
                cin >> N;
                cout << "Inserisci il punto: ";
                cin >> x;
                algoritmo2(N,x);
                cout << "\n-------------------------------------"<< endl;
                break;

            default:
                break;
        }
    }
    return 0;

}
