#include <iostream>
#include <cmath>
#include <fstream>

#define masa_sol 1.99e30
#define c 1.496e11  //Distancia Sol-Tierra
#define G 6.67191e-11  //Cte gravitacion
#define h 0.01  //Intervalo de tiempo
#define n 10
#define T 10000 //Iteraciones



using namespace std;

int main(){

    int i, j, k;
    double masa[n], posicion[2][n], velocidad[2][n], aceleracion[2][n]; //Se definen los vectores que portaran los datos de la masa, posicion y velocidad
    double w[2][n];
    double t, norma, energia, momento, periodo;
    bool condicion[9];



    //Abrimos los archivos con los datos de las masas, posiciones y velocidades iniciales
    ifstream m;
    ifstream posicion_ini;
    ifstream velocidad_ini;
    m.open("masas.txt");
    posicion_ini.open("posiciones.txt");
    velocidad_ini.open("velocidades.txt");


    //Definimos el fichero donde vamos a escribir las posiciones de los planetas
    ofstream r_pos;
    ofstream energia_momento;
    ofstream periodos;
    ofstream geocentrico;
    r_pos.open("resultados_posiciones.txt");
    energia_momento.open("Energia_Momento.txt");
    periodos.open("periodos.txt");
    geocentrico.open("geocentrico.txt");



    //Inicializamos el booleano a true
    for(i=0;i<9;i++){
        condicion[i]=true;
    }



    //Volcamos los datos de las masas en el vector correspondiente y la reescalamos dividiendo por la masa del sol
    for(i=0;i<n;i++){
        m>>masa[i];
        masa[i]=masa[i]/masa_sol;
    }

    //Volcamos los datos de las posiciones en el vector correspondiente y reescalamos diviendo por la distancia sol-tierra (c)

    posicion[0][3]=149.6e9/c;
    posicion[1][3]=0.0;
    for(i=0;i<n;i++){
        posicion_ini>>posicion[0][i];
        posicion[0][i]=posicion[0][i]/c;
        posicion[1][i]=0.0; //Posiciones del eje Y iniciales
        r_pos<<posicion[0][i]<<"\t"<<posicion[1][i]<<"\t";

        geocentrico<<posicion[0][i]-posicion[0][3]<<"\t"<<posicion[1][i]-posicion[1][3]<<"\t";
    }
    r_pos<<endl;
    geocentrico<<endl;

    //Volcamos los datos de las velocidades en el vector correspondiente 
    for(i=0;i<n;i++){
        velocidad_ini>>velocidad[1][i]; //Copiamos solo las velocidades del eje Y ya que seran las no nulas inicialmente
        velocidad[0][i]=0.0; //Hacemos que la componente horizontal de la velocidad inicialmente sea 0
        velocidad[1][i]=velocidad[1][i]*sqrt(c/(G*masa_sol));  //Reescalamos la velocidad
    }


    //Le damos un valor inicial a t de 0
    t=0.0;



    //Calculamos la aceleracion 
    for(i=0;i<n;i++){
        //Inicilizamos a 0 todas las componentes de la aceleracion
        aceleracion[0][i]=0.0;
        aceleracion[1][i]=0.0;

        for(j=0;j<n;j++){
            if(i!=j){
                //Calculamos la norma a parte para simplificar expresiones posteriores
                norma=sqrt(((posicion[0][i]-posicion[0][j])*(posicion[0][i]-posicion[0][j]))+((posicion[1][i]-posicion[1][j])*(posicion[1][i]-posicion[1][j])));

                //Se calculan las componentes de la aceleracion
                aceleracion[0][i]=aceleracion[0][i]-((masa[j]*(posicion[0][i]-posicion[0][j]))/(norma*norma*norma));
                aceleracion[1][i]=aceleracion[1][i]-((masa[j]*(posicion[1][i]-posicion[1][j]))/(norma*norma*norma));

            }
        }

    }




    for(k=0;k<=T;k++){
            
        //Calculamos las nuevas posiciones y omegas
        posicion[0][0]=0.0;
        posicion[1][0]=0.0;  //Posicion del sol

        for(i=1;i<n;i++){
            posicion[0][i]=posicion[0][i]+h*velocidad[0][i]+(0.5*h*h*aceleracion[0][i]);
            posicion[1][i]=posicion[1][i]+h*velocidad[1][i]+(0.5*h*h*aceleracion[1][i]);

            w[0][i]=velocidad[0][i]+(0.5*h*aceleracion[0][i]);
            w[1][i]=velocidad[1][i]+(0.5*h*aceleracion[1][i]);
        }



        energia=0.0;
        momento=0.0;

        //Calculamos las nuevas aceleraciones en base a las posiciones obtenidas
        for(i=0;i<n;i++){

            aceleracion[0][i]=0.0;
            aceleracion[1][i]=0.0;

            for(j=0;j<n;j++){
                if(i!=j){

                    //Calculamos la norma a parte para simplificar expresiones posteriores
                    norma=sqrt(((posicion[0][i]-posicion[0][j])*(posicion[0][i]-posicion[0][j]))+((posicion[1][i]-posicion[1][j])*(posicion[1][i]-posicion[1][j])));

                    //Se calculan las componentes de la aceleracion
                    aceleracion[0][i]=aceleracion[0][i]-((masa[j]*(posicion[0][i]-posicion[0][j]))/(norma*norma*norma));
                    aceleracion[1][i]=aceleracion[1][i]-((masa[j]*(posicion[1][i]-posicion[1][j]))/(norma*norma*norma));


                    //Calculamos la energia
                    energia=energia+0.5*masa[i]*sqrt((pow(velocidad[0][i], 2) + pow(velocidad[1][i], 2)))-((masa[i]*masa[j])/sqrt(pow((posicion[0][i]-posicion[0][j]), 2)+pow((posicion[1][i]-posicion[1][j]), 2)));

                    //Calculamos el momento angular
                    momento=momento+sqrt(pow((posicion[0][i]-posicion[0][j]), 2)+pow((posicion[1][i]-posicion[1][j]), 2))*masa[i]*sqrt(pow(velocidad[0][i]-velocidad[0][j], 2)+pow(velocidad[1][i]-velocidad[1][j], 2));

                }
            }
        }


        energia_momento<<t<<"\t"<<energia<<"\t"<<momento<<endl;





        //Se calculan las nuevas velocidades
        velocidad[0][0]=0.0;
        velocidad[1][0]=0.0;
        for(i=1;i<n;i++){
            velocidad[0][i]=w[0][i]+(0.5*h*aceleracion[0][i]);
            velocidad[1][i]=w[1][i]+(0.5*h*aceleracion[1][i]);

        }


        //Se suma el intervalo de tiempo
        t=t+h;


        //Copiamos las nuevas posiciones al fichero
        for(i=0;i<n;i++){
            r_pos<<posicion[0][i]<<"\t"<<posicion[1][i]<<"\t";


            geocentrico<<posicion[0][i]-posicion[0][3]<<"\t"<<posicion[1][i]-posicion[1][3]<<"\t";

        }
        r_pos<<endl;
        geocentrico<<endl;



        periodo=0.0;
        //Comprobamos si cada planeta ha dado media vuelta, y en caso afirmativo calculamos su periodo multiplicando por dos el periodo que lleve
        for(i=1;i<9;i++){
            if((posicion[1][i]<0.0)&&(condicion[i]==true)){
                condicion[i]=false;
                //Calculamos periodo en DIAS, lo mostramos en pantalla e imprimimos en fichero
                periodo=2*t/(60*60*24*sqrt(G*masa_sol/(c*c*c)));
                cout<<"El periodo del planeta "<<i<<" es: "<<periodo<<" dias."<<endl;
                periodos<<"Periodo planeta "<<i<<": "<<periodo<<" dias"<<endl;

            }

        }
        

    }


    r_pos.close();
    energia_momento.close();
    periodos.close();

    return 0;
}
