/*Ising Model in Two Dimensions -- Metropolis Algorithm*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
# define NSAMPLES 1000   // number of samples required for averaging
#define THSWEEPS 100000  // number of thermal sweeps 
# define N 10       // N by N lattice
int lattice[N][N];  




/* Initialize spins in the two dimensional lattice*/
void initialize(){

	int i, j;

	for (i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			lattice[i][j]=1;
		}
	}

}


/*Function for calculating the energy of a given configuration*/
float new_config()
{
	int i,j;
	float E=0.0;

	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++)
		{
			E=E-(1.0/2.0)*(lattice[i][j]*(lattice[(i+1)%N][j]+lattice[(i+N-1)%N][j]+lattice[i][(j+1)%N]+lattice[i][(j+N-1)%N]));
		}
	}
	return(E);
}



/* Function for calcuating Magnetization*/
float mag()
{
        int i,j,spin_up=0,spin_down=0;
        
        for(i=0;i<N;i++)
        {
                for(j=0;j<N;j++)
                {
                        if (lattice[i][j]==1)
                        {
                                spin_up++;
                        }
                        else
                        {
                                spin_down++;
                        }
                }
        }
        return((float)(spin_up-spin_down));
        
}


/* Function for generating Random Number */
double randomnumber()
{
        return((double) rand()/(double) (RAND_MAX+1.0));
        
}


/* Function for computing the square of a number*/

float sqr(float z)
{
        return(z*z);
}



/*One sweep of the Metropolis Algorithm*/
float sweep(float Beta)
{
        int i,j;
        float delta_E, E_initial;
	
	E_initial = new_config();

        for (i=0;i<N;i++)
        {
                for(j=0;j<N;j++)
                {
                        delta_E= 2*(lattice[i][j]*(lattice[(i+N-1)%N][j]+lattice[(i+1)%N][j]+lattice[i][(j+N-1)%N]+lattice[i][(j+1)%N]));
                        if (delta_E<=0 ||  exp(-(double)(Beta*delta_E))>randomnumber() )
                        {
                                E_initial=E_initial+delta_E;
                                lattice[i][j]=-lattice[i][j];
                        }
                }
        }
        return(E_initial);
}




int main(){

	float Beta;            // Inverse of Temperature
	int tau;   	       // iterator for thermalization		
	int sample_size;      // interator for averaging process
	float Mag_initial, E_initial;  // initial magnetization and energy
	float E_sweep;                // value of energy after one sweep of the metropolis algorithm
	float Avg_E, Avg_Esquared, Avg_Mag, Avg_Magsquared;
	float Total_E = 0.0,Total_Esquared = 0.0,Total_Mag = 0.0,Total_Magsquared=0.0;
	double Chi;      // Susceptiblity
	double Cv;	// Specific Heat
	double sigma;    // standard deviation of the Average Magnetization per site 
	FILE *ofp1, *ofp2, *ofp3;
	
	ofp1 = fopen("Avg_Magnetization_vs_Beta_metro", "w");
	ofp2 = fopen("Chi_vs_Temp_metro", "w");
	ofp3 = fopen("Cv_vs_Temp_metro", "w");

	for (Beta=0.04;Beta<2.5;Beta+=0.002)
	{
		
		srand(time(NULL));

		initialize();             // initialize spins in the 2D lattice

		E_initial=new_config();   // Calculate the initial energy 
		Mag_initial=mag();        // Calculate the initial magnetization

		/*Thermalization*/

		for (tau=0;tau<THSWEEPS;tau++)    
		{
			sweep(Beta);
		}

		/* Calculation of Average Magnetisation (Avg_Mag), Susceptibility (Chi) and Specific Heat (Cv) */

		for (sample_size=0;sample_size<NSAMPLES;sample_size++)
		{
			E_sweep=sweep(Beta);
			Total_E=Total_E+E_sweep;
			Total_Esquared=Total_Esquared+(E_sweep*E_sweep);
			Total_Mag=Total_Mag+fabs(mag());
			Total_Magsquared=Total_Magsquared+(mag()*mag());
		}

		Avg_E=Total_E/NSAMPLES;
		Avg_Esquared=Total_Esquared/NSAMPLES;
		Avg_Mag=Total_Mag/NSAMPLES;
		Avg_Magsquared=Total_Magsquared/NSAMPLES;

		Total_E=0.0;
		Total_Esquared=0.0;
		Total_Mag=0.0;
		Total_Magsquared=0.0;

		Chi = (double)(Beta)*(Avg_Magsquared-(Avg_Mag*Avg_Mag));
                Cv = (Beta*Beta)*((double)(Avg_Esquared-(Avg_E*Avg_E)));

		
                Avg_E=Avg_E/(N*N);
                Avg_Esquared=Avg_Esquared/(N*N*N*N);
		Avg_Mag=Avg_Mag/(N*N);
		Avg_Magsquared=Avg_Magsquared/(N*N*N*N);
		sigma = sqrt((double)(Avg_Magsquared-(Avg_Mag*Avg_Mag)));

		fprintf(ofp1, "%4f   %.4f   %.4f \n",Beta,Avg_Mag,sigma/sqrt(NSAMPLES-1));
		fprintf(ofp2, "%.4f  %.4f\n", 1/Beta, Chi);
		fprintf(ofp3, "%.4f  %.4f\n", 1/Beta, Cv);
	}
	return(0);
}












