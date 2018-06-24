/*Ising Model in Two Dimensions -- Wolff Algorithm*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define NSAMPLES 1000 
#define N 10 
#define BLOCKS 100
#define THSWEEPS 100000 

int lattice[N*N];  // 2D lattice with N*N lattice sites


/* Function for generating Random Number */

double randomnumber()
{
	return((double) rand()/(double) (RAND_MAX+1.0));
}



/* Randomize the spins in the two dimensional lattice*/
void initialize(){

	int i;

	for (i=0; i<N*N; i++)
	{
		if (randomnumber()>0.5)
			lattice[i]=1;
		else
			lattice[i]=-1;

	}

}


/* Function for calcuating Magnetization*/
double mag()
{
	int i,spin_up=0,spin_down=0;

	for(i=0;i<N*N;i++)
	{
		if (lattice[i]==1)
			spin_up++;
		else
			spin_down++;
	}

	return((double)(spin_up-spin_down));

}


/*Function for calculating the energy of a given configuration*/
float new_config()
{
	int i,j;
	float E=0.0;
	for (i=0;i<N*N;i++)
	{
		if((i-1)%N!=0 && i%N!=0)
			E=E-(1.0/2.0)*(lattice[i]*(lattice[(i+1)%(N*N)]+lattice[(i-1)%(N*N)]+lattice[(i-N+N*N)%(N*N)]+lattice[(i+N)%(N*N)]));
		else if( (i+1)%N==0 )
			E=E-(1.0/2.0)*(lattice[i]*(lattice[(i-1)%(N*N)]+lattice[(i-N+1)%(N*N)]+lattice[(i-N+N*N)%(N*N)]+lattice[(i+N)%(N*N)]));
		else if (i%N==0)
			E=E-(1.0/2.0)*(lattice[i]*(lattice[(i+N-1)%(N*N)]+lattice[(i+1)%(N*N)]+lattice[(i-N+N*N)%(N*N)]+lattice[(i+N)%(N*N)]));
	}
	return(E);
}




/*One sweep of the Wolff algorithm*/

double sweep(float Beta)
{
	double p_bond,y1,y2,y3,y4;
	int i1,j1,k1,n=0,a=0,b=0,c=0,flag1=1;
	int *x1bond,*x2bond,*y1bond,*y2bond;
	int DUMMY[N*N];
	int *cluster,*perimeter;
	int Nelements=0;

	p_bond=1.0-exp((double)(-2*Beta));                






	for(j1=0;j1<N*N;j1++){                                               /*  Invoke the DUMMY matrix */
		DUMMY[j1]=j1;
	}



	cluster=malloc(N*N*sizeof(int));perimeter=malloc(N*N*sizeof(int));
	x1bond=malloc(N*N*sizeof(int));x2bond=malloc(N*N*sizeof(int));y1bond=malloc(N*N*sizeof(int));y2bond=malloc(N*N*sizeof(int));

	i1=(N*N)*randomnumber();                                             /* Place the seed of a cluster */

	cluster[0]=i1;
	perimeter[0]=i1;
	k1=0;

	while (flag1 !=0 ){


		i1=cluster[k1];


		y2=randomnumber();
		y1=randomnumber();

		y3=randomnumber();
		y4=randomnumber();



		if ( (i1+1)%N!=0 && lattice[i1]*lattice[(i1+1)%(N*N)]==1 && DUMMY[(i1+1)%(N*N)] != -1 && y2 <= p_bond){         /* Form right bonds and assign -1 to DUMMY matrix */
			DUMMY[(i1+1)%(N*N)]=-1;
			x1bond[i1]=1;
		}
		else if ((i1+1)%N==0  && lattice[i1]*lattice[i1+1-N]==1 && DUMMY[i1+1-N] != -1 && y2 <= p_bond){
			DUMMY[i1+1-N]=-1;
			x1bond[i1]=1;
		}
		else 
			x1bond[i1]=0;


		if (i1%N!=0 && lattice[i1]*lattice[(i1-1)%(N*N)]==1 && DUMMY[(i1-1)%(N*N)] != -1 && y1 <= p_bond)           /*  Form left bonds and assign -1 to DUMMY matrix*/
		{
			DUMMY[(i1-1)%(N*N)]=-1;
			x2bond[i1]=1;
		}
		else if (i1%N==0 && lattice[i1]*lattice[i1-1+N]==1 && DUMMY[i1-1+N] != -1 && y1 <= p_bond){
			DUMMY[i1-1+N]=-1;
			x2bond[i1]=1;
		}
		else
			x2bond[i1]=0;


		if ( lattice[i1]*lattice[(i1+N)%(N*N)]==1 && DUMMY[(i1+N)%(N*N)] != -1 && y3 <= p_bond){                        /*  Form  down bonds and assign -1 to DUMMY matrix */
			DUMMY[(i1+N)%(N*N)]=-1;
			y1bond[i1]=1;
		}
		else
			y1bond[i1]=0;

		if ( lattice[i1]*lattice[(i1-N+N*N)%(N*N)]==1 && DUMMY[(i1-N+N*N)%(N*N)] != -1 && y4 <= p_bond){                         /* Form up bonds and assign -1 to DUMMY matrix  */             
			DUMMY[(i1-N+N*N)%(N*N)]=-1;
			y2bond[i1]=1;
		}
		else
			y2bond[i1]=0;


		/*  Build the Cluster */

		if(x1bond[i1]==1 && (i1+1)%N!=0)
		{
			b=b+1;
			cluster[b]=i1+1;
			perimeter[a]=i1+1;
			a=a+1;
			Nelements=Nelements+1;

		}
		else if (x1bond[i1]==1 && (i1+1)%N==0)
		{
			b=b+1;
			cluster[b]=i1+1-N;
			perimeter[a]=i1+1-N;
			a=a+1;
			Nelements=Nelements+1;
		}


		if (x2bond[i1]==1 && i1%N!=0)
		{
			b=b+1;
			cluster[b]=i1-1;
			perimeter[a]=i1-1;
			a=a+1;
			Nelements=Nelements+1;
		}
		else if (x2bond[i1]==1 && i1%N==0)
		{
			b=b+1;
			cluster[b]=i1-1+N;
			perimeter[a]=i1-1+N;
			a=a+1;
			Nelements=Nelements+1;
		}
		if(y1bond[i1]==1) 
		{
			b=b+1;
			cluster[b]=(i1+N)%(N*N);
			perimeter[a]=(i1+N)%(N*N);
			a=a+1;
			Nelements=Nelements+1;
		}
		if(y2bond[i1]==1)
		{
			b=b+1;
			cluster[b]=(i1-N+N*N)%(N*N);
			perimeter[a]=(i1-N+N*N)%(N*N);
			a=a+1;
			Nelements=Nelements+1;
		}
		DUMMY[i1]=-1;



		if(k1==Nelements)
			flag1=0;    // flag checks whether all elements have been included in the cluster 


		c=c+a;
		k1=k1+1;
		a=0;




	}

	Nelements=Nelements+1;b=b+1;                    // Original seed has to be included





	/* Flip all spins of cluster */

	for(j1=0;j1<b;j1++){              
		lattice[cluster[j1]]=-lattice[cluster[j1]];
	} 





	free(x1bond);free(x2bond);
	free(y1bond);free(y2bond);
	free(perimeter);
	free(cluster);

	return(mag());
}




int main()
{

	float Beta;      // Inverse of Temperature
	int block,tau,sample_size;   // iterators for thermalization and averaging process
	double Mag_initial;      // initial magnetization of the 2D lattice
	float E_initial;         // initial energy of the 2D lattice
	double Avg_Mag,Avg_Magsquared,Avg_Block_Mag,Avg_Block_Magsquared;
	double Total_Mag=0.0,Total_Block_Mag=0.0;
	double Total_Magsquared=0.0,Total_Block_Magsquared=0.0;
	double Avg_E,Avg_Esquared,Avg_Block_E,Avg_Block_Esquared;
	double Total_E=0.0,Total_Esquared=0.0,Total_Block_E=0.0,Total_Block_Esquared=0.0;
	double Cv;    // Specific Heat
	double Chi;   // Susceptibility
	double sigma;   // standard deviation for average magnetization per site
	FILE *ofp1, *ofp2, *ofp3;

	ofp1 = fopen("2D_Avg_Magnetization_vs_Beta_Wolff", "w");
	ofp2 = fopen("2D_Chi_vs_Temp_Wolff", "w");
	ofp3 = fopen("2D_Cv_vs_Temp_Wolff", "w");



	for (Beta=0.04;Beta<2.5;Beta+=0.002)
	{


		srand(time(NULL));

		initialize();   // initial configuration of spins in the 2D lattice

		E_initial=new_config();       // initial energy of the 2D lattice
		Mag_initial=mag();        // initial mag. of 2D lattice




		/*Thermalization*/

		for (tau=0;tau<THSWEEPS;tau++)
		{
			sweep(Beta);
		}



		/*Calculation of Average Energy and Average Magnetization */ 
		for(block=0;block<BLOCKS;block++){

			for (sample_size=0;sample_size<NSAMPLES;sample_size++)
			{
				Mag_initial=fabs(sweep(Beta));
				Total_Block_Mag=Total_Block_Mag+Mag_initial;
				Total_Block_Magsquared=Total_Block_Magsquared+(Mag_initial*Mag_initial);

				E_initial = new_config();
				Total_Block_E=Total_Block_E+E_initial;
				Total_Block_Esquared=Total_Block_Esquared+(E_initial*E_initial);

			}
			Avg_Block_Mag=(Total_Block_Mag/(NSAMPLES));
			Avg_Block_Magsquared=(Total_Block_Magsquared/(NSAMPLES));

			Avg_Block_E=(Total_Block_E/(NSAMPLES));
			Avg_Block_Esquared=(Total_Block_Esquared/(NSAMPLES));


			Total_Block_Mag=0.0;
			Total_Block_Magsquared=0.0;


			Total_Block_E=0.0;
			Total_Block_Esquared=0.0;


			Total_Mag=Total_Mag+Avg_Block_Mag;
			Total_Magsquared=Total_Magsquared+Avg_Block_Magsquared;

			Total_E=Total_E+Avg_Block_E;
			Total_Esquared=Total_Esquared+Avg_Block_Esquared;



		}
		Avg_Mag=Total_Mag/BLOCKS;
		Avg_Magsquared=Total_Magsquared/BLOCKS;

		Avg_E=Total_E/BLOCKS;
		Avg_Esquared=Total_Esquared/BLOCKS;


		Total_Mag=0.0;
		Total_Magsquared=0.0;

		Total_E=0.0;
		Total_Esquared=0.0;


		Chi = (double)(Beta)*(Avg_Magsquared-(Avg_Mag*Avg_Mag));
		Cv = (Beta*Beta)*((double)(Avg_Esquared-(Avg_E*Avg_E)));


		Avg_Mag=Avg_Mag/(N*N);
		Avg_Magsquared=Avg_Magsquared/(N*N*N*N);


		sigma = sqrt((double)(Avg_Magsquared-(Avg_Mag*Avg_Mag)));

		fprintf(ofp1, "%4f   %.4f   %.4f \n",Beta,Avg_Mag,sigma/sqrt(NSAMPLES-1));
		fprintf(ofp2, "%.4f  %.4f\n", 1/Beta, Chi);
		fprintf(ofp3, "%.4f  %.4f\n", 1/Beta, Cv);



	}
	return(0);
}

















