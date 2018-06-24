/*Ising Model in Three Dimensions - Wolff Algorithm*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define NSAMPLES 1000    // number of samples for averaging thermodynamic qtys 
#define THSWEEPS 100000  // number of thermalization steps
#define N 16             // Dimensions of the cubic lattice
int lattice[N*N*N];      // N by N by N lattice




/* Function for generating Random Number */

double randomnumber()
{
	return((double) rand()/(double) (RAND_MAX+1.0));
}




/* Initialize the spins in the two dimensional lattice*/
void initialize(){

	int i;

	for (i=0; i<N*N*N; i++)
	{
		lattice[i]=1;

	}

}




/* Function for calcuating Magnetization*/
double mag()
{
	int i,spin_up=0,spin_down=0;

	for(i=0;i<N*N*N;i++)
	{
		if (lattice[i]==1)
			spin_up++;
		else
			spin_down++;
	}

	return((double)(spin_up-spin_down));

}




/*Function for calculating the energy of a given configuration*/
double new_config()
{
	int dummy_lattice[N][N][N];
	int i, j, k;
	double E=0.0;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			for(k=0; k<N; k++){
				dummy_lattice[i][j][k] = lattice[k+j*N+i*N*N];
			}
		}
	}
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++)
		{
			for(k=0; k<N; k++)
			{
				E += -(1.0/2.0)*(dummy_lattice[i][j][k]*(dummy_lattice[(i+1)%N][j][k]+dummy_lattice[(i+N-1)%N][j][k]+dummy_lattice[i][(j+1)%N][k]+dummy_lattice[i][(j+N-1)%N][k]+dummy_lattice[i][j][(k+1)%N]+dummy_lattice[i][j][(k+N-1)%N]));
			}
		}
	}
	return(E);
}








/*One sweep of the Wolff algorithm */

double sweep(float Beta)
{
	double p_bond,y1,y2,y3,y4,y5;
	int i1,j1,k1,n=0,a=0,b=0,c=0,flag1=1;
	short *x1bond,*x2bond,*y1bond,*y2bond,*z1bond,*z2bond;
	short DUMMY[N*N*N];
	short *cluster,*perimeter;
	int Nelements=0;

	p_bond=1.0-exp((double)(-2*Beta));                






	for(j1=0;j1<N*N*N;j1++){                                              //   Invoke the DUMMY matrix 
		DUMMY[j1]=j1;
	}


	cluster=malloc(N*N*N*sizeof(short));perimeter=malloc(N*N*N*sizeof(short));
	x1bond=malloc(N*N*N*sizeof(short));x2bond=malloc(N*N*N*sizeof(short));y1bond=malloc(N*N*N*sizeof(short));y2bond=malloc(N*N*N*sizeof(short));z1bond=malloc(N*N*N*sizeof(short));z2bond=malloc(N*N*N*sizeof(short));

	i1=(N*N*N)*randomnumber();                                             // Place the seed of a cluster 

	cluster[0]=i1;
	perimeter[0]=i1;
	k1=0;

	while (flag1 !=0 ){


		i1=cluster[k1];


		y2=randomnumber();
		y1=randomnumber();

		y3=randomnumber();
		y4=randomnumber();

		y5=randomnumber();



		if ( (i1+1)%N!=0 && lattice[i1]*lattice[(i1+1)%(N*N*N)]==1 && y2<=p_bond && DUMMY[(i1+1)%(N*N*N)]!=-1){                                    // right bonds 
			DUMMY[(i1+1)%(N*N*N)]=-1;
			x1bond[i1]=1;
		}
		else if ((i1+1)%N==0  && lattice[i1]*lattice[i1+1-N]==1 &&  y2<=p_bond && DUMMY[i1+1-N]!=-1){
			DUMMY[i1+1-N]=-1;
			x1bond[i1]=1;
		}
		else
			x1bond[i1]=0;


		if (i1%N!=0 && lattice[i1]*lattice[(i1-1)%(N*N*N)]==1  && y1<=p_bond && DUMMY[(i1-1)%(N*N*N)]!=-1)           //  left bonds 
		{
			DUMMY[(i1-1)%(N*N*N)]=-1;
			x2bond[i1]=1;
		}
		else if (i1%N==0 && lattice[i1]*lattice[i1-1+N]==1 && y1<=p_bond && DUMMY[i1-1+N]!=-1){
			DUMMY[i1-1+N]=-1;
			x2bond[i1]=1;
		}
		else
			x2bond[i1]=0;



		if ( (i1+N)%(N*N)>N-1 &&  lattice[i1]*lattice[(i1+N)]==1 && y3<=p_bond && DUMMY[i1+N]!=-1 ){                        //  y  down bonds  
			DUMMY[i1+N]=-1;
			y1bond[i1]=1;
		}
		else if ( (i1+N)%(N*N)<=N-1 &&  lattice[i1]*lattice[i1-(N-1)*N]==1 && y3<=p_bond && DUMMY[i1-N*(N-1)]!=-1 ){                        
			DUMMY[i1-N*(N-1)]=-1;
			y1bond[i1]=1;
		}
		else
			y1bond[i1]=0;




		if ( (i1-N+N*N)%(N*N)<(N-1)*N &&  lattice[i1]*lattice[i1-N]==1 && y4<=p_bond && DUMMY[i1-N]!=-1 ){                         // y  up bonds                
			DUMMY[i1-N]=-1;
			y2bond[i1]=1;
		}
		else if ( (i1-N+N*N)%(N*N)>=(N-1)*N &&  lattice[i1]*lattice[i1+(N-1)*N]==1  && y4<=p_bond && DUMMY[i1+N*(N-1)]!=-1 ){                                    
			DUMMY[i1+N*(N-1)]=-1;
			y2bond[i1]=1;
		}
		else
			y2bond[i1]=0;



		if (  lattice[i1]*lattice[(i1+N*N)%(N*N*N)]==1 && y5<=p_bond && DUMMY[(i1+N*N)%(N*N*N)]!=-1){                                    // z down bonds   
			DUMMY[(i1+N*N)%(N*N*N)]=-1;
			z1bond[i1]=1;
		}
		else 
			z1bond[i1]=0;


		if (  lattice[i1]*lattice[(i1-N*N+N*N*N)%(N*N*N)]==1 && y5<=p_bond && DUMMY[(i1-N*N+N*N*N)%(N*N*N)]!=-1){                     // z up bonds   
			DUMMY[(i1-N*N+N*N*N)%(N*N*N)]=-1;
			z2bond[i1]=1;
		}
		else 
			z2bond[i1]=0;




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

		if(y1bond[i1]==1 && (i1+N)%(N*N)>N-1) 
		{
			b=b+1;
			cluster[b]=i1+N;
			perimeter[a]=i1+N;
			a=a+1;
			Nelements=Nelements+1;
		}
		else if(y1bond[i1]==1 && (i1+N)%(N*N)<=N-1)
		{
			b=b+1;
			cluster[b]=i1-N*(N-1);
			perimeter[a]=i1-N*(N-1);
			a=a+1;
			Nelements=Nelements+1;
		}

		if(y2bond[i1]==1 && (i1-N+N*N)%(N*N)<(N-1)*N)
		{
			b=b+1;
			cluster[b]=i1-N;
			perimeter[a]=i1-N;
			a=a+1;
			Nelements=Nelements+1;
		}
		if(y2bond[i1]==1 && (i1-N+N*N)%(N*N)>=(N-1)*N)
		{
			b=b+1;
			cluster[b]=i1+(N-1)*N;
			perimeter[a]=i1+(N-1)*N;
			a=a+1;
			Nelements=Nelements+1;
		}

		if(z1bond[i1]==1)
		{
			b=b+1;
			cluster[b]=(i1+N*N)%(N*N*N);
			perimeter[a]=(i1+N*N)%(N*N*N);
			a=a+1;
			Nelements=Nelements+1;
		}
		if(z2bond[i1]==1)
		{
			b=b+1;
			cluster[b]=(i1-N*N+N*N*N)%(N*N*N);
			perimeter[a]=(i1-N*N+N*N*N)%(N*N*N);
			a=a+1;
			Nelements=Nelements+1;
		}



		DUMMY[i1]=-1;




		if(k1==Nelements)
			flag1=0;


		c=c+a;
		k1=k1+1;
		a=0;


	}

	Nelements=Nelements+1;b=b+1;                    /* Original seed has to be included */


	for(j1=0;j1<b;j1++){              // Flip all spins of cluster
		lattice[cluster[j1]]=-lattice[cluster[j1]];
	}

	free(x1bond);free(x2bond);
	free(y1bond);free(y2bond);
	free(z1bond);free(z2bond);
	free(perimeter);
	free(cluster);
	b=0;Nelements=0;
	c=0;
	flag1=1;


	return(mag());
}





int main()
{

	float Beta;      // Inverse of Temperature
	int tau,sample_size;   // iterators for thermalization and averaging process
	double Mag_initial;      // initial magnetization of the 3D lattice
	float E_initial;         // initial energy of the 3D lattice
	double Avg_Mag,Avg_Magsquared;
	double Total_Mag=0.0;
	double Total_Magsquared=0.0;
	double Avg_E,Avg_Esquared;
	double Total_E=0.0,Total_Esquared=0.0;
	double Cv;    // Specific Heat
	double Chi;   // Susceptibility
	double sigma;   // standard deviation for average magnetization per site
	FILE *ofp1, *ofp2, *ofp3;   

	ofp1 = fopen("3D_Avg_Magnetization_vs_Beta_Wolff", "w");
	ofp2 = fopen("3D_Chi_vs_Temp_Wolff", "w");
	ofp3 = fopen("3D_Cv_vs_Temp_Wolff", "w");




	for (Beta=0.04;Beta<2.5;Beta+=0.002)
	{

		/* Initialize spins in the one dimensional lattice and calculate initial energy and magnetisation*/
		srand(time(NULL));


		initialize();

		E_initial=new_config();
		Mag_initial=mag();

		/*Thermalization and Calculation of Average Energy and Average Magnetisation*/

		for (tau=0;tau<THSWEEPS;tau++)
		{
			sweep(Beta);
		}


		for (sample_size=0;sample_size<NSAMPLES;sample_size++)
		{
			Mag_initial = fabs(sweep(Beta));
			Total_Mag += Mag_initial;
			Total_Magsquared += Mag_initial*Mag_initial;


			E_initial = new_config();
			Total_E += E_initial;
			Total_Esquared += E_initial*E_initial;



		}
		Avg_Mag=(Total_Mag/(NSAMPLES));
		Avg_Magsquared=(Total_Magsquared/(NSAMPLES));


                Avg_E=(Total_E/(NSAMPLES));
                Avg_Esquared=(Total_Esquared/(NSAMPLES));

                Chi = (double)(Beta)*(Avg_Magsquared-(Avg_Mag*Avg_Mag));
                Cv = (Beta*Beta)*((double)(Avg_Esquared-(Avg_E*Avg_E)));


		Total_Mag=0.0;
		Total_Magsquared=0.0;

                Total_E=0.0;
                Total_Esquared=0.0;


		Avg_Mag=Avg_Mag/(N*N*N);
		Avg_Magsquared=Avg_Magsquared/(N*N*N*N*N*N);
		sigma = sqrt((Avg_Magsquared-Avg_Mag*Avg_Mag));


		fprintf(ofp1, "%4f   %.4f   %.4f \n",Beta,Avg_Mag,sigma/sqrt(NSAMPLES-1));
		fprintf(ofp2, "%.4f  %.4f\n", 1/Beta, Chi);
		fprintf(ofp3, "%.4f  %.4f\n", 1/Beta, Cv);



	}

	return(0);
}


















