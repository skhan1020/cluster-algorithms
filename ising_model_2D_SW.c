/*Ising Model in Two Dimensions -- Wolff Algorithm*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define NSAMPLES 1000 
#define N 10 
#define BLOCKS 100
#define THSWEEPS 100000 
int lattice[N*N];      // 2D lattice with N*N lattice sites



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


/*Function for evaluating Integer Minimum */

int min(int a,int b)
{        
	if(a>b)
		return(b);
	else
		return(a);

}



/*One sweep of the SW algorithm*/

double sweep(float Beta)
{
	double p_bond,y2;
	int i1,j1,k1,l1,n=0,a=0,b=0,flag=0,tot_cluster;
	int xbond[N*N],ybond[N*N],number[N*N],DUMMY[N*N],CHECK[N*N];
	int cluster[N*N][N*N];
	int Nelements[N*N]={0};


	p_bond=1-exp((double)(-2*Beta));                

	for (i1=0;i1<N*N;i1++){

		y2=randomnumber();


		if ( (i1+1)%N!=0 && lattice[i1]*lattice[(i1+1)%(N*N)]==1  && y2 <= p_bond)                                    /* Form Bonds (Right) and Place Markers */
			xbond[i1]=1;
		else if ((i1+1)%N==0  && lattice[i1]*lattice[i1+1-N]==1 && y2 <= p_bond)
			xbond[i1]=1;
		else
			xbond[i1]=0;

	} 

	for(i1=0;i1<N*N;i1++){

		if ( lattice[i1]*lattice[(i1+N)%(N*N)]==1  && randomnumber()<=p_bond )                          /* Form Bonds (Down) and Place Markers */ 
			ybond[i1]=1;
		else
			ybond[i1]=0;
	}

	for(i1=0;i1<N*N;i1++){                                               /*  Invoke the DUMMY and CHECK matrices */
		DUMMY[i1]=i1;CHECK[i1]=DUMMY[i1];
	}


	while (flag != N*N){
		flag=0;
		for(i1=0;i1<N*N;i1++){                                            //  Assign new values to the DUMMY matrix based on the bonds. 
			if(xbond[i1]==1 && (i1+1)%N!=0){                          // Each value of the site represents the cluster number to which it belongs 
				DUMMY[i1]=min(DUMMY[i1],DUMMY[i1+1]);
				DUMMY[i1+1]=DUMMY[i1];
			}



			else if (xbond[i1]==1 && (i1+1)%N==0){
				DUMMY[i1]=min(DUMMY[i1],DUMMY[i1+1-N]);
				DUMMY[i1+1-N]=DUMMY[i1];
			}
		}

		for(i1=0;i1<N*N;i1++){
			if(ybond[i1]==1){
				DUMMY[i1]=min(DUMMY[i1],DUMMY[(i1+N)%(N*N)]);
				DUMMY[(i1+N)%(N*N)]=DUMMY[i1];
			}
		}


		for(i1=0;i1<N*N;i1++){                                         // CHECK matrix is used to determine whether all sites have the correct labels before building the clusters   
			if(CHECK[i1]==DUMMY[i1])
			{
				flag=flag+1;
			}
			else
			{
				CHECK[i1]=DUMMY[i1];
			}
		}



	}



	/* Form Clusters */

	for(i1=0;i1<N*N;i1++){                                           
		if(DUMMY[i1]>=0){
			k1=DUMMY[i1];
			for(j1=i1+1;j1<N*N;j1++){
				if(DUMMY[j1]>=0 && DUMMY[i1]==DUMMY[j1]){
					cluster[n][b]=j1;
					Nelements[n]=Nelements[n]+1;
					b=b+1;
					DUMMY[j1]=-1;
				}
			}
			cluster[n][b]=k1;
			Nelements[n]=Nelements[n]+1;
			b=0;
			n=n+1;
		}
		DUMMY[i1]=-1;
	}

	/*  Flip Spins of Each Cluster */

	for(i1=0;i1<n+1;i1++){                                   
		if(randomnumber()<0.5){
			for(j1=0;j1<Nelements[i1];j1++){
				lattice[cluster[i1][j1]]=-lattice[cluster[i1][j1]];
			}
		}
	} 



	return(mag());
}


int main()
{


	float Beta;   // Inverse of Temperature
	int block,tau,sample_size; // iterators for thermalization and averaging thermodynamic qtys.
	double Mag_initial;   // initial magnetization of the 2D lattice
	double E_initial;    // initial energy of the 2D lattice
	double Avg_Mag,Avg_Magsquared,Avg_Block_Mag,Avg_Block_Magsquared;
	double Total_Mag=0.0,Total_Magsquared=0.0,Total_Block_Mag=0.0,Total_Block_Magsquared=0.0;
	double Avg_E,Avg_Esquared,Avg_Block_E,Avg_Block_Esquared;
	double Total_E=0.0,Total_Esquared=0.0,Total_Block_E=0.0,Total_Block_Esquared=0.0;
	double Chi;   // Susceptiblity
	double Cv;    // Specific Heat
	double sigma;  // standard deviation of the avgerage magnetization per site
	FILE *ofp1, *ofp2, *ofp3;

	ofp1 = fopen("2D_Avg_Magnetization_vs_Beta_SW", "w");
	ofp2 = fopen("2D_Chi_vs_Temp_SW", "w");
	ofp3 = fopen("2D_Cv_vs_Temp_SW", "w");


	for (Beta=0.04;Beta<2.5;Beta+=0.002)
	{
		srand(time(NULL));

		initialize();  // initial configuration of spins in the 2D lattice


		E_initial=new_config();    // initial energy of the 2D lattice  
		Mag_initial=mag();        // initial magnetization of the 2D lattice


		/*Thermalization*/

		for (tau=0;tau<THSWEEPS;tau++)
		{
			sweep(Beta);
		}

		/* Calculation of Average Magnetisation (Avg_Mag) , Susceptibility (Chi) and Specific Heat (Cv) */
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




















