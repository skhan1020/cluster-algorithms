/* Wolff Algorithm for the XY Model in Two Dimensions*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define N 4    // dimensions along each side
#define NSAMPLES 1000    // number of samples for block averaging  
#define THSWEEPS 100000   // number of thermalization steps
#define BLOCKS 500    //  number of blocks
#define Pi 4.0*atan(1.0)   // define tan inverse
double lattice[N*N];    // N by N 2D square lattice
float Beta;    // Inverse of Temperature
int flag = 1;   // Flag value determines whether Avg. Stiffness (flag = 1) or Avg. Energy (flag = 0) is calculated

/* Randomize the spins in the two dimensional lattice*/
void initialize(){

	int i;

	for (i=0; i<N*N; i++)
	{
		lattice[i]=2*Pi;   // spins are pointing along +ve x-axis

	}

}


/*Function for evaluating Minimum */

double min(double a,double b)
{
	if(a>b)
		return(b);
	else
		return(a);

}


/* Function for generating Random Number */

double randomnumber()
{

	return(rand()/(double) (RAND_MAX+1.0));

}


/* Function for calcuating Magnetization*/

double mag()
{
	int i;
	double Mag=0.0;
	for(i=0;i<N*N;i++){

		Mag+=cos(lattice[i]);

	}

	return(Mag/(N*N));

}



/*Function for calculating N_vortices */

double NVORTICES(int a)
{
	double N_vor,Nplus=0.0,Nminus=0.0,N1,N2,N3,N4;
	double dummy_lattice[N][N];
	int l,k;

	for(l=0;l<N;l++){
		for(k=0;k<N;k++){
			if(lattice[k+l*N]>Pi)
				dummy_lattice[k][l]=-(2*Pi-lattice[k+l*N]);
			else
				dummy_lattice[k][l]=lattice[k+l*N];
		}
	}


	for(k=0;k<N;k++){
		for(l=0;l<N;l++){

			N1=(dummy_lattice[k][(l+1)%N]-dummy_lattice[k][l]);
			if(N1>Pi)
				N1=-(2*Pi-N1);
			if(N1<-Pi)
				N1=(2*Pi+N1);

			N2=(dummy_lattice[(k+1)%N][(l+1)%N]-dummy_lattice[k][(l+1)%N]);
			if(N2>Pi)
				N2=-(2*Pi-N2);
			if(N2<-Pi)
				N2=(2*Pi+N2);


			N3=(dummy_lattice[(k+1)%N][l]-dummy_lattice[(k+1)%N][(l+1)%N]);
			if(N3>Pi)
				N3=-(2*Pi-N3);
			if(N3<-Pi)
				N3=(2*Pi+N3);


			N4=(dummy_lattice[k][l]-dummy_lattice[(k+1)%N][l]);
			if(N4>Pi)
				N4=-(2*Pi-N4);
			if(N4<-Pi)
				N4=(2*Pi+N4);


			if(fabs(N1+N2+N3+N4+2*Pi)<0.000001)
				Nplus=Nplus+1.0;
			if(fabs(N1+N2+N3+N4-2*Pi)<0.000001)
				Nminus=Nminus+1.0;

		}
	}


	N_vor=fabs(Nplus-Nminus);
	if (a)
		return (Nplus);
	else
		return(Nminus);
}



/*Function for calculating the energy of a given configuration*/
double new_config()
{
	double dummy_lattice[N][N];
	int i2,j2;
	double E=0.0,Ix=0.0,Iy=0.0,Ixsquared=0.0,Iysquared=0.0;
	double rho;

	for(i2=0;i2<N;i2++){
		for(j2=0;j2<N;j2++){
			dummy_lattice[i2][j2]=lattice[j2+i2*N];
		}
	}

	for (i2=0;i2<N;i2++)
	{
		for (j2=0;j2<N;j2++)
		{
			E += -(1.0/2.0)*(cos(dummy_lattice[i2][j2]-dummy_lattice[(i2+1)%N][j2])+cos(dummy_lattice[i2][j2]-dummy_lattice[(i2-1+N)%N][j2])+cos(dummy_lattice[i2][j2]-dummy_lattice[i2][(j2+1)%N])+cos(dummy_lattice[i2][j2]-dummy_lattice[i2][(j2-1+N)%N]));


			Ix += (1.0/2.0)*(sin(dummy_lattice[i2][j2]-dummy_lattice[(i2+1)%N][j2])+sin(dummy_lattice[i2][j2]-dummy_lattice[(i2-1+N)%N][j2]));


			Iy += (1.0/2.0)*(sin(dummy_lattice[i2][j2]-dummy_lattice[i2][(j2+1)%N])+sin(dummy_lattice[i2][j2]-dummy_lattice[i2][(j2-1+N)%N]));

		}
	}


	rho=(-E-Beta*(Ix*Ix+Iy*Iy))/(2*N*N);

	E=E/(N*N);

	if (flag)
		//return(E);
		return(rho);
	else
		return(E);
}



/* Function for calculating the bond probability between n-n sites */
double p_bond( int p, int q, double r)
{

	return(1-exp(min(0,-2*Beta*cos(lattice[p]-r)*cos(lattice[q]-r))));

}




/*One sweep of the Wolff algorithm*/

double sweep()
{
	double y1,y2,y3,y4,R_angle,R_x,R_y,R_theta;
	int i1,j1,k1,list=0,n=0,a=0,b=0,c=0,flag1=1,flag2=0,p1_element,p2_element,tot_cluster;
	int DUMMY[N*N];
	int *cluster,*perimeter;
	int Nelements=0;

	for(j1=0;j1<N*N;j1++){                                               //  Invoke the DUMMY matrix 
		DUMMY[j1]=j1;
	}


	cluster=malloc(N*N*sizeof(int));perimeter=malloc(N*N*sizeof(int));

	if(cluster==NULL){printf("NO MEMORY CLUSTER! \n"); exit(0);}
	if(perimeter==NULL){printf("NO MEMORY PERIMETER! \n");exit(0);}


	i1=(N*N)*randomnumber();                                             // Place the seed of a cluster 
	cluster[0]=i1;
	perimeter[0]=i1;
	k1=0;

	R_angle=(2*Pi)*randomnumber();                                       // Angle made by r-vector (perp to the hyperplane) w.r.t. spin at i1

	while (flag1 !=0 ){




		i1=cluster[k1];


		y2=randomnumber();
		y1=randomnumber();

		y3=randomnumber();
		y4=randomnumber();


		if ( (i1+1)%N!=0 && y2<=p_bond(i1,(i1+1)%(N*N),R_angle) && DUMMY[(i1+1)%(N*N)]!=-1){                                    // right bonds 
			j1=(i1+1)%(N*N);
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
			a=a+1;
			Nelements=Nelements+1;
		}
		else if ((i1+1)%N==0  && y2<=p_bond(i1,i1+1-N,R_angle) && DUMMY[i1+1-N]!=-1){
			j1=i1+1-N;
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
			a=a+1;
			Nelements=Nelements+1;
		}


		if (i1%N!=0 && y1<=p_bond(i1,(i1-1)%(N*N),R_angle) && DUMMY[(i1-1)%(N*N)]!=-1)           //  left bonds 
		{
			j1=(i1-1)%(N*N);
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
			a=a+1;
			Nelements=Nelements+1;
		}
		else if (i1%N==0 && y1<=p_bond(i1,i1-1+N,R_angle) && DUMMY[i1-1+N]!=-1){
			j1=i1-1+N;
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
			a=a+1;
			Nelements=Nelements+1;
		}



		if ( y3<=p_bond(i1,(i1+N)%(N*N),R_angle) && DUMMY[(i1+N)%(N*N)]!=-1){                        //   down bonds  
			j1=(i1+N)%(N*N);
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
			a=a+1;
			Nelements=Nelements+1;
		}



		if ( y4<=p_bond(i1,(i1-N+N*N)%(N*N),R_angle) && DUMMY[(i1-N+N*N)%(N*N)]!=-1){                         //  up bonds                
			j1=(i1-N+N*N)%(N*N);
			DUMMY[j1]=-1;
			b=b+1;
			cluster[b]=j1;                                           //Add neighbour to cluster
			perimeter[a]=j1;                                        // Add neighbour to perimeter
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

	Nelements=Nelements+1;b=b+1;                    // Original seed has to be included 


	for(j1=0;j1<b;j1++){                                         // Flip all spins of cluster 
		R_x=cos(lattice[cluster[j1]])-2*cos(R_angle)*cos(lattice[cluster[j1]]-R_angle);
		R_y=sin(lattice[cluster[j1]])-2*sin(R_angle)*cos(lattice[cluster[j1]]-R_angle);
		if(R_x >=0.0 && R_y>=0.0)
			R_theta=atan(R_y/R_x);
		if(R_x <0.0 && R_y >=0.0)
			R_theta=Pi-atan(fabs(R_y/R_x));
		if(R_x<0.0 && R_y<0.0)
			R_theta=Pi+atan(R_y/R_x);
		if(R_x>=0.0 && R_y<0.0)
			R_theta=2*Pi-atan(fabs(R_y/R_x));

		lattice[cluster[j1]]=fmod(R_theta,2*Pi);
		lattice[cluster[j1]]=R_theta;
	}



	free(perimeter);
	free(cluster);
	b=0;Nelements=0;
	c=0;
	flag1=1;

	return(new_config());
}




int main()
{


	float TEMP;    // Temperature
	int block,tau,sample_size;   // iterators for thermalization and averaging process
	double rho_initial;   // Initial spin stiffness
	double E_vor;   // energy of a single vortex
	double Avg_rho,Avg_rhosquared,Avg_Block_rho,Avg_Block_rhosquared;
	double Total_rho=0.0,Total_Block_rho=0.0;
	double Total_rhosquared=0.0,Total_Block_rhosquared=0.0;
	double Avg_Nplus=0.0, Avg_Nminus=0.0;
	double Avg_Mag=0.0;  // Average Magnetization
	double Cv;   // Specific Heat
	double sigma;  // standard deviation for rho/energy
	FILE *ofp1,*ofp2,*ofp3,*ofp4,*ofp5;


	if (flag)
		ofp1=fopen("2D_xy_Wolff_rho","w");
	else
		ofp1=fopen("2D_xy_Wolff_energy","w");

	ofp2=fopen("2D_xy_Wolff_mag", "w");
	ofp3=fopen("2D_xy_Wolff_cv","w");
	ofp4=fopen("2D_xy_Wolff_Nplus","w");
	ofp5=fopen("2D_xy_Wolff_Nminus", "w");



	for (TEMP=0.04;TEMP<4.5;TEMP+=0.002)
	{Beta=1/TEMP;


		
		srand(time(NULL));

		initialize();      // Initialize the spins in the 2D lattice.

		E_vor=(Pi-2*TEMP)*log(N);   // Vortex energy


		/*Thermalization and Calculation of Average Energy/Average Spin Stiffness : rho and Average Magnetization*/


		for (tau=0;tau<THSWEEPS;tau++)
		{
			sweep();
		}

		for(block=0;block<BLOCKS;block++){

			for (sample_size=0;sample_size<NSAMPLES;sample_size++)
			{
				rho_initial=sweep();
				Total_Block_rho += rho_initial;
				Total_Block_rhosquared += rho_initial*rho_initial;
			}

			Avg_Block_rho=(Total_Block_rho/(NSAMPLES));
			Avg_Block_rhosquared=(Total_Block_rhosquared/(NSAMPLES));


			Total_Block_rho=0.0;
			Total_Block_rhosquared=0.0;


			Total_rho += Avg_Block_rho;
			Total_rhosquared += Avg_Block_rhosquared;


			Avg_Nplus += NVORTICES(1);
			Avg_Nminus += NVORTICES(0);

			Avg_Mag += mag();


		}

		Avg_rho=Total_rho/BLOCKS;
		Avg_rhosquared=Total_rhosquared/BLOCKS;


		Total_rho=0.0;
		Total_rhosquared=0.0;

		sigma=sqrt((Avg_rhosquared-(Avg_rho*Avg_rho))/(BLOCKS-1));

		Avg_Nplus = Avg_Nplus/BLOCKS;
		Avg_Nminus = Avg_Nminus/BLOCKS;



		fprintf(ofp1, "%4f   %.4f   %.4f \n",TEMP,Avg_rho,sigma);

		fprintf(ofp2, "%4f   %.4f   %.4f \n",TEMP,Avg_Mag/(BLOCKS));

		fprintf(ofp4, "%.4f  %.4f\n",TEMP,Avg_Nplus);
		fprintf(ofp5, "%.4f  %.4f\n",TEMP,Avg_Nminus);


		if(flag == 0)
		{
			Cv = (Beta*Beta)*(Avg_rhosquared-(Avg_rho*Avg_rho));
			fprintf(ofp3, "%.4f  %.4f\n",TEMP,Cv);
		}


	}

	return(0);
}















