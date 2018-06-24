/*Ising Model in Pyrochole Lattice*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define NSAMPLES 1000               // number of samples for averaging 
#define N 3                         // Size of lattice
#define THSWEEPS 100000             // number of thermalization steps
int lattice[4*N*N*N];               // 4N by N by N lattice
short *cluster,*perimeter;          // arrays for cluster and perimeter lattice sites



/* Function for generating Random Number */

double randomnumber()
{
	return( rand()/(double) (RAND_MAX+1.0) );
}



/*Initialize the spins*/

void initialize(){

	int i;

	for (i=0; i<4*N*N*N; i++)
	{
		lattice[i]=1;
	}


}



/* Function for calculating Magnetization*/
double mag()
{
	int l,m;
	int spin_up=0,spin_down=0;

	for(l=0;l<4*N*N*N;l++)
	{
		if (lattice[l]==1)
		{
			spin_up ++;
		}
		else
		{
			spin_down ++;
		}
	}

	return((double)(spin_up-spin_down));
}





/*One sweep of the Wolff algorithm*/

double sweep(Beta)
{

	int dummyA[N][N][N],dummyB[N][N][N],dummyC[N][N][N],dummyD[N][N][N];
	double p_bond;
	int i1,j1,k1,i2,j2,k2,n=0,a=0,b=0,c=0,flag1=1;
	short DUMMY[4*N*N*N];
	int Nelements=0;

	p_bond=1.0-exp((double)(-2*Beta));                



	for(j1=0;j1<4*N*N*N;j1++){                                               //  Invoke the DUMMY matrix 
		DUMMY[j1]=j1;
	}


	i1=(4*N*N*N)*randomnumber();                                             // Place the seed of a cluster 
	cluster[0]=i1;
	perimeter[0]=i1;
	k1=0;


	for(i2=0;i2<N;i2++){c=0;
		for(j2=0;j2<N;j2++){
			for(k2=0;k2<N;k2++){
				dummyA[i2][j2][k2]=lattice[k2+j2*(4*N)+i2*(4*N*N)+c];
				dummyB[i2][j2][k2]=lattice[k2+j2*(4*N)+i2*(4*N*N)+1+c];
				dummyC[i2][j2][k2]=lattice[k2+j2*(4*N)+i2*(4*N*N)+2+c];
				dummyD[i2][j2][k2]=lattice[k2+j2*(4*N)+i2*(4*N*N)+3+c];

				c=c+3;
			}
			c=0;
		}
	}



	while (flag1 !=0 ){


		i1=cluster[k1];

		if (i1%4==0)
		{
			i2=i1/(4*N*N);
			j2=(i1-i2*(4*N*N))/(4*N);
			k2=(i1-i2*(4*N*N)-j2*(4*N)+1)/4;


			if(dummyA[i2][j2][k2]*dummyB[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+1]!=-1){
				//bond formed ---- include in the cluster
				j1=i1+1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(dummyA[i2][j2][k2]*dummyC[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+2]!=-1){
				//bond formed --- include in the cluster
				j1=i1+2;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(dummyA[i2][j2][k2]*dummyD[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+3]!=-1){
				//bond formed --- include in the cluster
				j1=i1+3;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}


			if(k2!=0){
				if(dummyA[i2][j2][k2]*dummyB[i2][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1-3]!=-1){
					//bond formed ---include in the cluster
					j1=i1-3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else
			{
				if(dummyA[i2][j2][k2]*dummyB[i2][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N-3]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N-3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			if(i2!=0 && k2!=N-1){
				if(dummyA[i2][j2][k2]*dummyD[i2-1][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N+7]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N+7;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2!=N-1)
			{
				if(dummyA[i2][j2][k2]*dummyD[(i2-1+N)%N][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)+7]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)+7;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=0 && k2==N-1)
			{
				if(dummyA[i2][j2][k2]*dummyD[i2-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N-4*(N-1)+3]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N-4*(N-1)+3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2==N-1)
			{
				if(dummyA[i2][j2][k2]*dummyD[i2+N-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)-4*(N-1)+3]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)-4*(N-1)+3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}



			if(j2!=0){
				if(dummyA[i2][j2][k2]*dummyC[i2][(j2+N-1)%N][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N+2]!=-1){
					//bond formed ---include in the cluster
					j1=i1-4*N+2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else
			{
				if(dummyA[i2][j2][k2]*dummyC[i2][(j2+N-1)%N][k2]==1 && randomnumber()<=p_bond && DUMMY[(i1+8*N+2)]!=-1){
					//bond formed ---include in the cluster
					j1=i1+8*N+2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}


		}

		if((i1-1)%4==0)  // i1 belongs to B sublattice
		{

			i2=i1/(4*N*N);
			j2=(i1-i2*(4*N*N))/(4*N);
			k2=(i1-i2*(4*N*N)-j2*(4*N)+1)/4;

			if(dummyB[i2][j2][k2]*dummyA[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-1]!=-1){
				//bond formed --- include in the cluster
				j1=i1-1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(k2!=N-1){
				if(dummyB[i2][j2][k2]*dummyA[i2][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1+3]!=-1){
					//bond formed--- include in the cluster
					j1=i1+3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(k2==N-1)
			{
				if(dummyB[i2][j2][k2]*dummyA[i2][j2][k2+1-N]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N+3]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N+3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			if(dummyB[i2][j2][k2]*dummyC[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+1]!=-1){
				//bond formed -- include in the cluster
				j1=i1+1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}


			if(j2!=0 && k2!=N-1){
				if (dummyB[i2][j2][k2]*dummyC[i2][(j2-1+N)%N][(k2+1)%N]==1 && randomnumber() <=p_bond && DUMMY[(i1-4*N+5)]!=-1){
					// bond formed --- include in the cluster
					j1=(i1-4*N+5);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(j2==0 && k2!=N-1)
			{
				if (dummyB[i2][j2][k2]*dummyC[i2][(j2-1+N)%N][(k2+1)%N]==1 && randomnumber() <=p_bond && DUMMY[(i1+8*N+5)]!=-1){
					// bond formed --- include in the cluster
					j1=i1+8*N+5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(j2==0 && k2==N-1)
			{
				if (dummyB[i2][j2][k2]*dummyC[i2][(j2-1+N)%N][(k2+1)%N]==1 && randomnumber() <=p_bond && DUMMY[(i1+4*N+5)]!=-1){
					// bond formed --- include in the cluster
					j1=i1+4*N+5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if (j2!=0 && k2==N-1)
			{
				if (dummyB[i2][j2][k2]*dummyC[i2][(j2-1+N)%N][(k2+1)%N]==1 && randomnumber() <=p_bond && DUMMY[(i1-4*N-4*(N-1)+1)]!=-1){
					// bond formed --- include in the cluster
					j1=i1-4*N-4*(N-1)+1;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			if(dummyB[i2][j2][k2]*dummyD[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+2]!=-1){
				//bond formed -- include in the cluster
				j1=i1+2;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(i2!=0 && k2!=N-1){
				if(dummyB[i2][j2][k2]*dummyD[i2-1][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N+6]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N+6;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2!=N-1)
			{
				if(dummyB[i2][j2][k2]*dummyD[(i2-1+N)%N][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)+6]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)+6;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=0 && k2==N-1)
			{
				if(dummyB[i2][j2][k2]*dummyD[i2-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N-4*(N-1)+2]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N-4*(N-1)+2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2==N-1)
			{
				if(dummyB[i2][j2][k2]*dummyD[i2+N-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)-4*(N-1)+2]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)-4*(N-1)+2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}


		}


		if((i1+2)%4==0)   // i1 belongs to C sublattice
		{
			i2=i1/(4*N*N);
			j2=(i1-i2*(4*N*N))/(4*N);
			k2=(i1-i2*(4*N*N)-j2*(4*N)+1)/4;


			if(dummyC[i2][j2][k2]*dummyA[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-2]!=-1){
				//bond formed --- include in the cluster
				j1=i1-2;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(j2==N-1){
				if( dummyC[i2][j2][k2]*dummyA[i2][j2-N+1][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*(N-1)-2]!=-1){
					//bond formed --- include in the cluster
					j1=(i1-4*N*(N-1)-2);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else
			{
				if( dummyC[i2][j2][k2]*dummyA[i2][j2+1][j2]==1 && randomnumber()<=p_bond && DUMMY[(i1+4*N-2)%(4*N*N*N)]!=-1){
					//bond formed --- include in the cluster
					j1=(i1+4*N-2)%(4*N*N*N);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			if(dummyC[i2][j2][k2]*dummyB[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-1]!=-1){
				//bond formed --- include in the cluster
				j1=i1-1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(j2==N-1 && k2!=0){
				if(dummyC[i2][j2][k2]*dummyB[i2][(j2+1-N)%N][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*(N-1)-5]!=-1){
					// bond formed --- include in the cluster
					j1=(i1-4*N*(N-1)-5);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if (j2==N-1 && k2==0){
				if(dummyC[i2][j2][k2]*dummyB[i2][0][N-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N-5]!=-1){
					j1=i1-4*N-5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			else if(j2!=N-1 && k2!=0)
			{
				if(dummyC[i2][j2][k2]*dummyB[i2][(j2+1)%N][(k2-1+N)%N]==1 && randomnumber()<=p_bond && DUMMY[(i1+4*N-5)%(4*N*N*N)]!=-1){
					// bond formed --- include in the cluster
					j1=(i1+4*N-5)%(4*N*N*N);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(j2!=N-1 && k2==0) 
			{
				if(dummyC[i2][j2][k2]*dummyB[i2][j2+1][(k2-1+N)%N]==1 && randomnumber()<=p_bond && DUMMY[i1+8*N-5]!=-1){
					// bond formed --- include in the cluster
					j1=(i1+8*N-5);
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

			if(dummyC[i2][j2][k2]*dummyD[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1+1]!=-1){
				//bond formed --- include in the cluster
				j1=i1+1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(i2!=0 && k2!=N-1){
				if(dummyC[i2][j2][k2]*dummyD[i2-1][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N+5]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N+5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2!=N-1)
			{
				if(dummyC[i2][j2][k2]*dummyD[(i2-1+N)%N][j2][k2+1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)+5]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)+5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=0 && k2==N-1)
			{
				if(dummyC[i2][j2][k2]*dummyD[i2-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N-4*(N-1)+1]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N-4*(N-1)+1;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==0 && k2==N-1)
			{
				if(dummyC[i2][j2][k2]*dummyD[i2+N-1][j2][0]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N*(N-1)-4*(N-1)+1]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N*(N-1)-4*(N-1)+1;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}

		}

		if((i1+1)%4==0)   // i1 belongs to D sublattice
		{
			i2=i1/(4*N*N);
			j2=(i1-i2*(4*N*N))/(4*N);
			k2=(i1-i2*(4*N*N)-j2*(4*N)+1)/4;
			k2=k2-1;

			if(dummyD[i2][j2][k2]*dummyA[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-3]!=-1){
				//bond formed --- include in the cluster
				j1=i1-3;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(dummyD[i2][j2][k2]*dummyB[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-2]!=-1){
				//bond formed --- include in the cluster
				j1=i1-2;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(dummyD[i2][j2][k2]*dummyC[i2][j2][k2]==1 && randomnumber()<=p_bond && DUMMY[i1-1]!=-1){
				//bond formed --- include in the cluster
				j1=i1-1;
				DUMMY[j1]=-1;
				b=b+1;
				cluster[b]=j1;                                           //Add neighbour to cluster
				perimeter[a]=j1;                                        // Add neighbour to perimeter
				a=a+1;
				Nelements=Nelements+1;
			}

			if(i2!=N-1 && k2!=0){
				if(dummyD[i2][j2][k2]*dummyA[i2+1][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N-7]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N-7;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2!=0)
			{
				if(dummyD[i2][j2][k2]*dummyA[(i2+1)%N][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)-7]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)-7;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyA[i2+1][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N+4*(N-1)-3]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N+4*(N-1)-3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyA[(i2+1)%N][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)+4*(N-1)-3]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)+4*(N-1)-3;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}


			if(i2!=N-1 && k2!=0){
				if(dummyD[i2][j2][k2]*dummyB[i2+1][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N-6]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N-6;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2!=0)
			{
				if(dummyD[i2][j2][k2]*dummyB[(i2+1)%N][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)-6]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)-6;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyB[i2+1][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N+4*(N-1)-2]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N+4*(N-1)-2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyB[(i2+1)%N][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)+4*(N-1)-2]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)+4*(N-1)-2;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}


			if(i2!=N-1 && k2!=0){
				if(dummyD[i2][j2][k2]*dummyC[i2+1][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N-5]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N-5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2!=0)
			{
				if(dummyD[i2][j2][k2]*dummyC[(i2+1)%N][j2][k2-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)-5]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)-5;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2!=N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyC[i2+1][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1+4*N*N+4*(N-1)-1]!=-1){
					//bond formed --- include in the cluster
					j1=i1+4*N*N+4*(N-1)-1;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}
			else if(i2==N-1 && k2==0)
			{
				if(dummyD[i2][j2][k2]*dummyC[(i2+1)%N][j2][k2+N-1]==1 && randomnumber()<=p_bond && DUMMY[i1-4*N*N*(N-1)+4*(N-1)-1]!=-1){
					//bond formed --- include in the cluster
					j1=i1-4*N*N*(N-1)+4*(N-1)-1;
					DUMMY[j1]=-1;
					b=b+1;
					cluster[b]=j1;                                           //Add neighbour to cluster
					perimeter[a]=j1;                                        // Add neighbour to perimeter
					a=a+1;
					Nelements=Nelements+1;
				}
			}


		}

		DUMMY[i1]=-1;



		if(k1==Nelements)
			flag1=0;


		k1=k1+1;         
		a=0;



	}      

	Nelements=Nelements+1;b=b+1;                     // Original seed has to be included 



	for(j1=0;j1<b;j1++){              // Flip all spins of cluster 
		lattice[cluster[j1]]=-lattice[cluster[j1]];
	}


	b=0;Nelements=0;
	c=0;
	flag1=1;
	return(mag());
}








int main()
{

	float Beta;                               // Inverse of Temperature
	int tau,sample_size;               // iterators for thermalization and averaging thermodynamic qtys.
	double Mag_initial;                       // Initial Magnetization
	double Avg_Mag,Avg_Magsquared;
	double Total_Mag=0.0;
	double Total_Magsquared=0.0; 
	double sigma;                             // standard deviation for average  magnetization per site



	cluster=malloc(4*N*N*N*sizeof(short));perimeter=malloc(4*N*N*N*sizeof(short));

	for (Beta=0.04;Beta<2.5;Beta+=0.002)
	{

		srand(time(NULL));

		initialize();


		/*Thermalization*/

		for (tau=0;tau<THSWEEPS;tau++)
		{
			sweep(Beta);
		}


		/* Calculation of Average Magnetization */

		for (sample_size=0;sample_size<NSAMPLES;sample_size++)
		{
			Mag_initial=fabs(sweep(Beta));
			Total_Mag=Total_Mag+Mag_initial;
			Total_Magsquared=Total_Magsquared+(Mag_initial*Mag_initial );
		}
		Avg_Mag=(Total_Mag/(NSAMPLES));
		Avg_Magsquared=(Total_Magsquared/(NSAMPLES));


		Total_Mag=0.0;
		Total_Magsquared=0.0;
		Avg_Mag=Avg_Mag/(4*N*N*N);
		Avg_Magsquared=Avg_Magsquared/((4*N*N*N)*(4*N*N*N));
		sigma = sqrt((Avg_Magsquared-Avg_Mag*Avg_Mag));
		printf("%4f   %.4f   %.4f \n",Beta,Avg_Mag,sigma/sqrt(NSAMPLES-1));


	}


	free(perimeter);
	free(cluster);


	return(0);
}







