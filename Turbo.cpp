#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <random>

#define debug 0

int frame = 0;
const int block_cols = 8;
const int block_rows = 8;
const int length     = 64;
const int error_target = 200;
int times;
long double stddev;
long double Eb_N0_dB;
int max_iteration = 10;
double code_rate = 1./3.;
int Srand[length];
int Srand_inv[length];
int S = 8;

using namespace::std;

long double rand_normal(long double mean, long double stddev)
{//Box muller method
    static long double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        long double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            long double d = sqrt(-2.0*log(r)/r);
            long double n1 = x*d;
            n2 = y*d;
            long double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

long double max_star(long double a)
{
    return a;
}
long double max_star(long double a, long double b)
{
	if (a>b)
		swap(a,b);

	if (b-a > 500)	return b;
	
	return a+log(1+exp(b-a));
}
long double max_star(long double a, long double b,long double c)
{
	if (a>b)
		swap(a,b);
	if (a>c)
		swap(a,c);
	if (b>c)
		swap(b,c);
	
	if (c-b > 500)	return c;
	
	if (c-a > 500)	return max_star(b,c);
	
    return a+log(1+exp(b-a)+exp(c-a));
}
long double max_star(long double a, long double b,long double c,long double d)
{
    if (a>b)
		swap(a,b);
	if (a>c)
		swap(a,c);
	if (a>d)
		swap(a,d);
	if (b>c)
		swap(b,c);
	if (b>d)
		swap(b,d);
	if (c>d)
		swap(c,d);
	
	if (d-c > 500)	return d;
	
	if (d-b > 500)	return max_star(c,d); 
	
	if (d-a > 500)	return max_star(b,c,d);
	
    return a+log(1+exp(b-a)+exp(c-a)+exp(d-a));
}
int sign(long double a)
{
    return a>0?1:0;
}
long double max(long double a)
{
    return a;
}
long double max(long double a,long double b)
{
	return max_star(a,b);
}
long double max(long double a, long double b,long double c,long double d)
{
    //return max(max(a,b),max(c,d));
    return max_star(a,b,c,d);
}
bool interleave(int* source, int* to)				//slow, not use
{
	int** rotate;
	rotate = new int* [length];
	for (int i=0;i<length;i++)
	{
		rotate[i] = new int [length];
		for (int j=0;j<length;j++)
		{
			rotate[i][j] = 0;
		}
		rotate[i][(i%block_rows)*block_cols+(i/block_rows)] = 1;
	}
	
	for (int i=0;i<length;i++)
	{
		for (int j=0;j<length;j++)
		{
			if (rotate[i][j])
			{
				to[i] = source[j];
			}
		}
	}
	
	for (int i=0;i<length;i++)
		delete [] rotate[i];
	delete [] rotate;
	return true;
}
bool interleave_iterate(int* source, int* to)		//fast, not use
{
	for (int i=0;i<length;i++)
	{
		to[i] = source[(i%block_rows)*block_cols+(i/block_rows)];
	}
	return true;
}
int interleave_P(int pos)
{
	//return (pos%block_rows)*block_cols+(pos/block_rows);
	return Srand[pos];
}
int interleave_Pinv(int pos)
{
	//return (pos%block_rows)*block_cols+(pos/block_rows);
	return Srand_inv[pos];
}
void generates()
{
	while (1)
	{
		vector <int> vect;
		for (int i=0;i<length;i++)
		{
			vect.push_back(i);
		}
		int times;
		for (int i=0;i<length;i++)
		{
			for (times=0;times<length*10;times++)
			{
				int target = rand()%vect.size();
				bool successone = true;
				for (int j=i-S;j<i;j++)
				{
					if (j<0)	continue;
					if (abs(vect[target]-Srand[j]) <= S)
					{
						successone = false;
						break;
					}
				}
				if (successone)
				{
					Srand[i] = vect[target];
					vect.erase(vect.begin()+target);
					break;
				}	
			}
			if (times == length*10)	break;
		}
		if (times == length*10)	continue;
		//check;
		bool success = true;
		for (int i=0;i<length;i++)
		{
			for (int j=i-S;j<=i+S;j++)
			{
				if (j<0|| i==j || j>=length)	continue;
				if (abs(Srand[i]-Srand[j])<=S)
				{
					printf ("%d %d %d\n",i,j,abs(Srand[i]-Srand[j]));
					success = false;
					break;
				}
			}
			if (!success)	break;
		}
		if (success)	break;
		else
		{
			printf ("fail\n");
		}
	}
	system("cls");
	for (int i=0;i<length;i++)
	{
		printf ("%4d ",Srand[i]);
		Srand_inv[Srand[i]] = i;
	}
	printf ("\n");
	for (int i=0;i<length;i++)
	{
		printf ("%4d ",Srand_inv[i]);
	}
	printf ("\n");
	
	int mindist = INT_MAX;
	for (int i=0;i<length;i++)
	{
		for (int j=i-S;j<=i+S;j++)
		{
			if (j<0|| i==j || j>=length)	continue;
			if (abs(Srand[i]-Srand[j])<mindist)
				mindist = abs(Srand[i]-Srand[j]);
		}
		if (i != Srand_inv[Srand[i]])	printf ("Alert!\n");
	}
	printf ("min dist = %d\n",mindist);
	
	return;
}

int main()
{
    srand(time(NULL));
    generates();
    int* message_array;
    int* interleaved_message_array;
    long double* coded_array;
    long double* noise_array;
    long double* receive_array;
    long double** alpha_D1;
    long double** beta_D1;
    long double** gamma_D1;
//    long double** alpha_D2;
//    long double** beta_D2;
//    long double** gamma_D2;
    long double* L21;
    long double* L12;
    long double* L;
    
    printf ("length = %d\n",length);

	for (max_iteration = 1; max_iteration<15;max_iteration++)
	{
		char filename[1024];
		char filename_bit[1024];
		sprintf (filename,"Turbo_%d_%d.txt",length,max_iteration);
		FILE *f = fopen(filename,"w");
		fclose(f);
		sprintf (filename_bit,"Turbo_%d_%d_bit.txt",length,max_iteration);
		f = fopen(filename,"w");
		fclose(f);
		
		for (Eb_N0_dB=-3.;Eb_N0_dB<7.41;Eb_N0_dB+=0.1)
		{
			stddev = sqrt(pow(10,-Eb_N0_dB/10)/2/code_rate);
			times = 1;																				// change here
			cout <<"SNR = "<<Eb_N0_dB<<", stddev = "<<stddev<<", max_iteration = "<<max_iteration<<endl;
		    while (times)
		    {
		    	int error_count = 0;
		    	int error_count_bit = 0;
		    	frame = 0;
		    	while(1)
		    	{
		    		frame++;
			        message_array = new int[length];
			        interleaved_message_array = new int[length];
			        coded_array = new long double[length*3];
			        noise_array = new long double[length*3];
			        receive_array = new long double[length*3];
			
			        int reg0=0;
			        int reg1=0;
			        int reg2=0;
			        int reg3=0;
			        int reg4=0;
			        int reg5=0;
			
			        for (int i=0;i<length;i++)      //construct initial message
			        {
			            if (i>length-3) message_array[i]=0;
			            else message_array[i] = rand()/0.5 > RAND_MAX ? 1 : 0;
			            if (debug) printf ("%d ",message_array[i]);
			        }
			        if (debug) printf ("\n");
			        
			        interleave(message_array,interleaved_message_array);
			        
			        for (int i=0;i<length;i++)
			        {
			        	reg2 = reg1;
			            reg1 = reg0;
			            reg0 = message_array[i];
			            coded_array[3*i  ] = (reg0          )%2 ? 1 : -1;
			            coded_array[3*i+1] = (reg0+reg1+reg2)%2 ? 1 : -1;
			            reg5 = reg4;
			            reg4 = reg3;
			            reg3 = message_array[interleave_P(i)];
			            coded_array[3*i+2] = (reg3+reg4+reg5)%2 ? 1 : -1;
					}
			        
			        for (int i=0;i<3*length;i++)
			        {
			            noise_array[i] = rand_normal(0,stddev);
			            receive_array[i] = coded_array[i] + noise_array[i];
			            if (debug) cout << receive_array[i] << " ";
			            if (debug && i%3==2) printf ("\n");
			        }
			        if (debug) printf ("\n");
			        delete [] coded_array;
			        delete [] noise_array;
			
					//Parallel-Concatenated Convolutional Codes Turbo Codes
					//init:
					int route[4][4] = {{0,-1,1,-1},{2,-1,3,-1},{-1,4,-1,5},{-1,6,-1,7}};
					int iteration = 0;
					long double u,p,q;
			        alpha_D1 = new long double* [4];
			        beta_D1 = new long double* [4];
			        gamma_D1 = new long double* [8];		        
			        for (int i=0;i<4;i++)
			        {
			            alpha_D1[i] = new long double [length];
			            beta_D1[i] = new long double [length];
			        }
			        for (int i=0;i<8;i++)
			            gamma_D1[i] = new long double [length];
//			        alpha_D2 = new long double* [4];
//			        beta_D2 = new long double* [4];
//			        gamma_D2 = new long double* [8];		        
//			        for (int i=0;i<4;i++)
//			        {
//			            alpha_D2[i] = new long double [length];
//			            beta_D2[i] = new long double [length];
//			        }
//			        for (int i=0;i<8;i++)
//			            gamma_D2[i] = new long double [length];
			        
			        L21 = new long double [length];
			        L12 = new long double [length];
			        L   = new long double [length];
			        for (int i=0;i<length;i++)
			        {
			        	L21[i] = 0;
					}
					
					while (1)
					{
						// D1 (i)
				        for (int i=0;i<length;i++)
				        {
				            for (int j=0;j<8;j++)
				            {
				            	u = j==0||j==2||j==4||j==6?-1.:1.;	
				            	p = j==0||j==3||j==5||j==6?-1.:1.;
				                gamma_D1[j][i] = u*L21[interleave_Pinv(i)]/2+u*receive_array[3*i]/(stddev*stddev)+p*receive_array[3*i+1]/(stddev*stddev);				//check
				                if (debug) cout << gamma_D1[j][i] << " ";
				            }
				            if (debug) printf ("\n");
				            for (int j=0;j<4;j++)
				            {
				                if (i==0)   alpha_D1[j][i] = (j==0||j==2?max(0+gamma_D1[route[0][j]][i]):-100.);
				                else
				                {
				                    int s1 = (j==0||j==2)?0:2;
				                    int s2 = (j==0||j==2)?1:3;
				                    int r1 = route[s1][j];
				                    int r2 = route[s2][j];
				                    alpha_D1[j][i] = max(alpha_D1[s1][i-1]+gamma_D1[r1][i],alpha_D1[s2][i-1]+gamma_D1[r2][i]);
				                }
				                if (debug) cout << alpha_D1[j][i] << " ";
				            }
				            if (debug) printf ("\n");
				        }
				        //D1 (ii)
				        for (int i=length-1;i>0;i--)
				        {
				        	if (i==length-1)
				        	{
				        		for (int j=0;j<4;j++)
				        		{
				        			beta_D1[j][i] = (j==0||j==1?0:-100.);
								}
							}
				            for (int j=0;j<4;j++)
				            {
								int s1 = (j==0||j==1)?0:1;
				                int s2 = (j==0||j==1)?2:3;
				                int r1 = route[j][s1];
				                int r2 = route[j][s2];
				                beta_D1[j][i-1] = max(beta_D1[s1][i]+gamma_D1[r1][i],beta_D1[s2][i]+gamma_D1[r2][i]);
				            }
				        }
				        
				        //D1 (iii)
				        for (int i=0;i<length;i++)
				        {
				            if (i==0) 				L12[i] = max(0               +(+1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[2][i])  
													        -max(0               +(-1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[0][i]);
				        //	else if (i==length-1)   L12[i] = max(alpha_D1[0][i-1]+(-100)                                   +0)    
						//								    -max(alpha_D1[0][i-1]+(-1)*receive_array[3*i+1]/(stddev*stddev)+0            , alpha_D1[1][i-1]+(+1)*receive_array[3*i+1]/(stddev*stddev)+0);
							else    				L12[i] = max(alpha_D1[0][i-1]+(+1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[2][i], alpha_D1[1][i-1]+(-1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[2][i], alpha_D1[2][i-1]+(-1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[3][i], alpha_D1[3][i-1]+(+1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[3][i])
				                                            -max(alpha_D1[0][i-1]+(-1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[0][i], alpha_D1[1][i-1]+(+1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[0][i], alpha_D1[2][i-1]+(+1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[1][i], alpha_D1[3][i-1]+(-1)*receive_array[3*i+1]/(stddev*stddev)+beta_D1[1][i]);
						}
						
						iteration++;
						if (iteration>=max_iteration)	break;
						
						// D2 (i)
				        for (int i=0;i<length;i++)
				        {
				            for (int j=0;j<8;j++)
				            {
				            	u = j==0||j==2||j==4||j==6?-1.:1.;	
				            	q = j==0||j==3||j==5||j==6?-1.:1.;
				                gamma_D1[j][i] = u*L12[interleave_P(i)]/2+u*receive_array[3*interleave_P(i)]/(stddev*stddev)+q*receive_array[3*i+2]/(stddev*stddev);				//check
				                if (debug) cout << gamma_D1[j][i] << " ";
				            }
				            if (debug) printf ("\n");
				            for (int j=0;j<4;j++)
				            {
				                if (i==0)   alpha_D1[j][i] = (j==0||j==2?max(0+gamma_D1[route[0][j]][i]):-100.);
				                else
				                {
				                    int s1 = (j==0||j==2)?0:2;
				                    int s2 = (j==0||j==2)?1:3;
				                    int r1 = route[s1][j];
				                    int r2 = route[s2][j];
				                    alpha_D1[j][i] = max(alpha_D1[s1][i-1]+gamma_D1[r1][i],alpha_D1[s2][i-1]+gamma_D1[r2][i]);
				                }
				                if (debug) cout << alpha_D1[j][i] << " ";
				            }
				            if (debug) printf ("\n");
				        }
				        //D2 (ii)
				        for (int i=length-1;i>0;i--)
				        {
				        	if (i==length-1)
				        	{
				        		for (int j=0;j<4;j++)
				        		{
				        			beta_D1[j][i] = alpha_D1[j][i];
								}
							}
				            for (int j=0;j<4;j++)
				            {
								int s1 = (j==0||j==1)?0:1;
				                int s2 = (j==0||j==1)?2:3;
				                int r1 = route[j][s1];
				                int r2 = route[j][s2];
				                beta_D1[j][i-1] = max(beta_D1[s1][i]+gamma_D1[r1][i],beta_D1[s2][i]+gamma_D1[r2][i]);
				            }
				        }
				        
				        //D2 (iii)
				        for (int i=0;i<length;i++)
				        {
				            if (i==0) 				L21[i] = max(0               +(+1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[2][i])  
													        -max(0               +(-1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[0][i]);
				            //else if (i==length-1)   L21[i] = max(alpha_D1[0][i-1]+(-100)                                   +0)    
							//							    -max(alpha_D1[0][i-1]+(-1)*receive_array[3*i+2]/(stddev*stddev)+0            , alpha_D1[1][i-1]+(+1)*receive_array[3*i+2]/(stddev*stddev)+0);
							else    				L21[i] = max(alpha_D1[0][i-1]+(+1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[2][i], alpha_D1[1][i-1]+(-1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[2][i], alpha_D1[2][i-1]+(-1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[3][i], alpha_D1[3][i-1]+(+1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[3][i])
				                                            -max(alpha_D1[0][i-1]+(-1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[0][i], alpha_D1[1][i-1]+(+1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[0][i], alpha_D1[2][i-1]+(+1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[1][i], alpha_D1[3][i-1]+(-1)*receive_array[3*i+2]/(stddev*stddev)+beta_D1[1][i]);
						}
						
						iteration++;
						if (iteration>=max_iteration)	break;
					}
			        
					
			        //Decision
			        int* decoded_array;
			        decoded_array = new int[length];
			        for (int i=0;i<length;i++)
			        {
			        	L[i] = 2*receive_array[3*i]/(stddev*stddev)+L21[interleave_Pinv(i)]+L12[i];
						decoded_array[i] = sign (L[i]);
					}
			        
			
					//Summary
			        for (int i=0;i<length;i++)
			            if (decoded_array[i]!=message_array[i])
						{
							error_count++;
							break;					//Frame error rate
						} 
						
					for (int i=0;i<length;i++)
			            if (decoded_array[i]!=message_array[i])
						{
							error_count_bit++;
						} 
			
			        delete [] message_array;
			        delete [] interleaved_message_array;
			        delete [] receive_array;
			        for (int i=0;i<4;i++)
			        {
			            delete [] alpha_D1[i];
			            delete [] beta_D1[i];
			        }
			        for (int i=0;i<8;i++)
			            delete [] gamma_D1[i];
			        delete [] alpha_D1;
			        delete [] beta_D1;
			        delete [] gamma_D1;
//			        for (int i=0;i<4;i++)
//			        {
//			            delete [] alpha_D2[i];
//			            delete [] beta_D2[i];
//			        }
//			        for (int i=0;i<8;i++)
//			            delete [] gamma_D2[i];
//			        delete [] alpha_D2;
//			        delete [] beta_D2;
//			        delete [] gamma_D2;
			        delete [] decoded_array;
			        delete [] L21;
			        delete [] L12;
			        delete [] L;
			        //scanf ("%c",&temp);
			        if (error_count>=error_target)	break;   
				}
				printf ("frame count = %d\n",frame);
					        
			    FILE *f = fopen(filename,"a");
			    fprintf (f,"%.1lf %d %d\n",(double) Eb_N0_dB,error_count,frame);					//frame error rate
			    fclose(f);
			    
			    f = fopen(filename_bit,"a");
			    fprintf (f,"%.1lf %d %d\n",(double) Eb_N0_dB,error_count_bit,frame*length);					//frame error rate
			    fclose(f);
			    
			    char temp;
				if (debug) scanf ("%c",&temp);
				
				times--;
		    }
		}	
	}

    return 0;
}

