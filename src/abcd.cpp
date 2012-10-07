#include <math.h>
#include <vector>
#include <iostream>
#include <R.h>
using namespace std;

extern "C" {
	void abcd(double *x,double *y,double *d,double *c,double *lambda1,double *lambda2,
	int *nrowa, int *ncola){
	    int i, j;
		int n=(*nrowa);
		int p=(*ncola);
	    // Convert vector parameter a from R to matrices in C:a to A
	    double **X=new double*[*nrowa];
		for(i=0;i<*nrowa;i++)
			X[i]=new double[*ncola];
	    for (i=0; i<*nrowa; i++){
	        for (j=0; j<*ncola; j++){
	            X[i][j] = x[j * (*nrowa) + i];
	        }
	    }

	    double **Y=new double*[*nrowa];
		for(i=0;i<*nrowa;i++)
			Y[i]=new double[*ncola];

	    for (i=0; i<*nrowa; i++){
	        for (j=0; j<*ncola; j++){
	            Y[i][j] = y[j * (*nrowa) + i];
	        }
	    }

		double **D=new double*[n];
		for(i=0;i<n;i++)
			D[i]=new double[p];
		for (i=0; i<*nrowa; i++){
			        for (j=0; j<*ncola; j++){
			            D[i][j] = d[j * (*nrowa) + i];
			        }
	    }


        double t1=0;
		for(i=0;i<n;i++){
				(t1)+=(X[i][0]+D[i][0])*(X[i][0]+D[i][0]);
		}
		t1=sqrt(t1);

		double t2=0;
		for(i=0;i<n;i++){
				(t2)+=(X[i][0]+D[i][0]-X[i][1]-D[i][1])*(X[i][0]+D[i][0]-X[i][1]-D[i][1]);
		}
		t2=sqrt(t2);

		double t6=0;
		for(i=0;i<n;i++)
			(t6)+=(X[i][p-2]+D[i][p-2]-X[i][p-1]-D[i][p-1])*(X[i][p-2]+D[i][p-2]-X[i][p-1]-D[i][p-1]);
			t6=sqrt(t6);

	    double t7=0;
		for(i=0;i<n;i++)
			(t7)+=(X[i][p-1]+D[i][p-1])*(X[i][p-1]+D[i][p-1]);
	    t7=sqrt(t7);
	     //chase matrix A:
		    double **A=new double*[p];
		    for(i=0;i<p;i++)
				A[i]=new double[3];
		    A[0][0]=2+(*lambda1)/t1+(*lambda2)/t2;
		    A[0][1]=-2*(*lambda2)/t2;
		    A[0][2]=0;
		    A[p-1][0]=-(*lambda2)/t6;
		    A[p-1][1]=2+(*lambda1)/t7+(*lambda2)/t6;
		    A[p-1][2]=0;

			double *t3=new double[p-2];
			for(j=0;j<p-2;j++){
				t3[j]=0;
				for(i=0;i<n;i++)
					t3[j]+=(X[i][j]+D[i][j]-X[i][j+1]-D[i][j+1])*(X[i][j]+D[i][j]-X[i][j+1]-D[i][j+1]);
				t3[j]=-(*lambda2)/sqrt(t3[j]);
			}

			double *t5=new double[p-2];
		    for(j=0;j<p-2;j++){
				t5[j]=0;
				for(i=0;i<n;i++)
					t5[j]+=(X[i][j+1]+D[i][j+1]-X[i][j+2]-D[i][j+2])*(X[i][j+1]+D[i][j+1]-X[i][j+2]-D[i][j+2]);
				t5[j]=-(*lambda2)/sqrt(t5[j]);
			}


			double *t4=new double[p-2];
		    for(j=0;j<p-2;j++){
				t4[j]=0;
				for(i=0;i<n;i++)
					t4[j]+=(X[i][j+1]+D[i][j+1])*(X[i][j+1]+D[i][j+1]);
				t4[j]=(*lambda1)/sqrt(t4[j]);
			}

			for(i=1;i<p-1;i++){
				A[i][0]=t3[i-1];
				A[i][1]=2+t4[i-1]-t3[i-1]-t5[i-1];
				A[i][2]=t5[i-1];
			}
	//chase matrix A:



	//chase matrix B:
	    double **f=new double*[*nrowa];
		for(i=0;i<*nrowa;i++)
			f[i]=new double[*ncola];
	    for (i=0; i<*nrowa; i++){
	        for (j=0; j<*ncola; j++){
	            f[i][j] = 2*(X[i][j]-Y[i][j]);
	        }
	    }

		double **B=new double*[p];
		for(i=0;i<p;i++)
			B[i]=new double[n];
		for(j=0;j<n;j++){
			B[0][j]=-f[j][0]-((*lambda1)/t1+(*lambda2)/t2)*X[j][0]+(*lambda2)/t2*X[j][1];
			B[p-1][j]=-f[j][p-1]-(*lambda1)/t7*X[j][p-1]+(*lambda2)/t6*(X[j][p-2]-X[j][p-1]);
		}

		double *aj=new double[p-2];
		for(i=0;i<p-2;i++)
			aj[i]=t3[i]-t4[i]+t5[i];

		for(i=1;i<p-1;i++){
			for(j=0;j<n;j++){
				B[i][j]=-f[j][i]+aj[i-1]*X[j][i]-t3[i-1]*X[j][i-1]-t5[i-1]*X[j][i+1];
			}
	    }

//	    for(i=0;i<p;i++){
//			for(j=0;j<3;j++){
//				a[j*(p)+i]=A[i][j];
//			}
//		}

//		for(i=0;i<p;i++){
//			for(j=0;j<n;j++){
//				b[j*p+i]=B[i][j];
//			}
//		}


		double **C = new double*[p];
		for (i=0; i<p; i++){
		     C[i] = new double[n];
        }

        double *beta=new double[p];
		beta[0]=(A[0][1])/(A[0][0]);
		for(i=1;i<p;i++){
		 	beta[i]=A[i][2]/(A[i][1]-A[i][0]*beta[i-1]);
        }

        double *y1=new double[p];
        for(j=0;j<n;j++){
			y1[0]=B[0][j]/A[0][0];
			for(i=1;i<p;i++)
				y1[i]=(B[i][j]-A[i][0]*y1[i-1])/(A[i][1]-A[i][0]*beta[i-1]);
//
			C[p-1][j]=y1[p-1];
			for(i=p-2;i>=0;i--)
		       C[i][j]=y1[i]-beta[i]*C[i+1][j];
	 }

	 // Convert C matrix to parameter c from R
	     for (i=0; i<p; i++){
	         for (j=0; j<n; j++){
	             c[j * (p) + i] = C[i][j];
	         }
     }


		delete []t3;
		delete []t4;
		delete []t5;
		delete []aj;
		delete []beta;
		delete []y1;

        for(i=0;i<n;i++){
			delete [] X[i];
			delete [] Y[i];
			delete [] D[i];
			delete [] f[i];
		}
		delete [] X;
		delete [] Y;
		delete [] D;
		delete [] f;
		for(i=0;i<p;i++){
			delete [] A[i];
			delete [] B[i];
			delete [] C[i];
		}
		delete [] A;
		delete [] B;
		delete [] C;
	}
}

