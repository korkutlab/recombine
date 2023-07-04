#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "flsa.h"

#define EPSILON 1e-10

// This code has been adapted from the SLEP package (Liu, et al, 2011)

/*

  Some mathematical background.

  In this file, we discuss how to solve the following subproblem,

        min_x  1/2 \|x-v\|^2  + lambda \|A x\|_1,                 (1)

  which is a key problem used in the Fused Lasso Signal Approximator (FLSA).

  Also, note that, FLSA is a building block for solving the optimation problmes with fused Lasso penalty.
  
  In (1), x and v are n-dimensional vectors, 
        and A is a matrix with size (n-1) x n, and is defined as follows (e.g., n=4):
		                                         A= [ -1  1  0  0;
												       0  -1 1  0;
													   0  0  -1 1]

  The above problem can be reformulated as the following equivalent min-max optimization problem

        min_x  max_z  1/2 \|x-v\|^2  + <A x, z>
		subject to   \|z\|_{infty} \leq lambda                     (2)


  It is easy to get that, at the optimal point

                         x = v - AT z,                             (3)

  where z is the optimal solution to the following optimization problem

        min_z  1/2  z^T A AT z - < z, A v>,
		subject to  \|z\|_{infty} \leq lambda                      (4)


  
  Let B=A A^T. It is easy to get that B is a (n-1) x (n-1) tridiagonal matrix.
  When n=5, B is defined as:
                                                B= [ 2  -1   0    0;
												     -1  2   -1   0;
													 0  -1   2    -1;
													 0   0   -1   2]

  Let z0 be the solution to the linear system:

                               A A^T * z0 = A * v                  (5)

  The problem (5) can be solve by the Thomas Algorithm, in about 5n multiplications and 4n additions.

  It can also be solved by the Rose's Algorithm, in about 2n multiplications and 2n additions.

  Moreover, considering the special structure of the matrix A (and B), 
       it can be solved in about n multiplications and 3n additions

  If lambda \geq \|z0\|_{infty}, x_i= mean(v), for all i, 
                the problem (1) admits near analytical solution


  We have also added the restart technique, please refer to our paper for detail!

*/


/*
///////////////    Solving the Linear System via Thomas's Algorithm \\\\\\\\\\\\\\\\\\
*/

void Thomas(double *zMax, double *z0, double * Av, int nn){

	/*

	We apply the Tomas algorithm for solving the following linear system
	                   B * z0 = Av


    Thomas algorithm is also called the tridiagonal matrix algorithm

  B=[ 2  -1   0    0;
	  -1  2   -1   0;
	  0  -1   2    -1;
	  0   0   -1   2]

    z0 is the result,  Av is unchanged after the computation


    c is a precomputed nn dimensional vector
	c=[-1/2, -2/3, -3/4, -4/5, ..., -nn/(nn+1)]

    c[i]=- (i+1) / (i+2)
	c[i-1]=- i / (i+1)

    z0 is an nn dimensional vector
    
	*/

	int i;
	double tt, z_max;

	/*
	Modify the coefficients in Av (copy to z0)
	*/
	z0[0]=Av[0]/2;
	for (i=1;i < nn; i++){
		tt=Av[i] + z0[i-1];
		z0[i]=tt - tt / (i+2);
	}

	/*z0[i]=(Av[i] + z0[i-1]) * (i+1) / (i+2);*/
		
	/*z0[i]=(Av[i] + z0[i-1])/ ( 2 - i / (i+1));*/

	
	/*
	Back substitute (obtain the result in z0)
	*/
	z_max= fabs(z0[nn-1]);

	for (i=nn-2; i>=0; i--){

		z0[i]+=  z0[i+1] -  z0[i+1]/ (i+2);

		/*z0[i]+=  z0[i+1] * (i+1) / (i+2);*/

		tt=fabs(z0[i]);

		if (tt > z_max)
			z_max=tt;

	}
	*zMax=z_max;
	
}

			
/*
///////////////    Solving the Linear System via Rose's Algorithm \\\\\\\\\\\\\\\\\\
*/

void Rose(double *zMax, double *z0,	double * Av, int nn){

	/*
	We use the Rose algorithm for solving the following linear system
	                   B * z0 = Av


  B=[ 2  -1   0    0;
	  -1  2   -1   0;
	  0  -1   2    -1;
	  0   0   -1   2]

    z0 is the result,  Av is unchanged after the computation

    z0 is an nn dimensional vector
    
	*/

	int i, m;
	double s=0, z_max;


	/*
	We follow the style in CLAPACK
	*/
	m= nn % 5;
	if (m!=0){
		for (i=0;i<m; i++)
			s+=Av[i] * (i+1);
	}
	for(i=m;i<nn;i+=5)
		s+=   Av[i]   * (i+1) 
		    + Av[i+1] * (i+2) 
		    + Av[i+2] * (i+3) 
			+ Av[i+3] * (i+4) 
			+ Av[i+4] * (i+5);
	s/=(nn+1);


	/*
    from nn-1 to 0
	*/
	z0[nn-1]=Av[nn-1]- s;
	for (i=nn-2;i >=0; i--){
		z0[i]=Av[i] + z0[i+1];
	}

	/*
    from 0 to nn-1
	*/
	z_max= fabs(z0[0]);
	for (i=0; i<nn; i++){

		z0[i]+=  z0[i-1];

		s=fabs(z0[i]);

		if (s > z_max)
			z_max=s;

	}
	*zMax=z_max;
	
}


/*
////////////////    compute x for restarting \\\\\\\\\\\\\\\\\\\\\\\\\

x=omega(z)

v: the vector to be projected
z: the approximate solution
g: the gradient at z (g should be computed before calling this function

nn: the length of z, g, and S (maximal length for S)

n:  the length of x and v

S: records the indices of the elements in the support set
*/

int supportSet(double *x, double *v, double *z, double *g, int * S, double lambda, int nn){

	int i, j, n=nn+1, numS=0;
	double temp;


	/*
	we first scan z and g to obtain the support set S
	*/

	/*numS: number of the elements in the support set S*/
	for(i=0;i<nn; i++){
		if ( ( (z[i]==lambda) && (g[i] < EPSILON) ) || ( (z[i]==-lambda) && (g[i] >EPSILON) )){
			S[numS]=i;
			numS++;
		}
	}


	if (numS==0){ /*this shows that S is empty*/
		temp=0;
		for (i=0;i<n;i++)
			temp+=v[i];

		temp=temp/n;
		for(i=0;i<n;i++)
			x[i]=temp;

		return numS;
	}


    /*
	 Next, we deal with numS >=1
     */

	/*process the first block
	   j=0
	*/
	temp=0;
	for (i=0;i<=S[0]; i++)
		temp+=v[i];
	/*temp =sum (v [0: s[0] ]*/
	temp=( temp + z[ S[0] ] ) / (S[0] +1);
	for (i=0;i<=S[0]; i++)
		x[i]=temp;


	/*process the middle blocks
	  If numS=1, it belongs the last block
	*/
	for (j=1; j < numS; j++){
		temp=0;
		for (i= S[j-1] +1; i<= S[j]; i++){
			temp+=v[i];
		}

		/*temp =sum (v [ S[j-1] +1: s[j] ]*/

		temp=(temp - z[ S[j-1] ] + z[ S[j] ])/ (S[j]- S[j-1]);

		for (i= S[j-1] +1; i<= S[j]; i++){
			x[i]=temp;
		}
	}

	/*process the last block
	j=numS-1;
	*/
	temp=0;
	for (i=S[numS-1] +1 ;i< n; i++)
		temp+=v[i];
	/*temp =sum (v [  (S[numS-1] +1): (n-1) ]*/

	temp=( temp - z[ S[numS-1] ] ) / (nn - S[numS-1]); /*S[numS-1] <= nn-1*/

	for (i=S[numS-1] +1 ;i< n; i++)
		x[i]=temp;
	
	return numS;

}


/*
////////////  Computing the duality gap \\\\\\\\\\\\\\\\\\\\\\\\\\

we compute the duality corresponding the solution z

z: the approximate solution
g: the gradient at z (we recompute the gradient)
s: an auxiliary variable
Av: A*v

nn: the lenght for z, g, s, and Av

The variables g and s shall be revised.

The variables z and Av remain unchanged.
*/

void dualityGap(double *gap, double *z, double *g, double *s, double *Av, double lambda, int nn){

	int i, m;
	double temp;

	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];
	
	for (i=0;i<nn;i++)
		if (g[i]>0)
			s[i]=lambda + z[i];
		else
			s[i]=-lambda + z[i];
	
	temp=0;					
	m=nn%5;
	
	if (m!=0){
		for(i=0;i<m;i++)
			temp+=s[i]*g[i];
	}
	
	for(i=m;i<nn;i+=5)
		temp=temp + s[i]  *g[i]
		          + s[i+1]*g[i+1]
		          + s[i+2]*g[i+2]
		          + s[i+3]*g[i+3]
		          + s[i+4]*g[i+4];
	*gap=temp;
}


/*
Similar to dualityGap,

  The difference is that, we assume that g has been computed.
*/

void dualityGap2(double *gap, double *z, double *g, double *s, double *Av, double lambda, int nn){

	int i, m;
	double temp;

	/*
	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];
    */
	
	for (i=0;i<nn;i++)
		if (g[i]>0)
			s[i]=lambda + z[i];
		else
			s[i]=-lambda + z[i];

	temp=0;					
	m=nn%5;
	
	if (m!=0){
		for(i=0;i<m;i++)
			temp+=s[i]*g[i];
	}
	
	for(i=m;i<nn;i+=5)
		temp=temp + s[i]  *g[i]
		          + s[i+1]*g[i+1]
		          + s[i+2]*g[i+2]
		          + s[i+3]*g[i+3]
		          + s[i+4]*g[i+4];
	*gap=temp;
}


/*
generateSolution:

  generate the solution x based on the information of z and g 
  (!!!!we assume that g has been computed as the gradient of z!!!!)

*/

int generateSolution(double *x, double *z, double *gap,
					  double *v, double *Av,
					  double *g, double *s, int *S,
					  double lambda, int nn){

	int i, m, numS, n=nn+1;
	double temp, funVal1, funVal2;
	    
	/*
	z is the appropriate solution,
	and g contains its gradient
	*/


   /*
		We assume that n>=3, and thus nn>=2
		
		  We have two ways for recovering x. 
		  The first way is x = v - A^T z
		  The second way is x =omega(z)
  */

	temp=0;
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=z[i]*(g[i] + Av[i]);
    }
    for (i=m;i<nn;i+=5)
        temp=temp + z[i]  *(g[i]   + Av[i])
                  + z[i+1]*(g[i+1] + Av[i+1])
                  + z[i+2]*(g[i+2] + Av[i+2])
                  + z[i+3]*(g[i+3] + Av[i+3])
                  + z[i+4]*(g[i+4] + Av[i+4]);
    funVal1=temp /2;
    
    temp=0;
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=fabs(g[i]);
    }
    for (i=m;i<nn;i+=5)
        temp=temp + fabs(g[i])
        + fabs(g[i+1])
        + fabs(g[i+2])
        + fabs(g[i+3])
        + fabs(g[i+4]);
    funVal1=funVal1+ temp*lambda;
    
    
    /*
        we compute the solution by the second way
    */

    numS= supportSet(x, v, z, g, S, lambda, nn);
    
	/*
        we compute the objective function of x computed in the second way
    */
    
    temp=0;
    m=n%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=(x[i]-v[i]) * (x[i]-v[i]);
    }
    for (i=m;i<n;i+=5)
        temp=temp + (x[i]-  v[i]) * (  x[i]-  v[i])
                  + (x[i+1]-v[i+1]) * (x[i+1]-v[i+1])
                  + (x[i+2]-v[i+2]) * (x[i+2]-v[i+2])
                  + (x[i+3]-v[i+3]) * (x[i+3]-v[i+3])
                  + (x[i+4]-v[i+4]) * (x[i+4]-v[i+4]);
    funVal2=temp/2;
    
    temp=0;
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=fabs( x[i+1]-x[i] );
    }
    for (i=m;i<nn;i+=5)
        temp=temp + fabs( x[i+1]-x[i] )
                  + fabs( x[i+2]-x[i+1] )
                  + fabs( x[i+3]-x[i+2] )
                  + fabs( x[i+4]-x[i+3] )
                  + fabs( x[i+5]-x[i+4] );
    funVal2=funVal2 + lambda * temp;
    
    if (funVal2 > funVal1){  /*
                                  we compute the solution by the first way
                              */
        x[0]=v[0] + z[0];
        for(i=1;i<n-1;i++)
            x[i]= v[i] - z[i-1] + z[i];
        x[n-1]=v[n-1] - z[n-2];
    }
    else{
        
        /*
        the solution x is computed in the second way
        the gap can be further reduced
        (note that, there might be numerical error)
         */
        
        *gap=*gap - (funVal1- funVal2);
        if (*gap <0)
            *gap=0;
	}

	return (numS);
}


	/*

     /////////////////////////////////////// Explanation for the function sfa \\\\\\\\\\\\\\\\\\\\\\\\\\\\

     Our objective is to solve the fused Lasso signal approximator (flsa) problem:

             min_x  g(x) 1/2 \|x-v\|^2  + lambda \|A x\|_1,                      (1)

     Let x* be the solution (which is unique), it satisfies
	    
		              0 in  x* - v +  A^T * lambda *SGN(Ax*)                     (2)

     To solve x*, it suffices to find
	                
					  y*  in A^T * lambda *SGN(Ax*)                              (3)
	 that satisfies
	    
		              x* - v + y* =0                                             (4)
	 which leads to
	                  x*= v - y*                                                 (5)

     Due to the uniqueness of x*, we conclude that y* is unique. 
	 
	 As y* is a subgradient of lambda \|A x*\|_1, 
	         we name our method as Subgradient Finding Algorithm (sfa).

     y* in (3) can be further written as
	                  
					  y*= A^T * z*                                               (6)
	 where

	                  z* in lambda* SGN (Ax*)                                    (7)

     From (6), we have
	                  z* = (A A^T)^{-1} A * y*                                   (8)

     Therefore, from the uqniueness of y*, we conclude that z* is also unique.
	 Next, we discuss how to solve this unique z*.

     The problem (1) can reformulated as the following equivalent problem:	 
	   
		 min_x  max_z  f(x, z)= 1/2 \|x-v\|^2  + <A x, z>
		 subject to   \|z\|_{infty} \leq lambda                                  (9)

     At the saddle point, we have
              
				        x = v - AT z,                                            (10)

	 which somehow concides with (5) and (6)

     Plugging (10) into (9), we obtain the problem
	        
			  min_z  1/2  z^T A AT z - < z, A v>,
		      subject to  \|z\|_{infty} \leq lambda,                             (11)

     In this program, we apply the Nesterov's method for solving (11).


    Duality gap:
	
	At a given point z0, we compute x0= v - A^T z0.
	It is easy to show that
	                  min_x f(x, z0) = f(x0, z0) <= max_z f(x0, z)               (12)

    Moreover, we have
	                  max_z f(x0, z) - min_x f(x, z0) 
					     <= lambda * \|A x0\|_1 - < z0, Av - A A^T z0>           (13)

    It is also to get that
	     
		              f(x0, z0) <= f(x*, z*) <= max_z f(x0, z)                   (14)

					  g(x*)=f(x*, z*)                                            (15)

                      g(x0)=max_z f(x0, z)                                       (17)

    Therefore, we have

	                  g(x0)-g(x*) <= lambda * \|A x0\|_1 - < z0, Av - A A^T z0>  (18)


    We have applied a restarting technique, which is quite involved; and thus, we do not explain here.

     /////////////////////////////////////// Explanation for the function sfa \\\\\\\\\\\\\\\\\\\\\\\\\\\\
	*/


/*
////////////               sfa              \\\\\\\\\\\\\\\\\\\\\

For sfa_one, We do one gradient descent, and then use the restart technique at the beginning of each iteration.

The stepsize of the graident descent method is fixed to 1/4, so that no line search is needed.

    Explanation of the parameters:
    
	Output parameters
	x:    the solution to the primal problem
	gap:  the duality gap (pointer)

    Input parameters
	z:    the solution to the dual problem (before calling this function, z contains a starting point)
               !!!!we assume that the starting point has been successfully initialized in z !!!!
	z0:   a variable used for multiple purposes:
			  1) the previous solution z0
			  	 z0 may be assigned from one of the three condictions:
    			 	(i) 0;
					(ii)) the restarted value from the restart technique;
					(iii) warm start from other calculations.
			  2) the difference between z and z0, i.e., z0=z- z0

	lambda:   the regularization parameter (and the radius of the infity ball, see (11)).
	nn:       the length of z, z0, Av, g, and s
	maxStep:  the maximal number of iterations

	v:    the point to be projected (not changed after the program)
	Av:   A*v (not changed after the program)

	s:        the search point (used for multiple purposes)
	g:        the gradient at g (and it is also used for multiple purposes)

	tol:      the tolerance of the gap
	tau:  the duality gap or the restarting technique is done every tau steps

     We would like to emphasis that the following assumptions 
	       have been checked in the functions that call this function:
		   1) 0< lambda < z_max
		   2) nn >=2
		   3) z has been initialized with a starting point
		   4) z0 has been initialized with all zeros
		   
The termination condition is checked every tau iterations.

  	For the duality gap, please refer to (12-18)

*/

int sfa_one(double *x,     double *gap, int * activeS,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau){

	int i, iterStep, m, tFlag=0; //, n=nn+1;
	double temp;
	int* S=(int *) malloc(sizeof(int)*nn);
	double gapp=-1, gappp=-2;	/*gapp denotes the previous gap*/
	int numS=-100, numSp=-200, numSpp=-300;    
	/*
	numS denotes the number of elements in the Support Set S
	numSp denotes the number of elements in the previous Support Set S
	*/

	*gap=-1; /*initialize *gap a value*/

	/*
	we first do a gradient step based on z
	*/

	/*
	---------------------------------------------------
		A gradient step  begins
	*/
	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

	/* 
	do a gradient step based on z to get the new z
	*/
	m=nn%5;
	if (m!=0){
		for(i=0;i<m; i++)
			z[i]=z[i] - g[i]/4;
	}
	for (i=m;i<nn; i+=5){			
		z[i]   = z[i]   -  g[i]  /4;
		z[i+1] = z[i+1] -  g[i+1]/4;
		z[i+2] = z[i+2] -  g[i+2]/4;
		z[i+3] = z[i+3] -  g[i+3]/4;
		z[i+4] = z[i+4] -  g[i+4]/4;
	}

	/*
	project z onto the L_{infty} ball with radius lambda
	z is the new approximate solution
	*/			
	for (i=0;i<nn; i++){
		if (z[i]>lambda)
			z[i]=lambda;
		else
			if (z[i]<-lambda)
				z[i]=-lambda;
	}

	/*
	---------------------------------------------------
		A gradient descent step ends
	*/


	/*compute the gradient at z*/
	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

	/*
	main loop
	*/
	for (iterStep=1; iterStep<=maxStep; iterStep++){

		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)
		*/
				
		numSpp=numSp;
		numSp=numS; /*record the previous numS*/
		numS = supportSet(x, v, z, g, S, lambda, nn);
		
		/* With x, we compute z via
		AA^T z = Av - Ax
		*/
	
		/*
		compute s= Av -Ax
		*/
		for (i=0;i<nn; i++)
			s[i]=Av[i] - x[i+1] + x[i];
					
		/*
		Apply Thomas or Rose Algorithm for solving z
		*/
		Thomas(&temp, z, s, nn);

	    /*
	    project z to [-lambda, lambda]
	    */
		for(i=0;i<nn;i++){		
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)

        we have computed a new z, based on the above relationship
		*/


		/*
		---------------------------------------------------
		  A gradient step  begins
		*/
		g[0]=z[0] + z[0] - z[1] - Av[0];
		for (i=1;i<nn-1;i++){
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];

		/* 
		do a gradient step based on z to get the new z
		*/
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda
        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++){
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/


		/*compute the gradient at z*/
		g[0]=z[0] + z[0] - z[1] - Av[0];
		for (i=1;i<nn-1;i++){
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}	
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];


		/*
		checked convergence every tau iterations
		*/
		if (iterStep % tau==0){
			gappp=gapp;
			gapp=*gap;  /*record the previous gap*/

			dualityGap2(gap, z, g, s, Av, lambda, nn);
			/*g, the gradient of z should be computed before calling this function*/
			
			if (*gap <=tol){
				tFlag=1;
				break;
			}

			m=1;
			if (nn > 1000000)
				m=5;
			else
				if (nn > 100000)
					m=3;

			if ( abs( numS-numSp) <m ){
				m=generateSolution(x, z, gap, v, Av,
					                  g, s, S, lambda, nn);
				/*g, the gradient of z should be computed before calling this function*/

				if (*gap < tol){
					numS=m;
					tFlag=2;					
					break;
				}

				if ( (*gap ==gappp) && (numS==numSpp) ){			
					tFlag=2;
					break;					
				}
				
            } /*end of if*/

		}/*end of if tau*/

	} /*end of for*/

	if (tFlag!=2){
		numS=generateSolution(x, z, gap, v, Av, g, s, S, lambda, nn);
       /*g, the gradient of z should be computed before calling this function*/
    }

	free(S);

	*activeS=numS;
	return(iterStep);
}


/*

  Now, we solve the Fused Lasso Signal Approximator (FLSA) problem:

     min_x  1/2 \|x-v\|^2  + lambda1 * \|x\|_1 + lambda2 * \|A x\|_1,      (1)

  It can be shown that, if x* is the solution to

     min_x  1/2 \|x-v\|^2  + lambda2 \|A x\|_1,                            (2)

  then 
     x**= sgn(x*) max(|x*|-lambda_1, 0)                                    (3)

  is the solution to (1).

  By some derivation (see the description in sfa.h), (2) can be solved by

     x*= v - A^T z*,

  where z* is the optimal solution to

     min_z  1/2  z^T A AT z - < z, A v>,
		subject to  \|z\|_{infty} \leq lambda2                             (4)
*/



/*

  In flsa, we solve (1) corresponding to a given (lambda1, lambda2)

  void flsa(double *x, double *z, double *gap,
		  double * v, double *z0, 
		  double lambda1, double lambda2, int n, 
		  int maxStep, double tol)

  Output parameters:
      x:        the solution to problem (1)
	  z:        the solution to problem (4)
	  infor:    the information about running the subgradient finding algorithm
	                 infor[0] = gap:         the computed gap (either the duality gap
	                                            or the summation of the absolute change of the adjacent solutions)
					 infor[1] = steps:       the number of iterations
					 infor[2] = lambad2_max: the maximal value of lambda2_max
					 infor[3] = numS:        the number of elements in the support set
								
  Input parameters:
      v:        the input vector to be projected
	  z0:       a guess of the solution of z

	  lambad1:  the regularization parameter
	  labmda2:  the regularization parameter
	  n:        the length of v and x

      maxStep:  the maximal allowed iteration steps
	  tol:      the tolerance parameter
	  tau:      the program sfa is checked every tau iterations for termination

  The input variable z0 is not modified after calling sfa.
*/


void flsa(double *x, double *z, double *infor,
		  double *v, double *z0, 
		  double lambda1, double lambda2, int n, 
		  int maxStep, double tol, int tau){

	int i, nn=n-1, m;
	double zMax, temp;
	double *Av, *g, *s;
	int iterStep, numS;
	double gap;

    Av=(double *) malloc(sizeof(double)*nn);

	/*
	Compute Av= A*v                  (n=4, nn=3)
			                                         A= [ -1  1  0  0;
												          0  -1 1  0;
													      0  0  -1 1]
	*/

	for (i=0;i<nn; i++)
		Av[i]=v[i+1]-v[i];

	/*
	Sovlve the linear system via Thomas's algorithm or Rose's algorithm
        B * z0 = Av
	*/

    Thomas(&zMax, z, Av, nn);

	/*
	Rose(&zMax, z, Av, nn);
	*/

	/*
	We consider two cases:
	   1) lambda2 >= zMax, which leads to a solution with same entry values
	   2) lambda2 < zMax, which needs to first run sfa, and then perform soft thresholding
	*/

	/*
	First case: lambda2 >= zMax
	*/
	if (lambda2 >= zMax){
		
		temp=0;
		m=n%5;
		if (m!=0){
			for (i=0;i<m;i++)
				temp+=v[i];
		}		
		for (i=m;i<n;i+=5){
			temp += v[i] + v[i+1] + v[i+2] + v[i+3] + v[i+4];
		}
		temp/=n; 
		/* temp is the mean value of v*/


		/*
		soft thresholding by lambda1
		*/
		if (temp> lambda1)
			temp= temp-lambda1;
		else
			if (temp < -lambda1)
				temp= temp+lambda1;
			else
				temp=0;

		m=n%7;
		if (m!=0){
			for (i=0;i<m;i++)
				x[i]=temp;
		}
		for (i=m;i<n;i+=7){
			x[i]   =temp;
			x[i+1] =temp;
			x[i+2] =temp;
			x[i+3] =temp;
			x[i+4] =temp;
			x[i+5] =temp;
			x[i+6] =temp;
		}
		
		gap=0;

		free(Av);

		infor[0]= gap;
		infor[1]= 0;
		infor[2]=zMax;
		infor[3]=0;

		return;
	}

	/*
	Second case: lambda2 < zMax

    We need to call sfa for computing x, and then do soft thresholding

    Before calling sfa, we need to allocate memory for g and s, 
	           and initialize z and z0.
	*/


	/*
	Allocate memory for g and s
	*/
	g = (double *) malloc(sizeof(double)*nn),
	s = (double *) malloc(sizeof(double)*nn);
	
	for (i=0;i<nn;i++){
		if (z0[i] > lambda2)
			z[i]=lambda2;
		else
			if (z0[i]<-lambda2)
				z[i]=-lambda2;
			else
				z[i]=z0[i];	
	}
	
	/*
	// call sfa_one to compute z, for finding the subgradient and x
	*/
	iterStep = sfa_one(x, &gap, &numS,
					   z,  v,   Av, 
					   lambda2, nn,  maxStep,
					   s, g,
					   tol, tau);
		
	/*
	soft thresholding by lambda1
	*/
	for(i=0;i<n;i++)
		if (x[i] > lambda1)
			x[i]-=lambda1;
		else
			if (x[i]<-lambda1)
				x[i]+=lambda1;
			else
				x[i]=0;

	
	free(Av);
	free(g);
	free(s);

	infor[0]=gap;
	infor[1]=iterStep;
	infor[2]=zMax;
	infor[3]=numS;
}
