// Asymmetric channel steady-flow with sudden expansion on one side
// using streamfunction-vorticity formulation.
// Implicit first order upwind finite difference approach 
// by constructing the successive substitution (Gauss-Seidel) formula.

// Code by Shreya Agrawal (160662)
// Written in C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 #define N 406 /* N denotes number of x nodes (columns) */
 #define M 41 /* M denotes number of y nodes (rows) */
 #define n 10000 /* n denotes number of iterations in time */
 double omega_0[M*N]; /* omega array for previous time */
 double omega_old[M*N];
 double omega_new[M*N];/* omega array for the current time */
 double psi_0[M*N];  /* stream-function for previous time*/
 double psi_old[M*N];
 double psi_new[M*N];
 double u[M*N];  /* u velocity */
 double v[M*N]; /* v velocity */

 #define PI 3.14159265359
 float H = 1; /* half-height of the channel */
 float N1 = N;
 float M1 = M;
 float Re = 200;
 float del_t = 0.01; /* time-step size */

// float rho = 1.225; /* density of fluid */
// float mu = 0.000015; /* viscosity of fluid */

 int i,j,k,p,q;
 int c1,c2;
 int ie,iw,in,is;
 float ue,uw,vn,vs; /* velocities */
 float Fo_e,Fo_w,Fo_n,Fo_s;
 float Co_e,Co_w,Co_n,Co_s;
 float aE,aW,aN,aS,aP,bP;
 float psi_aE,psi_aW,psi_aN,psi_aS,psi_aP,psi_bP;


int main(){

// FILE *fileP  /*use file handling to store all n iterations as they are getting evaluated*/
 float del_x = (15*H)/(N1-1);
 float del_y = (2*H)/(M1-1);
 float nu = 1/Re; /* kinematic viscosity */
 float U = 1; /* Inlet velocity x-component */
 long double e1 = pow(10,-3);
 long double e2 = pow(10,-4);


 for(j=0;j<M;j=j+1){
    for(i=0;i<N;i=i+1){
 /* Initial Condition at t=0- */
        omega_0[i+j*N] = 0;
        omega_old[i+j*N] = omega_0[i+j*N];
        psi_0[i+j*N] = 0;
        psi_old[i+j*N] = psi_0[i+j*N];
        u[i+j*N] = 0;
        v[i+j*N] = 0;

    }

 }
 k = 0;
 Co_e = 0; Co_w = 0; Co_s = 0; Co_n = 0; /* Initializing values */

 do{

 p = 1;

    do{

      for(j=0;j<M;j=j+1){
       for(i=0;i<N;i=i+1){

         /* B.C.s */
         /* Not sure whether the condition on bottom wall should depend on omega_old */
         if(j==0){

                psi_new[i+j*N] = 0;

         } /* At bottom wall AB*/

         else if( (j==M-1) && (i > (4*(N-1))/15) ){

                psi_new[i+j*N] = U*H;

         } /* At top wall EF*/

         else if( (i==0) && (j <= (M-1)/2) && (j!=0) ){

                psi_new[i+j*N] = U*j*del_y;

         } /* At inlet AC*/

         else if( (i==N-1) && (j!=0) ){

                psi_new[i+j*N] = psi_new[i-1+j*N];

         } /* At outlet BF*/

         else if( ( i==(4*(N-1))/15 ) && ( j > (M-1)/2 ) ){

                psi_new[i+j*N] = U*H;

         } /* At the expansion vertical wall (x=4H) DE */

         else if( (j==(M-1)/2) && (i <= (4*(N-1))/15) && (i!=0) ){

                psi_new[i+j*N] = U*H;

         } /* At y=H, wall CD */

         else if ( (j > (M-1)/2) && (i < (4*(N-1))/15) ){

               psi_new[i+j*N] = 0;

         }   /* Imaginary nodes have zero conditions */

         else{

/* Finding psi coefficients */

           psi_aE = ((-1)*del_t)/(del_x*del_x);
           psi_aW = ((-1)*del_t)/(del_x*del_x);
           psi_aN = ((-1)*del_t)/(del_y*del_y);
           psi_aS = ((-1)*del_t)/(del_y*del_y);
           psi_aP = 1 - psi_aE - psi_aN - psi_aS - psi_aW;
           psi_bP = psi_0[i+j*N] + del_t*omega_old[i+j*N];

/* Solving psi equation for interior nodes */

   psi_new[i+j*N] = ( psi_bP - ( psi_aE*psi_old[i+1+j*N] + psi_aW*psi_new[i-1+j*N] + psi_aN*psi_old[i+(j+1)*N] + psi_aS*psi_new[i+(j-1)*N] ) )/( psi_aP );


         }

       }
      }


      for(j=0;j<M;j=j+1){
       for(i=0;i<N;i=i+1){

         /* B.C.s */
         /* Not sure whether the condition on bottom wall should depend on omega_old */
         if(j==0){

                u[i+j*N] = 0;
                v[i+j*N] = 0;

         } /* At bottom wall AB*/

         else if( (j==M-1) && (i > (4*(N-1))/15) ){

                u[i+j*N] = 0;
                v[i+j*N] = 0;

         } /* At top wall EF*/

         else if( (i==0) && (j <= (M-1)/2) && (j!=0) ){

                u[i+j*N] = U;
                v[i+j*N] = 0;

         } /* At inlet AC*/

         else if( (i==N-1) && (j!=0) ){

                u[i+j*N] = ( psi_new[i+(j+1)*N] - psi_new[i+(j-1)*N] )/(2*del_y);
                v[i+j*N] = 0;

         } /* At outlet BF*/

         else if( ( i==(4*(N-1))/15 ) && ( j > (M-1)/2 ) ){

                u[i+j*N] = 0;
                v[i+j*N] = 0;

         } /* At the expansion vertical wall (x=4H) DE */

         else if( (j==(M-1)/2) && (i <= (4*(N-1))/15) && (i!=0) ){

                u[i+j*N] = 0;
                v[i+j*N] = 0;

         } /* At y=H, wall CD */

         else if ( (j > (M-1)/2) && (i < (4*(N-1))/15) ){

               u[i+j*N] = 0;
               v[i+j*N] = 0;

         }   /* Imaginary nodes have zero conditions */

         else{


  u[i+j*N] = ( psi_new[i+(j+1)*N] - psi_new[i+(j-1)*N] )/(2*del_y);
  v[i+j*N] = ( -psi_new[i+1+j*N] + psi_new[i-1+j*N] )/(2*del_x);

         }

       }
      }


      for(j=0;j<M;j=j+1){
       for(i=0;i<N;i=i+1){

         /* B.C.s */
         /* Not sure whether the condition on bottom wall should depend on omega_old */
         if(j==0){

                omega_new[i+j*N] = ( (-2)*psi_new[i+(j+1)*N] )/(del_y*del_y);

         } /* At bottom wall AB*/

         else if( (j==M-1) && (i > (4*(N-1))/15) ){

                omega_new[i+j*N] = ( 2*U*H - 2*psi_new[i+(j-1)*N] )/( del_y*del_y );

         } /* At top wall EF*/

         else if( (i==0) && (j <= (M-1)/2) && (j!=0) ){

                omega_new[i+j*N] = 0;

         } /* At inlet AC*/

         else if( (i==N-1) && (j!=0) ){

                omega_new[i+j*N] = omega_new[i-1+j*N];

         } /* At outlet BF*/

         else if( ( i==(4*(N-1))/15 ) && ( j > (M-1)/2 ) ){

                omega_new[i+j*N] = ( 2*U*H - 2*psi_new[i+1+j*N] )/( del_x*del_x );

         } /* At the expansion vertical wall (x=4H) DE */

         else if( (j==(M-1)/2) && (i <= (4*(N-1))/15) && (i!=0) ){

                omega_new[i+j*N] = ( 2*U*H - 2*psi_new[i+(j-1)*N] )/( del_y*del_y );

         } /* At y=H, wall CD */

         else if ( (j > (M-1)/2) && (i < (4*(N-1))/15) ){

               omega_new[i+j*N] = 0;

         }   /* Imaginary nodes have zero conditions */

         else{


  ue = ( u[i+j*N] + u[i+1+j*N] )/2;
  uw = ( u[i+j*N] + u[i-1+j*N] )/2;
  vn = ( v[i+j*N] + v[i+(j+1)*N] )/2;
  vs = ( v[i+j*N] + v[i+(j-1)*N] )/2;;

  /* specify ie,iw,in,is */
  if(u[i+j*N]>0){
        ie=0; iw=1;

        if(uw>=0){
        Co_w = (uw*del_t)/del_x;
        }else{
        Co_w = ((-1)*uw*del_t)/del_x;
        }
  }

  else if(u[i+j*N]<0){
        ie=1; iw=0;

        if(ue>=0){
        Co_e = (ue*del_t)/del_x;
        }else{
        Co_e = ((-1)*ue*del_t)/del_x;
        }
  }

  else if(u[i+j*N]==0){
        ie=0; iw=0;
  }

  if(v[i+j*N]>0){
        in=0; is=1;

        if(vs>=0){
        Co_s = (vs*del_t)/del_y;
        }else{
        Co_s = ((-1)*vs*del_t)/del_y;
        }
  }

  else if(v[i+j*N]<0){
        in=1; is=0;

        if(vn>=0){
        Co_n = (vn*del_t)/del_y;
        }else{
        Co_n = ((-1)*vn*del_t)/del_y;
        }
  }

  else if(v[i+j*N]==0){
        in=0; is=0;
  }

/* Solving omega equation */
            Fo_e = ( nu*del_t )/(del_x*del_x);
            Fo_w = ( nu*del_t )/(del_x*del_x);
            Fo_n = ( nu*del_t )/(del_y*del_y);
            Fo_s = ( nu*del_t )/(del_y*del_y);

            aE = -Fo_e-(ie*Co_e);
            aW = -Fo_w-(iw*Co_w);
            aN = -Fo_n-(in*Co_n);
            aS = -Fo_s-(is*Co_s);
            aP = 1-aE-aW-aS-aN;
            bP = omega_0[i+j*N];

            omega_new[i+j*N] = ( bP-(aE*omega_old[i+1+j*N])-(aW*omega_new[i-1+j*N])-(aN*omega_old[i+(j+1)*N])-(aS*omega_new[i+(j-1)*N]) )/aP;

//            printf("omega = %lf \n",omega_new[i+j*N]);
         }

       }
      }


        c1 = 0;

       for(j=0;j<M;j=j+1){          /* stopping criterion: convergence */
       for(i=0;i<N;i=i+1){
       if( (((-1)*e1)<(omega_new[i+j*N]-omega_old[i+j*N])) && (e1>(omega_new[i+j*N]-omega_old[i+j*N])) && (((-1)*e1)<(psi_new[i+j*N]-psi_old[i+j*N])) && (e1>(psi_new[i+j*N]-psi_old[i+j*N])) ){
        c1 = c1 + 1;
       }
       }
       }

    if(c1==M*N){
        break;
    }else{
        for(j=0;j<M;j=j+1){
        for(i=0;i<N;i=i+1){
        psi_old[i+j*N] = psi_new[i+j*N];
        omega_old[i+j*N] = omega_new[i+j*N];
        }
    }
    }

    p = p+1;
    }
    while(p<5000);

  if(p==5000){ printf("Convergence not reached for time iteration number (n) = %d !\n",k);}


 c2 = 0;

       for(j=0;j<M;j=j+1){          /* stopping criterion: steady-state achieved */
       for(i=0;i<N;i=i+1){
       if( (((-1)*e2)<(psi_new[i+j*N]-psi_0[i+j*N])) && (e2>(psi_new[i+j*N]-psi_0[i+j*N])) && (((-1)*e2)<(omega_new[i+j*N]-omega_0[i+j*N])) && (e2>(omega_new[i+j*N]-omega_0[i+j*N])) ){
        c2 = c2 + 1;

       }
       }
       }
   /* store temp to permanent */
    for(j=0;j<M;j=j+1){
    for(i=0;i<N;i=i+1){

 omega_0[i+j*N] = omega_new[i+j*N];
 omega_old[i+j*N] = omega_new[i+j*N];
 psi_0[i+j*N] = psi_new[i+j*N];
 psi_old[i+j*N] = psi_new[i+j*N];

 }
    }



  k = k+1;

   if(c2==M*N){
        break;
    }

 }
 while(k<n);
 if(k==n){ printf("Steady state not reached\n");}

printf("u velocity: \n");

for(j=M-1;j>=0;j=j-1){
for(i=0;i<N;i=i+1){
printf("%lf \t", u[i+j*N] );
}
printf("\n");
}

printf("\n");
printf("psi: \n");

for(j=M-1;j>=0;j=j-1){
for(i=0;i<N;i=i+1){
printf("%lf \t", psi_new[i+j*N] );
}
printf("\n");
}

printf("\n");
printf("omega: \n");

for(j=M-1;j>=0;j=j-1){
for(i=0;i<N;i=i+1){
printf("%lf \t", omega_new[i+j*N] );
}
printf("\n");
}
 return 0;
}
