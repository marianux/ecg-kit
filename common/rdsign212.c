#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
/* #include <malloc.h> */
#include "mex.h"

/*  Gateway Routine */

void mexFunction(
		 int	nlhs,
		 mxArray	*plhs[],
		 int	nrhs,
		 const mxArray	*prhs[]
		 )
{
   double  *ecg;
   short  nleads,flag=0,flag2=0,j;
   long   n_inicio, n_final, Nc,i,low,high;
   FILE   *fp;
   char   *nombre;
   unsigned char  buf[3];
   unsigned int aux_ui;

   
   /* Compruebo el numero OK de los argumentos */

   if (nrhs != 4){
      mexErrMsgTxt("Write four input arguments");
      mexErrMsgTxt("Filename, number of leads, initial sample, end sample");
      return;
   } else if (nlhs > 1){
      mexErrMsgTxt("Write one output argument");
      return;
   }

   /* Debo crear un matriz para el argumento de salida */
  
   j        = (short) mxGetN(prhs[0])+1;
   nombre   = mxCalloc((int)j, sizeof(char));  
   mxGetString(prhs[0], nombre, j);
   nleads = (short) mxGetScalar(prhs[1]);
   n_inicio = (long) mxGetScalar(prhs[2]);
   n_final  = (long) mxGetScalar(prhs[3]);



   Nc = n_final - n_inicio + 1;    /* number of samples to read */

   plhs[0]  = mxCreateDoubleMatrix(Nc, nleads, mxREAL);

   if ((n_inicio*nleads) & 0x1) {flag2=1;;}

   if (plhs[0] == NULL){
       printf("Insufficient memory in rdsign212\n");
       exit(1);
   }

   /*
    * Referenciar argumentos y
    * llamar a la Computational Routine
    */

    ecg       = mxGetPr(plhs[0]);

/* Computational Routine */ 
  
   if((fp=fopen(nombre, "rb")) == NULL){
      	mexErrMsgTxt("Error: Unable to open input file.\n");     }
   else { 


        if (flag2){
           fseek(fp, (int) (((nleads*n_inicio)-1)*1.5) , 0);
	   aux_ui = fread(buf, sizeof(char), 3, fp);
   	   low=buf[1]&0X0F;
           high=buf[1]&0XF0;
	   flag=1;}
 	else
 	   fseek(fp, (int) nleads*1.5*n_inicio, 0);
 	
	
	for (i=0; i<Nc; i++){
	  for (j=0;j<nleads;j++){
 		switch (flag){
		  case 0:
		 	 aux_ui = fread(buf, sizeof(char), 3, fp);  
           		 low=buf[1]&0X0F;
          		 high=buf[1]&0XF0;
          		 if (low>7)
  	 		     ecg[i+Nc*j] = (buf[0]+(low<<8)-4096);
          		 else
  	  		     ecg[i+Nc*j] = (buf[0]+(low<<8));
			  flag++;
		          break;
		  case 1:
          		 if (high>127)
    	 		     ecg[i+Nc*j] = (buf[2]+(high<<4)-4096);
          		 else
    	   		     ecg[i+Nc*j] = (buf[2]+(high<<4));
			 flag=0;
			 break;
         	  }
	   }

        }

	fclose(fp);
   }
   mxFree(nombre);
}
