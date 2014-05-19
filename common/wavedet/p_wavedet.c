#include 	"mex.h"
#include        <stdio.h>
/*#define NUMBER_OF_STRUCTS (sizeof(regularity)/sizeof(struct reg))
#define NUMBER_OF_FIELDS (sizeof(field_reg)/sizeof(*field_reg))*/
#define NUMBER_OF_STRUCTS_POS (sizeof(position)/sizeof(struct pos))
#define NUMBER_OF_FIELDS_POS (sizeof(field_pos)/sizeof(*field_pos))
/*struct reg
{      
  double *time,*amp2,*amp3,*amp4,*alphalinha,*alpha1,*alpha2,*alpha3,*alpha4;
  double *thalpha1,*thalpha2,*thalpha3,*th2,*th3,*th4,*thlinha,*escolhidos;
};*/

struct pos
{
 double Pon, *P, *Poff, *QRSon, *Q, *R, *Fiducial, *qrs, *Rprima, *S, *QRSoff;
 double *Ton, *T, *Tprima, *Toff, *Ttipo;
};

void mexFunction(
int	nlhs,
mxArray	*plhs[],
int	nrhs,
const mxArray	*prhs[]
)
{
    /* Input variables declaration */
    char *sigdir, *headir, *ecgnr, *matdir, *anot;
    double ft, lead, qrs_flag;
    double *time;
    
    /* Struct array "regularity" declaration */
    /*const char *field_reg[] = {"time", "amp2", "amp3", "amp4", "alphalinha",
    "alpha1", "alpha2", "alpha3", "alpha4", "thalpha1", "thalpha2", "thalpha3",
    "th2", "th3", "th4", "thlinha", "escolhidos"};
    struct reg regularity[] = {{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},
    {NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},
    {NULL}};    
    int dims[2] = {1, NUMBER_OF_STRUCTS};
    double *time_field, *amp2_field, *amp3_field, *amp4_field, *alphalinha_field;
    double *alpha1_field, *alpha2_field, *alpha3_field, *alpha4_field, *thalpha1_field, *thalpha2_field, *thalpha3_field;
    double *th2_field, *th3_field, *th4_field, *thlinha_field, *escolhidos_field;*/
    
    /* Struct array "position" output declaration */
    const char *field_pos[] = {"Pon", "P", "Poff", "QRSon", "Q", "R", "Fiducial", "qrs", "Rprima",
    "S", "QRSoff", "Ton", "T", "Tprima", "Toff", "Ttipo"};
    struct pos position[] = {{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},
    {NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL}};
    int dims_pos[2] = {1, NUMBER_OF_STRUCTS_POS};    
    double *Pon_field, *P_field, *Poff_field, *QRSon_field, *Q_field, *R_field;
    double *Fiducial_field, *qrs_field, *Rprima_field, *S_field, *QRSoff_field;
    double *Ton_field, *T_field, *Tprima_field, *Toff_field, *Ttipo_field;
    
    
    /* Checking inputs */
  if (nrhs < 9)
    {
        mexErrMsgTxt ("Wavedet: nine input arguments required.");
    } 
  
  if (!mxIsChar(prhs[0]))
      mexErrMsgTxt("Signal directory must be a string.");
  else if (mxGetM(prhs[0]) != 1)
      mexErrMsgTxt("Signal directory must be a string in a row vector.");
  
  if (!mxIsChar(prhs[1]))
      mexErrMsgTxt("Header signal directory must be a string.");
  else if (mxGetM(prhs[1]) != 1)
      mexErrMsgTxt("Header signal directory must be a string in a row vector.");
  
  if (!mxIsChar(prhs[2]))
      mexErrMsgTxt("Output directory must be a string.");
  else if (mxGetM(prhs[2]) != 1)
      mexErrMsgTxt("Output directory must be a string in a row vector.");
  
  if (!mxIsChar(prhs[3]))
      mexErrMsgTxt("Name of the signal must be a string.");
  else if (mxGetM(prhs[3]) != 1)
      mexErrMsgTxt("Name of the signal must be a string in a row vector.");
  
  if (mxIsChar(prhs[4]))
      mexErrMsgTxt("Format file must be an integer: 0 (MIT), 1 (LUND) or 2 (mat files).");
  else if ((mxGetM(prhs[4]) != 1)||(mxGetN(prhs[4]) != 1))
      mexErrMsgTxt("Format file must be an integer!!");
  else
  {
      ft = mxGetScalar(prhs[4]);
      if ((ft != 1)&&(ft != 0)&&(ft != 2))
          mexErrMsgTxt("Format file must be an integer: 0, 1 or 2.");      
      
  }
  if (!mxIsChar(prhs[5]))
      mexErrMsgTxt("Annotation output file must be a string.");
  else if (mxGetM(prhs[5]) != 1)
      mexErrMsgTxt("Annotation output file must be a string in a row vector.");
  
  if (mxIsChar(prhs[6]))
      mexErrMsgTxt("Number of lead to analyze must be an integer.");
  else if ((mxGetM(prhs[6]) != 1)||(mxGetN(prhs[6]) != 1))
      mexErrMsgTxt("Number of lead to analyze must be an integer.");
  
    lead = mxGetScalar(prhs[6]);
    
  if (mxIsChar(prhs[7]))
      mexErrMsgTxt("Time vector must be a 1x2 matrix.");
  else if ((mxGetM(prhs[7]) != 1)||(mxGetN(prhs[7]) != 2))
      mexErrMsgTxt("Time vector must be a 1x2 matrix.");
    
    time = mxGetPr(prhs[7]);
    if (time[1]<=time[0])
        mexErrMsgTxt("Bad time vector, first element greater or equal than second one");          
    
  if (mxIsChar(prhs[8]))
      mexErrMsgTxt("Qrs_flag must be an integer: 1 (External), 0 (Internal) QRS Detector.");
    else
    {
        qrs_flag = mxGetScalar(prhs[8]);
        if ((qrs_flag != 1)&&(qrs_flag != 0))
            mexErrMsgTxt("Qrs_flag must be a number: 0 or 1.");
    }
  /* End of checking inputs */
    
    sigdir = mxArrayToString(prhs[0]);
    headir = mxArrayToString(prhs[1]);
    matdir = mxArrayToString(prhs[2]);
    ecgnr = mxArrayToString(prhs[3]);
    anot = mxArrayToString(prhs[5]);
//    plhs[0] = mxCreateString(anot);
    //plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    //lead = mxGetPr(plhs[0]);
    //lead = mxGetPr(prhs[6]);
    //mexPrintf("%f & %f",time[0],time[1]);
    
    /* Create an structure output array */
    
    plhs[0] = mxCreateStructArray(2, dims_pos, NUMBER_OF_FIELDS_POS, NULL);
    wavedet(sigdir,headir,matdir,ecgnr,anot,ft,lead,qrs_flag,time);
    mexPrintf("%f & %f",time[0],time[1]);
    
}

void wavedet(char *sigdir, char *headir, char *matdir, char *ecgnr, char *anot, 
double ft, double lead, double qrs_flag, double *time)
{
      
    
    
}









