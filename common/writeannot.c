/* MEX-file for writting MIT format annotation files for ECG Signals */
/* from a MATLAB structure similar to the one returned by readdannot*/
/* (c) Juan Pablo Martinez Cortes.  University of Zaragoza.
   e-mail: juanpabl@tsc1.cps.unizar.es */
/*
/	INCLUDE FILES
 */

#include 	"mex.h"
#include        <stdio.h>
#include        <string.h>
#define NCODIGOS 41
#define MAXINT10 0x03FF

int codeanntyp(char letter);
int getstringfield(char** string, long i, mxChar *auxarray, long dim, int auxlength);


void mexFunction(
int	nlhs,
mxArray	*plhs[],
int	nrhs,
const mxArray	*prhs[]
)
{
/* Variable definition*/
    char *filename, *fieldname[6], zero='\0', *string;
    short nfields, subtypflag=0, chanflag=0, numflag=0, auxflag=0;
    mxArray *mxtime, *mxanntyp, *mxsubtyp, *mxchan, *mxnum, *mxaux;
    int  buflen, A, II, I, nbytes, auxlength;
    long dim=-1,  N, M;
    double *time, currtime;
    long i, pos;
    mxChar *anntyp;
    mxChar *subtyp;
    mxChar *chan, *num, *auxarray;
    char **aux;
    unsigned char defaultnumfield, currnumfield, defaultchanfield, currchanfield, defaultsubtyp, *defaultauxfield;
    unsigned char data[2], skip[4];
    
    FILE *fid;
    
/*---------------------------------------------------------------------*/
/*----------------------------------------------------------------------/
/	Check for proper number of arguments				/
/----------------------------------------------------------------------*/
    
    if (nrhs!=2)
    {
        mexErrMsgTxt ( "writeannot: two input arguments required" );
    }
    if ( nlhs > 1 )
    {
        mexErrMsgTxt ( "writeannot: No output arguments required" );
    }
    
    if (!mxIsChar(prhs[0]))
        mexErrMsgTxt("First argument must be a string.");
    else if (mxGetM(prhs[0])!=1)
        mexErrMsgTxt("First argument must be a string in a row vector.");
    
/* Get the length of the input string. */
    buflen = mxGetN(prhs[0])+1;
/* Allocate memory for input filename */
    filename = mxCalloc(buflen, sizeof(char));
/* Copy argument to variable */
    mxGetString(prhs[0], filename, buflen);
    
/* Checking annotation struct format*/
    if (!mxIsStruct(prhs[1]))
        mexErrMsgTxt("Second argument must be a struct.");
    if ((mxGetM(prhs[1])!=1)||(mxGetN(prhs[1])!=1))
        mexErrMsgTxt("Second argument must be a 1x1 struct");
    if  ( (nfields=mxGetNumberOfFields(prhs[1])) <2)
        mexErrMsgTxt("Annotation struct must have at least two fields");
    if  (nfields>6)
        mexErrMsgTxt("Annotation struct cannot have more than six fields");
    for (i=0;i<nfields;i++)
    {if (mxGetNumberOfDimensions(mxGetFieldByNumber(prhs[1],0,i))>2)
        mexErrMsgTxt("Fields should have less than 3 dimensions");
    if ((N=mxGetN(mxGetFieldByNumber(prhs[1],0,i)))!=1)
        if (strcmp(mxGetFieldNameByNumber(prhs[1],i),"aux")!=0)
            mexErrMsgTxt("Only aux field can have more than one column");
    M = mxGetM(mxGetFieldByNumber(prhs[1],0,i));
    if ((dim==-1)||(dim==M))
        dim = M;
    else
        mexErrMsgTxt("All fields must have the same number of Columns.\n If you want to pass an empty aux field, you must create a Nx1 empty matrix. ");
    }
    if (dim < 1) mexErrMsgTxt ("Struct must contain at least one annotation");
    mxtime = mxGetField(prhs[1],0,"time");
    mxanntyp = mxGetField(prhs[1],0,"anntyp");
    if ((mxtime == NULL)||(mxanntyp == NULL))
        mexErrMsgTxt("Annotation struct must have at least time and anntyp fields");
    
    mxsubtyp = mxGetField(prhs[1],0,"subtyp");
    
    if (mxsubtyp==NULL)
        mexWarnMsgTxt("subtyp field not found. Default: 0");
    
    mxchan = mxGetField(prhs[1],0,"chan");
    if (mxchan==NULL)
        mexWarnMsgTxt("chan field not found. Default: 0");
    mxnum = mxGetField(prhs[1],0,"num");
    if (mxnum==NULL)
        mexWarnMsgTxt("num field not found. Default: 0");
    mxaux = mxGetField(prhs[1],0,"aux");
    if (mxaux==NULL)
        mexWarnMsgTxt("aux field not found. Default: No auxiliary text");
    else auxlength = mxGetN(mxaux);
    
    
       /* Check data types and get data from fields */
    if (!mxIsNumeric(mxtime)) mexErrMsgTxt ("field \"time\" elements must be numeric");
    if (!mxIsChar(mxanntyp)) mexErrMsgTxt ("field \"anntyp\" elements must be characters");
    if ((mxsubtyp!=NULL) && (!mxIsChar(mxsubtyp))) mexErrMsgTxt ("field \"subtyp\" elements must be characters");
    if ((mxnum!=NULL)&&(!mxIsChar(mxnum))) mexErrMsgTxt ("field \"num\" elements must be characters");
    if ((mxchan!=NULL)&&(!mxIsChar(mxchan))) mexErrMsgTxt ("field \"chan\" elements must be characters");
    if ((mxaux!=NULL)&&(!mxIsChar(mxaux))) mexErrMsgTxt ("field \"aux\" elements must be characters");
    
 /* time = (double*) mxCalloc(dim, sizeof(double));  */               /* Allocates memory and initializes to zeros */
 /* anntyp = (mxChar*) mxCalloc(dim, sizeof(mxChar));
  subtyp = (mxChar*) mxCalloc(dim, sizeof(mxChar));
  chan = (mxChar*) mxCalloc(dim, sizeof(mxChar));
  num = (mxChar*) mxCalloc(dim, sizeof(mxChar));
  aux = (char**) mxCalloc(dim, sizeof(char*));*/
  /* Think about if these last lines are really necessary  */
    
    time = mxGetPr(mxtime);
    anntyp = (mxChar*) mxGetData(mxanntyp);
    if (mxsubtyp!=NULL)
        subtyp = (mxChar*) mxGetData(mxsubtyp);
    else subtypflag = 1;
    if (mxnum!=NULL)
        num = (mxChar*) mxGetData(mxnum);
    else  numflag=1;
    if (mxchan!=NULL)
        chan = (mxChar*) mxGetData(mxchan);
    else  chanflag=1;
    
  
    aux = (char**) mxCalloc(dim, sizeof(char*));
    auxflag = (mxaux==NULL);
    
    
/*-----------------------------------------------------------------------------*/
/*------------------------- Opening the file for writting            ----------*/
/*-----------------------------------------------------------------------------*/
    
    if ((fid = fopen(filename, "wb"))==NULL)
        mexErrMsgTxt("Unable to open annotation file");
    rewind(fid);
    
    
/***************** Variable Initialization**********************/
    i=0;           /* Current annotation number (index) */
    pos = 0;       /* Current position */
    defaultnumfield ='0';
    currnumfield = '0';
    defaultchanfield = '0';
    currchanfield = '0';
    defaultsubtyp = '0';
    defaultauxfield = (unsigned char *) &zero;
    
/**************** Writing of first annotation ************************/
    for (i=0;i<dim;i++) {
        if ((A = codeanntyp((char) anntyp[i]))==-1)
            mexErrMsgTxt("An unknown anntyp code has been found. \n If you use new annotation codes, please revise this program");
        if ((II=(long) time[i])<0) mexErrMsgTxt("A negative time value has been found");
        if ((i>0)&&(II<=pos)) mexErrMsgTxt("Two annotations have been found at the same time or badly ordered annotation");
        II = II - pos;
        pos +=II;
        if  (II > MAXINT10)
        /*SKIP */
        {
            data[1]= (char) (59<<2); data[0]= 0;
            fwrite(data,sizeof(char),2,fid);
            skip[0]=(char)((II&0x00ff0000)>>16); skip[1]= (char) ((II&0xff000000)>>24);
            skip[2]=(char)(II&0x000000ff); skip[3]=(char)((II&0x0000ff00)>>8);
            fwrite(skip, sizeof(char),4,fid);
            II=0;
        }
        
        data[1]=(char) ((A<<2)+((II&0x0300)>>8)); data[0]=(char) (II&0x00ff);
        fwrite(data,sizeof(char),2,fid);
        if (!subtypflag)
            if (subtyp[i]!=defaultsubtyp) {
                I = subtyp[i]-defaultsubtyp;
                data[1]= (char) (61<<2); data[0]= (char) I;
                fwrite(data, sizeof(char),2,fid);
            }
        if (!numflag)
            if (num[i]!=currnumfield) {
                I = num[i]-defaultnumfield;
                data[1]= (char) (60<<2); data[0]= (char) I;
                fwrite(data, sizeof(char),2,fid);
                currnumfield = num[i];
            }
        if (!chanflag)
            if (chan[i]!=currchanfield) {
                I = chan[i]-defaultchanfield;
                data[1]= (char) (62<<2); data[0]= (char) I;
                fwrite(data, sizeof(char),2,fid);
                currchanfield = chan[i];
            }
        
        if ((!auxflag)&&(auxlength)) 
        {
            auxarray = (mxChar*) mxGetData(mxaux);
            nbytes = getstringfield(&string, i, auxarray, dim, auxlength);
            if (nbytes)
            {
                if (nbytes > MAXINT10) mexErrMsgTxt ("Too long an auxiliary field");
                data[1] = (char) ((63<<2)+((nbytes&0x0300)>>8)); data[0]=(char)(nbytes&0x00ff);
                fwrite(data, sizeof(char),2,fid);
                if (nbytes & 1)  string[nbytes++]='\0';
                fwrite(string,sizeof(char),nbytes,fid);
            }
            mxFree(string);
        }
    }
/* Finish file with two zero bytes */
    data[0]=0;data[1]=0;
    fwrite(data,sizeof(char),2,fid);
    fclose(fid); 
                  
                  
                  
/****************** Free Memory   ******************/
    
    mxFree(filename);
    mxFree(aux);
    
}

int codeanntyp(char letter)
{const char* typecode = "NLRaVFJASEj/Q~ | sT*D\"=pB^t+u?![]en@xf()r";
short i=0, found =-1;
while ((found==-1)&&(i<NCODIGOS))
{if (typecode[i]==letter) found=++i;
i++;}
return found;
}

int getstringfield(char** string, long i, mxChar *array,long dim, int auxlength)
{int k=0;
*string = (char*) mxCalloc(auxlength+1, sizeof(char));
while ((k<auxlength)&&(array[k*dim+i]!='\0')&&(array[k*dim+i]!=' '))
{
    (*string)[k] = (char) array[k*dim+i];
    k++;
}
return k;
}
