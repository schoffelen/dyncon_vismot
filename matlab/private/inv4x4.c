#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*declare variables*/
  const mwSize *dims;
  mwSize  *dimsout;
  mwIndex indx;
  int i, numdims;
  int numelin;
  mxClassID classid;
  double *input1r, *input1i, *output1r, *output1i;
  double x11,x12,x13,x14,x21,x22,x23,x24,x31,x32,x33,x34,x41,x42,x43,x44;
  double x11i,x12i,x13i,x14i,x21i,x22i,x23i,x24i,x31i,x32i,x33i,x34i,x41i,x42i,x43i,x44i;
  double tmp1,tmp2,tmp3,tmp4,tmp1i,tmp2i,tmp3i,tmp4i;
  double D,Di,Dabs;

  /*figure out the classid*/
  classid = mxGetClassID(prhs[0]);
     
  /*check inputs*/
  if (nrhs>1)
    mexErrMsgTxt("Too many input arguments");
  
  /*associate inputs*/
  input1r = mxGetData(prhs[0]);
  input1i = mxGetImagData(prhs[0]);
  
  /*figure out dimension info and number of elements*/
  dims    = mxGetDimensions(prhs[0]);
  numdims = mxGetNumberOfDimensions(prhs[0]);
  numelin = mxGetNumberOfElements(prhs[0]);
  
  dimsout    = mxMalloc(numdims * sizeof(mwSize));
  for (i=0; i<numdims; i++)
  {
    dimsout[i] = dims[i];
  }
  
  /*associate output*/
  if (input1i == NULL)
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxREAL);
    output1r = mxGetData(plhs[0]);
  }
  else
  {
    plhs[0]  = mxCreateNumericArray(numdims, dimsout, classid, mxCOMPLEX);
    output1r = mxGetData(plhs[0]);
    output1i = mxGetImagData(plhs[0]);
  }
  
  /* do the computation*/
  if (input1i == NULL)
  {  
    for (i=0; i<numelin/16; i++)
    {
      x11 = input1r[i*16   ];
      x21 = input1r[i*16+1 ];
      x31 = input1r[i*16+2 ];
      x41 = input1r[i*16+3 ];
      x12 = input1r[i*16+4 ];
      x22 = input1r[i*16+5 ];
      x32 = input1r[i*16+6 ];
      x42 = input1r[i*16+7 ];
      x13 = input1r[i*16+8 ];
      x23 = input1r[i*16+9 ];
      x33 = input1r[i*16+10];
      x43 = input1r[i*16+11];
      x14 = input1r[i*16+12];
      x24 = input1r[i*16+13];
      x34 = input1r[i*16+14];
      x44 = input1r[i*16+15];
      
      tmp1 =  ( (x22*(x33*x44-x34*x43)-x23*(x32*x44-x42*x34)+x24*(x32*x43-x33*x42)) );
      tmp2 = -( (x21*(x33*x44-x34*x43)-x23*(x31*x44-x41*x34)+x24*(x31*x43-x33*x41)) );
      tmp3 =  ( (x21*(x32*x44-x34*x42)-x22*(x31*x44-x41*x34)+x24*(x31*x42-x32*x41)) );
      tmp4 = -( (x21*(x32*x43-x33*x42)-x22*(x31*x43-x41*x33)+x23*(x31*x42-x32*x41)) );

      D    = tmp1*x11+tmp2*x12+tmp3*x13+tmp4*x14;

      output1r[i*16   ] = tmp1/D;
      output1r[i*16+1 ] = tmp2/D;
      output1r[i*16+2 ] = tmp3/D;
      output1r[i*16+3 ] = tmp4/D;
      
      output1r[i*16+4 ] = -(x12*x33*x44-x12*x34*x43-x13*x32*x44+x13*x42*x34+x14*x32*x43-x14*x33*x42)/D;
      output1r[i*16+5 ] =  (x11*x33*x44-x11*x34*x43-x13*x31*x44+x13*x41*x34+x14*x31*x43-x14*x33*x41)/D;
      output1r[i*16+6 ] = -(x11*x32*x44-x11*x34*x42-x12*x31*x44+x12*x41*x34+x14*x31*x42-x14*x32*x41)/D;
      output1r[i*16+7 ] =  (x11*x32*x43-x11*x33*x42-x12*x31*x43+x12*x41*x33+x13*x31*x42-x13*x32*x41)/D;

      output1r[i*16+8 ] =  (x12*x23*x44-x12*x24*x43-x13*x22*x44+x13*x42*x24+x14*x22*x43-x14*x23*x42)/D;
      output1r[i*16+9 ] = -(x11*x23*x44-x11*x24*x43-x13*x21*x44+x13*x41*x24+x14*x21*x43-x14*x23*x41)/D;
      output1r[i*16+10] =  (x11*x22*x44-x11*x24*x42-x12*x21*x44+x12*x41*x24+x14*x21*x42-x14*x22*x41)/D;
      output1r[i*16+11] = -(x11*x22*x43-x11*x23*x42-x12*x21*x43+x12*x41*x23+x13*x21*x42-x13*x22*x41)/D;
      
      output1r[i*16+12] = -(x12*x23*x34-x12*x24*x33-x13*x22*x34+x13*x32*x24+x14*x22*x33-x14*x23*x32)/D;
      output1r[i*16+13] =  (x11*x23*x34-x11*x24*x33-x13*x21*x34+x13*x31*x24+x14*x21*x33-x14*x23*x31)/D;
      output1r[i*16+14] = -(x11*x22*x34-x11*x24*x32-x12*x21*x34+x12*x31*x24+x14*x21*x32-x14*x22*x31)/D;
      output1r[i*16+15] =  (x11*x22*x33-x11*x23*x32-x12*x21*x33+x12*x31*x23+x13*x21*x32-x13*x22*x31)/D;
    }
    return;
  }
  else
  {  
    for (i=0; i<numelin/16; i++)
    {
      x11 = input1r[i*16   ];
      x21 = input1r[i*16+1 ];
      x31 = input1r[i*16+2 ];
      x41 = input1r[i*16+3 ];
      x12 = input1r[i*16+4 ];
      x22 = input1r[i*16+5 ];
      x32 = input1r[i*16+6 ];
      x42 = input1r[i*16+7 ];
      x13 = input1r[i*16+8 ];
      x23 = input1r[i*16+9 ];
      x33 = input1r[i*16+10];
      x43 = input1r[i*16+11];
      x14 = input1r[i*16+12];
      x24 = input1r[i*16+13];
      x34 = input1r[i*16+14];
      x44 = input1r[i*16+15];
      
      x11i = input1i[i*16   ];
      x21i = input1i[i*16+1 ];
      x31i = input1i[i*16+2 ];
      x41i = input1i[i*16+3 ];
      x12i = input1i[i*16+4 ];
      x22i = input1i[i*16+5 ];
      x32i = input1i[i*16+6 ];
      x42i = input1i[i*16+7 ];
      x13i = input1i[i*16+8 ];
      x23i = input1i[i*16+9 ];
      x33i = input1i[i*16+10];
      x43i = input1i[i*16+11];
      x14i = input1i[i*16+12];
      x24i = input1i[i*16+13];
      x34i = input1i[i*16+14];
      x44i = input1i[i*16+15];
      
      tmp1 =  ( (x22*(x33*x44-x34*x43)-x23*(x32*x44-x42*x34)+x24*(x32*x43-x33*x42)) - (x22i*(x33i*x44-x34i*x43)-x23i*(x32i*x44-x42i*x34)+x24i*(x32i*x43-x33i*x42)) - (x22i*(x33*x44i-x34*x43i)-x23i*(x32*x44i-x42*x34i)+x24i*(x32*x43i-x33*x42i)) - (x22*(x33i*x44i-x34i*x43i)-x23*(x32i*x44i-x42i*x34i)+x24*(x32i*x43i-x33i*x42i)));
      tmp2 = -( (x21*(x33*x44-x34*x43)-x23*(x31*x44-x41*x34)+x24*(x31*x43-x33*x41)) - (x21i*(x33i*x44-x34i*x43)-x23i*(x31i*x44-x41i*x34)+x24i*(x31i*x43-x33i*x41)) - (x21i*(x33*x44i-x34*x43i)-x23i*(x31*x44i-x41*x34i)+x24i*(x31*x43i-x33*x41i)) - (x21*(x33i*x44i-x34i*x43i)-x23*(x31i*x44i-x41i*x34i)+x24*(x31i*x43i-x33i*x41i)));
      tmp3 =  ( (x21*(x32*x44-x34*x42)-x22*(x31*x44-x41*x34)+x24*(x31*x42-x32*x41)) - (x21i*(x32i*x44-x34i*x42)-x22i*(x31i*x44-x41i*x34)+x24i*(x31i*x42-x32i*x41)) - (x21i*(x32*x44i-x34*x42i)-x22i*(x31*x44i-x41*x34i)+x24i*(x31*x42i-x32*x41i)) - (x21*(x32i*x44i-x34i*x42i)-x22*(x31i*x44i-x41i*x34i)+x24*(x31i*x42i-x32i*x41i)));
      tmp4 = -( (x21*(x32*x43-x33*x42)-x22*(x31*x43-x41*x33)+x23*(x31*x42-x32*x41)) - (x21i*(x32i*x43-x33i*x42)-x22i*(x31i*x43-x41i*x33)+x23i*(x31i*x42-x32i*x41)) - (x21i*(x32*x43i-x33*x42i)-x22i*(x31*x43i-x41*x33i)+x23i*(x31*x42i-x32*x41i)) - (x21*(x32i*x43i-x33i*x42i)-x22*(x31i*x43i-x41i*x33i)+x23*(x31i*x42i-x32i*x41i)));
      
      tmp1i =  (-(x22i*(x33i*x44i-x34i*x43i)-x23i*(x32i*x44i-x42i*x34i)+x24i*(x32i*x43i-x33i*x42i)) + (x22*(x33*x44i-x34*x43i)-x23*(x32*x44i-x42*x34i)+x24*(x32*x43i-x33*x42i)) + (x22*(x33i*x44-x34i*x43)-x23*(x32i*x44-x42i*x34)+x24*(x32i*x43-x33i*x42)) + (x22i*(x33*x44-x34*x43)-x23i*(x32*x44-x42*x34)+x24i*(x32*x43-x33*x42)));
      tmp2i = -(-(x21i*(x33i*x44i-x34i*x43i)-x23i*(x31i*x44i-x41i*x34i)+x24i*(x31i*x43i-x33i*x41i)) + (x21*(x33*x44i-x34*x43i)-x23*(x31*x44i-x41*x34i)+x24*(x31*x43i-x33*x41i)) + (x21*(x33i*x44-x34i*x43)-x23*(x31i*x44-x41i*x34)+x24*(x31i*x43-x33i*x41)) + (x21i*(x33*x44-x34*x43)-x23i*(x31*x44-x41*x34)+x24i*(x31*x43-x33*x41)));
      tmp3i =  (-(x21i*(x32i*x44i-x34i*x42i)-x22i*(x31i*x44i-x41i*x34i)+x24i*(x31i*x42i-x32i*x41i)) + (x21*(x32*x44i-x34*x42i)-x22*(x31*x44i-x41*x34i)+x24*(x31*x42i-x32*x41i)) + (x21*(x32i*x44-x34i*x42)-x22*(x31i*x44-x41i*x34)+x24*(x31i*x42-x32i*x41)) + (x21i*(x32*x44-x34*x42)-x22i*(x31*x44-x41*x34)+x24i*(x31*x42-x32*x41)));
      tmp4i = -(-(x21i*(x32i*x43i-x33i*x42i)-x22i*(x31i*x43i-x41i*x33i)+x23i*(x31i*x42i-x32i*x41i)) + (x21*(x32*x43i-x33*x42i)-x22*(x31*x43i-x41*x33i)+x23*(x31*x42i-x32*x41i)) + (x21*(x32i*x43-x33i*x42)-x22*(x31i*x43-x41i*x33)+x23*(x31i*x42-x32i*x41)) + (x21i*(x32*x43-x33*x42)-x22i*(x31*x43-x41*x33)+x23i*(x31*x42-x32*x41)));

      D    = tmp1*x11+tmp2*x12+tmp3*x13+tmp4*x14-tmp1i*x11i-tmp2i*x12i-tmp3i*x13i-tmp4i*x14i;
      Di   = tmp1*x11i+tmp2*x12i+tmp3*x13i+tmp4*x14i+tmp1i*x11+tmp2i*x12+tmp3i*x13+tmp4i*x14;
      Dabs = D*D+Di*Di;

      output1r[i*16   ] = (tmp1*D+tmp1i*Di)/Dabs;
      output1r[i*16+1 ] = (tmp2*D+tmp2i*Di)/Dabs;
      output1r[i*16+2 ] = (tmp3*D+tmp3i*Di)/Dabs;
      output1r[i*16+3 ] = (tmp4*D+tmp4i*Di)/Dabs;
      
      output1i[i*16   ] = (tmp1i*D-tmp1*Di)/Dabs;
      output1i[i*16+1 ] = (tmp2i*D-tmp2*Di)/Dabs;
      output1i[i*16+2 ] = (tmp3i*D-tmp3*Di)/Dabs;
      output1i[i*16+3 ] = (tmp4i*D-tmp4*Di)/Dabs;
      
      tmp1 = -( (x12*(x33*x44-x34*x43)-x13*(x32*x44-x42*x34)+x14*(x32*x43-x33*x42)) - (x12i*(x33i*x44-x34i*x43)-x13i*(x32i*x44-x42i*x34)+x14i*(x32i*x43-x33i*x42)) - (x12i*(x33*x44i-x34*x43i)-x13i*(x32*x44i-x42*x34i)+x14i*(x32*x43i-x33*x42i)) - (x12*(x33i*x44i-x34i*x43i)-x13*(x32i*x44i-x42i*x34i)+x14*(x32i*x43i-x33i*x42i)));
      tmp2 =  ( (x11*(x33*x44-x34*x43)-x13*(x31*x44-x41*x34)+x14*(x31*x43-x33*x41)) - (x11i*(x33i*x44-x34i*x43)-x13i*(x31i*x44-x41i*x34)+x14i*(x31i*x43-x33i*x41)) - (x11i*(x33*x44i-x34*x43i)-x13i*(x31*x44i-x41*x34i)+x14i*(x31*x43i-x33*x41i)) - (x11*(x33i*x44i-x34i*x43i)-x13*(x31i*x44i-x41i*x34i)+x14*(x31i*x43i-x33i*x41i)));
      tmp3 = -( (x11*(x32*x44-x34*x42)-x12*(x31*x44-x41*x34)+x14*(x31*x42-x32*x41)) - (x11i*(x32i*x44-x34i*x42)-x12i*(x31i*x44-x41i*x34)+x14i*(x31i*x42-x32i*x41)) - (x11i*(x32*x44i-x34*x42i)-x12i*(x31*x44i-x41*x34i)+x14i*(x31*x42i-x32*x41i)) - (x11*(x32i*x44i-x34i*x42i)-x12*(x31i*x44i-x41i*x34i)+x14*(x31i*x42i-x32i*x41i)));
      tmp4 =  ( (x11*(x32*x43-x33*x42)-x12*(x31*x43-x41*x33)+x13*(x31*x42-x32*x41)) - (x11i*(x32i*x43-x33i*x42)-x12i*(x31i*x43-x41i*x33)+x13i*(x31i*x42-x32i*x41)) - (x11i*(x32*x43i-x33*x42i)-x12i*(x31*x43i-x41*x33i)+x13i*(x31*x42i-x32*x41i)) - (x11*(x32i*x43i-x33i*x42i)-x12*(x31i*x43i-x41i*x33i)+x13*(x31i*x42i-x32i*x41i)));
      
      tmp1i = -(-(x12i*(x33i*x44i-x34i*x43i)-x13i*(x32i*x44i-x42i*x34i)+x14i*(x32i*x43i-x33i*x42i)) + (x12*(x33*x44i-x34*x43i)-x13*(x32*x44i-x42*x34i)+x14*(x32*x43i-x33*x42i)) + (x12*(x33i*x44-x34i*x43)-x13*(x32i*x44-x42i*x34)+x14*(x32i*x43-x33i*x42)) + (x12i*(x33*x44-x34*x43)-x13i*(x32*x44-x42*x34)+x14i*(x32*x43-x33*x42)));
      tmp2i =  (-(x11i*(x33i*x44i-x34i*x43i)-x13i*(x31i*x44i-x41i*x34i)+x14i*(x31i*x43i-x33i*x41i)) + (x11*(x33*x44i-x34*x43i)-x13*(x31*x44i-x41*x34i)+x14*(x31*x43i-x33*x41i)) + (x11*(x33i*x44-x34i*x43)-x13*(x31i*x44-x41i*x34)+x14*(x31i*x43-x33i*x41)) + (x11i*(x33*x44-x34*x43)-x13i*(x31*x44-x41*x34)+x14i*(x31*x43-x33*x41)));
      tmp3i = -(-(x11i*(x32i*x44i-x34i*x42i)-x12i*(x31i*x44i-x41i*x34i)+x14i*(x31i*x42i-x32i*x41i)) + (x11*(x32*x44i-x34*x42i)-x12*(x31*x44i-x41*x34i)+x14*(x31*x42i-x32*x41i)) + (x11*(x32i*x44-x34i*x42)-x12*(x31i*x44-x41i*x34)+x14*(x31i*x42-x32i*x41)) + (x11i*(x32*x44-x34*x42)-x12i*(x31*x44-x41*x34)+x14i*(x31*x42-x32*x41)));
      tmp4i =  (-(x11i*(x32i*x43i-x33i*x42i)-x12i*(x31i*x43i-x41i*x33i)+x13i*(x31i*x42i-x32i*x41i)) + (x11*(x32*x43i-x33*x42i)-x12*(x31*x43i-x41*x33i)+x13*(x31*x42i-x32*x41i)) + (x11*(x32i*x43-x33i*x42)-x12*(x31i*x43-x41i*x33)+x13*(x31i*x42-x32i*x41)) + (x11i*(x32*x43-x33*x42)-x12i*(x31*x43-x41*x33)+x13i*(x31*x42-x32*x41)));
      
      output1r[i*16+4 ] = (tmp1*D+tmp1i*Di)/Dabs;
      output1r[i*16+5 ] = (tmp2*D+tmp2i*Di)/Dabs;
      output1r[i*16+6 ] = (tmp3*D+tmp3i*Di)/Dabs;
      output1r[i*16+7 ] = (tmp4*D+tmp4i*Di)/Dabs;
      
      output1i[i*16+4 ] = (tmp1i*D-tmp1*Di)/Dabs;
      output1i[i*16+5 ] = (tmp2i*D-tmp2*Di)/Dabs;
      output1i[i*16+6 ] = (tmp3i*D-tmp3*Di)/Dabs;
      output1i[i*16+7 ] = (tmp4i*D-tmp4*Di)/Dabs;
      
      tmp1 =  ( (x12*(x23*x44-x24*x43)-x13*(x22*x44-x42*x24)+x14*(x22*x43-x23*x42)) - (x12i*(x23i*x44-x24i*x43)-x13i*(x22i*x44-x42i*x24)+x14i*(x22i*x43-x23i*x42)) - (x12i*(x23*x44i-x24*x43i)-x13i*(x22*x44i-x42*x24i)+x14i*(x22*x43i-x23*x42i)) - (x12*(x23i*x44i-x24i*x43i)-x13*(x22i*x44i-x42i*x24i)+x14*(x22i*x43i-x23i*x42i)));
      tmp2 = -( (x11*(x23*x44-x24*x43)-x13*(x21*x44-x41*x24)+x14*(x21*x43-x23*x41)) - (x11i*(x23i*x44-x24i*x43)-x13i*(x21i*x44-x41i*x24)+x14i*(x21i*x43-x23i*x41)) - (x11i*(x23*x44i-x24*x43i)-x13i*(x21*x44i-x41*x24i)+x14i*(x21*x43i-x23*x41i)) - (x11*(x23i*x44i-x24i*x43i)-x13*(x21i*x44i-x41i*x24i)+x14*(x21i*x43i-x23i*x41i)));
      tmp3 =  ( (x11*(x22*x44-x24*x42)-x12*(x21*x44-x41*x24)+x14*(x21*x42-x22*x41)) - (x11i*(x22i*x44-x24i*x42)-x12i*(x21i*x44-x41i*x24)+x14i*(x21i*x42-x22i*x41)) - (x11i*(x22*x44i-x24*x42i)-x12i*(x21*x44i-x41*x24i)+x14i*(x21*x42i-x22*x41i)) - (x11*(x22i*x44i-x24i*x42i)-x12*(x21i*x44i-x41i*x24i)+x14*(x21i*x42i-x22i*x41i)));
      tmp4 = -( (x11*(x22*x43-x23*x42)-x12*(x21*x43-x41*x23)+x13*(x21*x42-x22*x41)) - (x11i*(x22i*x43-x23i*x42)-x12i*(x21i*x43-x41i*x23)+x13i*(x21i*x42-x22i*x41)) - (x11i*(x22*x43i-x23*x42i)-x12i*(x21*x43i-x41*x23i)+x13i*(x21*x42i-x22*x41i)) - (x11*(x22i*x43i-x23i*x42i)-x12*(x21i*x43i-x41i*x23i)+x13*(x21i*x42i-x22i*x41i)));
      
      tmp1i =  (-(x12i*(x23i*x44i-x24i*x43i)-x13i*(x22i*x44i-x42i*x24i)+x14i*(x22i*x43i-x23i*x42i)) + (x12*(x23*x44i-x24*x43i)-x13*(x22*x44i-x42*x24i)+x14*(x22*x43i-x23*x42i)) + (x12*(x23i*x44-x24i*x43)-x13*(x22i*x44-x42i*x24)+x14*(x22i*x43-x23i*x42)) + (x12i*(x23*x44-x24*x43)-x13i*(x22*x44-x42*x24)+x14i*(x22*x43-x23*x42)));
      tmp2i = -(-(x11i*(x23i*x44i-x24i*x43i)-x13i*(x21i*x44i-x41i*x24i)+x14i*(x21i*x43i-x23i*x41i)) + (x11*(x23*x44i-x24*x43i)-x13*(x21*x44i-x41*x24i)+x14*(x21*x43i-x23*x41i)) + (x11*(x23i*x44-x24i*x43)-x13*(x21i*x44-x41i*x24)+x14*(x21i*x43-x23i*x41)) + (x11i*(x23*x44-x24*x43)-x13i*(x21*x44-x41*x24)+x14i*(x21*x43-x23*x41)));
      tmp3i =  (-(x11i*(x22i*x44i-x24i*x42i)-x12i*(x21i*x44i-x41i*x24i)+x14i*(x21i*x42i-x22i*x41i)) + (x11*(x22*x44i-x24*x42i)-x12*(x21*x44i-x41*x24i)+x14*(x21*x42i-x22*x41i)) + (x11*(x22i*x44-x24i*x42)-x12*(x21i*x44-x41i*x24)+x14*(x21i*x42-x22i*x41)) + (x11i*(x22*x44-x24*x42)-x12i*(x21*x44-x41*x24)+x14i*(x21*x42-x22*x41)));
      tmp4i = -(-(x11i*(x22i*x43i-x23i*x42i)-x12i*(x21i*x43i-x41i*x23i)+x13i*(x21i*x42i-x22i*x41i)) + (x11*(x22*x43i-x23*x42i)-x12*(x21*x43i-x41*x23i)+x13*(x21*x42i-x22*x41i)) + (x11*(x22i*x43-x23i*x42)-x12*(x21i*x43-x41i*x23)+x13*(x21i*x42-x22i*x41)) + (x11i*(x22*x43-x23*x42)-x12i*(x21*x43-x41*x23)+x13i*(x21*x42-x22*x41)));
      
      output1r[i*16+8 ] = (tmp1*D+tmp1i*Di)/Dabs;
      output1r[i*16+9 ] = (tmp2*D+tmp2i*Di)/Dabs;
      output1r[i*16+10] = (tmp3*D+tmp3i*Di)/Dabs;
      output1r[i*16+11] = (tmp4*D+tmp4i*Di)/Dabs;
      
      output1i[i*16+8 ] = (tmp1i*D-tmp1*Di)/Dabs;
      output1i[i*16+9 ] = (tmp2i*D-tmp2*Di)/Dabs;
      output1i[i*16+10] = (tmp3i*D-tmp3*Di)/Dabs;
      output1i[i*16+11] = (tmp4i*D-tmp4*Di)/Dabs;
      
      tmp1 = -( (x12*(x23*x34-x24*x33)-x13*(x22*x34-x32*x24)+x14*(x22*x33-x23*x32)) - (x12i*(x23i*x34-x24i*x33)-x13i*(x22i*x34-x32i*x24)+x14i*(x22i*x33-x23i*x32)) - (x12i*(x23*x34i-x24*x33i)-x13i*(x22*x34i-x32*x24i)+x14i*(x22*x33i-x23*x32i)) - (x12*(x23i*x34i-x24i*x33i)-x13*(x22i*x34i-x32i*x24i)+x14*(x22i*x33i-x23i*x32i)));
      tmp2 =  ( (x11*(x23*x34-x24*x33)-x13*(x21*x34-x31*x24)+x14*(x21*x33-x23*x31)) - (x11i*(x23i*x34-x24i*x33)-x13i*(x21i*x34-x31i*x24)+x14i*(x21i*x33-x23i*x31)) - (x11i*(x23*x34i-x24*x33i)-x13i*(x21*x34i-x31*x24i)+x14i*(x21*x33i-x23*x31i)) - (x11*(x23i*x34i-x24i*x33i)-x13*(x21i*x34i-x31i*x24i)+x14*(x21i*x33i-x23i*x31i)));
      tmp3 = -( (x11*(x22*x34-x24*x32)-x12*(x21*x34-x31*x24)+x14*(x21*x32-x22*x31)) - (x11i*(x22i*x34-x24i*x32)-x12i*(x21i*x34-x31i*x24)+x14i*(x21i*x32-x22i*x31)) - (x11i*(x22*x34i-x24*x32i)-x12i*(x21*x34i-x31*x24i)+x14i*(x21*x32i-x22*x31i)) - (x11*(x22i*x34i-x24i*x32i)-x12*(x21i*x34i-x31i*x24i)+x14*(x21i*x32i-x22i*x31i)));
      tmp4 =  ( (x11*(x22*x33-x23*x32)-x12*(x21*x33-x31*x23)+x13*(x21*x32-x22*x31)) - (x11i*(x22i*x33-x23i*x32)-x12i*(x21i*x33-x31i*x23)+x13i*(x21i*x32-x22i*x31)) - (x11i*(x22*x33i-x23*x32i)-x12i*(x21*x33i-x31*x23i)+x13i*(x21*x32i-x22*x31i)) - (x11*(x22i*x33i-x23i*x32i)-x12*(x21i*x33i-x31i*x23i)+x13*(x21i*x32i-x22i*x31i)));
      
      tmp1i = -(-(x12i*(x23i*x34i-x24i*x33i)-x13i*(x22i*x34i-x32i*x24i)+x14i*(x22i*x33i-x23i*x32i)) + (x12*(x23*x34i-x24*x33i)-x13*(x22*x34i-x32*x24i)+x14*(x22*x33i-x23*x32i)) + (x12*(x23i*x34-x24i*x33)-x13*(x22i*x34-x32i*x24)+x14*(x22i*x33-x23i*x32)) + (x12i*(x23*x34-x24*x33)-x13i*(x22*x34-x32*x24)+x14i*(x22*x33-x23*x32)));
      tmp2i =  (-(x11i*(x23i*x34i-x24i*x33i)-x13i*(x21i*x34i-x31i*x24i)+x14i*(x21i*x33i-x23i*x31i)) + (x11*(x23*x34i-x24*x33i)-x13*(x21*x34i-x31*x24i)+x14*(x21*x33i-x23*x31i)) + (x11*(x23i*x34-x24i*x33)-x13*(x21i*x34-x31i*x24)+x14*(x21i*x33-x23i*x31)) + (x11i*(x23*x34-x24*x33)-x13i*(x21*x34-x31*x24)+x14i*(x21*x33-x23*x31)));
      tmp3i = -(-(x11i*(x22i*x34i-x24i*x32i)-x12i*(x21i*x34i-x31i*x24i)+x14i*(x21i*x32i-x22i*x31i)) + (x11*(x22*x34i-x24*x32i)-x12*(x21*x34i-x31*x24i)+x14*(x21*x32i-x22*x31i)) + (x11*(x22i*x34-x24i*x32)-x12*(x21i*x34-x31i*x24)+x14*(x21i*x32-x22i*x31)) + (x11i*(x22*x34-x24*x32)-x12i*(x21*x34-x31*x24)+x14i*(x21*x32-x22*x31)));
      tmp4i =  (-(x11i*(x22i*x33i-x23i*x32i)-x12i*(x21i*x33i-x31i*x23i)+x13i*(x21i*x32i-x22i*x31i)) + (x11*(x22*x33i-x23*x32i)-x12*(x21*x33i-x31*x23i)+x13*(x21*x32i-x22*x31i)) + (x11*(x22i*x33-x23i*x32)-x12*(x21i*x33-x31i*x23)+x13*(x21i*x32-x22i*x31)) + (x11i*(x22*x33-x23*x32)-x12i*(x21*x33-x31*x23)+x13i*(x21*x32-x22*x31)));
      
      output1r[i*16+12] = (tmp1*D+tmp1i*Di)/Dabs;
      output1r[i*16+13] = (tmp2*D+tmp2i*Di)/Dabs;
      output1r[i*16+14] = (tmp3*D+tmp3i*Di)/Dabs;
      output1r[i*16+15] = (tmp4*D+tmp4i*Di)/Dabs;
      
      output1i[i*16+12] = (tmp1i*D-tmp1*Di)/Dabs;
      output1i[i*16+13] = (tmp2i*D-tmp2*Di)/Dabs;
      output1i[i*16+14] = (tmp3i*D-tmp3*Di)/Dabs;
      output1i[i*16+15] = (tmp4i*D-tmp4*Di)/Dabs;
    }
    return;
  }
}
