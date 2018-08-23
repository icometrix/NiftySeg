#include "_seg_common.h"
#include "_seg_matrix.h"
#include "_seg_tools.h"

void Gaussian_Filter_Short_4D_old(segPrecisionTYPE *ShortData,
                                 int *S2L,
                                 int *L2S,
                                 segPrecisionTYPE gauss_std,
                                 ImageSize *CurrSizes,
                                 int class_with_CSF){

    long kernelsize=0;
    long kernelsizemin=(int)floorf(gauss_std*6.0);
    long kernelsizemax=(int)ceilf(gauss_std*6.0);

    if((kernelsizemin/2.0)==(double)(kernelsizemin/2) && kernelsizemin!=kernelsizemax){ // Which one is odd? kernelsizemin or kernelsizemax?
        kernelsize=kernelsizemax;}
    else if((kernelsizemax/2.0)==(double)(kernelsizemax/2) && kernelsizemin!=kernelsizemax){
        kernelsize=kernelsizemin;}
    else{
        kernelsize=kernelsizemin+1;}

    if(kernelsize<3){
        kernelsize=3;
    }

    long kernelshift=(int)floorf(kernelsize/2);
    segPrecisionTYPE GaussKernel [100]= {0};

    for(long i=0; i<kernelsize; i++){
        float kernelvalue=expf((float)(-0.5*powf((i-kernelshift)/gauss_std, 2)))/(sqrtf(2*3.14159265*powf(gauss_std, 2)));
        GaussKernel[i]=kernelvalue;
    }

    segPrecisionTYPE * Buffer= new segPrecisionTYPE [CurrSizes->numel]();
    segPrecisionTYPE * LongData= new segPrecisionTYPE [CurrSizes->numel]();


    long shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=(int)CurrSizes->xsize;
    shiftdirection[2]=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
    long dim_array[3];
    dim_array[0]=(int)CurrSizes->xsize;
    dim_array[1]=(int)CurrSizes->ysize;
    dim_array[2]=(int)CurrSizes->zsize;


    //Outside the mask is considered Pure CSF
    int outsiderangevalue[10];
    for(long i=0; i<10; i++){
        outsiderangevalue[i]=0;
    }
    outsiderangevalue[class_with_CSF]=1;

    for(long curr4d=0; curr4d<CurrSizes->numclass; curr4d++){ //For Each Class
        int current_4dShift_short=curr4d*CurrSizes->numelmasked;
        for(long index=0;index<(long)CurrSizes->numelmasked;index++){ //Copy Class to Buffer in LongFormat
            Buffer[S2L[index]]=ShortData[index+current_4dShift_short];
        }


        long xyzpos[3];
        for(long currentdirection=0;currentdirection<3;currentdirection++){ //Blur Buffer along each direction
            int index=0;
            for(xyzpos[2]=0;xyzpos[2]<(long)CurrSizes->zsize;xyzpos[2]++){
                for(xyzpos[1]=0;xyzpos[1]<(long)CurrSizes->ysize;xyzpos[1]++){
                    for(xyzpos[0]=0;xyzpos[0]<(long)CurrSizes->xsize;xyzpos[0]++){
                        segPrecisionTYPE tmpvalue=0.0f;
                        segPrecisionTYPE tmpkernelsum=0.0f;
                        LongData[index]=0.0f;
                        if(L2S[index]>=0){
                            for(long shift=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);shift<=((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))?(int)dim_array[currentdirection]-xyzpos[currentdirection]-1:kernelshift) ; shift++){
                                tmpvalue+=(L2S[index+shift*shiftdirection[currentdirection]]==-1)?GaussKernel[shift+kernelshift]*outsiderangevalue[curr4d]:GaussKernel[shift+kernelshift]*Buffer[index+shift*shiftdirection[currentdirection]];
                                tmpkernelsum+=GaussKernel[shift+kernelshift];
                            }
                            LongData[index]=tmpvalue/tmpkernelsum;
                        }
                        index++;
                    }
                }
            }
            if(currentdirection<2){
                for(long index2=0;index2<(long)CurrSizes->numel;index2++){
                    Buffer[index2]=LongData[index2];
                }
            }
        }

        for(long index=0;index<(long)CurrSizes->numelmasked;index++){ //Copy Class to Buffer in LongFormat
            ShortData[index+current_4dShift_short]=LongData[S2L[index]];
        }


    }
    delete [] LongData;
    delete [] Buffer;
    return;
}
