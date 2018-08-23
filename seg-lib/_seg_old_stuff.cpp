#include "_seg_common.h"
#include "_seg_matrix.h"
#include "_seg_tools.h"

#define numclass(I)         ((I)->nt)

void get_xyz_pow_int(segPrecisionTYPE xpos,
                     segPrecisionTYPE ypos,
                     segPrecisionTYPE zpos,
                     segPrecisionTYPE *currxpower,
                     segPrecisionTYPE *currypower,
                     segPrecisionTYPE *currzpower,
                     int maxorder){
    int order=1;
    currxpower[0]=1;
    currypower[0]=1;
    currzpower[0]=1;
    int orderminusone=0;
    int maxorderplusone=maxorder+1;
    while (order<maxorderplusone){
        currxpower[order]=currxpower[orderminusone]*xpos;
        currypower[order]=currypower[orderminusone]*ypos;
        currzpower[order]=currzpower[orderminusone]*zpos;
        order++;
        orderminusone++;
    }
    return;
}

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

nifti_image * Get_Bias_Corrected(float * BiasField,
                                 nifti_image * T1,
                                 char * filename,
                                 ImageSize * CurrSizes){

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);
    BiasCorrected->scl_inter=0;
    BiasCorrected->scl_slope=1;
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);



    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);

    for(long multispec=0;multispec<CurrSizes->usize;multispec++){

        BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*BiasCorrected->nvox];

        T1data = static_cast<segPrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*BiasCorrected->nvox];

        for( unsigned int i=0; i<BiasCorrected->nvox; i++){
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;

        for(long i=0; i<(long)CurrSizes->numel; i++){
            to_resize=exp((BiasField[i]+T1data[i])*0.693147181)-1;
            BiasCorrected_PTR[i]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
        }
    }
    return BiasCorrected;
}


int Gaussian_Filter_4D(segPrecisionTYPE * LongData,
                       segPrecisionTYPE gauss_std,
                       ImageSize * CurrSizes){

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

    segPrecisionTYPE GaussKernel [500]= {0};

    if(kernelsize>500){
        cout << "ERROR: Std is to large for Gaussian filter"<<endl;
        flush(cout);
        exit(9);
    }
    for(long i=0; i<kernelsize; i++){
        GaussKernel[i]=expf((float)(-0.5*powf((i-kernelshift)/gauss_std, 2)))/(sqrtf(2*3.14159265*powf(gauss_std, 2)));
    }


    segPrecisionTYPE * Buffer= new segPrecisionTYPE [CurrSizes->numel]();



    long shiftdirection[3];
    shiftdirection[0]=1;
    shiftdirection[1]=(int)CurrSizes->xsize;
    shiftdirection[2]=(int)CurrSizes->xsize*(int)CurrSizes->ysize;
    long dim_array[3];
    dim_array[0]=(int)CurrSizes->xsize;
    dim_array[1]=(int)CurrSizes->ysize;
    dim_array[2]=(int)CurrSizes->zsize;


    for(long curr4d=0; curr4d<CurrSizes->tsize; curr4d++){ //For Each Class

        long current_4dShift_short=curr4d*CurrSizes->numel;
        segPrecisionTYPE * longdataptr=&LongData[current_4dShift_short];
        for(long index=0;index<(long)CurrSizes->numel;index++){ //Copy Class to Buffer in LongFormat
            Buffer[index]=longdataptr[index];
        }


        long xyzpos[3];
        for(long currentdirection=0;currentdirection<3;currentdirection++){ //Blur Buffer along each direction
            long index=0;
            segPrecisionTYPE * longdataptr2=&longdataptr[0];
            for(xyzpos[2]=0;xyzpos[2]<(long)CurrSizes->zsize;xyzpos[2]++){
                for(xyzpos[1]=0;xyzpos[1]<(long)CurrSizes->ysize;xyzpos[1]++){
                    xyzpos[0]=0;
                    int kernelshiftminus=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);
                    int kernelshiftmplus=kernelshift;
                    for(xyzpos[0]=0;xyzpos[0]<CurrSizes->xsize;xyzpos[0]++){

                        if(xyzpos[currentdirection]<=kernelshift){
                            kernelshiftminus=((xyzpos[currentdirection]<kernelshift)?-xyzpos[currentdirection]:-kernelshift);
                        }
                        if((xyzpos[currentdirection]>=(dim_array[currentdirection]-kernelshift))){
                            kernelshiftmplus=dim_array[currentdirection]-xyzpos[currentdirection]-1;
                        }

                        segPrecisionTYPE * GaussKernelptr=&GaussKernel[kernelshiftminus+kernelshift];
                        segPrecisionTYPE * Bufferptr=&Buffer[index+kernelshiftminus*shiftdirection[currentdirection]];
                        segPrecisionTYPE tmpvalue=0.0f;
                        segPrecisionTYPE tmpkernelsum=0.0f;

                        for(int shift=kernelshiftminus;shift<=kernelshiftmplus ; shift++, GaussKernelptr++,Bufferptr+=shiftdirection[currentdirection]){
                            segPrecisionTYPE GaussKernelvar=(*GaussKernelptr);
                            tmpvalue+=GaussKernelvar*(*Bufferptr);
                            tmpkernelsum+=GaussKernelvar;
                        }

                        (*longdataptr2)=tmpvalue/tmpkernelsum;

                        index++;
                        longdataptr2++;
                    }
                }
            }
            if(currentdirection<2){
                for(long index2=0;index2<(long)CurrSizes->numel;index2++){
                    Buffer[index2]=longdataptr[index2];
                }
            }
        }
    }

    delete [] Buffer;
    return 0;
}

nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs,
                                      nifti_image * T1,
                                      char * filename,
                                      ImageSize * CurrSizes,
                                      int biasOrder){

    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * T1data = static_cast<segPrecisionTYPE *>(T1->data);

    nifti_image * BiasCorrected = nifti_copy_nim_info(T1);
    BiasCorrected->dim[0]=4;
    BiasCorrected->dim[4]=CurrSizes->usize;
    BiasCorrected->datatype=DT_FLOAT32;
    BiasCorrected->cal_max=(CurrSizes->rescale_max[0]);
    BiasCorrected->scl_inter=0;
    BiasCorrected->scl_slope=1;

    float * brainmask= new float [CurrSizes->numel];
    for(long i=0;i<(long)CurrSizes->numel;i++){
        brainmask[i]=T1data[i];
    }
    otsu(brainmask,NULL,CurrSizes);
    Dillate(brainmask,5,CurrSizes);
    Erosion(brainmask,4,CurrSizes);
    Gaussian_Filter_4D(brainmask, 3.0f, CurrSizes);


    nifti_set_filenames(BiasCorrected,filename,0,0);
    nifti_update_dims_from_array(BiasCorrected);
    nifti_datatype_sizes(BiasCorrected->datatype,&BiasCorrected->nbyper,&BiasCorrected->swapsize);
    BiasCorrected->data = (void *) calloc(BiasCorrected->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);

    float BiasField=0;
    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)CurrSizes->xsize);
    segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)CurrSizes->ysize);
    segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)CurrSizes->zsize);
    segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->xsize);
    segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->ysize);
    segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)CurrSizes->zsize);
    int ind=0;

    for(long multispec=0;multispec<CurrSizes->usize;multispec++){

        BiasCorrected_PTR = static_cast<segPrecisionTYPE *>(BiasCorrected->data);
        BiasCorrected_PTR = &BiasCorrected_PTR[multispec*CurrSizes->numel];
        T1data = static_cast<segPrecisionTYPE *>(T1->data);
        T1data = &T1data[multispec*CurrSizes->numel];

        float * BiasFieldCoefs_multispec = &BiasFieldCoefs[multispec*UsedBasisFunctions];


        for(long i=0; i<(long)CurrSizes->numel; i++){
            BiasCorrected_PTR[i]=0;
        }

        float to_resize=0;
        int index_full=0;
        for (int iz=0; iz<CurrSizes->zsize; iz++) {
            for (int iy=0; iy<CurrSizes->ysize; iy++) {
                for (int ix=0; ix<CurrSizes->xsize; ix++) {
                    BiasField=0.0f;
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                    ind=0;
                    for(long order=0; order<=biasOrder; order++){
                        for(long xorder=0; xorder<=order; xorder++){
                            for(long yorder=0; yorder<=(order-xorder); yorder++){
                                int zorder=order-yorder-xorder;
                                BiasField-=BiasFieldCoefs_multispec[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField*=brainmask[index_full];

                    to_resize=exp((BiasField+T1data[index_full])*0.693147181)-1;
                    BiasCorrected_PTR[index_full]=(to_resize*(CurrSizes->rescale_max[multispec]-CurrSizes->rescale_min[multispec])+CurrSizes->rescale_min[multispec]);
                    index_full++;
                }
            }
        }
    }
    return BiasCorrected;
}
