#include "_seg_common.h"
#include "_seg_matrix.h"
#include "_seg_tools.h"

#define colsize(I)          ((I)->nx)
#define rowsize(I)          ((I)->ny)
#define depth(I)            ((I)->nz)
#define numclass(I)         ((I)->nt)
#define numelem(I)          ((I)->nx*(I)->ny*(I)->nz)
#define numelemtotal(I)     ((I)->nvox)

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

void get_xy_pow_int(segPrecisionTYPE xpos,
                    segPrecisionTYPE ypos,
                    segPrecisionTYPE *currxpower,
                    segPrecisionTYPE *currypower,
                    int maxorder){
    int order=1;
    currxpower[0]=1;
    currypower[0]=1;
    int orderminusone=0;
    int maxorderplusone=maxorder+1;
    while (order<maxorderplusone){
        currxpower[order]=currxpower[orderminusone]*xpos;
        currypower[order]=currypower[orderminusone]*ypos;
        order++;
        orderminusone++;
    }
    return;
}


segPrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                                 nifti_image * Priors,
                                                 long numclass,
                                                 bool PV_ON)
{
    register long numel=(int)(rowsize(Mask)*colsize(Mask)*depth(Mask));
    register long numel_masked=0;

    bool * Maskptrtmp = static_cast<bool *> (Mask->data);;
    for (long i=0; i<numel; i++, Maskptrtmp++) {
        *Maskptrtmp?numel_masked++:0;
    }
    int pluspv=(int)(PV_ON)*2;

    segPrecisionTYPE * Expec = new segPrecisionTYPE [numel_masked*(numclass+pluspv)] ();
    segPrecisionTYPE * tempExpec= (segPrecisionTYPE *) Expec;
    segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<numclass;cl++){
        Maskptrtmp = static_cast<bool *> (Mask->data);;
        for (int i=numel; i--; Maskptrtmp++,PriorPTR++) {
            if(*Maskptrtmp){
                *tempExpec = *PriorPTR;
                tempExpec++;
            }
        }
    }

    return Expec;
}

segPrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                            long numclass,
                                            bool PV_ON)
{
    register long numel=(int)(rowsize(Priors)*colsize(Priors)*depth(Priors));
    long pluspv=(int)(PV_ON)*2;
    segPrecisionTYPE * Expec = new segPrecisionTYPE [numel*(numclass+pluspv)] ();
    segPrecisionTYPE * Expec_PTR= Expec;
    segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(Priors->data);
    for(long cl=0; cl<numclass;cl++){
        for (int i=numel; i--; PriorPTR++,Expec_PTR++) {
            *Expec_PTR = *PriorPTR;
        }
    }
    return Expec;
}


int Gaussian_Filter_Short_4D(segPrecisionTYPE * ShortData,
                             int * S2L,
                             int * L2S,
                             segPrecisionTYPE gauss_std,
                             ImageSize * CurrSizes,
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
    return 1;
}

int PriorWeight_mask(float * ShortPrior,nifti_image * Priors, float * Expec,float GaussKernelSize,float RelaxFactor, int * S2L, int * L2S,ImageSize * CurrSizes,int verbose_level){

    for(long i=0; i<(CurrSizes->numclass*CurrSizes->numelmasked);i++)ShortPrior[i]=Expec[i];
    Gaussian_Filter_Short_4D(ShortPrior,S2L,L2S,GaussKernelSize,CurrSizes,CSFclass);
    float * PriorsPtr = static_cast<float *>(Priors->data);
    long currindex=0;
    for(long k=0; k<CurrSizes->numclass; k++){
        currindex=k*CurrSizes->numelmasked;
        for(long i=0; i<(CurrSizes->numelmasked);i++){
            ShortPrior[currindex]*=(1-RelaxFactor);
            ShortPrior[currindex]+=(RelaxFactor)*PriorsPtr[S2L[i]+k*CurrSizes->numel];
            currindex++;
        }
    }

    return 1;
}



int Normalize_NaN_Priors_mask(nifti_image * Priors,
                              nifti_image * Mask,
                              bool verbose)
{
    register int numel = Mask->nvox;
    register int ups=0;
    register int good=0;
    if(verbose>0){
        cout<< "Normalizing Priors" << endl;
    }
    if(Mask->datatype==DT_BINARY){
        if(Priors->datatype==NIFTI_TYPE_FLOAT32){
            segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(Priors->data);
            bool * brainmaskptr = static_cast<bool *> (Mask->data);

            for (int i=0; i<numel; i++) {
                if(brainmaskptr[i]){
                    float tempsum=0;
                    for (int j=0; j<Priors->nt; j++) {
                        int tempind=i+numel*j;
                        if( priorsptr[tempind]<0.0 || priorsptr[tempind]!=priorsptr[tempind] || priorsptr[tempind]>1000 ){
                            priorsptr[tempind]=0.0;
                        }
                        tempsum+=priorsptr[tempind];
                    }
                    if (tempsum>0 && tempsum<1000) {
                        for (int j=0; j<Priors->nt; j++) {
                            int tempind=i+numel*j;
                            priorsptr[tempind]=priorsptr[tempind]/tempsum;
                        }
                        good++;
                    }
                    else{
                        for (int j=0; j<Priors->nt; j++) {
                            int tempind=i+numel*j;
                            priorsptr[tempind]=1.0f/(Priors->nt);
                        }
                        ups++;

                    }
                }
                else{
                    for (int j=0; j<Priors->nt; j++) {
                        int tempind=i+numel*j;
                        priorsptr[tempind]=0;
                    }
                }
            }
        }
        else{

            printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

        }
    }
    else{

        printf("err\tNormalize_NaN_Priors\tWrong mask datatype\n");

    }

    if(verbose>0){
        cout<<"Priors: "<< good<<" good voxels and "<<ups<<" bad voxels" << endl;
        flush(cout);
    }

    return 1;
}


int Normalize_Image_mask(nifti_image * input,
                         nifti_image * Mask,
                         ImageSize * CurrSizes,
                         bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<float>(input);
    }
    if(verbose>0){
        cout<< "Normalizing Input Image" << endl;
    }
    int numel=(int)(rowsize(input)*colsize(input)*depth(input));
    if(Mask->datatype!=DT_BINARY){
        seg_convert2binary(Mask,0.0f);
    }
    if(input->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<segPrecisionTYPE>(input);
    }
    for(long udir=0; udir<CurrSizes->usize;udir++){ // Per Multispectral Image
        bool * brainmaskptr = static_cast<bool *> (Mask->data);
        segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(input->data);
        segPrecisionTYPE * Inputptr=&Inputptrtmp[numel*udir];

        float tempmax=-(1e32);
        float tempmin=1e32;

        for (int i=0; i<numel; i++) {
            if(brainmaskptr[i]){
                if (Inputptr[i]<tempmin) {
                    tempmin=Inputptr[i];
                }
                if (Inputptr[i]>tempmax) {
                    tempmax=Inputptr[i];
                }
            }
        }
        CurrSizes->rescale_max[udir]=tempmax;
        CurrSizes->rescale_min[udir]=tempmin;
        if(verbose>0){
            cout << "Normalization["<<udir<<"] = ["<<tempmin<<","<<tempmax<<"]"<<endl;
        }
        Inputptr=&Inputptrtmp[numel*udir];
        brainmaskptr = static_cast<bool *> (Mask->data);
        for (int i=0; i<numel; i++) {
            //if(brainmaskptr[i]>0){
            //log(number_between_0_and_1 + 1)/log(2)
            Inputptr[i]=logf((((Inputptr[i])-tempmin)/(tempmax-tempmin))+1)/0.693147181;
            if(Inputptr[i]!=Inputptr[i]){
                cout<< "Image has NaNs" << endl;
            }
            /*}
            else{
                Inputptr[i]=0;
            }*/
        }
    }
    return 1;
}


int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          long * shortsize){
    int numel_masked=0;
    int numel = Mask->nvox;
    if(Mask->datatype==DT_BINARY){
        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        for (int i=0; i<numel; i++, Maskptrtmp++) {
            (*Maskptrtmp)>0?numel_masked++:0;
        }
        shortsize[0]=numel_masked;

        int * Short_2_Long_Indices= new int [numel_masked]();
        int * Short_2_Long_Indices_PTR = (int *)(Short_2_Long_Indices);

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++) {
            if ((*Maskptrtmp)>0) {
                Short_2_Long_Indices_PTR[tempindex]=i;
                tempindex++;
            }
            Maskptrtmp++;
        }
        return Short_2_Long_Indices;
    }
    else {
        printf("err\tCreate_Short_2_Long_Matrix\tWrong Mask datatype\n");
        return NULL;

    }
}

int *  Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask){
    int numel = Mask->nvox;
    int * Long_2_Short_Indices= new int [numel]();
    if(Mask->datatype==DT_BINARY){
        bool * Maskptr = static_cast<bool *> (Mask->data);
        bool * Maskptrtmp = Maskptr;
        int * Long_2_Short_Indices_PTR = (int *) Long_2_Short_Indices;

        Maskptrtmp = Maskptr;
        int tempindex=0;
        for (int i=0; i<numel; i++,Maskptrtmp++,Long_2_Short_Indices_PTR++) {
            if ((*Maskptrtmp)>0) {
                (*Long_2_Short_Indices_PTR)=tempindex;
                tempindex++;
            }
            else{
                (*Long_2_Short_Indices_PTR)=-1;
            }
        }
        return Long_2_Short_Indices;
    }
    else {
        cout<< "err\tCreate_Correspondace_Matrices\tWrong Mask datatype\n" << endl;
    }
    return Long_2_Short_Indices;
}

int Normalize_NaN_Priors(nifti_image * Priors,
                         bool verbose)
{
    register int numel = Priors->nx*Priors->ny*Priors->nz;
    register int ups=0;
    register int good=0;
    if(verbose>0){
        cout<< "Normalizing Priors" << endl;
    }
    if(Priors->datatype==NIFTI_TYPE_FLOAT32){
        segPrecisionTYPE * priorsptr = static_cast<segPrecisionTYPE *>(Priors->data);
        for (int i=0; i<numel; i++) {
            float tempsum=0;
            for (int j=0; j<Priors->nt; j++) {
                int tempind=i+numel*j;
                if( priorsptr[tempind]<0.0 || priorsptr[tempind]!=priorsptr[tempind] || priorsptr[tempind]>1000 ){
                    priorsptr[tempind]=0.0;
                }
                tempsum+=priorsptr[tempind];
            }
            if (tempsum>0 && tempsum<1000) {
                for (int j=0; j<Priors->nt; j++) {
                    int tempind=i+numel*j;
                    priorsptr[tempind]=priorsptr[tempind]/tempsum;
                }
                good++;
            }
            else{
                for (int j=0; j<Priors->nt; j++) {
                    int tempind=i+numel*j;
                    priorsptr[tempind]=0.2f;
                }
                ups++;
            }
        }
    }
    else{

        printf("err\tNormalize_NaN_Priors\tWrong Image datatype\n");

    }

    if(verbose>0){
        cout<<"Priors: "<< good<<" good voxels and "<<ups<<" bad voxels" << endl;
        flush(cout);
    }

    return 1;
}


int Normalize_Image(nifti_image * input,
                    ImageSize * CurrSizes,
                    bool verbose)
{
    if(input->datatype!=NIFTI_TYPE_FLOAT32){
        seg_changeDatatype<float>(input);
    }
    if(verbose>0){
        cout<< "Normalizing Input Image" << endl;
    }
    if(input->datatype==NIFTI_TYPE_FLOAT32){
        // if mask is not set up
        for(long udir=0; udir<CurrSizes->usize;udir++){ // Per Multispectral Image
            int numel=(int)(rowsize(input)*colsize(input)*depth(input));
            segPrecisionTYPE * Inputptrtmp = static_cast<segPrecisionTYPE *>(input->data);
            segPrecisionTYPE * Inputptr=&Inputptrtmp[numel*udir];

            float tempmax=0;
            float tempmin=1000000.f;
            for (int i=0; i<numel; i++) {
                if (*Inputptr<tempmin) {
                    tempmin=*Inputptr;
                }
                if (*Inputptr>tempmax) {
                    tempmax=*Inputptr;
                }
                Inputptr++;
            }
            CurrSizes->rescale_max[udir]=tempmax;
            CurrSizes->rescale_min[udir]=tempmin;
            Inputptr=&Inputptrtmp[numel*udir];
            for (int i=0; i<numel; i++) {
                //log(number_between_0_and_1 + 1)/log(2)
                *Inputptr=logf((((*Inputptr)-tempmin)/(tempmax-tempmin))+1)/0.693147181;
                if(*Inputptr!=*Inputptr){
                    cout<< "Image has NaNs" << endl;
                }
                Inputptr++;
            }
        }
    }

    return 1;
}

nifti_image * Copy_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->dim[5]=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++){Resultdata[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(long currclass=0; currclass<CurrSizes->numclass;currclass++){

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++){
            Resultdata_class[Short_2_Long_Indices[i]]=*Expec_PTR;
            if( i == 1444500){
                //cout << "In Copy_Expec_to_Result_mask : " << (*Expec_PTR) << endl; 
            }
        }
    }
    return Result;
}
//saurabh added
nifti_image * Copy_Correct_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                                segPrecisionTYPE * Outlierness,
                                                int * Short_2_Long_Indices,
                                                nifti_image * T1,
                                                char * filename,
                                                ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->dim[5]=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++){Resultdata[i]=0;}

    int * Short_2_Long_Indices_PRT = (int *) Short_2_Long_Indices;

    int class_nvox=Result->nx*Result->ny*Result->nz;

    Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;
    for(long currclass=0; currclass<CurrSizes->numclass;currclass++){

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numelmasked];
        segPrecisionTYPE * Outlier_PTR = &Outlierness[(currclass)*CurrSizes->numelmasked];
        Short_2_Long_Indices_PRT= (int *) Short_2_Long_Indices;

        for(long i=0; i<(long)CurrSizes->numelmasked; i++,Short_2_Long_Indices_PRT++,Expec_PTR++,Outlier_PTR++){
            Resultdata_class[Short_2_Long_Indices[i]]=(*Expec_PTR) * (*Outlier_PTR);
            if( i == 1444500){
                //cout << "In Copy_Correct_Expec_to_Result_mask : " <<(*Expec_PTR) << ' ' << (*Outlier_PTR) << ' '<< (*Expec_PTR) * (*Outlier_PTR) <<endl;  
            }
        }
    }
    //cout << "I am in Copy_Correct_Expec_to_Result_mask function... !!!!" << endl;
    return Result;
}


nifti_image * Copy_Expec_to_Result(segPrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes){

    nifti_image * Result = nifti_copy_nim_info(T1);
    Result->dim[0]=4;
    Result->dim[4]=CurrSizes->numclass;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    Result->scl_inter=0;
    Result->scl_slope=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++){Resultdata[i]=0;}
    int class_nvox=Result->nx*Result->ny*Result->nz;

    for(long currclass=0; currclass<CurrSizes->numclass;currclass++){

        segPrecisionTYPE * Resultdata_class = &Resultdata[(currclass)*class_nvox];
        segPrecisionTYPE * Expec_PTR = &Expec[(currclass)*CurrSizes->numel];

        for(long i=0; i<(long)CurrSizes->numel; i++,Expec_PTR++){
            Resultdata_class[i]=*Expec_PTR;
        }
    }
    return Result;
}


void BiasCorrection(segPrecisionTYPE * BiasField,
                    segPrecisionTYPE * BiasFieldCoefs,
                    nifti_image * T1,
                    segPrecisionTYPE * Expec,
                    segPrecisionTYPE * Outlierness,
                    segPrecisionTYPE * M,
                    segPrecisionTYPE * V,
                    int biasOrder,
                    ImageSize * CurrSizes,
                    bool flag_Bias,
                    int verbose_level){

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};

//#ifdef _OPENMP
//        segPrecisionTYPE * Athread =new segPrecisionTYPE [omp_get_max_threads()*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]();
//#endif
// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0){
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0){
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }

// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;

        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    currindex=iz*plane_size+iy*col_size+ix;
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++){
                        Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }

// Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)Dims3d[2]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[2]);



        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        segPrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        segPrecisionTYPE current_Tempvar=0.0f;

// Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

                    Basisptr1= (segPrecisionTYPE *) Basis;
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    z_bias_index_shift=2;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++){
                        // Because Powerorder is alwais int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        // Instead, although slower, one can use
                        // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
                    }
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    //#ifdef _OPENMP
                    //                    Aptr=&Athread[omp_get_thread_num()*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)];
                    //                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                    //                        Basisptr2= &Basis[j2];
                    //                        segPrecisionTYPE * Aptr2= &Aptr[j2+j2*UsedBasisFunctions];
                    //                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr2++, Basisptr2++){
                    //                            (*Aptr2)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                    //                        }
                    //                    }
                    //#else
                    Aptr= (segPrecisionTYPE *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }


        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
        }

// CALC MATRIX B

//Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;

                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++){
                        segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[linearindexes+TotalLength*j];
                        Wij=tmpexpec*(invV[j]);
                        Wi+=Wij;
                        Yest+=Wij*(currM[j]);
                        Ysum+=Wij;
                    }
                    Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                    tempvarindex++;

                }
            }
        }

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                        B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                               pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                               pow_int((((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
                               Tempvar[tempvarindex];
                        if(B[i2]!=B[i2]){
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }



        for (int iz=0; iz<maxiz; iz++) {
            for (int iy=0; iy<maxiy; iy++) {
                for (int ix=0; ix<maxix; ix++) {
                    segPrecisionTYPE tmpbiasfield=0.0f;
                    segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
                    segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
                    segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    tmpbiasfield=0.0f;
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            for(int yorder=0; yorder<=(order-xorder); yorder++){
                                int zorder=order-yorder-xorder;
                                tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField[linearindexes+multispec*CurrSizes->numel]=tmpbiasfield;
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }

        delete [] Basis;
        delete [] Tempvar;
    }


}


void BiasCorrection_SPARCS(float * BiasField,
                           float * T1,
                           float * Expec,
                           float * Mask,
                           float * M,
                           float * V,
                           int biasOrder,
                           int nrOfClasses,
                           int aceletation_factor,
                           int xyzsize[3]) {


    if(aceletation_factor<1){
        cout<<"ERROR: The acceleration factor has to be above or equal to 1"<<endl;
        flush(cout);
        return;
    }
    //int aceletation_factor=redux_factor_for_bias;
    int TotalLength = xyzsize[0]*xyzsize[1]*xyzsize[2];
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }

    float invV[maxNumbClass];
    float currM[maxNumbClass];

    float * sampledData = &T1[0];
    // Precompute the M and V  inverses

    for(int i=0; i<nrOfClasses; i++){
        invV[i]=1.0f/V[i];

        currM[i]=M[i];
    }

    float A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
    float B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
    float C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};

    // Precompute sizes
    int maxix = (int)(xyzsize[0]);
    int maxiy = (int)(xyzsize[1]);
    int maxiz = (int)(xyzsize[2]);
    int Dims3d[3]={0};
    Dims3d[0]=maxix;
    Dims3d[1]=maxiy;
    Dims3d[2]=maxiz;

    // Precompute number of samples as it was never computed
    int samplecount=0;
    int linearindexes=0;

    for (int iz=0; iz<maxiz; iz+=aceletation_factor) {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor) {
            for (int ix=0; ix<maxix; ix+=aceletation_factor) {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0){
                    samplecount++;
                }
            }
        }
    }


    // CALC MATRIX A

    // Calc W (Van Leemput 1999 eq 7)

    float * Tempvar= new float [samplecount] ();
    float Tempvar_tmp=0;
    int tempvarindex=0;
    for (int iz=0; iz<maxiz; iz+=aceletation_factor) {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor) {
            for (int ix=0; ix<maxix; ix+=aceletation_factor) {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0){
                    Tempvar_tmp=0;
                    for(int j=0; j<nrOfClasses; j++){
                        Tempvar_tmp+=Expec[linearindexes+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }
    }

    // Precompute shifts
    float not_point_five_times_dims_x=(0.5f*(float)Dims3d[0]);
    float not_point_five_times_dims_y=(0.5f*(float)Dims3d[1]);
    float not_point_five_times_dims_z=(0.5f*(float)Dims3d[2]);
    float inv_not_point_five_times_dims_x=1.0f/(0.5f*(float)Dims3d[0]);
    float inv_not_point_five_times_dims_y=1.0f/(0.5f*(float)Dims3d[1]);
    float inv_not_point_five_times_dims_z=1.0f/(0.5f*(float)Dims3d[2]);
    float * Basis= new float[UsedBasisFunctions]();
    float xpos=0.0f;
    float ypos=0.0f;
    float zpos=0.0f;
    int x_bias_index_shift=0;
    int y_bias_index_shift=1;
    int z_bias_index_shift=2;
    float current_Tempvar=0.0f;

    // Calc A'WA (Van Leemput 1999 eq 7)
    tempvarindex=0;
    float * Basisptr1= (float *) Basis;
    float * Basisptr2= (float *) Basis;
    float * Aptr= (float *) A;
    for (int iz=0; iz<maxiz; iz+=aceletation_factor) {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor) {
            for (int ix=0; ix<maxix; ix+=aceletation_factor) {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0){
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    z_bias_index_shift=2;
                    Basisptr1= (float *) Basis;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++){
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                    }
                    Basisptr1= (float *) Basis;
                    Aptr= (float *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }
    }

    seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

    for(int j2=0; j2<UsedBasisFunctions; j2++){
        for(int i2=j2; i2<UsedBasisFunctions; i2++){
            RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
            RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
        }
    }

    seg_Matrix <double> RealA_inv(UsedBasisFunctions);
    RealA_inv.copymatrix(RealA);
    RealA_inv.invert();


    // CALC MATRIX B
    //Precompute WR (Van Leemput 1999 eq 7)
    float Wi;
    float Wij;
    float Yest;
    float Ysum;
    tempvarindex=0;

    for (int iz=0; iz<maxiz; iz+=aceletation_factor) {
        for (int iy=0; iy<maxiy; iy+=aceletation_factor) {
            for (int ix=0; ix<maxix; ix+=aceletation_factor) {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0){
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(int j=0; j<nrOfClasses; j++){
                        float tmpexpec = (float)Expec[linearindexes+TotalLength*j];
                        Wij=tmpexpec*(invV[j]);
                        Wi+=Wij;
                        Yest+=Wij*(currM[j]);
                        Ysum+=Wij;
                    }
                    Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                    tempvarindex++;
                }
            }
        }
    }


    for(int i2=0; i2<UsedBasisFunctions; i2++){
        tempvarindex=0;
        B[i2]=0;
        for (int iz=0; iz<maxiz; iz+=aceletation_factor) {
            for (int iy=0; iy<maxiy; iy+=aceletation_factor) {
                for (int ix=0; ix<maxix; ix+=aceletation_factor) {
                    linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                    if((Mask[linearindexes])>0){
                        linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                        B[i2]+=pow_int((((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                               pow_int((((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                               pow_int((((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+i2*3])*
                               Tempvar[tempvarindex];
                        if(B[i2]!=B[i2]){
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }
    }

    seg_Matrix <double> RealB(UsedBasisFunctions,1);

    for(int i2=0; i2<UsedBasisFunctions; i2++){
        RealB.setvalue(i2,0,(double)(B[i2]));
    }

    seg_Matrix <double> RealC(UsedBasisFunctions,1);

    RealC.settoproduct(RealA_inv,RealB);

    double cvalue=0.0f;
    bool success;
    for(int i2=0; i2<UsedBasisFunctions; i2++){
        RealC.getvalue(i2,0,cvalue,success);
        C[i2]=(float)(cvalue);
    }

    for(int j2=0; j2<UsedBasisFunctions; j2++){
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            double cvalue=0.0f;
            bool success;
            RealB.getvalue(i2,j2,cvalue,success);

        }
    }

    for (int iz=0; iz<maxiz; iz++) {
        for (int iy=0; iy<maxiy; iy++) {
            for (int ix=0; ix<maxix; ix++) {
                linearindexes=(iz)*(xyzsize[0])*(xyzsize[1])+(iy)*(xyzsize[0])+ix;
                if((Mask[linearindexes])>0){
                    float tmpbiasfield=0.0f;
                    float currxpower[maxAllowedBCPowerOrder];
                    float currypower[maxAllowedBCPowerOrder];
                    float currzpower[maxAllowedBCPowerOrder];
                    tmpbiasfield=0.0f;
                    xpos=(((float)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((float)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    zpos=(((float)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                    get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            for(int yorder=0; yorder<=(order-xorder); yorder++){
                                int zorder=order-yorder-xorder;
                                tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                ind++;
                            }
                        }
                    }
                    BiasField[linearindexes]=tmpbiasfield;

                }
                else{
                    BiasField[linearindexes]=1;
                }
            }
        }
    }

    delete [] Basis;
    delete [] Tempvar;
}


void BiasCorrection_mask(segPrecisionTYPE * BiasField,
                         segPrecisionTYPE * BiasFieldCoefs,
                         nifti_image * T1,
                         int * Long_2_Short_Indices,
                         segPrecisionTYPE * Expec,
                         segPrecisionTYPE * Outlierness,
                         segPrecisionTYPE * M,
                         segPrecisionTYPE * V,
                         int biasOrder,
                         ImageSize * CurrSizes,
                         bool flag_Bias,
                         int verbose_level) {

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    int TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2 *(biasOrder+3)/3);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3))]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            for(int yorder=0; yorder<=(order-xorder); yorder++){
                int zorder=order-yorder-xorder;
                PowerOrder[ind] =xorder;
                PowerOrder[ind+1] =yorder;
                PowerOrder[ind+2] =zorder;
                ind += 3;
            }
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = static_cast<segPrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2*(maxAllowedBCPowerOrder+3)/3)]={0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int plane_size = (int)(CurrSizes->xsize)*(CurrSizes->ysize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int maxiz = (int)(CurrSizes->zsize);
        int Dims3d[3]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;
        Dims3d[2]=maxiz;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0){
            for (int iz=0; iz<maxiz; iz+=reduxfactor) {
                for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                    for (int ix=0; ix<maxix; ix+=reduxfactor) {
                        linearindexes=iz*plane_size+iy*col_size+ix;
                        if(Long_2_Short_Indices[linearindexes]>=0){
                            samplecount++;
                        }
                    }
                }
            }
            if(verbose_level>0){
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)
//cout << "Calculating C = inv(A'WA) WR";
//flush(cout);

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    currshortindex=Long_2_Short_Indices[iz*plane_size+iy*col_size+ix];
                    if(currshortindex>=0){
                        Tempvar_tmp=0;
                        if(Outlierness==NULL){
                            for(int j=0; j<nrOfClasses; j++){
                                Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                            }
                        }
                        else{
                            for(int j=0; j<nrOfClasses; j++){
                                Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*Outlierness[currshortindex+TotalLength*j]*invV[j];
                            }

                        }
                        Tempvar[tempvarindex]=Tempvar_tmp;
                        tempvarindex++;
                    }
                }
            }
        }

// Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE not_point_five_times_dims_z=(0.5f*(segPrecisionTYPE)Dims3d[2]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);
        segPrecisionTYPE inv_not_point_five_times_dims_z=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[2]);

        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        segPrecisionTYPE zpos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        int z_bias_index_shift=2;
        segPrecisionTYPE current_Tempvar=0.0f;

// Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        Basisptr1= (segPrecisionTYPE *) Basis;
                        current_Tempvar=Tempvar[tempvarindex];
                        xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        x_bias_index_shift=0;
                        y_bias_index_shift=1;
                        z_bias_index_shift=2;
                        for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=3,y_bias_index_shift+=3,z_bias_index_shift+=3, Basisptr1++){
                            // Because Powerorder is always int, use a special power function (faster)
                            *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift])*pow_int(zpos,PowerOrder[z_bias_index_shift]));
                        }


                        Basisptr1= (segPrecisionTYPE *) Basis;
                        Aptr= (segPrecisionTYPE *) A;
                        for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                            Basisptr2= &Basis[j2];
                            Aptr= &A[j2+j2*UsedBasisFunctions];
                            for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                                (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                            }
                        }
                        tempvarindex++;
                    }
                }
            }
        }



        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



// CALC MATRIX B

//Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iz=0; iz<maxiz; iz+=reduxfactor) {
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        Wi=0;
                        Wij=0;
                        Yest=0;
                        Ysum=0;
                        for(int j=0; j<nrOfClasses; j++){
                            segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[currshortindex+TotalLength*j];
                            Wij=tmpexpec*(invV[j]);
                            Wi+=Wij;
                            Yest+=Wij*(currM[j]);
                            Ysum+=Wij;
                        }
                        Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                        tempvarindex++;
                    }
                }
            }
        }

//#ifdef _OPENMP
//#pragma omp parallel shared(Tempvar,PowerOrder,CurrSizes,B,Long_2_Short_Indices) private(UsedBasisFunctions,maxiz,maxiy,maxix, not_point_five_times_dims_x, not_point_five_times_dims_y, not_point_five_times_dims_z,reduxfactor,inv_not_point_five_times_dims_x, inv_not_point_five_times_dims_y, inv_not_point_five_times_dims_z)
//#endif
        for(int bfindex=0; bfindex<UsedBasisFunctions; bfindex++){
            int tempvarindex2=0;
            int linearindexes2=0;
            segPrecisionTYPE * BPTR=&B[bfindex];
            *BPTR=0;
            for (int iz2=0; iz2<maxiz; iz2+=reduxfactor) {
                for (int iy2=0; iy2<maxiy; iy2+=reduxfactor) {
                    for (int ix2=0; ix2<maxix; ix2+=reduxfactor) {
                        linearindexes2=(iz2)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy2)*(CurrSizes->xsize)+ix2;
                        if(Long_2_Short_Indices[linearindexes2]>=0){
                            *BPTR+=pow_int((((segPrecisionTYPE)ix2-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+bfindex*3])*
                                   pow_int((((segPrecisionTYPE)iy2-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+bfindex*3])*
                                   pow_int((((segPrecisionTYPE)iz2-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z),PowerOrder[2+bfindex*3])*
                                   Tempvar[tempvarindex2];


                            tempvarindex2++;
                        }
                    }
                }
            }
            if(*BPTR!=*BPTR){
                *BPTR=1;
            }
        }

//#ifdef _OPENMP
//#pragma omp barrier
//#endif
        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }

//#ifdef _OPENMP
//#pragma omp parallel
//#endif
        for (int iz=0; iz<maxiz; iz++) {
            segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
            segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
            segPrecisionTYPE currzpower[maxAllowedBCPowerOrder];
            segPrecisionTYPE tmpbiasfield=0.0f;
            for (int iy=0; iy<maxiy; iy++) {
                for (int ix=0; ix<maxix; ix++) {
                    linearindexes=(iz)*(CurrSizes->xsize)*(CurrSizes->ysize)+(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        tmpbiasfield=0.0f;
                        xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                        ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                        zpos=(((segPrecisionTYPE)iz-not_point_five_times_dims_z)*inv_not_point_five_times_dims_z);
                        get_xyz_pow_int(xpos, ypos, zpos, currxpower, currypower, currzpower, biasOrder);
                        int ind2=0;
                        for(int order=0; order<=biasOrder; order++){
                            for(int xorder=0; xorder<=order; xorder++){
                                for(int yorder=0; yorder<=(order-xorder); yorder++){
                                    int zorder=order-yorder-xorder;
                                    tmpbiasfield-=C[ind2]*currxpower[xorder]*currypower[yorder]*currzpower[zorder];
                                    ind2++;
                                }
                            }
                        }
                        BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                    }
                }
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
}




void BiasCorrection2D(segPrecisionTYPE * BiasField,
                      segPrecisionTYPE * BiasFieldCoefs,
                      nifti_image * T1,
                      segPrecisionTYPE * Expec,
                      segPrecisionTYPE * Outlierness,
                      segPrecisionTYPE * M,
                      segPrecisionTYPE * V,
                      int biasOrder,
                      ImageSize * CurrSizes,
                      bool flag_Bias,
                      int verbose_level) {

    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    int nrOfClasses = CurrSizes->numclass;

    int TotalLength = CurrSizes->numel;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            int yorder=order-xorder;
            PowerOrder[ind] =xorder;
            PowerOrder[ind+1] =yorder;
            ind += 3;
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses

        for(int i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }



        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currindex=0;
        if(CurrSizes->numelbias==0){
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    samplecount++;
                }
            }
            if(verbose_level>0){
                cout << "Samplecount = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)

        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                currindex=iy*col_size+ix;
                Tempvar_tmp=0;
                for(int j=0; j<nrOfClasses; j++){
                    Tempvar_tmp+=Expec[currindex+TotalLength*j]*invV[j];
                }
                Tempvar[tempvarindex]=Tempvar_tmp;
                tempvarindex++;
            }
        }


// Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);


        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        segPrecisionTYPE current_Tempvar=0.0f;

// Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                Basisptr1= (segPrecisionTYPE *) Basis;
                current_Tempvar=Tempvar[tempvarindex];
                xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                x_bias_index_shift=0;
                y_bias_index_shift=1;
                for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++){
                    // Because Powerorder is alwais int, use a special power function (faster)
                    *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    // Instead, although slower, one can use
                    // TmpA=(pow(xpos,PowerOrder[0+j2*3])*pow(ypos,PowerOrder[1+j2*3])*pow(zpos,PowerOrder[2+j2*3]));
                }
                Basisptr1= (segPrecisionTYPE *) Basis;
                Aptr= (segPrecisionTYPE *) A;
                for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                    Basisptr2= &Basis[j2];
                    Aptr= &A[j2+j2*UsedBasisFunctions];
                    for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                        (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                    }
                }
                tempvarindex++;
            }
        }




        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



// CALC MATRIX B

//Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;

                Wi=0;
                Wij=0;
                Yest=0;
                Ysum=0;
                for(int j=0; j<nrOfClasses; j++){
                    segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[linearindexes+TotalLength*j];
                    Wij=tmpexpec*(invV[j]);
                    Wi+=Wij;
                    Yest+=Wij*(currM[j]);
                    Ysum+=Wij;
                }
                Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                tempvarindex++;

            }
        }

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iy)*(CurrSizes->xsize)+ix;
                    B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*3])*
                           pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*3])*
                           Tempvar[tempvarindex];
                    if(B[i2]!=B[i2]){
                        B[i2]=1;
                    }
                    tempvarindex++;
                }
            }
        }

        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
        segPrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++) {
            for (int ix=0; ix<maxix; ix++) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                tmpbiasfield=0.0f;
                xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                get_xy_pow_int(xpos, ypos, currxpower, currypower, biasOrder);
                ind=0;
                for(int order=0; order<=biasOrder; order++){
                    for(int xorder=0; xorder<=order; xorder++){
                        int yorder=order-xorder;
                        tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                        ind++;
                    }
                }
                BiasField[linearindexes+multispec*CurrSizes->numel]=tmpbiasfield;
            }
        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }

        delete [] Basis;
        delete [] Tempvar;
    }


}





void BiasCorrection_mask2D(segPrecisionTYPE * BiasField,
                           segPrecisionTYPE * BiasFieldCoefs,
                           nifti_image * T1,
                           int * Long_2_Short_Indices,
                           segPrecisionTYPE * Expec,
                           segPrecisionTYPE * Outlierness,
                           segPrecisionTYPE * M,
                           segPrecisionTYPE * V,
                           int biasOrder,
                           ImageSize * CurrSizes,
                           bool flag_Bias,
                           int verbose_level)
{
    if(verbose_level>0){
        cout << "Optimising the Bias Field with order " << biasOrder<< endl;
        if(CurrSizes->usize>1){
            cout<< "Assuming fully decoupled bias-fields" << endl;
        }
        flush(cout);
    }
    int reduxfactor=reduxFactorForBias;
    long nrOfClasses = CurrSizes->numclass;
    //nrOfClasses = 1;
    //int nrOfClasses = non_PV_numclass;
    long TotalLength = CurrSizes->numelmasked;
    int UsedBasisFunctions=(int)((biasOrder+1) * (biasOrder+2)/2);
    segPrecisionTYPE * sampledData = static_cast<segPrecisionTYPE *>(T1->data);


    // Precompute Powers depending on the current BiasOrder
    int PowerOrder [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0};
    int ind=0;
    for(int order=0; order<=biasOrder; order++){
        for(int xorder=0; xorder<=order; xorder++){
            int yorder=order-xorder;
            PowerOrder[ind] =xorder;
            PowerOrder[ind+1] =yorder;
            ind += 2;
        }
    }
    segPrecisionTYPE invV[maxNumbClass];
    segPrecisionTYPE currM[maxNumbClass];

    for(long multispec=0; multispec<CurrSizes->usize; multispec++){
        sampledData = static_cast<segPrecisionTYPE *>(T1->data);
        sampledData = &sampledData[multispec*CurrSizes->numel];
        // Precompute the M and V  inverses
        for(long i=0; i<nrOfClasses; i++){
            invV[i]=1.0f/(V[i*CurrSizes->usize*CurrSizes->usize+multispec+multispec*CurrSizes->usize]);
            currM[i]=M[i*CurrSizes->usize+multispec];
        }
        segPrecisionTYPE A [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2*((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2)/2)]={0.0f};
        segPrecisionTYPE B [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]={0.0f};
        segPrecisionTYPE C [((maxAllowedBCPowerOrder+1)*(maxAllowedBCPowerOrder+2))/2]={0.0f};


// Precompute sizes
        int col_size = (int)(CurrSizes->xsize);
        int maxix = (int)(CurrSizes->xsize);
        int maxiy = (int)(CurrSizes->ysize);
        int Dims3d[2]={0};
        Dims3d[0]=maxix;
        Dims3d[1]=maxiy;

// Precompute number of samples as it was never computed
        int samplecount=0;
        int linearindexes=0;
        int currshortindex=0;
        if(CurrSizes->numelbias==0){
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {
                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=iy*col_size+ix;
                    if(Long_2_Short_Indices[linearindexes]>=0){
                        samplecount++;
                    }
                }
            }
            if(verbose_level>0){
                cout << "Number of samples for BiasField = " << samplecount<<"\n";
                flush(cout);
            }
            CurrSizes->numelbias=samplecount;
        }
        else{
            samplecount=CurrSizes->numelbias;
        }



// CALC MATRIX A

// Calc W (Van Leemput 1999 eq 7)
//cout << "Calculating C = inv(A'WA) WR";
//flush(cout);
        segPrecisionTYPE * Tempvar= new segPrecisionTYPE [samplecount] ();
        segPrecisionTYPE Tempvar_tmp=0;
        currshortindex=0;
        int tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                currshortindex=Long_2_Short_Indices[iy*col_size+ix];
                if(currshortindex>=0){
                    Tempvar_tmp=0;
                    for(long j=0; j<nrOfClasses; j++){
                        Tempvar_tmp+=Expec[currshortindex+TotalLength*j]*invV[j];
                    }
                    Tempvar[tempvarindex]=Tempvar_tmp;
                    tempvarindex++;
                }
            }
        }
// Precompute shifts
        segPrecisionTYPE not_point_five_times_dims_x=(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE not_point_five_times_dims_y=(0.5f*(segPrecisionTYPE)Dims3d[1]);

        segPrecisionTYPE inv_not_point_five_times_dims_x=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[0]);
        segPrecisionTYPE inv_not_point_five_times_dims_y=1.0f/(0.5f*(segPrecisionTYPE)Dims3d[1]);


        segPrecisionTYPE * Basis= new segPrecisionTYPE[UsedBasisFunctions]();
        segPrecisionTYPE xpos=0.0f;
        segPrecisionTYPE ypos=0.0f;
        int x_bias_index_shift=0;
        int y_bias_index_shift=1;
        segPrecisionTYPE current_Tempvar=0.0f;
// Calc A'WA (Van Leemput 1999 eq 7)
        tempvarindex=0;
        segPrecisionTYPE * Basisptr1= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Basisptr2= (segPrecisionTYPE *) Basis;
        segPrecisionTYPE * Aptr= (segPrecisionTYPE *) A;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    current_Tempvar=Tempvar[tempvarindex];
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    x_bias_index_shift=0;
                    y_bias_index_shift=1;
                    for(int j2=0; j2<UsedBasisFunctions; j2++,x_bias_index_shift+=2,y_bias_index_shift+=2, Basisptr1++){
                        // Because Powerorder is always int, use a special power function (faster)
                        *Basisptr1=(pow_int(xpos,PowerOrder[x_bias_index_shift])*pow_int(ypos,PowerOrder[y_bias_index_shift]));
                    }
                    Basisptr1= (segPrecisionTYPE *) Basis;
                    Aptr= (segPrecisionTYPE *) A;
                    for(int j2=0; j2<UsedBasisFunctions; j2++, Basisptr1++){
                        Basisptr2= &Basis[j2];
                        Aptr= &A[j2+j2*UsedBasisFunctions];
                        for(int i2=j2; i2<UsedBasisFunctions; i2++, Aptr++, Basisptr2++){
                            (*Aptr)+=(*Basisptr2)*(current_Tempvar)*(*Basisptr1);
                        }
                    }
                    tempvarindex++;
                }
            }
        }

        seg_Matrix <double> RealA(UsedBasisFunctions,UsedBasisFunctions);

        for(int j2=0; j2<UsedBasisFunctions; j2++){
            for(int i2=j2; i2<UsedBasisFunctions; i2++){
                RealA.setvalue(i2,j2,(double)(A[i2+j2*UsedBasisFunctions]));
                RealA.setvalue(j2,i2,(double)(A[i2+j2*UsedBasisFunctions]));
            }
        }

        seg_Matrix <double> RealA_inv(UsedBasisFunctions);
        RealA_inv.copymatrix(RealA);
        RealA_inv.invert();

        if(verbose_level>1){
            seg_Matrix <double> RealA_test(UsedBasisFunctions);
            RealA_test.settoproduct(RealA,RealA_inv);
            RealA_test.comparetoidentity();
            //RealA.dumpmatrix();
        }



// CALC MATRIX B

//Precompute WR (Van Leemput 1999 eq 7)
        segPrecisionTYPE Wi;
        segPrecisionTYPE Wij;
        segPrecisionTYPE Yest;
        segPrecisionTYPE Ysum;
        tempvarindex=0;
        for (int iy=0; iy<maxiy; iy+=reduxfactor) {
            for (int ix=0; ix<maxix; ix+=reduxfactor) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    Wi=0;
                    Wij=0;
                    Yest=0;
                    Ysum=0;
                    for(long j=0; j<nrOfClasses; j++){
                        segPrecisionTYPE tmpexpec = (segPrecisionTYPE)Expec[currshortindex+TotalLength*j];
                        Wij=tmpexpec*(invV[j]);
                        Wi+=Wij;
                        Yest+=Wij*(currM[j]);
                        Ysum+=Wij;
                    }
                    Tempvar[tempvarindex]=Wi*(sampledData[linearindexes]-(Yest/Ysum));
                    tempvarindex++;
                }
            }
        }

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            tempvarindex=0;
            B[i2]=0;
            for (int iy=0; iy<maxiy; iy+=reduxfactor) {

                for (int ix=0; ix<maxix; ix+=reduxfactor) {
                    linearindexes=(iy)*(CurrSizes->xsize)+ix;
                    currshortindex=Long_2_Short_Indices[linearindexes];
                    if(currshortindex>=0){
                        B[i2]+=pow_int((((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x),PowerOrder[0+i2*2])*
                               pow_int((((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y),PowerOrder[1+i2*2])*
                               Tempvar[tempvarindex];

                        if(B[i2]!=B[i2]){
                            B[i2]=1;
                        }
                        tempvarindex++;
                    }
                }
            }
        }

        seg_Matrix <double> RealB(UsedBasisFunctions,1);

        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealB.setvalue(i2,0,(double)(B[i2]));
        }

        seg_Matrix <double> RealC(UsedBasisFunctions,1);

        RealC.settoproduct(RealA_inv,RealB);
        if(verbose_level>1){
            cout << "C= " << endl;
            RealC.dumpmatrix();
        }

        double cvalue=0.0f;
        bool success;
        for(int i2=0; i2<UsedBasisFunctions; i2++){
            RealC.getvalue(i2,0,cvalue,success);
            C[i2]=(segPrecisionTYPE)(cvalue);
        }

        segPrecisionTYPE currxpower[maxAllowedBCPowerOrder];
        segPrecisionTYPE currypower[maxAllowedBCPowerOrder];
        segPrecisionTYPE tmpbiasfield=0.0f;
        for (int iy=0; iy<maxiy; iy++) {
            for (int ix=0; ix<maxix; ix++) {
                linearindexes=(iy)*(CurrSizes->xsize)+ix;
                currshortindex=Long_2_Short_Indices[linearindexes];
                if(currshortindex>=0){
                    tmpbiasfield=0.0f;
                    xpos=(((segPrecisionTYPE)ix-not_point_five_times_dims_x)*inv_not_point_five_times_dims_x);
                    ypos=(((segPrecisionTYPE)iy-not_point_five_times_dims_y)*inv_not_point_five_times_dims_y);
                    get_xy_pow_int(xpos, ypos, currxpower, currypower, biasOrder);
                    ind=0;
                    for(int order=0; order<=biasOrder; order++){
                        for(int xorder=0; xorder<=order; xorder++){
                            int yorder=order-xorder;
                            tmpbiasfield-=C[ind]*currxpower[xorder]*currypower[yorder];
                            ind++;
                        }
                    }
                    BiasField[currshortindex+multispec*CurrSizes->numelmasked]=tmpbiasfield;
                }
            }

        }

        for( int i=0; i<UsedBasisFunctions; i++){
            BiasFieldCoefs[i+multispec*UsedBasisFunctions]=C[i];
        }
        delete [] Basis;
        delete [] Tempvar;
    }
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

int calcE_mask(nifti_image * T1,
               segPrecisionTYPE * IterPrior,
               segPrecisionTYPE * Expec,
               double * loglik,
               segPrecisionTYPE * BiasField,
               segPrecisionTYPE * Outlierness,
               segPrecisionTYPE OutliernessThreshold,
               int * S2L,
               segPrecisionTYPE * M,
               segPrecisionTYPE * V,
               ImageSize * CurrSizes,
               int verbose)
{
    int numel_masked=CurrSizes->numelmasked;
    int num_class=CurrSizes->numclass;
    bool OutliernessFlag=(Outlierness==NULL)?0:1;
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]={0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]={0.0f};

    int Expec_offset [maxNumbClass]={0};

    for (int cl=0; cl<num_class; cl++) {
        Expec_offset[cl]=(int) cl*numel_masked;
        if(CurrSizes->usize>1){
            seg_Matrix <double> Vmat(CurrSizes->usize,CurrSizes->usize);

            for(long j2=0; j2<CurrSizes->usize; j2++){
                for(long i2=j2; i2<CurrSizes->usize; i2++){
                    Vmat.setvalue(i2,j2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                    Vmat.setvalue(j2,i2,(double)(V[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]));
                }
            }
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if(verbose>1){
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            Vmat.invert();
            double cvalue=0.0f;
            bool success;
            if(verbose>1){
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for(long j2=0; j2<CurrSizes->usize; j2++){
                if(verbose>1){
                    if(j2!=0){
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<CurrSizes->usize; i2++){
                    Vmat.getvalue(i2,j2,cvalue,success);
                    inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]=(segPrecisionTYPE)(cvalue);
                    if(verbose>1){
                        cout<<inv_v[i2+j2*CurrSizes->usize+cl*CurrSizes->usize*CurrSizes->usize]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(verbose>1){
                cout<< endl;
            }
        }
        else{
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*V[cl]));
            inv_v[cl]=1/V[cl];
        }
    }
    loglik[0]=0;

    //int * Expec_offset_PTR= (int *) Expec_offset;

    float logliktmp=0.0f;


#ifdef _OPENMP
    float * loglikthread = new float [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(Expec,loglikthread,T1,BiasField,Outlierness,IterPrior)
#endif
    for (int i=0; i<numel_masked;i++) {
        segPrecisionTYPE * T1_PTR = static_cast<segPrecisionTYPE *>(T1->data);
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize];
        segPrecisionTYPE SumExpec=0.0f;

        for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++)
            T1_Bias_corr[Multispec]=(BiasField!=NULL)?(T1_PTR[S2L[i]+Multispec*CurrSizes->numel] + BiasField[i+Multispec*numel_masked]):(T1_PTR[S2L[i]+Multispec*CurrSizes->numel]);



        //Expec_offset_PTR=Expec_offset;

        for (int cl=0; cl<num_class; cl++) {
            segPrecisionTYPE mahal=0.0f;
            for(long Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                segPrecisionTYPE tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - M[cl*(CurrSizes->usize)+Multispec]);
                for(long Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(CurrSizes->usize)+Multispec2])*inv_v[cl*CurrSizes->usize*CurrSizes->usize+Multispec+Multispec2*CurrSizes->usize]*tmpT1_BC_minusM;
                }
            }

            if(OutliernessFlag){
                float outvalue=(expf(mahal)+0.01)/(expf(mahal)+expf(-0.5*(OutliernessThreshold*OutliernessThreshold))+0.01);
                Outlierness[i+Expec_offset[cl]]=outvalue;
            }
            Expec[i+Expec_offset[cl]]=IterPrior[i+Expec_offset[cl]] * expf(mahal) * inv_sqrt_V_2pi[cl];
            SumExpec+=Expec[i+Expec_offset[cl]];
        }

        if (SumExpec<=0.0 || SumExpec!=SumExpec){
            for (int cl=0; cl<num_class; cl++) {
                Expec[i+Expec_offset[cl]]=(float)(1)/(float)(num_class);
            }

        }
        else{

            for (int cl=0; cl<num_class; cl++) {
                Expec[i+Expec_offset[cl]]=Expec[i+Expec_offset[cl]]/SumExpec;
            }
#ifdef _OPENMP
            loglikthread[omp_get_thread_num()]+=logf(SumExpec);
#else
            logliktmp+=logf(SumExpec);
#endif

        }
    }

#ifdef _OPENMP
    for(long i =0; i<(long)omp_get_max_threads(); i++)
        logliktmp+=loglikthread[i];
#endif

    loglik[0]=logliktmp;
    return 1;
}
