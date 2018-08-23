#ifndef _SEG_EM_CPP
#define _SEG_EM_CPP
#include "_seg_EM_old.h"
#include "_seg_matrix.h"

seg_EM_old::seg_EM_old(int _numb_classes, int _nu,int _nt)
{

    this->InputImage=NULL; // pointer to external

    this->FilenameOut="Segmentation.nii.gz";

    this->dimentions=1;
    this->nx=1;
    this->ny=0;
    this->nz=0;
    this->nu=_nu*_nt;
    this->nt=1;
    this->dx=0;
    this->dy=0;
    this->dz=0;
    this->numElements=0;
    this->aprox=true;
    this->iter=0;
    this->checkpoint_iter=0;
    this->ratio=1000;

    this->numberOfClasses=_numb_classes;
    this->M= new float [maxMultispectalSize*maxNumbClass];
    for(int i=0; i<(maxMultispectalSize*maxNumbClass); i++){
        this->M[i]=0.0f;
    }
    this->V= new float [maxMultispectalSize*maxMultispectalSize*maxNumbClass];
    for(int i=0; i<(maxMultispectalSize*maxMultispectalSize*maxNumbClass); i++){
        this->V[i]=0.0f;
    }

    this->Expec=NULL;
    this->ShortPrior=NULL;
    this->S2L=NULL;
    this->L2S=NULL;
    this->CurrSizes=NULL;
    this->reg_factor=1.1f;

    this->maxIteration=100;
    this->minIteration=4;
    this->verbose_level=0;
    this->loglik=2.0;
    this->oldloglik=1.0;


    this->maskImageStatus=false;
    this->Mask=NULL; // pointer to external
    this->numElementsMasked=0;

    this->priorsStatus=false;
    this->Priors=NULL;

    this->mrfStatus=false;
    this->MRF_strength=0.0f;
    this->MRF=NULL;
    this->MRFBeta=NULL;
    this->MRFTransitionMatrix=NULL;

    this->Outlierness=NULL;
    this->OutliernessUSE=NULL;
    this->outliernessStatus=false;
    this->outliernessThreshold=0.0f;
    this->outliernessRatio=0.01f;

    this->BiasField_status=false;
    this->BiasField_order=0;
    this->BiasField=NULL;
    this->BiasField_coeficients=NULL;
    this->BiasField_ratio=0;

    this->LoAd_phase=-1;
    this->pvModelStatus=false;
    this->SG_deli_status=false;
    this->Relax_status=false;
    this->Relax_factor=0;
    this->RelaxGaussKernelSize=0;
    this->mapStatus=0;
    this->MAP_M=NULL;
    this->MAP_V=NULL;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

seg_EM_old::~seg_EM_old()
{
    if(this->Expec!=NULL){
        delete[] this->Expec;
    }
    this->Expec=NULL;

    if(this->ShortPrior!=NULL){
        delete [] this->ShortPrior;
    }
    this->ShortPrior=NULL;

    if(this->BiasField!=NULL){
        delete [] this->BiasField;
    }
    this->BiasField=NULL;

    if(this->BiasField_coeficients!=NULL){
        delete [] this->BiasField_coeficients;
    }
    this->BiasField_coeficients=NULL;

    if(this->Outlierness!=NULL){
        delete [] this->Outlierness;
    }

    if(this->mrfStatus){
        delete [] this->MRF;
        delete [] this->MRFTransitionMatrix;
    }
    this->MRF=NULL;
    this->MRFTransitionMatrix=NULL;

    if(this->maskImageStatus){
        delete [] this->L2S;
        delete [] this->S2L;

    }

    if(this->CurrSizes!=NULL)
        delete [] this->CurrSizes;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetInputImage(nifti_image *r)
{
    this->InputImage = r;
    this->inputImage_status = true;
    // Size
    this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1)+(int)((r->nt)>1)+(int)((r->nu)>1);
    this->nx=r->nx;
    this->ny=r->ny;
    this->nz=r->nz;
    this->nt=r->nt;
    this->nu=r->nu;
    this->dx=r->dx;
    this->dy=r->dy;
    this->dz=r->dz;
    this->numElements=r->nz*r->ny*r->nx;
    if(this->nx==1 ||this->ny==1){
        cout<<"Error: The segmentation algorithm only takes 2D, 3D and 5D images. 2D images should be on the XY plane"<<endl;
        return 1;
    }

    return 0;
}


int seg_EM_old::SetMAP( float *M, float* V)
{
    this->mapStatus=true;
    this->MAP_M=new float[this->numberOfClasses];
    this->MAP_V=new float[this->numberOfClasses];


    for(int i=0;i<this->numberOfClasses;i++){
        this->MAP_M[i]=M[i];
        this->MAP_V[i]=V[i];
    }

    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetPriorImage(nifti_image *r)
{
    this->Priors = r;
    this->priorsStatus = true;
    // Size
    this->dimentions=(int)((r->nx)>1)+(int)((r->ny)>1)+(int)((r->nz)>1);
    if(this->nx==r->nx && this->ny==r->ny && this->nz==r->nz){
        return 0;
    }
    else{
        cout << "ERROR: Priors have wrong dimentions" << endl;
        return 1;
    }


}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetFilenameOut(char *f)
{
    this->FilenameOut = f;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetMaskImage(nifti_image *f)
{
    this->Mask = f;
    this->maskImageStatus = true;

    if(Mask->datatype!=DT_BINARY){
        seg_convert2binary(Mask,0.5f);
    }

    if(Mask->nt>1 ||Mask->nt>2){
        cout << "ERROR: Mask has wrong dimentions" << endl;
        this->Mask->nt=1;
        this->Mask->nu=1;
        this->Mask->ndim=3;
    }

    bool * MaskDataPtr = static_cast<bool *>(Mask->data);
    this->numElementsMasked=0;
    for(int i=0; i<this->numElements; i++,MaskDataPtr++){
        if((*MaskDataPtr)>0){
            this->numElementsMasked++;
        }
    }
    if(this->nx==f->nx && this->ny==f->ny && this->nz==f->nz){
        return 0;
    }
    else{
        cout << "ERROR: Mask has wrong dimentions" << endl;
        return 1;
    }

    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetVerbose(unsigned int verblevel)
{
    if(verblevel>2){
        this->verbose_level=2;
    }
    else{
        this->verbose_level=verblevel;
    }
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetRegValue(float reg)
{
    if(reg>=1){
        this->reg_factor=reg;
    }
    else{
        this->reg_factor=1;
        cout << "Non valid regularization value. Will assume -reg 1."<<endl;
    }
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetMaximalIterationNumber(unsigned int numberiter)
{
    this->maxIteration=numberiter;
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetMinIterationNumber(unsigned int numberiter)
{
    this->minIteration=numberiter;
    this->checkpoint_iter=this->minIteration;
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::SetAprox(bool aproxval)
{
    this->aprox=aproxval;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::Turn_MRF_ON(float strength)
{
    this->mrfStatus=true;
    this->MRF_strength=strength;
    this->MRFTransitionMatrix=new float[this->numberOfClasses*this->numberOfClasses]();
    Create_diagonal_MRF_transitionMatrix();


    if(this->maskImageStatus){
        this->MRF=new float[this->numElementsMasked*this->numberOfClasses*this->nu*this->nt]();
        for(int i=0; i<(this->numElementsMasked*this->numberOfClasses); i++){
            MRF[i]=(float)(1.0);
        }
    }
    else{
        this->MRF=new float[this->numElements*this->numberOfClasses*this->nu*this->nt]();
        for(int i=0; i<(this->numElements*this->numberOfClasses); i++){
            MRF[i]=(float)(1.0);
        }
    }

    Create_diagonal_MRF_transitionMatrix();
    return 0;
}


int seg_EM_old::Turn_Relaxation_ON(float relax_factor,float relax_gauss_kernel){
    this->Relax_status=true;
    this->Relax_factor=relax_factor;
    this->RelaxGaussKernelSize=relax_gauss_kernel;
    return 1;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::Turn_BiasField_ON(int powerOrder,float ratiothresh)
{
    this->BiasField_status=true;
    this->BiasField_order=powerOrder;
    this->BiasField_coeficients = new segPrecisionTYPE[((powerOrder+1)*(powerOrder+2)/2*(powerOrder+3)/3)*this->nu*this->nt]();
    this->BiasField_ratio=ratiothresh;
    if(this->maskImageStatus){
        this->BiasField = new segPrecisionTYPE[this->numElementsMasked*this->nu*this->nt]();
    }
    else{
        this->BiasField = new segPrecisionTYPE[this->numElements*this->nu*this->nt]();
    }
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::Create_diagonal_MRF_transitionMatrix()
{
    for(int i=0;i<this->numberOfClasses;i++){
        for(int j=0;j<this->numberOfClasses;j++){
            if(j==i){
                this->MRFTransitionMatrix[i+j*this->numberOfClasses]=0;
            }
            else{
                this->MRFTransitionMatrix[i+j*this->numberOfClasses]=this->MRF_strength;
            }
        }
    }
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::Create_CurrSizes()
{
    this->CurrSizes = new ImageSize [1]();
    CurrSizes->numel=(int)(this->nx*this->ny*this->nz);
    CurrSizes->xsize=this->nx;
    CurrSizes->ysize=this->ny;
    CurrSizes->zsize=this->nz;
    CurrSizes->usize=(this->nu>1)?this->nu:1;
    CurrSizes->tsize=(this->nt>1)?this->nt:1;
    CurrSizes->numclass=this->numberOfClasses;
    CurrSizes->numelmasked=0;
    CurrSizes->numelbias=0;
    return 0;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_EM_old::OutliernessON(float in_OutliernessThreshold, float ratio){

    this->outliernessStatus=true;
    this->outliernessThreshold=in_OutliernessThreshold;
    this->outliernessRatio=ratio;
    return 0;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
void seg_EM_old::RunMaximization()
{

    int verbose=this->verbose_level;
    if(this->verbose_level>0)
    {
        cout<< "Optimising Mixture Model Parameters" << endl;
        flush(cout);
    }
    bool OutliernessFlag=(Outlierness==NULL)?0:1;

    int numel_masked=this->numElementsMasked;
    int num_class=this->numberOfClasses;
    int Expec_offset[maxNumbClass];
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=cl*numel_masked;
    }

    segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;
    segPrecisionTYPE * OutliernessTmpPTR = (segPrecisionTYPE *) this->Outlierness;
    segPrecisionTYPE * InputImageTmpPtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
    segPrecisionTYPE * BiasFieldTmpPTR= (segPrecisionTYPE *) this->BiasField;

#ifdef _OPENMP
#pragma omp parallel for shared(InputImageTmpPtr,BiasFieldTmpPTR,OutliernessTmpPTR)
#endif
    // ***********
    // For each class, get all the temporary pointers
    for( int cl=0; cl<num_class; cl++)
    {
        int * S2L_PTR = (int *) this->S2L;
        segPrecisionTYPE * ExpectationPTR = (segPrecisionTYPE *) ExpectationTmpPTR;
        segPrecisionTYPE * OutliernessPTR = (segPrecisionTYPE *) OutliernessTmpPTR;
        segPrecisionTYPE * InputImagePtr = (segPrecisionTYPE *) InputImageTmpPtr;
        segPrecisionTYPE * BiasFieldPTR= (segPrecisionTYPE *) BiasFieldTmpPTR;
        segPrecisionTYPE * T1_PTR2= (segPrecisionTYPE *) InputImageTmpPtr;
        segPrecisionTYPE * BiasField_PTR2= (segPrecisionTYPE *) BiasFieldTmpPTR;

        // MEAN
        // For each multispectral data (or each time point), estimate the mean vector. This involves doing a weighted sum (tempsum/SumPriors) of the observed intensities, weighted by the responsabilities Expec_PTR, the OutliernessPTR. The Estimated mean uses the bias field corrected intensities (T1_PTR-BiasField_PTR)
        for(long Multispec=0; Multispec<this->nu; Multispec++)
        {
            ExpectationPTR=(segPrecisionTYPE *) &Expec[Expec_offset[cl]];
            OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];
            S2L_PTR = (int *) this->S2L;

            InputImagePtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
            InputImagePtr = &InputImagePtr[Multispec*this->numElements];
            segPrecisionTYPE tempsum=(segPrecisionTYPE)0.0;
            segPrecisionTYPE SumPriors=(segPrecisionTYPE)0.0;
            // First, it estimates the weights tempsum and SumPriors.
            if(OutliernessFlag)
            {
                if(BiasField!=NULL)
                {
                    BiasFieldPTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,OutliernessPTR++,BiasFieldPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR)*(*OutliernessPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR)*(*OutliernessPTR);
                        }
                    }
                }
            }
            else
            {
                if(BiasField!=NULL)
                {
                    BiasFieldPTR= &BiasField[Multispec*numel_masked];
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,BiasFieldPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR));
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<numel_masked; i++, ExpectationPTR++,S2L_PTR++)
                    {
                        segPrecisionTYPE current_value=(*ExpectationPTR)*(InputImagePtr[(*S2L_PTR)]);
                        if(current_value==current_value)
                        {
                            tempsum+=current_value;
                            SumPriors+=(*ExpectationPTR);
                        }
                    }
                }

            }

            // Second, estimates the actual mean M, the variance V
            if(SumPriors==SumPriors && SumPriors>0)
            {
                if(this->mapStatus && this->MAP_M!=NULL)
                {
                    // The mean M
                    M[cl*(this->nu)+Multispec]=(tempsum/SumPriors/powf(V[cl*(this->nu)+Multispec],2)+this->MAP_M[cl*(this->nu)+Multispec]/powf(this->MAP_V[cl*(this->nu)+Multispec],2))/(1/powf(V[cl*(this->nu)+Multispec],2)+1/powf(this->MAP_V[cl*(this->nu)+Multispec],2));
                }
                else
                {
                    M[cl*(this->nu)+Multispec]=tempsum/SumPriors;

                }

                // The Covariances V, which require looping through all the multispectral chanels again.
                for(long Multispec2=Multispec; Multispec2<this->nu; Multispec2++)
                {
                    S2L_PTR = (int *) this->S2L;

                    InputImagePtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    InputImagePtr =&InputImagePtr[Multispec*this->numElements];

                    T1_PTR2 = static_cast<segPrecisionTYPE *>(this->InputImage->data);
                    T1_PTR2 =&T1_PTR2[Multispec2*this->numElements];
                    segPrecisionTYPE tmpM=this->M[cl*this->nu+Multispec];
                    segPrecisionTYPE tmpM2=this->M[cl*this->nu+Multispec2];
                    tempsum=0;
                    ExpectationPTR=&Expec[Expec_offset[cl]];
                    OutliernessPTR=(segPrecisionTYPE *) &Outlierness[Expec_offset[cl]];

                    if(BiasField!=NULL)
                    {
                        BiasFieldPTR=&BiasField[Multispec*numel_masked];
                        BiasField_PTR2=&BiasField[Multispec2*numel_masked];
                        if(OutliernessFlag)
                        {
                            for (int i=0; i<numel_masked; i++,ExpectationPTR++,BiasFieldPTR++,OutliernessPTR++,BiasField_PTR2++,S2L_PTR++)
                            {

                                segPrecisionTYPE currentValue=(*ExpectationPTR) * (*OutliernessPTR)*(InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(currentValue==currentValue)
                                {
                                    tempsum+=currentValue;
                                }
                            }
                        }
                        else
                        {
                            for (int i=0; i<numel_masked; i++,ExpectationPTR++,BiasFieldPTR++,BiasField_PTR2++,S2L_PTR++)
                            {
                                segPrecisionTYPE currentValue=(*ExpectationPTR) * (InputImagePtr[(*S2L_PTR)]+(*BiasFieldPTR)-tmpM) * (T1_PTR2[(*S2L_PTR)]+(*BiasField_PTR2)-tmpM2);
                                if(currentValue==currentValue)
                                {
                                    tempsum+=currentValue;
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int i=0; i<numel_masked; i++,ExpectationPTR++,S2L_PTR++)
                        {
                            segPrecisionTYPE current_vaue=(*ExpectationPTR) * (InputImagePtr[(*S2L_PTR)]-tmpM) * (T1_PTR2[(*S2L_PTR)]-tmpM2);
                            if(current_vaue==current_vaue)
                            {
                                tempsum+=current_vaue;
                            }
                        }

                    }
                    if( (tempsum/SumPriors>0) && SumPriors>0  && (!isnan(tempsum/SumPriors)))
                    {
                        // assign tempsum/SumPriors to the uper triangular part
                        V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]=tempsum/SumPriors;
                        //if the value is not in a diagonal
                        if(Multispec2!=Multispec)
                        {
                            // then copy to the botom triangular part
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]=V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu];
                            // and regularise the value by dividing by the reg_factor.
                            V[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]/=reg_factor;
                            V[cl*this->nu*this->nu+Multispec2+Multispec*this->nu]/=reg_factor;
                        }
                    }
                }
            }
        }
    }

    // This section is only for printing.
    if(verbose>0)
    {
        for (int cl=0; cl<num_class; cl++)
        {
            if(this->nu==1)
            {
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else
            {
                cout<< "M["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<this->nu; Multispec++)
                {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*this->nu+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(long Multispec=0; Multispec<this->nu; Multispec++)
                {
                    if(Multispec>0)
                    {
                        cout<< "      ";
                    }
                    for(long Multispec2=0; Multispec2<this->nu; Multispec2++)
                    {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*this->nu*this->nu+Multispec*this->nu+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM_old::RunExpectation()
{
    if(this->ratio<(segPrecisionTYPE)(this->outliernessRatio) && this->outliernessStatus)
    {
        this->OutliernessUSE=this->Outlierness;
        if(this->verbose_level>0)
        {
            cout << "Updating Outlierness - LogRatio = "<<ratio<<endl;
        }
    }

    segPrecisionTYPE * IterPrior=NULL;
    if(this->mrfStatus)
    {
        IterPrior=this->MRF;
    }
    else
    {
        IterPrior=this->ShortPrior;
    }

    int numel_masked=this->numElementsMasked;
    int num_class=this->numberOfClasses;
    bool OutliernessFlag=(this->OutliernessUSE==NULL)?0:1;
    segPrecisionTYPE inv_v [maxNumbClass*maxMultispectalSize*maxMultispectalSize]= {0.0f};
    segPrecisionTYPE inv_sqrt_V_2pi [maxNumbClass]= {0.0f};

    int Expec_offset [maxNumbClass]= {0};

    // for each class, do
    for (int cl=0; cl<num_class; cl++)
    {
        Expec_offset[cl]=(int) cl*numel_masked;
        // if multimodal, do
        if(this->nu>1)
        {
            seg_Matrix <double> Vmat(this->nu,this->nu);
            for (long j2=0; j2<this->nu; j2++)
            {
                for (long i2=j2; i2<this->nu; i2++)
                {
                    Vmat.setvalue(i2,j2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                    Vmat.setvalue(j2,i2,(double)(this->V[i2+j2*this->nu+cl*this->nu*this->nu]));
                }
            }
            // Get the Gaussian normaliser
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*Vmat.determinant()));
            if (this->verbose_level>1)
            {
                cout<<endl<<"inv_sqrt_V_2pi["<< cl <<"]= "<< inv_sqrt_V_2pi[cl] << endl;
                flush(cout);
            }
            // Get the inverted covariance matrix
            Vmat.invert();
            double covarianceValue=0.0f;
            bool success;
            // Print if in debug, i.e. verbose_level==2
            if (this->verbose_level>1)
            {
                cout<<"inv_V["<< cl <<"]= ";
                flush(cout);
            }
            for (long j2=0; j2<this->nu; j2++)
            {
                if(this->verbose_level>1)
                {
                    if(j2!=0)
                    {
                        cout<< endl << "          ";
                    }
                }
                for(long i2=0; i2<this->nu; i2++)
                {
                    // copy data from Vmat to inv_v
                    Vmat.getvalue(i2,j2,covarianceValue,success);
                    inv_v[i2+j2*this->nu+cl*this->nu*this->nu]=(segPrecisionTYPE)(covarianceValue);
                    if(this->verbose_level>1)
                    {
                        cout<<inv_v[i2+j2*this->nu+cl*this->nu*this->nu]<< "\t";
                        flush(cout);
                    }
                }

            }
            if(this->verbose_level>1)
            {
                cout<< endl;
            }
        }
            // else, just get the gaussian normaliser and the inverse covariance determinant
        else
        {
            inv_sqrt_V_2pi[cl]=1/(sqrtf(2*M_PI*V[cl]));
            inv_v[cl]=1/V[cl];
        }
    }
    this->loglik=0;
    segPrecisionTYPE logliktmp=0.0f;


    segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;
    segPrecisionTYPE * OutliernessTmpPTR = (segPrecisionTYPE *) this->Outlierness;
    segPrecisionTYPE * InputImageTmpPtr = static_cast<segPrecisionTYPE *>(this->InputImage->data);
    segPrecisionTYPE * BiasFieldTmpPTR= (segPrecisionTYPE *) this->BiasField;

#ifdef _OPENMP
    segPrecisionTYPE * loglikthread = new segPrecisionTYPE [omp_get_max_threads()]();
    for(long i=0; i<(long)omp_get_max_threads(); i++)
        loglikthread[i]=0;

#pragma omp parallel for shared(ExpectationTmpPTR,loglikthread,InputImageTmpPtr,BiasFieldTmpPTR,OutliernessTmpPTR,IterPrior)
#endif
    // Now that we have the Gaussian normaliser and the inverted determinant of the covariance, estimate the Gaussian PDF. We do that independently per voxel.
    for (int i=0; i<numel_masked; i++)
    {
        segPrecisionTYPE * T1_PTR = (segPrecisionTYPE *)(InputImageTmpPtr);
        segPrecisionTYPE T1_Bias_corr[maxMultispectalSize];
        segPrecisionTYPE SumExpec=0.0f;

        // For each modality, estimate the bias field corrected intensity first
        for(long Multispec=0; Multispec<this->nu; Multispec++)
            T1_Bias_corr[Multispec]=(BiasFieldTmpPTR!=NULL)?(T1_PTR[this->S2L[i]+Multispec*this->numElements] + BiasFieldTmpPTR[i+Multispec*numel_masked]):(T1_PTR[this->S2L[i]+Multispec*this->numElements]);

        // Then, for each class and for each modality, iterate over all modalities and subtract their (Y-\Mu)/inv(|\Sigma|) in order to get the Mahalanobis distance.
        for (int cl=0; cl<num_class; cl++)
        {
            segPrecisionTYPE mahal=0.0f;
            for(long Multispec=0; Multispec<this->nu; Multispec++)
            {
                segPrecisionTYPE tmpT1_BC_minusM=(T1_Bias_corr[Multispec] - this->M[cl*(this->nu)+Multispec]);
                for(long Multispec2=0; Multispec2<this->nu; Multispec2++)
                {
                    mahal-=(0.5f)*(T1_Bias_corr[Multispec2] - M[cl*(this->nu)+Multispec2])*inv_v[cl*this->nu*this->nu+Multispec+Multispec2*this->nu]*tmpT1_BC_minusM;
                }
            }

            // If the outlierness threshold is used, estimate here the new outlierness values. Note that this value is different per class. The 0.01 are there for stability reasons.
            if(OutliernessFlag)
            {
                segPrecisionTYPE outvalue=(expf(mahal)+0.01)/(expf(mahal)+expf(-0.5*(this->outliernessThreshold*this->outliernessThreshold))+0.01);
                OutliernessTmpPTR[i+Expec_offset[cl]]=outvalue;
            }
            // Update the Expectation by calculating, Prior*GaussianPDF
            ExpectationTmpPTR[i+Expec_offset[cl]]=IterPrior[i+Expec_offset[cl]] * ( inv_sqrt_V_2pi[cl] * expf(mahal) ) ;
            // Update the normaliser
            SumExpec+=ExpectationTmpPTR[i+Expec_offset[cl]];
        }

        // If something went wrong, just set the expectation to 1/K
        if (SumExpec<=0.0 || SumExpec!=SumExpec)
        {
            for (int cl=0; cl<num_class; cl++)
            {
                ExpectationTmpPTR[i+Expec_offset[cl]]=(segPrecisionTYPE)(1)/(segPrecisionTYPE)(num_class);
            }

        }
            // If it worked, then normalise the expectations, thus obtaining a per voxel responsability
        else
        {

            for (int cl=0; cl<num_class; cl++)
            {
                ExpectationTmpPTR[i+Expec_offset[cl]]=ExpectationTmpPTR[i+Expec_offset[cl]]/SumExpec;
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

    // update the log likelihood
    loglik=logliktmp;

    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the MRF optimisation step
///
/// This function will run the MRF in 2D or 3D depending on the clique structure. A ND version is coming soon, providing a temporal MRF.
///
void seg_EM_old::RunMRF()
{
    if(this->nz>1)
    {
        RunMRF3D();
    }
    else if(this->nz==1)
    {
        RunMRF2D();
    }
    return;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the 2D MRF
///
/// Estimates the 2D MRF energy given the MRFTransitionMatrix. This is similar to this->RunMRF3D() but in 2D. \sa this->RunMRF3D()
///
void seg_EM_old::RunMRF2D()
{

    int numelmasked=this->numElementsMasked;
    int numclass=this->numberOfClasses;
    segPrecisionTYPE * G =this->MRFTransitionMatrix;

    // if the MRF optimisation is ON
    if(this->mrfStatus)
    {
        // Define all the pointers and varibles
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *) this->MRF;
        int * Long_2_Short_IndicesPtr = (int *) this->L2S;
        int col_size, indexCentre, indexWest, indexEast, indexSouth, indexNorth;
        int ix, iy,maxiy, maxix, neighbourclass;
        segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
        segPrecisionTYPE Gplane[maxNumbClass];
        segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
        register int currclass;
        unsigned int numelmasked_currclass_shift[maxNumbClass];
        col_size = (int)(this->nx);
        maxix = (int)(this->nx);
        maxiy = (int)(this->ny);


        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }

        // precompute the index shifts between modalities (performance reasons)
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;
        }

        // As it is 2D, iterate over all Y and X's
        int curr_short_centreindex;
        for (iy=1; iy<maxiy-1; iy++)
        {
            // This updates the indexCentre to the correct column
            indexCentre=(col_size*iy);
            for (ix=1; ix<maxix-1; ix++)
            {
                // indexcentre is incremented before because ix starts at 1
                indexCentre++;
                Sum_Temp_MRF_Class_Expect = 0;
                curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];

                if (curr_short_centreindex>=0)
                {
                    // Get the index of all the 4 neighbours
                    indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                    indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                    indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                    indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        // Get the sum of the expectations for each class over all the neighbours
                        Gplane[currclass] = 0.0;
                        Temp_MRF_Class_Expect[currclass] = 0.0;
                        Gplane[currclass]+=this->Expec[indexWest];
                        Gplane[currclass]+=this->Expec[indexEast];
                        Gplane[currclass]+=this->Expec[indexNorth];
                        Gplane[currclass]+=this->Expec[indexSouth];
                        // increment the indexes to shift them to the next class
                        if(currclass<numclass)
                        {
                            indexWest+=numelmasked;
                            indexEast+=numelmasked;
                            indexNorth+=numelmasked;
                            indexSouth+=numelmasked;
                        }
                    }
                    // After you have the probabilities, estimate exp(-Beta*U_MRF), with U_MRF beeing the sum over all classes of G*Gplane
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        for (neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                        {
                            Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass];
                        }

                        if(this->MRFBeta==NULL)
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                        }
                        else
                        {
                            Temp_MRF_Class_Expect[currclass] = exp(this->MRFBeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                        }
                        // Also estimate the normaliser
                        Sum_Temp_MRF_Class_Expect += Temp_MRF_Class_Expect[currclass];
                    }
                    // Normalise the MRF prob using the MRF_Exp/Sum_MRF_Exp
                    for (currclass=0; currclass<numclass; currclass++)
                    {
                        MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                    }
                }
            }
        }

    }
    return;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/// @brief Optimisation function taking care of the 3D MRF
///
/// Estimates the 3D MRF energy given the MRFTransitionMatrix. This is similar to this->RunMRF2D() but in 3D. \sa this->RunMRF2D()
///
void seg_EM_old::RunMRF3D()
{
    int numelmasked=this->numElementsMasked;
    int numclass=this->numberOfClasses;

    segPrecisionTYPE * G =this->MRFTransitionMatrix;
    segPrecisionTYPE * H =this->MRFTransitionMatrix;

    // if the MRF optimisation is ON
    if(this->mrfStatus)
    {
        segPrecisionTYPE * MRFpriorPtr = (segPrecisionTYPE *)this->MRF;
        int * Long_2_Short_IndicesPtr = (int *)this->L2S;
        int col_size, plane_size;
        int maxiy, maxix, maxiz;
        col_size = (int)(this->nx);
        plane_size = (int)(this->nx)*(this->ny);

        maxix = (int)(this->nx);
        maxiy = (int)(this->ny);
        maxiz = (int)(this->nz);

        if(verbose_level>0)
        {
            cout << "Optimising MRF"<<endl;
            flush(cout);
        }

        unsigned int numelmasked_currclass_shift[maxNumbClass];

        // precompute the index shifts between modalities (performance reasons)
        for(int i=0; i<numclass; i++)
        {
            numelmasked_currclass_shift[i]=i*numelmasked;

        }

        segPrecisionTYPE * ExpectationTmpPTR = (segPrecisionTYPE *) this->Expec;

#ifdef _OPENMP
#pragma omp parallel for shared(ExpectationTmpPTR,MRFpriorPtr,Long_2_Short_IndicesPtr)
#endif
        // As it is 3D, iterate over all Z, Y and X's
        for (int iz=1; iz<maxiz-1; iz++)
        {
            // Define all the pointers and varibles
            segPrecisionTYPE Sum_Temp_MRF_Class_Expect;
            register int currclass;
            segPrecisionTYPE Temp_MRF_Class_Expect[maxNumbClass];
            segPrecisionTYPE Gplane[maxNumbClass];
            segPrecisionTYPE Hplane[maxNumbClass];
            int indexCentre, indexWest, indexEast, indexSouth, indexNorth, indexTop, indexBottom;
            for (int iy=1; iy<maxiy-1; iy++)
            {
                indexCentre=(col_size*iy)+(plane_size*iz);
                for (int ix=1; ix<maxix-1; ix++)
                {
                    // indexcentre is incremented before because ix starts at 1
                    indexCentre++;
                    Sum_Temp_MRF_Class_Expect = 0;
                    int curr_short_centreindex=Long_2_Short_IndicesPtr[indexCentre];
                    if (curr_short_centreindex>=0)
                    {
                        // Get the index of all the 6 neighbours
                        indexWest=Long_2_Short_IndicesPtr[indexCentre-col_size]>-1?Long_2_Short_IndicesPtr[indexCentre-col_size]:0;
                        indexEast=Long_2_Short_IndicesPtr[indexCentre+col_size]>-1?Long_2_Short_IndicesPtr[indexCentre+col_size]:0;
                        indexNorth=Long_2_Short_IndicesPtr[indexCentre-1]>-1?Long_2_Short_IndicesPtr[indexCentre-1]:0;
                        indexSouth=Long_2_Short_IndicesPtr[indexCentre+1]>-1?Long_2_Short_IndicesPtr[indexCentre+1]:0;
                        indexBottom=Long_2_Short_IndicesPtr[indexCentre+plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre+plane_size]:0;
                        indexTop=Long_2_Short_IndicesPtr[indexCentre-plane_size]>-1?Long_2_Short_IndicesPtr[indexCentre-plane_size]:0;
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            // Get the sum of the expectations for each class over all the neighbours
                            Gplane[currclass] = 0.0;
                            Hplane[currclass] = 0.0;
                            Temp_MRF_Class_Expect[currclass] = 0.0;
                            Gplane[currclass]+=ExpectationTmpPTR[indexWest];
                            Gplane[currclass]+=ExpectationTmpPTR[indexEast];
                            Gplane[currclass]+=ExpectationTmpPTR[indexNorth];
                            Gplane[currclass]+=ExpectationTmpPTR[indexSouth];
                            Hplane[currclass]+=ExpectationTmpPTR[indexTop];
                            Hplane[currclass]+=ExpectationTmpPTR[indexBottom];
                            // increment the indexes to shift them to the next class

                            if(currclass<numclass)
                            {
                                indexWest+=numelmasked;
                                indexEast+=numelmasked;
                                indexNorth+=numelmasked;
                                indexSouth+=numelmasked;
                                indexTop+=numelmasked;
                                indexBottom+=numelmasked;
                            }
                        }
                        // After you have the probabilities, estimate exp(-Beta*U_MRF), with U_MRF beeing the sum over all classes of ( G*Gplane + H*Hplane )

                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            for (int neighbourclass=0; neighbourclass<numclass; neighbourclass++)
                            {
                                Temp_MRF_Class_Expect[currclass]-=G[currclass+(numclass)*neighbourclass]*Gplane[neighbourclass]+H[currclass+(numclass)*neighbourclass]*Hplane[neighbourclass];
                            }
                            if(this->MRFBeta==NULL)
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            else
                            {
                                Temp_MRF_Class_Expect[currclass] = exp(this->MRFBeta[curr_short_centreindex]*Temp_MRF_Class_Expect[currclass])*this->ShortPrior[curr_short_centreindex+numelmasked_currclass_shift[currclass]];
                            }
                            // Estimate the normaliser
                            Sum_Temp_MRF_Class_Expect+=Temp_MRF_Class_Expect[currclass];
                        }
                        // Normalise the MRF prob using the MRF_Exp/Sum_MRF_Exp
                        for (currclass=0; currclass<numclass; currclass++)
                        {
                            MRFpriorPtr[curr_short_centreindex+numelmasked_currclass_shift[currclass]]=(Temp_MRF_Class_Expect[currclass]/Sum_Temp_MRF_Class_Expect);
                        }
                    }
                }
            }
        }
    }
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::RunBiasField()
{
    if(this->BiasField_status){
        if((((this->loglik-this->oldloglik)/fabs(this->oldloglik))<(segPrecisionTYPE)(this->BiasField_ratio) && this->iter>3)||((segPrecisionTYPE)(this->BiasField_ratio)==0.0f)){
            if(this->maskImageStatus){
                if(this->nz>1){

                    BiasCorrection_mask(this->BiasField,this->BiasField_coeficients,this->InputImage,this->L2S,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
                }
                else{
                    BiasCorrection_mask2D(this->BiasField,this->BiasField_coeficients,this->InputImage,this->L2S,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
                }
            }
            else{
                BiasCorrection(this->BiasField,this->BiasField_coeficients,this->InputImage,Expec,this->OutliernessUSE,this->M,this->V,this->BiasField_order,CurrSizes,this->BiasField_status,this->verbose_level);
            }
        }
    }
    return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::RunPriorRelaxation()
{
    if(this->Relax_status){

        if(this->maskImageStatus){
            if((int)(this->verbose_level)>(int)(0)){
                cout << "Relaxing Priors"<< endl;
            }
            PriorWeight_mask(this->ShortPrior,this->Priors,this->Expec,this->RelaxGaussKernelSize,this->Relax_factor,this->S2L,this->L2S,CurrSizes,this->verbose_level);
        }
        else{
            if((int)(this->verbose_level)>(int)(0)){
                cout << "Relaxing Priors only available on masked images"<< endl;
            }
            // PriorWeight(this->ShortPrior,this->Priors,this->Expec,CurrSizes,this->verbose_level);
        }
    }
    return 1;
}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

int seg_EM_old::InitializeAndNormalizeImageAndPriors()
{
    if(this->maskImageStatus){
        if(this->priorsStatus)Normalize_NaN_Priors_mask(this->Priors,this->Mask,this->verbose_level);
        Normalize_Image_mask(this->InputImage,this->Mask,CurrSizes,this->verbose_level);

        this->S2L = Create_Short_2_Long_Matrix_from_NII(this->Mask,&(CurrSizes->numelmasked));
        this->L2S = Create_Long_2_Short_Matrix_from_NII(this->Mask);
        this->numElementsMasked=(CurrSizes->numelmasked);
    }
    else{
        if(this->priorsStatus)Normalize_NaN_Priors(this->Priors,this->verbose_level);

        Normalize_Image(this->InputImage,CurrSizes,this->verbose_level);
    }

    if(this->mapStatus){
        for(int i=0;i<this->numberOfClasses;i++){
            this->MAP_M[i]=logf(((this->MAP_M[i]-CurrSizes->rescale_min[0])/(CurrSizes->rescale_max[0]-CurrSizes->rescale_min[0]))+1)/0.693147181;;
            if(this->verbose_level>0){
                cout << "MAP_M["<<i<<"] = "<<this->MAP_M[i]<< endl;
            }
            this->M[i]=this->MAP_M[i];
            this->V[i]=1.0/this->numberOfClasses;
        }
    }
    return 0;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

void seg_EM_old::InitializeAndAllocate()
{
    if(this->priorsStatus)
    {
        register long numel=(int)(this->Mask->nx*this->Mask->ny*this->Mask->nz);
        register long numel_masked=0;

        bool * Maskptrtmp = static_cast<bool *> (this->Mask->data);
        for (long i=0; i<numel; i++, Maskptrtmp++)
        {
            *Maskptrtmp?numel_masked++:0;
        }
        int pluspv=(int)(this->pvModelStatus)*2;

        this->Expec = new segPrecisionTYPE [numel_masked*(this->numberOfClasses+pluspv)] ();
        segPrecisionTYPE * tempExpec= (segPrecisionTYPE *) Expec;
        this->ShortPrior = new segPrecisionTYPE [numel_masked*(this->numberOfClasses+pluspv)] ();
        segPrecisionTYPE * tempShortPrior= (segPrecisionTYPE *) ShortPrior;
        segPrecisionTYPE * PriorPTR = static_cast<segPrecisionTYPE *>(this->Priors->data);
        for(long cl=0; cl<this->numberOfClasses; cl++)
        {
            Maskptrtmp = static_cast<bool *> (this->Mask->data);;
            for (int i=numel; i--; Maskptrtmp++,PriorPTR++)
            {
                if(*Maskptrtmp)
                {
                    *tempExpec = *PriorPTR;
                    *tempShortPrior= *PriorPTR;
                    tempExpec++;
                    tempShortPrior++;
                }
            }
        }
    }
    else
    {
        this->InitializeMeansUsingIntensity();
        int tmpnumb_elem=0;
        if(this->maskImageStatus>0)
        {
            tmpnumb_elem=(this->numElementsMasked*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
        }
        else
        {
            tmpnumb_elem=(numElements*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
        }

        segPrecisionTYPE tmpnumb_clas=((this->numberOfClasses+(int)(this->pvModelStatus)*2));
        this->Expec=new segPrecisionTYPE [tmpnumb_elem] ();
        this->ShortPrior=new segPrecisionTYPE [tmpnumb_elem] ();
        for(int i=0; i<tmpnumb_elem; i++)
        {
            this->Expec[i]=1.0/tmpnumb_clas;
            this->ShortPrior[i]=1.0/tmpnumb_clas;
        }
        this->RunExpectation();

    }

    for (int cl=0; cl<this->numberOfClasses; cl++)
    {
        if(this->outliernessStatus)
        {
            int tmpnumb_elem=0;
            if(this->maskImageStatus>0)
            {
                tmpnumb_elem=(this->numElementsMasked*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
            }
            else
            {
                tmpnumb_elem=(numElements*(this->numberOfClasses+(int)(this->pvModelStatus)*2));
            }
            this->OutliernessUSE=NULL;
            this->Outlierness=new segPrecisionTYPE [tmpnumb_elem] ();
            for(int i=0; i<tmpnumb_elem; i++)
            {
                this->Outlierness[i]=1.0;
            }
        }
    }

    return;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
int seg_EM_old::InitializeMeansUsingIntensity()
{

    for(int ms=0; ms<this->nu; ms++){
        segPrecisionTYPE * Intensity_PTR = static_cast<segPrecisionTYPE *>(this->InputImage->data);
        bool * MaskDataPtr=NULL;
        if(this->maskImageStatus){
            MaskDataPtr = static_cast<bool *>(this->Mask->data);
        }
        int mycounter=0;
        float meanval=0.0;
        float variance=0.0;

        for(int i=0; i<this->numElements; i++){
            if(!this->maskImageStatus || MaskDataPtr[i]>0){
                mycounter++;
                meanval+=(Intensity_PTR[i+ms*this->numElements]);
            }
        }
        meanval=meanval/mycounter;

        for(int i=0; i<this->numElements; i++){
            if(!this->maskImageStatus || MaskDataPtr[i]>0 ){
                variance+=pow((meanval-Intensity_PTR[i+ms*this->numElements]),2);
            }
        }
        variance=variance/mycounter;
        int histogram[1001];
        for(int i=0;i<1000;i++){
            histogram[i]=0;
        }
        float tmpmax=-1e32f;
        float tmpmin=1e32f;



        for(int i=0; i<this->numElements; i++){
            if(!this->maskImageStatus || MaskDataPtr[i]>0){
                if(tmpmax<(int)(Intensity_PTR[i+ms*this->numElements])){
                    tmpmax=(int)(Intensity_PTR[i+ms*this->numElements]);
                }
                if(tmpmin>(int)(Intensity_PTR[i+ms*this->numElements])){
                    tmpmin=(int)(Intensity_PTR[i+ms*this->numElements]);
                }
            }
        }

        for(int i=0; i<this->numElements; i++){
            if(!this->maskImageStatus || MaskDataPtr[i]>0){

                int index4hist=(int)(1000.0f*(float)(Intensity_PTR[i+ms*this->numElements]-tmpmin)/(float)(tmpmax-tmpmin));
                if((index4hist>1000) & (index4hist<0)){
                    cout<< "error"<<endl;
                }
                histogram[(int)(1000.0*(float)(Intensity_PTR[i+ms*this->numElements]-tmpmin)/(float)(tmpmax-tmpmin))]++;
            }
        }


        for(int clas=0; clas<this->numberOfClasses; clas++){
            float tmpsum=0;
            int tmpindex=0;
            float percentile=((float)clas+1)/(this->numberOfClasses+1);
            for(int i=999;i>0;i--){
                tmpsum+=histogram[i];
                tmpindex=i;
                if((float)(tmpsum)>((1.0f-percentile)*(float)(mycounter))){
                    i=0;
                }
            }
            M[clas*CurrSizes->usize+ms]=float(tmpindex)*(tmpmax-tmpmin)/1000.0f+(tmpmin);
            V[clas*CurrSizes->usize*CurrSizes->usize+ms*CurrSizes->usize+ms]=variance/this->numberOfClasses/2;
        }
    }


    for (int cl=0; cl<this->numberOfClasses; cl++) {
        if(this->verbose_level>0){
            if(CurrSizes->usize==1){
                cout.fill('0');
                cout<< "M["<<(int)(cl)<<"]= "<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl])<<"\tV["<<(int)(cl)<<"]="<<setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl])<< endl;
                flush(cout);
            }
            else{

                cout<< "M["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(M[cl*CurrSizes->usize+Multispec])<<"\t";
                }
                cout<< endl;
                flush(cout);
                cout<< "V["<<(int)(cl)<<"]= ";
                for(int Multispec=0; Multispec<CurrSizes->usize; Multispec++) {
                    if(Multispec>0){
                        cout<< "      ";
                    }
                    for(int Multispec2=0; Multispec2<CurrSizes->usize; Multispec2++) {
                        cout<< setw(10)<<setprecision(7)<<left<<(segPrecisionTYPE)(V[cl*CurrSizes->usize*CurrSizes->usize+Multispec*CurrSizes->usize+Multispec2])<<"\t";
                    }
                    cout<< endl;
                }
                cout<< endl;
                flush(cout);
            }
        }
    }

    return 0;
}
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

float * seg_EM_old::GetMeans()
{
    float * OutM= new float [this->nu*this->numberOfClasses];
    int index=0;
    for(int j=0; j<(this->numberOfClasses); j++){
        for(int i=0; i<(this->nu); i++){
            //cout << "M["<<j<<"]="<<this->M[j+i*this->numberOfClasses]<<endl;
            float resize=exp((this->M[index++])*0.693147181)-1;
            OutM[j+i*this->numberOfClasses]=(resize*(CurrSizes->rescale_max[i]-CurrSizes->rescale_min[i])+CurrSizes->rescale_min[i]);
        }
    }

    return OutM;
}
/*
float * seg_EM_old::GetSTD()
{
  float * OutV= new float [this->nu *this->numberOfClasses*this->numberOfClasses];
  for(int i=0; i<(this->nu); i++){
      for(int j=0; j<(this->numberOfClasses*this->numberOfClasses); j++){
          float resize=exp((sqrt(this->V[j+i*this->numberOfClasses*this->numberOfClasses]))*0.693147181)-1;
          OutV[j+i*this->numberOfClasses*this->numberOfClasses]=(resize*(CurrSizes->rescale_max[i]-CurrSizes->rescale_min[i]));
          cout << (j+i*this->numberOfClasses*this->numberOfClasses) << " value is : " << OutV[j+i*this->numberOfClasses*this->numberOfClasses] << " resize :"<< resize<< " max:"<<CurrSizes->rescale_max[i]<<" min"<<CurrSizes->rescale_min[i]<<" diff: "<<CurrSizes->rescale_max[i]-CurrSizes->rescale_min[i]<<endl;
        }
    }
  return OutV;
}
*/

float * seg_EM_old::GetSTD()
{
    float * OutV= new float [this->nu *this->numberOfClasses*this->numberOfClasses];
    for(int i=0; i<(this->numberOfClasses); i++){
        for(int j=0; j<(this->nu*this->nu); j++){
            float resize=exp((sqrt(this->V[j+i*this->nu*this->nu]))*0.693147181)-1;
            int x = j / this->nu;//row
            int y = j - (this->nu*x);//col
            float rescale_1 = CurrSizes->rescale_max[x]-CurrSizes->rescale_min[x];
            float rescale_2 = CurrSizes->rescale_max[y]-CurrSizes->rescale_min[y];
            OutV[j+i*this->nu*this->nu]=(resize*(sqrt(rescale_1 * rescale_2)));
            //cout << (j+i*this->nu*this->nu) << " value is : " << OutV[j+i*this->nu*this->nu];
            //cout << " resize :"<< resize<< " x: "<< x;
            //cout <<" y: "<< y << endl;
        }

    }
    return OutV;
}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM_old::GetResult()
{
    nifti_image * Result=NULL;

    if(this->maskImageStatus){
        Result = Copy_Expec_to_Result_mask(this->Expec,this->S2L,this->InputImage,(char*)this->FilenameOut.c_str(),this->CurrSizes);
    }
    else{
        Result = Copy_Expec_to_Result(this->Expec,this->InputImage,(char*)this->FilenameOut.c_str(),this->CurrSizes);
    }
    return Result;

}

//saurabh added
nifti_image * seg_EM_old::GetResultNew(char * fileName)
{
    nifti_image * Result=NULL;

    if(this->maskImageStatus){
        Result = Copy_Correct_Expec_to_Result_mask(this->Expec,this->Outlierness,this->S2L,this->InputImage,fileName,this->CurrSizes);
    }
    else{
        Result = Copy_Expec_to_Result(this->Expec,this->InputImage,(char*)this->FilenameOut.c_str(),this->CurrSizes);
    }
    return Result;

}


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM_old::GetBiasCorrected(char * filename)
{
    nifti_image * Result=NULL;

    if(this->maskImageStatus){
        Result = Get_Bias_Corrected_mask(this->BiasField_coeficients,this->InputImage,filename,this->CurrSizes,this->BiasField_order);
    }
    else{
        Result = Get_Bias_Corrected(this->BiasField,this->InputImage,filename,this->CurrSizes);
    }
    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */

nifti_image * seg_EM_old::GetOutlierness(char * filename)
{
    nifti_image * Result = nifti_copy_nim_info(this->InputImage);
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    Result->datatype=DT_FLOAT32;
    Result->cal_max=1;
    nifti_set_filenames(Result,filename,0,0);
    nifti_update_dims_from_array(Result);
    nifti_datatype_sizes(Result->datatype,&Result->nbyper,&Result->swapsize);
    Result->data = (void *) calloc(Result->nvox, sizeof(segPrecisionTYPE));
    segPrecisionTYPE * Resultdata = static_cast<segPrecisionTYPE *>(Result->data);
    for(unsigned int i=0; i<Result->nvox; i++){Resultdata[i]=0;}



    if(this->maskImageStatus){
        for(int i=0; i<CurrSizes->numelmasked; i++){
            float currsum=0;
            for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
                //currsum+=this->Outlierness[i+(currclass)*CurrSizes->numElementsMasked]*Expec[i+(currclass)*CurrSizes->numElementsMasked];
                currsum+= this->Outlierness[i+(currclass)*CurrSizes->numelmasked] * this->Expec[i+(currclass)*CurrSizes->numelmasked];
            }
            Resultdata[S2L[i]]=1-currsum;
        }
        //cout << "I m in GetOutlierness with mask code..."<< endl;
    }
    else{
        int class_nvox=Result->nx*Result->ny*Result->nz;
        for(int i=0; i<CurrSizes->numel; i++){
            float currsum=0;
            for(int currclass=0; currclass<CurrSizes->numclass;currclass++){
                currsum+=this->Outlierness[i+(currclass)*class_nvox];
            }
            Resultdata[i]=1-currsum;
        }
    }


    return Result;

}

/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */


int *  seg_EM_old::Run_EM()
{
    time_t start,end;
    time(&start);
    if((int)(this->verbose_level)>(int)(0)){
        cout << "EM: Verbose level " << this->verbose_level << endl;
    }
    if(!priorsStatus){

    }

    if(this->CurrSizes==NULL) {
        Create_CurrSizes();
    }

    InitializeAndNormalizeImageAndPriors();
    InitializeAndAllocate();
    if((int)(this->verbose_level)>(int)(0)){
        cout << "Number of voxels inside the mask = " << this->numElementsMasked << endl;
    }
    //**************
    // EM Algorithm
    //**************
    //bool MRFreset=0;
    this->iter=0;
    bool out= true;

    while (out) {

        if(this->verbose_level>0){
            cout << endl << "*******************************" << endl;
            cout << "Iteration " << iter << endl;
        }

        // Iterative Components - EM, MRF, Bias Correction

        //RunMaximization
        /*
         * No Difference
         * */
        this->RunMaximization();
        //Expectation
        /*
         * Very Small Difference
         * MAX(ABS(X)): 0.000378907
         * SUM(ABS(X)): 16.314
         * SUM(X): 5.35529e-05
         * */
        this->RunExpectation();
        //MRF
        /*
         * No Additional Difference
         * */
        this->RunMRF();
        //Bias Correction
        this->RunBiasField();
        //Update Weight
        this->RunPriorRelaxation();

        // Print LogLik depending on the verbose level
        if(this->verbose_level>0 && this->iter>0)
        {
            if(iter>0)
            {
                if ((this->loglik-this->oldloglik)/fabs(this->oldloglik)>0 && (this->loglik-this->oldloglik)/fabs(this->oldloglik)<100)
                {
                    cout<< "Loglik = " << setprecision(7)<<this->loglik <<
                        " : Ratio = " << (this->loglik-this->oldloglik)/fabs(this->oldloglik) << endl;
                }
                else
                {
                    cout<< "Loglik = " << setprecision(7)<<this->loglik << endl;
                }
            }
            else
            {
                cout<< "Initial Loglik = " << setprecision(7) <<this->loglik << endl;
            }
        }
        // Preform Exit
        if((((this->loglik-this->oldloglik)/fabs(this->oldloglik))<(segPrecisionTYPE)(0.0005) && this->iter > this->checkpoint_iter) || iter>=this->maxIteration || (isinf(this->loglik) && this->iter>3)){
            out=false;
        }
        this->ratio=((this->loglik-this->oldloglik)/fabs(this->oldloglik));
        // Update LogLik
        this->oldloglik=this->loglik;
        iter++;
    }

    time(&end);

    if(this->verbose_level>0){
        int minutes = (int)floorf(float(end-start)/60.0f);
        int seconds = (int)(end-start - 60*minutes);
        cout << "Finished in "<<minutes<<"min "<<seconds<<"sec"<< endl;
    }
    return 0;
}

#endif