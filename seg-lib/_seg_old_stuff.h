
void Gaussian_Filter_Short_4D_old(segPrecisionTYPE *ShortData,
                                 int *S2L,
                                 int *L2S,
                                 segPrecisionTYPE gauss_std,
                                 ImageSize *CurrSizes,
                                 int class_with_CSF);

void get_xyz_pow_int(segPrecisionTYPE xpos,
                     segPrecisionTYPE ypos,
                     segPrecisionTYPE zpos,
                     segPrecisionTYPE currxpower[10],
                     segPrecisionTYPE currypower[10],
                     segPrecisionTYPE currzpower[10],
                     int maxorder);


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
                         int verbose_level);

nifti_image * Get_Bias_Corrected(float * BiasField,
                                 nifti_image * T1,
                                 char * filename,
                                 ImageSize * CurrSizes);

nifti_image * Get_Bias_Corrected_mask(float * BiasFieldCoefs,
                                      nifti_image * T1,
                                      char * filename,
                                      ImageSize * CurrSizes,
                                      int biasOrder);

int Gaussian_Filter_4D(segPrecisionTYPE * LongData,
                       segPrecisionTYPE gauss_std,
                       ImageSize * CurrSizes);