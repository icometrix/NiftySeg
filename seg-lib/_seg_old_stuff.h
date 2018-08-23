int calcE(nifti_image * T1,
          segPrecisionTYPE * MRF,
          segPrecisionTYPE * Expec,
          double * loglik,
          segPrecisionTYPE * BiasField,
          segPrecisionTYPE * Outlierness,
          segPrecisionTYPE OutliernessThreshold,
          segPrecisionTYPE * M,
          segPrecisionTYPE * V,
          ImageSize * CurrSizes,
          int verbose);

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
               int verbose);


void MRFregularization_mask(const segPrecisionTYPE * Expec,
                            const segPrecisionTYPE * G,
                            const segPrecisionTYPE * H,
                            segPrecisionTYPE * MRFbeta,
                            segPrecisionTYPE * MRFprior,
                            segPrecisionTYPE * AtlasPrior,
                            int * Long_2_Short_Indices,
                            int * Short_2_Long_Indices,
                            ImageSize * CurrSizes,
                            bool MRFflag,
                            int verbose_level);

void MRFregularization(const segPrecisionTYPE * Expec,
                       const segPrecisionTYPE * G,
                       const segPrecisionTYPE * H,
                       segPrecisionTYPE * MRFbeta,
                       segPrecisionTYPE * MRFprior,
                       segPrecisionTYPE * AtlasPrior,
                       ImageSize * CurrSizes,
                       bool MRFflag,
                       int verbose_level);

void MRFregularization_mask2D(const segPrecisionTYPE * Expec,
                              const segPrecisionTYPE * G,
                              const segPrecisionTYPE * H,
                              segPrecisionTYPE * MRFbeta,
                              segPrecisionTYPE * MRFprior,
                              segPrecisionTYPE * AtlasPrior,
                              int * Long_2_Short_Indices,
                              int * Short_2_Long_Indices,
                              ImageSize * CurrSizes,
                              bool MRFflag,
                              int verbose_level);

void MRFregularization2D(const segPrecisionTYPE * Expec,
                         const segPrecisionTYPE * G,
                         const segPrecisionTYPE * H,
                         segPrecisionTYPE * MRFbeta,
                         segPrecisionTYPE * MRFprior,
                         segPrecisionTYPE * AtlasPrior,
                         ImageSize * CurrSizes,
                         bool MRFflag,
                         int verbose_level);

segPrecisionTYPE * Create_cArray_from_Prior_mask(nifti_image * Mask,
                                                 nifti_image * Priors,
                                                 long numclass,
                                                 bool PV_ON);

segPrecisionTYPE * Create_cArray_from_Prior(nifti_image * Priors,
                                            long numclass,
                                            bool PV_ON);

int PriorWeight_mask(float * ShortPrior,
                     nifti_image * Priors,
                     float * Expec,
                     float GaussKernelSize,
                     float RelaxFactor,
                     int * S2L,
                     int * L2S,
                     ImageSize * CurrSizes,
                     int verbose_level);

int Gaussian_Filter_Short_4D(segPrecisionTYPE * ShortData,
                             int * S2L,
                             int * L2S,
                             segPrecisionTYPE gauss_std,
                             ImageSize * CurrSizes,
                             int class_with_CSF);


int Normalize_NaN_Priors_mask(nifti_image * Priors,
                              nifti_image * Mask,
                              bool verbose);

int Normalize_Image_mask(nifti_image * T1,
                         nifti_image * Mask,
                         ImageSize * CurrSizes,
                         bool verbose);

int * Create_Long_2_Short_Matrix_from_NII(nifti_image * Mask);

int * Create_Short_2_Long_Matrix_from_NII(nifti_image * Mask,
                                          long *shortsize);

int Normalize_NaN_Priors(nifti_image * Priors,
                         bool verbose);

int Normalize_Image(nifti_image * T1,
                    ImageSize * CurrSizes,
                    bool verbose);

nifti_image * Copy_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                        int * Short_2_Long_Indices,
                                        nifti_image * T1,
                                        char * filename,
                                        ImageSize * CurrSizes);

//added by saurabh
//save p_ik * t_ik and not just p_ik as done by fnc Copy_Expec_to_Result_mask
nifti_image * Copy_Correct_Expec_to_Result_mask(segPrecisionTYPE * Expec,
                                                segPrecisionTYPE * Outlierness,
                                                int * Short_2_Long_Indices,
                                                nifti_image * T1,
                                                char * filename,
                                                ImageSize * CurrSizes);

nifti_image * Copy_Expec_to_Result(segPrecisionTYPE * Expec,
                                   nifti_image * T1,
                                   char * filename,
                                   ImageSize * CurrSizes);


inline segPrecisionTYPE pow_int(const segPrecisionTYPE x,
                                int exp);

void BiasCorrection_SPARCS(float * BiasField,
                           float * T1,
                           float * Expec,
                           float * Mask,
                           float * M,
                           float * V,
                           int biasOrder,
                           int nrOfClasses,
                           int aceletation_factor,
                           int xyzsize[3]);

void get_xyz_pow_int(segPrecisionTYPE xpos,
                     segPrecisionTYPE ypos,
                     segPrecisionTYPE zpos,
                     segPrecisionTYPE currxpower[10],
                     segPrecisionTYPE currypower[10],
                     segPrecisionTYPE currzpower[10],
                     int maxorder);

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
                    int verbose_level);

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

void get_xy_pow_int(segPrecisionTYPE xpos,
                    segPrecisionTYPE ypos,
                    segPrecisionTYPE currxpower[10],
                    segPrecisionTYPE currypower[10],
                    int maxorder);

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
                      int verbose_level);

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