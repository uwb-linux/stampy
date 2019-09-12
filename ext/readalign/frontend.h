#ifdef __cplusplus
extern "C"
#endif
unsigned int sensitiveAlign( char* iSeq1, 
			     char* iSeq2, 
			     int iLen1, 
			     int iLen2, 
			     char* pQuals1, 
			     int iInsNucPrior, 
			     int iGapOpen, 
			     int iGapExtend, 
			     char* line1, 
			     char* line2, 
			     char* line1q,
			     char* baq1,
			     int maxlikelihood, 
			     char* line3, 
			     char* line4,
			     char* line3q,
			     char* baq3,
			     int iStartMean,
			     int iStartSD,
			     int band,
			     double* amapq);



