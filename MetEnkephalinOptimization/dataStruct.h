#if !defined(ZGJ_INFO_ZJUT_EDU_CN_STRUCT)
#define ZGJ_INFO_ZJUT_EDU_CN_STRUCT

typedef	double (*ObjectiveFunction)(int dim, double *x);

typedef enum _SeedStatus {DEUnknown = 0, DETrue = 1, DEFalse = 2} SeedStatus;

struct	FunctionStruct {
	ObjectiveFunction objectiveFunction;	
	double	optimum;							 
	double	lowerBound;
};

#endif

