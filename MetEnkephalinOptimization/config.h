#if !defined(ZGJ_INFO_ZJUT_EDU_CN_CONFIG)
#define ZGJ_INFO_ZJUT_EDU_CN_CONFIG

#define	 TCONDITION		 1				// 定义算法中止方式(0: 差异终止; 1: 代数终止)
#define	 EPSILON		 0.0001			// 差异中止方式中止误差
#define	 MAX_GENERATIONS 100			// 代数中止方式中指定运行总代数
#define  OPTIMUM		 -11.7073		// 定义问题最优解(如果已知，便于统计)
#define	 RUN_NUM		 1				//定义运行次数





#define	BIG_NUMBER		 1E25

#define	SINGLE_ENERGY_EVALUATE_SOURCE_FILE		"singleEnergyEvaluateSourceFile.inp"
#define	SINGLE_ENERGY_EVALUATE_SOURCE			"singleEnergyEvaluateSourceFile"
#define	SINGLE_ENERGY_EVALUATE_MAINOUT_FILE		"main_out.MET_SINGLE_EVALUATE"

#define	SINGLE_ENERGY_MINIMIZE_SOURCE_FILE		"singleEnergyMinimizeSourceFile.inp"
#define	SINGLE_ENERGY_MINIMIZE_SOURCE			"singleEnergyMinimizeSourceFile"
#define	SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE		"main_out.MET_SINGLE_MIN"
#define SINGLE_ENERGY_MINIMIZE_OUTO_FILE		"outo.MET_SINGLE_MIN"


#define	MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE	"metEnkephalinEnergyMultipleEvaluationSourceFile.inp"
#define	MULTIPLE_ENERGY_EVALUATION_INPUT_FILE	"outo.MultipleEnergyEvaluateInput"
#define	MULTIPLE_ENERGY_EVALUATION_RESULT_FILE	"main_out.MET_ENKEPHALIN_ENERGY_EVALUATION"			


#define	MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE	"metEnkephalinEnergyMultipleMinimizeSourceFile.inp"
#define MULTIPLE_ENERGY_MINIMIZE_RESULT_FILE	"main_out.MET_MULTI_MIN"
#define	MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE		"outo.MultipleEnergyMinimizeInput"


#define	LOG_FILE								"log"


#define	SINGLE_ENERGY_EVALUATION_SOURCE_FILE	"metEnkephalinEnergyEvaluationSourceFile.inp"
#define	SINGLE_ENERGY_EVALUATION_RESULT_FILE	"main_out.MET_ENKEPHALIN_ENERGY_EVALUATION"			


#endif
