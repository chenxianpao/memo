#if !defined(ZGJ_INFO_ZJUT_EDU_CN_CONFIG)
#define ZGJ_INFO_ZJUT_EDU_CN_CONFIG

#define	 TCONDITION		 1				// �����㷨��ֹ��ʽ(0: ������ֹ; 1: ������ֹ)
#define	 EPSILON		 0.0001			// ������ֹ��ʽ��ֹ���
#define	 MAX_GENERATIONS 100			// ������ֹ��ʽ��ָ�������ܴ���
#define  OPTIMUM		 -11.7073		// �����������Ž�(�����֪������ͳ��)
#define	 RUN_NUM		 1				//�������д���





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
