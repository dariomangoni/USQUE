/*
 * UKF filter; following Crassidis
 */

//#define ARM_MATH_CM4
#include <arm_math.h>
#include <math_helper.h>

#include "stm32f4xx_conf.h"

#define N 			15
#define N_obs 		8
#define N_sigma 	2*N+1
#define gravityZ	9.81f
#define pi			3.14159265359f

/*
 * Function prototypes
 */
void UKF_init(void);
void UKF_filterTimeUpdate (void);
void UKF_filterMeasurementUpdate (void);

void process_model( float32_t* quat, float32_t* sigmaX, int colquat, int colsigmaX);
void observation_model (float32_t* quat, float32_t* sigmaX, float32_t* sigmaY, int colquat, int colsigmaX, int colsigmaY);

void quat_multiply ( float32_t* quat1, float32_t* quat2, float32_t* quatdest, int col1, int col2, int coldest );
void cholesky ( float32_t* source, float32_t* dest, int n);
void setGPSinit ( float32_t* gps_init);
