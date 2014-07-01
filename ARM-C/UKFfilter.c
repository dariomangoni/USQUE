#include "UKFfilter.h"

/*
 * Global variables
 */
extern float32_t dt;
extern float32_t omega_meas[3];
extern float32_t acc_meas[3];
extern float32_t magn_meas[3];
extern float32_t gps_meas[6];
extern float32_t state[N+1];
extern float32_t sonic_altitude;
extern float32_t gps_init[3];


/*
 *	Filter variables
 */
float32_t lambda		= 1.0f;
float32_t a				= 1.0f;
float32_t earth_radius	= 6362697.0f;
float32_t magn_NORTH[3] = {0.48528, 0.01856, -0.87416}; // vector measured when QUAD points north


float32_t vel_abs;
float32_t vel_angle;
uint8_t row, col, row2, col2;
float32_t fpow2, f;
float32_t one_minus_apow2;
float32_t fadq, N_plus_lambda;
float32_t scarto_var;
float32_t normsigmaXpow2;

float32_t x_post		[N];
float32_t x_pre			[N];
float32_t y_pre			[N_obs];
float32_t y_post		[N_obs];
float32_t K_innov		[N];
float32_t sigmaYmagn	[3];
float32_t gps_diff		[3];
float32_t gps_vel		[2];
float32_t innovation	[N_obs];

float32_t dq			[4]		[N_sigma];
float32_t q				[4]		[N_sigma];
float32_t sigmaX		[N]		[N_sigma];
float32_t sigmaXscarto	[N]		[N_sigma];
float32_t sigmaY		[N_obs]	[N_sigma];
float32_t sigmaYscarto	[N_obs]	[N_sigma];
float32_t K				[N]		[N_obs];
float32_t K_trasp		[N_obs]	[N];
float32_t PQ			[N]		[N];
float32_t Px			[N]		[N];
float32_t Px_pre		[N]		[N];
float32_t Px_temp		[N]		[N];
float32_t Pxy			[N]		[N_obs];
float32_t Py			[N_obs]	[N_obs];
float32_t Py_inverse	[N_obs]	[N_obs];
float32_t Py_temp		[N_obs]	[N_obs];
float32_t Q				[N]		[N];
float32_t R				[N_obs]	[N_obs];
float32_t PQ_temp		[N]		[N];
float32_t temp			[N_obs]	[N];


//ARM matrix
arm_matrix_instance_f32 Py_arm;
arm_matrix_instance_f32 Py_inverse_arm;
arm_matrix_instance_f32 Pxy_arm;
arm_matrix_instance_f32 K_arm;
arm_matrix_instance_f32 y_post_arm;
arm_matrix_instance_f32 y_pre_arm;
arm_matrix_instance_f32 innovation_arm;
arm_matrix_instance_f32 K_innov_arm;
arm_matrix_instance_f32 x_pre_arm;
arm_matrix_instance_f32 x_post_arm;
arm_matrix_instance_f32 K_trasp_arm;
arm_matrix_instance_f32 Px_arm;
arm_matrix_instance_f32 Q_arm;
arm_matrix_instance_f32 PQ_arm;
arm_matrix_instance_f32 Px_pre_arm;
arm_matrix_instance_f32 Px_temp_arm;
arm_matrix_instance_f32 magn_NORTH_arm;
arm_matrix_instance_f32 sigmaYmagn_arm;
arm_matrix_instance_f32 temp_arm;
arm_matrix_instance_f32 rotation_matrix_arm;
arm_matrix_instance_f32 rotation_matrix_inverse_arm;


/*
 * Filter initalization
 */

void UKF_init(void)
{

	f = 2*(a+1);
	fpow2 = f*f;
	one_minus_apow2 = (1-a*a);
	N_plus_lambda = (float32_t)(N+lambda);
	
	// matrix names cannot be used as pointers; must use &matrixname[0][0]; on the other hand vector names are pointers...
	arm_mat_init_f32(&Py_inverse_arm,		N_obs,	N_obs, &Py_inverse[0][0]);
	arm_mat_init_f32(&Pxy_arm,				N,		N_obs, &Pxy[0][0]);
	arm_mat_init_f32(&K_arm,				N,		N_obs, &K[0][0]);
	arm_mat_init_f32(&y_post_arm,			N_obs,	1,		y_post);
	arm_mat_init_f32(&y_pre_arm,			N_obs,	1,		y_pre);
	arm_mat_init_f32(&innovation_arm,		N_obs,	1,		innovation);
	arm_mat_init_f32(&K_innov_arm,			N,		1,		K_innov);
	arm_mat_init_f32(&x_pre_arm,			N,		1,		x_pre);
	arm_mat_init_f32(&x_post_arm,			N,		1,		x_post);
	arm_mat_init_f32(&Px_arm,				N,		N,		&Px[0][0]);
	arm_mat_init_f32(&Q_arm,				N,		N,		&Q[0][0]);
	arm_mat_init_f32(&K_trasp_arm,			N_obs,	N,		&K_trasp[0][0]);
	arm_mat_init_f32(&Px_pre_arm,			N,		N,		&Px_pre[0][0]);
	arm_mat_init_f32(&PQ_arm,				N,		N,		&PQ[0][0]);
	arm_mat_init_f32(&Px_temp_arm,			N,		N,		&Px_temp[0][0]);
	arm_mat_init_f32(&sigmaYmagn_arm,		3,		1,		sigmaYmagn);
	arm_mat_init_f32(&magn_NORTH_arm,		3,		1,		magn_NORTH);
	arm_mat_init_f32(&temp_arm,				N_obs,	N,		&temp[0][0]);
	

	float32_t diagonalPx[N] = {0.1,0.1,0.1,0.1,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};
	float32_t diagonalQ[N] = {1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12};
	float32_t diagonalR[N_obs] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};

	// initialize Px and Q
	int i;
	for (i=0;i<N;i++)
	{
		Px[i][i] = diagonalPx[i];
		Q[i][i] = diagonalQ[i];
	}

	// initialize R
	for (i=0;i<N_obs;i++)
	{
		R[i][i] = diagonalR[i];
	}

/*
 **************************
 * WARNING: Put here the code for activating timer interrupt that calls UKF_filterTimeUpdate() function!!!!
 **************************
 */

}

/*
 * Filter
 */

void UKF_filterTimeUpdate(void)
{
	// cholesky factorization (n+lambda)(Px+Q)
	arm_mat_add_f32	   (&Px_arm, &Q_arm,				  &PQ_arm);
	arm_mat_scale_f32  (&PQ_arm, N_plus_lambda,			  &PQ_arm);
	
	for (row=0;row<N;row++){
		for (col=0;col<N;col++){
			PQ_temp[row][col] = 0;
		}
	}

	cholesky(&PQ[0][0], &PQ_temp[0][0], N);

	// Initialize sigmapointX (equation 5)
	// mean in first column

	for (row=0;row<3;row++)
	{
		sigmaX[row][0] = 0.0f;
	}
	for (row=3;row<N;row++)
	{
		sigmaX[row][0] = state[row+1];
	}

	// other column
	for (row=0;row<N;row++)
	{
		for(col=0;col<N;col++)
		{
			sigmaX[row][col+1]	 = sigmaX[row][0] + PQ_temp[row][col];
			sigmaX[row][col+N+1] = sigmaX[row][0] - PQ_temp[row][col];
		}
	}

	// error quaternion (eq.33)
	dq[0][0] = 1.0f;
	dq[1][0] = 0.0f;
	dq[2][0] = 0.0f;
	dq[3][0] = 0.0f;
	for (col=1;col<N_sigma;col++)
	{
		normsigmaXpow2 = sigmaX[0][col]*sigmaX[0][col]+sigmaX[1][col]*sigmaX[1][col]+sigmaX[2][col]*sigmaX[2][col];
		dq[0][col] = (-a*normsigmaXpow2 + f*sqrtf( fpow2 + one_minus_apow2*normsigmaXpow2 ) ) / ( fpow2 + normsigmaXpow2 ); // can be simplified
		fadq = ( a + dq[0][col] ) / f;
		dq[1][col] = fadq * sigmaX[0][col];
		dq[2][col] = fadq * sigmaX[1][col];
		dq[3][col] = fadq * sigmaX[2][col];
	}

	// quaternion to be propagated (eq.32)
	q[0][0] = state[0];
	q[1][0] = state[1];
	q[2][0] = state[2];
	q[3][0] = state[3];
	for (col=1;col<N_sigma;col++)
	{
		quat_multiply ( &dq[0][col], &q[0][0], &q[0][col], N_sigma, N_sigma, N_sigma );
	}

	// propagation (eq.34)
	for (col=0;col<N_sigma;col++)
	{
		process_model( &q[0][col], &sigmaX[0][col], N_sigma, N_sigma );
	}

	// propagated quaterion error (eq.36)
	float32_t q_conj[4];
	q_conj[0] =	 q[0][0];
	q_conj[1] = -q[1][0];
	q_conj[2] = -q[2][0];
	q_conj[3] = -q[3][0];
	for (col=0;col<N_sigma;col++)
	{
		quat_multiply ( &q[0][col], q_conj, &dq[0][col], N_sigma, 1, N_sigma );
	}

	// propagated sigmapoint (eq.37,38)
	sigmaX[0][0] = 0.0f;
	sigmaX[1][0] = 0.0f;
	sigmaX[2][0] = 0.0f;
	for (col=1;col<N_sigma;col++)
	{
		fadq = f / ( a + dq[0][col] );
		sigmaX[0][col] = fadq * dq[1][col];
		sigmaX[1][col] = fadq * dq[2][col];
		sigmaX[2][col] = fadq * dq[3][col];
	}


	// predicted mean (eq.7)
	for (row=0; row<N; row++)
	{
		x_pre[row] = 0;
	}
	
	for (row=0;row<N;row++)
	{
		for(col=1;col<N_sigma;col++)
		{
			x_pre[row] += sigmaX[row][col];
		}
	}

	for (row=0;row<N;row++)
	{
		x_pre[row] = ( x_pre[row] * 0.5f + lambda*sigmaX[row][0] ) / N_plus_lambda;
	}

	// predicted state covariance (eq.8);
	// zeros in the matrix
	for (row=0;row<N;row++)
	{
		for(col=0;col<N;col++)
		{
			Px_pre[row][col] = 0.0f;
		}
	}
	
	for (col=1;col<N_sigma;col++)
	{
		// calculate scarto
		for(row=0;row<N;row++)
		{
			sigmaXscarto[row][col] = sigmaX[row][col] - x_pre[row];
		}
		
		// compose matrix
		for (col2=0;col2<N;col2++)
		{
			for (row2=0;row2<=col2;row2++)
			{
				Px_pre[row2][col2] += sigmaXscarto[row2][col]*sigmaXscarto[col2][col];
			}
		}
	}
	
	// add mean variance
	for(row=0;row<N;row++)
	{
		sigmaXscarto[row][0] = sigmaX[row][0] - x_pre[row];
	}

	for (col2=0;col2<N;col2++)
	{
		for (row2=0;row2<=col2;row2++)
		{
			Px_pre[row2][col2] = (0.5f * Px_pre[row2][col2] + lambda * sigmaXscarto[row2][0]*sigmaXscarto[col2][0]) / N_plus_lambda + Q [row2][col2];
			if (row2!=col2)
			{
				Px_pre[col2][row2] = Px_pre[row2][col2];
			}
		}
	}
	
}

void UKF_filterMeasurementUpdate(void)
{
	// observation model (eq.10,43)
	for (col=0;col<N_sigma;col++)
	{
		observation_model ( &q[0][col], &sigmaX[0][col], &sigmaY[0][col], N_sigma, N_sigma, N_sigma);
	}
	
	// observation mean (eq.9)
	for (row=0;row<N_obs;row++)
	{
		y_pre[row] = 0.0f; //inglobabile nel ciclo qui appena dopo
	}
	
	for (row=0;row<N_obs;row++)
	{
		for(col=1;col<N_sigma;col++)
		{
			y_pre[row] += sigmaY[row][col];
		}
	}

	for (row=0;row<N_obs;row++)
	{
		y_pre[row] = ( y_pre[row] * 0.5f + lambda*sigmaY[row][0] ) / N_plus_lambda;
	}
	
	// measurement/output covariance (eq.11,12) Pvv = Pyy for memory
	for (row=0;row<N_obs;row++)
	{
		for(col=0;col<N_obs;col++)
		{
			Py[row][col] = 0.0f;
		}
	}
	
	for (col=1;col<N_sigma;col++) // cycle upon sigmapoints
	{
		// calculate scarto
		for(row=0;row<N_obs;row++)
		{
			sigmaYscarto[row][col] = sigmaY[row][col] - y_pre[row];
		}
		
		// compose matrix
		for (col2=0;col2<N_obs;col2++)
		{
			for (row2=0;row2<=col2;row2++)
			{
				Py[row2][col2] += sigmaYscarto[row2][col]*sigmaYscarto[col2][col];
			}
		}
	}
	// add mean variance
	for(row=0;row<N_obs;row++)
	{
		sigmaYscarto[row][0] = sigmaY[row][0] - y_pre[row];
	}

	for (col2=0;col2<N_obs;col2++)
	{
		for (row2=0;row2<=col2;row2++)
		{
			Py[row2][col2] = (0.5f*Py[row2][col2] + lambda * sigmaYscarto[row2][0]*sigmaYscarto[col2][0]) / N_plus_lambda + R [row2][col2];
			Py_temp[row2][col2] = Py[row2][col2]; // for future inversion
			if (row2!=col2)
			{
				Py[col2][row2] = Py[row2][col2];
				Py_temp[col2][row2] = Py[col2][row2];
			}

		}
	}
	
	// cross-covariance (13)
	//supposed N>=N_obs
	for (row=0;row<N;row++)
	{
		for(col=0;col<N_obs;col++)
		{
			Pxy[row][col] = 0.0f;
		}
	}
	
	
	for (col=1;col<N_sigma;col++)
	{		
		// compose matrix
		for (col2=0;(col2<N_obs);col2++)
		{
			for (row2=0;(row2<N);row2++) // the last N-N_obs rows are not symmetric
			{
				scarto_var = sigmaXscarto[row2][col]*sigmaYscarto[col2][col];
				Pxy[row2][col2] += scarto_var;
			}
		}
	}
	
	// add mean variance
	for (col2=0;(col2<N_obs);col2++)
	{
		for (row2=0;(row2<N);row2++)
		{
			scarto_var = lambda*sigmaXscarto[row2][0]*sigmaYscarto[col2][0];
			Pxy[row2][col2] = (0.5f*Pxy[row2][col2] + scarto_var) / N_plus_lambda;
		}
	}
	
	// update state vector x_post
	
	// GPS conversion: degree-->meters
	gps_diff[0] = (gps_meas[0] - gps_init[0]) *pi / 180.0f * earth_radius * (float32_t)cos((float64_t) (gps_meas[0] * pi /180.0f) );
	gps_diff[1] = (gps_meas[1] - gps_init[1]) *pi / 180.0f * earth_radius;
	gps_diff[2] = gps_meas[2] - gps_init[2]; // not implemented
	
	// GPS conversion: km/h & heading --> m/s xy
	vel_abs = gps_meas[5] / 3.6f;
	vel_angle = gps_meas[4] * pi / 180.0f;
	
	gps_vel[0] = vel_abs * (float32_t)sin((float64_t)vel_angle);
	gps_vel[1] = vel_abs * (float32_t)cos((float64_t)vel_angle);
	
	// update state
	arm_mat_init_f32(&Py_arm, N_obs, N_obs, &Py_temp[0][0]);	//needs to stay here!
	
	y_post[0] = magn_meas[0];
	y_post[1] = magn_meas[1];
	y_post[2] = magn_meas[2];
	y_post[3] = gps_diff[0];
	y_post[4] = gps_diff[1];
	y_post[5] = sonic_altitude; // alternative to GPS altitude gps_diff[2]
	y_post[6] = gps_vel[0];
	y_post[7] = gps_vel[1];



	
	arm_mat_inverse_f32(&Py_arm,	 &Py_inverse_arm);
	arm_mat_mult_f32   (&Pxy_arm,	 &Py_inverse_arm, &K_arm);
	arm_mat_scale_f32  (&y_pre_arm,	 -1.0f,			  &y_pre_arm);
	arm_mat_add_f32	   (&y_post_arm, &y_pre_arm,	  &innovation_arm);
	arm_mat_mult_f32   (&K_arm,		 &innovation_arm, &K_innov_arm);
	arm_mat_add_f32	   (&x_pre_arm,	 &K_innov_arm,	  &x_post_arm);
	
	// error quaternion (eq.33)
	float32_t normadppow2 = x_post[0]*x_post[0] + x_post[1]*x_post[1] + x_post[2]*x_post[2];
	dq[0][0] = (-a*normadppow2 + f*sqrtf( fpow2 + one_minus_apow2*normadppow2 ) ) / ( fpow2 + normadppow2 ); // can be simplified
	fadq = (a + dq[0][0]) / f;
	dq[1][0] = fadq * x_post[0];
	dq[2][0] = fadq * x_post[1];
	dq[3][0] = fadq * x_post[2];
	

	//update state
	quat_multiply( &dq[0][0], &q[0][0], &state[0], N_sigma, N_sigma, 1);
	for (row=3;row<N;row++)
	{
		state[row+1] = x_post[row];
	}
	
	// Px update
	for(row=0;row<N_obs;row++)
	{
		for(col=0;col<N;col++)
		{
			K_trasp[row][col] = K[col][row];
		}
	}

	arm_mat_init_f32   (&Py_arm,	  N_obs, N_obs,	 &Py[0][0]);	//needs to stay here!

	for (row=0;row<N_obs;row++)
		{
			for(col=0;col<N;col++)
			{
				temp[row][col] = 0.0f;
			}
		}


	arm_mat_mult_f32   (&Py_arm,	  &K_trasp_arm,	 &temp_arm);
	arm_mat_mult_f32   (&K_arm,		  &temp_arm,  &Px_temp_arm);
	arm_mat_scale_f32  (&Px_temp_arm,  -1.0f,		 &Px_temp_arm);
	arm_mat_add_f32	   (&Px_pre_arm,  &Px_temp_arm,	 &Px_arm);
}


/*************************************
 * Auxiliary functions
 *************************************/

void quat_multiply ( float32_t* quat1, float32_t* quat2, float32_t* quatdest, int col1, int col2, int coldest )
{
	float32_t quat_in1[4] =
	{
			*quat1,
			*(quat1+col1),
			*(quat1+2*col1),
			*(quat1+3*col1)
	};

	float32_t quat_in2[4] =
	{
			*quat2,
			*(quat2+col2),
			*(quat2+2*col2),
			*(quat2+3*col2)
	};

	float32_t quat_temp[4]; // to avoid overwriting destination quat while quat it's used as quat 1 or quat 2
	quat_temp[0] = quat_in2[0] * quat_in1[0] - quat_in2[1] * quat_in1[1] - quat_in2[2] * quat_in1[2] - quat_in2[3] * quat_in1[3];
	quat_temp[1] = quat_in2[1] * quat_in1[0] + quat_in2[0] * quat_in1[1] - quat_in2[3] * quat_in1[2] + quat_in2[2] * quat_in1[3];
	quat_temp[2] = quat_in2[2] * quat_in1[0] + quat_in2[3] * quat_in1[1] + quat_in2[0] * quat_in1[2] - quat_in2[1] * quat_in1[3];
	quat_temp[3] = quat_in2[3] * quat_in1[0] - quat_in2[2] * quat_in1[1] + quat_in2[1] * quat_in1[2] + quat_in2[0] * quat_in1[3];

	*quatdest				= quat_temp[0];
	*(quatdest + coldest)	= quat_temp[1];
	*(quatdest + 2*coldest) = quat_temp[2];
	*(quatdest + 3*coldest) = quat_temp[3];
}

void cholesky(float32_t* source, float32_t* dest, int n)
{
	float32_t s;
	int i,j,k;
	for (i = 0; i < n; i++)
		for (j = 0; j < (i+1); j++)
		{
			s = 0.0f;
			for (k = 0; k < j; k++)
				s += dest[i*n+k] * dest[j*n+k];
			dest[i*n+j] = (i==j) ? sqrtf(source[i*n+i]-s) : ( (source[i*n+j]-s)/ dest[j*n+j] );
		}
}

void process_model( float32_t* quat, float32_t* sigmaX, int colquat, int colsigmaX)
{
	int i;

	// create temp variables
	float32_t* quatern[4];
	quatern[0] = quat;
	quatern[1] = quat+colquat;
	quatern[2] = quat+2*colquat;
	quatern[3] = quat+3*colquat;

	float32_t* sigmaXtemp[15];
	for (i=0;i<15;i++)
		sigmaXtemp[i] = sigmaX+i*colsigmaX;

	float32_t rotation_matrix[3][3] =
	{
	{0.5f - *quatern[2] * *quatern[2] - *quatern[3] * *quatern[3],
	*quatern[1] * *quatern[2] - *quatern[0] * *quatern[3],
	*quatern[1] * *quatern[3] + *quatern[0] * *quatern[2]},
	{*quatern[1] * *quatern[2] + *quatern[0] * *quatern[3],
	0.5f - *quatern[1] * *quatern[1] - *quatern[3] * *quatern[3],
	*quatern[2] * *quatern[3] - *quatern[0] * *quatern[1]},
	{*quatern[1] * *quatern[3] - *quatern[0] * *quatern[2],
	*quatern[2] * *quatern[3] + *quatern[0] * *quatern[1],
	0.5f - *quatern[1] * *quatern[1] - *quatern[2] * *quatern[2]}
	};
	
	arm_mat_init_f32(&rotation_matrix_arm,	3,		3,		&rotation_matrix[0][0]);
	arm_mat_scale_f32  (&rotation_matrix_arm, 2.0f, &rotation_matrix_arm);

	// calculate data for update
	float32_t omega_truedt2[3];
	omega_truedt2[0] = (omega_meas[0] - *sigmaXtemp[3])*dt*0.5f;
	omega_truedt2[1] = (omega_meas[1] - *sigmaXtemp[4])*dt*0.5f;
	omega_truedt2[2] = (omega_meas[2] - *sigmaXtemp[5])*dt*0.5f;

	float32_t acc_true[3];
	acc_true[0] = acc_meas[0] - *sigmaXtemp[12];
	acc_true[1] = acc_meas[1] - *sigmaXtemp[13];
	acc_true[2] = acc_meas[2] - *sigmaXtemp[14];

	float32_t pos_1dt[3];
	pos_1dt[0] = *sigmaXtemp[9]	 * dt;
	pos_1dt[1] = *sigmaXtemp[10] * dt;
	pos_1dt[2] = *sigmaXtemp[11] * dt;

	float32_t vel_1dt[3];
	float32_t gravity[3] = {0.0f, 0.0f, -gravityZ};
	arm_matrix_instance_f32 vel_1dt_arm;
	arm_matrix_instance_f32 acc_true_arm;
	arm_matrix_instance_f32 gravity_arm;

	arm_mat_init_f32   (&vel_1dt_arm,  3, 1, vel_1dt);
	arm_mat_init_f32   (&acc_true_arm, 3, 1, acc_true);
	arm_mat_init_f32   (&gravity_arm,  3, 1, gravity);

	// calculate acceleration in inertial frame
	arm_mat_mult_f32   (&rotation_matrix_arm, &acc_true_arm, &vel_1dt_arm);
	arm_mat_add_f32	   (&vel_1dt_arm,		  &gravity_arm,	 &vel_1dt_arm);
	arm_mat_scale_f32  (&vel_1dt_arm, dt,	  &vel_1dt_arm);

	// executing update
		// update position
	*sigmaXtemp[6] += pos_1dt[0];
	*sigmaXtemp[7] += pos_1dt[1];
	*sigmaXtemp[8] += pos_1dt[2];

		// update speed
	*sigmaXtemp[9]	+= vel_1dt[0];
	*sigmaXtemp[10] += vel_1dt[1];
	*sigmaXtemp[11] += vel_1dt[2];

		// update quaternion
	float32_t normomegadt = sqrtf(omega_truedt2[0] * omega_truedt2[0] + omega_truedt2[1] * omega_truedt2[1] + omega_truedt2[2] * omega_truedt2[2]);
	if (normomegadt>0.0f)
	{
		//float32_t deltanormquat = 1.0f - sqrtf(*quatern[0] * *quatern[0] + *quatern[1] * *quatern[1] + *quatern[2] * *quatern[2] + *quatern[3] * *quatern[3]);
		//float32_t jdt = 0.99999f;
		float32_t normomegadt_sin = sin(normomegadt) / normomegadt;
		float32_t omegaquat[4] =
		{
		//((float32_t)cos((float64_t)normomegadt) + jdt * deltanormquat),
		(cos(normomegadt)),
		omega_truedt2[0]*normomegadt_sin,
		omega_truedt2[1]*normomegadt_sin,
		omega_truedt2[2]*normomegadt_sin
		};

		quat_multiply ( &omegaquat[0], quat, quat, 1, 31, 31 );
	};

}

void observation_model (float32_t* quat, float32_t* sigmaX, float32_t* sigmaY, int colquat, int colsigmaX, int colsigmaY)
{
	float32_t* quatern[4];
	quatern[0] = quat;
	quatern[1] = quat+colquat;
	quatern[2] = quat+2*colquat;
	quatern[3] = quat+3*colquat;

	// define matrix for coordinates transform
	float32_t rotation_matrix_inverse[3][3] =
	{
	{0.5f - *quatern[2] * *quatern[2] - *quatern[3] * *quatern[3],
	*quatern[1] * *quatern[2] + *quatern[0] * *quatern[3],
	*quatern[1] * *quatern[3] - *quatern[0] * *quatern[2]},
	{*quatern[1] * *quatern[2] - *quatern[0] * *quatern[3],
	0.5f - *quatern[1] * *quatern[1] - *quatern[3] * *quatern[3],
	*quatern[2] * *quatern[3] + *quatern[0] * *quatern[1]},
	{*quatern[1] * *quatern[3] + *quatern[0] * *quatern[2],
	*quatern[2] * *quatern[3] - *quatern[0] * *quatern[1],
	0.5f - *quatern[1] * *quatern[1] - *quatern[2] * *quatern[2]}
	};

	arm_mat_init_f32   (&rotation_matrix_inverse_arm, 3, 3, &rotation_matrix_inverse[0][0]);
	arm_mat_scale_f32  (&rotation_matrix_inverse_arm, 2.0f, &rotation_matrix_inverse_arm);

	// transport magnetometer vector from earth to body frame
	arm_mat_mult_f32   (&rotation_matrix_inverse_arm, &magn_NORTH_arm, &sigmaYmagn_arm);

	*(sigmaY) = sigmaYmagn[0];
	*(sigmaY+colsigmaY) = sigmaYmagn[1];
	*(sigmaY+2*colsigmaY) = sigmaYmagn[2];

	// PosX PosY PosZ
	*(sigmaY+3*colsigmaY) = *(sigmaX + 6*colsigmaX);
	*(sigmaY+4*colsigmaY) = *(sigmaX + 7*colsigmaX);
	*(sigmaY+5*colsigmaY) = *(sigmaX + 8*colsigmaX);

	// SpeedX SpeedY
	*(sigmaY+6*colsigmaY) = *(sigmaX + 9*colsigmaX);
	*(sigmaY+7*colsigmaY) = *(sigmaX + 10*colsigmaX);
}

void setGPSinit(float32_t* gps_init)
{
  int cont = 5;
  int cont_temp = cont;
  float32_t mean_gps[3] = {0,0,0};

  while (cont>0)
	{
	  if (gps_meas[4]>((float32_t) 4))
		{
		  cont--;
		  mean_gps[0] += gps_meas[0]; //lat
		  mean_gps[1] += gps_meas[1]; //long
		  mean_gps[2] += gps_meas[2]; //alt --> low precision
		  //delay(1000); wait, it's not critical. doesn't need interrupt or hardware timer
		}
	}

  mean_gps[0] = mean_gps[0]/((float32_t) cont_temp);
  mean_gps[1] = mean_gps[1]/((float32_t) cont_temp);
  mean_gps[2] = mean_gps[2]/((float32_t) cont_temp);

  *gps_init		= mean_gps[0];
  *(gps_init+1) = mean_gps[1];
  *(gps_init+2) = mean_gps[2];
}
