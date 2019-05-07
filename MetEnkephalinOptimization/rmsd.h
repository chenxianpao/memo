#ifndef	RMSD_ZGJ_H
#define RMSD_ZGJ_H

/*
 *  ***************************************************************************
 *
 *  (c) 2010 G.J.Zhang
 *
 *  Implementation of the Kabsch algorithm to find the RMSD, 
 *  and the least-squares rotation matrix for a superposition between
 *  two sets of vectors.
 *
 *  This implementation is completely self-contained. No other dependencies.
 *
 *  ***************************************************************************
 *
 */

#ifndef PI
#define PI 3.14159265358979323846
#endif

class RMSD  
{
public:
	RMSD();
	virtual ~RMSD();

	/**************************************************************************
	 ***input parameters***
	 * -ref_xlist	: reference list of x, y, z coordinates, ref_xlist[][3]
	 * -mov_xlist	: trial list of x, y, z coordinates, mov_xlist[][3]
	 * -n_list		: the length of two lists
	 *-------------------------------------------------------------------------
	 ***output parameters***
	 * -mov_com		: the centre of mass of the mov list
	 * -mov_to_ref	: vector between the com of mov and ref	
	 * U			: the rotation matrix for least-squares
	 *-------------------------------------------------------------------------
	 ***retrun values***
	 * rmsd			: measures similarity between the vectors
	**************************************************************************/
	double calculate_rotation_rmsd(	double **ref_xlist,	// ref_xlist[][3]
									double **mov_xlist,	// mov_xlist[][3]
									int    n_list,
									double mov_com[3],
									double mov_to_ref[3],
									double U[3][3]);

	double fast_rmsd(double **ref_xlist, double** mov_xlist, int n_list);

private:
	double raw_rmsd (double **ref_xlist, double **mov_xlist, int n_list);
	double raw_translation_rmsd(double **ref_xlist, double **mov_xlist, int n_list);

public:
	void setup_rotation( double **ref_xlist,			// ref_xlist[][3]
                         double **mov_xlist,			// mov_xlist[][3]
                         int    n_list,
                         double mov_com[3],
                         double mov_to_ref[3],
                         double R[3][3],
                         double *E0);

	int calculate_rotation_matrix(double R[3][3],
                                  double U[3][3],
                                  double E0,
                                  double *residual);
	
	int diagonalize_symmetric(double matrix[3][3],
							  double eigen_vec[3][3],
                              double eigenval[3]);
	
	
	void normalize(double a[3]);
	void cross(double a[3], double b[3], double c[3]);
	double dot(double a[3], double b[3]);
	int	jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot);
};


#endif 
