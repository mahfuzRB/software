//    rbAPPmit: An Android front-end for the Certified Reduced Basis Method
//    Copyright (C) 2010 David J. Knezevic and Phuong Huynh
//
//    This file is part of rbAPPmit
//
//    rbAPPmit is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rbAPPmit is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rbAPPmit.  If not, see <http://www.gnu.org/licenses/>. 

package rb.app;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;

import android.content.Context;
import android.util.Log;

// This class provides the Online reduced
// basis functionality for quadratically nonlinear
// time-dependent problems.
// This class is modeled on the QNTransientRBSystem
// class in rbOOmit

public class QNTransientRBSystem extends TransientRBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "QNTransientRBSystem";

	// The tolerance for the Newton solver
	private double nonlinear_tolerance;

	// The maximum number of Newton iterations
	private int n_newton_steps;

	// The nominal lower bound for the stability factor
	private double nominal_rho_min;
	private double nominal_rho_max;

	// The RB data for the trilinear form
	protected double[][][] RB_trilinear_form;
	
	/**
	 * Vectors storing the residual representor inner products
	 * to be used in computing the residuals online.
	 */
	double[][][]   Fq_C_representor_norms;
	double[][][][] Mq_C_representor_norms;
	double[][][][] Aq_C_representor_norms;
	double[][][][] C_C_representor_norms;
	
	/**
	 * Member variable that stores the exponential growth factor
	 * of the error bound
	 */
	private double tau_prod_k;

	// Private member that we need in calling the SCM
	private int _N_in_RB_solve;
	
	/**
	 * Constructor.
	 */
	public QNTransientRBSystem(Context context) {
		super(context);
	}

	/**
	 * Perform online solve for current_params
	 * with the N basis functions. Overload this
	 * to solve the nonlinear RB system using
	 * Newton's method.
	 */
	@Override
	public double RB_solve(int N) {
		
		current_N = N;
		_N_in_RB_solve = N;
		
		// Initialize tau_prod_k
		tau_prod_k = 1.;

		if (N > get_n_basis_functions()) {
			throw new RuntimeException(
					"ERROR: N cannot be larger than the number "
					+ "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException(
			"ERROR: N must be greater than 0 in RB_solve");
		}

		double dt = get_dt();

		// First assemble the mass matrix
		RealMatrix RB_mass_matrix_N = new Array2DRowRealMatrix(N, N);
		for(int q_m=0; q_m<get_Q_m(); q_m++)
		{
			RB_mass_matrix_N = RB_mass_matrix_N.add(
					RB_M_q_vector[q_m].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply(eval_theta_q_m(q_m)) );
		}

		RealMatrix RB_RHS_Aq_matrix = new Array2DRowRealMatrix(N, N);

		RealMatrix Base_RB_LHS_matrix = RB_mass_matrix_N.scalarMultiply(1./dt);

		for(int q_a=0; q_a<get_Q_a(); q_a++) {
			double cached_theta_q_a = eval_theta_q_a(q_a);
			Base_RB_LHS_matrix = Base_RB_LHS_matrix.add(
					RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply( get_euler_theta()*cached_theta_q_a ) );
			RB_RHS_Aq_matrix   = RB_RHS_Aq_matrix.add( 
					RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply( -cached_theta_q_a ) );
		}

		// Set system time level to 0
		set_time_level(0);

		// Initialize a vector to store our current Newton iterate
		RealVector RB_u_bar = new ArrayRealVector(N);

		// and load the _k=0 data
		RB_solution = RB_u_bar;
		RB_temporal_solution_data[_k] = RB_u_bar; // Should use .copy() here!

		double error_bound_sum = 0.;

		// Set error bound at _k=0
		error_bound_all_k[_k] = Math.sqrt(error_bound_sum);

		// Compute the outputs and associated error bounds at _k=0
		for(int i=0; i<get_n_outputs(); i++) {
			RB_outputs_all_k[i][_k] = 0.;
			RB_output_error_bounds_all_k[i][_k] = 0.;
			for(int q_l=0; q_l<get_Q_l(i); q_l++) {
				RB_outputs_all_k[i][_k] +=
					eval_theta_q_l(i, q_l) * 
					(RB_output_vectors[i][q_l].getSubVector(0,N).dotProduct(RB_solution) );
			}
			RB_output_error_bounds_all_k[i][_k] = compute_output_dual_norm(i) * error_bound_all_k[_k];
		}

		// Initialize a vector to store the solution from the old time-step
		RealVector RB_u_old = new ArrayRealVector(N);

		// Initialize a vector to store the Newton increment, RB_delta_u
		RealVector RB_delta_u = new ArrayRealVector(N);

		// Pre-compute eval_theta_c()
		double cached_theta_c = eval_theta_c();

		for(int time_level=1; time_level<=_K; time_level++) {
			
			set_time_level(time_level); // update the member variable _k

			// Set RB_u_old to be the result of the previous Newton loop
			RB_u_old = RB_u_bar.copy();

			// Now we begin the nonlinear loop
			for(int l=0; l<n_newton_steps; ++l) {
				// Get u_euler_theta = euler_theta*RB_u_bar + (1-euler_theta)*RB_u_old
				RealVector RB_u_euler_theta =
					RB_u_bar.mapMultiply(get_euler_theta()).add( RB_u_old.mapMultiply(1.-get_euler_theta()) );

				// Assemble the left-hand side for the RB linear system
				RealMatrix RB_LHS_matrix = Base_RB_LHS_matrix.copy();

				// Add the trilinear term
				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						double new_entry = RB_LHS_matrix.getEntry(i,j);
						for(int n=0; n<N; n++) {
							new_entry += cached_theta_c*get_euler_theta()*RB_u_euler_theta.getEntry(n)*
								(RB_trilinear_form[n][i][j] + RB_trilinear_form[j][i][n]);
						}
						RB_LHS_matrix.setEntry(i,j, new_entry);
					}
				}

				// Assemble the right-hand side for the RB linear system (the residual)
				// First add forcing terms
				RealVector RB_rhs = new ArrayRealVector(N);
				
				for(int q_f=0; q_f<get_Q_f(); q_f++) {
					RB_rhs = RB_rhs.add( RB_F_q_vector[q_f].getSubVector(0, N)
										.mapMultiply( eval_theta_q_f(q_f)) );
				}

				// Now add -1./dt * M * (RB_u_bar - RB_u_old)
				RB_rhs = RB_rhs.add( RB_mass_matrix_N.operate( RB_u_bar ).mapMultiply(-1./dt) );
				RB_rhs = RB_rhs.add( RB_mass_matrix_N.operate( RB_u_old ).mapMultiply(1./dt) );

				// Now add the Aq stuff
				RB_rhs = RB_rhs.add( RB_RHS_Aq_matrix.operate(RB_u_euler_theta) );
				
				// Finally add the trilinear term
				for(int i=0; i<N; i++) {
					double new_entry = RB_rhs.getEntry(i);
					
					for(int j=0; j<N; j++) {
						double RB_u_euler_theta_j = RB_u_euler_theta.getEntry(j);
						
						for(int n=0; n<N; n++) {
							new_entry -= cached_theta_c*RB_u_euler_theta.getEntry(n)
								*RB_u_euler_theta_j*RB_trilinear_form[n][i][j];
						}
					}
					RB_rhs.setEntry(i, new_entry);
				}

				DecompositionSolver solver = new LUDecompositionImpl(RB_LHS_matrix).getSolver();
				RB_delta_u = solver.solve(RB_rhs);

				// update the Newton iterate
				RB_u_bar = RB_u_bar.add( RB_delta_u );

				// Compute the l2 norm of RB_delta_u
				double RB_delta_u_norm = RB_delta_u.getNorm();

				if( RB_delta_u_norm < nonlinear_tolerance) {
					break;
				}

				if( (l==(n_newton_steps-1)) && (RB_delta_u_norm > nonlinear_tolerance) )
				{
					throw new RuntimeException("RB Newton loop did not converge");
				}
			}

			// Load RB_solution into RB_solution_vector for residual computation
			RB_solution = RB_u_bar;
			old_RB_solution = RB_u_old;
			RB_temporal_solution_data[_k] = RB_u_bar; //should use copy here!

			double rho_LB = (mRbScmSystem == null) ? get_nominal_rho_LB() : get_SCM_lower_bound();

			// Evaluate the dual norm of the residual for RB_solution_vector
			double epsilon_N = compute_residual_dual_norm(N);

			error_bound_sum += residual_scaling_numer(rho_LB) * Math.pow(epsilon_N, 2.);

			// store error bound at time-level _k
			error_bound_all_k[_k] = Math.sqrt(error_bound_sum/residual_scaling_denom(rho_LB));

			// Now compute the outputs and associated error bounds
			for(int i=0; i<get_n_outputs(); i++) {
				RB_outputs_all_k[i][_k] = 0.;
				RB_output_error_bounds_all_k[i][_k] = 0.;
				for(int q_l=0; q_l<get_Q_l(i); q_l++) {
					RB_outputs_all_k[i][_k] +=
						 eval_theta_q_l(i, q_l) * 
						(RB_output_vectors[i][q_l].getSubVector(0,N).dotProduct(RB_solution) );
				}
				RB_output_error_bounds_all_k[i][_k] = compute_output_dual_norm(i) * error_bound_all_k[_k];
			}
			Log.d(DEBUG_TAG, "output = " + RB_outputs_all_k[0][_k] + ", bound="+RB_output_error_bounds_all_k[0][_k]);
		}

		// Now compute the L2 norm of the RB solution at time-level _K
		// to normalize the error bound
		// We reuse RB_rhs here
		double final_RB_L2_norm = 
			Math.sqrt( RB_mass_matrix_N.operate(RB_solution).dotProduct(RB_solution) );

		return ( return_rel_error_bound ? error_bound_all_k[_K]/final_RB_L2_norm : error_bound_all_k[_K] );
	}
	
	/**
	 * Override appropriate for quadratically nonlinear error bounds.
	 */
	@Override
	protected double residual_scaling_numer(double rho_LB) {
		  double tau_LB = (0.5*rho_LB < 0.) ? 0.5*rho_LB : 0.;

		  return get_dt() * tau_prod_k / ( 1.- tau_LB*get_dt() );
	}

	/**
	 * Override appropriate for quadratically nonlinear error bounds.
	 */
	@Override
	protected double residual_scaling_denom(double rho_LB) {
		  double tau_LB = (0.5*rho_LB < 0.) ? 0.5*rho_LB : 0.;

		  // Update tau_prod_k
		  tau_prod_k *= (1.+tau_LB*get_dt())/(1.-tau_LB*get_dt());

		  // and return it
		  return tau_prod_k;
	}

	/**
	 * Set the nonlinear tolerance for Newton's
	 * method for both the truth and RB solves.
	 */
	public void set_nonlinear_tolerance(double nonlinear_tolerance_in) {
		nonlinear_tolerance = nonlinear_tolerance_in;
	}
	
	/**
	 * Get the nonlinear tolerance for Newton's method.
	 */
	public double get_nonlinear_tolerance() {
		return nonlinear_tolerance;
	}

	/**
	 * Set the maximum number of Newton steps
	 * for both the truth and RB solves.
	 */
	public void set_n_newton_steps(int n_newton_steps_in) {
		n_newton_steps = n_newton_steps_in;
	}
	
	/**
	 * Set the maximum number of Newton steps
	 * for both the truth and RB solves.
	 */
	public int get_n_newton_steps() {
		return n_newton_steps;
	}

	/**
	 * Evaluate theta_c (for the quadratic nonlinearity) at the current parameter.
	 */
	public double eval_theta_c() {
		
		Method meth;

		try {
			Class partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = mAffineFnsClass.getMethod("evaluateC", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateC failed", nsme);
		}

		Double theta_val;
		try {
            Object arglist[] = new Object[1];
            arglist[0] = current_parameters.getArray();

			Object theta_obj = meth.invoke(mTheta, arglist);
			theta_val = (Double) theta_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector.
	 */
	@Override
	protected double compute_residual_dual_norm(int N) {
		// Use the stored representor inner product values
		// to evaluate the residual dual norm
		double dt = get_dt();

		double residual_norm_sq = 0.;

		// Use TransientRBSystem to compute all the linear terms
		residual_norm_sq += Math.pow( super.uncached_compute_residual_dual_norm(N), 2.);

		// Now just need to add the terms involving the nonlinearity
		RealVector RB_u_euler_theta =
			RB_solution.mapMultiply(get_euler_theta()).add( old_RB_solution.mapMultiply(1.-get_euler_theta()) );
		RealVector mass_coeffs =
			RB_solution.subtract( old_RB_solution ).mapMultiply(-1./dt);

		// Pre-compute eval_theta_c()
		double cached_theta_c = eval_theta_c();

		// All residual terms can be treated as positive quantities...
		for(int q_f=0; q_f<get_Q_f(); q_f++) {
			double cached_theta_q_f = eval_theta_q_f(q_f);
			for(int n1=0; n1<N; n1++) {
				for(int j1=0; j1<N; j1++) {
					residual_norm_sq +=
						2.*cached_theta_q_f*cached_theta_c
						*RB_u_euler_theta.getEntry(n1)*RB_u_euler_theta.getEntry(j1)
						*Fq_C_representor_norms[q_f][n1][j1];
				}
			}
		}

		for(int q_m=0; q_m<get_Q_m(); q_m++) {
			double cached_theta_q_m = eval_theta_q_m(q_m);
			for(int i=0; i<N; i++) {
				for(int n1=0; n1<N; n1++) {
					for(int j1=0; j1<N; j1++) {
						residual_norm_sq += 2.*cached_theta_q_m*cached_theta_c*
						mass_coeffs.getEntry(i)*RB_u_euler_theta.getEntry(n1)*RB_u_euler_theta.getEntry(j1)*
						Mq_C_representor_norms[q_m][i][n1][j1];
					}
				}
			}
		}

		for(int q_a=0; q_a<get_Q_a(); q_a++) {
			double cached_theta_q_a = eval_theta_q_a(q_a);
			for(int i=0; i<N; i++) {
				for(int n1=0; n1<N; n1++) {
					for(int j1=0; j1<N; j1++) {
						residual_norm_sq += 2.*cached_theta_q_a*cached_theta_c*
							RB_u_euler_theta.getEntry(i)*RB_u_euler_theta.getEntry(n1)*RB_u_euler_theta.getEntry(j1)*
							Aq_C_representor_norms[q_a][i][n1][j1];
					}
				}
			}
		}

		for(int n1=0; n1<N; n1++) {
			for(int j1=0; j1<N; j1++) {
				double RB_u_euler_theta_1 = RB_u_euler_theta.getEntry(n1)*RB_u_euler_theta.getEntry(j1);

				for(int n2=n1; n2<N; n2++) {
					int init_j2_index = (n2 == n1) ? j1 : 0;
					for(int j2=init_j2_index; j2<N; j2++) {
						double RB_u_euler_theta_2 = RB_u_euler_theta.getEntry(n2)*RB_u_euler_theta.getEntry(j2);

						double delta = ( (n2 == n1) && (j2 == j1) ) ? 1. : 2.;

						residual_norm_sq += delta*
							cached_theta_c*cached_theta_c*
							RB_u_euler_theta_1*RB_u_euler_theta_2*
							C_C_representor_norms[n1][j1][n2][j2];
					}
				}

			}
		}

		if(residual_norm_sq < 0)
		{
			// Sometimes this is negative due to rounding error,
			// but this error shouldn't affect the error bound
			// too much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt( residual_norm_sq );
	}

	// Get/set the nominal rho min/max values, we typically read these
	// in from the input file
	public void set_nominal_rho_min(double nominal_rho_LB_in) {
		nominal_rho_min = nominal_rho_LB_in;
	}
	
	public void set_nominal_rho_max(double nominal_rho_LB_in) {
		nominal_rho_max = nominal_rho_LB_in;
	}
	
	public double get_nominal_rho_min() {
		return nominal_rho_min;
	}
	
	public double get_nominal_rho_max() {
		return nominal_rho_max;
	}
	
	/**
	 * Get the nominal stability factor lower bound.
	 * By default, this is a linear function of parameter 0.
	 */
	public double get_nominal_rho_LB() {
		double mu_min = getParameterMin(0);
		double mu_max = getParameterMax(0);
		double current_mu = getCurrentParameters().getEntry(0);
		return ( nominal_rho_max*(mu_max - current_mu) +
		         nominal_rho_min*(current_mu - mu_min) ) /
		       ( mu_max - mu_min );
	}
	
	/**
	 * @return the SCM lower bound for current_parameters
	 */
	public double get_SCM_lower_bound() {

		if( mRbScmSystem != null) {
			// Cast to a QNTransientSCMSystem
			QNTransientSCMSystem qnScmSystem = (QNTransientSCMSystem)mRbScmSystem;
			
			// Tell the SCM system the number of basis functions
			qnScmSystem.set_n_basis_functions(_N_in_RB_solve);

			// Create a parameter vector in which the current time-level
			// is appended to current_parameters.
			try {
				Parameter params = new Parameter(get_n_params()+1);
				for(int i=0; i<get_n_params(); i++)
					params.setEntry(i, getCurrentParameters().getEntry(i) );
				params.setEntry(params.getNEntries()-1, _k);
				
				// Set the parameter
				qnScmSystem.setCurrentParameters( params );
			}
			catch(InconsistentStateException ise) {
				Log.e(DEBUG_TAG, ise.getMessage());
			}

			// Also, construct a vector storing the RB coefficients
			RealVector RB_u_euler_theta =
				RB_solution.mapMultiply(get_euler_theta()).
				add( old_RB_solution.mapMultiply(1.-get_euler_theta()) );

			// Pass params and RB_u_euler_theta to the associated SCM system
			qnScmSystem.set_current_RB_coeffs( RB_u_euler_theta );

			return mRbScmSystem.get_SCM_LB();
		}
		else {
			return get_SCM_from_AffineFunction();
		}
	}
	
	/**
	 * @param parameters_filename
	 *            The name of the file to parse Parse the input file to
	 *            initialize this RBSystem.
	 */
	@Override
	public void parse_parameters_file(String parameters_filename, boolean isAssetFile)
	throws InconsistentStateException {
		super.parse_parameters_file(parameters_filename, isAssetFile);

		GetPot infile = new GetPot(context, parameters_filename, isAssetFile);

		double nonlinear_tolerance_in = infile.call("nonlinear_tolerance",1.e-8);
		set_nonlinear_tolerance(nonlinear_tolerance_in);
		
		int n_newton_steps_in = infile.call("n_newton_steps", 15);
		set_n_newton_steps(n_newton_steps_in);
		
		double nominal_rho_min_in = infile.call("nominal_rho_min", 0);
		set_nominal_rho_min(nominal_rho_min_in);
		
		double nominal_rho_max_in = infile.call("nominal_rho_max", 0);
		set_nominal_rho_max(nominal_rho_max_in);

		Log.d(DEBUG_TAG,"QNTransientRBSystem parameters from " + parameters_filename + ":");
		Log.d(DEBUG_TAG,"Tolerance for Newton's method: " + get_nonlinear_tolerance());
		Log.d(DEBUG_TAG,"Max number of Newton steps: " + get_n_newton_steps());
		Log.d(DEBUG_TAG,"Nominal rho min: " + get_nominal_rho_min());
		Log.d(DEBUG_TAG,"Nominal rho max: " + get_nominal_rho_max());
	}
	
	/**
	 * Override read_offline_data_from_files in order to
	 * read in the data for the nonlinear problem.
	 */
	@Override
	public void read_offline_data(String directory_name, boolean isAssetFile)
	throws IOException {

		super.read_offline_data(directory_name, isAssetFile);

		HttpClient client = new DefaultHttpClient();

		int buffer_size = 8192;

		int n_bfs = get_n_basis_functions();

		// Read in the trlinear form
		Log.d(DEBUG_TAG, "Starting read RB_trilinear_form.dat");
		{
			InputStreamReader isr;
			String dataString = directory_name + "/RB_trilinear_form.dat";

			if(!isAssetFile) { // Read from server
				HttpGet request = new HttpGet(dataString);
				HttpResponse response = client.execute(request);
				isr = new InputStreamReader(response.getEntity()
						.getContent());
			}
			else { // Read from assets
				isr = new InputStreamReader(
						context.getAssets().open(dataString));
			}

			BufferedReader reader = new BufferedReader(isr,buffer_size);

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			// Set the size of the inner product matrix
			RB_trilinear_form = new double[n_bfs][n_bfs][n_bfs];

			// Fill the array
			int count = 0;
			for(int i=0; i<n_bfs; i++) {
				for(int j=0; j<n_bfs; j++) {
					for(int l=0; l<n_bfs; l++) {
						RB_trilinear_form[i][j][l] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading RB_trilinear_form.dat");

		// Read in Fq_C representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Fq_C_norms.dat";

			if(!isAssetFile) { // Read from server
				HttpGet request = new HttpGet(dataString);
				HttpResponse response = client.execute(request);
				isr = new InputStreamReader(response.getEntity()
						.getContent());
			}
			else { // Read from assets
				isr = new InputStreamReader(
						context.getAssets().open(dataString));
			}

			BufferedReader reader = new BufferedReader(isr,buffer_size);

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			// Declare the array
			Fq_C_representor_norms = new double[get_Q_f()][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < get_Q_f(); q_f++)
				for (int i = 0; i < n_bfs; i++) {
					for (int j=0; j < n_bfs; j++) {
						Fq_C_representor_norms[q_f][i][j] = Double
							.parseDouble(tokens[count]);
						count++;
					}
				}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading Fq_C_norms.dat");

		// Read in M_M representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Mq_C_norms.dat";

			if(!isAssetFile) { // Read from server
				HttpGet request = new HttpGet(dataString);
				HttpResponse response = client.execute(request);
				isr = new InputStreamReader(response.getEntity()
						.getContent());
			}
			else { // Read from assets
				isr = new InputStreamReader(
						context.getAssets().open(dataString));
			}

			BufferedReader reader = new BufferedReader(isr,buffer_size);

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			// Declare the array
			Mq_C_representor_norms = new double[get_Q_m()][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int q_m=0; q_m<get_Q_m(); q_m++)
				for (int i = 0; i < n_bfs; i++)
					for (int j = 0; j < n_bfs; j++) {
						for (int k = 0; k < n_bfs; k++) {
							Mq_C_representor_norms[q_m][i][j][k] = Double
								.parseDouble(tokens[count]);
							count++;
						}
					}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading Mq_C_norms.dat");

		// Read in Aq_C representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Aq_C_norms.dat";

			if(!isAssetFile) { // Read from server
				HttpGet request = new HttpGet(dataString);
				HttpResponse response = client.execute(request);
				isr = new InputStreamReader(response.getEntity()
						.getContent());
			}
			else { // Read from assets
				isr = new InputStreamReader(
						context.getAssets().open(dataString));
			}

			BufferedReader reader = new BufferedReader(isr,buffer_size);

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			// Declare the array
			Aq_C_representor_norms = new double[get_Q_a()][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for(int q_a=0; q_a<get_Q_a(); q_a++) {
				for(int i=0; i<n_bfs; i++) {
					for(int n1=0; n1<n_bfs; n1++) {
						for(int j1=0; j1<n_bfs; j1++) {
							Aq_C_representor_norms[q_a][i][n1][j1] = Double
							.parseDouble(tokens[count]);
							count++;
						}
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading Aq_C_norms.dat");
		
		// Read in C_C representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/C_C_norms.dat";

			if(!isAssetFile) { // Read from server
				HttpGet request = new HttpGet(dataString);
				HttpResponse response = client.execute(request);
				isr = new InputStreamReader(response.getEntity()
						.getContent());
			}
			else { // Read from assets
				isr = new InputStreamReader(
						context.getAssets().open(dataString));
			}

			BufferedReader reader = new BufferedReader(isr,buffer_size);

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			// Declare the array
			C_C_representor_norms = new double[n_bfs][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for(int ii=0; ii<n_bfs; ii++) {
				for(int i=0; i<n_bfs; i++) {
					for(int n1=0; n1<n_bfs; n1++) {
						for(int j1=0; j1<n_bfs; j1++) {
							C_C_representor_norms[ii][i][n1][j1] = Double
							.parseDouble(tokens[count]);
							count++;
						}
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading C_C_norms.dat");

	}

}

