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
// basis functionality for linear parabolic
// problems.
// This class is modeled on the TransientRBSystem
// class in rbOOmit

public class TransientRBSystem extends RBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "TransientRBSystem";

	/**
	 * A secondary SCM object since we might need a lower bound
	 * for the mass matrix and the stiffness matrix.
	 */
	private RBSCMSystem mSecondRbScmSystem;
	
	/**
	 * Time step size.
	 */
	private double dt;

	/**
	 * Parameter that defines the temporal discretization:
	 * euler_theta = 0   ---> Forward Euler
	 * euler_theta = 0.5 ---> Crank-Nicolson
	 * euler_theta = 1   ---> Backward Euler
	 */
	private double euler_theta;
	
	/**
	 * boolean flag that determines whether or not we impose a temporal
	 * filter on each output
	 */
	private boolean apply_temporal_filter_flag;
	
	/**
	 * The standard deviation of the filter (in terms of time steps).
	 */
	private double filter_width;

	/**
	 * The current time-level, 0 <= _k <= _K.
	 */
	protected int _k;
	
	/**
	 * Total number of time-steps.
	 */
	protected int _K;
	
	/**
	 * The number of terms in the affine expansion of the mass matrix
	 */
	private int mQ_m;

	/**
	 * RBSystem has a member RB_solution, we also need
	 * old_RB_solution here.
	 */
	protected RealVector old_RB_solution;

	/**
	 * The error bound for the field variable at each
	 * time level.
	 */
	protected double[] error_bound_all_k;

	/**
	 * The solution coefficients at each time level
	 * from the most recent RB_solve.
	 */
	protected RealVector[] RB_temporal_solution_data;

	/**
	 * Dense RB mass matrix.
	 */
	protected RealMatrix RB_L2_matrix;
	
	/**
	 * Dense matrices for the RB computations.
	 */
	protected RealMatrix[] RB_M_q_vector;

	/**
	 * Vectors storing the residual representor inner products
	 * to be used in computing the residuals online.
	 */
	protected double[][][]   Fq_Mq_representor_norms;
	protected double[][][]   Mq_Mq_representor_norms;
	protected double[][][][] Aq_Mq_representor_norms;

	// basis vector
//	protected float[][][] Z_vector;
	
	/**
	 * Private residual caching data.
	 */
	private double     cached_Fq_term;
	private double[]   cached_Fq_Aq_vector;
	private double[][] cached_Aq_Aq_matrix;
	private double[]   cached_Fq_Mq_vector;
	private double[][] cached_Aq_Mq_matrix;
	private double[][] cached_Mq_Mq_matrix;
	
	/**
	 * Constructor.
	 */
	public TransientRBSystem(Context context) {
		super(context);
	}

	/**
	 * Set the secondary SCM system
	 */
	public void setSecondarySCM(RBSCMSystem second_scm_system) {
		mSecondRbScmSystem = second_scm_system;
	}

	/**
	 * Evaluate theta_q_m (for the q^th mass matrix term) at the current parameter.
	 */
	public double eval_theta_q_m(int q) {		
        Method meth;

        try {
            // Get a reference to get_n_M_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[2];
        	partypes[0] = Integer.TYPE;
            partypes[1] = double[].class;
        	
            meth = mAffineFnsClass.getMethod("evaluateM", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateM failed", nsme);
        }
        
        Double theta_val;
        try {
            Object arglist[] = new Object[2];
            arglist[0] = new Integer(q);
            arglist[1] = current_parameters.getArray();
        	
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
	 * @return Q_m, the number of terms in the affine expansion
	 * of the mass matrix
	 */
	public int get_Q_m() {
		return mQ_m;
	}
	
	/**
	 * Perform online solve with the N RB basis functions, for the set of
	 * parameters in current_params, where 1 <= N <= RB_size.
	 */
	@Override
	public double RB_solve(int N) {
		current_N = N;

		if (N > get_n_basis_functions()) {
			throw new RuntimeException(
					"ERROR: N cannot be larger than the number "
					+ "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException(
			"ERROR: N must be greater than 0 in RB_solve");
		}

		// First assemble the mass matrix
		RealMatrix RB_mass_matrix_N = new Array2DRowRealMatrix(N, N);
		for(int q_m=0; q_m<get_Q_m(); q_m++)
		{
			RB_mass_matrix_N = RB_mass_matrix_N.add(
					RB_M_q_vector[q_m].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply(eval_theta_q_m(q_m)) );
		}

		RealMatrix RB_LHS_matrix = new Array2DRowRealMatrix(N, N);
		RealMatrix RB_RHS_matrix = new Array2DRowRealMatrix(N, N);

		RB_LHS_matrix = RB_LHS_matrix.add( RB_mass_matrix_N.scalarMultiply(1./get_dt()) );
		RB_RHS_matrix = RB_RHS_matrix.add( RB_mass_matrix_N.scalarMultiply(1./get_dt()) );

		for(int q_a=0; q_a<get_Q_a(); q_a++) {
			RB_LHS_matrix = RB_LHS_matrix.add( 
					RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply(get_euler_theta()*eval_theta_q_a(q_a)) );
			RB_RHS_matrix = RB_RHS_matrix.add( 
					RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1)
					.scalarMultiply(-(1.-get_euler_theta())*eval_theta_q_a(q_a)) );
		}

		set_time_level(0); // Sets the member variable _k to zero

		// Zero initial condition
		RealVector RB_solution_N = new ArrayRealVector(N);
		RealVector old_RB_solution_N = new ArrayRealVector(N);

		// Initialize rhs
		RealVector RB_rhs_N = new ArrayRealVector(N);

		// Load the initial condition into RB_temporal_solution_data
		RB_temporal_solution_data[_k] = RB_solution;

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
					(RB_output_vectors[i][q_l].getSubVector(0,N).dotProduct(RB_solution_N) );
			}
			RB_output_error_bounds_all_k[i][_k] = compute_output_dual_norm(i) * error_bound_all_k[_k];
		}

		double alpha_LB = get_SCM_lower_bound();
		
		// Precompute time-invariant parts of the dual norm of the residual.
		cache_online_residual_terms(N);

		for(int time_level=1; time_level<=_K; time_level++) {
			
			set_time_level(time_level); // This updates the member variable _k
			old_RB_solution_N = RB_solution_N;

			// Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
			RB_rhs_N = RB_RHS_matrix.operate(old_RB_solution_N);

			// Add forcing terms
			for(int q_f=0; q_f<get_Q_f(); q_f++) {
				RB_rhs_N = RB_rhs_N.add( 
						RB_F_q_vector[q_f].getSubVector(0,N).mapMultiply(eval_theta_q_f(q_f)) );
			}

			// Solve the linear system
			DecompositionSolver solver = new LUDecompositionImpl(RB_LHS_matrix).getSolver();
			RB_solution_N = solver.solve(RB_rhs_N);

			// Save RB_solution for current time level
			RB_temporal_solution_data[_k] = RB_solution_N;

			// Evaluate the dual norm of the residual for RB_solution_vector
			RB_solution     = RB_solution_N;
			old_RB_solution = old_RB_solution_N;
			double epsilon_N = compute_residual_dual_norm(N);

			error_bound_sum += residual_scaling_numer(alpha_LB) * Math.pow(epsilon_N, 2.);

			// store error bound at time-level _k
			error_bound_all_k[_k] = Math.sqrt(error_bound_sum/residual_scaling_denom(alpha_LB));

			// Now compute the outputs and associated errors
			for(int i=0; i<get_n_outputs(); i++) {
				RB_outputs_all_k[i][_k] = 0.;
				RB_output_error_bounds_all_k[i][_k] = 0.;
				for(int q_l=0; q_l<get_Q_l(i); q_l++) {
					RB_outputs_all_k[i][_k] +=
						eval_theta_q_l(i, q_l) * 
						(RB_output_vectors[i][q_l].getSubVector(0,N).dotProduct(RB_solution_N) );
				}
				RB_output_error_bounds_all_k[i][_k] = compute_output_dual_norm(i) * error_bound_all_k[_k];
			}
		}

		// Now compute the L2 norm of the RB solution at time-level _K
		// to normalize the error bound
		// We reuse RB_rhs here
		RealMatrix RB_L2_matrix_N = RB_L2_matrix.getSubMatrix(0, N-1, 0, N-1);
		double final_RB_L2_norm = 
			Math.sqrt( RB_L2_matrix_N.operate(RB_solution_N).dotProduct(RB_solution_N) );
		
		if(apply_temporal_filter_flag) {
			apply_temporal_filter();
		}

		return ( return_rel_error_bound ? error_bound_all_k[_K]/final_RB_L2_norm :
			error_bound_all_k[_K] );
	}
	
	/**
	 * Apply the temporal filter to the outputs
	 */
	public void apply_temporal_filter() {
			double[][] RB_convolved_outputs_all_k = new double[get_n_outputs()][get_K()+1];
			double[][] RB_convolved_error_bounds_all_k = new double[get_n_outputs()][get_K()+1];
			
	        for (int n=0; n<get_n_outputs(); n++)
	        {
	        	double output_dual_norm = compute_output_dual_norm(n);
	
	            for (int time_level=0; time_level<=get_K(); time_level++)
	            {
	                double conv_weight_integral_sq = 0.;
	                RB_convolved_outputs_all_k[n][time_level] = 0.;
	
	                for (int k_prime=0; k_prime<=get_K(); k_prime++)
	                {
	                    double time_diff = get_dt()*(time_level-k_prime);
	                    RB_convolved_outputs_all_k[n][time_level] += get_dt() *
	                            conv_weight(time_diff) *
	                            RB_outputs_all_k[n][k_prime];
	
	                    conv_weight_integral_sq += get_dt() * Math.pow( conv_weight(time_diff), 2.);
	                }
	                
	                RB_convolved_error_bounds_all_k[n][time_level] = error_bound_all_k[get_K()] *
	                output_dual_norm*Math.sqrt(conv_weight_integral_sq);
	            }
	        }
	        
	        RB_outputs_all_k = RB_convolved_outputs_all_k;
	        RB_output_error_bounds_all_k = RB_convolved_error_bounds_all_k;
	}
	
	protected double conv_weight(double x) {
        // Specify a Gaussian with standard deviation sigma
        
        double sigma = filter_width*get_dt();

        return 1./Math.sqrt( 2.*Math.PI * sigma*sigma ) * Math.exp( -x*x/(2.*sigma*sigma) );
	}

	/**
	 * Get/set dt, the time-step size.
	 */
	public double get_dt()           { return dt; }
	public void set_dt(double dt_in) { this.dt = dt_in; }

	/**
	 * Get/set euler_theta, parameter that determines
	 * the temporal discretization.
	 * euler_theta = 0   ---> Forward Euler
	 * euler_theta = 0.5 ---> Crank-Nicolson
	 * euler_theta = 1   ---> Backward Euler
	 */
	public double get_euler_theta()              { return euler_theta; }
	public void set_euler_theta(double euler_theta_in) { this.euler_theta = euler_theta_in; }

	/**
	 * Get/set the current time-level.
	 */
	public int get_time_level()          { return _k; }
	public void set_time_level(int k_in) { this._k = k_in; }

	/**
	 * Get/set K, the total number of time-steps.
	 */
	public int get_K()          { return _K; }
	public void set_K(int K_in) { this._K = K_in; }

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
		infile.isAssetFile = isAssetFile;

		double dt_in = infile.call("dt",0.);
		set_dt(dt_in);

		int K_in = infile.call("K",0);
		set_K(K_in);

		double euler_theta_in = infile.call("euler_theta",1.);
		set_euler_theta(euler_theta_in);
		
		int apply_temporal_filter_flag_in = infile.call("apply_temporal_filter_flag", 0);
		apply_temporal_filter_flag = (apply_temporal_filter_flag_in != 0);
		
		double filter_width_in = infile.call("filter_width", 2.);
		filter_width = filter_width_in;
		
		int n_plotting_steps_in = infile.call("n_plotting_steps", get_K()+1);
		n_plotting_steps = n_plotting_steps_in;

		Log.d(DEBUG_TAG,"TransientRBSystem parameters from " + parameters_filename + ":");
		Log.d(DEBUG_TAG,"dt: " + get_dt());
		Log.d(DEBUG_TAG,"Number of time steps: " + get_K());
		Log.d(DEBUG_TAG,"euler_theta (for generalized Euler): " + get_euler_theta());
		Log.d(DEBUG_TAG,"Apply a temporal filter? " + apply_temporal_filter_flag);
		if(apply_temporal_filter_flag) {
			Log.d(DEBUG_TAG,"Temporal filter std. dev. " + filter_width);
			Log.d(DEBUG_TAG,"Number of timesteps to be plotted" + n_plotting_steps);
		}
		
	}
	
	/**
	 * Set the Q_a variable from the mTheta object.
	 */
	protected void read_in_Q_m() {
		Method meth;

		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = mAffineFnsClass.getMethod("get_n_M_functions", null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_M_functions failed", nsme);
		}

		Integer Q_m;
		try {
			Object Q_m_obj = meth.invoke(mTheta, null);
			Q_m = (Integer) Q_m_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		mQ_m = Q_m.intValue();
	}

	/**
	 * Override read_offline_data_from_files in order to
	 * read in the mass matrix and initial condition
	 * data as well.
	 */
	@Override
	public void read_offline_data(String directory_name, boolean isAssetFile)
		throws IOException {

		super.read_offline_data(directory_name, isAssetFile);
		
		// Initialize the residual caching data storage
		cached_Fq_Aq_vector  = new double[get_n_basis_functions()];
		cached_Aq_Aq_matrix  = new double[get_n_basis_functions()][get_n_basis_functions()];
		cached_Fq_Mq_vector  = new double[get_n_basis_functions()];
		cached_Aq_Mq_matrix  = new double[get_n_basis_functions()][get_n_basis_functions()];
		cached_Mq_Mq_matrix  = new double[get_n_basis_functions()][get_n_basis_functions()];

		HttpClient client = new DefaultHttpClient();

		int buffer_size = 8192;

		// Read in the L2 matrix
		{
			InputStreamReader isr;
			String dataString = directory_name + "/RB_L2_matrix.dat";
			
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
			RB_L2_matrix = new Array2DRowRealMatrix(get_n_basis_functions(),
					get_n_basis_functions());

			// Fill the matrix
			int count = 0;
			for (int i = 0; i < get_n_basis_functions(); i++)
				for (int j = 0; j < get_n_basis_functions(); j++) {
					RB_L2_matrix.setEntry(i, j, Double
							.parseDouble(tokens[count]));
					count++;
				}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading RB_L2_matrix.dat");
		
		// Read in the M_q matrices
		RB_M_q_vector = new RealMatrix[get_Q_m()];
		for (int q_m = 0; q_m < get_Q_m(); q_m++) {
			InputStreamReader isr;
			String dataString = directory_name + "/RB_M_"
				+ String.format("%03d", q_m) + ".dat";
			
			if(!isAssetFile) {
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
			RB_M_q_vector[q_m] = new Array2DRowRealMatrix(get_n_basis_functions(),
					                                      get_n_basis_functions());

			// Fill the vector
			int count = 0;
			for (int i = 0; i < get_n_basis_functions(); i++)
				for (int j = 0; j < get_n_basis_functions(); j++) {
					RB_M_q_vector[q_m].setEntry(i, j, Double
							.parseDouble(tokens[count]));
					count++;
				}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading RB_M_q data");

		// Read in Fq_Mq representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Fq_Mq_norms.dat";
			
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
			Fq_Mq_representor_norms = new double[get_Q_f()][get_Q_m()][get_n_basis_functions()];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < get_Q_f(); q_f++)
				for (int q_m = 0; q_m < get_Q_m(); q_m++)
					for (int i = 0; i < get_n_basis_functions(); i++) {
						Fq_Mq_representor_norms[q_f][q_m][i] = Double
							.parseDouble(tokens[count]);
						count++;
					}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading Fq_Mq_norms.dat");

		// Read in M_M representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Mq_Mq_norms.dat";
			
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
		    int Q_m_hat = get_Q_m()*(get_Q_m()+1)/2;
			Mq_Mq_representor_norms = new double[Q_m_hat][get_n_basis_functions()][get_n_basis_functions()];

			// Fill it
			int count = 0;
			for(int q=0; q<Q_m_hat; q++)
				for (int i = 0; i < get_n_basis_functions(); i++)
					for (int j = 0; j < get_n_basis_functions(); j++) {
						Mq_Mq_representor_norms[q][i][j] = Double
							.parseDouble(tokens[count]);
						count++;
					}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading Mq_Mq_norms.dat");

		// Read in Aq_M representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Aq_Mq_norms.dat";
			
			try{
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
			Aq_Mq_representor_norms = new double[get_Q_a()][get_Q_m()][get_n_basis_functions()][get_n_basis_functions()];

			// Fill it
			int count = 0;
			for (int q_a = 0; q_a < get_Q_a(); q_a++)
				for (int q_m = 0; q_m < get_Q_m(); q_m++)
					for (int i = 0; i < get_n_basis_functions(); i++)
						for (int j = 0; j < get_n_basis_functions(); j++) {
							Aq_Mq_representor_norms[q_a][q_m][i][j] = Double
								.parseDouble(tokens[count]);
							count++;
						}
			reader.close();
			} catch (IOException iae) {
				// Declare the array
				Aq_Mq_representor_norms = new double[get_Q_a()][get_Q_m()][get_n_basis_functions()][get_n_basis_functions()];
				
				int count = 0;
				for (int i = 0; i < get_Q_a(); i++)
					for (int j = 0; j < get_Q_m(); j++){					
						BinaryReader dis;
						String bindataString = directory_name + "/Aq_Mq_" 
						   + String.format("%03d", i) + "_" + String.format("%03d", j) + "_norms.bin";					
						if(!isAssetFile) {
							HttpGet request = new HttpGet(bindataString);
							HttpResponse response = client.execute(request);
							dis = new BinaryReader(response.getEntity()
									.getContent());
						}
						else { // Read from assets
							dis = new BinaryReader(
									context.getAssets().open(bindataString));
						}
						
						for (int k = 0; k < get_n_basis_functions(); k++)
							for (int l = 0; l < get_n_basis_functions(); l++)
								Aq_Mq_representor_norms[i][j][k][l] = dis.ReadDouble();
						count++;
						dis.close();
					}
			}
		}
		
	
		Log.d(DEBUG_TAG, "Finished reading Aq_Mq_norms.dat");

// This is read in in RBSystem, no?
//		// Read calN number
//		if (get_mfield() > 0)
//		{
//			InputStreamReader isr;
//			String dataString = directory_name + "/calN.dat";
//			
//			if(!isAssetFile) {
//				HttpGet request = new HttpGet(dataString);
//				HttpResponse response = client.execute(request);
//				isr = new InputStreamReader(response.getEntity()
//						.getContent());
//			}
//			else { // Read from assets
//				isr = new InputStreamReader(
//						context.getAssets().open(dataString));
//			}
//			BufferedReader reader = new BufferedReader(isr,buffer_size);
//
//			String line = reader.readLine();
//						
//			set_calN(Integer.parseInt(line));
//			reader.close();
//		}
//		
//		Log.d(DEBUG_TAG, "Finished reading calN.dat");
		/*
		// Read in Z data
		if (get_mfield() > 0){
			Z_vector = new float[get_mfield()][get_n_basis_functions()][get_calN()];
			for (int imf = 0; imf < get_mfield(); imf++)
				for (int inbfs = 0; inbfs < get_n_basis_functions(); inbfs++){					
					BinaryReader dis;
					String dataString = directory_name + "/Z_"
					+ String.format("%03d", imf) + "_" + String.format("%03d", inbfs) + ".bin";
					
					if(!isAssetFile) {
						HttpGet request = new HttpGet(dataString);
						HttpResponse response = client.execute(request);
						dis = new BinaryReader(response.getEntity()
								.getContent());
					}
					else { // Read from assets
						dis = new BinaryReader(
								context.getAssets().open(dataString));
					}				
								
					//Z_vector[imf][inbfs] = dis.ReadFloat(get_calN());
					for (int i = 0; i < get_calN(); i++)
						Z_vector[imf][inbfs][i] = dis.ReadFloat();
								
					dis.close();
				}
		}
				
		Log.d(DEBUG_TAG, "Finished reading Z.dat");
*/
	}

	// PROTECTED FUNCTIONS

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector. This assumes that the time-independent
	 * quantities were cached in RB_solve.
	 */
	@Override
	protected double compute_residual_dual_norm(int N) {
		// This assembly assumes we have already called cache_online_residual_terms
		// and that the RB_solve parameter is constant in time

		RealVector RB_u_euler_theta = RB_solution.mapMultiply(get_euler_theta()).add(
				old_RB_solution.mapMultiply(1.-get_euler_theta()) );
		RealVector mass_coeffs = RB_solution.subtract(old_RB_solution).mapMultiply(-1./get_dt());

		double residual_norm_sq = cached_Fq_term;

		for(int i=0; i<N; i++) {
			residual_norm_sq += RB_u_euler_theta.getEntry(i)*cached_Fq_Aq_vector[i];
			residual_norm_sq += mass_coeffs.getEntry(i)*cached_Fq_Mq_vector[i];
		}

		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
			{
				residual_norm_sq += RB_u_euler_theta.getEntry(i)*RB_u_euler_theta.getEntry(j)*cached_Aq_Aq_matrix[i][j];
				residual_norm_sq += mass_coeffs.getEntry(i)*mass_coeffs.getEntry(j)*cached_Mq_Mq_matrix[i][j];
				residual_norm_sq += RB_u_euler_theta.getEntry(i)*mass_coeffs.getEntry(j)*cached_Aq_Mq_matrix[i][j];
			}


		if(residual_norm_sq < 0)
		{
			Log.d(DEBUG_TAG, "Warning: Square of residual norm is negative " +
			"in TransientRBSystem::compute_residual_dual_norm()");

			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect result
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt( residual_norm_sq );
	}
	
	/**
	 * Helper function that caches the time-independent residual quantities.
	 */
	protected void cache_online_residual_terms(int N) {

		cached_Fq_term = 0.;
		int q=0;
		for(int q_f1=0; q_f1<get_Q_f(); q_f1++) {

			double cached_theta_q_f1 = eval_theta_q_f(q_f1);
			for(int q_f2=q_f1; q_f2<get_Q_f(); q_f2++) {
				double delta = (q_f1==q_f2) ? 1. : 2.;
				cached_Fq_term += delta*cached_theta_q_f1*eval_theta_q_f(q_f2) * Fq_representor_norms[q];

				q++;
			}
		}

		for(int q_f=0; q_f<get_Q_f(); q_f++) {
			double cached_theta_q_f = eval_theta_q_f(q_f);
			for(int q_a=0; q_a<get_Q_a(); q_a++) {
				double cached_theta_q_a = eval_theta_q_a(q_a);
				for(int i=0; i<N; i++) {
					// Clear the entries on the first pass
					if( (q_f==0) && (q_a==0) )
						cached_Fq_Aq_vector[i] = 0.;

					cached_Fq_Aq_vector[i] += 2.*cached_theta_q_f*cached_theta_q_a*
					Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q=0;
		for(int q_a1=0; q_a1<get_Q_a(); q_a1++) {
			double cached_theta_q_a1 = eval_theta_q_a(q_a1);
			for(int q_a2=q_a1; q_a2<get_Q_a(); q_a2++) {
				double cached_theta_q_a2 = eval_theta_q_a(q_a2);
				double delta = (q_a1==q_a2) ? 1. : 2.;

				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						// Clear the entries on the first pass
						if(q==0)
							cached_Aq_Aq_matrix[i][j] = 0.;

						cached_Aq_Aq_matrix[i][j] += delta*
						cached_theta_q_a1*cached_theta_q_a2*
						Aq_Aq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

		for(int q_f=0; q_f<get_Q_f(); q_f++) {
			double cached_theta_q_f = eval_theta_q_f(q_f);

	    	for(int q_m=0; q_m<get_Q_m(); q_m++) {
	    		double cached_theta_q_m = eval_theta_q_m(q_m);
		
	    		for(int i=0; i<N; i++) {
	    			// Clear the entries on the first pass
					if((q_f==0) && (q_m==0))
						cached_Fq_Mq_vector[i] = 0.;
	
					cached_Fq_Mq_vector[i] += 2.*cached_theta_q_f * cached_theta_q_m * 
						Fq_Mq_representor_norms[q_f][q_m][i];
	    		}
	    	}
		}

		for(int q_a=0; q_a<get_Q_a(); q_a++) {
			double cached_theta_q_a = eval_theta_q_a(q_a);
			
	    	for(int q_m=0; q_m<get_Q_m(); q_m++) {
	    		double cached_theta_q_m = eval_theta_q_m(q_m);
	    		
				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						// Clear the entries on the first pass
						if( (q_a==0) && (q_m == 0) )
							cached_Aq_Mq_matrix[i][j] = 0.;
	
						cached_Aq_Mq_matrix[i][j] += 2.*cached_theta_q_a*cached_theta_q_m*
							Aq_Mq_representor_norms[q_a][q_m][i][j];
					}
				}
	    	}
		}
		
		q=0;
		for(int q_m1=0; q_m1<get_Q_m(); q_m1++) {
			double cached_theta_q_m1 = eval_theta_q_m(q_m1);
			for(int q_m2=q_m1; q_m2<get_Q_m(); q_m2++) {
				double cached_theta_q_m2 = eval_theta_q_m(q_m2);
				double delta = (q_m1==q_m2) ? 1. : 2.;

				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						if( q==0 )
							cached_Mq_Mq_matrix[i][j] = 0.;

						cached_Mq_Mq_matrix[i][j] += delta*
							cached_theta_q_m1*cached_theta_q_m2*
						Mq_Mq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

	}
	
	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector. This does not assume cached data hence
	 * works for parameters that change as a function of time.
	 */
	protected double uncached_compute_residual_dual_norm(int N) {

		RealVector RB_u_euler_theta = RB_solution.mapMultiply(get_euler_theta()).add(
				old_RB_solution.mapMultiply(1.-get_euler_theta()) );
		RealVector mass_coeffs = RB_solution.subtract(old_RB_solution).mapMultiply(-1./get_dt());

		double residual_norm_sq = 0.;

		int q=0;
		for(int q_f1=0; q_f1<get_Q_f(); q_f1++) {
			double cached_theta_q_f1 = eval_theta_q_f(q_f1);
			for(int q_f2=q_f1; q_f2<get_Q_f(); q_f2++) {
				double delta = (q_f1==q_f2) ? 1. : 2.;
				residual_norm_sq += delta*cached_theta_q_f1*eval_theta_q_f(q_f2) * Fq_representor_norms[q];

				q++;
			}
		}

		for(int q_f=0; q_f<get_Q_f(); q_f++) {
			double cached_theta_q_f = eval_theta_q_f(q_f);
			for(int q_a=0; q_a<get_Q_a(); q_a++) {
				double cached_theta_q_a = eval_theta_q_a(q_a);
				for(int i=0; i<N; i++) {
					residual_norm_sq += 2.*RB_u_euler_theta.getEntry(i)*cached_theta_q_f*cached_theta_q_a*
					Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q=0;
		for(int q_a1=0; q_a1<get_Q_a(); q_a1++) {
			double cached_theta_q_a1 = eval_theta_q_a(q_a1);
			for(int q_a2=q_a1; q_a2<get_Q_a(); q_a2++) {
				double cached_theta_q_a2 = eval_theta_q_a(q_a2);
				double delta = (q_a1==q_a2) ? 1. : 2.;

				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						residual_norm_sq += delta*RB_u_euler_theta.getEntry(i)*RB_u_euler_theta.getEntry(j)*
						cached_theta_q_a1*cached_theta_q_a2*
						Aq_Aq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

		// Now add the terms due to the time-derivative
		q=0;
		  for(int q_m1=0; q_m1<get_Q_m(); q_m1++)
		  {
		    double cached_theta_q_m1 = eval_theta_q_m(q_m1);
		    for(int q_m2=q_m1; q_m2<get_Q_m(); q_m2++)
		    {
		      double cached_theta_q_m2 = eval_theta_q_m(q_m2);
		      double delta = (q_m1==q_m2) ? 1. : 2.;

		      for(int i=0; i<N; i++)
		      {
		        for(int j=0; j<N; j++)
		        {
		          residual_norm_sq += delta*mass_coeffs.getEntry(i)*mass_coeffs.getEntry(j)*
		                              cached_theta_q_m1*cached_theta_q_m2*
		                              Mq_Mq_representor_norms[q][i][j];
		        }
		      }
		      q++;
		    }
		  }

		for(int q_f=0; q_f<get_Q_f(); q_f++) {
			double cached_theta_q_f = eval_theta_q_f(q_f);
			
			for(int q_m=0; q_m<get_Q_m(); q_m++) {
				double cached_theta_q_m = eval_theta_q_m(q_m);
				
				for(int i=0; i<N; i++) {
					residual_norm_sq += 2.*mass_coeffs.getEntry(i)*cached_theta_q_f*cached_theta_q_m
					* Fq_Mq_representor_norms[q_f][q_m][i];
				}
			}
		}

		for(int q_a=0; q_a<get_Q_a(); q_a++) {
			double cached_theta_q_a = eval_theta_q_a(q_a);
			for(int q_m=0; q_m<get_Q_m(); q_m++) {
				double cached_theta_q_m = eval_theta_q_m(q_m);
				
				for(int i=0; i<N; i++) {
					for(int j=0; j<N; j++) {
						residual_norm_sq += 2.*RB_u_euler_theta.getEntry(i)*mass_coeffs.getEntry(j)
						*cached_theta_q_a*cached_theta_q_m*Aq_Mq_representor_norms[q_a][q_m][i][j];
					}
				}
			}
		}

		if(residual_norm_sq < 0)
		{
			Log.d(DEBUG_TAG, "Warning: Square of residual norm is negative " +
					"in TransientRBSystem::compute_residual_dual_norm()");

			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect result
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt( residual_norm_sq );
	}

	/**
	 * Specifies the residual scaling on the numerator to
	 * be used in the a posteriori error bound. Overload
	 * in subclass in order to obtain the desired error bound.
	 */
	protected double residual_scaling_numer(double alpha_LB) {
		return get_dt();
	}

	/**
	 * Specifies the residual scaling on the denominator to
	 * be used in the a posteriori error bound. Overload
	 * in subclass in order to obtain the desired error bound.
	 */
	@Override
	protected double residual_scaling_denom(double alpha_LB) {
		return alpha_LB;
	}
	
	/**
	 * Resize the output vectors according to n_outputs.
	 */
	@Override
	protected void initialize_data_vectors() {
		super.initialize_data_vectors();
		
		RB_outputs_all_k = new double[get_n_outputs()][get_K()+1];
		RB_output_error_bounds_all_k = new double[get_n_outputs()][get_K()+1];
		
		// Resize the error bound vector
		error_bound_all_k = new double[get_K()+1];
		
		// Resize the array that stores the solution data at all time levels
		RB_temporal_solution_data = new RealVector[get_K()+1];
	}

	
	// return truth solution 
	public float[][][] get_truth_sol(){
		int N = RB_temporal_solution_data[1].getDimension();
		int nt = get_nt();
		float[][][] truth_sol = new float[get_mfield()][1][get_calN()*nt];		
		for (int ifn = 0; ifn < get_mfield(); ifn++){
			double tmpval;
			for (_k = 1; _k <= nt; _k++)
				for (int i = 0; i < get_calN(); i++){
					tmpval = 0;
					for (int j = 0; j < N; j++)
						tmpval += Z_vector[ifn][j][i]*RB_temporal_solution_data[(int)Math.round(Math.floor(_k*_K/nt))].getEntry(j);
						//tmpval += Z_vector[ifn][j][i]*RB_temporal_solution_data[_k].getEntry(j);
					truth_sol[ifn][0][(_k-1)*get_calN()+i] = (float) tmpval;
				}
		}
		return truth_sol;
	}
	/*
	public double[] get_RBsolution(int nt){
		int N = get_N();
		double[] tmpsol = new double[N*nt];
		for (_k = 1; _k <= nt; _k++)
			for (int j = 0; j < N; j++)
				tmpsol[(_k-1)*N+j] = RB_temporal_solution_data[(int)Math.round(Math.floor(_k*_K/nt))].getEntry(j);
		return tmpsol;
	}
	*/
	public int get_nt(){
		/*
		int nt = Math.round(50000/get_calN());
		nt = nt>_K?_K:nt;
		return nt;
		*/
		int nt = (int) Math.round(75000/get_calN()/(1+0.4*(get_mfield()-1))); // can go up to 150000
		nt = nt>25?25:nt; // cap nt at 25
		return nt>_K?_K:nt;
	}
}

