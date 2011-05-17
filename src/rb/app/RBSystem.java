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

import org.apache.commons.math.complex.Complex;
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

// This class provides the Online stage for
// the reduced basis method for elliptic
// steady state problems.
// This class is modeled on RBSystem from rbOOmit

public class RBSystem extends RBBase {

	// Logging tag
	private static final String DEBUG_TAG = "RBSystem";

	// PUBLIC MEMBER VARIABLES
	
	/**
	 * A reference to the SCM system.
	 */
	public RBSCMSystem mRbScmSystem;
	
	/**
	 * The RB solution vector. Stored as a Vector so that we can easily resize
	 * it during an RB_solve.
	 */
	protected RealVector RB_solution;

	/**
	 * The vector storing the RB output values in the steady-state case.
	 */
	public double[] RB_outputs;

	/**
	 * The vector storing the RB error bounds in the steady-state case.
	 */
	public double[] RB_output_error_bounds;

	/**
	 * Boolean flag to indicate whether RB_solve returns an absolute or relative
	 * error bound. True => relative, false => absolute.
	 */
	public boolean return_rel_error_bound;

	// PRIVATE MEMBER VARIABLES

	/**
	 * The number of basis functions in the RB space.
	 */
	private int n_bfs;
	
	private int calN;
	
	// current N
	protected int current_N;
	
	/**
	 * The number of terms in the affine expansion of the rhs
	 */	
	private int mQ_f;
	
	private int mQ_uL;
	
	/**
	 * The number of output functionals
	 */
	private int mN_outputs;
	
	/**
	 * The number of time steps we actually plot in the output plotter.
	 * This is sometimes less than K so that we don't plot the effect
	 * of the filter near the final time.
	 */
	protected int n_plotting_steps;

	/**
	 * Dense matrices for the RB computations.
	 */
	protected RealMatrix[] RB_A_q_vector;

	/**
	 * Dense vector for the RHS.
	 */
	protected RealVector[] RB_F_q_vector;

	/**
	 * The vectors storing the RB output vectors
	 */
	protected RealVector[][] RB_output_vectors;

	protected float[][][] Z_vector;

	/* The number of fields: 0 means no visualization */
	private int mfield; 
	
	protected float[][] uL_vector;

	/**
	 * This array stores the dual norms for each output.
	 * Row n stores the Q_l dual norms for the expansion
	 * of the n^th output.
	 */
	public double[][] output_dual_norms;

	/**
	 * Arrays storing the residual representor inner products to be used in
	 * computing the residuals in the Online stage. These are resized by reading
	 * in the Offline data written out by rbOOmit.
	 */
	protected double[] Fq_representor_norms;
	protected double[][][] Fq_Aq_representor_norms;
	protected double[][][] Aq_Aq_representor_norms;

	// (those two below are moved from TransientRBSystem
	/**
	 * The outputs at all time levels.
	 */
	protected double[][] RB_outputs_all_k;
	/**
	 * The output error bounds at all time levels.
	 */
	protected double[][] RB_output_error_bounds_all_k;

	protected double[][][] RB_sweep_solution;
	
	// PUBLIC FUNCTIONS

	/**
	 * Constructor.
	 */
	public RBSystem(Context context) {
		super(context);
		
		// Initialize n_bfs to 0
		n_bfs = 0;
	}
	
	/**
	 * Set the primary SCM
	 */
	public void setPrimarySCM(RBSCMSystem scm_system) {
		mRbScmSystem = scm_system;
	}
	
	/**
	 * Static builder function.
	 */
	public static RBSystem buildRBSystem(Context context,
			RBActivityEnums.SystemTypeEnum type) {
		
		switch(type) {
		case NONE:
			return null;
		case LINEAR_STEADY:
			return new RBSystem(context);
		case LINEAR_COMPLEX_STEADY:
			return new ComplexRBSystem(context);			
		case LINEAR_UNSTEADY:
			return new TransientRBSystem(context);
		case QN_UNSTEADY:
			return new QNTransientRBSystem(context);
		default:
			return null;
		}
		
	}

	/**
	 * Perform online solve with the N RB basis functions, for the set of
	 * parameters in current_params, where 1 <= N <= RB_size.
	 */
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

		// Assemble the RB system
		RealMatrix RB_system_matrix_N = new Array2DRowRealMatrix(N, N);

		for (int q_a = 0; q_a < get_Q_a(); q_a++) {
			RB_system_matrix_N = RB_system_matrix_N.add (
				RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1).
				scalarMultiply(eval_theta_q_a(q_a) ) );
		}

		// Assemble the RB rhs
		RealVector RB_rhs_N = new ArrayRealVector(N);

		for (int q_f = 0; q_f < get_Q_f(); q_f++) {
			// Note getSubVector takes an initial index and the number of entries
			// i.e. the interface is a bit different to getSubMatrix
			RB_rhs_N = RB_rhs_N.add(RB_F_q_vector[q_f].getSubVector(0, N)
					.mapMultiply(eval_theta_q_f(q_f)));
		}

		// Solve the linear system
		DecompositionSolver solver = new LUDecompositionImpl(RB_system_matrix_N).getSolver();
		RB_solution = solver.solve(RB_rhs_N);

		// Evaluate the dual norm of the residual for RB_solution_vector
		double epsilon_N = compute_residual_dual_norm(N);

		// Get lower bound for coercivity constant
		double alpha_LB = get_SCM_lower_bound();
		
		// If SCM lower bound is negative
		if (alpha_LB<=0){ // Get an upper bound instead
			alpha_LB = get_SCM_upper_bound();
		}

		// Store (absolute) error bound
		double abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

		// Compute the norm of RB_solution
		double RB_solution_norm = RB_solution.getNorm();

		// Now compute the outputs and associated errors
		RealVector RB_output_vector_N = new ArrayRealVector(N);
		for (int i = 0; i < get_n_outputs(); i++) {
			RB_outputs[i] = 0.;
			
			for(int q_l=0; q_l < get_Q_l(i); q_l++) {
				RB_output_vector_N = RB_output_vectors[i][q_l].getSubVector(0, N);
				RB_outputs[i] += eval_theta_q_l(i, q_l) * RB_solution.dotProduct(RB_output_vector_N);
			}
			RB_output_error_bounds[i] = compute_output_dual_norm(i) * abs_error_bound;
		}

		return (return_rel_error_bound ? abs_error_bound / RB_solution_norm
				: abs_error_bound);
	}
	
	public double get_dt(){
		return 0.;
	}
	
	public int get_K(){
		return 1;
	}
	
	public float[][][] get_truth_sol(){
		int N = RB_solution.getDimension();
		float[][][] truth_sol = new float[get_mfield()][1][calN];
		for (int ifn = 0; ifn < get_mfield(); ifn++){
			double tmpval;
			for (int i = 0; i < calN; i++){
				tmpval = 0;
				for (int j = 0; j < N; j++)
					tmpval += Z_vector[ifn][j][i]*get_soln_coeff(j);
				truth_sol[ifn][0][i] = (float) tmpval;
			}
		}
		return truth_sol;
	}
	
	void set_sweep_sol(double[][][] _sweep_sol){
		RB_sweep_solution = _sweep_sol; 
	}
	
	public float[][][] get_sweep_truth_sol(){
		int N = RB_sweep_solution[0][0].length;
		int numSweep = RB_sweep_solution.length;
		float[][][] truth_sol = new float[get_mfield()][1][calN*numSweep];
		for (int ifn = 0; ifn < get_mfield(); ifn++){
			double tmpval;
			for (int iSweep = 0; iSweep < numSweep; iSweep++)
				for (int i = 0; i < calN; i++){
					tmpval = 0;
					for (int j = 0; j < N; j++)
						tmpval += Z_vector[ifn][j][i]*RB_sweep_solution[iSweep][0][j];
					truth_sol[ifn][0][iSweep*calN+i] = (float) tmpval;
				}
		}
		return truth_sol;
	}
		
	/**
	 * @return the SCM lower bound for current_parameters
	 */
	public double get_SCM_lower_bound() {
		if( mRbScmSystem != null) {
			mRbScmSystem.setCurrentParameters(getCurrentParameters());
			return mRbScmSystem.get_SCM_LB();
		}
		else {
			return get_SCM_from_AffineFunction();
		}
	}
	
	/**
	 * @return the SCM upper bound for current_parameters
	 */
	double get_SCM_upper_bound() {
		
		if( mRbScmSystem != null) {
			mRbScmSystem.setCurrentParameters(getCurrentParameters());
			return mRbScmSystem.get_SCM_UB();
		}
		else {
			return get_SCM_from_AffineFunction();
		}
	}
	
	/**
	 * A private helper function to get the SCM from AffineFunctions
	 * in the case that SCM_TYPE = NONE
	 */
	protected double get_SCM_from_AffineFunction() {
		// we assume that an SCM LB function has been specified
		// in AffineFunctions.jar
        Method meth;

        try {
        	Class partypes[] = new Class[1];
            partypes[0] = double[].class;
        	
            meth = mAffineFnsClass.getMethod("get_SCM_LB", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for get_SCM_LB failed", nsme);
        }
        
        Double SCM_val;
        try {
            Object arglist[] = new Object[1];
            arglist[0] = current_parameters.getArray();
        	
        	Object SCM_obj = meth.invoke(mTheta, arglist);
        	SCM_val = (Double) SCM_obj;
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
        
        return SCM_val.doubleValue();
	}	

	/**
	 * @return coefficient i of RB_solution from the most recent RB_solve.
	 */
	double get_soln_coeff(int i) {
		return RB_solution.getEntry(i);
	}

	double[][] get_RBsolution(){
		double[][] RBsol = new double[1][];
		RBsol[0] = RB_solution.toArray();
		return RBsol;
	}

	/**
	 * @return The number of basis functions in the system.
	 */
	public int get_n_basis_functions() {
		return n_bfs;
	}
	
	public void set_n_basis_functions(int _N) {
		n_bfs = _N;
	}
	
	// Return calN
	public int get_calN() {
		return calN;
	}
	
	// Set calN
	public void set_calN(int _calN) {
		calN = _calN;
	}
	
	/**
	 * @return Q_f, the number of term in the affine expansion
	 * of the right-hand side
	 */
	public int get_Q_f() {
		return mQ_f;
	}
	
	/**
	 * @return the number of output functionals
	 */
	public int get_n_outputs() {
		return mN_outputs;
	}
	
	/**
	 * Set the Q_f variable from the mTheta object.
	 */
	protected void read_in_Q_f() {
		Method meth;

		try {
			// Get a reference to get_n_F_functions, which does not
			// take any arguments
			meth = mAffineFnsClass.getMethod("get_n_F_functions", null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_F_functions failed", nsme);
		}

		Integer Q_f;
		try {
			Object Q_f_obj = meth.invoke(mTheta, null);
			Q_f = (Integer) Q_f_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		mQ_f = Q_f.intValue();
	}

	/**
	 * Set the n_outputs variable from the mTheta object.
	 */
	protected void read_in_n_outputs() {
		Method meth;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments
			meth = mAffineFnsClass.getMethod("get_n_outputs", null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_outputs failed", nsme);
		}

		Integer n_outputs;
		try {
			Object n_outputs_obj = meth.invoke(mTheta, null);
			n_outputs = (Integer) n_outputs_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		mN_outputs = n_outputs.intValue();
	}
	
	/**
	 * @param output_index The index of the output we are interested in
	 * @return the number of terms in the affine expansion of the specified output
	 */
	protected int get_Q_l(int output_index) {
		Method meth;

		try {
			// Get a reference to get_Q_l, which takes an int argument
        	Class partypes[] = new Class[1];
        	partypes[0] = Integer.TYPE;
			meth = mAffineFnsClass.getMethod("get_Q_l", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_Q_l failed", nsme);
		}

		Integer Q_l;
		try {
            Object arglist[] = new Object[1];
            arglist[0] = new Integer(output_index);
			Object Q_l_obj = meth.invoke(mTheta, arglist);
			Q_l = (Integer) Q_l_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return Q_l.intValue();
	}

	protected void read_in_Q_uL() {
		Method meth;

		boolean noQ_uLdefined = false;
		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = mAffineFnsClass.getMethod("get_n_uL_functions", null);
		} catch (NoSuchMethodException nsme) {
			//throw new RuntimeException("getMethod for get_n_uL_functions failed", nsme);
			noQ_uLdefined = true;
			meth = null;
		}

		if (noQ_uLdefined)
			mQ_uL = 0;
		else {
			Integer Q_uL;
			try {
				Object Q_uL_obj = meth.invoke(mTheta, null);
				Q_uL = (Integer) Q_uL_obj;
			}
			catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			}
			catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}
	
			mQ_uL = Q_uL.intValue();
		}
	}
	
	public int get_Q_uL(){
		return mQ_uL;
	}
	
	/**
	 * Evaluate theta_q_f (for the q^th rhs function) at the current parameter.
	 */
	public double eval_theta_q_f(int q) {
        Method meth;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[2];
        	partypes[0] = Integer.TYPE;
            partypes[1] = double[].class;
        	
            meth = mAffineFnsClass.getMethod("evaluateF", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateF failed", nsme);
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
	 * Evaluate theta_q_l (for the n^th output) at the current parameter.
	 */
	public double eval_theta_q_l(int n, int q_l) {
        Method meth;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[3];
        	partypes[0] = Integer.TYPE;
        	partypes[1] = Integer.TYPE;
            partypes[2] = double[].class;
        	
            meth = mAffineFnsClass.getMethod("evaluateL", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateL failed", nsme);
        }
        
        Double theta_val;
        try {
            Object arglist[] = new Object[3];
            arglist[0] = new Integer(n);
            arglist[1] = new Integer(q_l);
            arglist[2] = current_parameters.getArray();
        	
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
	
	public float[][] get_tranformation_data() {
        Method meth = null;
        boolean isOK = true;
        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[1];
        	partypes[0] = double[].class;
        	
            meth = mAffineFnsClass.getMethod("get_local_transformation", partypes);
        } catch (NoSuchMethodException nsme) {
            //throw new RuntimeException("getMethod for evaluateF failed", nsme);
        	isOK = false;
        }
        
        float[][] T_vector;
        if (isOK){	        
	        try {
	            Object arglist[] = new Object[1];
	            arglist[0] = current_parameters.getArray();
	        	
	        	Object theta_obj = meth.invoke(mTheta, arglist);
	        	T_vector = (float[][]) theta_obj;
	        }
	        catch (IllegalAccessException iae) {
	            throw new RuntimeException(iae);
	        }
	        catch (InvocationTargetException ite) {
	            throw new RuntimeException(ite.getCause());
	        }
        } else{
        	T_vector = new float[1][12];
        	T_vector[0][0] = 1f; T_vector[0][1] = 0f; T_vector[0][2] = 0f;
        	T_vector[0][3] = 0f; T_vector[0][4] = 1f; T_vector[0][5] = 0f;
        	T_vector[0][6] = 0f; T_vector[0][7] = 0f; T_vector[0][8] = 1f;
        	T_vector[0][9] = 0f; T_vector[0][10] = 0f; T_vector[0][11] = 0f;
        }	        
        
        return T_vector;
	}
	
	/**
	 * @param parameters_filename
	 *            The name of the file to parse Parse the input file to
	 *            initialize this RBSystem.
	 * @param isAssetFile indicates whether the file we're reading from is in
	 *        Android's asset directory
	 */
	public void parse_parameters_file(String parameters_filename, boolean isAssetFile)
	throws InconsistentStateException {

		Log.d(DEBUG_TAG, "Entered parse_parameters_file, filename = " + parameters_filename);

		GetPot infile = new GetPot(context, parameters_filename, isAssetFile);

		Log.d(DEBUG_TAG, "Created GetPot object");
		
		int n_parameters = infile.call("n_parameters",1);
		Log.d(DEBUG_TAG, "n_parameters = " + n_parameters);

		min_parameter = new Parameter(n_parameters);
		max_parameter = new Parameter(n_parameters);
		for(int i=0; i<n_parameters; i++) {
			// Read in the min/max for the i^th parameter
			String min_string = new String("mu" + i + "_min");
			double mu_i_min = infile.call(min_string, 0.);
			min_parameter.setEntry(i, mu_i_min);

			String max_string = new String("mu" + i + "_max");
			double mu_i_max = infile.call(max_string, 0.);
			max_parameter.setEntry(i, mu_i_max);
		}

		Log.d(DEBUG_TAG,"RBBase parameters from " + parameters_filename + ":");
		for(int i=0; i<n_parameters; i++) {
			Log.d(DEBUG_TAG,"Parameter " + i +
					": Min = " + getParameterMin(i) +
					", Max = " + getParameterMax(i) );
		}

		mfield = infile.call("n_field",1);
		Log.d(DEBUG_TAG,"n_field = " + mfield);
		

		int return_rel_error_bound_in =
			infile.call("return_rel_error_bound",1);
		return_rel_error_bound =
			(return_rel_error_bound_in != 0);

		Log.d(DEBUG_TAG,"RBSystem parameters from " + parameters_filename + ":");
		Log.d(DEBUG_TAG,"return a relative error bound from RB_solve? " +
				return_rel_error_bound);
	}

	/**
	 * @param directory_name
	 *            The URL of the directory containing the Offline data Read in
	 *            the Offline data to initialize this RBSystem.
	 */
	public void read_offline_data(String directory_name, boolean isAssetFile)
			throws IOException {

		HttpClient client = new DefaultHttpClient();
		
		
		int buffer_size = 8192;

		// Find out dimension of the RB space
		{
			InputStreamReader isr;
			String dataString = directory_name + "/n_bfs.dat";
			
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

			n_bfs = Integer.parseInt(line);

			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");

		// Read in output data
		if (get_n_outputs() > 0) {
			// Get output dual norms
			{

				RB_output_vectors = new RealVector[get_n_outputs()][];
				output_dual_norms = new double[get_n_outputs()][];
				for (int i = 0; i < get_n_outputs(); i++) {
					
					InputStreamReader isr1;
					String dataString1 = directory_name + "/output_" + String.format("%03d", i) + 
						"_dual_norms.dat";
					
					if(!isAssetFile) { // Read from server
						HttpGet request1 = new HttpGet(dataString1);
						HttpResponse response1 = client.execute(request1);
						isr1 = new InputStreamReader(response1
								.getEntity().getContent());
					}
					else { // Read from assets
						isr1 = new InputStreamReader(
								context.getAssets().open(dataString1));
					}
					BufferedReader reader1 = new BufferedReader(isr1,buffer_size);

					String line1 = reader1.readLine();
					String[] dual_norms_tokens = line1.split(" ");
					
					int Q_l_hat = get_Q_l(i) * (get_Q_l(i) + 1) / 2;
					output_dual_norms[i] = new double[Q_l_hat];
					for(int q=0; q<Q_l_hat; q++) {
						output_dual_norms[i][q] = Double
							.parseDouble(dual_norms_tokens[q]);
					}
					
					
					RB_output_vectors[i] = new RealVector[get_Q_l(i)];
					for (int q_l=0; q_l < get_Q_l(i); q_l++) {
						// Now read in the RB output vectors
						InputStreamReader isr_i;
						String dataString_i = directory_name + "/output_"
							+ String.format("%03d", i) + "_" + String.format("%03d", q_l) + ".dat";
						if(!isAssetFile) { // Read from server
							HttpGet request_i = new HttpGet(dataString_i);
							HttpResponse response_i = client.execute(request_i);
							isr_i = new InputStreamReader(response_i
									.getEntity().getContent());
						}
						else { // Read from assets
							isr_i = new InputStreamReader(
									context.getAssets().open(dataString_i));
						}
						BufferedReader reader_i = new BufferedReader(isr_i,buffer_size);

						String line_i = reader_i.readLine();
						String[] output_i_tokens = line_i.split(" ");

						RB_output_vectors[i][q_l] = new ArrayRealVector(n_bfs);
						for (int j = 0; j < n_bfs; j++) {
							RB_output_vectors[i][q_l].setEntry(j, Double
									.parseDouble(output_i_tokens[j]));
						}
						reader_i.close();
						
					}

					reader1.close();
				}
			}
		}
		
		Log.d(DEBUG_TAG, "Finished reading output data");

		// Read in the F_q vectors
		RB_F_q_vector = new RealVector[get_Q_f()];
		for (int q_f = 0; q_f < get_Q_f(); q_f++) {
			InputStreamReader isr;
			String dataString = directory_name + "/RB_F_"
				+ String.format("%03d", q_f) + ".dat";
			
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
			RB_F_q_vector[q_f] = new ArrayRealVector(n_bfs);

			// Fill the vector
			for (int i = 0; i < n_bfs; i++) {
				RB_F_q_vector[q_f].setEntry(i, Double.parseDouble(tokens[i]));
			}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading RB_F_q data");

		// Read in the A_q matrices
		RB_A_q_vector = new RealMatrix[get_Q_a()];
		for (int q_a = 0; q_a < get_Q_a(); q_a++) {
			InputStreamReader isr;
			String dataString = directory_name + "/RB_A_"
				+ String.format("%03d", q_a) + ".dat";
			
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
			RB_A_q_vector[q_a] = new Array2DRowRealMatrix(n_bfs, n_bfs);

			// Fill the vector
			int count = 0;
			for (int i = 0; i < n_bfs; i++)
				for (int j = 0; j < n_bfs; j++) {
					RB_A_q_vector[q_a].setEntry(i, j, Double
							.parseDouble(tokens[count]));
					count++;
				}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading RB_A_q data");

		// Read in F_q representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Fq_norms.dat";
			
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

			// Declare the array
			int Q_f_hat = get_Q_f() * (get_Q_f() + 1) / 2;
			Fq_representor_norms = new double[Q_f_hat];

			// Fill it
			for (int i = 0; i < Q_f_hat; i++) {
				Fq_representor_norms[i] = Double.parseDouble(tokens[i]);
			}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading Fq_norms.dat");

		// Read in Fq_Aq representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Fq_Aq_norms.dat";
			
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

			// Declare the array
			Fq_Aq_representor_norms = new double[get_Q_f()][get_Q_a()][n_bfs];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < get_Q_f(); q_f++)
				for (int q_a = 0; q_a < get_Q_a(); q_a++)
					for (int i = 0; i < n_bfs; i++) {
						Fq_Aq_representor_norms[q_f][q_a][i] = Double
								.parseDouble(tokens[count]);
						count++;
					}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading Fq_Aq_norms.dat");
/*
		// Read in Aq_Aq representor norm data
		{
			InputStreamReader isr;
			String dataString = directory_name + "/Aq_Aq_norms.dat";
			
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

			// Declare the array
			int Q_a_hat = get_Q_a() * (get_Q_a() + 1) / 2;
			Aq_Aq_representor_norms = new double[Q_a_hat][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int i = 0; i < Q_a_hat; i++)
				for (int j = 0; j < n_bfs; j++)
					for (int l = 0; l < n_bfs; l++) {
						Aq_Aq_representor_norms[i][j][l] = Double
								.parseDouble(tokens[count]);
						count++;
					}

			reader.close();
		}
*/
		// Read in Aq_Aq representor norm data
		{
			// Declare the array
			int Q_a_hat = get_Q_a() * (get_Q_a() + 1) / 2;
			Aq_Aq_representor_norms = new double[Q_a_hat][n_bfs][n_bfs];

			int count = 0;
			for (int i = 0; i < get_Q_a(); i++)
				for (int j = i; j < get_Q_a(); j++){					
					BinaryReader dis;
					String dataString = directory_name + "/Aq_Aq_" 
					   + String.format("%03d", i) + "_" + String.format("%03d", j) + "_norms.bin";					
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
					
					for (int k = 0; k < n_bfs; k++)
						for (int l = 0; l < n_bfs; l++)
							Aq_Aq_representor_norms[count][k][l] = dis.ReadDouble();
					count++;
					dis.close();
				}
		}
		Log.d(DEBUG_TAG, "Finished reading Aq_Aq_norms.dat");
		
		// Read calN number
		if (get_mfield() > 0)
		{
			InputStreamReader isr;
			String dataString = directory_name + "/calN.dat";
			
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
						
			calN = Integer.parseInt(line);
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading calN.dat");
		
		// Reading uL data
		if (get_Q_uL()>0){
			uL_vector = new float[get_Q_uL()][get_calN()];
			for (int q_uL = 0; q_uL < get_Q_uL(); q_uL++)
			{
				BinaryReader dis;
				String dataString = directory_name + "/uL_" + String.format("%03d", q_uL) + ".bin";
				
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
				
				uL_vector[q_uL] = dis.ReadFloat(get_calN());				
				
				dis.close();
			}
		}
		Log.d(DEBUG_TAG, "Finished reading uL.dat");
		
		// Read in Z data
		if (get_mfield() > 0){
			Z_vector = new float[get_mfield()][n_bfs][calN];
			for (int imf = 0; imf < get_mfield(); imf++)
				for (int inbfs = 0; inbfs < n_bfs; inbfs++){					
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
								
					Z_vector[imf][inbfs] = dis.ReadFloat(calN);
								
					dis.close();
				}
		}
				
		Log.d(DEBUG_TAG, "Finished reading Z.dat");
		
		initialize_data_vectors();
	}

	// PROTECTED FUNCTIONS

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector.
	 */
	protected double compute_residual_dual_norm(int N) {

		// Use the stored representor inner product values
		// to evaluate the residual norm
		double residual_norm_sq = 0.;

		int q = 0;
		for (int q_f1 = 0; q_f1 < get_Q_f(); q_f1++) {
			for (int q_f2 = q_f1; q_f2 < get_Q_f(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				residual_norm_sq += delta * eval_theta_q_f(q_f1)
						* eval_theta_q_f(q_f2) * Fq_representor_norms[q];

				q++;
			}
		}

		for (int q_f = 0; q_f < get_Q_f(); q_f++) {
			for (int q_a = 0; q_a < get_Q_a(); q_a++) {
				for (int i = 0; i < N; i++) {
					double delta = 2.;
					residual_norm_sq += get_soln_coeff(i) * delta
					* eval_theta_q_f(q_f) * eval_theta_q_a(q_a)
					* Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q = 0;
		for (int q_a1 = 0; q_a1 < get_Q_a(); q_a1++) {
			for (int q_a2 = q_a1; q_a2 < get_Q_a(); q_a2++) {
				double delta = (q_a1 == q_a2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						residual_norm_sq += get_soln_coeff(i)
								* get_soln_coeff(j) * delta
								* eval_theta_q_a(q_a1) * eval_theta_q_a(q_a2)
								* Aq_Aq_representor_norms[q][i][j];
					}
				}

				q++;
			}
		}

		if (residual_norm_sq < 0.) {
			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect error bound much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}
	
	/**
	 * Compute the dual norm of the i^th output function at the current
	 * parameter value
	 */
	protected double compute_output_dual_norm(int i) {
		
		// Use the stored representor inner product values
		// to evaluate the output dual norm
		double output_norm_sq = 0.;
		
		int q = 0;
		for (int q_l1 = 0; q_l1 < get_Q_l(i); q_l1++) {
			for (int q_l2 = q_l1; q_l2 < get_Q_l(i); q_l2++) {
				double delta = (q_l1 == q_l2) ? 1. : 2.;
				output_norm_sq += delta * eval_theta_q_l(i, q_l1)
						* eval_theta_q_l(i, q_l2) * output_dual_norms[i][q];

				q++;
			}
		}
		
		return Math.sqrt(output_norm_sq);
	}

	/**
	 * Specifies the residual scaling on the denominator to be used in the a
	 * posteriori error bound. Overload in subclass in order to obtain the
	 * desired error bound.
	 */
	protected double residual_scaling_denom(double alpha_LB) {
		return Math.sqrt(alpha_LB);
	}
	
	/**
	 * Resize the vectors that store solution data and output data.
	 */
	protected void initialize_data_vectors() {
		// Also, resize RB_outputs and RB_output_error_error_bounds arrays
		RB_outputs = new double[get_n_outputs()];
		RB_output_error_bounds = new double[get_n_outputs()];
	}

	public int get_nt(){
		return 1;
	}
	
	public int get_N(){
		return current_N;
	}
	
	public int get_mfield() {
		return mfield;
	}
	
	double get_RB_output(int n_output, boolean Rpart){
		return RB_outputs[n_output];
	}
	
	double get_RB_output_error_bound(int n_output, boolean Rpart){
		return RB_output_error_bounds[n_output];
	}

}

