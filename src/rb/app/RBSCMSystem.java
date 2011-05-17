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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.linear.LinearConstraint;
import org.apache.commons.math.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math.optimization.linear.Relationship;
import org.apache.commons.math.optimization.linear.SimplexSolver;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;

import android.content.Context;
import android.util.Log;

// This class implements the Online stage
// of the Successive Constraint Method
// for coercive problems.
// This class is modeled on the RBSCMSystem
// class in rbOOmit

public class RBSCMSystem extends RBBase {

	// Logging tag
	private static final String DEBUG_TAG = "RBSCMSystem";

	/**
	 * The maximum number of constraints we impose.
	 */
	private int SCM_M;

	/**
	 * The bounding box values
	 */
	private double[] B_min;
	private double[] B_max;

	/**
	 * The Greedily selected parameters.
	 */
	private Vector<Parameter> C_J;

	/**
	 * The values of the stability factor at the
	 * greedily selected parameters.
	 */
	protected double[] C_J_stability_vector;

	/**
	 * This matrix stores the infimizing vectors
	 * y_1(\mu),...,y_Q_a(\mu), for each \mu in
	 * C_J, which are used in computing the SCM
	 * upper bounds.
	 */
	private double[][] SCM_UB_vectors;

	/**
	 * A private Parameter used to temporarily store
	 * current_parameters during the SCM calculation
	 */
	private Parameter saved_parameters;

	/**
	 * Constructor.
	 */
	public RBSCMSystem(Context context) {
		super(context);
	}
	
	/**
	 * Static builder function.
	 */
	public static RBSCMSystem buildSCMSystem(Context context, RBActivityEnums.SCMTypeEnum type) {
		
		switch(type) {
		case NONE:
			return null;
		case COERCIVE:
			return new RBSCMSystem(context);
		case QN_TRANSIENT_SCM:
			return new QNTransientSCMSystem(context);
		case COMPLEX_NONCOERCIVE:
			return new RBnSCMCSystem(context);
		default:
			return null;
		}
		
	}

	/**
	 * @return the SCM lower bound for the current parameters.
	 */
	public double get_SCM_LB()
	{
		double min_J_obj = 0.;

		try {

			// First, declare the constraints
			Collection constraints = new ArrayList();

			// Add bounding box constraints for the get_Q_a() variables
			for(int q=0; q<get_Q_a(); q++)
			{
				double[] index = new double[get_Q_a()];
				index[q] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ, B_min[q]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ, B_max[q]));
			}

			// Sort the indices of C_J based on distance from current_parameters
			List<Integer> sortedIndices = getSorted_CJ_Indices();

			// Save the current_parameters since we'll change them in the loop below
			save_current_parameters();

			// Add the constraint rows
			int n_rows = Math.min( SCM_M, C_J.size() );
			int count = 1;
			
			if(n_rows > 0) {
				for (Iterator it = sortedIndices.iterator(); it.hasNext();) {
					Integer mu_index = (Integer)it.next();
	
					get_current_parameters_from_C_J(mu_index);
					// current_parameters = C_J.get(mu_index);
	
					double[] constraint_row = new double[get_Q_a()];
					for(int q=0; q<get_Q_a(); q++)
					{
						constraint_row[q] = eval_theta_q_a(q);
					}
	
					constraints.add(new LinearConstraint(constraint_row, 
							Relationship.GEQ, C_J_stability_vector[mu_index]));
	
					if(count >= n_rows)
						break;
	
					count++;
	
				}
			}


			// Now load the original parameters back into current_parameters
			// in order to set the coefficients of the objective function
			reload_current_parameters();
			
			// Create objective function object
			double[] objectiveFn = new double[get_Q_a()];
			for(int q=0; q<get_Q_a(); q++)
			{
				objectiveFn[q] = eval_theta_q_a(q);
			}
			LinearObjectiveFunction f = new LinearObjectiveFunction(objectiveFn, 0.);
			
			SimplexSolver solver = new SimplexSolver();			
			RealPointValuePair opt_pair = solver.optimize(f, constraints, GoalType.MINIMIZE, false);
			min_J_obj = opt_pair.getValue();
		}
		catch(OptimizationException e) {
			Log.e("DEBUG_TAG", "Optimal solution not found");
			e.printStackTrace();
		}
		catch(Exception e) {
			Log.e("DEBUG_TAG", "Exception occurred during SCM_LB calculation");
			e.printStackTrace();
		}

		Log.d(DEBUG_TAG, "SCM val = " + min_J_obj);
		return min_J_obj;
	}

	/**
	 * Evaluate the SCM upper bound for current_parameters.
	 */
	public double get_SCM_UB() {

		// Sort the indices of C_J based on distance from current_parameters
		List<Integer> sortedIndices = getSorted_CJ_Indices();

		// For each mu, we just find the minimum of J_obj over
		// the subset of vectors in SCM_UB_vectors corresponding
		// to C_J_M (SCM_UB_vectors contains vectors for all of
		// C_J).
		double min_J_obj = 0.;
		int n_rows = Math.min( SCM_M, C_J.size() );
		int count = 1;
		for (Iterator it = sortedIndices.iterator(); it.hasNext();) {
			Integer mu_index = (Integer)it.next();

			current_parameters = C_J.get(mu_index);

			double[] UB_vector = SCM_UB_vectors[mu_index];

			double J_obj = 0.;
			for(int q=0; q<get_Q_a(); q++) {
				J_obj += eval_theta_q_a(q)*UB_vector[q];
			}

			if( (count==1) || (J_obj < min_J_obj) ) {
				min_J_obj = J_obj;
			}

			if(count >= n_rows)
				break;

			count++;
		}

		return min_J_obj;
	}

	/**
	 * @return Euclidean distance between two parameters
	 */
	public static double param_dist(Parameter mu_1, Parameter mu_2)
	{
		// Default distance is Euclidean norm
		double sum = 0.;

		for(int i=0; i<mu_1.getNEntries(); i++)
		{
			sum += Math.pow(mu_1.getEntry(i) - mu_2.getEntry(i),2.);
		}

		return Math.sqrt(sum);
	}
	
	/**
	 * Load the current_parameters from the set C_J.
	 */
	protected void get_current_parameters_from_C_J(int index) {
		current_parameters = C_J.get(index);
	}

	/**
	 * Save current_parameters in saved_parameters.
	 */
	protected void save_current_parameters() {
		saved_parameters = current_parameters.clone();
	}

	/**
	 * Reload from saved_parameters
	 */
	protected void reload_current_parameters() {

		// Just point the current_parameters reference to
		// saved_parameters. No need to clone in this case.
		current_parameters = saved_parameters;
	}

	/**
	 * @return the current parameters
	 */
	public Parameter get_current_parameters() {
		return current_parameters;
	}

	/**
	 * @param parameters_filename
	 *            The name of the file to parse Parse the input file to
	 *            initialize this RBSCMSystem.
	 */
	public void parse_parameters_file(String parameters_filename, boolean isAssetFile)
	throws InconsistentStateException {

		GetPot infile = new GetPot(context, parameters_filename, isAssetFile);
		
		//int n_SCM_parameters = infile.call("n_SCM_parameters",1);
		int n_SCM_parameters = infile.call("n_SCM_parameters",infile.call("n_parameters",1));
		Log.d(DEBUG_TAG, "n_parameters = " + n_SCM_parameters);
		
		SCM_M = infile.call("SCM_M",0);

		min_parameter = new Parameter(n_SCM_parameters);
		max_parameter = new Parameter(n_SCM_parameters);
		for(int i=0; i<n_SCM_parameters; i++) {
			// Read in the min/max for the i^th parameter
			String min_string = new String("mu" + i + "_min");
			double mu_i_min = infile.call(min_string, 0.);
			min_parameter.setEntry(i, mu_i_min);

			String max_string = new String("mu" + i + "_max");
			double mu_i_max = infile.call(max_string, 0.);
			max_parameter.setEntry(i, mu_i_max);
		}

		Log.d(DEBUG_TAG,"RBBase parameters from " + parameters_filename + ":");
		for(int i=0; i<n_SCM_parameters; i++) {
			Log.d(DEBUG_TAG,"Parameter " + i +
					": Min = " + getParameterMin(i) +
					", Max = " + getParameterMax(i) );
		}

		Log.d(DEBUG_TAG,"RBSCMSystem parameters from " + parameters_filename + ":");
		Log.d(DEBUG_TAG,"SCM_M: " + SCM_M);
	}
	
	/**
	 * Read in the stored data from the specified URL in order to initialize
	 * the SCM.
	 * 
	 * @param directory_name The URL of the directory containing the Offline data Read in
	 *                       the Offline data to initialize this RBSystem.
	 */
	public void read_offline_data(String directory_name, boolean isAssetFile)
	throws IOException, InconsistentStateException {

		HttpClient client = new DefaultHttpClient();
		
		int buffer_size = 8192;

		// Read in the bounding box minimum values
		{
			
			InputStreamReader isr;
			String dataString = directory_name + "/B_min.dat";
			
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

			B_min = new double[get_Q_a()];
			for(int i=0; i<B_min.length; i++) {
				B_min[i] = Double.parseDouble(tokens[i]);
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading B_min.dat");


		// Read in the bounding box maximum values
		{
			InputStreamReader isr;
			String dataString = directory_name + "/B_max.dat";
			
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

			B_max = new double[get_Q_a()];
			for(int i=0; i<B_max.length; i++) {
				B_max[i] = Double.parseDouble(tokens[i]);
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading B_max.dat");


		// Read in the stability constant values 
		{
			InputStreamReader isr;
			String dataString = directory_name + "/C_J_stability_vector.dat";
			
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
			
			try {
				String[] tokens = line.split(" ");
				
				if( (tokens.length == 1) && (tokens[0] == "") ) {
					C_J_stability_vector = null;
				}
				else { 
					C_J_stability_vector = new double[tokens.length];
					for(int i=0; i<C_J_stability_vector.length; i++) {
						C_J_stability_vector[i] = Double.parseDouble(tokens[i]);
					}
				}
			} catch (Exception e) {
				Log.d(DEBUG_TAG, "Exception occurred when splitting string, " +
						               "setting C_J_stability_vector to null");
				C_J_stability_vector = null;
			}

			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading C_J_stability_vector.dat");

		// Read in C_J, the Greedily selected parameters
		{
			InputStreamReader isr;
			String dataString = directory_name + "/C_J.dat";
			
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

			C_J = new Vector<Parameter>(0);
			if(C_J_stability_vector != null) {
				
				String line = reader.readLine();
				String[] tokens = line.split(" ");				
				
				int count = 0;
				for(int i=0; i<C_J_stability_vector.length; i++) {
					C_J.add( new Parameter(get_n_params()) );
					for(int j=0; j<get_n_params(); j++)
					{
						C_J.get(i).setEntry(j, Double.parseDouble(tokens[count]) );
						count++;
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading C_J.dat");


		// Read in SCM_UB_vectors
		{
			InputStreamReader isr;
			String dataString = directory_name + "/SCM_UB_vectors.dat";
			
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

			if(C_J_stability_vector != null) {
				
				String line = reader.readLine();
				String[] tokens = line.split(" ");
				
				int count = 0;
				// Resize SCM_UB_vectors based on C_J_stability_vector and Q_a
				SCM_UB_vectors = new double[C_J_stability_vector.length][get_Q_a()];
				for(int i=0; i<SCM_UB_vectors.length; i++)
				{
					for(int j=0; j<get_Q_a(); j++) {
						SCM_UB_vectors[i][j] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading SCM_UB_vectors.dat");

	}


	/**
	 * Private helper function to sort the indices of C_J
	 * based on distance from current_parameters
	 */
	private List<Integer> getSorted_CJ_Indices() {

		int J = C_J.size();

		LinkedHashMap<Double,Integer> dist_from_mu =
			new LinkedHashMap<Double,Integer>(J);

		for(int j=0; j<J; j++)
		{
			double dist = param_dist(get_current_parameters(), C_J.get(j) );
			dist_from_mu.put(dist, j);
		}		

		List< Map.Entry<Double,Integer> > list =
			new LinkedList< Map.Entry<Double,Integer> >(dist_from_mu.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
				return ((Comparable) ((Map.Entry) (o1)).getKey())
				.compareTo(((Map.Entry) (o2)).getKey());
			}
		});

		// Create a sorted list of values to return
		List<Integer> result = new LinkedList<Integer>();
		for (Iterator it = list.iterator(); it.hasNext();) {
			Map.Entry<Double,Integer> entry = (Map.Entry<Double,Integer>)it.next();
			result.add(entry.getValue());
		}

		return result;
	}

}

