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

import org.apache.commons.math.complex.Complex;
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

public class RBnSCMCSystem extends RBSCMSystem{
	
	// Logging tag
	private static final String DEBUG_TAG = "RBnSCMCSystem";
	
	private double[] B_min;
	private double[] B_max;
	private int n_mubar;
	private int[] n_muhat;
	private Vector<Parameter> mu_bar;
	private Vector<Parameter>[] mu_hat;
	private double[] beta_bar;
	private double[][] beta_hat;
	private double[][][] zval;

	public RBnSCMCSystem(Context context) {
		super(context);
	}
	
	public void read_offline_data(String directory_name, boolean isAssetFile)
	throws IOException, InconsistentStateException {
		HttpClient client = new DefaultHttpClient();
		
		int buffer_size = 8192;
		
		InputStreamReader isr;
		String dataString = directory_name + "/SCMdata.dat";
		
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

		int count = 0;
		
		B_min = new double[get_Q_a()];
		B_max = new double[get_Q_a()];
		for(int i=0; i<B_min.length; i++) {
			B_max[i] = Double.parseDouble(tokens[count]);
			B_min[i] = -B_max[i];
			count++; 
		}
		
		n_mubar = Integer.parseInt(tokens[count]);
		count++;
		
		mu_bar = new Vector<Parameter>(0);			
		for(int i=0; i<n_mubar; i++) {
			mu_bar.add( new Parameter(get_n_params()) );
			for(int j=0; j<get_n_params(); j++){
				mu_bar.get(i).setEntry(j, Double.parseDouble(tokens[count]) );
				count++;
			}
		}
		
		beta_bar = new double[n_mubar];
		for (int i = 0; i<n_mubar; i++){
			beta_bar[i] = Double.parseDouble(tokens[count]);
			count++;
		}
		
		mu_hat = (Vector<Parameter>[]) new Vector<?>[n_mubar];
		n_muhat = new int[n_mubar];
		beta_hat = new double[n_mubar][];
		zval = new double[n_mubar][][];
		for (int i = 0; i < n_mubar; i++){
			n_muhat[i] = Integer.parseInt(tokens[count]);
			count++;
			beta_hat[i] = new double[n_muhat[i]];
			zval[i] = new double[n_muhat[i]][get_Q_a()*2];		
			
			mu_hat[i] = new Vector<Parameter>(0);
			for(int j=0; j<n_muhat[i]; j++) {
				mu_hat[i].add( new Parameter(get_n_params()) );
				for(int k=0; k<get_n_params(); k++){
					mu_hat[i].get(j).setEntry(k, Double.parseDouble(tokens[count]) );
					count++;
				}
			}
			
			for(int j=0; j<n_muhat[i]; j++) {
				beta_hat[i][j] = Double.parseDouble(tokens[count]);
				count++;
			}
			
			for(int k=0; k<get_Q_a()*2; k++)
				for(int j=0; j<n_muhat[i]; j++){				
					zval[i][j][k] = Double.parseDouble(tokens[count]);
					count++;
			}
		}
		
		reader.close();

		Log.d(DEBUG_TAG, "Finished reading SCMdata.dat");
	
	}
	
	public double get_SCM_LB() {
		//return 0.01;
		
		double min_J_obj = 0.;
		double[] min_Jlocal_obj = new double[n_mubar];
		
		// Sort the indices of mu_bar based on distance from current_parameters
		List<Integer> sortedmubarIndices = getSorted_CJ_Indices(mu_bar);
		int icount = 0;
		//while ((min_J_obj<=0) && (icount < sortedmubarIndices.size())){
		while ((min_J_obj<=0) && (icount < sortedmubarIndices.size())){	
			int imubar = sortedmubarIndices.get(icount);
		
			// First, declare the constraints
			Collection constraints = new ArrayList();

			// Add bounding box constraints for the get_Q_a() variables
			for(int q=0; q<get_Q_a(); q++)
			{
				double[] index = new double[get_Q_a()*2];
				index[q] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ, B_min[q]/beta_bar[imubar]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ, B_max[q]/beta_bar[imubar]));
				
				index[q] = 0.; index[q+get_Q_a()] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ, B_min[q]/beta_bar[imubar]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ, B_max[q]/beta_bar[imubar]));
			}
				
			// Save the current_parameters since we'll change them in the loop below
			save_current_parameters();

			// Add the constraint rows
			if(n_muhat[imubar] > 0) {
				for (int imuhat = 0; imuhat < n_muhat[imubar]; imuhat++) {
					current_parameters = mu_hat[imubar].get(imuhat);
		
					double[] constraint_row = new double[get_Q_a()*2];
					for(int q=0; q<get_Q_a(); q++)
					{
						Complex theta_q_a = complex_eval_theta_q_a(q);
						constraint_row[q] 			= theta_q_a.getReal()*beta_bar[imubar];
						constraint_row[q+get_Q_a()] = theta_q_a.getImaginary()*beta_bar[imubar];
					}
		
					constraints.add(new LinearConstraint(constraint_row, Relationship.GEQ, beta_hat[imubar][imuhat]));
 				}
			}


			// Now load the original parameters back into current_parameters
			// in order to set the coefficients of the objective function
			reload_current_parameters();
				
			// Create objective function object
			double[] objectiveFn = new double[get_Q_a()*2];
			for(int q=0; q<get_Q_a(); q++){
				Complex theta_q_a = complex_eval_theta_q_a(q);
				objectiveFn[q] 			 = theta_q_a.getReal()*beta_bar[imubar];
				objectiveFn[q+get_Q_a()] = theta_q_a.getImaginary()*beta_bar[imubar];
			}
			LinearObjectiveFunction f = new LinearObjectiveFunction(objectiveFn, 0.);
			
			try {
				SimplexSolver solver = new SimplexSolver(1e-6);		
				RealPointValuePair opt_pair = solver.optimize(f, constraints, GoalType.MINIMIZE, false);
				min_Jlocal_obj[icount] = opt_pair.getValue();		
			}
			catch(OptimizationException e) {
				Log.e("RBSCMSYSTEM_TAG", "Optimal solution not found");
				e.printStackTrace();
			}
			catch(Exception e) {
				Log.e("RBSCMSYSTEM_TAG", "Exception occurred during SCM_LB calculation");
				e.printStackTrace();
			}			 
			
			min_J_obj = min_J_obj > min_Jlocal_obj[icount] ? min_J_obj : min_Jlocal_obj[icount];
			icount++;
		}
		return min_J_obj;
	}
	
	public double get_SCM_UB() {
		
		// cheating, since betaUB is rarely good
		double maxbetabar = beta_bar[0];
		for (int i = 0; i < n_mubar; i++)
			maxbetabar = maxbetabar < beta_bar[i] ? beta_bar[i] : maxbetabar;
		return maxbetabar;
		
		/*
		double[] theta_a = new double [get_Q_a()*2]; 
		for(int q=0; q<get_Q_a(); q++){
			Complex theta_q_a = complex_eval_theta_q_a(q);
			theta_a[q] = theta_q_a.getReal();
			theta_a[q+get_Q_a()] = theta_q_a.getImaginary();
		}
		double betaUB = -1e8;
		for (int i = 0; i < n_mubar; i++){
			double localbetaUB = 1e8;
			for (int j = 0; j < n_muhat[i]; j++){
				double calbetaUB = 0.;
				for (int k = 0; k < get_Q_a()*2; k++)
					calbetaUB += zval[i][j][k]*theta_a[k];
				localbetaUB = localbetaUB > calbetaUB ? calbetaUB : localbetaUB;
			}
			betaUB = betaUB < localbetaUB ? localbetaUB : betaUB;
		}
		return betaUB;
		*/
		
	}
	
	private List<Integer> getSorted_CJ_Indices(Vector<Parameter> C_J) {

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

