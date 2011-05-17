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

import org.apache.commons.math.linear.RealVector;
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

public class QNTransientSCMSystem extends RBSCMSystem {

	// Logging tag
	private static final String DEBUG_TAG = "QNTransientSCMSystem";

	/**
	 * The current number of basis functions.
	 */
	private int n_bfs;

	/**
	 * The RB coefficients at the C_J parameters.
	 */
	private double[][] C_J_RB_coeffs;
	
	/**
	 * The current RB coefficientss.
	 */
	private RealVector current_RB_coeffs;
	
	// We also need to save the RB coefficients during LB calculations
	private RealVector saved_RB_coeffs;
	
	// We may need to scale the theta_c function for the sake of the SCM!
	private double SCM_theta_c_scaling;
	
	/**
	 * Constructor.
	 */
	public QNTransientSCMSystem(Context context) {
		super(context);
	}
	
	/**
	 * Get/set the number of basis functions
	 */
	int get_n_basis_functions() {
		return n_bfs;
	}
	
	void set_n_basis_functions(int n_bfs_in) {
		n_bfs = n_bfs_in;
	}
	
	/**
	 * Set the current RB coefficients.
	 */
	void set_current_RB_coeffs(RealVector RB_coeffs) {
		current_RB_coeffs = RB_coeffs;
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

		return SCM_theta_c_scaling * theta_val.doubleValue();
	}
	
	/**
	 * Override get_Q_a since we have n_bfs extra terms
	 */
	@Override
	public int get_Q_a() {
		return super.get_Q_a() + get_n_basis_functions();
	}
	
	/**
	 * Override eval_theta_q_a in order to account for the affine terms
	 * related to basis functions
	 */
	@Override
	public double eval_theta_q_a(int q) {

		if(q < get_n_basis_functions()) {
			double theta_c_value = eval_theta_c();
			return current_RB_coeffs.getEntry(q)*theta_c_value;
		}
		else {

			Method meth;

			try {        	
				Class partypes[] = new Class[2];
				partypes[0] = Integer.TYPE;
				partypes[1] = double[].class;

				meth = mAffineFnsClass.getMethod("evaluateA", partypes);
			} catch (NoSuchMethodException nsme) {
				throw new RuntimeException("getMethod for evaluateA failed", nsme);
			}

			Double theta_val;
			try {
				Object arglist[] = new Object[2];
				arglist[0] = new Integer(q-get_n_basis_functions());
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
	}
	
	/**
	 * Override to also load the RB coefficients.
	 */
	@Override
	protected void get_current_parameters_from_C_J(int index) {
		super.get_current_parameters_from_C_J(index);
		
		for(int i=0; i<get_n_basis_functions(); i++)
			current_RB_coeffs.setEntry(i, C_J_RB_coeffs[index][i]);
	}
	
	/**
	 * Override to also save the RB coefficients.
	 */
	@Override
	protected void save_current_parameters() {
		super.save_current_parameters();
		
		saved_RB_coeffs = current_RB_coeffs.copy();
	}

	/**
	 * Override to also load the RB coefficients.
	 */
	@Override
	protected void reload_current_parameters() {
		super.reload_current_parameters();
		
		set_current_RB_coeffs(saved_RB_coeffs);
	}
	
	/**
	 * Override read_offline_data in order to read in the extra data
	 * in the QNTransient case.
	 */
	@Override
	public void read_offline_data(String directory_name, boolean isAssetFile)
	throws IOException, InconsistentStateException {
		
		HttpClient client = new DefaultHttpClient();
		
		int buffer_size = 8192;
		
		// Initially set number of basis functions from n_bfs.dat
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
		
		
		super.read_offline_data(directory_name, isAssetFile);



		// Read C_J_RB_coeffs
		{
			InputStreamReader isr;
			String dataString = directory_name + "/C_J_RB_coeffs.dat";
			
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

			C_J_RB_coeffs = new double[C_J_stability_vector.length][get_n_basis_functions()];
			if(C_J_stability_vector != null) {
				
				String line = reader.readLine();
				String[] tokens = line.split(" ");
				
				int count = 0;
				for(int i=0; i<C_J_stability_vector.length; i++) {
					for(int j=0; j<get_n_basis_functions(); j++)
					{
						C_J_RB_coeffs[i][j] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading C_J_RB_coeffs.dat");
	}
	
	/**
	 * @param parameters_filename
	 *            The name of the file to parse Parse the input file to
	 *            initialize this RBSCMSystem.
	 */
	@Override
	public void parse_parameters_file(String parameters_filename, boolean isAssetFile)
	throws InconsistentStateException {
		super.parse_parameters_file(parameters_filename, isAssetFile);
		
		GetPot infile = new GetPot(context, parameters_filename, isAssetFile);
		
		SCM_theta_c_scaling = infile.call("SCM_theta_c_scaling", 1.);
		Log.d(DEBUG_TAG, "SCM_theta_c_scaling = " + SCM_theta_c_scaling);
	}

}

