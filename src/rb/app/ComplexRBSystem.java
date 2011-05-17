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
import org.apache.commons.math.linear.Array2DRowFieldMatrix;
import org.apache.commons.math.linear.ArrayFieldVector;
import org.apache.commons.math.linear.FieldMatrix;
import org.apache.commons.math.linear.FieldVector;
import org.apache.commons.math.linear.RealVector;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;

import android.content.Context;
import android.util.Log;

public class ComplexRBSystem extends RBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "ComplexRBSystem";
	
	public Complex[] RB_outputs;
	public Complex[] RB_output_error_bounds;
	
	protected FieldVector<Complex> RB_solution;
		
	public Complex[][] output_dual_norms;
	protected FieldVector<Complex>[][] RB_output_vectors;
	protected FieldMatrix<Complex>[] RB_A_q_vector;
	protected FieldVector<Complex>[] RB_F_q_vector;
	
	protected FieldVector<Complex> theta_a;
	
	protected Complex[] Fq_representor_norms;
	protected Complex[][][] Fq_Aq_representor_norms;
	protected Complex[][][] Aq_Aq_representor_norms;
	
	protected Complex[][][] Z_vector;
	protected Complex[][] uL_vector;
	
	public ComplexRBSystem(Context context) {
		super(context);
	}	

	@Override
	public void read_offline_data(String directory_name, boolean isAssetFile)
	throws IOException {

		isReal = false;
		
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
		
			set_n_basis_functions(Integer.parseInt(line));
		
			reader.close();
		}
		Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");
		
		// Read in output data
		if (get_n_outputs() > 0) {
			// Get output dual norms
			{
		
				RB_output_vectors = (FieldVector<Complex>[][]) new FieldVector<?>[get_n_outputs()][];
				output_dual_norms = new Complex[get_n_outputs()][];
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
					output_dual_norms[i] = new Complex[Q_l_hat];
					for(int q=0; q<Q_l_hat; q++) {
						output_dual_norms[i][q] = new Complex(Double.parseDouble(dual_norms_tokens[q]),
															  Double.parseDouble(dual_norms_tokens[Q_l_hat+q]));
					}
					
					
					RB_output_vectors[i] = (FieldVector<Complex>[]) new FieldVector<?>[get_Q_l(i)];					
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
		
						RB_output_vectors[i][q_l] = new ArrayFieldVector<Complex>(get_n_basis_functions(), new Complex(0d,0d));
						for (int j = 0; j < get_n_basis_functions(); j++) {
							RB_output_vectors[i][q_l].setEntry(j, new Complex(Double.parseDouble(output_i_tokens[j]),
									Double.parseDouble(output_i_tokens[get_n_basis_functions()+j])));
						}
						reader_i.close();
						
					}
		
					reader1.close();
				}
			}
		}
		
		Log.d(DEBUG_TAG, "Finished reading output data");
		/*
		// Read in the inner product matrix
		{
			InputStreamReader isr;
			String dataString = directory_name + "/RB_inner_product_matrix.dat";
			
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
			RB_inner_product_matrix = new Array2DRowRealMatrix(n_bfs, n_bfs);
		
			// Fill the matrix
			int count = 0;
			for (int i = 0; i < n_bfs; i++)
				for (int j = 0; j < n_bfs; j++) {
					RB_inner_product_matrix.setEntry(i, j, Double
							.parseDouble(tokens[count]));
					count++;
				}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading RB_inner_product_matrix.dat");
		*/
		
		// Read in the F_q vectors
		RB_F_q_vector = (FieldVector<Complex>[]) new FieldVector<?>[get_Q_f()];
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
			RB_F_q_vector[q_f] = new ArrayFieldVector<Complex>(get_n_basis_functions(), new Complex(0d,0d));
		
			// Fill the vector
			for (int i = 0; i < get_n_basis_functions(); i++) {
				RB_F_q_vector[q_f].setEntry(i, new Complex(Double.parseDouble(tokens[i]),Double.parseDouble(tokens[get_n_basis_functions()+i])));
			}
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading RB_F_q data");
		
		// Read in the A_q matrices
		RB_A_q_vector = (FieldMatrix<Complex>[]) new FieldMatrix<?>[get_Q_a()];
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
			RB_A_q_vector[q_a] = new Array2DRowFieldMatrix<Complex>((new Complex(0,0)).getField(), get_n_basis_functions(),get_n_basis_functions());
					
			// Fill the vector
			int count = 0;
			for (int i = 0; i < get_n_basis_functions(); i++)
				for (int j = 0; j < get_n_basis_functions(); j++) {
					RB_A_q_vector[q_a].setEntry(i, j, new Complex(Double.parseDouble(tokens[count]),
															      Double.parseDouble(tokens[count+get_n_basis_functions()*get_n_basis_functions()])));
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
			Fq_representor_norms = new Complex[Q_f_hat];
		
			// Fill it
			for (int i = 0; i < Q_f_hat; i++) {
				Fq_representor_norms[i] = new Complex(Double.parseDouble(tokens[i*2+0]),
													  Double.parseDouble(tokens[i*2+1]));
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
			Fq_Aq_representor_norms = new Complex[get_Q_f()][get_Q_a()][get_n_basis_functions()];
		
			double[][][] Rdata = new double[get_Q_f()][get_Q_a()][get_n_basis_functions()];
			double[][][] Idata = new double[get_Q_f()][get_Q_a()][get_n_basis_functions()];
			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < get_Q_f(); q_f++)
				for (int q_a = 0; q_a < get_Q_a(); q_a++){
					for (int i = 0; i < get_n_basis_functions(); i++){
						Rdata[q_f][q_a][i] = Double.parseDouble(tokens[count]);
						count++;
					}
					for (int i = 0; i < get_n_basis_functions(); i++){
						Idata[q_f][q_a][i] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			
			for (int q_f = 0; q_f < get_Q_f(); q_f++)
				for (int q_a = 0; q_a < get_Q_a(); q_a++)
					for (int i = 0; i < get_n_basis_functions(); i++)
						Fq_Aq_representor_norms[q_f][q_a][i] = new Complex(Rdata[q_f][q_a][i],Idata[q_f][q_a][i]);
				
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading Fq_Aq_norms.dat");
		
		// Read in Aq_Aq representor norm data
		{
			// Declare the array
			int Q_a_hat = get_Q_a() * (get_Q_a() + 1) / 2;
			Aq_Aq_representor_norms = new Complex[Q_a_hat][get_n_basis_functions()][get_n_basis_functions()];
		
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
		
					double[][] Rdata = new double[get_n_basis_functions()][get_n_basis_functions()];
					double[][] Idata = new double[get_n_basis_functions()][get_n_basis_functions()];
					for (int k = 0; k < get_n_basis_functions(); k++)
						for (int l = 0; l < get_n_basis_functions(); l++)
							Rdata[k][l] = dis.ReadDouble();
					for (int k = 0; k < get_n_basis_functions(); k++)
						for (int l = 0; l < get_n_basis_functions(); l++)
							Idata[k][l] = dis.ReadDouble();
					for (int k = 0; k < get_n_basis_functions(); k++)
						for (int l = 0; l < get_n_basis_functions(); l++)
							Aq_Aq_representor_norms[count][k][l] = new Complex(Rdata[k][l],Idata[k][l]);
					
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
						
			set_calN(Integer.parseInt(line));
			reader.close();
		}
		
		Log.d(DEBUG_TAG, "Finished reading calN.dat");
		
		// Reading uL data
		if (get_Q_uL()>0){
			uL_vector = new Complex[get_Q_uL()][get_calN()];
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
				
				float[] Rdata = new float[get_calN()];
				float[] Idata = new float[get_calN()];
				Rdata = dis.ReadFloat(get_calN()); // read in the real part
				Idata = dis.ReadFloat(get_calN()); // read in the imagine part
				
				for (int i = 0; i < get_calN(); i++)
					uL_vector[q_uL][i] = new Complex((double)Rdata[i],(double)Idata[i]);
				
				dis.close();
			}
		}
		Log.d(DEBUG_TAG, "Finished reading uL.dat");
		
		// Read in Z data
		if (get_mfield() > 0){ 
			Z_vector = new Complex[get_mfield()][get_n_basis_functions()][get_calN()];
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
					
					float[] Rdata = new float[get_calN()];
					float[] Idata = new float[get_calN()];
					Rdata = dis.ReadFloat(get_calN()); // read in the real part
					Idata = dis.ReadFloat(get_calN()); // read in the imagine part
					
					for (int i = 0; i < get_calN(); i++)
						Z_vector[imf][inbfs][i] = new Complex((double)Rdata[i],(double)Idata[i]);
					
					dis.close();
				}
		}
				
		Log.d(DEBUG_TAG, "Finished reading Z.dat");
		
		initialize_data_vectors();
	}
	
	protected void initialize_data_vectors() {
		// Also, resize RB_outputs and RB_output_error_error_bounds arrays
		RB_outputs = new Complex[get_n_outputs()];
		RB_output_error_bounds = new Complex[get_n_outputs()];
	}

	@Override
	public double RB_solve(int N) {
		
		current_N = N;
		
		theta_a = complex_eval_theta_q_a();

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
		FieldMatrix<Complex> RB_system_matrix_N = new Array2DRowFieldMatrix<Complex>((new Complex(0,0)).getField(), N, N);

		for (int q_a = 0; q_a < get_Q_a(); q_a++) {
			RB_system_matrix_N = RB_system_matrix_N.add (
				RB_A_q_vector[q_a].getSubMatrix(0, N-1, 0, N-1).
				scalarMultiply(theta_a.getEntry(q_a)));
//				scalarMultiply(complex_eval_theta_q_a(q_a) ) );
		}

		// Assemble the RB rhs
		FieldVector<Complex> RB_rhs_N = new ArrayFieldVector<Complex>(N, new Complex(0d,0d));

		for (int q_f = 0; q_f < get_Q_f(); q_f++) {
			// Note getSubVector takes an initial index and the number of entries
			// i.e. the interface is a bit different to getSubMatrix
			RB_rhs_N = RB_rhs_N.add(RB_F_q_vector[q_f].getSubVector(0, N)
					.mapMultiply(complex_eval_theta_q_f(q_f)));
		}

		// Solve the linear system by Gaussian elimination

		RB_solution = new ArrayFieldVector<Complex>(N, new Complex(0., 0.));
		for (int j = 1; j < N; j++)
			for (int i = j; i < N; i++){
				Complex m = RB_system_matrix_N.getEntry(i, j-1).divide(RB_system_matrix_N.getEntry(j-1, j-1));
				for (int k = 0; k < N; k++)
					RB_system_matrix_N.setEntry(i, k, RB_system_matrix_N.getEntry(i,k)
							.subtract(RB_system_matrix_N.getEntry(j-1,k).multiply(m)));
				RB_rhs_N.setEntry(i, RB_rhs_N.getEntry(i).subtract(m.multiply(RB_rhs_N.getEntry(j-1))));
				}
		RB_solution.setEntry(N-1, RB_rhs_N.getEntry(N-1).divide(RB_system_matrix_N.getEntry(N-1,N-1)));
		for (int j = N-2; j >= 0; j--){
			Complex m = new Complex(0., 0.);
			for (int i = j+1; i<N; i++)
				m = m.add(RB_system_matrix_N.getEntry(j, i).multiply(RB_solution.getEntry(i)));
			RB_solution.setEntry(j, (RB_rhs_N.getEntry(j).subtract(m)).divide(RB_system_matrix_N.getEntry(j, j)));
		}

		// Evaluate the dual norm of the residual for RB_solution_vector
		double epsilon_N = compute_residual_dual_norm(N);

		// Get lower bound for coercivity constant
		double alpha_LB = get_SCM_lower_bound();
				
		// If SCM lower bound is negative
		if (alpha_LB<0){ // Get an upper bound instead
			alpha_LB = get_SCM_upper_bound();
		}

		// Store (absolute) error bound
		double abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

		// Compute the norm of RB_solution
		/*
		RealMatrix RB_inner_product_matrix_N =
			RB_inner_product_matrix.getSubMatrix(0, N-1, 0, N-1);
		*/

		double RB_solution_norm = 0.0d;
		for (int i = 0; i < N; i++)
			RB_solution_norm += ((RB_solution.getEntry(i)).multiply((RB_solution.getEntry(i)).conjugate())).getReal();
		RB_solution_norm = Math.sqrt(RB_solution_norm);

		// Now compute the outputs and associated errors
		FieldVector<Complex> RB_output_vector_N = new ArrayFieldVector<Complex>(N, new Complex(0d,0d));
		for (int i = 0; i < get_n_outputs(); i++) {
			RB_outputs[i] = new Complex(0.,0.);
			
			RB_output_vector_N = (RB_output_vectors[i][0].getSubVector(0, N)).mapMultiply(complex_eval_theta_q_l(i, 0));
			for(int q_l = 1; q_l < get_Q_l(i); q_l++)
				RB_output_vector_N = RB_output_vector_N.add((RB_output_vectors[i][q_l].getSubVector(0, N)).mapMultiply(complex_eval_theta_q_l(i, q_l)));
			for (int j = 0; j < N; j++)
				RB_outputs[i] = RB_outputs[i].add((RB_output_vector_N.getEntry(j).conjugate()).multiply((RB_solution.getEntry(j))));
						
			RB_output_error_bounds[i] = new Complex(compute_output_dual_norm(i) * abs_error_bound, compute_output_dual_norm(i) * abs_error_bound);
		}

		cal_derived_output();
		
		return (return_rel_error_bound ? abs_error_bound / RB_solution_norm
				: abs_error_bound);
	}

	@Override
	protected double compute_residual_dual_norm(int N) {

		// Use the stored representor inner product values
		// to evaluate the residual norm
		double res_ff = 0;
		double res_af = 0;
		double res_aa = 0;

		int q = 0;
		for (int q_f1 = 0; q_f1 < get_Q_f(); q_f1++) {
			for (int q_f2 = q_f1; q_f2 < get_Q_f(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				res_ff += delta *( (complex_eval_theta_q_f(q_f1).multiply(complex_eval_theta_q_f(q_f2).conjugate()))
					                .multiply(Fq_representor_norms[q]) ).getReal();				
				q++;
			}
		}

		for (int q_f = 0; q_f < get_Q_f(); q_f++) {
			for (int q_a = 0; q_a < get_Q_a(); q_a++) {
				for (int i = 0; i < N; i++) {
					res_af += 2.*( get_complex_soln_coeff(i).conjugate().multiply(
					//complex_eval_theta_q_f(q_f).multiply(complex_eval_theta_q_a(q_a).conjugate())
					complex_eval_theta_q_f(q_f).multiply(theta_a.getEntry(q_a).conjugate())
					).multiply(Fq_Aq_representor_norms[q_f][q_a][i]) ).getReal();
				}
			}
		}
		
		q = 0;
		for (int q_a1 = 0; q_a1 < get_Q_a(); q_a1++) {
			for (int q_a2 = q_a1; q_a2 < get_Q_a(); q_a2++) {
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						double delta = (q_a1 == q_a2) ? 1. : 2.;
						res_aa += delta*( (get_complex_soln_coeff(i).conjugate().multiply(get_complex_soln_coeff(j))).multiply(
									      //(complex_eval_theta_q_a(q_a1).conjugate().multiply(complex_eval_theta_q_a(q_a2)))
								           (theta_a.getEntry(q_a1).conjugate().multiply(theta_a.getEntry(q_a2)))
										   ).multiply(Aq_Aq_representor_norms[q][i][j]) ).getReal();
						}
				}
				q++;
			}
		}

		double residual_norm_sq = res_ff + res_af + res_aa;
		
		if (residual_norm_sq < 0.) {
			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect error bound much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}
	
	@Override
	protected double compute_output_dual_norm(int i) {
		
		// Use the stored representor inner product values
		// to evaluate the output dual norm
		double output_norm_sq = 0.;
		
		int q = 0;
		for (int q_l1 = 0; q_l1 < get_Q_l(i); q_l1++) {
			for (int q_l2 = q_l1; q_l2 < get_Q_l(i); q_l2++) {
				if (q_l1 == q_l2)
					output_norm_sq += 1.*( (complex_eval_theta_q_l(i, q_l1).multiply(complex_eval_theta_q_l(i, q_l2).conjugate()))
									       .multiply(output_dual_norms[i][q]) ).getReal();
				else
					output_norm_sq += 2.*( (complex_eval_theta_q_l(i, q_l1).multiply(complex_eval_theta_q_l(i, q_l2).conjugate()))
						                   .multiply(output_dual_norms[i][q]) ).getReal();
				q++;
			}
		}
		
		return Math.sqrt(output_norm_sq);
	}

	Complex get_complex_soln_coeff(int i) {
		return RB_solution.getEntry(i);
	}
	
	@Override
	double[][] get_RBsolution(){
		Complex[] c = RB_solution.toArray();
		double[][] d = new double[2][c.length];
		for (int i = 0; i < c.length; i++){
			d[0][i] = c[i].getReal();
			d[1][i] = c[i].getImaginary();
		}
		return d;
	}
	
	@Override
	double get_RB_output(int n_output, boolean Rpart){
		if (Rpart)
			return RB_outputs[n_output].getReal();
		else
			return RB_outputs[n_output].getImaginary();
	}
	
	@Override
	double get_RB_output_error_bound(int n_output, boolean Rpart){
		if (Rpart)
			return RB_output_error_bounds[n_output].getReal();
		else
			return RB_output_error_bounds[n_output].getImaginary();
	}

	@Override
	public float[][][] get_truth_sol(){
		int N = RB_solution.getDimension();
		float[][][] truth_sol = new float[get_mfield()][3][get_calN()];
		for (int ifn = 0; ifn < get_mfield(); ifn++){
			Complex tmpval;
			for (int i = 0; i < get_calN(); i++){
				tmpval = new Complex(0., 0.);
				for (int j = 0; j < N; j++)
					tmpval = tmpval.add(Z_vector[ifn][j][i].multiply(get_complex_soln_coeff(j)));
				truth_sol[ifn][0][i] = (float) tmpval.getReal();
				truth_sol[ifn][1][i] = (float) tmpval.getImaginary();
			}
			if (get_Q_uL()>0){
				for (int q_uL = 0; q_uL < get_Q_uL(); q_uL++)
					for (int i = 0; i < get_calN(); i++){
						truth_sol[ifn][0][i] += uL_vector[q_uL][i].getReal();
						truth_sol[ifn][1][i] += uL_vector[q_uL][i].getImaginary();
					}
			}
			for (int i = 0; i < get_calN(); i++)
				truth_sol[ifn][2][i] = (float) Math.sqrt(truth_sol[ifn][0][i]*truth_sol[ifn][0][i] + truth_sol[ifn][1][i]*truth_sol[ifn][1][i]);							
		}
		return truth_sol;
	}
	
	@Override
	public float[][][] get_sweep_truth_sol(){
		int N = RB_sweep_solution[0][0].length;
		int numSweep = RB_sweep_solution.length;
		
		Complex[][] RB_sweep_sol = new Complex[numSweep][N];
		for (int i = 0; i < numSweep; i++)
			for (int j = 0; j < N; j++)
				RB_sweep_sol[i][j] = new Complex(RB_sweep_solution[i][0][j],RB_sweep_solution[i][1][j]);
		
		float[][][] truth_sol = new float[get_mfield()][3][get_calN()*numSweep];
		for (int ifn = 0; ifn < get_mfield(); ifn++){
			Complex tmpval;
			for (int iSweep = 0; iSweep < numSweep; iSweep++){
				for (int i = 0; i < get_calN(); i++){
					tmpval = new Complex(0., 0.);
					for (int j = 0; j < N; j++)
						tmpval = tmpval.add(Z_vector[ifn][j][i].multiply(RB_sweep_sol[iSweep][j]));
					truth_sol[ifn][0][iSweep*get_calN()+i] = (float) tmpval.getReal();
					truth_sol[ifn][1][iSweep*get_calN()+i] = (float) tmpval.getImaginary();
				}
			}
			if (get_Q_uL()>0){
				for (int q_uL = 0; q_uL < get_Q_uL(); q_uL++)
					for (int iSweep = 0; iSweep < numSweep; iSweep++)
						for (int i = 0; i < get_calN(); i++){
							truth_sol[ifn][0][iSweep*get_calN()+i] += uL_vector[q_uL][i].getReal();
							truth_sol[ifn][1][iSweep*get_calN()+i] += uL_vector[q_uL][i].getImaginary();
						}
			}
			for (int i = 0; i < get_calN()*numSweep; i++)
				truth_sol[ifn][2][i] = (float) Math.sqrt(truth_sol[ifn][0][i]*truth_sol[ifn][0][i] + truth_sol[ifn][1][i]*truth_sol[ifn][1][i]);							
		}
		return truth_sol;
	}

	public boolean is_derived_output() {		
        Method meth;       

        try {
        	Class partypes[] = null;
            meth = mAffineFnsClass.getMethod("is_derived_output", partypes);
        } catch (NoSuchMethodException nsme) {
            return false;
        }
        
        try {
        	Object arglist[] = null;
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	boolean val = (Boolean) theta_obj;
        	return val;        	
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }        
	}
	
	public void cal_derived_output(){
		if (is_derived_output()){
			Method meth;       

	        try {        	
	        	
	        	Class partypes[] = new Class[1];
	        	partypes[0] = double[].class;
	        	
	            meth = mAffineFnsClass.getMethod("cal_derived_output", partypes);
	        } catch (NoSuchMethodException nsme) {
	        	throw new RuntimeException("getMethod for cal_derived_output failed", nsme);
	        }
	        
	        try {
	        	
	        	Object arglist[] = new Object[1];
	        	double[] input = new double[4];
	        	for (int i = 0; i < get_n_outputs(); i++){
		        	input[0] = RB_outputs[i].getReal();
		        	input[1] = RB_outputs[i].getImaginary();
		        	input[2] = RB_output_error_bounds[i].getReal();
		        	input[3] = RB_output_error_bounds[i].getImaginary();		        	
		        	
		        	arglist[0] = input;            	            
		        	
		        	Object theta_obj = meth.invoke(mTheta, arglist);
		        	double[] output = (double[]) theta_obj;
		        	RB_outputs[i] = new Complex(output[0],output[1]);
		        	RB_output_error_bounds[i] = new Complex(output[2],output[3]);
	        	}	                	
	        }
	        catch (IllegalAccessException iae) {
	            throw new RuntimeException(iae);
	        }
	        catch (InvocationTargetException ite) {
	            throw new RuntimeException(ite.getCause());
	        }       
			
		}		
	}
}

