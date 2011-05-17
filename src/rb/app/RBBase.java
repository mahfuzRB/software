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

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.linear.ArrayFieldVector;
import org.apache.commons.math.linear.FieldVector;

import android.content.Context;
import android.util.Log;

// Base class for RB and SCM systems, stores
// the current parameter value, parameter ranges,
// as well as other data common to RBSystems and
// SCMSystems.
// This class is modeled on the RBBase class in rbOOmit

public class RBBase {

	// Logging tag
	private static final String DEBUG_TAG = "RBBase";
	
	/**
	 * The Android Context of RBActivity.
	 */
	public Context context;
	
	/**
	 * Vector storing the current parameters.
	 */
	protected Parameter current_parameters;
	
	/**
	 * Parameters that define the corners of \calD.
	 */
	protected Parameter min_parameter;
	protected Parameter max_parameter;
	
	protected boolean isReal = true;

	/**
	 * The Class object for the affine functions class.
	 */
	protected Class mAffineFnsClass;
	
	/**
	 * The member object that defines the parameter-dependent functions for the
	 * affine expansion of the left-hand, right-hand side and outputs.
	 * We need to access this object at runtime, hence we use a ClassLoader and
	 * reflection in order to call the relevant functions.
	 */
	protected Object mTheta;
	
	/**
	 * The number of terms in the affine expansion of the bilinear form
	 */
	private int mQ_a;
	
	/**
	 * Constructor. Sets the Android context member.
	 */
	public RBBase(Context context_in) {
		context = context_in;
	}
	

	/**
	 * Set the current_parameters member variable
	 */
	public void setCurrentParameters(Parameter p) {
		current_parameters = p;
	}

	/**
	 * @return the current parameters
	 */
	public Parameter getCurrentParameters() {
		return current_parameters;
	}
	
	/**
	 * Print out the parameters stored in current_parameters.
	 */
	public void printCurrentParameters() {
		Log.d(DEBUG_TAG, "Current parameters: "
				+ current_parameters.toString());
	}

	/**
	 * Get the number of parameters. Value is determined by specifying the
	 * parameter ranges.
	 */
	public int get_n_params() throws InconsistentStateException {
		int mu_min_size = min_parameter.getNEntries();
		int mu_max_size = max_parameter.getNEntries();

		if (mu_min_size != mu_max_size) {
			throw new InconsistentStateException(
					"Inconsistent Parameter sizes in get_n_params()");
		} else {
			return mu_min_size;
		}
	}
	

	/**
	 * Get minimum parameter value for i^th parameter.
	 */
	public double getParameterMin(int i) {
		return min_parameter.getEntry(i);
	}

	/**
	 * Get maximum parameter value for i^th parameter.
	 */
	public double getParameterMax(int i) {
		return max_parameter.getEntry(i);
	}
	

	/**
	 * @return Q_a, the number of term in the affine expansion
	 * of the bilinear form
	 */
	public int get_Q_a() {
		return mQ_a;
	}
	
	/**
	 * Evaluate theta_q_a (for the q^th bilinear form) at the current parameter.
	 */
	public double eval_theta_q_a(int q) {		
        Method meth;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
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
	
	public Complex complex_eval_theta_q_a(int q) {		
        Method meth;
        //Complex c;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[3];
        	partypes[0] = Integer.TYPE;
            partypes[1] = double[].class;
            partypes[2] = boolean.class;
        	
            meth = mAffineFnsClass.getMethod("evaluateA", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateA failed", nsme);
        }
        
        Double theta_val_r, theta_val_i;
        try {
            Object arglist[] = new Object[3];
            arglist[0] = new Integer(q);
            arglist[1] = current_parameters.getArray();
            arglist[2] = true;
        	
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	theta_val_r = (Double) theta_obj;
        	
        	arglist[2] = false;
        	theta_val_i = (Double) meth.invoke(mTheta, arglist);
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
        Complex c = new Complex(theta_val_r.doubleValue(),theta_val_i.doubleValue());
        
        return c;
	}
	
	public FieldVector<Complex> complex_eval_theta_q_a() {		
        Method meth;
        //Complex c;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[1];        	
            partypes[0] = double[].class;            
        	
            meth = mAffineFnsClass.getMethod("evaluateA_array", partypes);
        } catch (NoSuchMethodException nsme) {
        	FieldVector<Complex> c = new ArrayFieldVector<Complex>(mQ_a, new Complex(0d,0d));
        	for (int i = 0; i < mQ_a; i++)
            	c.setEntry(i,complex_eval_theta_q_a(i));
        	return c;
            //throw new RuntimeException("getMethod for evaluateA failed", nsme);
        }
        
        double[][] theta_val;
        try {
            Object arglist[] = new Object[1];            
            arglist[0] = current_parameters.getArray();
                    	
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	theta_val = (double[][]) theta_obj;
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
        FieldVector<Complex> c = new ArrayFieldVector<Complex>(mQ_a, new Complex(0d,0d));
        for (int i = 0; i < mQ_a; i++)
        	c.setEntry(i,new Complex(theta_val[i][0],theta_val[i][1]));
        
        return c;
	}
	
	public Complex complex_eval_theta_q_f(int q) {		
        Method meth;
        //Complex c;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[3];
        	partypes[0] = Integer.TYPE;
            partypes[1] = double[].class;
            partypes[2] = boolean.class;
        	
            meth = mAffineFnsClass.getMethod("evaluateF", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateF failed", nsme);
        }
        
        Double theta_val_r, theta_val_i;
        try {
            Object arglist[] = new Object[3];
            arglist[0] = new Integer(q);
            arglist[1] = current_parameters.getArray();
            arglist[2] = true;
        	
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	theta_val_r = (Double) theta_obj;
        	
        	arglist[2] = false;
        	theta_val_i = (Double) meth.invoke(mTheta, arglist);
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
        Complex c = new Complex(theta_val_r.doubleValue(),theta_val_i.doubleValue());
        
        return c;
	}
	
	public Complex complex_eval_theta_q_l(int n, int q) {		
        Method meth;
        //Complex c;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[4];
        	partypes[0] = Integer.TYPE;
        	partypes[1] = Integer.TYPE;
            partypes[2] = double[].class;
            partypes[3] = boolean.class;
        	
            meth = mAffineFnsClass.getMethod("evaluateL", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for evaluateL failed", nsme);
        }
        
        Double theta_val_r, theta_val_i;
        try {
            Object arglist[] = new Object[4];
            arglist[0] = new Integer(n);
            arglist[1] = new Integer(q);
            arglist[2] = current_parameters.getArray();
            arglist[3] = true;
        	
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	theta_val_r = (Double) theta_obj;
        	
        	arglist[3] = false;
        	theta_val_i = (Double) meth.invoke(mTheta, arglist);
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
        Complex c = new Complex(theta_val_r.doubleValue(),theta_val_i.doubleValue());
        
        return c;
	}

	/**
	 * Set the Q_a variable from the mTheta object.
	 */
	protected void read_in_Q_a() {
		Method meth;

		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = mAffineFnsClass.getMethod("get_n_A_functions", null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_A_functions failed", nsme);
		}

		Integer Q_a;
		try {
			Object Q_a_obj = meth.invoke(mTheta, null);
			Q_a = (Integer) Q_a_obj;
		}
		catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		}
		catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		mQ_a = Q_a.intValue();
	}
	
	public boolean is_custom_mesh_transform() {		
        Method meth;       

        try {
        	Class partypes[] = null;
            meth = mAffineFnsClass.getMethod("is_custom_mesh_transform", partypes);
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

	public float[] mesh_transform(double[] mu, float[] x) {		
        Method meth;

        try {
            // Get a reference to get_n_L_functions, which does not
            // take any arguments
        	
        	Class partypes[] = new Class[2];        	
            partypes[0] = double[].class;
            partypes[1] = float[].class;
        	
            meth = mAffineFnsClass.getMethod("mesh_transform", partypes);
        } catch (NoSuchMethodException nsme) {
            throw new RuntimeException("getMethod for mesh_transform failed", nsme);
        }
        
        float[] xt;
        try {
            Object arglist[] = new Object[2];            
            arglist[0] = mu;
            arglist[1] = x;
        	
        	Object theta_obj = meth.invoke(mTheta, arglist);
        	xt = (float[]) theta_obj;
        }
        catch (IllegalAccessException iae) {
            throw new RuntimeException(iae);
        }
        catch (InvocationTargetException ite) {
            throw new RuntimeException(ite.getCause());
        }
                
        return xt;
	}
	
}

