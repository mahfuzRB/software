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

// A simple class to define a parameter vector

public class Parameter {

	/**
	 * The actual parameter values, stored in a Vector.
	 */
	private double[] mValues;
	
	/**
	 * Constructor, clone the input array.
	 */
	public Parameter(double[] param_in) {
		mValues = (double[])param_in.clone();
	}
	
	/**
	 * Constructor, set the size of mValues based on input
	 * argument and set mValues entries to zero.
	 */
	public Parameter(int size) {
		mValues = new double[size];
		for(int i=0; i<mValues.length; i++)
			mValues[i] = 0.;
	}
	
	/**
	 * @param i
	 * @return the i^th entry of this parameter.
	 */
	public double getEntry(int i) {
		return mValues[i];
	}
	
	/**
	 * @param i Index to set
	 * @param value The new value
	 */
	public void setEntry(int i, double value) {
		mValues[i] = value;
	}
	
	/**
	 * @return the number of entries contained in this parameter.
	 */
	public int getNEntries() {
		return mValues.length;
	}
	
	/**
	 * Clone this Parameter object.
	 */
	public Parameter clone() {
		Parameter new_p = new Parameter(getNEntries());
		for(int i=0; i<getNEntries(); i++) {
			new_p.setEntry(i, getEntry(i) );
		}
		
		return new_p;
	}
	
	/**
	 * Return a String representation of this parameter.
	 */
	public String toString() {
		String s = "";
		for(int i=0; i<this.getNEntries(); i++)
			s += this.getEntry(i) + " ";
		return s;
	}
	
	/**
	 * Return a reference to the internal array. 
	 */
	protected double[] getArray() {
		return mValues;
	}
	
}

