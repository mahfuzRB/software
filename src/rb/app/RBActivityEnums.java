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

// This class defines some enums that are useful in
// determining which type of RB and SCM systems to
// initialize for an Online solve

public class RBActivityEnums {

	// Enum for the RB system type
	public enum SystemTypeEnum { NONE, LINEAR_STEADY, LINEAR_UNSTEADY, QN_UNSTEADY, LINEAR_COMPLEX_STEADY }
	
	// Enum for the SCM system type
	public enum SCMTypeEnum { NONE, COERCIVE, COERCIVE_ALPHASIGMA, QN_TRANSIENT_SCM, COMPLEX_NONCOERCIVE }

}

