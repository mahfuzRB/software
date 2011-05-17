
public class AffineFunctions {

	public double evaluateM(int i, double[] p) {

		if( (i<0) || (i>get_n_M_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateM()," +
				" i = " + i + " but get_n_M_functions() = " + get_n_M_functions());
		}

		return 1.;
	}
	
	public double[] evaluateDiffM(int i, double[] p) {

		if( (i<0) || (i>get_n_M_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffM()," +
				" i = " + i + " but get_n_M_functions() = " + get_n_M_functions());
		}

		double[] diff = {0., 0.};
		return diff;
	}
	
	public double evaluateA(int i, double[] p) {

		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}

		switch(i) {
			case 0:
				return  p[0];
			case 1:
				return  p[1];
			case 2:
				return  1.;
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double[] evaluateDiffA(int i, double[] p) {

		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}

		double[] diff = {0., 0.};
		switch(i) {
			case 0:
				diff[0] = 1.;
				return  diff;
			case 1:
				diff[1] = 1.;
				return  diff;
			case 2:
				return  diff;
			default:
				throw new Error("Should not reach here");
		}
	}

	public double evaluateF(int i, double[] p) {

		if( (i<0) || (i>get_n_F_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateF()," +
				" i = " + i + " but get_n_F_functions() = " + get_n_F_functions());
		}

		switch(i) {
			case 0:
				return 1.;
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double[] evaluateDiffF(int i, double[] p) {

		if( (i<0) || (i>get_n_F_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffF()," +
				" i = " + i + " but get_n_F_functions() = " + get_n_F_functions());
		}

		switch(i) {
			case 0:
				double[] diff = {0., 0.};
				return  diff;
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double evaluateL(int i, int Q_l, double[] p) {

		if( (i<0) || (i>get_n_outputs()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateL()," +
				" i = " + i + " but get_n_outputs() = " + get_n_outputs());
		}

		if( (Q_l<0) || (Q_l>get_Q_l(i)-1) ) {
		  throw new RuntimeException(
		    "Input parameter is invalid in evaluateL()," +
				" Q_l = " + Q_l + " but get_Q_l(i) = " + get_Q_l(i));
		}

		return 1.;
	}

	public double[] evaluateDiffL(int i, int Q_l, double[] p) {

		if( (i<0) || (i>get_n_outputs()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffL()," +
				" i = " + i + " but get_n_outputs() = " + get_n_outputs());
		}

		if( (Q_l<0) || (Q_l>get_Q_l(i)-1) ) {
		  throw new RuntimeException(
		    "Input parameter is invalid in evaluateL()," +
				" Q_l = " + Q_l + " but get_Q_l(i) = " + get_Q_l(i));
		}

		double[] diff = {0., 0.};
		return  diff;
	}

  public double get_SCM_LB(double[] p) {
		return Math.min(p[0],p[1]);
	}

	public int get_n_M_functions() {
		return 1;
	}
  
	public int get_n_A_functions() {
		return 3;
	}

	public int get_n_F_functions() {
		return 1;
	}

	public int get_n_outputs() {
		return 1;
	}

	public int get_Q_l(int output_index) {
	  return 1;
	}

}

