
public class AffineFunctions {

	public double evaluateM(int i, double[] p) {

		if( (i<0) || (i>get_n_M_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateM()," +
				" i = " + i + " but get_n_M_functions() = " + get_n_M_functions());
		}

		return 1.;
	}
	
	public double evaluateA(int i, double[] p) {

		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}

		switch(i) {
			case 0:
				return  0.05;
			case 1:
				return  p[0];
			case 2:
				return  p[1];
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

  public double get_SCM_LB(double[] p) {
		return 0.05;
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
		return 4;
	}

	public int get_Q_l(int output_index) {
	  return 1;
	}

}

