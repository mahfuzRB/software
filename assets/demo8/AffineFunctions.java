
public class AffineFunctions {

	public double evaluateM(int i, double[] p) {
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
				return  1.;
			case 1:
				return  p[0];
			case 2:
				return  p[1];
			case 3:
				return  p[2];
			case 4:
				return  p[3];
			case 5:
				return  p[4];
			case 6:
				return  p[5];
			case 7:
				return  p[6];
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
	
	public double evaluateL(int i, int q_l, double[] p) {
		
		if( (i<0) || (i>get_n_outputs()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateL()," +
				" i = " + i + " but get_n_outputs() = " + get_n_outputs());
		}
		
		if( (q_l<0) || (q_l>get_Q_l(i)-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateL()," +
				" q_l = " + q_l + " but get_Q_l(i) = " + get_Q_l(i));
		}
		
		return 1.;
	}

	public double get_SCM_LB(double[] p) {
		double p_min = p[0];

		for(int i=1; i<p.length; i++) {
			if( p[i] < p_min ) {
				p_min = p[i];
			}
		}

		return p_min;
	}

	public int get_n_M_functions() {
		return 1;
	}
	
	public int get_n_A_functions() {
		return 8;
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

