
public class AffineFunctions {

	public double evaluateM(int i, double[] p) {
		return 1.;
	}

	public double evaluateA(int i, double[] p) {
		return 1.;
	}
	
	public double evaluateF(int i, double[] p) {
		return p[0];
	}
	
	public double evaluateL(int i, int Q_l, double[] p) {
		return 1.;
	}

	public double evaluateC(double[] p) {		
		return 1.;
	}

	public int get_n_M_functions() {
		return 1;
	}

	public int get_n_A_functions() {
		return 1;
	}
	
	public int get_n_F_functions() {
		return 1;
	}
	
	public int get_n_outputs() {
		return 2;
	}

	public int get_Q_l(int output_index){
		return 1;
	}

	public double get_SCM_LB(double[] p) {
		return 1.;
	}

}

