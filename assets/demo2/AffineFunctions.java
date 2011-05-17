
public class AffineFunctions {

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
				return  1.0;			
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
				return  1.0;			
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double evaluateL(int i, int Q_l, double[] p) {
		
		if( (Q_l<0) || (Q_l>get_Q_l(i)-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateL()," +
				" i = " + i + " but get_Q_l() = " + get_Q_l(i));
		}
		
		switch(Q_l) {
			case 0:
				return  1.0;			
			default:
				throw new Error("Should not reach here");
		}
	}

    public double get_SCM_LB(double[] p) {
		return 1.f;
	}

    
	public int get_n_A_functions() {
		return 2;
	}
	
	public int get_n_F_functions() {
		return 1;
	}
	
    public int get_n_outputs(){
        return 1;
    }
    
    public int get_Q_l(int output_index){
        return 1;
    }
    
	public float[][] get_local_transformation(double[] p){
            float[][] LTfunc = new float[1][12];
            for (int i = 0; i < 1; i++){
                  LTfunc[0][0] = (float)( 1.0);
                  LTfunc[0][1] = (float)( 0.0);
                  LTfunc[0][2] = (float)(0);
                  LTfunc[0][3] = (float)( 0.0);
                  LTfunc[0][4] = (float)( 1.0);
                  LTfunc[0][5] = (float)(0);
                  LTfunc[0][6] = (float)(0);
                  LTfunc[0][7] = (float)(0);
                  LTfunc[0][8] = (float)(1);
                  LTfunc[0][9] = (float)( 0.0);
                  LTfunc[0][10] = (float)( 0.0);
                  LTfunc[0][11] = (float)(0);
            }			
		return LTfunc;
    }

}
