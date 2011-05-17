
public class AffineFunctions {

	public double evaluateA(int i, double[] p, boolean is_real) {
		
		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}
		
		switch(i) {
			case 0:
				return (is_real) ? ( 1.75*(4.0/1.75E2)+(9.1E1/5.0E1)/1.75-2.0/2.5E1) : ( 0.0);
			case 1:
				return (is_real) ? ( 1.75*(p[0]*p[0])*(-4.0/7.0)) : ( 0.0);
			case 2:
				return (is_real) ? ( -p[0]*p[0]) : ( 0.0);
			case 3:
				return (is_real) ? ( (7.0/4.0)/1.75) : ( 0.0);
			case 4:
				return (is_real) ? ( 1.0) : ( 0.0);
			case 5:
				return (is_real) ? ( 1.75) : ( 0.0);
			case 6:
				return (is_real) ? ( 1.75*1.75) : ( 0.0);
			case 7:
				return (is_real) ? ( 0.0) : ( p[0]);			
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double evaluateF(int i, double[] p, boolean is_real) {
		
		if( (i<0) || (i>get_n_F_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateF()," +
				" i = " + i + " but get_n_F_functions() = " + get_n_F_functions());
		}
		
		switch(i) {
			case 0:
				return (is_real) ? ( 0.0) : ( p[0]*2.0);			
			default:
				throw new Error("Should not reach here");
		}
	}
	
	public double evaluateL(int i, int Q_l, double[] p, boolean is_real) {
		
		if( (Q_l<0) || (Q_l>get_Q_l(i)-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateL()," +
				" i = " + i + " but get_Q_l() = " + get_Q_l(i));
		}
		
		switch(Q_l) {
			case 0:
				return (is_real) ? ( 1.0) : ( 0.0);			
			default:
				throw new Error("Should not reach here");
		}
	}

    public double get_SCM_LB(double[] p) { 
  return 5e-2; 
	}
    
	public int get_n_A_functions() {
		return 8;
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
    
    public int get_n_uL_functions() {
		return 0;
	}
    
	public float[][] get_local_transformation(double[] p){
            float[][] LTfunc = new float[8][12];
            for (int i = 0; i < 8; i++){
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
                  LTfunc[1][0] = (float)( 1.0);
                  LTfunc[1][1] = (float)( 0.0);
                  LTfunc[1][2] = (float)(0);
                  LTfunc[1][3] = (float)( 0.0);
                  LTfunc[1][4] = (float)( 1.0);
                  LTfunc[1][5] = (float)(0);
                  LTfunc[1][6] = (float)(0);
                  LTfunc[1][7] = (float)(0);
                  LTfunc[1][8] = (float)(1);
                  LTfunc[1][9] = (float)( 0.0);
                  LTfunc[1][10] = (float)( 0.0);
                  LTfunc[1][11] = (float)(0);
                  LTfunc[2][0] = (float)( 1.0);
                  LTfunc[2][1] = (float)( 0.0);
                  LTfunc[2][2] = (float)(0);
                  LTfunc[2][3] = (float)( 1.75*(4.0/3.5E1)-1.0/5.0);
                  LTfunc[2][4] = (float)( 1.75*(4.0/7.0));
                  LTfunc[2][5] = (float)(0);
                  LTfunc[2][6] = (float)(0);
                  LTfunc[2][7] = (float)(0);
                  LTfunc[2][8] = (float)(1);
                  LTfunc[2][9] = (float)( 0.0);
                  LTfunc[2][10] = (float)( 1.75*(-6.0/7.0)+3.0/2.0);
                  LTfunc[2][11] = (float)(0);
                  LTfunc[3][0] = (float)( 1.0);
                  LTfunc[3][1] = (float)( 0.0);
                  LTfunc[3][2] = (float)(0);
                  LTfunc[3][3] = (float)( 1.75*(-4.0/3.5E1)+1.0/5.0);
                  LTfunc[3][4] = (float)( 1.75*(4.0/7.0));
                  LTfunc[3][5] = (float)(0);
                  LTfunc[3][6] = (float)(0);
                  LTfunc[3][7] = (float)(0);
                  LTfunc[3][8] = (float)(1);
                  LTfunc[3][9] = (float)( 0.0);
                  LTfunc[3][10] = (float)( 1.75*(6.0/7.0)-3.0/2.0);
                  LTfunc[3][11] = (float)(0);
                  LTfunc[4][0] = (float)( 1.0);
                  LTfunc[4][1] = (float)( 0.0);
                  LTfunc[4][2] = (float)(0);
                  LTfunc[4][3] = (float)( 0.0);
                  LTfunc[4][4] = (float)( 1.75*(4.0/7.0));
                  LTfunc[4][5] = (float)(0);
                  LTfunc[4][6] = (float)(0);
                  LTfunc[4][7] = (float)(0);
                  LTfunc[4][8] = (float)(1);
                  LTfunc[4][9] = (float)( 0.0);
                  LTfunc[4][10] = (float)( 0.0);
                  LTfunc[4][11] = (float)(0);
                  LTfunc[5][0] = (float)( 1.0);
                  LTfunc[5][1] = (float)( 0.0);
                  LTfunc[5][2] = (float)(0);
                  LTfunc[5][3] = (float)( 0.0);
                  LTfunc[5][4] = (float)( 1.75*(4.0/7.0));
                  LTfunc[5][5] = (float)(0);
                  LTfunc[5][6] = (float)(0);
                  LTfunc[5][7] = (float)(0);
                  LTfunc[5][8] = (float)(1);
                  LTfunc[5][9] = (float)( 0.0);
                  LTfunc[5][10] = (float)( 0.0);
                  LTfunc[5][11] = (float)(0);
                  LTfunc[6][0] = (float)( 1.0);
                  LTfunc[6][1] = (float)( 0.0);
                  LTfunc[6][2] = (float)(0);
                  LTfunc[6][3] = (float)( 1.75*(-2.0/5.0)+7.0/1.0E1);
                  LTfunc[6][4] = (float)( 1.0);
                  LTfunc[6][5] = (float)(0);
                  LTfunc[6][6] = (float)(0);
                  LTfunc[6][7] = (float)(0);
                  LTfunc[6][8] = (float)(1);
                  LTfunc[6][9] = (float)( 0.0);
                  LTfunc[6][10] = (float)( 1.75*4.0-7.0);
                  LTfunc[6][11] = (float)(0);
                  LTfunc[7][0] = (float)( 1.0);
                  LTfunc[7][1] = (float)( 0.0);
                  LTfunc[7][2] = (float)(0);
                  LTfunc[7][3] = (float)( 1.75*(2.0/5.0)-7.0/1.0E1);
                  LTfunc[7][4] = (float)( 1.0);
                  LTfunc[7][5] = (float)(0);
                  LTfunc[7][6] = (float)(0);
                  LTfunc[7][7] = (float)(0);
                  LTfunc[7][8] = (float)(1);
                  LTfunc[7][9] = (float)( 0.0);
                  LTfunc[7][10] = (float)( 1.75*(-4.0)+7.0);
                  LTfunc[7][11] = (float)(0);
            }			
		return LTfunc;
    }
    
    public boolean is_derived_output(){
        return true;
    }
    
    public double[] cal_derived_output(double[] input){
        double s = input[0]*input[0] + input[1]*input[1] - 2*input[0] + input[2]*input[2] + input[3]*input[3] + 1;
        double ds = 2*input[2]*(Math.signum(input[0])*input[0] + Math.signum(input[1])*input[1]) + 2*input[2];
        
        double s1 = Math.sqrt(Math.abs(s + ds));
        double s2 = Math.sqrt(Math.abs(s - ds));
        
        double[] output = new double[4];
        output[0] = 0.5f*(s1+s2);
        output[1] = 0;
        output[2] = 0.5f*(s1-s2);
        output[3] = 0;
        
        return output;
    }

}
