
public class AffineFunctions {

	public double evaluateA(int i, double[] p, boolean is_real) {
		
		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}
		
		switch(i) {
			case 0:
				return (is_real) ? ( 1.0) : ( 0.0);
			case 1:
				return (is_real) ? ( -p[2]*p[2]) : ( 0.0);
			case 2:
				return (is_real) ? ( 0.0) : ( p[2]);
			case 3:
				return (is_real) ? ( (p[1]*p[2])/(p[0]*p[0]+p[1]*p[1])) : ( (p[0]*p[2])/(p[0]*p[0]+p[1]*p[1]));			
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
				return (is_real) ? ( 0.0) : ( -p[2]);			
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

    public double[] evaluateDiffA(int i, double[] p, boolean is_real) {

		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}
		double[] diff = new double[3];
		switch (i) {
			case 0:
			    diff[0] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[1] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[2] = (is_real) ? ( 0.0) : ( 0.0);
			    return diff;
			case 1:
			    diff[0] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[1] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[2] = (is_real) ? ( p[2]*(-2.0)) : ( 0.0);
			    return diff;
			case 2:
			    diff[0] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[1] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[2] = (is_real) ? ( 0.0) : ( 1.0);
			    return diff;
			case 3:
			    diff[0] = (is_real) ? (-p[0]*p[1]*p[2]*1/Math.pow(p[0]*p[0]+p[1]*p[1],2.0)*2.0) : (-p[2]*(p[0]*p[0]-p[1]*p[1])*1/Math.pow(p[0]*p[0]+p[1]*p[1],2.0));
			    diff[1] = (is_real) ? ( p[2]*(p[0]*p[0]-p[1]*p[1])*1/Math.pow(p[0]*p[0]+p[1]*p[1],2.0)) : (-p[0]*p[1]*p[2]*1/Math.pow(p[0]*p[0]+p[1]*p[1],2.0)*2.0);
			    diff[2] = (is_real) ? ( p[1]/(p[0]*p[0]+p[1]*p[1])) : ( p[0]/(p[0]*p[0]+p[1]*p[1]));
			    return diff;
		}
		return diff;
	}
    
    public double[] evaluateDiffF(int i, double[] p, boolean is_real) {

		if( (i<0) || (i>get_n_F_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateDiffF()," +
				" i = " + i + " but get_n_F_functions() = " + get_n_F_functions());
		}
		double[] diff = new double[3];
		switch (i) {
			case 0:
			    diff[0] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[1] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[2] = (is_real) ? ( 0.0) : ( -1.0);
			    return diff;
		}
		return diff;
	}	
	
	public double[] evaluateDiffL(int i, int Q_l, double[] p, boolean is_real) {

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
		double[] diff = new double[3];
		switch (i) {
			case 0:
			    diff[0] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[1] = (is_real) ? ( 0.0) : ( 0.0);
			    diff[2] = (is_real) ? ( 0.0) : ( 0.0);
			    return diff;
		}
		return diff;
	}
    
    public double get_SCM_LB(double[] p) { 
  return 5e-3; 
	}
    
	public int get_n_A_functions() {
		return 4;
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
    
    public boolean is_derived_output(){
        return false;
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
