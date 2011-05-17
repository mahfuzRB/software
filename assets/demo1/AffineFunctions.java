
public class AffineFunctions {

	public double evaluateA(int i, double[] p) {
		
		if( (i<0) || (i>get_n_A_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateA()," +
				" i = " + i + " but get_n_A_functions() = " + get_n_A_functions());
		}
		
		switch(i) {
			case 0:
				return  -1000.0/273.0*p[0]/(-5.0+p[1]);
			case 1:
				return  -200.0/91.0*(-1.0+p[0])/p[1];
			case 2:
				return  -3.0/26.0*(-5.0+p[1])/p[0];
			case 3:
				return  -5.0/26.0*p[1]/(-1.0+p[0]);
			case 4:
				return  10.0/13.0*(-5.0+p[1])/(-1.0+p[0]);
			case 5:
				return  1000.0/637.0*(-1.0+p[0])/(-5.0+p[1]);
			case 6:
				return  15.0/182.0*p[1]/p[0];
			case 7:
				return  200.0/39.0*p[0]/p[1];
			case 8:
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
				return  10.0/3.0*Math.abs(p[0]);
			case 1:
				return  1.0;
			case 2:
				return  p[0];			
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
				return  10.0/3.0*Math.abs(p[0]);
			case 1:
				return  10.0/7.0-10.0/7.0*p[0];			
			default:
				throw new Error("Should not reach here");
		}
	}

	public int get_n_A_functions() {
		return 9;
	}
	
	public int get_n_F_functions() {
		return 3;
	}
	
    public int get_n_outputs(){
        return 1;
    }
    
    public int get_Q_l(int output_index){
        return 2;
    }
    
	public float[][] get_local_transformation(double[] p){
            float[][] LTfunc = new float[16][12];
            for (int i = 0; i < 16; i++){
                  LTfunc[0][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[0][1] = (float)( 0.0);
                  LTfunc[0][2] = (float)(0);
                  LTfunc[0][3] = (float)( 0.0);
                  LTfunc[0][4] = (float)( p[1]/4.0);
                  LTfunc[0][5] = (float)(0);
                  LTfunc[0][6] = (float)(0);
                  LTfunc[0][7] = (float)(0);
                  LTfunc[0][8] = (float)(0);
                  LTfunc[0][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[0][10] = (float)( 0.0);
                  LTfunc[0][11] = (float)(0);
                  LTfunc[1][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[1][1] = (float)( 0.0);
                  LTfunc[1][2] = (float)(0);
                  LTfunc[1][3] = (float)( 0.0);
                  LTfunc[1][4] = (float)( p[1]/4.0);
                  LTfunc[1][5] = (float)(0);
                  LTfunc[1][6] = (float)(0);
                  LTfunc[1][7] = (float)(0);
                  LTfunc[1][8] = (float)(0);
                  LTfunc[1][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[1][10] = (float)( 0.0);
                  LTfunc[1][11] = (float)(0);
                  LTfunc[2][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[2][1] = (float)( 0.0);
                  LTfunc[2][2] = (float)(0);
                  LTfunc[2][3] = (float)( 0.0);
                  LTfunc[2][4] = (float)( p[1]/4.0);
                  LTfunc[2][5] = (float)(0);
                  LTfunc[2][6] = (float)(0);
                  LTfunc[2][7] = (float)(0);
                  LTfunc[2][8] = (float)(0);
                  LTfunc[2][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[2][10] = (float)( 0.0);
                  LTfunc[2][11] = (float)(0);
                  LTfunc[3][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[3][1] = (float)( 0.0);
                  LTfunc[3][2] = (float)(0);
                  LTfunc[3][3] = (float)( 0.0);
                  LTfunc[3][4] = (float)( p[1]/4.0);
                  LTfunc[3][5] = (float)(0);
                  LTfunc[3][6] = (float)(0);
                  LTfunc[3][7] = (float)(0);
                  LTfunc[3][8] = (float)(0);
                  LTfunc[3][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[3][10] = (float)( 0.0);
                  LTfunc[3][11] = (float)(0);
                  LTfunc[4][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[4][1] = (float)( 0.0);
                  LTfunc[4][2] = (float)(0);
                  LTfunc[4][3] = (float)( 0.0);
                  LTfunc[4][4] = (float)( 5.0-p[1]);
                  LTfunc[4][5] = (float)(0);
                  LTfunc[4][6] = (float)(0);
                  LTfunc[4][7] = (float)(0);
                  LTfunc[4][8] = (float)(0);
                  LTfunc[4][9] = (float)( 0.0);
                  LTfunc[4][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[4][11] = (float)(0);
                  LTfunc[5][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[5][1] = (float)( 0.0);
                  LTfunc[5][2] = (float)(0);
                  LTfunc[5][3] = (float)( 0.0);
                  LTfunc[5][4] = (float)( 5.0-p[1]);
                  LTfunc[5][5] = (float)(0);
                  LTfunc[5][6] = (float)(0);
                  LTfunc[5][7] = (float)(0);
                  LTfunc[5][8] = (float)(0);
                  LTfunc[5][9] = (float)( 0.0);
                  LTfunc[5][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[5][11] = (float)(0);
                  LTfunc[6][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[6][1] = (float)( 0.0);
                  LTfunc[6][2] = (float)(0);
                  LTfunc[6][3] = (float)( 0.0);
                  LTfunc[6][4] = (float)( 5.0-p[1]);
                  LTfunc[6][5] = (float)(0);
                  LTfunc[6][6] = (float)(0);
                  LTfunc[6][7] = (float)(0);
                  LTfunc[6][8] = (float)(0);
                  LTfunc[6][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[6][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[6][11] = (float)(0);
                  LTfunc[7][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[7][1] = (float)( 0.0);
                  LTfunc[7][2] = (float)(0);
                  LTfunc[7][3] = (float)( 0.0);
                  LTfunc[7][4] = (float)( 5.0-p[1]);
                  LTfunc[7][5] = (float)(0);
                  LTfunc[7][6] = (float)(0);
                  LTfunc[7][7] = (float)(0);
                  LTfunc[7][8] = (float)(0);
                  LTfunc[7][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[7][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[7][11] = (float)(0);
                  LTfunc[8][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[8][1] = (float)( 0.0);
                  LTfunc[8][2] = (float)(0);
                  LTfunc[8][3] = (float)( 0.0);
                  LTfunc[8][4] = (float)( p[1]/4.0);
                  LTfunc[8][5] = (float)(0);
                  LTfunc[8][6] = (float)(0);
                  LTfunc[8][7] = (float)(0);
                  LTfunc[8][8] = (float)(0);
                  LTfunc[8][9] = (float)( 0.0);
                  LTfunc[8][10] = (float)( 0.0);
                  LTfunc[8][11] = (float)(0);
                  LTfunc[9][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[9][1] = (float)( 0.0);
                  LTfunc[9][2] = (float)(0);
                  LTfunc[9][3] = (float)( 0.0);
                  LTfunc[9][4] = (float)( p[1]/4.0);
                  LTfunc[9][5] = (float)(0);
                  LTfunc[9][6] = (float)(0);
                  LTfunc[9][7] = (float)(0);
                  LTfunc[9][8] = (float)(0);
                  LTfunc[9][9] = (float)( 0.0);
                  LTfunc[9][10] = (float)( 0.0);
                  LTfunc[9][11] = (float)(0);
                  LTfunc[10][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[10][1] = (float)( 0.0);
                  LTfunc[10][2] = (float)(0);
                  LTfunc[10][3] = (float)( 0.0);
                  LTfunc[10][4] = (float)( 5.0-p[1]);
                  LTfunc[10][5] = (float)(0);
                  LTfunc[10][6] = (float)(0);
                  LTfunc[10][7] = (float)(0);
                  LTfunc[10][8] = (float)(0);
                  LTfunc[10][9] = (float)( -3.0/7.0+10.0/7.0*p[0]);
                  LTfunc[10][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[10][11] = (float)(0);
                  LTfunc[11][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[11][1] = (float)( 0.0);
                  LTfunc[11][2] = (float)(0);
                  LTfunc[11][3] = (float)( 0.0);
                  LTfunc[11][4] = (float)( 5.0-p[1]);
                  LTfunc[11][5] = (float)(0);
                  LTfunc[11][6] = (float)(0);
                  LTfunc[11][7] = (float)(0);
                  LTfunc[11][8] = (float)(0);
                  LTfunc[11][9] = (float)( -3.0/7.0+10.0/7.0*p[0]);
                  LTfunc[11][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[11][11] = (float)(0);
                  LTfunc[12][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[12][1] = (float)( 0.0);
                  LTfunc[12][2] = (float)(0);
                  LTfunc[12][3] = (float)( 0.0);
                  LTfunc[12][4] = (float)( 5.0-p[1]);
                  LTfunc[12][5] = (float)(0);
                  LTfunc[12][6] = (float)(0);
                  LTfunc[12][7] = (float)(0);
                  LTfunc[12][8] = (float)(0);
                  LTfunc[12][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[12][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[12][11] = (float)(0);
                  LTfunc[13][0] = (float)( 10.0/3.0*p[0]);
                  LTfunc[13][1] = (float)( 0.0);
                  LTfunc[13][2] = (float)(0);
                  LTfunc[13][3] = (float)( 0.0);
                  LTfunc[13][4] = (float)( 5.0-p[1]);
                  LTfunc[13][5] = (float)(0);
                  LTfunc[13][6] = (float)(0);
                  LTfunc[13][7] = (float)(0);
                  LTfunc[13][8] = (float)(0);
                  LTfunc[13][9] = (float)( 1.0/2.0-5.0/3.0*p[0]);
                  LTfunc[13][10] = (float)( -20.0+5.0*p[1]);
                  LTfunc[13][11] = (float)(0);
                  LTfunc[14][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[14][1] = (float)( 0.0);
                  LTfunc[14][2] = (float)(0);
                  LTfunc[14][3] = (float)( 0.0);
                  LTfunc[14][4] = (float)( p[1]/4.0);
                  LTfunc[14][5] = (float)(0);
                  LTfunc[14][6] = (float)(0);
                  LTfunc[14][7] = (float)(0);
                  LTfunc[14][8] = (float)(0);
                  LTfunc[14][9] = (float)( -3.0/7.0+10.0/7.0*p[0]);
                  LTfunc[14][10] = (float)( 0.0);
                  LTfunc[14][11] = (float)(0);
                  LTfunc[15][0] = (float)( 10.0/7.0-10.0/7.0*p[0]);
                  LTfunc[15][1] = (float)( 0.0);
                  LTfunc[15][2] = (float)(0);
                  LTfunc[15][3] = (float)( 0.0);
                  LTfunc[15][4] = (float)( p[1]/4.0);
                  LTfunc[15][5] = (float)(0);
                  LTfunc[15][6] = (float)(0);
                  LTfunc[15][7] = (float)(0);
                  LTfunc[15][8] = (float)(0);
                  LTfunc[15][9] = (float)( -3.0/7.0+10.0/7.0*p[0]);
                  LTfunc[15][10] = (float)( 0.0);
                  LTfunc[15][11] = (float)(0);
            }			
		return LTfunc;
    }

}
