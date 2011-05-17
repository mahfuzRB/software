
public class AffineFunctions {

	public double[][] evaluateA_array(double[] p) {
		
        double[][] theta_a = new double[37][2];
        
        for (int i = 0; i<37;i++){
            theta_a[i][0] = 0.0;
            theta_a[i][1] = 0.0;
        }
        
        theta_a[0][0] = 1.; theta_a[0][1] = 0.;
        theta_a[1][0] = 0.; theta_a[1][1] = p[1];
     
        int[] M = {9,7,7,9,3};
        
        double[][][] B = {{{ 1.0000,         0,         0,         0,         0,         0,         0,         0,         0},
                           { 0.4590,    1.0000,         0,         0,         0,         0,         0,         0,         0},
                           { 0.6726,    0.3165,    1.0000,         0,         0,    0.0000,         0,         0,         0},
                           { 0.5372,    0.8292,    0.0808,    1.0000,    0.0000,   -0.0000,         0,   -0.0000,    0.0000},
                           { 0.5544,    0.6873,    0.4653,    0.8065,    1.0000,    0.0000,         0,         0,   -0.0000},
                           { 0.7360,    0.3551,    0.4542,    0.5857,   -0.2405,    1.0000,         0,         0,         0},
                           { 0.6325,    0.6640,    0.0568,   -0.3469,    0.4672,   -0.6491,    1.0000,   -0.0000,         0},
                           { 0.5069,    0.9280,   -0.0660,    0.7929,   -0.4710,    0.2684,    0.4111,    1.0000,   -0.0000},
                           { 0.5496,    0.7459,    0.2916,    0.9286,    0.5826,   -0.0127,   -0.1861,   -0.5628,    1.0000}},
                          {{ 1.0000,         0,         0,         0,         0,         0,         0},
                           {-0.5053,    1.0000,    0.0000,         0,         0,    0.0000,         0},
                           {-0.7224,    0.2812,    1.0000,         0,         0,         0,         0},
                           {-0.6226,   -0.7142,   -0.3301,    1.0000,    0.0000,   -0.0000,         0},
                           { 0.4063,    0.3991,    0.2123,   -0.2277,    1.0000,         0,    0.0000},
                           {-0.4297,    0.6803,   -0.0123,   -0.1558,    0.9075,    1.0000,         0},
                           { 0.6924,    0.1602,   -0.1638,   -0.3214,   -0.1203,    0.3026,    1.0000}},
                          {{ 1.0000,         0,         0,         0,         0,         0,         0},
                           {-0.5053,    1.0000,    0.0000,         0,         0,    0.0000,         0},
                           {-0.7224,    0.2812,    1.0000,         0,         0,         0,         0},
                           {-0.6226,   -0.7142,   -0.3301,    1.0000,    0.0000,   -0.0000,         0},
                           { 0.4063,    0.3991,    0.2123,   -0.2277,    1.0000,         0,    0.0000},
                           {-0.4297,    0.6803,   -0.0123,   -0.1558,    0.9075,    1.0000,         0},
                           { 0.6924,    0.1602,   -0.1638,   -0.3214,   -0.1203,    0.3026,    1.0000}},
                          {{ 1.0000,         0,         0,         0,         0,         0,         0,         0,         0},
                           { 0.4523,    1.0000,         0,         0,         0,         0,         0,    0.0000,    0.0000},
                           { 0.5933,    0.4621,    1.0000,         0,         0,         0,         0,         0,         0},
                           { 0.8413,    0.3711,   -0.3006,    1.0000,         0,         0,         0,   -0.0000,         0},
                           { 0.5637,    0.7162,    0.2954,   -0.6951,    1.0000,         0,         0,         0,   -0.0000},
                           { 0.7586,    0.2422,    0.7021,    0.5239,   -0.0294,    1.0000,         0,    0.0000,         0},
                           { 0.7792,    0.4204,   -0.0650,    0.3166,    0.1683,    0.4210,    1.0000,   -0.0000,         0},
                           { 0.6480,    0.7869,   -0.5141,    0.0403,    0.3448,   -0.3798,    0.5584,    1.0000,         0},
                           { 0.5034,    0.8000,    0.3846,   -0.3438,    0.7693,    0.6718,    0.3794,   -0.5221,    1.0000}},
                          {{ 1.0000,         0,         0},
                           { 0.6668,    1.0000,         0},
                           { 0.9840,    0.0481,    1.0000}}};

        double[][][] x = {{{ -0.4696,    0.0000},
                           { -0.6976,    0.7122},
                           { -0.3522,   -0.3512},
                           { -0.4797,   -0.4756},
                           { -0.3563,   -0.3554},
                           { -0.4346,    0.2483},
                           { -0.3475,    0.4494},
                           { -0.5849,    0.5812},
                           { -0.4020,   -0.4011}},
                          {{ -0.0915,   -0.9925},
                           { -0.3563,   -0.3554},
                           { -0.3522,   -0.3512},
                           { -0.3517,    0.3563},
                           { -0.4728,   -0.4891},
                           { -0.4387,   -0.4367},
                           { -0.1371,   -0.4823}},
                          {{ -0.0915,   -0.9925},
                           { -0.3563,   -0.3554},
                           { -0.3522,   -0.3512},
                           { -0.3517,    0.3563},
                           { -0.4728,   -0.4891},
                           { -0.4387,   -0.4367},
                           { -0.1371,   -0.4823}},
                          {{ -0.6976,    0.7122},
                           { -0.4696,    0.0000},
                           { -0.1491,   -0.4785},
                           { -0.3563,   -0.3554},
                           {  0.2538,    0.9639},
                           { -0.3517,    0.3563},
                           { -0.4434,   -0.3494},
                           { -0.3534,    0.2585},
                           { -0.6421,   -0.0836}},
                          {{ -0.2987,   -0.2217},
                           { -0.3126,    0.3909},
                           { -0.5303,    0.0435}}};
                                   
        int icount = 2;
        for (int i = 0; i<5; i++){
            double[] gvec = new double[M[i]];
            for (int j = 0; j<M[i]; j++){
                double[] fval = eval_func(p, x[i][j]);
                gvec[j] = fval[i];
            }
            double[] theta_mp = blackslash(B[i],gvec);
            for (int j = 0; j<M[i]; j++){
                if (i==4)
                    theta_a[icount][0] = -p[1]*p[1]*theta_mp[j];
                else
                    theta_a[icount][0] = theta_mp[j];
                icount++;
            }                
        }
                           
        return theta_a;
	}
	
	public double evaluateF(int i, double[] p, boolean is_real) {
		
		if( (i<0) || (i>get_n_F_functions()-1) ) {
			throw new RuntimeException(
				"Input parameter is invalid in evaluateF()," +
				" i = " + i + " but get_n_F_functions() = " + get_n_F_functions());
		}
		
		switch(i) {
			case 0:
				return (is_real) ? ( 1.0) : ( 0.0);
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
  return 0.0197; 
	}
    
	public int get_n_A_functions() {
		return 37;
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
            float[][] LTfunc = new float[4][12];
            for (int i = 0; i < 4; i++){
                  LTfunc[i][0] = (float)( 1.0);
                  LTfunc[i][1] = (float)( 0.0);
                  LTfunc[i][2] = (float)(0);
                  LTfunc[i][3] = (float)( 0.0);
                  LTfunc[i][4] = (float)( 1.0);
                  LTfunc[i][5] = (float)(0);
                  LTfunc[i][6] = (float)(0);
                  LTfunc[i][7] = (float)(0);
                  LTfunc[i][8] = (float)(1);
                  LTfunc[i][9] = (float)( 0.0);
                  LTfunc[i][10] = (float)( 0.0);
                  LTfunc[i][11] = (float)(0);
            }			
		return LTfunc;
    }
    
    public boolean is_derived_output(){
        return false;
    }

    public double[][] matrix22_transpose(double[][] a){
        double[][] b = new double[2][2];
        b[0][0] = a[0][0]; b[0][1] = a[1][0]; b[1][0] = a[0][1]; b[1][1] = a[1][1];
        return b;
    }
    
    public double[][] matrix22_mult(double[][] a, double[][] b){
        double[][] c = new double[2][2];
        c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
        c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];
        c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
        c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
        return c;
    }
    
    public double matrix22_det(double[][] a){
        return a[0][0]*a[1][1] - a[0][1]*a[1][0];
    }
        
    public double[][] matrix22_inv(double[][] a){
        double[][] b = new double[2][2];
        double deta_inv = 1./matrix22_det(a);
        b[0][0] = deta_inv*( a[1][1]);
        b[0][1] = deta_inv*(-a[0][1]);
        b[1][0] = deta_inv*(-a[1][0]);
        b[1][1] = deta_inv*( a[0][0]);
        return b;
    }
        
    public double[] eval_func(double[] p, double[] x){
        double pi = Math.PI;
        double mu_r = pi/4;
        double R2 = 1.0; 
        double R1 = 0.5;
        
        double mu = p[0];
        
        double[] f = new double[5];
        
        double r = Math.sqrt(x[0]*x[0]+x[1]*x[1]) + 0.0001;
        double theta = Math.atan2(x[1],x[0]);
        
        double[][] J = {{1,0},{0,1}};
        int ri = 0;
        if (r>R2){ri = 1;}
        else { if (r<R1){
            ri = 2;
                   theta = Math.atan2(x[1],-x[0]);
                   double[][] J1 = {{-Math.cos(theta), Math.sin(theta)}, {r*Math.sin(theta), r*Math.cos(theta)}};
                   double[][] J2 = {{-Math.cos((mu*theta)/mu_r),             Math.sin((mu*theta)/mu_r)},
                                    {(mu*r*Math.sin((mu*theta)/mu_r))/mu_r, (mu*r*Math.cos((mu*theta)/mu_r))/mu_r}};
                   J = matrix22_mult(matrix22_inv(J1),J2);}
               else { if (Math.abs(theta)<pi-mu_r){
                   ri = 3;
                          theta = Math.atan2(x[1],x[0]);
                          double[][] J1 = {{Math.cos(theta), Math.sin(theta)},{-r*Math.sin(theta), r*Math.cos(theta)}};
                          double[][] J2 = {{ Math.cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) + r*theta*Math.sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r))), Math.sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) - r*theta*Math.cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r)))},
                                           {-r*Math.sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))),                                                                      r*Math.cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))}};
                          J = matrix22_mult(matrix22_inv(J1),J2);}
                      else {
                   ri = 4;
                          theta = Math.atan2(x[1],-x[0]);
                          double[][] J1 = {{-Math.cos(theta), Math.sin(theta)},{r*Math.sin(theta), r*Math.cos(theta)}};
                          double[][] J2 = {{-Math.cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*Math.sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2))), Math.sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*Math.cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2)))},
                                           { r*Math.sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))),                                                        r*Math.cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))}};
                          J = matrix22_mult(matrix22_inv(J1),J2);}
                    }
             }
        double detJ = matrix22_det(J);
        double[][] invJ = matrix22_inv(J);
        double[][] T = matrix22_mult(invJ,matrix22_transpose(invJ));
        f[0] = T[0][0]*detJ;
        f[1] = T[0][1]*detJ;
        f[2] = T[1][0]*detJ;
        f[3] = T[1][1]*detJ;
        f[4] = detJ;
        
        return f;
    }
    
    public double[] blackslash(double[][] A, double[] b){
        double[] x = new double[b.length];
        for (int i = 0; i<b.length; i++){
            double tmp = b[i];
            for (int j = 0; j<=i-1; j++)
                tmp-= A[i][j]*x[j];
            x[i] = tmp/A[i][i];
        }
        return x;
    }
    
    public boolean is_custom_mesh_transform(){
        return true;
    }
    
    public float[] mesh_transform(double[] p, float[] x){
        double pi = Math.PI;
        double mu_r = pi/4;
        double R2 = 1.0f; 
        double R1 = 0.5f;
        
        double mu = p[0];
        
        double[] krat = {(pi-mu)/(pi-mu_r),mu/mu_r};
        
        int node_num = x.length/3;
        for (int i = 0; i< node_num; i++){
            double nR = Math.sqrt(x[i*3]*x[i*3] + x[i*3+1]*x[i*3+1]);
            double theta = Math.atan2(x[i*3+1],x[i*3]);
            if (nR>R2){}
            else { if (nR<R1-0.0001){                   
                   double ntheta = Math.atan2(x[i*3+1],-x[i*3]);
                   double rtheta = krat[1];
                   x[i*3+0] = (float)(nR*(-Math.cos(rtheta*ntheta)));
                   x[i*3+1] = (float)(nR*( Math.sin(rtheta*ntheta)));}
               else { if (Math.abs(theta)<(pi-mu_r)){                   
                          double ntheta = Math.atan2(x[i*3+1],x[i*3]);
                          double Rrat = (R2-nR)/(R2-R1);
                          double rtheta = Rrat*krat[0] + (1-Rrat)*1;
                          x[i*3+0] = (float)(nR*(Math.cos(rtheta*ntheta)));
                          x[i*3+1] = (float)(nR*(Math.sin(rtheta*ntheta)));}
                      else {
                          double ntheta = Math.atan2(x[i*3+1],-x[i*3]);
                          double Rrat = (R2-nR)/(R2-R1);
                          double rtheta = Rrat*krat[1] + (1-Rrat)*1;
                          x[i*3+0] = (float)(nR*(-Math.cos(rtheta*ntheta)));
                          x[i*3+1] = (float)(nR*( Math.sin(rtheta*ntheta)));}                          
                    }
             }
        }
        return x;
    }
    
}

/*

for i = 1:size(x,1)    
    r = norm(x(i,:));
%     switch pgnum(i)
%         case 1
%             J = eye(2);
%         case 2
%             theta = atan2(x(i,2),-x(i,1));
%             J1 = [-cos(theta), sin(theta); r*sin(theta), r*cos(theta)];
%             J2 = [-cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2))), sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2)));...
%                    r*sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))),                                                        r*cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))];
%             J = J1\J2;
%         case 3       
%             theta = atan2(x(i,2),-x(i,1));            
%             J1 = [-cos(theta), sin(theta); r*sin(theta), r*cos(theta)];
%             J2 = [-cos((mu*theta)/mu_r),             sin((mu*theta)/mu_r);...
%                    (mu*r*sin((mu*theta)/mu_r))/mu_r, (mu*r*cos((mu*theta)/mu_r))/mu_r];
%             J = J1\J2;
%         case 4
%             theta = atan2(x(i,2),x(i,1));
%             J1 = [cos(theta), sin(theta);-r*sin(theta), r*cos(theta)];
%             J2 = [ cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) + r*theta*sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r))), sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) - r*theta*cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r)));...
%                   -r*sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))),                                                                      r*cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))];            
%             J = J1\J2;
%     end
    theta = atan2(x(i,2),x(i,1));
    if r > R2        
        J = eye(2);
    elseif (r<R1)
        theta = atan2(x(i,2),-x(i,1));            
        J1 = [-cos(theta), sin(theta); r*sin(theta), r*cos(theta)];
        J2 = [-cos((mu*theta)/mu_r),             sin((mu*theta)/mu_r);...
            (mu*r*sin((mu*theta)/mu_r))/mu_r, (mu*r*cos((mu*theta)/mu_r))/mu_r];
        J = J1\J2;
    elseif abs(theta)<pi-mu_r
        theta = atan2(x(i,2),x(i,1));
        J1 = [cos(theta), sin(theta);-r*sin(theta), r*cos(theta)];
        J2 = [ cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) + r*theta*sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r))), sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))) - r*theta*cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*(1/(R1 - R2) - (pi - mu)/((R1 - R2)*(pi - mu_r)));...
            -r*sin(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))),                                                                      r*cos(theta*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r))))*((R1 - r)/(R1 - R2) - ((pi - mu)*(R2 - r))/((R1 - R2)*(pi - mu_r)))];
        J = J1\J2;
    else
        theta = atan2(x(i,2),-x(i,1));
        J1 = [-cos(theta), sin(theta); r*sin(theta), r*cos(theta)];
        J2 = [-cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2))), sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))) - r*theta*cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*(1/(R1 - R2) - mu/(mu_r*(R1 - R2)));...
            r*sin(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))),                                                        r*cos(theta*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2))))*((R1 - r)/(R1 - R2) - (mu*(R2 - r))/(mu_r*(R1 - R2)))];
        J = J1\J2;
    end
    detJ = det(J);
    invJ = inv(J);
    T = invJ*invJ'*detJ;
    T = T'; T = T(:); T(5) = detJ;
    f(i,:) = T;
end
if fnum ~= 0
    f = f(:,fnum);
end
*/