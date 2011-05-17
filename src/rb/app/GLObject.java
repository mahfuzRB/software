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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.ShortBuffer;

import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;

import android.content.Context;
import android.util.Log;

public class GLObject extends Object{
	
    private int reg_num; // number of subdomains
    private int[] node_reg; // tell us which subdomain our vertices belong to
    private int[] face_dom; // tell us which subdomain our faces belong to
    
    // the bounding box (xyz range) of the model
    public float[] nminmax = {1e9f, 1e9f, 1e9f, -1e9f, -1e9f, -1e9f};    
    public int node_num; // number of vertices
    public int face_num; // number of faces
    public float[] node = null; // vertex data
    public float[] normal; // vertex normal data
    public float[] fnormal; // face normal data
    public float[] ref_node; // the original vertex data
    public short[] face; // face data
    public short[] face_wf; // edge data
    public float[][] ucolor; // color data
    public float boxsize; // bounding box size
    public int[] frame_num; // number of animation frame for solution field
    public float[][] sol; // field solution data
    public int field_num = 1; // number of solution field
    
    boolean is2D = true; // is our model 2D?
    boolean isconstant = false; // is the field solution constant?
    
    public boolean isgeoani = false;
    public int vframe_num = 1; // number of animation frame for vertices
    public float[][][] vLTfunc = null;
    public float[] vnode = null; // animation node data    
    
    public ShortBuffer _shortBuffer;
    public FloatBuffer _floatBuffer;
    
    private Context context;
    
    public void allocateBuffer(){
    	int SHORT_MAX = 250000;
    	int FLOAT_MAX = 1000000;
    	
    	Log.d("GLRenderer","Allocate (short):"+SHORT_MAX*2+" bytes");
        ByteBuffer vbb = ByteBuffer.allocateDirect(SHORT_MAX*2);
        vbb.order(ByteOrder.nativeOrder());
        _shortBuffer = vbb.asShortBuffer();
        _shortBuffer.position(0);
        
        Log.d("GLRenderer","Allocate (float):"+FLOAT_MAX*4+" bytes");
        ByteBuffer fbb = ByteBuffer.allocateDirect(FLOAT_MAX*4);
        fbb.order(ByteOrder.nativeOrder());
        _floatBuffer = fbb.asFloatBuffer();
        _floatBuffer.position(0);                
    }
    
    public GLObject(){
    	// void constructor
    }
    
    public GLObject(Context _context){
    	context = _context;
    }

    // get the bounding box data
    public float[] get_minmax(){
    	return nminmax;
    }
    
    // calculate the color data (red green blue alpha) from the solution field
    public void cal_color_data(){
    	float calphaval = 0.8f; // default alpha value, use 1.0f for nonblend rendering
    	
    	ucolor = null;
    	// all field must have the same length (for now)
    	ucolor = new float[field_num][sol[0].length*4];
    	
    	isconstant = true;
    	for (int ifn = 0; ifn < field_num; ifn++)
    	{
    		// range of the current solution field
	    	float val_min = sol[ifn][0];
	    	float val_max = sol[ifn][0];
	    	for (int i = 0; i < sol[ifn].length; i++){
	    		val_min = (val_min>sol[ifn][i])?sol[ifn][i]:val_min;
	    		val_max = (val_max<sol[ifn][i])?sol[ifn][i]:val_max;
	    	}
	    	if (Math.abs(val_max-val_min)>1e-8){
	    		isconstant = false;
	    		
	    		// calculate color data
		    	for (int i = 0; i < sol[0].length; i++){    		
			    	float tmpvar = (sol[ifn][i]-val_min)/(val_max-val_min);
			    	if (tmpvar<=0.125f)
			    	{
			    		ucolor[ifn][i*4+0] = 0.0f;
			    		ucolor[ifn][i*4+1] = 0.0f;
			    		ucolor[ifn][i*4+2] = 0.5f + tmpvar/0.25f;
			    		ucolor[ifn][i*4+3] = calphaval; 
			    	}
			    	if ((tmpvar>0.125f)&&(tmpvar<=0.375f))
			    	{
			    		ucolor[ifn][i*4+0] = 0.0f;
			    		ucolor[ifn][i*4+1] = 0.0f + (tmpvar-0.125f)/0.25f;
			    		ucolor[ifn][i*4+2] = 1.0f;
			    		ucolor[ifn][i*4+3] = calphaval;
			    	}
			    	if ((tmpvar>0.375f)&&(tmpvar<=0.625f))
			    	{
			    		ucolor[ifn][i*4+0] = 0.0f + (tmpvar-0.375f)/0.25f;
			    		ucolor[ifn][i*4+1] = 1.0f;
			    		ucolor[ifn][i*4+2] = 1.0f - (tmpvar-0.375f)/0.25f;
			    		ucolor[ifn][i*4+3] = calphaval;
			    	}
			    	if ((tmpvar>0.625f)&&(tmpvar<=0.875f))
			    	{
			    		ucolor[ifn][i*4+0] = 1.0f;
			    		ucolor[ifn][i*4+1] = 1.0f - (tmpvar-0.625f)/0.25f;
			    		ucolor[ifn][i*4+2] = 0.0f;
			    		ucolor[ifn][i*4+3] = calphaval;
			    	}
			    	if (tmpvar>0.875f)
			    	{
			    		ucolor[ifn][i*4+0] = 1.0f - (tmpvar-0.875f)/0.25f;
			    		ucolor[ifn][i*4+1] = 0.0f;
			    		ucolor[ifn][i*4+2] = 0.0f;
			    		ucolor[ifn][i*4+3] = calphaval;
			    	}
		    	}
	    	} else {
	    		for (int i = 0; i < sol[0].length; i++){
	    			ucolor[ifn][i*4+0] = 0.0f;
		    		ucolor[ifn][i*4+1] = 0.0f;
		    		ucolor[ifn][i*4+2] = 1.0f;
		    		ucolor[ifn][i*4+3] = calphaval;
	    		}
	    	}
	    	
    	}
    }
    
    // calculate normal data for the current model
    public void cal_normal_data(){
    	normal = new float[node_num*3];
    	fnormal = new float[face_num*3];
    	float[] vecAB = new float[3];
    	float[] vecAC = new float[3];
    	int i, j, k;
    	float length;
    	// Initialize normal data and contribution flag
    	int[] icount = new int[node_num];
    	for (i = 0; i < node_num; i++){
    		normal[i*3+0] = 0.0f; normal[i*3+1] = 0.0f; normal[i*3+2] = 0.0f;
    		icount[i] = 0;
    	}
    	// calculate local face normal
    	for (i = 0; i < face_num; i++){
    		for (j = 0; j < 3; j++){
    			vecAB[j] = node[face[i*3+1]*3+j] - node[face[i*3+0]*3+j];
    			vecAC[j] = node[face[i*3+2]*3+j] - node[face[i*3+0]*3+j];    			
    		}
    		// normal of the face is the cross product of AB and AC
    		fnormal[i*3+0] = vecAB[1]*vecAC[2] - vecAB[2]*vecAC[1];
    		fnormal[i*3+1] = vecAB[2]*vecAC[0] - vecAB[0]*vecAC[2];
    		fnormal[i*3+2] = vecAB[0]*vecAC[1] - vecAB[1]*vecAC[0];
    		// normalize
    		length = (float)Math.sqrt((fnormal[i*3+0]*fnormal[i*3+0] + fnormal[i*3+1]*fnormal[i*3+1] + fnormal[i*3+2]*fnormal[i*3+2]));
    		for (j = 0; j < 3; j++) 
    			fnormal[i*3+j] = fnormal[i*3+j]/length;
    		// add in contribution to all three vertices 
    		for (j = 0; j < 3; j++){
    			icount[face[i*3+j]]++;
    			for (k = 0; k < 3; k++)
    				normal[face[i*3+j]*3+k] += fnormal[i*3+k];
    		}    		
    	}
    	// average and normalize all normal vectors
    	for (i = 0; i < node_num; i++){
    		for (j = 0; j < 3; j++) 
    			normal[i*3+j] = normal[i*3+j]/icount[i];
    		length = (float)Math.sqrt((normal[i*3+0]*normal[i*3+0] + normal[i*3+1]*normal[i*3+1] + normal[i*3+2]*normal[i*3+2]));
    		for (j = 0; j < 3; j++) 
    			normal[i*3+j] = normal[i*3+j]/length;
    	}    	
    }

    // move the model center to (0,0,0)
    public void model_centerize(){
    	cal_boxsize();
    	float xcen = 0.5f*(nminmax[0]+nminmax[3]);
    	float ycen = 0.5f*(nminmax[1]+nminmax[4]);
    	float zcen = 0.5f*(nminmax[2]+nminmax[5]);
    	for (int i = 0; i < node_num; i++){
    		node[i*3+0] -= xcen;
    		node[i*3+1] -= ycen;
    		node[i*3+2] -= zcen;
    	}
    	// recalculating minmax box
    	nminmax[0] -= xcen; nminmax[3] += xcen;
    	nminmax[1] -= ycen; nminmax[4] += ycen;
    	nminmax[2] -= zcen; nminmax[5] += zcen;
    }
    
    // calculate bounding box data
    public void cal_boxsize(){
    	nminmax[0] =  1e9f; nminmax[1] =  1e9f; nminmax[2] = 1e9f;
    	nminmax[3] = -1e9f; nminmax[4] = -1e9f; nminmax[5] = -1e9f;
    	for (int i = 0; i < node_num; i++){
			for (int j = 0; j < 3; j++){
				nminmax[0+j] = (nminmax[0+j]>node[i*3+j])?node[i*3+j]:nminmax[0+j];
				nminmax[3+j] = (nminmax[3+j]<node[i*3+j])?node[i*3+j]:nminmax[3+j];
			}				
		}
    	is2D = false;
    	if (Math.abs(nminmax[5]-nminmax[2])<1e-8)
    		is2D = true;
    	    	
		boxsize = 0.0f;
		boxsize = (nminmax[3]-nminmax[0])>boxsize?(nminmax[3]-nminmax[0]):boxsize;
		boxsize = (nminmax[4]-nminmax[1])>boxsize?(nminmax[4]-nminmax[1]):boxsize;
		boxsize = (nminmax[5]-nminmax[2])>boxsize?(nminmax[5]-nminmax[2]):boxsize;    
    }
    
    // read geometry data
    public void read_offline_data(String directory_name, boolean isAssetFile)
		throws IOException {
    	
    	HttpClient client = new DefaultHttpClient();
		int buffer_size = 8192;
		
		// Read in output data
		InputStreamReader isr;
		String dataString = directory_name + "/geometry.dat";
		
		if(!isAssetFile) {
			HttpGet request = new HttpGet(dataString);
			HttpResponse response = client.execute(request);
			isr = new InputStreamReader(response.getEntity()
					.getContent());
		}
		else { // Read from assets
			isr = new InputStreamReader(
					context.getAssets().open(dataString));
		}
		BufferedReader reader = new BufferedReader(isr,buffer_size);

		String line = reader.readLine();
		String[] tokens = line.split(" ");

		node_num = Integer.parseInt(tokens[0]);
		int count = 1;
		ref_node = new float[node_num*3];
		node = new float[node_num*3];
		sol = new float[field_num][node_num];
		ucolor = new float[field_num][node_num*4];
		for (int i = 0; i < node_num; i++){
			ref_node[i*3+0] = Float.parseFloat(tokens[count]);
			ref_node[i*3+1] = Float.parseFloat(tokens[count+1]);
			ref_node[i*3+2] = Float.parseFloat(tokens[count+2]);
			count += 3;
			node[i*3+0] = ref_node[i*3+0];
			node[i*3+1] = ref_node[i*3+1];
			node[i*3+2] = ref_node[i*3+2];
		}		
		cal_boxsize();
		
		reg_num = Integer.parseInt(tokens[count]);
		face_num = Integer.parseInt(tokens[count+1]);
		count += 2;
		face = new short[face_num*3];
		for (int i = 0; i < face_num; i++){
			face[i*3+0] = Short.parseShort(tokens[count]);
			face[i*3+1] = Short.parseShort(tokens[count+1]);
			face[i*3+2] = Short.parseShort(tokens[count+2]);
			count += 3;
		}
    	node_reg = new int[node_num];
    	for (int i = 0; i < node_num; i++){
    		node_reg[i] = Integer.parseInt(tokens[count]);
    		count++;
    	}
    	face_dom = new int[face_num];
    	for (int i = 0; i < face_num; i++){
    		face_dom[i] = Integer.parseInt(tokens[count]);
    		count++;
    	}
    	isr.close();
    	
    	// Create a wireframe list (edge data)
        face_wf = new short[face_num*3*2]; 
        for (int i = 0; i < face_num; i++){
        	face_wf[i*6+0*2+0] = face[i*3+0]; face_wf[i*6+0*2+1] = face[i*3+1];
        	face_wf[i*6+1*2+0] = face[i*3+1]; face_wf[i*6+1*2+1] = face[i*3+2];
        	face_wf[i*6+2*2+0] = face[i*3+2]; face_wf[i*6+2*2+1] = face[i*3+0];     	
        }
        
        if (!is2D())
        	cal_normal_data();
    }
    
    // Affine transformation
    // LTfunc is of size [number of subdomain, 12]
    // a row of LTfunc is the rowwise flatten of the [3,3] transformation matrix and the [3,1] translation vector
    public void nodal_transform(float[][] LTfunc){
    	float[] old_node = new float[3];
        for (int i = 0; i < node_num; i++){
        	old_node[0] = ref_node[i*3+0];
        	old_node[1] = ref_node[i*3+1];
        	old_node[2] = ref_node[i*3+2];	
        	node[i*3+0] = LTfunc[node_reg[i]][0]*old_node[0] + LTfunc[node_reg[i]][1]*old_node[1] + LTfunc[node_reg[i]][2]*old_node[2]
        				  + LTfunc[node_reg[i]][9];
        	node[i*3+1] = LTfunc[node_reg[i]][3]*old_node[0] + LTfunc[node_reg[i]][4]*old_node[1] + LTfunc[node_reg[i]][5]*old_node[2]
        	              + LTfunc[node_reg[i]][10];
        	node[i*3+2] = LTfunc[node_reg[i]][6]*old_node[0] + LTfunc[node_reg[i]][7]*old_node[1] + LTfunc[node_reg[i]][8]*old_node[2]
        	              + LTfunc[node_reg[i]][11];
        }
        cal_boxsize();
        model_centerize();
    }
    
    // assign 1 field solution
    public void set_field_data(float[] _val){
    	sol = null;
    	field_num = 1;
    	sol = new float[1][_val.length];
    	sol[0] = _val;
    	frame_num = new int[1];
    	frame_num[0] = sol[0].length/node_num;
    	cal_color_data();
    }

    public void set_field_data(float[] _val1, float[] _val2){
    	set_field_data(_val1, _val2, true);
    }
    
    // assign 2 field solutions
    public void set_field_data(float[] _val1, float[] _val2, boolean isDeformed){
    	if (isDeformed){
	    	sol = null;    	
	    	// the solution field is the displacement field
	    	// merge displacement field into current vertex data
	    	field_num = 1;
	    	float val_min = _val1[0];
	    	float val_max = _val1[0];
	    	for (int i = 0; i < _val1.length; i++){
	    		val_min = (val_min>_val1[i])?_val1[i]:val_min;
	    		val_max = (val_max<_val1[i])?_val1[i]:val_max;
	    	}
	    	for (int i = 0; i < _val2.length; i++){
	    		val_min = (val_min>_val2[i])?_val2[i]:val_min;
	    		val_max = (val_max<_val2[i])?_val2[i]:val_max;
	    	}
	    	float sval = (val_max-val_min)/boxsize*5;
	    	if ((_val1.length/node_num == 1) && (_val2.length/node_num == 1)){
		    	sol = new float[1][node_num];
		    	for (int i = 0; i<node_num; i++){
		    		node[i*3+0] = node[i*3+0] + _val1[i]/sval;
		    		node[i*3+1] = node[i*3+1] + _val2[i]/sval;
		    		node[i*3+2] = node[i*3+2];
		    		sol[0][i] = 0.0f;
		    	}
		    	frame_num = new int[1];
	        	frame_num[0] = 1;
	    	} else {
	    		int _vframe_num = (_val1.length/node_num);
		    	sol = new float[1][node_num*_vframe_num];
		    	for (int i = 0; i < node_num; i++)
		    		for (int j = 0; j < _vframe_num; j++){
			    		vnode[j*node_num*3+i*3+0] = vnode[j*node_num*3+i*3+0] + _val1[j*node_num+i]/sval;
			    		vnode[j*node_num*3+i*3+1] = vnode[j*node_num*3+i*3+1] + _val2[j*node_num+i]/sval;
			    		vnode[j*node_num*3+i*3+2] = vnode[j*node_num*3+i*3+2];
			    		sol[0][j*node_num+i] = 0.0f;
			    	}
		    	frame_num = new int[1];
		    	frame_num[0] = _vframe_num;
	    	}    	
	    	
    	}
    	else{
    		sol = null;    	
	    	field_num = 2;
	    	if ((_val1.length/node_num == 1) && (_val2.length/node_num == 1)){
		    	sol = new float[2][node_num];
			    for (int i = 0; i<node_num; i++){
			    	sol[0][i] = _val1[i];
			    	sol[1][i] = _val2[i];
			    }
			    frame_num = new int[2];
		       	frame_num[0] = 1;
		       	frame_num[1] = 1;
	    	} else {
	    		int _vframe_num = (_val1.length/node_num);
		    	sol = new float[2][node_num*_vframe_num];
		    	for (int i = 0; i < node_num; i++)
		    		for (int j = 0; j < _vframe_num; j++){			    		
			    		sol[0][j*node_num+i] = _val1[j*node_num+i];
			    		sol[1][j*node_num+i] = _val2[j*node_num+i];
			    	}
		    	frame_num = new int[2];
		    	frame_num[0] = _vframe_num;
		    	frame_num[1] = _vframe_num;
	    	}
	    		
    	}
    	
    	cal_color_data();
    }
    
    public void set_field_data(float[] _val1, float[] _val2, float[] _val3){
    	set_field_data(_val1, _val2, _val3, true);
    }
    // assign 3 field solutions
    public void set_field_data(float[] _val1, float[] _val2, float[] _val3, boolean isDeformed){
    	if (isDeformed){
	    	sol = null;    	
	    	// the solution field is the displacement field
	    	// merge displacement field into current vertex data
	    	field_num = 1;
	    	float val_min = _val1[0];
	    	float val_max = _val1[0];
	    	for (int i = 0; i < _val1.length; i++){
	    		val_min = (val_min>_val1[i])?_val1[i]:val_min;
	    		val_max = (val_max<_val1[i])?_val1[i]:val_max;
	    	}
	    	for (int i = 0; i < _val2.length; i++){
	    		val_min = (val_min>_val2[i])?_val2[i]:val_min;
	    		val_max = (val_max<_val2[i])?_val2[i]:val_max;
	    	}
	    	for (int i = 0; i < _val3.length; i++){
	    		val_min = (val_min>_val3[i])?_val3[i]:val_min;
	    		val_max = (val_max<_val3[i])?_val3[i]:val_max;
	    	}
	    	float sval = (val_max-val_min)/boxsize*5;
	    	if ((_val1.length/node_num == 1) && (_val2.length/node_num == 1) && (_val3.length/node_num == 1)){
		    	sol = new float[1][node_num];
		    	for (int i = 0; i<node_num; i++){
		    		node[i*3+0] = node[i*3+0] + _val1[i]/sval;
		    		node[i*3+1] = node[i*3+1] + _val2[i]/sval;
		    		node[i*3+2] = node[i*3+2] + _val3[i]/sval;;
		    		sol[0][i] = 0.0f;
		    		frame_num = new int[1];
		        	frame_num[0] = 1;
		    	}
	    	} else {
	    		int _vframe_num = (_val1.length/node_num);
		    	sol = new float[1][node_num*_vframe_num];
		    	for (int i = 0; i < node_num; i++)
		    		for (int j = 0; j < _vframe_num; j++){
			    		vnode[j*node_num*3+i*3+0] = vnode[j*node_num*3+i*3+0] + _val1[j*node_num+i]/sval;
			    		vnode[j*node_num*3+i*3+1] = vnode[j*node_num*3+i*3+1] + _val2[j*node_num+i]/sval;
			    		vnode[j*node_num*3+i*3+2] = vnode[j*node_num*3+i*3+2] + _val3[j*node_num+i]/sval;
			    		sol[0][j*node_num+i] = 0.0f;
			    	}
		    	frame_num = new int[1];
		    	frame_num[0] = _vframe_num;
	    	}
    	} else {
    		sol = null;    	
	    	field_num = 3;
	    	if ((_val1.length/node_num == 1) && (_val2.length/node_num == 1) && (_val3.length/node_num == 1)){
		    	sol = new float[3][node_num];
			    for (int i = 0; i<node_num; i++){
			    	sol[0][i] = _val1[i];
			    	sol[1][i] = _val2[i];
			    	sol[2][i] = _val3[i];
			    }
			    frame_num = new int[3];
		       	frame_num[0] = 1;
		       	frame_num[1] = 1;
		       	frame_num[2] = 1;
	    	} else {
	    		int _vframe_num = (_val1.length/node_num);
		    	sol = new float[3][node_num*_vframe_num];
		    	for (int i = 0; i < node_num; i++)
		    		for (int j = 0; j < _vframe_num; j++){			    		
			    		sol[0][j*node_num+i] = _val1[j*node_num+i];
			    		sol[1][j*node_num+i] = _val2[j*node_num+i];
			    		sol[2][j*node_num+i] = _val3[j*node_num+i];
			    	}
		    	frame_num = new int[3];
		    	frame_num[0] = _vframe_num;
		    	frame_num[1] = _vframe_num;
		    	frame_num[2] = _vframe_num;
	    	}
    	}
    	cal_color_data();
    }
    
    // assign 4 field solutions (4th field is sol col)
    public void set_field_data(float[] _val1, float[] _val2, float[] _val3, float[] _val4){
    	sol = null;
    	
    	// the solution field is the displacement field
    	// merge displacement field into current vertex data
    	field_num = 1;
    	float val_min = _val1[0];
    	float val_max = _val1[0];
    	for (int i = 0; i < _val1.length; i++){
    		val_min = (val_min>_val1[i])?_val1[i]:val_min;
    		val_max = (val_max<_val1[i])?_val1[i]:val_max;
    	}
    	for (int i = 0; i < _val2.length; i++){
    		val_min = (val_min>_val2[i])?_val2[i]:val_min;
    		val_max = (val_max<_val2[i])?_val2[i]:val_max;
    	}
    	for (int i = 0; i < _val3.length; i++){
    		val_min = (val_min>_val3[i])?_val3[i]:val_min;
    		val_max = (val_max<_val3[i])?_val3[i]:val_max;
    	}
    	float sval = (val_max-val_min)/boxsize*5.0f;
    	if ((_val1.length/node_num == 1) && (_val2.length/node_num == 1) && (_val3.length/node_num == 1)){
	    	sol = new float[1][node_num];
	    	for (int i = 0; i<node_num; i++){
	    		node[i*3+0] = node[i*3+0] + _val1[i]/sval;
	    		node[i*3+1] = node[i*3+1] + _val2[i]/sval;
	    		node[i*3+2] = node[i*3+2] + _val3[i]/sval;
	    		sol[0][i] = _val4[i];
	    		frame_num = new int[1];
	        	frame_num[0] = 1;
	    	}
    	} else {
    		int _vframe_num = (_val1.length/node_num);
	    	sol = new float[1][node_num*_vframe_num];
	    	for (int i = 0; i < node_num; i++)
	    		for (int j = 0; j < _vframe_num; j++){
		    		vnode[j*node_num*3+i*3+0] = vnode[j*node_num*3+i*3+0] + _val1[j*node_num+i]/sval;
		    		vnode[j*node_num*3+i*3+1] = vnode[j*node_num*3+i*3+1] + _val2[j*node_num+i]/sval;
		    		vnode[j*node_num*3+i*3+2] = vnode[j*node_num*3+i*3+2] + _val3[j*node_num+i]/sval;
		    		sol[0][j*node_num+i] = _val4[j*node_num+i];
		    	}
	    	frame_num = new int[1];
	    	frame_num[0] = _vframe_num;
    	}    	
    	cal_color_data();
    }
    
    // assign vertex data
    public void set_node_data(float[] _node){
    	node_num = _node.length/3;
    	node = null;
    	node = _node;
    	cal_boxsize();
    	model_centerize();
    }
    
    // assign node region data
    public void set_node_reg_data(int[] _node_reg){
    	node_reg = null;
    	node_reg = _node_reg;    	
    }
    
    // assign ref_node data
    public void set_ref_node_data(float[] _node){
    	node_num = _node.length/3;
    	ref_node = null;
    	ref_node = _node;
    	// copy ref_node to node
    	if (node == null){
    		node = new float[node_num*3];
    		for (int i = 0; i < node_num; i++)
    			for (int j = 0; j < 3; j++)
    				node[i*3+j] = ref_node[i*3+j];
    	}
    	cal_boxsize();
    	model_centerize();
    	//node = null;
    }
    
    // assign face data
    public void set_face_data(short[] _face){
    	face_num = _face.length/3;
    	face = null;
    	face = _face;
    	
    	// Create a wireframe list
    	face_wf = null;
        face_wf = new short[face_num*3*2]; 
        for (int i = 0; i < face_num; i++){
        	face_wf[i*6+0*2+0] = face[i*3+0]; face_wf[i*6+0*2+1] = face[i*3+1];
        	face_wf[i*6+1*2+0] = face[i*3+1]; face_wf[i*6+1*2+1] = face[i*3+2];
        	face_wf[i*6+2*2+0] = face[i*3+2]; face_wf[i*6+2*2+1] = face[i*3+0];     	
        } 
        
        // Calculate normal data for 3D object
        if (!is2D())
        	cal_normal_data();
    }    
    
    // how many field we currently have?
    public float[] get_field_data(int ifn){
    	return sol[ifn];
    }
    
    // the vertex data
    public float[] get_node_data(){
    	return node;
    }
    
    // the reference vertex data
    public float[] get_ref_node_data(){
    	return ref_node;
    }
    
    // the face data
    public short[] get_face_data(){
    	return face;
    }
    
    // the number of subdomain
    public int get_reg_num(){
    	return reg_num;
    }
 
    // is our model flat (2D)?
    public boolean is2D(){
    	return is2D;
    }
    
    // get node_reg data
    public int[] get_node_reg(){
    	return node_reg;
    }
    
    // get LTfunc data
    // this also get us vertex animation data
    public void set_LTfunc(float[][][] _LTfunc, int _reg_num, int _vframe_num){    	
    	reg_num = _reg_num;
    	vframe_num = _vframe_num;
    	if (vframe_num == 0)
    		isgeoani = false;
    	else
    		isgeoani = true;
    	vLTfunc = _LTfunc;
    	vnode = new float[vframe_num*node_num*3];
    	for (int i = 0; i < vframe_num; i++){
    		// get current nodal data
    		nodal_transform(vLTfunc[i]);
    		// copy current nodal data into animation list
    		for (int j = 0; j < node_num; j++)
    			for (int k = 0; k < 3; k++)
    				vnode[i*node_num*3 + j*3 + k] = node[j*3+k];
    	}    	
    }
    
    // get LTfunc data
    // this also get us vertex animation data
    public void set_LTfunc(float[][][] _LTfunc){    	
    	vframe_num = _LTfunc.length;
    	if (vframe_num == 1)
    		isgeoani = false;
    	else
    		isgeoani = true;
    	vLTfunc = _LTfunc;
    	vnode = new float[vframe_num*node_num*3];
    	for (int i = 0; i < vframe_num; i++){
    		// get current nodal data
    		nodal_transform(vLTfunc[i]);
    		// copy current nodal data into animation list
    		for (int j = 0; j < node_num; j++)
    			for (int k = 0; k < 3; k++)
    				vnode[i*node_num*3 + j*3 + k] = node[j*3+k];
    	}    	
    }
    
    public void mesh_transform_custom(Parameter[] mu){    	
    	vframe_num = mu.length;
    	if (vframe_num == 1)
    		isgeoani = false;
    	else
    		isgeoani = true;    	
    	vnode = new float[vframe_num*node_num*3];
    	for (int i = 0; i < vframe_num; i++){
    		// get current nodal data
    		float[] tmpnode = RBActivity.mRbSystem.mesh_transform(mu[i].getArray(), ref_node.clone());
    		Log.d("GLRenderer",mu[i].getEntry(0)+" "+mu[i].getEntry(1));
    		Log.d("GLRenderer",tmpnode[4]+" "+node[4]);
    		node = tmpnode.clone();
    		// copy current nodal data into animation list
    		for (int j = 0; j < node_num; j++)
    			for (int k = 0; k < 3; k++){
    				vnode[i*node_num*3 + j*3 + k] = tmpnode[j*3+k];    				
    			}
    	}
    }
}
