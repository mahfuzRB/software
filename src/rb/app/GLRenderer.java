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

import java.nio.FloatBuffer;
import java.nio.ShortBuffer;

import javax.microedition.khronos.egl.EGLConfig;
import javax.microedition.khronos.opengles.GL10;

import android.opengl.GLSurfaceView;
import android.util.Log;
 
public class GLRenderer implements GLSurfaceView.Renderer {
    private static final String LOG_TAG = GLRenderer.class.getSimpleName();
 
    GLObject _object;
    
    // Camera control
    private Camera camera;
    
    private ShortBuffer _shortBuffer;
    private FloatBuffer _floatBuffer;
    
    private int _vertex_off = 0;
    private int _index_off = 0;
    private int _indexwf_off = 0;
    private int _normal_off = 0;
    private int[] _color_off;
    private int _anivertex_off = 0;
    
    private float _width = 480f;
    private float _height = 800f;
    private float[] AR = {1f,1f}; // aspect ratio
    
    private float pos[] = {0f, 0f, 0f}; // touchscreeen control data
    public float scale_rat = 1.0f; // scaling ratio (for zooming)
    
    private int current_field = 0;
    private int current_frame = 0;
    private float current_framef = 0f;
    
    public boolean ispaused = false;
    public boolean isconstant = false;
    
    public boolean isFrontFace = true;
    
    public boolean isContinuousRotation = true;
 
    public void setOrientation(boolean pmode){
    	if (pmode) { // portrait mode
    		_width = 480f; _height = 762f;
    	    AR[0] = 1.0f; AR[1] = _width/_height;
    	}
    	else { // landscape mode
    		_width = 762f; _height = 480f;
    	    AR[0] = _height/_width; AR[1] = 1f;
    	}    	
    }
     /*
    public GLRenderer(){
    	_object = new GLObject();
    }
    */
    public GLRenderer(GLObject _input_object){
    	_object = _input_object;
    }
    
    //@Override
    public void onSurfaceCreated(GL10 gl, EGLConfig config) {    	
    	
    	initRendering();
    	camera = new Camera();
    	// set initial position away from the model in the y-direction
    	// looking toward the center of the model (0,0,0) horizontally
    	camera.setCamera(0f, -_object.boxsize, 0f, 0f, 1f, 0f, 0f, 0f, 1f);
    	
    	// define the color we want to be displayed as the "clipping wall"
        gl.glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        gl.glMatrixMode(GL10.GL_PROJECTION);
        
        float exrat; // marginal extension ratio 
        if (_object.is2D())
        	exrat = 0.65f;
        else
        	exrat = 0.95f;
        // orthographic view
        gl.glOrthof(-exrat*_object.boxsize/AR[0], exrat*_object.boxsize/AR[0],
        		    -exrat*_object.boxsize/AR[1], exrat*_object.boxsize/AR[1],
        		    -100,100);
        	
        gl.glViewport(0, 0, (int) _width, (int) _height);
        gl.glMatrixMode(GL10.GL_MODELVIEW);
        
        // define the color we want to be displayed as the "clipping wall"
        //gl.glClearColor(0f, 0f, 0f, 1.0f);
        gl.glClearColor(1f, 1f, 1f, 1.0f);
    
	    // enable the differentiation of which side may be visible 
	    gl.glEnable(GL10.GL_CULL_FACE);
	    // which is the front? the one which is drawn counter clockwise
	    gl.glFrontFace(GL10.GL_CCW);
	    // which one should NOT be drawn
	    gl.glCullFace(GL10.GL_BACK);
    
    	gl.glEnableClientState(GL10.GL_VERTEX_ARRAY);    	
    	gl.glEnableClientState(GL10.GL_COLOR_ARRAY);
    	
    	// Enable normal for 3D object
    	if (!_object.is2D())
    		gl.glEnableClientState(GL10.GL_NORMAL_ARRAY);
    }
 
    //@Override
    public void onSurfaceChanged(GL10 gl, int w, int h) {
        gl.glViewport(0, 0, w, h);
    }
 
    private void initRendering() {
    	int short_count = 0;
    	int float_count = 0;
    	
    	_shortBuffer = RBActivity.mRbModel._shortBuffer;
    	_shortBuffer.clear();
    	//_shortBuffer.position(0);    	
    	_floatBuffer = RBActivity.mRbModel._floatBuffer;
    	_floatBuffer.clear();
    	//_floatBuffer.position(0);
    	
    	// Init vertex, face, color and normal openGL buffers
    	
        _vertex_off = float_count;
        _floatBuffer.put(_object.node);
        float_count += _object.node.length;
        Log.d("GLRenderer","float_count (vertex) = "+float_count+"/"+_floatBuffer.capacity());
        
        _index_off = short_count;
        _shortBuffer.put(_object.face);
        short_count += _object.face.length;
        Log.d("GLRenderer","short_count (index) = "+short_count+"/"+_shortBuffer.capacity());
        
        _indexwf_off = short_count;
        _shortBuffer.put(_object.face_wf);
        short_count += _object.face_wf.length;
        Log.d("GLRenderer","short_count (indexwf) = "+short_count+"/"+_shortBuffer.capacity());
        
        // animation color buffer has (frame_num) the time of data
        _color_off = new int[_object.field_num];
        for (int i = 0; i < _object.field_num; i++){        	
        	_color_off[i] = float_count;
            _floatBuffer.put(_object.ucolor[i]);
            float_count += _object.ucolor[i].length;
            Log.d("GLRenderer","float_count (color["+i+"]) = "+float_count+"/"+_floatBuffer.capacity());
        }
        
        // init array for 3D object
        if (!_object.is2D()){        	
        	_normal_off = float_count;
            _floatBuffer.put(_object.normal);
            float_count += _object.normal.length;
            Log.d("GLRenderer","float_count (normal) = "+float_count+"/"+_floatBuffer.capacity());
        }
        
        // init vertex animation buffer
        if (_object.isgeoani){        	
        	_anivertex_off = float_count;
	        _floatBuffer.put(_object.vnode);
	        float_count += _object.vnode.length;
	        Log.d("GLRenderer","float_count (anivertex) = "+float_count+"/"+_floatBuffer.capacity());
        } else {        	
        	_anivertex_off = float_count;
        	_floatBuffer.put(_object.node);
        	float_count += _object.node.length;
        	Log.d("GLRenderer","float_count (anivertex) = "+float_count+"/"+_floatBuffer.capacity());
        }
	     
    }
    
    // pass touchscreeen control data
    public void setPos(boolean iCR, float posx, float posy, float posz){
    	pos[0] += posx;    	pos[1] += posy;    	pos[2] += posz;
    	isContinuousRotation = iCR; 
    }
    
    // set the current rendering field
    public void setcField(int cField){
    	if (cField > (_object.field_num - 1))
    		cField = 0;
    	current_field = cField;
    }
    
    // get the current rendering field
    public int getcField(){
    	return current_field;
    }
    
    //@Override
    public void onDrawFrame(GL10 gl) {
    	// clear the screen to black (0,0,0) color
    	//gl.glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    	gl.glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glClear(GL10.GL_COLOR_BUFFER_BIT | GL10.GL_DEPTH_BUFFER_BIT);
        if (!_object.is2D()) // enable depth test for 3D rendering
        	gl.glEnable(GL10.GL_DEPTH_TEST);
        
        
        //gl.glEnable(GL10.GL_CULL_FACE);
    	if ((isFrontFace) || (_object.is2D())){
    		// enable blending (for rendering wireframe)
    		if (!_object.is2D())
    			gl.glDisable(GL10.GL_CULL_FACE);
    		gl.glFrontFace(GL10.GL_CCW);
    	}
	    else{
	    	gl.glEnable(GL10.GL_CULL_FACE);	
    		gl.glFrontFace(GL10.GL_CW);
	    }
	    
        // reset transformation matrix
        gl.glLoadIdentity();
        
        // setup Light
        if (!_object.is2D()){
        	gl.glEnable(GL10.GL_LIGHTING); // Enable light
        	gl.glEnable(GL10.GL_LIGHT0); // turn on the light
        	gl.glEnable(GL10.GL_COLOR_MATERIAL); // turn on color lighting
        	        	
        	// material shininess
        	gl.glMaterialx(GL10.GL_FRONT_AND_BACK, GL10.GL_SHININESS, 128);
        	// ambient light
        	float lightAmbient[] = {0.5f, 0.5f, 0.5f, 1.0f};
        	gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_AMBIENT, lightAmbient, 0);
        	// diffuse light
        	float lightDiffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
        	gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_DIFFUSE, lightDiffuse, 0);
        	// specular light
        	float[] lightSpecular = {0.7f, 0.7f, 0.7f, 1.0f};
        	gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_SPECULAR, lightSpecular, 0);
        	// light position
        	float[] lightPosition = {-_object.boxsize, -_object.boxsize, 0.0f, 0.0f};
        	gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_POSITION, lightPosition, 0);
        	// light direction
        	float[] lightDirection = {_object.boxsize, _object.boxsize, _object.boxsize, 0.0f};
        	gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_SPOT_DIRECTION, lightDirection, 0);
        	// 90 degree FOV
        	gl.glLightf(GL10.GL_LIGHT0, GL10.GL_SPOT_CUTOFF, 45.0f);        	
            
        	// using our normal data
        	_floatBuffer.position(_normal_off);
          	gl.glNormalPointer(GL10.GL_FLOAT, 0, _floatBuffer);
        }
        
        // zoom in/out the model
        gl.glScalef(scale_rat, scale_rat, scale_rat);
        
        // touchscreen control
        if (_object.is2D()) // we just move the object around in 2D cases 
        	gl.glTranslatef(pos[0]*_object.boxsize/20f, pos[1]*_object.boxsize/20f, pos[2]*_object.boxsize/20f);
        else{ // but we rotate the object in 3D cases
        	// set yawing/pitching rotation angles and update camera
        	camera.SetRotation(-pos[0]*8f, -pos[1]*8f);
        	// update rotation matrix
        	gl.glMultMatrixf(camera.M,0);
        	// update rotation parameters
        	if (isContinuousRotation){
        		float minrot = 0.02f/scale_rat;
	        	// delay the rotation parameters...
	        	pos[0] = pos[0]*(1-(float)Math.exp(-Math.abs(pos[0])));
	        	pos[1] = pos[1]*(1-(float)Math.exp(-Math.abs(pos[1])));        	
	        	pos[0] = Math.abs(pos[0])>3.00f?Math.signum(pos[0])*3.00f:pos[0]; 
	        	pos[1] = Math.abs(pos[1])>3.00f?Math.signum(pos[1])*3.00f:pos[1]; 
	        	pos[0] = Math.abs(pos[0])>minrot?pos[0]:Math.signum(pos[0])*minrot;
	        	pos[1] = Math.abs(pos[1])>minrot?pos[1]:Math.signum(pos[1])*minrot;
        	} else {
        		// reset the rotation parameters
        		pos[0] = 0.0f;
        		pos[1] = 0.0f;
        	}
	        	
        	//gl.glTranslatef(-camera.Position[0],-camera.Position[1],-camera.Position[2]);  	
        }        
        
        // using our vertex data
        if (!_object.isgeoani){
        	_floatBuffer.position(_vertex_off);
        	gl.glVertexPointer(3, GL10.GL_FLOAT, 0, _floatBuffer);}
        else{
        	_floatBuffer.position(_anivertex_off+current_frame*(_object.node_num*3));
        	gl.glVertexPointer(3, GL10.GL_FLOAT, 0, _floatBuffer);
        }
       
        // specify the color data for the current frame
        _floatBuffer.position(_color_off[current_field]+current_frame*(_object.node_num*4));
    	gl.glColorPointer(4, GL10.GL_FLOAT, 0, _floatBuffer);
        _shortBuffer.position(_index_off);
        gl.glDrawElements(GL10.GL_TRIANGLES, _object.face_num*3, GL10.GL_UNSIGNED_SHORT, _shortBuffer);
        
        // Draw the wireframe for a n field object
        if ((_object.isconstant)|(!_object.is2D())){
        	// Draw the wireframe mesh
        	gl.glColor4f(0.1f, 0.1f, 0.1f, 0.5f);
        	_shortBuffer.position(_indexwf_off);
            gl.glDrawElements(GL10.GL_LINES, _object.face_num*6, GL10.GL_UNSIGNED_SHORT, _shortBuffer);            
        }

        // Draw next animation frame
        if (!ispaused) 
        	increase_frame(0.01f);
    }

    public void zoom(float pzoom){
    	pzoom = (pzoom<1)?1:pzoom;
    	scale_rat = pzoom;    	
    }
    
    // zoom in
    public void zoomin(){
    	scale_rat += 0.1f;    	
    }

    // zoom out
    public void zoomout(){
    	scale_rat -= 0.1f;
    	if (scale_rat < 0.25f)
    		scale_rat = 0.25f;
    }

    // reset zoom parameter
    public void zoomreset(){
    	scale_rat = 1.0f;
    }
    
    // delayed frame increasing, only update animation after 5 frames
    public void increase_frame(float fdelay){
    	current_framef += fdelay*_object.frame_num[current_field];
    	current_frame = (int) Math.round(current_framef);
        if (current_frame >= _object.frame_num[current_field]){
        	current_frame = 0;
        	current_framef = 0;
        }
        if (current_frame < 0){
        	current_frame = _object.frame_num[current_field]-1;
        	current_framef = _object.frame_num[current_field]-1;
        }
    }
    
 // nondelayed frame increasing
    public void increase_ndframe(float fdelay){
    	current_framef += fdelay;
    	current_frame = (int) Math.round(current_framef);
        if (current_frame >= _object.frame_num[current_field]){
        	current_frame = 0;
        	current_framef = 0;
        }
        if (current_frame < 0){
        	current_frame = _object.frame_num[current_field]-1;
        	current_framef = _object.frame_num[current_field]-1;
        }
    }
    
    // pause the animation if there is any
    public void pause(){
    	ispaused = true;
    }
    
    // resume animation
    public void unpause(){
    	ispaused = false;
    }
    
    /*--------------------------------------------------*/
    /* Camera class using quaternion */
    private class Camera{
    	// View position, focus points and view, up and right vectors
    	// view vector pointing from view position to the focus point
    	// up vector perpendicular to the view vector and current horizontal plane
    	// right vector lies on the current horizontal plane
    	public float[] Position, Center, View, Up, Right;
    	public float[] M = new float[16]; // resulted rotation matrix
    	// initialize data
    	public Camera(){
    		Position = new float[3];
    		Center = new float[3];
    		View = new float[3];
    		Up = new float[3];
    		Right = new float[3];
    	}
    	// set position, focus points and the up vector
    	public void setCamera(float px, float py, float pz, float vx, float vy, float vz, float ux, float uy, float uz){    		
    		Position[0] = px; Position[1] = py; Position[2] = pz;
    		Center[0]   = vx; Center[1]   = vy; Center[2]   = vz;
    		Up[0]       = ux; Up[1]       = uy; Up[2]       = uz;
    		// calculate View and Right vectors
    		View = MinusVec(Center, Position); // pointing from the looking point toward the focus point
    		Right = NormVec(CrossVec(Up, View)); // right vector perpendicular to the the up and view vectors
    	}
    	
    	// calculate the new position if we rotate rot1 degree about the current up vectors (yaw) 
    	// and rot2 degree about the right vector (pitch)
    	public void SetRotation(float rot1, float rot2){    		
    		// rotate rot1 about Up
    		RotateCamera(rot1, Up[0], Up[1], Up[2]);
    		// recalculae Right
    		Right = NormVec(CrossVec(Up,View));
    		// rotate rot2 about Right
    		RotateCamera(rot2, Right[0], Right[1], Right[2]);
    		// recalculate Up
    		Up = NormVec(CrossVec(View,Right));
    		// recalculate Position
    		Position = MinusVec(Center,View);
    		// Get transformation matrix
    		cal_M();
    	}

    	// rotate current view vector an "angle" degree about an aribitrary vector [x,y,z] 
    	public void RotateCamera(float angle, float x, float y, float z){
    		float[] temp = new float[4];
    		float[] conjtemp = new float[4];
    		float[] quat_view = new float[4];
    		float[] result = new float[4];
    		// temp is the rotation quaternion
    		float deg = (float)((angle/180f)*Math.PI);
			float sinhtheta = (float)Math.sin(deg/2f);			
			temp[0] = x*sinhtheta;
			temp[1] = y*sinhtheta;
			temp[2] = z*sinhtheta;
			temp[3] = (float)Math.cos(deg/2f);
			// the conjugate rotation quaternion
			conjtemp[0] = -temp[0];
			conjtemp[1] = -temp[1];
			conjtemp[2] = -temp[2];
			conjtemp[3] =  temp[3];
			// convert view vector into quaternion
			quat_view[0] = View[0];
			quat_view[1] = View[1];
			quat_view[2] = View[2];
			quat_view[3] = 0.0f;
			// rotate by quaternion temp'*quat_view*temp
			result = MultQuat(MultQuat(temp, quat_view),conjtemp);
			// retrieve the new view vector from the resulted quaternion
			View[0] = result[0];
			View[1] = result[1];
			View[2] = result[2];
    	}
    	// calculate the rotation matrix from the current information
    	// equivalent to gluLookat
    	public void cal_M(){
    		float[] s, u, f;
    		// Calculate M for gluLookat
    		f = NormVec(MinusVec(Center,Position));    		    		
    		s = CrossVec(f, NormVec(Up));
    		u = CrossVec(s, f);
    		// Get transformation Matrix
    		M[0 ] =  s[0]; M[4 ] =  s[1]; M[8 ] =  s[2]; M[12] =  0.0f;
    		M[1 ] =  u[0]; M[5 ] =  u[1]; M[9 ] =  u[2]; M[13] =  0.0f;
    		M[2 ] = -f[0]; M[6 ] = -f[1]; M[10] = -f[2]; M[14] =  0.0f;
    		M[3 ] =  0.0f; M[7 ] =  0.0f; M[11] =  0.0f; M[15] =  1.0f;
    	}
    	// return cross product of two vectors
    	public float[] CrossVec(float[] A, float[] B){
    		float[] C = new float[3];
    		C[0] = A[1]*B[2] - A[2]*B[1];
    		C[1] = A[2]*B[0] - A[0]*B[2];
    		C[2] = A[0]*B[1] - A[1]*B[0];
    		return C;
    	}    	
    	// return the subtraction from two vectors
    	public float[] MinusVec(float[] A, float[] B){
    		float[] C = new float[3];
    		C[0] = A[0] - B[0];	C[1] = A[1] - B[1];	C[2] = A[2] - B[2];
    		return C;
    	}
    	// normalize a vector
    	public float[] NormVec(float[] A){
    		float[] C = new float[3];
    		float length = (float)Math.sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    		C[0] = A[0]/length; C[1] = A[1]/length; C[2] = A[2]/length;
    		return C;
    	}
    	// quaternion multiplication
    	public float[] MultQuat(float[] A, float[] B){
    		float[] C = new float[4];    		
    		C[0] = A[3]*B[0] + A[0]*B[3] + A[1]*B[2] - A[2]*B[1];
    		C[1] = A[3]*B[1] - A[0]*B[2] + A[1]*B[3] + A[2]*B[0];
    		C[2] = A[3]*B[2] + A[0]*B[1] - A[1]*B[0] + A[2]*B[3];
    		C[3] = A[3]*B[3] - A[0]*B[0] - A[1]*B[1] - A[2]*B[2];    		  						
    		return C;
    	}
    }
}

