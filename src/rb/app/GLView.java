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

import android.content.Context;
import android.content.res.Configuration;
import android.opengl.GLSurfaceView;
import android.view.KeyEvent;
import android.view.MotionEvent;
 
public class GLView extends GLSurfaceView {
    private static final String LOG_TAG = GLView.class.getSimpleName();
    private GLRenderer _renderer;
 
    private float _x = 0;
    private float _y = 0;
    private float _dist = 1.0f;
    private float old_zoom = 1.0f;
    
    boolean ismTouch = false;
    
    boolean isSensorCtrl = false;
    boolean current_paused = true;
    
    public GLView(Context context, GLObject _object) {
        super(context);
        setFocusableInTouchMode(true);
        _renderer = new GLRenderer(_object);
        
        Configuration c = getResources().getConfiguration(); 
        if(c.orientation == Configuration.ORIENTATION_PORTRAIT ) {
        	//_renderer.setcField(0);
            _renderer.setOrientation(true);
        } else if(c.orientation == Configuration.ORIENTATION_LANDSCAPE ){ 
        	//_renderer.setcField(1);
        	_renderer.setOrientation(false);
        }
        
        setRenderer(_renderer);
    }
    
    public boolean onTouchEvent(final MotionEvent event) {
    	switch (event.getAction() & MotionEvent.ACTION_MASK) {
     	case MotionEvent.ACTION_DOWN:
    		ismTouch = false;
            _x = event.getX();
            _y = event.getY();
            current_paused = _renderer.ispaused;
            _renderer.pause();
            break;
        
    	case MotionEvent.ACTION_UP:
    	case MotionEvent.ACTION_POINTER_UP:
    		_renderer.ispaused = current_paused;
    		break;
    	case MotionEvent.ACTION_MOVE:
    		// pass touchscreen data to the renderer    	
    		if (!ismTouch){
	            final float xdiff = (_x - event.getX());
	            final float ydiff = (_y - event.getY());
	            queueEvent(new Runnable() {
	                public void run() {
	                    _renderer.setPos(true,-xdiff/20.0f, ydiff/20.0f, 0.0f);
	                }
	            });
	            _x = event.getX();
	            _y = event.getY();
    		} else 
    		{
    			_x = event.getX(0)-event.getX(0);
                _y = event.getY(1)-event.getY(0);
                final float dist = (float)Math.sqrt(_x*_x + _y*_y);
                if (dist > 10f){
                queueEvent(new Runnable() {
	                public void run() {
	                	_renderer.zoom(dist/_dist*old_zoom);
	                }
                });}
    		}
    		break;
    	case MotionEvent.ACTION_POINTER_DOWN:
    		ismTouch = true;
    		_x = event.getX(0)-event.getX(0);
            _y = event.getY(1)-event.getY(0);
            _dist = (float) Math.sqrt(_x*_x + _y*_y);
            old_zoom = _renderer.scale_rat;
            current_paused = _renderer.ispaused;
            _renderer.pause();
            break;
    	}
        return true;
    }
    
    public boolean onKeyDown(int keyCode, KeyEvent event){
    	if (event.getAction() == MotionEvent.ACTION_DOWN) {
            switch (keyCode){
            case 24: // KEYCODE_VOLUME_UP
                //isSensorCtrl = !isSensorCtrl;
                _renderer.ispaused = true;
                _renderer.isContinuousRotation = false;
                _renderer.increase_ndframe(1f);
            	return true;
            case 25: // KEYCODE_VOLUME_DOWN
                //isSensorCtrl = !isSensorCtrl;
            	_renderer.ispaused = true;
            	_renderer.isContinuousRotation = false;
            	_renderer.increase_ndframe(-1f);
                return true;
            case 82: // KEYCODE_MENU
            	if ((_renderer._object.is2D()) || (_renderer._object.field_num>0))
            		_renderer.setcField(_renderer.getcField()+1);
            	else
            		// enable tilting
            		//isSensorCtrl = !isSensorCtrl;
            		// swap face rendering
            		_renderer.isFrontFace = !_renderer.isFrontFace;
            	return true;
            case 83: // KEYCODE_HOME
            	// do nothing!
            	return true;
            case 84: // KEYCODE_SEARCH
            	_renderer.ispaused = !_renderer.ispaused;
            	if (!_renderer._object.is2D())
            		_renderer.isFrontFace = !_renderer.isFrontFace;
            	return true;
            default:
            	return super.onKeyDown(keyCode, event);
            }          
        }
    	return true;
    }
    
    
    public boolean onTrackballEvent(MotionEvent event){ 
        float TBx = event.getX(); 
        float TBy = event.getY();
        if (event.getAction() == MotionEvent.ACTION_MOVE){ 
        	if ((TBx>=0)&(TBy<=0)) // zoom in if trackball is moving in the 2D "positive" direction        		
        		_renderer.zoomin();
        	else
        		_renderer.zoomout(); // and zoom out if not
        }
        if (event.getAction() == MotionEvent.ACTION_DOWN){ 
        	_renderer.zoomreset(); // reset to original status when users push the "pearl"
        	_renderer.setPos(true, 0.0f, 0.0f, 0.0f);
        	_renderer.ispaused = false;
        	_renderer.isContinuousRotation = true;
        }
        return true; 
    }
    
    public void setSensorParam(float x, float y, float z){
    	if (isSensorCtrl){
    		_renderer.setPos(false,-x/1.50f, -y/1.50f, 0.0f);
    	}
    }   
}


/*
//--------------------------------------------------------//
// This is the v1.5 compatible codes

package rb.app;

import android.content.Context;
import android.content.res.Configuration;
import android.opengl.GLSurfaceView;
import android.view.KeyEvent;
import android.view.MotionEvent;
 
public class GLView extends GLSurfaceView {
    private static final String LOG_TAG = GLView.class.getSimpleName();
    private GLRenderer _renderer;
 
    private float _x = 0;
    private float _y = 0;
    private float _dist = 1.0f;
    private float old_zoom = 1.0f;
    
    boolean ismTouch = false;
    
    boolean isSensorCtrl = false;
    boolean current_paused = true;
    
    public GLView(Context context, GLObject _object) {
        super(context);
        setFocusableInTouchMode(true);
        _renderer = new GLRenderer(_object);
        
        Configuration c = getResources().getConfiguration(); 
        if(c.orientation == Configuration.ORIENTATION_PORTRAIT ) {
        	//_renderer.setcField(0);
            _renderer.setOrientation(true);
        } else if(c.orientation == Configuration.ORIENTATION_LANDSCAPE ){ 
        	//_renderer.setcField(1);
        	_renderer.setOrientation(false);
        }
        
        setRenderer(_renderer);
    }
    
    public boolean onTouchEvent(final MotionEvent event) {
    	switch (event.getAction()){
     	case MotionEvent.ACTION_DOWN:
            _x = event.getX();
            _y = event.getY();
            current_paused = _renderer.ispaused;
            _renderer.pause();
            break;
        
    	case MotionEvent.ACTION_UP:    	
    		_renderer.ispaused = current_paused;
    		break;
    	case MotionEvent.ACTION_MOVE:
    		// pass touchscreen data to the renderer    	
    		final float xdiff = (_x - event.getX());
    		final float ydiff = (_y - event.getY());
    		queueEvent(new Runnable() {
    			public void run() {
    				_renderer.setPos(true,-xdiff/20.0f, ydiff/20.0f, 0.0f);
    			}
    		});
    		_x = event.getX();
    		_y = event.getY();    		
    		break;    	
    	}
        return true;
    }
    
    public boolean onKeyDown(int keyCode, KeyEvent event){
    	if (event.getAction() == MotionEvent.ACTION_DOWN) {
            switch (keyCode){
            case 24: // KEYCODE_VOLUME_UP
                //isSensorCtrl = !isSensorCtrl;
                _renderer.ispaused = true;
                _renderer.isContinuousRotation = false;
                _renderer.increase_ndframe(1f);
            	return true;
            case 25: // KEYCODE_VOLUME_DOWN
                //isSensorCtrl = !isSensorCtrl;
            	_renderer.ispaused = true;
            	_renderer.isContinuousRotation = false;
            	_renderer.increase_ndframe(-1f);
                return true;
            case 82: // KEYCODE_MENU
            	if ((_renderer._object.is2D()) || (_renderer._object.field_num>0))
            		_renderer.setcField(_renderer.getcField()+1);
            	else
            		// enable tilting
            		//isSensorCtrl = !isSensorCtrl;
            		// swap face rendering
            		_renderer.isFrontFace = !_renderer.isFrontFace;
            	return true;
            case 83: // KEYCODE_HOME
            	// do nothing!
            	return true;
            case 84: // KEYCODE_SEARCH
            	_renderer.ispaused = !_renderer.ispaused;
            	if (!_renderer._object.is2D())
            		_renderer.isFrontFace = !_renderer.isFrontFace;
            	return true;
            default:
            	return super.onKeyDown(keyCode, event);
            }          
        }
    	return true;
    }
    
    
    public boolean onTrackballEvent(MotionEvent event){ 
        float TBx = event.getX(); 
        float TBy = event.getY();
        if (event.getAction() == MotionEvent.ACTION_MOVE){ 
        	if ((TBx>=0)&(TBy<=0)) // zoom in if trackball is moving in the 2D "positive" direction        		
        		_renderer.zoomin();
        	else
        		_renderer.zoomout(); // and zoom out if not
        }
        if (event.getAction() == MotionEvent.ACTION_DOWN){ 
        	_renderer.zoomreset(); // reset to original status when users push the "pearl"
        	_renderer.setPos(true, 0.0f, 0.0f, 0.0f);
        	_renderer.ispaused = false;
        	_renderer.isContinuousRotation = true;
        }
        return true; 
    }
    
    public void setSensorParam(float x, float y, float z){
    	if (isSensorCtrl){
    		_renderer.setPos(false,-x/1.50f, -y/1.50f, 0.0f);
    	}
    }   
}
*/

