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

import java.util.List;

import android.app.Activity;
import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.util.Log;

public class RBVisualization extends Activity {
    private static String DEBUG_TAG = "RBVISUALIZATION";
    private GLView _glView;
    private GLObject _object;

    private SensorManager myManager; 
    private List<Sensor> sensors; 
    private Sensor accSensor; 
    private float oldX, oldY, oldZ = 0f; 
    
    @Override
    protected void onCreate(Bundle savedInstanceState){
        super.onCreate(savedInstanceState);        
        
        //_object = new GLObject();
        _object = RBActivity.mRbModel;        
        
        Bundle extras = getIntent().getExtras();
        
        if (!extras.getBoolean("isSweep")){			
        	if (!RBActivity.mRbSystem.is_custom_mesh_transform()){
	        	float[][] LT = RBActivity.mRbSystem.get_tranformation_data();
	        	float[][][] LTfunc_array = new float[1][LT.length][LT[0].length];
	        	LTfunc_array[0] = LT;        	
	        	_object.set_LTfunc(LTfunc_array);
        	} else {
        		Parameter[] p = new Parameter[1];
        		p[0] = RBActivity.mRbSystem.current_parameters.clone();
        		_object.mesh_transform_custom(p);
        	}

        	float[][][] truth_sol = RBActivity.mRbSystem.get_truth_sol();
        	        	
			if (RBActivity.mRbSystem.isReal)
				switch (RBActivity.mRbSystem.get_mfield()){
					case 1:
						_object.set_field_data(truth_sol[0][0]);
						break;
					case 2:
						if (truth_sol[0][0].length == RBActivity.mRbModel.node_num)					
							_object.set_field_data(truth_sol[0][0],truth_sol[1][0]);
						else
							_object.set_field_data(truth_sol[0][0],truth_sol[1][0],false);
						break;
					case 3:
						if (truth_sol[0][0].length == RBActivity.mRbModel.node_num)	
							_object.set_field_data(truth_sol[0][0],truth_sol[1][0],truth_sol[2][0]);
						else
							_object.set_field_data(truth_sol[0][0],truth_sol[1][0],truth_sol[2][0], false);
						break;
					case 4:
						_object.set_field_data(truth_sol[0][0],truth_sol[1][0],truth_sol[2][0],truth_sol[3][0]);
						break;
				}
			else
				switch (RBActivity.mRbSystem.get_mfield()){
				case 1:
					_object.set_field_data(truth_sol[0][0], truth_sol[0][1], truth_sol[0][2], false);
					break;
				}
        }else{
        	if (!RBActivity.mRbSystem.is_custom_mesh_transform()){
        		_object.set_LTfunc(_object.vLTfunc);    			
        	} else {
        		_object.mesh_transform_custom(RBActivity.mSweepParam);
        	}
        	float[][][] truth_sol = RBActivity.mRbSystem.get_sweep_truth_sol();
			
			if (RBActivity.mRbSystem.isReal)
				switch (RBActivity.mRbSystem.get_mfield()){
				case 1:
					_object.set_field_data(truth_sol[0][0]);
					break;
				case 2:
					_object.set_field_data(truth_sol[0][0],truth_sol[1][0]);
					break;
				case 3:
					_object.set_field_data(truth_sol[0][0],truth_sol[1][0],truth_sol[2][0]);
					break;
				case 4:
					_object.set_field_data(truth_sol[0][0],truth_sol[1][0],truth_sol[2][0],truth_sol[3][0]);
					break;
				}
			else {
				switch (RBActivity.mRbSystem.get_mfield()){
				case 1:
					_object.set_field_data(truth_sol[0][0], truth_sol[0][1], truth_sol[0][2], false);
					break;
				}
			}
        }
        
        // Set Sensor + Manager 
        myManager = (SensorManager)getSystemService(Context.SENSOR_SERVICE); 
        sensors = myManager.getSensorList(Sensor.TYPE_ACCELEROMETER);
        if(sensors.size() > 0){ 
        	accSensor = sensors.get(0); 
        }        
        
        _glView = new GLView(this, _object);
        setContentView(_glView);
    }
    
    private final SensorEventListener mySensorListener = new SensorEventListener(){ 
    	public void onSensorChanged(SensorEvent event){ 
    		// send data
    		_glView.setSensorParam(event.values[0], event.values[1], event.values[2]);
    		// update
    		oldX = event.values[0];
    		oldY = event.values[1];
    		oldZ = event.values[2];
    	} 
      
    	public void onAccuracyChanged(Sensor sensor, int accuracy) {} 
    }; 
    
    @Override 
    protected void onResume(){ 
    	super.onResume(); 
    	myManager.registerListener(mySensorListener, accSensor, SensorManager.SENSOR_DELAY_GAME);       
    } 
    
    @Override 
    protected void onStop(){      
    	myManager.unregisterListener(mySensorListener); 
    	super.onStop(); 
    }
    
}

