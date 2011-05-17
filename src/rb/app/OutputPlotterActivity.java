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

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.achartengine.ChartFactory;
import org.achartengine.GraphicalView;
import org.achartengine.chart.LineChart;
import org.achartengine.chart.PointStyle;
import org.achartengine.chart.XYChart;
import org.achartengine.model.XYMultipleSeriesDataset;
import org.achartengine.model.XYSeries;
import org.achartengine.renderer.XYMultipleSeriesRenderer;
import org.achartengine.renderer.XYSeriesRenderer;

import android.app.Activity;
import android.app.AlertDialog;
import android.app.Dialog;
import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.graphics.Color;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.text.InputType;
import android.util.Log;
import android.view.KeyEvent;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.view.MotionEvent;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.LinearLayout;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.TableRow.LayoutParams;

/**
 * Class that plots outputs as functions of time
 * using the AChartEngine library.
 * This class is modeled on the demo provided with
 * AChartEngine library.
 */
public class OutputPlotterActivity extends Activity {
	// String for log printing
	private static final String DEBUG_TAG = "OutputPlotter";
	
	//dialog IDs for graphing options
	static final int AXES_DIALOG_ID = 0;
	static final int PLOT_DIALOG_ID = 1;
	static final int LABEL_DIALOG_ID = 2;
	static final int WARNING_DIALOG_ID = 3;
	static final int INFO_DIALOG_ID = 4;
	
	// The Chart!
	private GraphicalView mChartView; 
	
	// An array that stores the colors used in plotting the outputs
	private int[] mColors;
	
	// An array storing which outputs will not be shown
	private boolean[] plotRemoved;
	
	//int to determine whether touch events trigger new sweep points
	private int add_sweep_pts_toggle; //0 for adding points, 1 for showing values, 2 for nothing
	
	//stores which series to show data for
	private int labeled_series;
	
	//stores which time point to show data for
	private int labeled_time_point;
	
	//stores the value at which a new sweep point is to be graphed
	private double progressXVal;
	
	//stores the interpolated y-value indicating where a new sweep point is to be graphed
	private double progressYVal;
	
	// A boolean to indicate whether we are doing a sweep over parameters
	public boolean isSweep;
	
	// Strings for labeling the plot
	public String title;
	public String xLabel;
	
	// The time step size
	public double dt;
	
	// The x-range
	public double xMin;
	public double xMax;
	
	// The y-range (determined dynamically by dataset)
	private double yMin;
	private double yMax;
	
	// Textviews for displaying x-range
	private TextView xMinView;
	private TextView xMaxView;
	
	// The number of time steps
	public int n_time_steps;
	
	// The number of outputs we will plot
	public int n_outputs;

	// The array of time step data
	double[] time_step_array;
	
	// The array of output data
	public double[][] RB_outputs_all_k;
	
	// The upper and lower bounds for each output
	public double[][] RB_outputs_LB;
	public double[][] RB_outputs_UB;
	
	// saved information from last screen
	private Bundle extras;
	
	private int current_output;
	
	/**
	 * Constructor.
	 */
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.output_plotter_layout);
		
		extras = getIntent().getExtras();
		
		isSweep = extras.getBoolean("isSweep");
		title = extras.getString("title");
		dt = extras.getDouble("dt");
		xMin = extras.getDouble("xMin");
		xMax = extras.getDouble("xMax");
		xLabel = extras.getString("xLabel");
		n_time_steps = extras.getInt("n_time_steps");
		n_outputs = extras.getInt("n_outputs");
		//MW
		//yMin = yMax = 0;
		add_sweep_pts_toggle = 2;
		labeled_time_point = -1;
		progressXVal = xMin-1;
		progressYVal = 0;
		labeled_series = 0;
		
		time_step_array = new double[n_time_steps];
		for(int time_step=0; time_step<n_time_steps; time_step++) {
			time_step_array[time_step] = xMin + time_step*dt;
		}
		
		// Get the output and error bound data
		RB_outputs_all_k = new double[n_outputs][];
		double[][] RB_output_error_bounds_all_k = new double[n_outputs][];
		for(int i=0; i<n_outputs; i++) {
			RB_outputs_all_k[i] = extras.getDoubleArray("output_data_" + i);
			RB_output_error_bounds_all_k[i] = extras.getDoubleArray("output_bound_" + i);
		}
		
		// Create the UB and LB array based on
		// RB_outputs_all_k and RB_output_error_bounds_all_k
		
		RB_outputs_LB = new double[n_outputs][n_time_steps];
		RB_outputs_UB = new double[n_outputs][n_time_steps];
		
		for(int i=0; i<n_outputs; i++) {
			for(int time_step=0; time_step<n_time_steps; time_step++) {
				RB_outputs_LB[i][time_step] =
					RB_outputs_all_k[i][time_step] - RB_output_error_bounds_all_k[i][time_step];
				RB_outputs_UB[i][time_step] =
					RB_outputs_all_k[i][time_step] + RB_output_error_bounds_all_k[i][time_step];
			}
		}
		
		plotRemoved = new boolean[n_outputs];
		for(int i=0; i<n_outputs; i++){ plotRemoved[i]=true;}
		plotRemoved[0] = false;
		current_output = 0;
		
	    LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
	    mChartView = execute(this);
	    layout.addView(mChartView, new LayoutParams
	    		(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
	    
		// Add the legend
	    TableLayout legendLayout = (TableLayout) findViewById(R.id.legendLayout);
	    
	    int n_labels_per_row = 5; // The number of legend labels on each row
	    
    	TableRow row = new TableRow(this); // Build an initial row
		row.setLayoutParams(new LayoutParams(
				LayoutParams.FILL_PARENT,
				LayoutParams.FILL_PARENT) );
	    
	    int count = 0;
	    while(count<n_outputs) {
	    	TextView legend_i = new TextView(this);
	    	if(count%n_labels_per_row==0) legend_i.setText("   " + "output " + (count+1) + "   ");
	    	else legend_i.setText("output " + (count+1) + "   ");
	    	legend_i.setTextColor(mColors[count]);
	    	legend_i.setLayoutParams(new LayoutParams(
                    LayoutParams.FILL_PARENT,
                    LayoutParams.WRAP_CONTENT));
	    	
	    	// set background color of the label
	    	legend_i.setBackgroundColor(Color.WHITE);
	    	
	    	row.addView(legend_i);
	    	
	    	// set background color of the row
	    	row.setBackgroundColor(Color.WHITE);
	    	
	    	count++;
	    	
	    	// Add the current row and construct a new row
	    	if((count+1)%n_labels_per_row == 0) {
			    legendLayout.addView(row, new TableLayout.LayoutParams(
		                LayoutParams.FILL_PARENT,
		                LayoutParams.WRAP_CONTENT));
			    
		    	row = new TableRow(this);
		    	row.setLayoutParams(new LayoutParams(
						LayoutParams.FILL_PARENT,
						LayoutParams.FILL_PARENT) );
	    	}
	    }
    	legendLayout.addView(row, new TableLayout.LayoutParams(
            	LayoutParams.FILL_PARENT,
            	LayoutParams.WRAP_CONTENT));
	   
		// Add listener to the Visualize button
	    if(!isSweep) {
			Button visButton = (Button) findViewById(R.id.unsteadyVisButton);
			visButton.setOnClickListener(new View.OnClickListener() {
	
				public void onClick(View view) {
					// Next create the bundle and initialize it
					Bundle bundle = new Bundle();	
					bundle.putBoolean("isSweep", false);
									
					Intent intent = new Intent(OutputPlotterActivity.this, RBVisualization.class);
					intent.putExtras(bundle);
					OutputPlotterActivity.this.startActivity(intent);
							
				}
			});
	    } else{ // sweep case
	    	Button visButton = (Button) findViewById(R.id.unsteadyVisButton);
			visButton.setOnClickListener(new View.OnClickListener() {
	
				public void onClick(View view) {
					// Next create the bundle and initialize it
					Bundle bundle = new Bundle();
					bundle.putBoolean("isSweep", true);
					
					Intent intent = new Intent(OutputPlotterActivity.this, RBVisualization.class);
					intent.putExtras(bundle);
					OutputPlotterActivity.this.startActivity(intent);
					
				}
			});	    	
	    }

	}
	
	/**
	 * Builds an XY multiple dataset using the provided values.
	 * @param titles the series titles
	 * @param xValues the values for the X axis
	 * @param yValues the values for the Y axis
	 * @return the XY multiple dataset
	 */
	protected XYMultipleSeriesDataset buildDataset(String[] titles,
			                                       List<double[]> xValues,
			                                       List<double[]> yValues) {
		XYMultipleSeriesDataset dataset = new XYMultipleSeriesDataset();
		int num_series = 3 * n_outputs;
		
		yMin = yMax = (yValues.get(0))[0];
		for (int i = 0; i < num_series; i++) {
			XYSeries series = new XYSeries(titles[i]);
			double[] xV = xValues.get(i);
			double[] yV = yValues.get(i);
			int seriesLength = xV.length;
			
			for (int k = 0; k < seriesLength; k++) {
				//determine which values not to plot based on given xMin & xMax value
				//do not add values for hidden plots
				if((xV[k]>=xMin)&&(xV[k]<=xMax)&& !plotRemoved[i/3]) series.add(xV[k], yV[k]);
				if(yV[k]<yMin && !plotRemoved[i/3]) yMin = yV[k];
				if(yV[k]>yMax && !plotRemoved[i/3]) yMax = yV[k];
			}
			dataset.addSeries(series);
		}
		
		//an additional series for displaying what sweep values are to be computed with a vertical line
		XYSeries series = new XYSeries("");
		series.add(progressXVal, yMin);
		series.add(progressXVal, yMax);
		if(!isSweep || progressYVal>yMax || progressYVal < yMin) progressYVal = yMin; //don't want to stretch graph
		series.add(progressXVal, progressYVal);
		dataset.addSeries(series);
		
		return dataset;
	}

	/**
	 * Builds an XY multiple series renderer.
	 * @param colors the series rendering colors
	 * @param styles the series point styles
	 * @return the XY multiple series renderers
	 */
	protected XYMultipleSeriesRenderer buildRenderer() {
		XYMultipleSeriesRenderer renderer = new XYMultipleSeriesRenderer();		
		
		mColors = new int[n_outputs];
		for (int i = 0; i < n_outputs; i++) {
			XYSeriesRenderer r_out = new XYSeriesRenderer();
			XYSeriesRenderer r_LB  = new XYSeriesRenderer();
			XYSeriesRenderer r_UB  = new XYSeriesRenderer();
			
			// Set colors
			switch(i) {
			case 0:
				mColors[0] = Color.BLUE;
				r_out.setColor(Color.BLUE); // Output
				r_LB.setColor(Color.BLUE); // LB
				r_UB.setColor(Color.BLUE); // UB
				break;
			case 1:
				mColors[1] = Color.RED;
				r_out.setColor(Color.RED);
				r_LB.setColor(Color.RED);
				r_UB.setColor(Color.RED);
				break;
			case 2:
				mColors[2] = Color.GREEN;
				r_out.setColor(Color.GREEN);
				r_LB.setColor(Color.GREEN);
				r_UB.setColor(Color.GREEN);
				break;
			case 3:
				mColors[3] = Color.CYAN;
				r_out.setColor(Color.CYAN);
				r_LB.setColor(Color.CYAN);
				r_UB.setColor(Color.CYAN);
				break;
			case 4:
				mColors[4] = Color.YELLOW;
				r_out.setColor(Color.YELLOW);
				r_LB.setColor(Color.YELLOW);
				r_UB.setColor(Color.YELLOW);
				break;
			case 5:
				mColors[5] = Color.MAGENTA;
				r_out.setColor(Color.MAGENTA);
				r_LB.setColor(Color.MAGENTA);
				r_UB.setColor(Color.MAGENTA);
				break;
			case 6:
				mColors[6] = Color.GRAY;
				r_out.setColor(Color.GRAY);
				r_LB.setColor(Color.GRAY);
				r_UB.setColor(Color.GRAY);
				break;
			case 7:
				mColors[7] = Color.BLACK;
				r_out.setColor(Color.BLACK);
				r_LB.setColor(Color.BLACK);
				r_UB.setColor(Color.BLACK);
				break;
			default:
				Random numGen = new Random();
				mColors[i] = 
					Color.argb(255, numGen.nextInt(256), numGen.nextInt(256), numGen.nextInt(256));
				r_out.setColor(mColors[i]);
				r_LB.setColor(mColors[i]);
				r_UB.setColor(mColors[i]);
				break;
			}
			
			r_out.setPointStyle(PointStyle.X);    // output
			r_LB.setPointStyle(PointStyle.POINT); // LB
			r_UB.setPointStyle(PointStyle.POINT); // UB

			renderer.addSeriesRenderer(r_out);
			renderer.addSeriesRenderer(r_LB);
			renderer.addSeriesRenderer(r_UB);
		}
		
		XYSeriesRenderer r_add = new XYSeriesRenderer();
		r_add.setColor(Color.BLACK);
		r_add.setPointStyle(PointStyle.POINT);
		renderer.addSeriesRenderer(r_add);
		
		return renderer;
	}

	/**
	 * Sets a few of the series renderer settings.
	 * @param renderer the renderer to set the properties to
	 * @param title the chart title
	 * @param xTitle the title for the X axis
	 * @param yTitle the title for the Y axis
	 * @param xMin the minimum value on the X axis
	 * @param xMax the maximum value on the X axis
	 * @param yMin the minimum value on the Y axis
	 * @param yMax the maximum value on the Y axis
	 * @param axesColor the axes color
	 * @param labelsColor the labels color
	 */
	protected void setChartSettings(XYMultipleSeriesRenderer renderer, String title, String xLabel,
			double xMin, double xMax) {
		renderer.setChartTitle(title);
		renderer.setXTitle(xLabel);
		renderer.setXAxisMin(xMin);
		renderer.setXAxisMax(xMax);
		renderer.setAxesColor(Color.BLACK);
		renderer.setLabelsColor(Color.BLACK);
		
		//MW
		renderer.setDisplayChartValues(true);
		
		// Set the approximate number of "ticks" on the x and y axes
		renderer.setXLabels(10);
		renderer.setYLabels(10);
		renderer.setChartValuesTextSize(16);
		
		// Set the background color
		renderer.setApplyBackgroundColor(true);
		renderer.setBackgroundColor(Color.WHITE);
		
		// Turn off the legend
		renderer.setShowLegend(false);
	}

	/**
	 * Executes the chart demo.
	 * @param context the context
	 * @return the built intent
	 */
	public GraphicalView execute(Context context) {
		
		String[] labels = new String[3*n_outputs];
		for(int i=0; i<n_outputs; i++) {
			labels[3*i]   = "output " + (i+1);
			labels[3*i+1] = "LB " + (i+1);
			labels[3*i+2] = "UB " + (i+1);
		}
		
		List<double[]> t = new ArrayList<double[]>();
		List<double[]> output_values = new ArrayList<double[]>();		
		
		// Add each output and corresponding UB and LB 
		for(int i=0; i<n_outputs; i++) {
			t.add(time_step_array);
			output_values.add(RB_outputs_all_k[i]);
			
			t.add(time_step_array);
			output_values.add(RB_outputs_LB[i]);
			
			t.add(time_step_array);
			output_values.add(RB_outputs_UB[i]);
		}

		XYMultipleSeriesRenderer renderer = buildRenderer();
		setChartSettings(renderer, "", "", xMin, xMax+(xMax-xMin)*.05);

		//return ChartFactory.getLineChartView(context, buildDataset(labels, t, output_values), renderer);
		boolean sweepInProgress = !(progressYVal==yMin);
		labeled_series = sweepInProgress?n_outputs*3:labeled_series;
		
		SingleLabelChart theChart = new SingleLabelChart(buildDataset(labels, t, output_values),
				renderer, labeled_series, labeled_time_point, sweepInProgress);
		return new GraphicalView(context, theChart);
	}

	
	public boolean onTouchEvent(final MotionEvent event) {
		if (event.getAction() == MotionEvent.ACTION_DOWN) {
			float X = event.getX();
			float Y = event.getY();		
			if(add_sweep_pts_toggle==0)add_sweep_point(X/mChartView.getWidth());
			else if(add_sweep_pts_toggle==1)showPoint(X/mChartView.getWidth());
			Log.d(DEBUG_TAG,"touch point (x) = " + (event.getX()-15)/290);
		}
		if (event.getAction() == MotionEvent.ACTION_UP) {
		}
		return true;
	}
	
	/**
	 * Adds a new point in non-time-dependent problems on where a touch event has occurred
	 */
	public void add_sweep_point(float xpos){
		if (isSweep) {
			//make sure xpos is within the range (0.01, 0.99)
			xpos = xpos<0.01f?0.01f:xpos;
			xpos = xpos>0.99f?0.99f:xpos;
			
			
			//xrange is the current range of the graph's x-values (determined by xMin and xMax)
			//xcurrent is the x-value of the touch point
			double xrange = (time_step_array[n_time_steps-1]-time_step_array[0])*
							((xMax-xMin)/(extras.getDouble("xMax")-extras.getDouble("xMin")));
			
			double xcurrent = time_step_array[0] + 
							 (time_step_array[n_time_steps-1]-time_step_array[0])*
							 (xMin-extras.getDouble("xMin"))/(extras.getDouble("xMax")-extras.getDouble("xMin"))
							 + xrange*(xpos);
			
			
			//find time step value closest to y-axis
			int leftmostStep = 0;
			while(time_step_array[leftmostStep]<xMin){
				leftmostStep++;
			}
			double xAxisValue = ((xMax-xMin)*xpos)+xMin;
			//now, adjusting for the current leftmost point, find the smallest time step greater than the touch point
			int closestTimestep = 0; 
			while(time_step_array[closestTimestep+leftmostStep]<xAxisValue){
				closestTimestep++;
			}
			
			//take the average of the values of the first visible dataset to get the y axis value for showing
			//progress label
			int firstVisible = 0;
			while(plotRemoved[firstVisible]) firstVisible++;
			if(firstVisible < n_outputs){ //there is a visible output
				
				int ts = closestTimestep+leftmostStep;
				//before computing the value, set progress values to indicate the computation is taking place
				progressXVal = xcurrent;
				progressYVal = RB_outputs_all_k[firstVisible][ts-1]+
				              ((RB_outputs_all_k[firstVisible][ts]
				              -RB_outputs_all_k[firstVisible][ts-1])/
				              (time_step_array[ts]-time_step_array[ts-1]))
				              *(xcurrent-time_step_array[ts-1]);
				
				//disable new sweep points until computation is complete
				add_sweep_pts_toggle = 2;
				
				//repaint graph
				LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
				mChartView = execute(this);
				layout.removeAllViews();
				layout.addView(mChartView, new LayoutParams
						(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
		    
				//start separate thread to actually solve for sweep point and graph it
				SolveThread st = new SolveThread();
				st.start();
			}
		    
		}
	}
	
	/**
	 * Displays information on an individual user-selected point on the graph
	 */
	public void showPoint(float xpos)
	{
		//make sure xpos is within the range (0.01, 0.99)
		xpos = xpos<0.01f?0.01f:xpos;
		xpos = xpos>0.99f?0.99f:xpos;
		
		double xAxisValue = ((xMax-xMin)*xpos)+xMin;

		//find time step value closest to y-axis
		int leftmostStep = 0;
		while(time_step_array[leftmostStep]<xMin){
			leftmostStep++;
		}

		//now, adjusting for the current leftmost point, find the smallest time step greater than the touch point
		int closestTimestep = 0; 
		while(time_step_array[closestTimestep+leftmostStep]<xAxisValue){
			closestTimestep++;
		}

		//again adjusting for the current leftmost point, take either the previous time step or the one before it
		if(closestTimestep>0) 
			closestTimestep = time_step_array[closestTimestep+leftmostStep]-xAxisValue>
							  xAxisValue-time_step_array[closestTimestep+leftmostStep-1]?
							  closestTimestep-1:closestTimestep;
		
		//labeled_time_point is given to the constructor of SingleLabelChart to determine
		//which values get labeled instead of all of them (in the default XYChart)
		labeled_time_point = closestTimestep*2; 
		Log.d(DEBUG_TAG, "Closest timestep: "+closestTimestep);
		
		
		//repaint graph
		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
		mChartView = execute(this);
		layout.removeAllViews();
		layout.addView(mChartView, new LayoutParams
				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
		
		
		
	}
	
	
	/**
	 * This function takes care of constructing the dialogs that pop up.
	 */
	protected Dialog onCreateDialog(int id) {

		Dialog dialog;

		switch (id)
		{
		case AXES_DIALOG_ID:
			dialog = new Dialog(this);
			dialog.setContentView(R.layout.axes_dialog);
			dialog.setTitle("Coordinate Axes Options");
			dialog.setCancelable(false);
			
			//initialize text fields
			TextView xMinLabel = (TextView) dialog.findViewById(R.id.x_min_label);
			xMinLabel.setText("Minimum x-value");
			xMinView = (EditText) dialog.findViewById(R.id.x_min_entry);
			xMinView.setText(String.valueOf(xMin));
			xMinView.setInputType(InputType.TYPE_CLASS_NUMBER | InputType.TYPE_NUMBER_FLAG_DECIMAL
					| InputType.TYPE_NUMBER_FLAG_SIGNED);
			
			TextView xMaxLabel = (TextView) dialog.findViewById(R.id.x_max_label);
			xMaxLabel.setText("Maximum x-value");
			xMaxView = (EditText) dialog.findViewById(R.id.x_max_entry);
			xMaxView.setText(String.valueOf(xMax));
			xMaxView.setInputType(InputType.TYPE_CLASS_NUMBER | InputType.TYPE_NUMBER_FLAG_DECIMAL
					| InputType.TYPE_NUMBER_FLAG_SIGNED);
			
			//button for resetting back to default xMin and xMax
			Button resetButton = (Button) dialog.findViewById(R.id.reset_default);
			resetButton.setOnClickListener(new View.OnClickListener() {
				
				public void onClick(View view)
				{
					//change local values and re-plot
					xMin = extras.getDouble("xMin");
					xMax = extras.getDouble("xMax");
					LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
	    		    mChartView = execute(getApplicationContext());
	    		    layout.removeAllViews();
	    		    layout.addView(mChartView, new LayoutParams
	    		    		(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
	    		    
	    		    dismissDialog(AXES_DIALOG_ID);
					removeDialog(AXES_DIALOG_ID);
				}
			});
			
			//handle changing local values to user-submitted xMin and xMax values on clicking the ok button
			Button okButton = (Button) dialog.findViewById(R.id.axes_okButton);
			okButton.setOnClickListener(new View.OnClickListener() {

				public void onClick(View view){
					//determine if value in minimum input field is within acceptable range
					String userMinString = xMinView.getText().toString();
					double userMin;
					try{
						userMin = Double.parseDouble(userMinString);
					}
					catch(NumberFormatException e){ 
						//if user submits non-double, default value is out of bounds to trigger toast
						userMin = extras.getDouble("xMin")-1;
					}
					//determine if value in maximum input field is within acceptable range
					String userMaxString = xMaxView.getText().toString();
					double userMax;
					try{
						userMax = Double.parseDouble(userMaxString);
					}
					catch(NumberFormatException e){ 
						//if user submits non-double, default value is out of bounds to trigger toast
						userMax = extras.getDouble("xMax")+1;
					}
					//make sure the given max value is greater than the given min, and that they are both
					//within original values
					if((userMax>userMin)&&(userMax<=extras.getDouble("xMax"))&&(userMin>=extras.getDouble("xMin")))
					{
						//change local values
						xMin = userMin;
						xMax = userMax;
						//add data points at new endpoints
						//add_sweep_point(0.011f);
						//add_sweep_point(1f);
						//re-plot
						LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
		    		    mChartView = execute(getApplicationContext());
		    		    layout.removeAllViews();
		    		    layout.addView(mChartView, new LayoutParams
		    		    		(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));	
					}
					else{ Toast.makeText(getApplicationContext(), "Value out of bounds", Toast.LENGTH_SHORT).show();}
					dismissDialog(AXES_DIALOG_ID);
					removeDialog(AXES_DIALOG_ID);
					
				}
			});
			break;
		case PLOT_DIALOG_ID:
			String[] plotStrings = new String[n_outputs];
			
			for (int i = 0; i < n_outputs; i++) {
					plotStrings[i] = "output " +(i+1);
			}
			
			boolean[] plotBooleans = new boolean[n_outputs];
			for(int i = 0; i < n_outputs; i++){
				if(plotRemoved[i] == true) plotBooleans[i] = false;
				else plotBooleans[i] = true;
			}
			
			
			AlertDialog.Builder builder = new AlertDialog.Builder(OutputPlotterActivity.this);
			builder.setTitle("Enable/Disable Plots");
			builder.setMultiChoiceItems(plotStrings, plotBooleans, new DialogInterface.OnMultiChoiceClickListener() {
			    public void onClick(DialogInterface dialog, int item, boolean checked) {
			    	//reset colors appropriately
			    	if (checked){
			    		plotRemoved[item] = false;
			    	}
			    	else{
			    		plotRemoved[item] = true;
			    	}	
			    }
			});
			builder.setNeutralButton("OK", new DialogInterface.OnClickListener(){
				public void onClick(DialogInterface dialog, int item)
				{	
					//re-plot the graph
		    		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
		    		mChartView = execute(getApplicationContext());
		    		layout.removeAllViews();
		    		layout.addView(mChartView, new LayoutParams
		    				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
		    	}
			});
			dialog = builder.create();
			break;
		case LABEL_DIALOG_ID:
			AlertDialog.Builder builder2 = new AlertDialog.Builder(OutputPlotterActivity.this);
			if(isSweep){
				String[] labelStrings = new String[3];
				labelStrings[0]="Add sweep points";
				labelStrings[1]="Show data labels";
				labelStrings[2]="Do nothing";
				
				boolean[] labelBooleans = new boolean[3];
				labelBooleans[add_sweep_pts_toggle]=true;
				
				builder2.setTitle("Set on touch behavior");
				builder2.setMultiChoiceItems(labelStrings, labelBooleans, new DialogInterface.OnMultiChoiceClickListener() {
				    public void onClick(DialogInterface dialog, int item, boolean checked) {
			    		if(item==0){
				    		add_sweep_pts_toggle = 0; //now on touch, sweep points will be added
			    		}
			    		if(item==1){
				    		add_sweep_pts_toggle = 1; //now on touch, data values will show
				    		
				    		//determine whether more than one output is showing: if this is the case,
				    		//recommend that the user disable all but one output
							int visibleCount = 0;
							for(int i=0; i<n_outputs; i++)
							{
								if(plotRemoved[i]==false) visibleCount++;
							}
							if(visibleCount>1)
							{
								showDialog(WARNING_DIALOG_ID);
							}
				    	}
			    		if(item==2){
			    			add_sweep_pts_toggle = 2;
			    		}
			    		labeled_time_point = -1; //remove the label
				    	
				    	//re-plot the graph
			    		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
			    		mChartView = execute(getApplicationContext());
			    		layout.removeAllViews();
			    		layout.addView(mChartView, new LayoutParams
			    				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
			    		//have to remove dialog because this is not automatic for a multi-choice dialog
			    		dismissDialog(LABEL_DIALOG_ID);
			    		removeDialog(LABEL_DIALOG_ID);
				    }
				});
			}
			else{
				String[] labelStrings = new String[2];
				labelStrings[0]="Show data labels";
				labelStrings[1]="Do nothing";
				
				boolean[] labelBooleans = new boolean[2];
				if(add_sweep_pts_toggle==1)labelBooleans[0]=true;
				else labelBooleans[1]=true;
				
				builder2.setTitle("Set on touch behavior");
				builder2.setMultiChoiceItems(labelStrings, labelBooleans, new DialogInterface.OnMultiChoiceClickListener() {
				    public void onClick(DialogInterface dialog, int item, boolean checked) {
			    		if(item==1){
				    		add_sweep_pts_toggle = 2; //now on touch, sweep points will be added
			    		}
			    		if(item==0){
				    		add_sweep_pts_toggle = 1; //now on touch, data values will show
				    		
				    		//determine whether more than one output is showing: if this is the case,
				    		//recommend that the user disable all but one output
							int visibleCount = 0;
							for(int i=0; i<n_outputs; i++)
							{
								if(plotRemoved[i]==false) visibleCount++;
							}
							if(visibleCount>1)
							{
								showDialog(WARNING_DIALOG_ID);
							}
				    	}
			    		labeled_time_point = -1; //remove the label
				    	
				    	//re-plot the graph
			    		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
			    		mChartView = execute(getApplicationContext());
			    		layout.removeAllViews();
			    		layout.addView(mChartView, new LayoutParams
			    				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
			    		//have to remove dialog because this is not automatic for a multi-choice dialog
			    		dismissDialog(LABEL_DIALOG_ID);
			    		removeDialog(LABEL_DIALOG_ID);
				    }
				});
			}
			
			dialog = builder2.create();
			break;
		case WARNING_DIALOG_ID:
			String[] plotStrings2 = new String[n_outputs];
			
			for (int i = 0; i < n_outputs; i++) {
					plotStrings2[i] = "output " +(i+1);
			}
			
			AlertDialog.Builder builder3 = new AlertDialog.Builder(OutputPlotterActivity.this);
			builder3.setTitle("Show labels for which output?");
			builder3.setItems(plotStrings2, new DialogInterface.OnClickListener() {
			    public void onClick(DialogInterface dialog, int item) {

			    	labeled_series = item;
			    		
			    	//re-plot the graph
			    	LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
			    	mChartView = execute(getApplicationContext());
			    	layout.removeAllViews();
			    	layout.addView(mChartView, new LayoutParams
			    			(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
			    }
			});
			dialog = builder3.create();
			break;
		case INFO_DIALOG_ID:
			AlertDialog.Builder infoBuilder = new AlertDialog.Builder(OutputPlotterActivity.this);
			infoBuilder.setTitle("Current Parameters");
			
			String message = "Online N = " + RBActivity.mOnlineNForGui + "\n\n" + "Parameters: \n\n";
			for(int i=0; i<RBActivity.mCurrentParamForGUI.getNEntries(); i++)
			{
				if(isSweep && Integer.parseInt(xLabel) == (i+1)) message = message + (i+1) + ": " + "Sweep\n";
				else message = message + (i+1) + ": " + RBActivity.mCurrentParamForGUI.getEntry(i) + "\n";
			}
			message = message + "\nChange values?";
			infoBuilder.setMessage(message);
			infoBuilder.setPositiveButton("Yes", new DialogInterface.OnClickListener() {
			           public void onClick(DialogInterface dialog, int id) {
			        	   onBackPressed();
			           }
			       });
			infoBuilder.setNegativeButton("No", null);
			dialog = infoBuilder.create();
			break;
		default:
			dialog = null;
			break;
		}
		return dialog;

	}
	

	
	//populate main menu
	public boolean onPrepareOptionsMenu(Menu menu)
	{
		menu.clear();
		MenuInflater inflater = getMenuInflater();
		inflater.inflate(R.menu.plotter_main_menu, menu);
		return true;
	}
	
	//handle item selection in main menu
	public boolean onOptionsItemSelected(MenuItem item)
	{
		switch (item.getItemId())
		{
		case R.id.plot_info:
			showDialog(INFO_DIALOG_ID);
			return true;
		case R.id.change_axes:
			showDialog(AXES_DIALOG_ID);
			return true;
		case R.id.en_dis_plots:
			showDialog(PLOT_DIALOG_ID);
			return true;
		case R.id.on_touch:
			showDialog(LABEL_DIALOG_ID);
			return true;
		}
		return false;
	}
	
	//solves for new sweep values separately from UI thread, allowing UI to indicate that a solve
	//is being performed
	private class SolveThread extends Thread
	{
		
		public void run()
		{
			int icount = 0; 
			double xcurrent = progressXVal;
			while ((icount < n_time_steps) && (xcurrent > time_step_array[icount])){
				icount++;
			}
			//xleft and xright are the numbers of the time steps immediately to the left and right
			//of the touch point, respectively
			int xleft = icount-1;
			int xright = icount;
			
			//double max_time_step = time_step_array[xleft+1]-time_step_array[xleft+1];
			double max_time_step = 0;
			int i_max_step = xleft+1;
			for (int i = xleft; i < xright; i++){
				double current_time_step = time_step_array[i+1]-time_step_array[i];
				if (max_time_step<current_time_step){
					max_time_step = current_time_step;
					i_max_step = i+1;
				}					
			}				
			
			RBActivity.mCurrentParamForGUI.setEntry(RBActivity.mSweepIndex, xcurrent);
			RBActivity.mRbSystem.setCurrentParameters(RBActivity.mCurrentParamForGUI);				
			RBActivity.mRbSystem.RB_solve(RBActivity.mOnlineNForGui);
			
			double[][] new_RB_outputs_all_k = new double[n_outputs][n_time_steps+1];
			double[][] new_RB_outputs_LB = new double[n_outputs][n_time_steps+1];
			double[][] new_RB_outputs_UB = new double[n_outputs][n_time_steps+1];
			double[] new_time_step_array = new double[n_time_steps+1];				

			new_time_step_array[i_max_step] = xcurrent;
			
			for (int j = 0; j < n_time_steps; j++){
				int k = (j<i_max_step)?j:j+1;
				new_time_step_array[k] = time_step_array[j];
				for (int i = 0; i < n_outputs; i++){
					new_RB_outputs_all_k[i][k] = RB_outputs_all_k[i][j];
					new_RB_outputs_LB[i][k] = RB_outputs_LB[i][j];
					new_RB_outputs_UB[i][k] = RB_outputs_UB[i][j];
				}
			}

			
			if (RBActivity.mRbSystem.isReal)
				for(int n=0; n<n_outputs; n++) {
					new_RB_outputs_all_k[n][i_max_step] = RBActivity.mRbSystem.RB_outputs[n];
					double OutputBound = RBActivity.mRbSystem.RB_output_error_bounds[n];
					new_RB_outputs_LB[n][i_max_step] = new_RB_outputs_all_k[n][i_max_step] - OutputBound;
					new_RB_outputs_UB[n][i_max_step] = new_RB_outputs_all_k[n][i_max_step] + OutputBound;
				}
			else
				for(int n=0; n<n_outputs/2; n++) {
					new_RB_outputs_all_k[n][i_max_step] = RBActivity.mRbSystem.get_RB_output(n,true);
					double OutputBound = RBActivity.mRbSystem.get_RB_output_error_bound(n, true);
					new_RB_outputs_LB[n][i_max_step] = new_RB_outputs_all_k[n][i_max_step] - OutputBound;
					new_RB_outputs_UB[n][i_max_step] = new_RB_outputs_all_k[n][i_max_step] + OutputBound;
					
					new_RB_outputs_all_k[n+ n_outputs/2][i_max_step] = RBActivity.mRbSystem.get_RB_output(n,false);
					OutputBound = RBActivity.mRbSystem.get_RB_output_error_bound(n, false);
					new_RB_outputs_LB[n+ n_outputs/2][i_max_step] = new_RB_outputs_all_k[n + n_outputs/2][i_max_step] - OutputBound;
					new_RB_outputs_UB[n+ n_outputs/2][i_max_step] = new_RB_outputs_all_k[n + n_outputs/2][i_max_step] + OutputBound;					
				}
			
			RB_outputs_all_k = new_RB_outputs_all_k;
			RB_outputs_LB = new_RB_outputs_LB;
			RB_outputs_UB = new_RB_outputs_UB;
			time_step_array = new_time_step_array;
			
			n_time_steps++;
			
			//hide progress vertical line
			progressXVal = xMin-1;
			progressYVal = yMin;
			
			//re-enable adding sweep points
			add_sweep_pts_toggle = 0;
			
			//update graph
			handler.sendEmptyMessage(0);
		}
		
		private Handler handler = new Handler()
		{
			public void handleMessage(Message msg){
				LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
			    mChartView = execute(getApplicationContext());
			    layout.removeAllViews();
			    layout.addView(mChartView, new LayoutParams
			    		(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
			}
		};
	}

	public boolean onKeyDown(int keyCode, KeyEvent event){
    	if (event.getAction() == MotionEvent.ACTION_DOWN) {
            switch (keyCode){
            case 24: // KEYCODE_VOLUME_UP
            {
            	current_output += 1;
            	if (current_output == n_outputs)
            		current_output = 0;
            	for(int i=0; i<n_outputs; i++){ plotRemoved[i]=true;}
        		plotRemoved[current_output] = false;
        		
        		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
	    		mChartView = execute(getApplicationContext());
	    		layout.removeAllViews();
	    		layout.addView(mChartView, new LayoutParams
	    				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
            
	    		return true;
            }
            case 25: // KEYCODE_VOLUME_DOWN
            {
            	current_output -= 1;
            	if (current_output == -1)
            		current_output = n_outputs-1;
            	for(int i=0; i<n_outputs; i++){ plotRemoved[i]=true;}
        		plotRemoved[current_output] = false;

        		LinearLayout layout = (LinearLayout) findViewById(R.id.chart);
	    		mChartView = execute(getApplicationContext());
	    		layout.removeAllViews();
	    		layout.addView(mChartView, new LayoutParams
	    				(LayoutParams.FILL_PARENT, LayoutParams.FILL_PARENT));
                return true;
            }
            default:
            	return super.onKeyDown(keyCode, event);
            }          
        }
    	return true;
    }
}