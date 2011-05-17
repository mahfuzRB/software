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

import java.math.BigDecimal;
import java.text.DecimalFormat;

import org.achartengine.chart.LineChart;
import org.achartengine.model.XYMultipleSeriesDataset;
import org.achartengine.model.XYSeries;
import org.achartengine.renderer.XYMultipleSeriesRenderer;

import android.util.Log;

/**
 * SingleLabelChart is a subclass of LineChart identical to it, except it can only show
 * data label at a time.
 */
public class SingleLabelChart extends LineChart {

	//which series to show labels for
	private int seriesNum;
	//indicates which data point in the series should be labeled
	private int pointInSeries;
	//whether to show a progress label
	private boolean progressLabel;

	
	public SingleLabelChart(XYMultipleSeriesDataset dataset, XYMultipleSeriesRenderer renderer,
			int s, int p, boolean pr)
	{
		super(dataset, renderer);
		seriesNum = s;
		pointInSeries = p;
		progressLabel = pr;
	}
	

	@Override
	public void drawChartValuesText(android.graphics.Canvas canvas,
            XYSeries series,
            android.graphics.Paint paint,
            float[] points,
            int seriesIndex) //seriesIndex is 0 (mod 3) for output, 1 for LB, 2 for UB
	{
		 
		//showing numbered labels case
		if(!progressLabel){
			if((pointInSeries >= 0)&&(pointInSeries < points.length)){
		
				double accuracyPlace = (series.getMaxY()-series.getMinY())/100;
				//need to get multiple of ten
				double decimalPlace = 0.000001;
				while(accuracyPlace >= decimalPlace*10) decimalPlace *=10;
				accuracyPlace = decimalPlace;
				
				float verticalAdjustment;
				switch(seriesIndex%3){
				default:verticalAdjustment=0; break;
				case 1:verticalAdjustment=8f; break;
				case 2:verticalAdjustment=-8f; break;
				}

			
				//this is essentially identical to LineChart's drawChartValuesText, but it rounds the data
				//points depending on the range of values present, and only plots for one series and one
				//point in that series
				for(int k = 0; k < points.length; k += 2){
					if(k==pointInSeries && seriesIndex/3 == seriesNum)
						super.drawText(canvas, roundToNthPlace(series.getY(k/2), accuracyPlace),
								points[k], (points[k+1]-3.5f+verticalAdjustment), paint, 0);
				}
			}
		}
		//showing progress label case
		else{
			for(int k=0; k< points.length; k+=2){
				if(k>3 && seriesIndex == seriesNum)
					super.drawText(canvas, "working", points[k]+1f, points[k+1]-1f, paint, 0);
			}
		}
	}
	
	//rounds toRound to an arbitrary place, using the multiple of ten accPlace
	//(e.g. 0.01 to round to the hundredths place)
	public static String roundToNthPlace(double toRound, double accPlace)
	{	
		//find leftmost place in double
		double count = 1;
		while(count<=Math.abs(toRound))count*=10;
		int placesBeforeDecimal = (int)Math.log10(count);
		
		int placesAfterDecimal = -(int)Math.log10(accPlace);
		
		String forDecimalFormat = "";
		//if accPlace is greater than or equal to 1, round to whole number
		if(placesAfterDecimal<=0) placesAfterDecimal = 0;

		if(toRound<0)forDecimalFormat += "-";
		for(int i=0; i<placesBeforeDecimal; i++) forDecimalFormat += "0";
		if(placesAfterDecimal>0) forDecimalFormat += ".";
		for(int i=0; i<placesAfterDecimal; i++) forDecimalFormat +="0";
		
		DecimalFormat decimal_format = new DecimalFormat(forDecimalFormat);
		return decimal_format.format(toRound);
				
	}

}

