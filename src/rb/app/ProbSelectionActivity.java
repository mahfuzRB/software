//rbAPPmit: An Android front-end for the Certified Reduced Basis Method
//Copyright (C) 2010 David J. Knezevic and Phuong Huynh
//
//This file is part of rbAPPmit
//
//rbAPPmit is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//rbAPPmit is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with rbAPPmit.  If not, see <http://www.gnu.org/licenses/>. 

package rb.app;


import android.app.Activity;
import android.app.Dialog;
import android.content.Context;
import android.content.Intent;
import android.os.Bundle;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.AdapterView;
import android.widget.BaseAdapter;
import android.widget.Button;
import android.widget.EditText;
import android.widget.GridView;
import android.widget.ImageView;
import android.widget.LinearLayout;
import android.widget.TextView;


//This activity generates a simple ListView
//to allow selection of demo problems or
//to download data
public class ProbSelectionActivity extends Activity {

// Activity ID
public static final int DEMO_1 = 0;
public static final int DEMO_2 = 1;
public static final int DEMO_3 = 2;
public static final int DEMO_4 = 3;
public static final int DEMO_5 = 4;
public static final int DEMO_6 = 5;
public static final int DEMO_7 = 6;
public static final int DEMO_8 = 7;
public static final int DEMO_9 = 8;
public static final int DEMO_10 = 9;
public static final int DEMO_11 = 10;
public static final int DEMO_12 = 11;
public static final int DEMO_13 = 12;             // Check (Modified for Demo 13)
public static final int DOWNLOAD_FROM_SERVER = 13;
//public static final int DOWNLOAD_FROM_SERVER = 12;
//public static final int DEMO_13 = 13;             // Check (Modifications)


private static int selected_problem = 0;
private static String directory_name;

private GridView gridView;
private Dialog downloadDialog;


public void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);
	
	setContentView(R.layout.prob_selection);
	gridView = (GridView)findViewById(R.id.gridview);
	
	gridView.setStretchMode(GridView.STRETCH_COLUMN_WIDTH);
	gridView.setAdapter(new CaptionImageAdapter(this));
	
	gridView.setOnItemClickListener(new GridView.OnItemClickListener() {

		public void onItemClick(AdapterView<?> av, View view,
				int position, long id) {
			selected_problem = position;
			if(position == DOWNLOAD_FROM_SERVER){
				showDialog(1); 
			}
			else{
				directory_name = "file:///android_asset/demo" + String.valueOf(position+1);
				Intent intent = new Intent(ProbSelectionActivity.this, RBGroupActivity.class);
				ProbSelectionActivity.this.startActivityForResult(intent, 0);
			}
		}
	});
	
}


public static int getProblemSelected(){
	return selected_problem;
}

public static String getDirectory(){
	return directory_name;
}

@Override
protected void onActivityResult(int requestCode, int resultCode, Intent data)
{
	if(resultCode<0) finish();
}

protected Dialog onCreateDialog(int id)
{
	Dialog dialog;
	switch(id){
	case 1:
		downloadDialog = new Dialog(this);
		downloadDialog.setContentView(R.layout.download_dialog);
		downloadDialog.setTitle("Connect to server");
		downloadDialog.setCancelable(false);

		// When the download button is pressed, we read the offline_data
		// files from the specified URL
		Button downloadButton = (Button) downloadDialog
				.findViewById(R.id.downloadButton);
		downloadButton.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {

				// Get the directory_name from the EditText object
				EditText urlEntry = (EditText) downloadDialog
						.findViewById(R.id.urlEntry);
				directory_name = urlEntry.getText().toString();

				// Dismiss the URL specification dialog
				dismissDialog(1);
				Intent intent = new Intent(ProbSelectionActivity.this, RBGroupActivity.class);
				ProbSelectionActivity.this.startActivityForResult(intent, 0);
			}
		});
		
		// Add listener to cancel download
		Button quitDownloadButton = (Button) downloadDialog
				.findViewById(R.id.quitDownloadButton);
		quitDownloadButton.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {
				dismissDialog(1);
			}

		});

		dialog = downloadDialog;
		break;
		default:
			dialog = null;
			break;
	}
	return dialog;
}


//populates the gridview
public class CaptionImageAdapter extends BaseAdapter
{
	private Context mContext;
	
	public CaptionImageAdapter(Context c)
	{
		mContext = c;
	}
	
	@Override
	public int getCount() {
		return mThumbIds.length;
	}

	@Override
	public Object getItem(int arg0) {
		return null;
	}

	@Override
	public long getItemId(int position) {
		return 0;
	}


	@Override
	public View getView(int position, View convertView, ViewGroup parent) {
		View v;
		LayoutInflater li = getLayoutInflater();
		
		v = li.inflate(R.layout.grid_item, null);
		
		TextView tv = (TextView)v.findViewById(R.id.label);
		ImageView iv = (ImageView)v.findViewById(R.id.icon);
		
		tv.setText(mThumbCaptions[position]);
		iv.setImageResource(mThumbIds[position]);
		//necessary for uniform size of icons
		iv.setLayoutParams(new LinearLayout.LayoutParams
				(GridView.LayoutParams.FILL_PARENT, gridView.getHeight()/6));
		Log.d("adapter", "getview called pos: "+position);
		
		
		v.setLayoutParams(new GridView.LayoutParams((int)(gridView.getWidth()/2.5), (int)(gridView.getHeight()/3.5)));
		v.setPadding(2, 4, 2, 8);
		return v;
	}
	
	private final Integer[] mThumbIds = {
			R.drawable.sample_1, 
			R.drawable.sample_2,
			R.drawable.sample_3, 
			R.drawable.sample_4,
			R.drawable.sample_5, 
			R.drawable.sample_6,
			R.drawable.sample_7, 
			R.drawable.sample_8,
			R.drawable.sample_9, 
			R.drawable.sample_10,
			R.drawable.sample_11,
			R.drawable.sample_12,
    		R.drawable.sample_13,        //Check(Modified for Demo 13)
			R.drawable.server_icon 
	};		

	private final String[] mThumbCaptions = {
			"Crack in a pillar",
			"Composite unit cell",
			"MIT logo",
			"Horn",
			"Acoustic impedance",
			"Convection-diffusion",
			"Two-layer conduction",
			"Swiss cheese",
			"FRP",
			"MICA",
			"Navier-Stokes channel flow",
			"Pac-man",
		//	"Uni-Stuttgart",                           //Check(Modified for Demo 13)
			"TinyRB",
			"Download from server"
		};		
}
         
}


