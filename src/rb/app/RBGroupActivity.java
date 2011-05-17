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

import android.app.AlertDialog;
import android.app.Dialog;
import android.app.TabActivity;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.res.ColorStateList;
import android.graphics.Color;
import android.os.Bundle;
import android.text.SpannableString;
import android.text.method.LinkMovementMethod;
import android.text.util.Linkify;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.widget.TabHost;
import android.widget.TextView;

public class RBGroupActivity extends TabActivity {

	static final int ABOUT_DIALOG_ID = 0;
	static final int HELP_DIALOG_ID = 1;
	
	public void onCreate(Bundle savedInstanceState)
	{
		super.onCreate(savedInstanceState);
		
		final TabHost tabHost = getTabHost();
		
		tabHost.setup();
		
		//set up formatting for text labels
		TextView label1 = new TextView(this);
		label1.setGravity(android.view.Gravity.CENTER); label1.setTextSize(16); label1.setLines(4);
        label1.setTextColor(new ColorStateList(
                          new int[][] {
                                  new int[] { android.R.attr.state_selected },
                                  new int[0]}
                              , new int[] {
                                  Color.WHITE,
                                  Color.GRAY,
                              }
                          ));
        label1.setSelected(true);

		TextView label2 = new TextView(this);
		label2.setGravity(android.view.Gravity.CENTER); label2.setTextSize(16);
		label2.setTextColor(new ColorStateList(
                new int[][] {
                        new int[] { android.R.attr.state_selected },
                        new int[0]}
                    , new int[] {
                        Color.WHITE,
                        Color.GRAY,
                    }
                ));
		label2.setSelected(false);
		
		
		//add tabs
		TabHost.TabSpec one = tabHost.newTabSpec("tab 1");
		label1.setText("About This Problem");
		one.setIndicator(label1);
		Intent intentOne = new Intent(RBGroupActivity.this, SimpleBrowser.class);
		one.setContent(intentOne);
		tabHost.addTab(one);
		
		TabHost.TabSpec two = tabHost.newTabSpec("tab 2");
		label2.setText("Reduced Basis Solve");
		two.setIndicator(label2);
		Intent intentTwo = new Intent(RBGroupActivity.this, RBActivity.class);
		two.setContent(intentTwo);
		tabHost.addTab(two);
		
		//set which tab has initial focus
		tabHost.setCurrentTab(0);

	}
	
	
	public void onBackPressed(){
		setResult(0);
		getLocalActivityManager().removeAllActivities();
		finish();
	}
	
	
	/**
	 * Re-populates the main menu every time it is brought up depending on
	 * whether help should be displayed
	 */
	public boolean onPrepareOptionsMenu(Menu menu) {

		menu.clear();
		MenuInflater inflater = getMenuInflater();
		if(getTabHost().getCurrentTab() == 0) inflater.inflate(R.menu.info_main_menu, menu);
		else inflater.inflate(R.menu.main_menu, menu); //don't want help text for info screen
		return true;
	}

	/** Handles item selection in the main menu */
	public boolean onOptionsItemSelected(MenuItem item) {
		switch (item.getItemId()) {
		case R.id.about:
			showDialog(ABOUT_DIALOG_ID);
			return true;
		case R.id.help:
			showDialog(HELP_DIALOG_ID);
			return true;
		case R.id.new_problem:
			setResult(0);
			getLocalActivityManager().removeAllActivities();
			finish();
			return true;
		case R.id.quit:
			setResult(-1);
			getLocalActivityManager().removeAllActivities();
			finish();
			return true;
		}
		return false;
	}
	
	protected Dialog onCreateDialog(int id)
	{
		Dialog dialog;
		switch(id){
		case ABOUT_DIALOG_ID:
			AlertDialog.Builder about_builder = new AlertDialog.Builder(
					RBGroupActivity.this);
			about_builder.setTitle("rbAPPmit version 0.1")
			       .setMessage("rbAPPmit:\nAn Android front-end for the Reduced Basis Method\n\n"
			    		    + "Copyright (C) 2010\nDavid Knezevic\nPhuong Huynh\n\n"
			    		    + "Implementation by:\n" +
			    		      "David Knezevic\n" +
			    		      "Phuong Huynh\n" +
			    		      "Mark Wittels\n\n"
			    		    + "This is free software released under the GPLv3 license")
			       .setPositiveButton("OK", new DialogInterface.OnClickListener() {
			           public void onClick(DialogInterface dialog, int id) {
			        	   dismissDialog(ABOUT_DIALOG_ID);
			           }
			       });
			dialog = about_builder.create();

			break;			
		case HELP_DIALOG_ID:			
			final TextView message = new TextView(RBGroupActivity.this);

			final SpannableString s = 
				new SpannableString("Set RB dimension slider to specify the number of basis functions in the approximation.\n\n" +
						            "Set the parameters at which the solve is performed via the parameter sliders or by pressing the adjacent parameter buttons.\n\n" +
						            "Hit solve to perform the reduced basis solve and visualize output(s) and field variable(s).\n\n" +
				                    "For more information see http://augustine.mit.edu/methodology.htm");
			Linkify.addLinks(s, Linkify.WEB_URLS);
			message.setText(s);
			message.setTextSize(18); // 18 sp
			message.setPadding(5,5,5,5); // 5 pixels of padding on all sides
			message.setMovementMethod(LinkMovementMethod.getInstance());

			AlertDialog.Builder help_builder = new AlertDialog.Builder(
					RBGroupActivity.this);
			help_builder.setTitle("rbAPPmit help")
			            .setView(message)
			            .setPositiveButton("OK", new DialogInterface.OnClickListener() {
			           public void onClick(DialogInterface dialog, int id) {
			        	   dismissDialog(HELP_DIALOG_ID);
			           }
			       });

			dialog = help_builder.create();

			break;
			default:
				dialog = null;
		}
		return dialog;
	}
	
}

