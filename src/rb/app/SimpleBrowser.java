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

import android.app.Activity;
import android.os.Bundle;
import android.util.Log;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.Toast;

public class SimpleBrowser extends Activity
{
	private WebView browser;
	private String url;
	
	public void onCreate(Bundle icicle)
	{
		super.onCreate(icicle);
		
		browser = new WebView(this);
		setContentView(browser);
		
		//load correct url based on demo
		//url = "file:///android_asset/demo" + String.valueOf(ProbSelectionActivity.getProblemSelected()+1)
		url = ProbSelectionActivity.getDirectory() + "/site_info.html";

		
		//if info page was not supplied/another error occurs while loading, redirect user to home site
		browser.setWebViewClient(new WebViewClient(){
			public void onReceivedError(WebView view, int errorCode, String description, String failingUrl) {
				
				// Turn off the toast if we are downloading from server
				if(ProbSelectionActivity.getProblemSelected() != ProbSelectionActivity.DOWNLOAD_FROM_SERVER) {
					Toast.makeText(SimpleBrowser.this, "Sorry, no informational page was supplied with this problem",
						Toast.LENGTH_LONG).show();
					Toast.makeText(SimpleBrowser.this, " You are being redirected to the main RBAppmit site.",
			    		 Toast.LENGTH_LONG).show();
				}
				
			     url = "http://augustine.mit.edu/methodology.htm";
			     if(!failingUrl.equals("http://augustine.mit.edu/methodology.htm"))browser.loadUrl(url);
			   }
		});
		
		//removes white scrollbar from view while keeping scrolling funcionality
		browser.setScrollBarStyle(WebView.SCROLLBARS_OUTSIDE_OVERLAY);
		browser.loadUrl(url);
		
		while(browser.getUrl()==null){
			
		}
		Log.d("browser", browser.getUrl());
	}
	
	public void onBackPressed(){
		//need to tell parent activity to close all activities
		getParent().onBackPressed();
	}
	
}
