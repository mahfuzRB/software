package rb.app;

import java.io.FileOutputStream;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.text.DecimalFormat;
import android.app.Activity;
import android.app.AlertDialog;
import android.app.Dialog;
import android.app.ProgressDialog;
import android.content.DialogInterface;
import android.content.Intent;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.text.Html;
import android.text.InputType;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.LinearLayout;
import android.widget.SeekBar;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.TableRow.LayoutParams;
import dalvik.system.DexClassLoader;

// This is the main Activity class for the app.
// This Activity handles downloading the stored
// data, initializing the systems and performing
// the RB solve.

public class RBActivity extends Activity {

// Dialog IDs
static final int DOWNLOAD_DIALOG_ID = 0;
static final int PROGRESS_DIALOG_ID = 1;
static final int DOWNLOAD_FAILED_DIALOG_ID = 2;
static final int RB_SOLVE_DIALOG_ID = 3;
static final int LOAD_DEMO_DIALOG_ID = 4;             //check
static final int SWEEP_DIALOG_ID = 5;
static final int PARAM_DIALOG_ID = 6;


// Activity ID
static final int SELECT_PROBLEM_TYPE = 0;

// String for log printing
private static final String DEBUG_TAG = "RBActivity";

// Member variables that indicate the RB system and SCM types
public RBActivityEnums.SystemTypeEnum mSystemType;
public RBActivityEnums.SCMTypeEnum mSCMType;

/**
* Descriptive member variables for the problem title, variable names, and
* general info.
*/
public String problemTitle;
public String[] paramLabels;
public String problemDescription;
public String descriptionURL;

/**
* The RBSystem object
*/
public static RBSystem mRbSystem;

/**
* The main RBSCMSystem object
*/
public RBSCMSystem mRbScmSystem;

/**
* The second RBSCMSystem object, needed in some time-dependent problems
*/
public RBSCMSystem mSecondRbScmSystem;

public static GLObject mRbModel;

/**
* String to indicate the (URL of the) directory for the offline data.                          //check
*/
private String directory_name;

/**
* Thread object to handle data downloading.
*/
private DownloadThread downloadThread;

/**
* Custom dialog to allow specification of URL to download Offline data
* from.
*/
private Dialog downloadDialog;

/**
* ProgressDialog to display while downloading data.
*/
private ProgressDialog pd;

/**
* Array of TextViews and SeekBars for constructing the parameter selection
*/
private TextView[] mParamLabels;
private SeekBar[] mParamBars;
private Button[] mParamButtons;

/**
* Member variable to store index number of currently selected parameter
* button
*/
private int paramButtonIndex;
private TextView paramInputField;

/**
* The online N constructed by the GUI. We use this value when we call
* RB_solve.
*/
public static int mOnlineNForGui;

/**
* The current parameter constructed by the GUI. We set the RBSystem's
* current parameter to this before performing a solve.
*/
public static Parameter mCurrentParamForGUI;
public static Parameter[] mSweepParam;

/**
* The name of the input file that stores the parameters for the current
* problem.
*/
private String parameters_filename = "input.in";  //check

/**
* The name of the jar file (containing compiled files in .dex form) that we
* download from the server.
*/
private String jarFileName = "AffineFunctions.jar";     //Check

/**
* The corresponding dex file we create locally.
*/
private String dexFileName = "AffineFunctions.dex";   //Check

/**
* The index for the parameter sweep, -1 implies no sweep.
*/
public static int mSweepIndex;


/** Called when the activity is first created. */
@Override
public void onCreate(Bundle savedInstanceState) {
	super.onCreate(savedInstanceState);
	setContentView(R.layout.main);

	// Initialize the RB system and SCM types to NONE
	mSystemType = RBActivityEnums.SystemTypeEnum.NONE;
	mSCMType = RBActivityEnums.SCMTypeEnum.NONE;

	// Set mRbScmSystem and mRbSystem to null initially
	mRbSystem = null;
	mRbScmSystem = null;
	mRbModel = null;

	mRbModel = new GLObject(RBActivity.this);
	mRbModel.allocateBuffer();

	// Set the secondary SCM system to null also
	mSecondRbScmSystem = null;

	// initialize sweep index to -1
	mSweepIndex = -1;
	
	// Add listener to the Solve button
	Button solveButton = (Button) findViewById(R.id.solveButton);
	solveButton.setOnClickListener(new View.OnClickListener() {	
		public void onClick(View view) {
			pd = ProgressDialog.show(RBActivity.this, "", "Solving...");
			SolveThread st = new SolveThread();
			st.start();
		}

	});

	// Attach a listener to onlineNSeekBar
	SeekBar onlineNSeekBar = (SeekBar) findViewById(R.id.onlineNSeekbar);
	onlineNSeekBar
			.setOnSeekBarChangeListener(new SeekBar.OnSeekBarChangeListener() {

				public void onProgressChanged(SeekBar seekBar,
						int progress, boolean fromUser) {
					mOnlineNForGui = (progress + 1);
					TextView currentOnlineNView = (TextView) findViewById(R.id.currentOnlineN);
					currentOnlineNView.setText("Online N =  "
							+ mOnlineNForGui);
				}

				public void onStartTrackingTouch(SeekBar seekBar) {
				}

				public void onStopTrackingTouch(SeekBar seekBar) {
				}
			});

	// Now, call the ListActivity to select problem
	//Intent intent = new Intent(RBActivity.this, ProbSelectionActivity.class);
	//RBActivity.this.startActivityForResult(intent, SELECT_PROBLEM_TYPE);
	
	switch(ProbSelectionActivity.getProblemSelected()) {   //  Check                     
	case ProbSelectionActivity.DOWNLOAD_FROM_SERVER:
		directory_name = ProbSelectionActivity.getDirectory();
		// Open a progress dialog
		pd = ProgressDialog.show(RBActivity.this, "",
				"Reading Offline data files from server...", true);

		// Start the download
		downloadThread = new DownloadThread(downloadHandler,
				directory_name, true);
		downloadThread.start();
		break;
	case ProbSelectionActivity.DEMO_1:
	case ProbSelectionActivity.DEMO_2:
	case ProbSelectionActivity.DEMO_3:
	case ProbSelectionActivity.DEMO_4:
	case ProbSelectionActivity.DEMO_5:
	case ProbSelectionActivity.DEMO_6:
	case ProbSelectionActivity.DEMO_7:
	case ProbSelectionActivity.DEMO_8:
	case ProbSelectionActivity.DEMO_9:
	case ProbSelectionActivity.DEMO_10:
	case ProbSelectionActivity.DEMO_11:
	case ProbSelectionActivity.DEMO_12:  
 	case ProbSelectionActivity.DEMO_13:  //check (Modified for Demo 13)
		
		// Open a progress dialog while we're loading in the data
		//int demo_number = resultCode - ProbSelectionActivity.DEMO_1 + 1;
		int demo_number = ProbSelectionActivity.getProblemSelected() + 1; // Check(ProbselectionActivity.java )
		// it starts from 0 so we need to add 1. 

		Log.d("DEBUG_TAG", "demo_number = " + demo_number);

		pd = ProgressDialog.show(RBActivity.this, "", "Loading Demo "
				+ demo_number + "...", true);
		downloadThread = new DownloadThread(downloadHandler,
				directory_name, false);
		downloadThread.mDemoIndex = demo_number;
		downloadThread.start();
		break;
	default:
		setResult(0);
		finish();
	
	}
	
}

public void onBackPressed(){
	//need to tell parent activity to close all activities
	getParent().onBackPressed();
}

/** Called when the activity is destroyed. */
@Override
protected void onDestroy() {
	super.onDestroy();
	Log.d(DEBUG_TAG, "onDestroy() called");
	// Clean up the files that were downloaded
	delete_downloaded_files();
}

/**
* This function takes care of constructing the dialogs that pop up.
*/
protected Dialog onCreateDialog(int id) {

	Dialog dialog;

	switch (id) {

	
	case DOWNLOAD_FAILED_DIALOG_ID:
		downloadDialog.setTitle("Invalid URL, try again");
		dialog = downloadDialog;
		break;

	case RB_SOLVE_DIALOG_ID:

		dialog = new Dialog(this);
		dialog.setContentView(R.layout.rb_result_dialog);
		dialog.setTitle("RB Solve Results");
		dialog.setCancelable(false);

		Button okButton = (Button) dialog.findViewById(R.id.okButton);
		okButton.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {
				dismissDialog(RB_SOLVE_DIALOG_ID);
				removeDialog(RB_SOLVE_DIALOG_ID);
			}

		});

 		Button visButton = (Button) dialog
				.findViewById(R.id.steadyVisButton);
		visButton.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {

				// mRbModel.nodal_transform(mRbSystem.get_tranformation_data());
				if (mRbSystem.get_mfield() > 0) {
					// Next create the bundle and initialize it
					Bundle bundle = new Bundle();
					/*
					 * bundle.putFloatArray("node",
					 * mRbModel.get_node_data());
					 * bundle.putShortArray("face",
					 * mRbModel.get_face_data()); bundle.putInt("nField",
					 * mRbSystem.get_mfield()); bundle.putBoolean("isReal",
					 * mRbSystem.isReal);
					 */
					/*
					 * if (mRbSystem.isReal) for (int i = 0;
					 * i<mRbSystem.get_mfield(); i++)
					 * bundle.putFloatArray("field"+String.valueOf(i),
					 * mRbSystem.get_truth_sol(i)); else for (int i = 0;
					 * i<mRbSystem.get_mfield(); i++){ float[][] truth_sol =
					 * mRbSystem.get_complex_truth_sol(i);
					 * bundle.putFloatArray("field"+String.valueOf(i)+"R",
					 * truth_sol[0]);
					 * bundle.putFloatArray("field"+String.valueOf(i)+"I",
					 * truth_sol[1]); }
					 */
					Intent intent = new Intent(RBActivity.this,
							RBVisualization.class);
					intent.putExtras(bundle);
					RBActivity.this.startActivity(intent);
				}
			}
		});

		// Create the output string
		String rb_solve_message = "Online N = " + mOnlineNForGui
				+ "\n\u00B5 = [ " + mRbSystem.getCurrentParameters()
				+ "]\n\n";

		DecimalFormat twoPlaces = new DecimalFormat("0.###E0");

		// Create a string that shows each output and error bound
		if (mRbSystem.isReal)
			for (int i = 0; i < mRbSystem.get_n_outputs(); i++) {

				double output_i = mRbSystem.RB_outputs[i];
				double output_bound_i = mRbSystem.RB_output_error_bounds[i];

				rb_solve_message += "Output " + (i + 1) + ":\n"
						+ "Value = " + twoPlaces.format(output_i) + "\n"
						+ "Error bound = "
						+ twoPlaces.format(output_bound_i) + "\n\n";
			}
		else
			for (int i = 0; i < mRbSystem.get_n_outputs(); i++) {

				double output_i_r = mRbSystem.get_RB_output(i, true);
				double output_bound_i_r = mRbSystem
						.get_RB_output_error_bound(i, true);
				double output_i_i = mRbSystem.get_RB_output(i, false);
				double output_bound_i_i = mRbSystem
						.get_RB_output_error_bound(i, false);

				rb_solve_message += "Output " + (i + 1) + ":\n"
						+ "Value = " + twoPlaces.format(output_i_r) + " + "
						+ twoPlaces.format(output_i_i) + "i\n"
						+ "Error bound = "
						+ twoPlaces.format(output_bound_i_r) + " + "
						+ twoPlaces.format(output_bound_i_i) + "i\n\n";
			}

		TextView outputView = (TextView) dialog
				.findViewById(R.id.rb_solve_output);
		outputView.setText(rb_solve_message);

		break;

	case SWEEP_DIALOG_ID:

		try {
			final String[] paramStrings = new String[mRbSystem
					.get_n_params() + 1];

			paramStrings[0] = "No Sweep";
			for (int i = 0; i < paramStrings.length; i++) {
				if (i > 0) {
					paramStrings[i] = "Parameter " + i;
				}
			}

			AlertDialog.Builder builder = new AlertDialog.Builder(
					RBActivity.this);
			builder.setTitle("Pick sweep parameter");
			builder.setItems(paramStrings,
					new DialogInterface.OnClickListener() {
						public void onClick(DialogInterface dialog, int item) {
							// Show a Toast for the selected item
							Toast.makeText(getApplicationContext(),
									paramStrings[item], Toast.LENGTH_SHORT)
									.show();
							// Send a message that indicates which parameter
							// was chosen
							sweepHandler.sendEmptyMessage(item - 1);

							// disable selected slider, enable all others
							// set disabled slider's progress to 0, all
							// others to old values
							try {
								for (int i = 0; i < mRbSystem
										.get_n_params(); i++) {
									mParamBars[i].setEnabled(true);
									double slopeVal = (100 / (mRbSystem
											.getParameterMax(i) - mRbSystem
											.getParameterMin(i)));
									Double progressVal = Double
											.valueOf((slopeVal * mCurrentParamForGUI
													.getEntry(i))
													- (mRbSystem
															.getParameterMin(i) * slopeVal));
									mParamBars[i].setProgress(progressVal
											.intValue());
								}
							} catch (Exception e) {
							}
							if (item >= 1) {
								mParamBars[item - 1].setProgress(0);
								mParamBars[item - 1].setEnabled(false);
							}

							// disable selected parameter button, enable all
							// others
							// set disabled button to "sweep", all others to
							// old values
							try {
								for (int i = 0; i < mRbSystem
										.get_n_params(); i++) {
									displayParamValue(i,
											mCurrentParamForGUI.getEntry(i));
									mParamButtons[i].setEnabled(true);
								}
							} catch (Exception e) {
							}
							if (item >= 1) {
								mParamButtons[item - 1].setText("Sweep");
								mParamButtons[item - 1].setEnabled(false);
							}
						}
					});
			dialog = builder.create();
		} catch (Exception e) {
			Log.e(DEBUG_TAG,
					"Exception thrown during creation of Sweep dialog");
			dialog = null;
		}

		break;

		// Check
		
	case PARAM_DIALOG_ID:
		dialog = new Dialog(this);
		dialog.setContentView(R.layout.param_dialog);
		dialog.setTitle("Minimum: "
				+ mRbSystem.getParameterMin(paramButtonIndex)
				+ " Maximum: "
				+ mRbSystem.getParameterMax(paramButtonIndex));
		dialog.setCancelable(false);

		paramInputField = (EditText) dialog
				.findViewById(R.id.param_input_textview);

		// field should accept signed doubles only
		paramInputField.setInputType(InputType.TYPE_CLASS_NUMBER
				| InputType.TYPE_NUMBER_FLAG_DECIMAL
				| InputType.TYPE_NUMBER_FLAG_SIGNED);

		// user-submitted parameter value will be handled when the ok button
		// is pressed
		Button okButton2 = (Button) dialog
				.findViewById(R.id.param_okButton);
		okButton2.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {
				// determine if value in input field is within acceptable
				// range
				String userParamString = paramInputField.getText()
						.toString();
				double userParam;
				try {
					userParam = Double.parseDouble(userParamString);
				} catch (NumberFormatException e) {
					// if user submits non-double, default value is out of
					// bounds to trigger toast
					userParam = mRbSystem.getParameterMin(paramButtonIndex) - 1;
				}

				if (userParam <= mRbSystem
						.getParameterMax(paramButtonIndex)
						&& userParam >= mRbSystem
								.getParameterMin(paramButtonIndex)) {
					// update parameter bars
					double slopeVal = (100 / (mRbSystem
							.getParameterMax(paramButtonIndex) - mRbSystem
							.getParameterMin(paramButtonIndex)));
					Double progressVal = Double
							.valueOf((slopeVal * userParam)
									- (mRbSystem
											.getParameterMin(paramButtonIndex) * slopeVal));
					mParamBars[paramButtonIndex].setProgress(progressVal
							.intValue());

					// call displayParamValue to update parameter value
					displayParamValue(paramButtonIndex, userParam);
				} else {
					Toast.makeText(getApplicationContext(),
							"Invalid Value", Toast.LENGTH_SHORT).show();
				}

				dismissDialog(PARAM_DIALOG_ID);
				removeDialog(PARAM_DIALOG_ID);
			}

		});

		break;
	default:
		dialog = null;
	}
	return dialog;
}

/**
* A Helper function to display the value of the currently selected
* parameter in the TextView.
* 
* @param current_param
*            The parameter value to display
*/
private void displayParamValue(int index, double current_param) {
	String current_param_str;
	double abs = Math.abs(current_param);
	if ((abs < 0.1) && (current_param != 0.)) {
		DecimalFormat decimal_format = new DecimalFormat("0.###E0");
		current_param_str = decimal_format.format(current_param);
	} 
	else if((abs<1)&&(abs>=0)){
		DecimalFormat decimal_format = new DecimalFormat("@@@");
		current_param_str = decimal_format.format(current_param);
	}
	else{
		DecimalFormat decimal_format = new DecimalFormat("@@@@");
		current_param_str = decimal_format.format(current_param);
	}

	// Make sure we set the parameter to be the same as what the TextView
	// shows
	mCurrentParamForGUI.setEntry(index, Double
			.parseDouble(current_param_str));

	mParamLabels[index].setText(Html.fromHtml(paramLabels[index]));
	mParamButtons[index].setText(Html.fromHtml(current_param_str));

	// Set title
	TextView problemTitleView = (TextView) findViewById(R.id.problemTitle);
	problemTitleView.setText(problemTitle);
}

/**
* Initialize the SeekBar that allows us to select Online N
*/
private void initializeOnlineNBar() {
	// Set max/min of online N seekbar
	SeekBar onlineNSeekBar = (SeekBar) findViewById(R.id.onlineNSeekbar);
	onlineNSeekBar.setMax(mRbSystem.get_n_basis_functions() - 1);

	// Change the progress state so that online N gets initialized
	// to 1
	// If we don't change away from 0, then onProgressChanged
	// doesn't get called
	onlineNSeekBar.setProgress(1);
	onlineNSeekBar.setProgress(0);
}

/**
* Initialize the list view depending on the number of parameters in the
* system and on the parameter ranges.
*/
private void initializeParamBars() {

	// Create String array of parameters to store in the ListView
	try {
		TableLayout paramLayout = (TableLayout) findViewById(R.id.paramLayout);

		// Clear the paramLayout in case we're doing a new problem
		paramLayout.removeAllViews();

		mParamLabels = new TextView[mRbSystem.get_n_params()];
		mParamBars = new SeekBar[mRbSystem.get_n_params()];
		mParamButtons = new Button[mRbSystem.get_n_params()];

		for (int i = 0; i < mRbSystem.get_n_params(); i++) {
			TableRow row = new TableRow(this);
			row.setLayoutParams(new LayoutParams(LayoutParams.FILL_PARENT,
					LayoutParams.FILL_PARENT));

			// First add the text label
			mParamLabels[i] = new TextView(this);
			mParamLabels[i].setTextSize(15); // Size is in scaled pixels
			mParamLabels[i].setLayoutParams(new TableRow.LayoutParams(
					TableRow.LayoutParams.WRAP_CONTENT, TableRow.LayoutParams.WRAP_CONTENT));
			mParamLabels[i].setPadding(0, 0, 4, 0);
			row.addView(mParamLabels[i]);

			// Next add the SeekBar
			mParamBars[i] = new IndexedSeekBar(this);
			((IndexedSeekBar) mParamBars[i]).setIndex(i);
			mParamBars[i].setLayoutParams(new LayoutParams(
					LayoutParams.FILL_PARENT, LayoutParams.WRAP_CONTENT));
			mParamBars[i].setPadding(10, 10, 10, 0); // Set 10px padding on
			// left and right
			row.addView(mParamBars[i]);

			// Finally add the parameter value text
			mParamButtons[i] = new IndexedButton(this);
			((IndexedButton) mParamButtons[i]).setIndex(i);
			row.addView(mParamButtons[i]);

			paramLayout.addView(row, new TableLayout.LayoutParams(
					LayoutParams.FILL_PARENT, LayoutParams.WRAP_CONTENT));
		}

		// Initialize mCurrentParamForGUI to min_parameter
		mCurrentParamForGUI = new Parameter(mRbSystem.get_n_params());
		for (int i = 0; i < mRbSystem.get_n_params(); i++) {
			double min_param = mRbSystem.getParameterMin(i);

			mCurrentParamForGUI.setEntry(i, min_param);
			displayParamValue(i, min_param);

			// Also set param bars to zero to match min_param
			mParamBars[i].setProgress(0);
		}
	} catch (InconsistentStateException e) {
		Log.e(DEBUG_TAG, e.getMessage());
	} catch (Exception e) {
		e.printStackTrace();
	}

	try {
		addParamBarListeners();
	} catch (InconsistentStateException e) {
		Log
				.e(DEBUG_TAG,
						"Exception occurred when adding listeners to the parameter SeekBars");
	}

	try {
		addParamButtonListeners();
	} catch (InconsistentStateException e) {
	}

}

/**
* Add a new button to perform a parameter sweep
*/
private void initializeParamSweep() {

	try {

		LinearLayout sweepLayout = (LinearLayout) findViewById(R.id.sweepButtonHolder);

		Button sweepButton = new Button(this);
		sweepButton.setText("\u00B5 Sweep");
		sweepButton.setTextSize(22);
		sweepButton.setOnClickListener(new View.OnClickListener() {

			public void onClick(View view) {

				// Create an alert dialog with radio buttons for
				// selecting the sweep parameter
				showDialog(SWEEP_DIALOG_ID);
			}
		});

		sweepLayout.addView(sweepButton, new LinearLayout.LayoutParams(
				LayoutParams.FILL_PARENT, LayoutParams.WRAP_CONTENT));

	} catch (Exception e) {
		e.printStackTrace();
	}

}


// This helper function adds listeners to
// the parameter value buttons
private void addParamButtonListeners() throws InconsistentStateException {

	for (int i = 0; i < mRbSystem.get_n_params(); i++) {
		mParamButtons[i].setOnClickListener(new View.OnClickListener() {
			public void onClick(View v) {
				paramButtonIndex = ((IndexedButton) v).getIndex();
				showDialog(PARAM_DIALOG_ID);
			}
		});
	}
}

// This helper function adds the listeners
// to the newly built parameter SeekBars
private void addParamBarListeners() throws InconsistentStateException {

	// Add a listener to each SeekBar
	for (int i = 0; i < mRbSystem.get_n_params(); i++) {
		mParamBars[i]
				.setOnSeekBarChangeListener(new SeekBar.OnSeekBarChangeListener() {

					public void onProgressChanged(SeekBar seekBar,
							int progress, boolean fromUser) {

						if (fromUser) {
							IndexedSeekBar isb = (IndexedSeekBar) seekBar;
							int index = isb.getIndex();

							if (mRbSystem != null) {
								double param_range = mRbSystem
										.getParameterMax(index)
										- mRbSystem.getParameterMin(index);

								double current_param = mRbSystem
										.getParameterMin(index)
										+ param_range
										* progress
										/ seekBar.getMax();

								displayParamValue(index, current_param);
							}

						}
					}

					public void onStartTrackingTouch(SeekBar seekBar) {
					}

					public void onStopTrackingTouch(SeekBar seekBar) {
					}
				});
	}

}

// Private helper function that reads in which type of systems we
// need for this problem and sets the member flags appropriately
private void read_system_types_from_input_file(String parameters_filename,
		boolean isAssetFile) {
	GetPot infile = new GetPot(this, parameters_filename, isAssetFile);

	problemTitle = infile.call("title", "RB Online");
	int parameter_number = infile.call("n_parameters", 0);
	paramLabels = new String[parameter_number];
	for (int n = 0; n < parameter_number; n++) {
		paramLabels[n] = infile.call("param" + Integer.toString(n)
				+ "_label", "DEFAULT");

		// If DEFAULT was read in, then replace with a default
		// mu label
		if (paramLabels[n] == "DEFAULT") {
			paramLabels[n] = "\u00B5" + (n + 1);
		}
	}

	descriptionURL = infile.call("descriptionURL", "");

	String SystemTypeEnum_in = infile.call("system_type", "NONE");
	mSystemType = getSystemEnumFromString(SystemTypeEnum_in);

	String SCMTypeEnum_in = infile.call("scm_type", "NONE");
	mSCMType = getSCMEnumFromString(SCMTypeEnum_in);

	Log.d(DEBUG_TAG, "RB system type = " + mSystemType);
	Log.d(DEBUG_TAG, "SCM type = " + mSCMType);
}

// Private helper function that reads in which type of systems we
// need for this problem and sets the member flags appropriately
private void attach_affine_functions(InputStream in) throws Exception {

	// Create a local copy of AffineFunctions.jar
	FileOutputStream f = openFileOutput("AffineFunctions.jar",
			MODE_WORLD_READABLE);

	byte[] buffer = new byte[1024];
	int len1 = 0;
	while ((len1 = in.read(buffer)) > 0) {
		f.write(buffer, 0, len1);
	}
	in.close();

	Log.d(DEBUG_TAG, "Finished copying jar file");

	DexClassLoader cl = new DexClassLoader("/data/data/rb.app/files/"
			+ jarFileName, "/data/data/rb.app/files/", null, ClassLoader
			.getSystemClassLoader());

	Log.d(DEBUG_TAG, "Created local dex file");

	if (mRbSystem != null) {
		mRbSystem.mAffineFnsClass = cl.loadClass("AffineFunctions");            // Check
		mRbSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();             // Check

		// Set Q_a, Q_f and n_outputs from the loaded class
		mRbSystem.read_in_Q_a();
		Log.d(DEBUG_TAG, "Q_a = " + mRbSystem.get_Q_a());

		mRbSystem.read_in_Q_f();
		Log.d(DEBUG_TAG, "Q_f = " + mRbSystem.get_Q_f());

		mRbSystem.read_in_n_outputs();
		Log.d(DEBUG_TAG, "n_outputs = " + mRbSystem.get_n_outputs());

		mRbSystem.read_in_Q_uL();

		if (mSystemType == RBActivityEnums.SystemTypeEnum.LINEAR_UNSTEADY
				|| mSystemType == RBActivityEnums.SystemTypeEnum.QN_UNSTEADY) {
			TransientRBSystem trans_rb = (TransientRBSystem) mRbSystem;
			trans_rb.read_in_Q_m();
			Log.d(DEBUG_TAG, "Q_m = " + trans_rb.get_Q_m());
		}
	}

	if (mRbScmSystem != null) {
		mRbScmSystem.mAffineFnsClass = cl.loadClass("AffineFunctions");
		mRbScmSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

		// set Q_a
		mRbScmSystem.read_in_Q_a();
	}
	if (mSecondRbScmSystem != null) {
		mSecondRbScmSystem.mAffineFnsClass = cl
				.loadClass("AffineFunctions");
		mSecondRbScmSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

		// set Q_a
		mSecondRbScmSystem.read_in_Q_a();
	}

}

// Private helper function that initializes the RB and SCM systems
private void initialize_systems(String path, boolean asset_file_bool)
		throws Exception {

	// First, clear all the systems in case we're doing a new (different)
	// problem
	mRbScmSystem = null;
	mSecondRbScmSystem = null;
	mRbSystem = null;

	// Find out which type of systems we'll be initializing
	RBActivity.this
			.read_system_types_from_input_file(path, asset_file_bool);

	// Initialize the SCM systems
	if (mSCMType == RBActivityEnums.SCMTypeEnum.COERCIVE_ALPHASIGMA) {
		mRbScmSystem = RBSCMSystem.buildSCMSystem(RBActivity.this,
				RBActivityEnums.SCMTypeEnum.COERCIVE);
		mSecondRbScmSystem = RBSCMSystem.buildSCMSystem(RBActivity.this,
				RBActivityEnums.SCMTypeEnum.COERCIVE);
	} else if (mSCMType == RBActivityEnums.SCMTypeEnum.COERCIVE
			|| mSCMType == RBActivityEnums.SCMTypeEnum.QN_TRANSIENT_SCM) {
		mRbScmSystem = RBSCMSystem
				.buildSCMSystem(RBActivity.this, mSCMType);
	} else if (mSCMType == RBActivityEnums.SCMTypeEnum.COMPLEX_NONCOERCIVE) {
		mRbScmSystem = RBSCMSystem
				.buildSCMSystem(RBActivity.this, mSCMType);
	} else {
		mRbScmSystem = null;
	}

	// Read parameters into SCM systems
	if (mRbScmSystem != null) {
		mRbScmSystem.parse_parameters_file(path, asset_file_bool);
	}
	if (mSecondRbScmSystem != null) {
		mSecondRbScmSystem.parse_parameters_file(path, asset_file_bool);
	}

	mRbSystem = RBSystem.buildRBSystem(RBActivity.this, mSystemType);
	mRbSystem.setPrimarySCM(mRbScmSystem);
	if (mSystemType == RBActivityEnums.SystemTypeEnum.LINEAR_UNSTEADY) {
		TransientRBSystem trans_rb = (TransientRBSystem) mRbSystem;
		trans_rb.setSecondarySCM(mSecondRbScmSystem);
	}

	// Read parameters into RB systems
	if (mRbSystem != null) {
		mRbSystem.parse_parameters_file(path, asset_file_bool);
		/*
		 * if(mRbSystem.get_mfield() > 0) mRbModel = new
		 * GLObject(RBActivity.this);
		 */
	}
}

// Private helper function to delete the downloaded files (if they exist)
private void delete_downloaded_files() {
	deleteFile(parameters_filename);
	deleteFile(jarFileName);
	deleteFile(dexFileName);
}

// Private helper function to return a SystemTypeEnum based on a String
private static RBActivityEnums.SystemTypeEnum getSystemEnumFromString(
		String s) {

	for (RBActivityEnums.SystemTypeEnum type : RBActivityEnums.SystemTypeEnum
			.values()) {
		if (type.toString().equals(s)) {
			return type;
		}
	}

	return RBActivityEnums.SystemTypeEnum.NONE;
}

// Private helper function to return a SystemTypeEnum based on a String
private static RBActivityEnums.SCMTypeEnum getSCMEnumFromString(String s) {

	for (RBActivityEnums.SCMTypeEnum type : RBActivityEnums.SCMTypeEnum
			.values()) {
		if (type.toString().equals(s))
			return type;
	}

	return RBActivityEnums.SCMTypeEnum.NONE;
}

final Handler sweepHandler = new Handler() {
	public void handleMessage(Message msg) {
		mSweepIndex = msg.what;
		Log.d(DEBUG_TAG, "Sweep index set to " + mSweepIndex);
	}
};

// Define the Handler that receives messages from the thread and updates the
// progress
// handler is a final member object of the RBActivity class
final Handler downloadHandler = new Handler() {
	public void handleMessage(Message msg) {
		pd.dismiss();

		// Now check if there was a problem or not
		boolean downloadSuccessful = msg.getData().getBoolean(
				"Download successful");
		Log.d(DEBUG_TAG, "downloadSuccessful = " + downloadSuccessful);

		if (!downloadSuccessful) {
			showDialog(DOWNLOAD_FAILED_DIALOG_ID);
			delete_downloaded_files();
		} else {

			// Initialize the SeekBar for Online N
			RBActivity.this.initializeOnlineNBar();

			// Initialize the ListView for the parameters
			RBActivity.this.initializeParamBars();

//			// Set link to problem info page
//			TextView linkView = (TextView) findViewById(R.id.link_view);
//			linkView.setText
//			("http://augustine.mit.edu/methodology.htm"); // re-write this to be
//											              // problem-specific

			if ((mSystemType == RBActivityEnums.SystemTypeEnum.LINEAR_STEADY)
					|| (mSystemType == RBActivityEnums.SystemTypeEnum.LINEAR_COMPLEX_STEADY)) {
				RBActivity.this.initializeParamSweep();
			}
		}

	}
};


/** Nested class that performs progress calculations (counting) */
private class DownloadThread extends Thread {

	// URL to download from
	String mUrl;

	// Handler to interact with this thread
	Handler mHandler;

	// Boolean to indicate if we download from server
	// or load from the assets directory
	boolean isDownload;

	// The index of the demo we have chosen
	int mDemoIndex;

	/**
	 * Constructor takes a handler and a String specifying the URL of the
	 * Offline data directory
	 */
	DownloadThread(Handler h, String directory_name, boolean isDownload_in) {
		mHandler = h;
		mUrl = directory_name;
		isDownload = isDownload_in;
	}

	/**
	 * This function gets called when we execute Thread.start()
	 */
	public void run() {

		boolean download_successful;

		if (isDownload) {
			download_successful = download_from_server();
		} else {
			download_successful = load_from_asset(mDemoIndex);
		}

		Message msg = mHandler.obtainMessage();
		Bundle b = new Bundle();
		b.putBoolean("Download successful", download_successful);
		msg.setData(b);
		mHandler.sendMessage(msg);
	}

	public boolean download_from_server() {

		boolean download_successful = true;

		final boolean asset_file_bool = false; // Do not load from Android
		// asset

		// Copy parameters filename to a local file
		try {
			// Access parameters_filename from the server
			URL u = new URL(directory_name + "/" + parameters_filename);
			HttpURLConnection c = (HttpURLConnection) u.openConnection();
			c.setRequestMethod("GET");
			c.setDoOutput(true);
			c.connect();

			// Create a local copy of parameters_filename so that we
			// can use GetPot
			FileOutputStream f = openFileOutput(parameters_filename,
					MODE_WORLD_READABLE);

			InputStream in = c.getInputStream();

			byte[] buffer = new byte[1024];
			int len1 = 0;
			while ((len1 = in.read(buffer)) > 0) {
				f.write(buffer, 0, len1);
			}
			in.close();

			Log.d(DEBUG_TAG, "Successfully downloaded "
					+ parameters_filename);

			initialize_systems(parameters_filename, asset_file_bool);

		} catch (InconsistentStateException e) {
			Log.e(DEBUG_TAG,
					"Inconsistent state exception occurred when parsing input file: "
							+ e.getMessage());
			download_successful = false;
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception thrown when accessing "
					+ directory_name + "/" + parameters_filename);
			download_successful = false;
		}

		try {
			URL u = new URL(directory_name + "/" + jarFileName);
			HttpURLConnection c = (HttpURLConnection) u.openConnection();
			c.setRequestMethod("GET");
			c.setDoOutput(true);
			c.connect();

			Log.d(DEBUG_TAG, "Connected to url of jar file");

			InputStream in = c.getInputStream();

			attach_affine_functions(in);
		} catch (Exception e) {
			Log.e(DEBUG_TAG,
					"Exception occurred while loading affine functions");
			download_successful = false;
		}

		try {
			if (mRbSystem != null) {
				mRbSystem.read_offline_data(mUrl, asset_file_bool);
			}

			if (mRbScmSystem != null) {
				mRbScmSystem.read_offline_data(mUrl, asset_file_bool);
			}

			if (mRbModel != null) {
				mRbModel.read_offline_data(mUrl, asset_file_bool);
			}

		} catch (Exception e) {
			Log
					.e(DEBUG_TAG,
							"Exception occurred while reading saved data from server");
			download_successful = false;
		}

		return download_successful;
	}

	// This function initializes the Online solver based on the demo
	// files stored in res/raw/
	private boolean load_from_asset(int demoIndex) {

		String directory_name = "demo" + demoIndex;
		final boolean asset_file_bool = true; // Load from Android asset

		try {
			initialize_systems(directory_name + "/" + parameters_filename,
					asset_file_bool);
		} catch (InconsistentStateException e) {
			Log.e(DEBUG_TAG,
					"Inconsistent state exception occurred when parsing input file: "
							+ e.getMessage());
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception thrown when accessing "
					+ directory_name + "/" + parameters_filename);
		}

		// Now copy AffineFunctions.jar to data/data/files so that we can
		// unpack it
		try {
			InputStream in = getAssets().open(
					directory_name + "/" + jarFileName);

			attach_affine_functions(in);

		} catch (Exception e) {
			Log.e(DEBUG_TAG,
					"Exception occurred while loading affine functions");
		}

		// Finally, initialize the RB and SCM systems
		try {
			if (mRbSystem != null) {
				mRbSystem
						.read_offline_data(directory_name, asset_file_bool);
			}

			if (mRbScmSystem != null) {
				mRbScmSystem.read_offline_data(directory_name,
						asset_file_bool);
			}

			if (mRbModel != null) {
				mRbModel.read_offline_data(directory_name, asset_file_bool);
			}

		} catch (Exception e) {
			Log
					.e(DEBUG_TAG,
							"Exception occurred while reading saved data from assets");
		}

		return true;
	}

}

private class SolveThread extends Thread
{
	public void run()
	{
		switch (mSystemType) {

		case LINEAR_STEADY:
		case LINEAR_COMPLEX_STEADY:

			if (mSweepIndex == -1) {

				mRbSystem.setCurrentParameters(mCurrentParamForGUI);
				mRbSystem.RB_solve(mOnlineNForGui);

				handler.sendEmptyMessage(0);
			} else { // We need to perform a sweep
				int numSweepPts = 10;
				numSweepPts = Math.round(100000 / (mRbSystem
						.get_mfield() * mRbSystem.get_calN()));
				if (!mRbSystem.isReal)
					numSweepPts /= 3;
				numSweepPts = numSweepPts > 10 ? 10 : numSweepPts;
				//numSweepPts = 50;
				int n_outputs = mRbSystem.get_n_outputs();

				mSweepParam = new Parameter[numSweepPts];
				
				double[][][] RB_sweep_sol = null;
				if (mRbSystem.isReal) {
					RB_sweep_sol = new double[numSweepPts][1][mRbSystem
							.get_N()];
				} else {
					RB_sweep_sol = new double[numSweepPts][2][mRbSystem
							.get_N()];
					n_outputs *= 2;
				}

				double[][] sweepOutputs = new double[n_outputs][numSweepPts];
				double[][] sweepOutputBounds = new double[n_outputs][numSweepPts];

				// Create the bundle and initialize it
				Bundle bundle = new Bundle();

				double sweepParamRange = mRbSystem
						.getParameterMax(mSweepIndex)
						- mRbSystem.getParameterMin(mSweepIndex);
				double sweepIncrement = sweepParamRange
						/ (numSweepPts - 1);

				float[][][] vLTfunc = new float[numSweepPts][][];

				for (int i = 0; i < numSweepPts; i++) {
					double new_param = mRbSystem
							.getParameterMin(mSweepIndex)
							+ i * sweepIncrement;
					mCurrentParamForGUI
							.setEntry(mSweepIndex, new_param);
					mRbSystem.setCurrentParameters(mCurrentParamForGUI);
					mSweepParam[i] = mCurrentParamForGUI.clone();
					Log.d(DEBUG_TAG, "Set new param "
							+ mCurrentParamForGUI);
					mRbSystem.RB_solve(mOnlineNForGui);

					if (mRbSystem.isReal)
						for (int n = 0; n < n_outputs; n++) {
							sweepOutputs[n][i] = mRbSystem.RB_outputs[n];
							sweepOutputBounds[n][i] = mRbSystem.RB_output_error_bounds[n];
						}
					else
						for (int n = 0; n < n_outputs / 2; n++) {
							sweepOutputs[n][i] = mRbSystem
									.get_RB_output(n, true);
							sweepOutputs[n + n_outputs / 2][i] = mRbSystem
									.get_RB_output(n, false);
							sweepOutputBounds[n][i] = mRbSystem
									.get_RB_output_error_bound(n, true);
							sweepOutputBounds[n + n_outputs / 2][i] = mRbSystem
									.get_RB_output_error_bound(n, false);
						}

					if (mRbSystem.get_mfield() > 0) {
						RB_sweep_sol[i] = mRbSystem.get_RBsolution();
						vLTfunc[i] = mRbSystem.get_tranformation_data();
					}
				}
				mRbModel.vLTfunc = vLTfunc;
				mRbSystem.set_sweep_sol(RB_sweep_sol);

				bundle.putBoolean("isSweep", true);
				bundle.putString("title", "Online N = "
						+ mOnlineNForGui);
				bundle.putDouble("dt", sweepIncrement);
				bundle.putDouble("xMin", mRbSystem
						.getParameterMin(mSweepIndex));
				bundle.putDouble("xMax", mRbSystem
						.getParameterMax(mSweepIndex));
				bundle.putString("xLabel", Integer.toString(mSweepIndex + 1));
				bundle.putInt("n_time_steps", numSweepPts);
				bundle.putInt("n_outputs", n_outputs);
				for (int i = 0; i < n_outputs; i++) {
					bundle.putDoubleArray("output_data_" + i,
							sweepOutputs[i]);
					bundle.putDoubleArray("output_bound_" + i,
							sweepOutputBounds[i]);
				}

				// Add this bundle to the intent and plot
				Intent intent = new Intent(RBActivity.this,
						OutputPlotterActivity.class);
				intent.putExtras(bundle);
				RBActivity.this.startActivity(intent);
			}

			break;
		case LINEAR_UNSTEADY:
		case QN_UNSTEADY:

			// Perform the solve
			mRbSystem.setCurrentParameters(mCurrentParamForGUI);
			mRbSystem.RB_solve(mOnlineNForGui);

			// Next create the bundle and initialize it
			Bundle bundle = new Bundle();
			bundle.putBoolean("isReal", mRbSystem.isReal);
			bundle.putBoolean("isSweep", false);
			bundle
					.putString("title", "Online N = " + mOnlineNForGui
							+ ", parameter = "
							+ mCurrentParamForGUI.toString());
			bundle.putDouble("dt", mRbSystem.get_dt());
			bundle.putDouble("xMin", 0);
			bundle.putDouble("xMax", mRbSystem.get_dt()
					* mRbSystem.get_K());
			bundle.putString("xLabel", "time");
			bundle.putInt("n_time_steps", mRbSystem.n_plotting_steps); //mRbSystem.get_K() + 1
			bundle.putInt("n_outputs", mRbSystem.get_n_outputs());
			for (int i = 0; i < mRbSystem.get_n_outputs(); i++) {
				bundle.putDoubleArray("output_data_" + i,
						mRbSystem.RB_outputs_all_k[i]);
				bundle.putDoubleArray("output_bound_" + i,
						mRbSystem.RB_output_error_bounds_all_k[i]);
			}

			// Add this bundle to the intent
			Intent intent = new Intent(RBActivity.this,
					OutputPlotterActivity.class);
			intent.putExtras(bundle);
			RBActivity.this.startActivity(intent);

			break;
		default:
			throw new RuntimeException(
					"Invalid RB system type for solve");
		}
		
		handler.sendEmptyMessage(-1);
	}
	
	private Handler handler = new Handler()
	{
		public void handleMessage(Message msg){
			pd.dismiss();
			if(msg.what==0) showDialog(RB_SOLVE_DIALOG_ID);
		}
	};
	
}

}
