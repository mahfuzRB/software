package rb.app; // Added line to specify the rb.app package

//  -*- java -*- 
//  GetPot Version 1.0                                        Sept/13/2002
//  
//  WEBSITE: http://getpot.sourceforge.net
//  
//  This library is  free software; you can redistribute  it and/or modify
//  it  under  the terms  of  the GNU  Lesser  General  Public License  as
//  published by the  Free Software Foundation; either version  2.1 of the
//  License, or (at your option) any later version.
//  
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//  
//  You  should have  received a  copy of  the GNU  Lesser  General Public
//  License along  with this library; if  not, write to  the Free Software
//  Foundation, Inc.,  59 Temple Place,  Suite 330, Boston,  MA 02111-1307
//  USA
//  
//  (C) 2001 Frank R. Schaefer  
//==========================================================================
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PushbackInputStream;
import java.util.Enumeration;
import java.util.StringTokenizer;
import java.util.Vector;

import android.content.Context;
import android.util.Log;

public class GetPot {
	
	// Added Android Context for file I/O
	public Context androidContext; 
	// Flag that indicates whether the file we're reading is in Android's asset directory
	public boolean isAssetFile;
	
    // variable structure
    public class variable {
	// (*) constructor
	variable(String Name, String Value) {
	    name = Name;
	    take(Value);
	}
	variable() {
	    name = "";
	    value = new String[1];
	    value[0] = "";
	    original = "";
	}
	// (*) get a specific element in the string vector
	String  get_element(int Idx) { 
	    if( Idx < 0 || Idx > value.length ) return new String("");
	    else                                return value[Idx];
	}

	void take(String Value) {
	    original = Value;
	    StringTokenizer  st = new StringTokenizer(Value);
	    value = new String[st.countTokens()];
	    for(int i =0; st.hasMoreTokens(); i++ )
		value[i] = st.nextToken();
	}
	// (*) data memebers
	String  name;      // identifier of variable
	String  value[];   // value of variable stored in vector
	String  original;  // value of variable as given on command line
    };
    String    prefix;          // prefix indicating section for search
    String    section;         // section for dollar bracket parsing
    Vector    argv;            // String[]
    Vector    section_list;
    int       cursor; 
    Vector    idx_nominus;     // int[]
    int       nominus_cursor;
    boolean   search_loop_f;
    boolean   search_failed_f;
    Vector    variables;       // class variable


    // (*) constructor    
    GetPot(String argv_[], String applicationName) {
	prefix = "";
	cursor = 0;
	nominus_cursor = -1;
	search_loop_f = true;
	search_failed_f = false;

	Vector __argv      = new Vector();
	idx_nominus = new Vector();
	variables   = new Vector();
	section_list = new Vector();
	
	// in java the argv[0] does not contain the application name
	// => append application name to argv vector by hand

	for(int i=0; i < argv_.length ; i++) {
	    String copied_str = new String(argv_[i]);
	    __argv.addElement(copied_str);
	}

	argv = new Vector();
	__argv.insertElementAt(applicationName, 0);
	__parse_argument_vector(__argv);
    }
     
    GetPot(Context c_in, String FileName, boolean isAssetFile_in) {
    	
    androidContext = c_in; // Added Android Context for file I/O
    isAssetFile    = isAssetFile_in; // Added boolean for app
    	
	prefix = "";
	cursor = 0;
	nominus_cursor = -1;
	search_loop_f = true;
	search_failed_f = false;

	Vector __argv      = new Vector();
	idx_nominus = new Vector();
	variables   = new Vector();
	section_list = new Vector();

	__argv = read_in_file(FileName);

	argv = new Vector();
	__argv.insertElementAt(FileName, 0);
	__parse_argument_vector(__argv);
    }

    void __parse_argument_vector(Vector _argv) {
	// build internal database
	//   1) array with no-minus arguments (usually used as filenames)
	//   2) section label treatment
	//   3) variable assignments:
	//             'variable name' '=' number | string
	int i=0;
	section = "";
	Vector section_stack = new Vector();
	if( _argv.size() < 1 ) return;

	for(Enumeration arge = _argv.elements(); arge.hasMoreElements() ;i++) {
	    String arg = (String)arge.nextElement();

	    if( i == 0 ) { argv.addElement(arg); continue; }
	    if( arg.length() == 0) continue;

	    // -- [section] labels
	    if( arg.length() > 1 && arg.charAt(0) == '[' && arg.charAt(arg.length()-1) == ']' ) {
		String name = DBE_expand_string(arg.substring(1, arg.length()-1));
		// String name = arg.substring(1, arg.length()-1);
		section = process_section_label(name, section_stack);
		// new section --> append to list of sections
		if( section_list.indexOf(section) == -1 )
		    if( section.length() != 0 ) section_list.addElement(section);
		argv.addElement(arg);
	    }
	    else {
		arg = section + DBE_expand_string(arg);
		argv.addElement(arg);
	    }

	    // 1) separate array for nominus arguments
	    if( arg.charAt(0) != '-' ) 
		idx_nominus.addElement(new Integer(i));
	    
	    // 2) variables: does _argv[i] contain a '=' operator ?
	    int k = arg.indexOf('=');
	    if( k > 0 ) { 
		variable Var = find_variable(arg.substring(0,k));
		if( Var.name.compareTo("") == 0 ) 
		    variables.addElement(new variable(arg.substring(0,k), 
						      arg.substring(k+1)));
		else                              
		    Var.take(arg.substring(k+1));
	    }
	}
    }
    // file parsing
    Vector  read_in_file(String filename) {
	Vector Empty = new Vector();
	try { 
	    // InputStream  i = new FileInputStream(filename);
		InputStream i;
		if(!isAssetFile) {  // Modified for Android
			i = androidContext.openFileInput(filename);
		}
		else {
			i = androidContext.getAssets().open(filename);
		}
	    if( i.available() == 0 ) return Empty;
	    // argv[0] == the filename of the file that was read in
	    return read_in_stream(new PushbackInputStream(i));
	} 
	catch(FileNotFoundException ID_10T) { return Empty; } 
	catch(IOException ID_20T) { return Empty; }
    }

    Vector read_in_stream(PushbackInputStream istr) {
	Vector  brute_tokens = new Vector();
	String  token;
	Vector  arglist = new Vector();

	while( 1 + 1 == 2 ) {
	    if( __skip_whitespace(istr) == -1 ) break;
	    token = __get_next_token(istr);
	    if( token.length() == 0 ) break;
	    brute_tokens.addElement(token);
	}

	// for(int i=0; i<brute_tokens.size() ; i++)
	// System.out.println((String)(brute_tokens.elementAt(i)));

	// -- reduce expressions of token1'='token2 to a single 
	//    string 'token1=token2'
	// -- copy everything into 'result'
	// -- arguments preceded by something like '[' name ']' (section)
	//    produce a second copy of each argument with a prefix '[name]argument'
	int i1 = 0;
	int i2 = 1;
	int i3 = 2;
	String   section = new String("");

	while( i1 < brute_tokens.size() ) {
	    String SRef = (String)brute_tokens.elementAt(i1);
	    // 1) concatinate 'abcdef' '=' 'efgasdef' to 'abcdef=efgasdef'
	    // note: java.lang.String: substring(a,b) = from a to b-1
	    //        C++ string:      substr(a,b)    = from a to a + b
	    if( i2 < brute_tokens.size() && ((String)brute_tokens.elementAt(i2)).charAt(0) == '=' ) {
		if( i3 >= brute_tokens.size() )
		    arglist.addElement((String)brute_tokens.elementAt(i1) 
				      + (String)brute_tokens.elementAt(i2));
		else
		    arglist.addElement((String)brute_tokens.elementAt(i1) 
				      + (String)brute_tokens.elementAt(i2) 
				      + (String)brute_tokens.elementAt(i3));
		
		i1 = i3+1; i2 = i3+2; i3 = i3+3;
		continue;
	    }
	    else {
		arglist.addElement(SRef);
		i1=i2; i2=i3; i3++;
	    }
	}
	return arglist;
    }


    int __skip_whitespace(PushbackInputStream istr) { 
	try { 
	    while( 1 + 1 == 2 ) {
		int tmp = istr.read();  // read first character out of stream
		if( tmp == -1 ) return -1;
		while( (char)tmp == ' ' || (char)tmp == '\t' || (char)tmp == '\n' )
		{ tmp = istr.read(); if( tmp == -1 ) return -1; }
	    
		// found a non whitespace 
		if( (char)tmp == '#' ) {
		    // comment => skip until end of line
		    do { tmp = istr.read(); if( tmp == -1 ) return -1; }
		    while( (char)tmp != '\n' );

		    istr.unread(tmp);
		    continue; 
		}
		else {
		    // go back to last mark
		    istr.unread(tmp);
		    return 1;
		}
	    }
	}
	catch(IOException ID_10T) {}
	return 1;
    }

    String __get_next_token(PushbackInputStream istr) {
	// get next concatinates string token. consider quotes that embrace
	// whitespaces
	try {
	    String  token = new String("");
	    int     tmp0 = 0;
	    char    tmp = 0; 
	    char    last_letter = 0;
	    while(1+1 == 2) {
		last_letter = tmp; tmp0 = istr.read(); tmp = (char)tmp0;
		if( tmp == -1
		    || ((tmp == ' ' || tmp == '\t' || tmp == '\n') && last_letter != '\\') ) {
		    return token;
		}
		else if( tmp == '\'' && last_letter != '\\' ) {
		    // QUOTES: un-backslashed quotes => it's a string
		    token += __get_string(istr);
		    continue;
		}
		else if( tmp == '{' && last_letter == '$') {
		    token += '{' + __get_until_closing_bracket(istr);
		    continue;
		}
		else if( tmp == '$' && last_letter == '\\') {
		    token += tmp; tmp = 0;  //  so that last_letter will become = 0, not '$';
		    continue;
		}
		else if( tmp == '\\' && last_letter != '\\') 
		    continue;              // don't append un-backslashed backslashes
		token += tmp;
	    }
	}
	catch(IOException ID_10T) {}
	return "";
    }

    String __get_string(PushbackInputStream istr) {
	// parse input until next matching '

	try {
	    String str = new String("");
	    int     tmp0 = 0;
	    char    tmp = 0;
	    char    last_letter = 0;
	    while(1 + 1 == 2) {
		last_letter = tmp; tmp0 = istr.read(); tmp = (char)tmp0;
		if( tmp0 == -1 )  return str;
		// un-backslashed quotes => it's the end of the string
		else if( tmp == '\'' && last_letter != '\\')  return str;
		else if( tmp == '\\' && last_letter != '\\')  continue; // don't append 

		str += tmp;
	    }
	    // return str;
	}
	catch(IOException ID_10T) {}
	return "";
    }

    String  __get_until_closing_bracket(PushbackInputStream istr) {
	try {
	    // parse input until next matching }
	    String str = new String(""); 
	    int    tmp0 = 0;
	    char   tmp = 0;
	    char   last_letter = 0;
	    int    brackets = 1;
	    while(1 + 1 == 2) {
		last_letter = tmp; tmp0 = istr.read(); tmp = (char)tmp0;
		if( tmp == -1 ) return str;
		else if( tmp == '{' && last_letter == '$') brackets += 1;
		else if( tmp == '}') {
		    brackets -= 1;
		    // un-backslashed brackets => it's the end of the string
		    if( brackets == 0) return str + '}';
		    else if( tmp == '\\' && last_letter != '\\') 
			continue;  // do not append an unbackslashed backslash
		}
		str += tmp;
	    }
	}
	catch(IOException ID_10T) {}
	return "";
    }

    String  process_section_label(String label, 
				  Vector section_stack) {
	String sname = label;
	//  1) subsection of actual section ('./' prefix)
	if( sname.length() >= 2 && sname.substring(0, 2).compareTo("./") == 0) {
	    sname = sname.substring(2, sname.length());
	}
	//  2) subsection of parent section ('../' prefix)
	else if( sname.length() >= 3 && sname.substring(0, 3).compareTo("../") == 0) {
	    do {
		if( section_stack.size() != 0 ) 
		    section_stack.removeElementAt(section_stack.size()-1);
		sname = sname.substring(3, sname.length());
	    } while( sname.substring(0, 3).compareTo("../") == 0 );
	}
	// 3) subsection of the root-section
	else {
	    // a root section
	    section_stack.removeAllElements();
	}

	if( sname.compareTo("") != 0 ) {
	    // parse section name for 'slashes'
	    int i=0;
	    while( i < sname.length() ) {
		if( sname.charAt(i) == '/' ) {
		    section_stack.addElement(sname.substring(0,i));
		    if( i+1 < sname.length()-1 )
			sname = sname.substring(i+1, sname.length());
		    i = 0;
		}
		else 
		    i++;
	    }	   
	    section_stack.addElement(sname);
	} 
	// else: "[]" => back to root section

	String section = new String("");
	if( section_stack.size() != 0 )
	    for(int i=0; i< section_stack.size() ; i++)
		section += (String)(section_stack.elementAt(i)) + "/";
	return section;
    }
    // (*) helper functions
    int __convert_to_type(String S, int Default) {
	int start_i = 0;
	if( S.length() >= 2 && S.substring(0,2).compareTo("0x") == 0 )     start_i = 2;
	else if( S.length() >=3 && S.substring(0,3).compareTo("-0x") == 0) start_i = 3;
	else {
	    // a normal integer, not a hex number
	    try { return Integer.parseInt(S); }
	    catch(NumberFormatException ID_10T) { 
		try { return (int)java.lang.Math.rint(Double.valueOf(S).doubleValue()); }
		catch(NumberFormatException ID_20T) { return Default; }
	    }
	}

	// a hexadecimal number
	int number = 0;
	for(int i=start_i; i<S.length(); i++) {
	    int c = (int)S.charAt(i);
	    int digit = 0;
	    if( c >= '0' && c <= '9')      digit = c - '0';
	    else if( c >= 'a' && c <= 'f') digit = c - 'a';
	    else if( c >= 'A' && c <= 'F') digit = c - 'a';
	    else 
		break;
	    number *= 16;
	    number += digit;
	}
	if( start_i == 2 ) return number;
	else               return -number;
    }

    double __convert_to_type(String S, double Default) {
	try { return Double.valueOf(S).doubleValue(); }
	catch(NumberFormatException ID_10T) { return Default; }
    }

    void   set_prefix(String Prefix) { prefix = Prefix; }

    String  __get_remaining_string(String String, String Start) {
	// Checks if 'String' begins with 'Start' and returns the remaining String.
	// Returns None if String does not begin with Start.
	
	// note: java.lang.String: substring(a,b) = from a to b-1
	//        C++ string:      substr(a,b)    = from a to a + b
	if( Start.compareTo("") == 0 )   return String;
	if( String.indexOf(Start) == 0 ) return String.substring(Start.length(), String.length());
	else                             return "";
    }


    ///////////////////////////////////////////////////////////////////////////////
    // (*) cursor oriented functions
    //.............................................................................
    //
    void disable_loop() { search_loop_f = false; }
    void enable_loop()  { search_loop_f = true; }
    void reset_cursor()	{ search_failed_f = false; cursor = 0; }
    
    void init_multiple_occurrence()
	{ disable_loop(); reset_cursor(); }

    //     -- search for a certain argument and set cursor to position
    boolean search(String Option) { 
	int     OldCursor  = cursor;
	String  SearchTerm = prefix + Option;

	if( OldCursor >= argv.size() ) OldCursor = argv.size() - 1;
	search_failed_f = true;

	// (*) first loop from cursor position until end
	for(int c = cursor; c < argv.size(); c++) {
	    if( ((String)argv.elementAt(c)).compareTo(SearchTerm) == 0 ) 
	    { cursor = c; search_failed_f = false; return true; }
	}
	if( ! search_loop_f ) return false;
	
	// (*) second loop from 0 to old cursor position
	for(int c = 1; c < OldCursor; c++) {
	    if( ((String)argv.elementAt(c)).compareTo(SearchTerm) == 0 ) 
	    { cursor = c; search_failed_f = false; return true; }
	}
	// in case nothing is found the cursor stays where it was
	return false;
    }


    boolean search(String P[]) {
	for(int i=0; i<P.length; i++)
	    if( search(P[i]) == true ) return true;
	return false;
    }

    boolean search(String P1, String P2) {
	if( search(P1) == true ) return true;
	else if( search(P2) == true ) return true;
	else return false;
    }

    boolean search(String P1, String P2, String P3) {
	if( search(P1, P2) == true ) return true;
	else if( search(P3) == true ) return true;
	else return false;
    }

    boolean search(String P1, String P2, String P3, String P4) {
	if( search(P1, P2, P3) == true ) return true;
	else if( search(P4) == true ) return true;
	else return false;
    }

    boolean search(String P1, String P2, String P3, String P4, 
		   String P5) {
	if( search(P1, P2, P3, P4) == true ) return true;
	else if( search(P5) == true ) return true;
	else return false;
    }

    boolean search(String P1, String P2, String P3, String P4, 
		   String P5, String P6) {
	if( search(P1, P2, P3, P4, P5) == true ) return true;
	else if( search(P6) == true ) return true;
	else return false;
    }
    ///////////////////////////////////////////////////////////////////////////////
    // (*) No operator overloading in java
    //     => give them the 'python names'
    String getitem(int idx) { 
	if( idx < argv.size() ) return (String)argv.elementAt(idx) ; 
	else                    return new String(""); 
    }

    int get(int Idx, int Default) { 
	if( Idx >= argv.size() ) return Default;
	return __convert_to_type((String)argv.elementAt(Idx), Default);
    }

    double get(int Idx, double Default) { 
	if( Idx >= argv.size() ) return Default;
	return __convert_to_type((String)argv.elementAt(Idx), Default);
    }

    String get(int Idx, String Default) { 
	if( Idx >= argv.size() ) return Default;
	else                     return (String)argv.elementAt(Idx);
    }

    int    size() { return argv.size(); }
    //     -- next() function group
    int next(int Default) { 
	if( search_failed_f ) return Default;
	cursor++; 
	if( cursor >= argv.size() ) 
	{ cursor = argv.size(); return Default; }

	if( prefix.compareTo("") == 0) 
	    return __convert_to_type((String)argv.elementAt(cursor), Default);
	
	String Remain = __get_remaining_string((String)argv.elementAt(cursor), prefix);
	
	return Remain.compareTo("") != 0 ? __convert_to_type(Remain, Default) : Default;
    }

    double next(double Default) { 
	if( search_failed_f ) return Default;
	cursor++;
	if( cursor >= argv.size() )
	{ cursor = argv.size(); return Default; }

	if( prefix.compareTo("") == 0) 
	    return __convert_to_type((String)argv.elementAt(cursor), Default);
	
	String Remain = __get_remaining_string((String)argv.elementAt(cursor), prefix);
	
	return Remain.compareTo("") != 0 ? __convert_to_type(Remain, Default) : Default;
    }

    String next(String Default) { 
	if( search_failed_f ) return Default;
	cursor++; 
	if( cursor >= argv.size() )
	{ cursor = argv.size(); return Default; }

	if( prefix.compareTo("") == 0) 
	    return (String)argv.elementAt(cursor);
	
	String Remain = __get_remaining_string((String)argv.elementAt(cursor), prefix);
	
	return Remain.compareTo("") != 0 ? Remain : Default;
    }
    //     -- follow() function group 
    //        distinct option to be searched for
    int follow(int Default, String A1) { 
	if( search(A1) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A[]) { 
	if( search(A) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A1, String A2) { 
	if( search(A1, A2) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A1, String A2, String A3) { 
	if( search(A1, A2, A3) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A1, String A2, String A3, String A4) { 
	if( search(A1, A2, A3, A4) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A1, String A2, String A3, String A4,
	       String A5) { 
	if( search(A1, A2, A3, A4, A5) == false ) return Default;
	return next(Default);
    }
    int follow(int Default, String A1, String A2, String A3, String A4,
	       String A5, String A6) { 
	if( search(A1, A2, A3, A4, A5, A6) == false ) return Default;
	return next(Default);
    }

    ///////// double
    double follow(double Default, String A1) { 
	if( search(A1) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A[]) { 
	if( search(A) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A1, String A2) { 
	if( search(A1, A2) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A1, String A2, String A3) { 
	if( search(A1, A2, A3) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A1, String A2, String A3, String A4) { 
	if( search(A1, A2, A3, A4) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A1, String A2, String A3, String A4,
	       String A5) { 
	if( search(A1, A2, A3, A4, A5) == false ) return Default;
	return next(Default);
    }
    double follow(double Default, String A1, String A2, String A3, String A4,
	       String A5, String A6) { 
	if( search(A1, A2, A3, A4, A5, A6) == false ) return Default;
	return next(Default);
    }
    ///////// String ...
    String follow(String Default, String A1) { 
	if( search(A1) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A[]) { 
	if( search(A) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A1, String A2) { 
	if( search(A1, A2) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A1, String A2, String A3) { 
	if( search(A1, A2, A3) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A1, String A2, String A3, String A4) { 
	if( search(A1, A2, A3, A4) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A1, String A2, String A3, String A4,
	       String A5) { 
	if( search(A1, A2, A3, A4, A5) == false ) return Default;
	return next(Default);
    }
    String follow(String Default, String A1, String A2, String A3, String A4,
	       String A5, String A6) { 
	if( search(A1, A2, A3, A4, A5, A6) == false ) return Default;
	return next(Default);
    }
    // (*) directly followed arguments
    int direct_follow(int Default, String Option) {
	String FollowStr = match_starting_string(Option);
	if( FollowStr.compareTo("") == 0 ) return Default;

	int Value = __convert_to_type(FollowStr, Default);
	if( ++cursor >= argv.size() ) cursor = argv.size();
	return Value;
    }

    double direct_follow(double Default, String Option) {
	String FollowStr = match_starting_string(Option);
	if( FollowStr.compareTo("") == 0 ) return Default;
	double Value = __convert_to_type(FollowStr, Default);
	if( ++cursor >= argv.size() ) cursor = argv.size();
	return Value;
    }
    
    String direct_follow(String Default, String Option)	{
	String FollowStr = match_starting_string(Option);
	if( FollowStr.compareTo("") == 0) return Default;
	if( ++cursor >= argv.size() ) cursor = argv.size();
	return FollowStr;
    }
    
    String match_starting_string(String StartString) {
	// pointer  to the place where the string after 
	//          the match inside the found argument starts.
	// 0        no argument matches the starting string.
	int N = StartString.length();
	int OldCursor = cursor;
	if( OldCursor >= argv.size() ) OldCursor = argv.size() - 1;
	search_failed_f = true;
	
	// (*) first loop from cursor position until end
	for(int c = cursor; c < argv.size(); c++) {
	    String Arg = (String)argv.elementAt(c);
	    if( Arg.length() >= N && Arg.substring(0,N).compareTo(StartString) == 0) { 
		cursor = c; search_failed_f = false; 
		return ((String)argv.elementAt(c)).substring(N); 
	    }
	}
	
	if( search_loop_f == false ) return "";
	
	// (*) second loop from 0 to old cursor position
	for(int c = 1; c < OldCursor; c++) {
	    String Arg = (String)argv.elementAt(c);
	    if( Arg.length() >= N && Arg.substring(0,N).compareTo(StartString)==0 ) { 
		cursor = c; search_failed_f = false; 
		return ((String)argv.elementAt(c)).substring(N); 
	    }
	}
	return "";
    }
    ///////////////////////////////////////////////////////////////////////////////
    // (*) search for flags
    //.............................................................................
    //
    boolean options_contain(String FlagList) {
	// go through all arguments that start with a '-', but not more than one
	for(int i = 0; i<argv.size() ; i++) {
	    String str;
	    if( prefix.compareTo("") != 0 ) str = __get_remaining_string((String)argv.elementAt(i), prefix);
	    else                            str = (String)argv.elementAt(i);

	    if( str.length() > 1 && str.charAt(0) == '-' && str.charAt(1) != '-' )
		if( check_flags(str, FlagList) ) return true;
	}
	return false;
    }

    boolean argument_contains(int Idx, String FlagList) {
	if( Idx >= argv.size() || Idx < 0) return false;

	if( prefix.compareTo("") == 0 )
	    return check_flags((String)argv.elementAt(Idx), FlagList);
	
	// if a prefix is set, then the argument index is the index
	//   inside the 'namespace'
	// => only check list of arguments that start with prefix
	int no_matches = 0;
	for(int i=0; i < argv.size(); i++) {
	    String Remain = __get_remaining_string((String)argv.elementAt(i), prefix);
	    if( Remain.compareTo("") != 0 ) {
		no_matches += 1;
		if( no_matches == Idx)
		    return check_flags(Remain, FlagList);
	    }
	}
	// no argument in this namespace
	return false;

    }

    boolean check_flags(String Str, String FlagList) {
	for(int i=0; i != FlagList.length() ; i++)
	    if( Str.indexOf(FlagList.charAt(i)) != -1 )
		// one of the flags was found in given string
		return true;
	return false;
    }
    ///////////////////////////////////////////////////////////////////////////////
    // (*) nominus arguments
    String[] nominus_vector() {
	String nv[] = new String[idx_nominus.size()];
	int i=0;
	for(Enumeration inm = idx_nominus.elements() ; inm.hasMoreElements(); i++) {
	    Integer e = (Integer)inm.nextElement();
	    nv[i] = (String)argv.elementAt(e.intValue());
	}
	return nv;
    }

    String next_nominus() {
	if( nominus_cursor >= idx_nominus.size()-1 ) return "";
	int C = ++nominus_cursor;
	return (String)argv.elementAt(((Integer)idx_nominus.elementAt(C)).intValue());
    }

    void  reset_nominus_cursor()
	{ nominus_cursor = -1; }

    boolean  search_failed() 
	{ return search_failed_f; }
    ///////////////////////////////////////////////////////////////////////////////
    // (*) variables
    //     (no operator overloading, => functions with 'python names' for 
    //     operators
    int call(String VarName, int Default) {
	variable sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0 ) return Default;
	return __convert_to_type(sv.original, Default);
    }

    double call(String VarName,  double Default) {
	variable  sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return Default;
	return __convert_to_type(sv.original, Default);
    }
    
    String call(String VarName, String Default) {
	variable  sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return Default;
	return sv.original;
    }
    
    int call(String VarName, int Default, int Idx) {
	variable sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return Default;
	String   element = sv.get_element(Idx);
	if( element.compareTo("") == 0) return Default;
	return __convert_to_type(element, Default);
    }

    double call(String VarName,  double Default, int Idx) {
	variable sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return Default;
	String  element = sv.get_element(Idx);
	if( element.compareTo("") == 0) return Default;
	return __convert_to_type(element, Default);
    }

    String call(String VarName, String Default, int Idx) {
	variable  sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return Default;
	String element = sv.get_element(Idx);     
	if( element.compareTo("") == 0) return Default;
	return element;
    }

    int vector_variable_size(String VarName) {
	variable  sv = find_variable(VarName);
	if( sv.name.compareTo("") == 0) return 0;
	return sv.value.length;
    }

    variable find_variable(String NameX) {
	String Name = prefix + NameX;
	for(Enumeration inm = variables.elements() ; inm.hasMoreElements(); ) {
	    variable v = (variable)inm.nextElement();
	    if( v.name.compareTo(Name) == 0 ) return v;
	}
	return new variable("","");
    }

    String[]  get_variable_names() {
	Vector result = new Vector();
	int    length = 0;
	for(Enumeration inm = variables.elements() ; inm.hasMoreElements(); ) {
	    String Name = ((variable)inm.nextElement()).name;
	    String Tmp = __get_remaining_string(Name, prefix);
	    if( Tmp != "") {
		result.addElement(new String(Tmp));
		length++;
	    }
	}
	String end[] = new String[length];
	int i = 0;
	for(Enumeration inm = result.elements() ; inm.hasMoreElements(); i++) {
	    String Name = (String)inm.nextElement();
	    end[i] = Name;
	}
	return end;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // (*) ouput (basically for debugging reasons)
    //.............................................................................
    //
    void print() {
	System.out.println("argc = " + argv.size() + "\n");	
	for(int i=0; i<argv.size() ; i++)
	    System.out.println(argv.elementAt(i));
	return;
    }

  // (*) dollar bracket expressions (DBEs) ------------------------------------
  //
  //     1) Entry Function: dbe_expand_string()
  //        Takes a string such as
  //
  //          "${+ ${x} ${y}}   Subject-${& ${section} ${subsection}}:   ${title}"
  //
  //        calls DBE_expand_string() for each of the expressions
  //
  //           ${+ ${x} ${y}}
  //           ${& ${section} ${subsection}}
  //           ${Title}
  //
  //        and returns the string
  //
  //          "4711 Subject-1.01:   Mit den Clowns kamen die Schwaene"
  //
  //        assuming that
  //            x          = "4699"
  //            y          = "12"
  //            section    = "1."
  //            subsection = "01"
  //            title      = "Mit den Clowns kamen die Schwaene"
  //
  //      2) DBE_expand():
  //
  //           checks for the command, i.e. the 'sign' that follows '${'
  //           divides the argument list into sub-expressions using
  //           DBE_get_expr_list()
  //
  //           ${+ ${x} ${y}}                 -> "${x}"  "${y}"
  //           ${& ${section} ${subsection}}  -> "${section}" "${subsection}"
  //           ${Title}                       -> Nothing, variable expansion
  //
  //      3) DBE_expression_list():
  //
  //           builds a vector of unbracketed whitespace separated strings, i.e.
  //
  //           "  ${Number}.a ${: Das Marmorbild} AB-${& Author= ${Eichendorf}-1870}"
  //
  //           is split into a vector
  //
  //              [0] ${Number}.a
  //              [1] ${: Das Marmorbild}
  //              [2] AB-${& Author= ${Eichendorf}}-1870
  //
  //           Each sub-expression is expanded using expand(). 
  //---------------------------------------------------------------------------    
    String DBE_expand_string(String str) {
	// Parses for closing operators '${ }' and expands them letting
	// white spaces and other letters as they are.
	String   new_string = "";
	int      open_brackets = 0;
	int      first = 0;
	for(int i = 0;  i<str.length(); i++) {
	    if( i < str.length() - 2 && str.substring(i, i+2).compareTo("${") == 0 ) {
		if( open_brackets == 0 ) first = i+2;
		open_brackets++;
	    }
	    else if( str.charAt(i) == '}' && open_brackets > 0) {
		open_brackets -= 1;
		if( open_brackets == 0 ) {
		    String Replacement = DBE_expand(str.substring(first, i));
		    new_string += Replacement;
		}
	    }
	    else if( open_brackets == 0 )
		new_string += str.charAt(i);
	}
	return new_string;
    }

    Vector DBE_get_expr_list(String str_, int ExpectedNumber) {
	// ensures that the resulting vector has the expected number
	// of arguments, but they may contain an error message

	String str = str_;
	// Separates expressions by non-bracketed whitespaces, expands them
	// and puts them into a list.
	int i=0;

	// (1) eat initial whitespaces
	for(; i < str.length(); i++) {
	    char tmp = str.charAt(i);
	    if( tmp != ' ' && tmp != '\t' && tmp != '\n' ) break;
	}

	Vector expr_list = new Vector();
	int    open_brackets = 0;
	Vector start_idx = new Vector();
	int    start_new_string = i;
	int    l = str.length();

	// (2) search for ${ } expressions ...
	while( i < l ) {
	    char letter = str.charAt(i);
	    // whitespace -> end of expression
	    if( (letter == ' ' || letter == '\t' || letter == '\n') && open_brackets == 0) {
		expr_list.addElement(str.substring(start_new_string, i));
		boolean no_breakout_f = true;
		for(i++; i < l ; i++) {
		    char tmp = str.charAt(i);
		    if( tmp != ' ' && tmp != '\t' && tmp != '\n' ) 
		    { no_breakout_f = false; start_new_string = i; break; }
		}
		if( no_breakout_f == true ) {
		    // end of expression list
		    if( expr_list.size() < ExpectedNumber ) {
			String   pre_tmp = "<< ${ }: missing arguments>>";
			for(int ni = expr_list.size(); ni < ExpectedNumber; ni++)
			    expr_list.addElement(pre_tmp);
		    }
		    return expr_list;
		}
	    }
	    
	    // dollar-bracket expression
	    if( str.length() >= i+2 && str.substring(i, i+2).compareTo("${") == 0 ) {
		open_brackets++;
		start_idx.addElement(new Integer(i+2));
	    }
	    else if( letter == '}' && open_brackets > 0) {
		int start = ((Integer)start_idx.elementAt(start_idx.size()-1)).intValue();
		start_idx.removeElementAt(start_idx.size()-1);
		String Replacement = DBE_expand(str.substring(start, i));
		if( start - 3 < 0)
		    str = Replacement + str.substring(i+1);
		else
		    str = str.substring(0, start-2) + Replacement + str.substring(i+1);
		l = str.length();
		i = start + Replacement.length() - 3;
		open_brackets--;
	    }
	    i++;
	}
	
	// end of expression list
	expr_list.addElement(str.substring(start_new_string, i));

	if( expr_list.size() < ExpectedNumber ) {
	    String   pre_tmp = "<< ${ }: missing arguments>>";
	    for(int ni = expr_list.size(); ni < ExpectedNumber; ni++)
		expr_list.addElement(pre_tmp);
	}
	return expr_list;
    }
    
    variable DBE_get_variable(String VarName) {
	variable  ev = new variable();
	String    secure_Prefix = prefix;
    
	prefix = section;
	// (1) first search in currently active section
	variable var = find_variable(VarName);
	if( var.name != "" ) { prefix = secure_Prefix; return var; }
    
	// (2) search in root name space
	prefix = "";
	var = find_variable(VarName);
	if( var.name != "" ) { prefix = secure_Prefix; return var; }

	prefix = secure_Prefix;
    
	// error occured => variable name == ""
	ev.name = "";
	ev.original = "<<${ } variable '" + VarName + "' undefined>>";
	return ev;
    }
  
    String DBE_expand(String expr) {
	// ${: } pure text
	if( expr.charAt(0) == ':' )
	    return expr.substring(1);
      
	// ${& expr expr ... } text concatination
	else if( expr.charAt(0) == '&' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 1);
	    String result = "";
	    for(Enumeration e = A.elements() ; e.hasMoreElements() ; ) {
		String add = (String)e.nextElement();
		result += add;
	    }       
	    return result;
	}
      
	// ${<-> expr expr expr} text replacement
	else if( expr.length() >= 3 && expr.substring(0, 3).compareTo("<->") == 0) {
	    Vector A = DBE_get_expr_list(expr.substring(3), 3);
	    String Arg0 = (String)A.elementAt(0);
	    String Arg1 = (String)A.elementAt(1);
	    String Arg2 = (String)A.elementAt(2);
	    int    L = Arg1.length();

	    while( 1 + 1 == 2 ) {
		int tmp = Arg0.indexOf(Arg1);
		if( tmp == -1 ) break;
		if( tmp == 0 ) Arg0 = Arg2 + Arg0.substring(L);
		else           Arg0 = Arg0.substring(0, tmp) + Arg2 + Arg0.substring(L+tmp);
	    }
	    return Arg0;
	}
	// ${+ ...}, ${- ...}, ${* ...}, ${/ ...} expressions
	else if( expr.charAt(0) == '+' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 2);

	    double result = 0.0;

	    for(Enumeration e = A.elements() ; e.hasMoreElements() ; ) {
		String add = (String)e.nextElement();
		result += __convert_to_type(add, 0.0);
	    }       
	    
	    return (new Double(result)).toString();
	}
	else if( expr.charAt(0) == '-' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 2);

	    Enumeration e = A.elements(); 
	    double result = __convert_to_type((String)e.nextElement(), 0.0);
	    for(; e.hasMoreElements() ; ) {
		String sub = (String)e.nextElement();
		result -= __convert_to_type(sub, 0.0);
	    }       
	    return (new Double(result)).toString();
	}
	else if( expr.charAt(0) == '*' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 2);

	    double result = 1.0;
	    for(Enumeration e = A.elements() ; e.hasMoreElements() ; ) {
		String mul = (String)e.nextElement();
		result *= __convert_to_type(mul, 0.0);
	    }       

	    return (new Double(result)).toString();
	}
	else if( expr.charAt(0) == '/' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 2);
	    Enumeration e = A.elements();
	    double result = __convert_to_type((String)e.nextElement(), 0.0);
	    if( result == 0 ) return "0.0";
	    for(; e.hasMoreElements() ; ) {
		double Q = __convert_to_type((String)e.nextElement(), 0.0);
		if( Q == 0.0) return "0.0";
		result /= Q;		
	    }
	    return (new Double(result)).toString();
	}

	// ${^ ... } power expressions
	else if( expr.charAt(0) == '^' ) {
	    Vector A = DBE_get_expr_list(expr.substring(1), 2);

	    Enumeration e = A.elements(); 
	    double result = __convert_to_type((String)e.nextElement(), 0.0);
	    for(; e.hasMoreElements() ; ) {
		double p = __convert_to_type((String)e.nextElement(), 0.0);
		result = java.lang.Math.pow(result, p);
	    }
	    
	    return (new Double(result)).toString();
	}
    
	// ${==  } ${<=  } ${>= } comparisons (return the number of the first 'match'
	else if(expr.length() >= 2  && 
		(expr.substring(0,2).compareTo("==") == 0 || 
		 expr.substring(0,2).compareTo(">=") == 0 || 
		 expr.substring(0,2).compareTo("<=") == 0 || 
		 expr.charAt(0) == '>' || 
		 expr.charAt(0) == '<') ) {
	    // differentiate between two and one sign operators
	    int op = 0;
	    // EQ = 0   GEQ = 1   LEQ = 2   GT = 3   LT = 4
	    if      ( expr.substring(0, 2).compareTo("==") == 0 ) op = 0;
	    else if ( expr.substring(0, 2).compareTo(">=") == 0 ) op = 1;
	    else if ( expr.substring(0, 2).compareTo("<=") == 0 ) op = 2;
	    else if ( expr.charAt(0) == '>' )                     op = 3;
	    else                                                  op = 4;
	
	    Vector a = new Vector();
	    if ( op >= 3 ) a = DBE_get_expr_list(expr.substring(1), 2);
	    else           a = DBE_get_expr_list(expr.substring(2), 2);

	    Enumeration e = a.elements();
	    String   x_orig = (String)e.nextElement();
	    double   x = __convert_to_type(x_orig, 1e37);
	    int i = 1;
	    for(; e.hasMoreElements() ; ) {
		String y_orig = (String)e.nextElement();
		double y = __convert_to_type(y_orig, 1e37);

		// set the strings as reference if something wasn't a number
		if ( x == 1e37 || y == 1e37 ) {
		    // it's a string comparison
		    int comparison = x_orig.compareTo(y_orig);

		    if( (   op == 0  && comparison == 0)
			|| (op == 1  && comparison >= 0)
			|| (op == 2  && comparison <= 0)
			|| (op == 3  && comparison > 0)
			|| (op == 4  && comparison < 0) ) {
			return (new Integer(i)).toString();
		    }
		}
		else {
		    // it's a number comparison
		    if( (op == 0  && x == y) || (op == 1 && x >= y) ||
			(op == 2  && x <= y) || (op == 3 && x >  y) || 
			(op == 4  && x <  y) ) {
			return (new Integer(i)).toString();
		    }
		}
		i++;
	    }
	    // nothing fulfills the condition => return 0
	    return "0";
	}
	// ${?? expr expr} select 
	else if( expr.length() >= 2 && expr.substring(0, 2).compareTo("??") == 0) {
	    Vector a = DBE_get_expr_list(expr.substring(2), 2);
	    double x = __convert_to_type((String)a.elementAt(0), 1e37);
	    // last element is always the default argument
	    if( x == 1e37 || x < 0 || x >= a.size() - 1 ) 
		return (String)a.elementAt(a.size()-1);

	    // round x to closest integer
	    int rx = (int)java.lang.Math.rint(x);
	    return (String)a.elementAt(rx);
	}
	// ${? expr expr expr} if then else conditions
	else if( expr.charAt(0) == '?' ) {
	    Vector a = DBE_get_expr_list(expr.substring(1), 2);
	    if( __convert_to_type((String)a.elementAt(0), 0.0) == 1.0 ) 
		return (String)a.elementAt(1);
	    else if( a.size() > 2 ) 
		return (String)a.elementAt(2);
	}
	// ${! expr} maxro expansion 
	else if( expr.charAt(0) == '!' ) {
	    variable Var = DBE_get_variable(expr.substring(1));
	    // error
	    if( Var.name == "" ) return Var.original;

	    Vector A = DBE_get_expr_list(Var.original, 2);
	    return (String)A.elementAt(0);
	}
	// ${@: } - string subscription
	else if( expr.length() >=2 && expr.substring(0,2).compareTo("@:") == 0 ) {
	    Vector A = DBE_get_expr_list(expr.substring(2), 2);
	    double x = __convert_to_type((String)A.elementAt(1), 1e37);
	
	    // last element is always the default argument
	    if( x == 1e37 || x < 0 || x >= (double)((String)A.elementAt(0)).length() - 1)
		return "<<1st index out of range>>";
	
	    if( A.size() > 2 ) {
		double y = __convert_to_type((String)A.elementAt(2), 1e37);
		if ( y != 1e37 && y > 0 && y <= (double)((String)A.elementAt(0)).length() - 1 && y > x ) {
		    // note: java.lang.String: substring(a,b) = from a to b - 1
		    //        C++ string:      substr(a,b)    = from a to a + b		    
		    int rx1 = (int)java.lang.Math.rint(x);
		    int rx2 = (int)java.lang.Math.rint(y + 1.0);
		    return ((String)A.elementAt(0)).substring(rx1, rx2);
		}
		else if( y == -1 ) {		    
		    int rx = (int)java.lang.Math.rint(x);
		    return ((String)A.elementAt(0)).substring(rx);
		}
		return "<<2nd index out of range>>";
	    }
	    else {
		String tmp = (String)A.elementAt(0);
		int rx = (int)java.lang.Math.rint(x);
		return (String)tmp.substring(rx, rx+1);
	    }
	}
	// ${@ } - vector subscription
	else if( expr.charAt(0) == '@' ) {
	    Vector    A   = DBE_get_expr_list(expr.substring(1), 2);
	    variable  Var = DBE_get_variable((String)A.elementAt(0));
	    // error
	    if( Var.name == "" ) {    
		// make a copy of the string if an error occured
		// (since the error variable is a static variable inside get_variable())
		return Var.original;
	    }

	    double x = __convert_to_type((String)A.elementAt(1), 1e37);
      
	    // last element is always the default argument
	    if (x == 1e37 || x < 0 || x >= Var.value.length )
		return "<<1st index out of range>>";
	    
	    if ( A.size() > 2) {
		double y = __convert_to_type((String)A.elementAt(2), 1e37);
		int    begin = (int)java.lang.Math.rint(x);
		int    end = 0;
		if ( y != 1e37 && y > 0 && y <= Var.value.length && y > x)
		    end = (int)java.lang.Math.rint(y + 1.0);
		else if( y == -1 )
		    end = Var.value.length;
		else
		    return "<<2nd index out of range>>";	    	    
		
		String result = Var.get_element(begin);
		for(int i = begin+1; i < end; i++) 
		    result += " " + Var.get_element(i);	    
		return result;
	    }
	    else {
		int tmp = (int)java.lang.Math.rint(x);
		return Var.get_element(tmp);
	    }
	}    
	Vector    A = DBE_get_expr_list(expr, 1);
	variable  B = DBE_get_variable((String)A.elementAt(0));
	
	// make a copy of the string if an error occured
	// (since the error variable is a static variable inside get_variable())
	if( B.name == "" ) return (String)B.original;
	else               return B.original;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // (*) unidentified flying objects
    //.............................................................................
    //
    String[] string_vector_to_string_array(Vector Vec) {
	String[] result = new String[Vec.size()];
	for(int i =0; i < Vec.size(); i++ )
	    result[i] = (String)Vec.elementAt(i);
	return result;
    }

    int __search_string_array(String[] Arr, String Str) {
	for(int i=0; i < Arr.length ; i++)
	    if( Str.compareTo(Arr[i]) == 0 ) return i;
	return -1;
    }

    String[] unidentified_arguments(String KnownArguments[]) {
	Vector ufos = new Vector();
	Enumeration arge = argv.elements();
	if( arge.hasMoreElements() == false ) 
	    return string_vector_to_string_array(ufos);

        // forget about first element (application name)
	String it = (String)arge.nextElement(); 
	for(; arge.hasMoreElements() ;) {	    
	    it = (String)arge.nextElement();
	    // -- argument belongs to prefixed section ?
	    String arg = __get_remaining_string(it, prefix);
	    if( arg == "" ) continue;

	    // -- check if argument is known
	    if( __search_string_array(KnownArguments, arg) == -1 )
		ufos.addElement(it);
	}    
	return string_vector_to_string_array(ufos);
    }


    String[]  unidentified_options(String KnownOptions[]) {
	Vector ufos = new Vector();

	Enumeration arge = argv.elements();
	if( arge.hasMoreElements() == false )
	    return string_vector_to_string_array(ufos);
        // forget about first element (application name)
	String it = (String)arge.nextElement(); 
	for(; arge.hasMoreElements() ;) {
	    it = (String)arge.nextElement();
	    // -- argument belongs to prefixed section ?
	    String arg = __get_remaining_string(it, prefix);
	    if( arg == "" ) continue;
	    // -- argument == option ?
	    if( arg.length() < 1 || arg.charAt(0) != '-' ) continue;
	    // -- check if argument is known 
	    if( __search_string_array(KnownOptions, arg) == -1 )
		ufos.addElement(it);
	}    
	return string_vector_to_string_array(ufos);
    }

    String  unidentified_flags(String KnownFlagList) {
	return unidentified_flags(KnownFlagList, -1);
    }

    String  unidentified_flags(String KFL, int ArgumentNumber) {
	// Two modes:
	//  ArgumentNumber >= 0 check specific argument 
	//  ArgumentNumber == -1 check all options starting with one '-' 
	//                       for flags

	String ufos = "";
	Vector known_arguments;

	// (2) iteration over '-' arguments (options)
	if( ArgumentNumber == -1 ) {
	    // iterate over all arguments 
	    Enumeration arge = argv.elements();
	    if( arge.hasMoreElements() == false ) return ufos;
	    // forget about first element (application name)
	    String it = (String)arge.nextElement();
	    for(; arge.hasMoreElements() ;) {
		it = (String)arge.nextElement();
		// -- argument belongs to prefixed section ?
		String arg = __get_remaining_string(it, prefix);
		if( arg == "" ) continue;

		// does arguments start with '-' (but not '--')
		if     ( arg.length() < 2 )     continue;
		else if( arg.charAt(0) != '-' ) continue;
		else if( arg.charAt(1) == '-' ) continue;
	    
		// -- check out if flags inside option are contained in KnownFlagList
		for(int i=1; i<arg.length() ; i++) {
		    char flag = arg.charAt(i);
		    if( KFL.indexOf(flag) == -1 ) ufos += flag;
		}
	    }
	}
	else {
	    int no_matches = 0;
	    // -- only check arguments that start with prefix
	    Enumeration arge = argv.elements();
	    if( arge.hasMoreElements() == false ) return ufos;
	    // forget about first element (application name)
	    String it = (String)arge.nextElement();
	    for(; arge.hasMoreElements() ;) {
		it = (String)arge.nextElement();
		String Remain = __get_remaining_string(it, prefix);
		if( Remain != "") {
		    no_matches += 1;
		    if( no_matches == ArgumentNumber) {
			// -- the right argument number inside the section is found
			// => check it for flags
			for(int i=0; i<Remain.length() ; i++) {
			    char flag = Remain.charAt(i);
			    if( KFL.indexOf(flag) == -1 ) ufos += flag;
			}
		    }
		}
	    }
	}
	return ufos;
    }


    String[]   unidentified_variables(String[] Knowns) {
	Vector ufos = new Vector();

	for(Enumeration inm = variables.elements() ; inm.hasMoreElements(); ) {
	    variable it = (variable)inm.nextElement();
	    // -- variable belongs to prefixed section ?
	    String var_name = __get_remaining_string(it.name, prefix);
	    if( var_name == "" ) continue;

	    // -- check if variable is known
	    if( __search_string_array(Knowns, var_name) == -1 )
		ufos.addElement(it.name);
	}
	return string_vector_to_string_array(ufos);
    }


    String[] unidentified_sections(String[] Knowns) {
	Vector ufos = new Vector();

	for(Enumeration inm = section_list.elements() ; inm.hasMoreElements(); ) {
	    String it = (String)inm.nextElement();
	    // -- section a subsection of prefix ?
	    String sec_name = __get_remaining_string(it, prefix);
	    if( sec_name == "" ) continue;

	    // -- section is known ?
	    if( __search_string_array(Knowns, sec_name) == -1 )
		ufos.addElement(it);
	}
	return string_vector_to_string_array(ufos);	
    }


    String[] unidentified_nominuses(String[] KnownNominuses) {
	Vector ufos = new Vector();

	// iterate over all arguments 
	Enumeration arge = argv.elements();
	if( arge.hasMoreElements() == false )
	    return string_vector_to_string_array(ufos);
        // forget about first element (application name)
	String it = (String)arge.nextElement(); 
	for(; arge.hasMoreElements() ;) {
	    it = (String)arge.nextElement();
	    // -- argument belongs to prefixed section ?
	    String arg = __get_remaining_string(it, prefix);
	    if( arg == "" ) continue;

	    if( arg.length() < 1 )                                                continue;
	    // option ? --> not a nomius 
	    if( arg.charAt(0) == '-' )                                            continue;
	    // section ? --> not a real nominus
	    if( arg.charAt(0) == '[' && arg.charAt(arg.length()-1) == ']' ) continue;
	    // variable definition ? --> not a real nominus
	    boolean continue_f = false;
	    for(int i=0; i<arg.length() ; i++)
		if( arg.charAt(i) == '=' ) { continue_f = true; break; }
	    if( continue_f )                   		                             continue;

	    // real nominuses are compared with the given list
	    if( __search_string_array(KnownNominuses, arg) == -1 )
		ufos.addElement(it);
	}    
	return string_vector_to_string_array(ufos);
    }   
}