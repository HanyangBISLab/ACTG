package ExonGraph;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;

class Result implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 113834286059755650L;
	String flatResult = null;
	ArrayList<String> GFFResult = null;
	ArrayList<String> logResult = new ArrayList<String>();
	String chr = null;
}

public class ResultMapping implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 2764314750344313208L;
	public static Hashtable <String, Result> ResultTable = new Hashtable<String, Result>();
	
}