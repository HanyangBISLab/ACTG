package data;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.TreeMap;

class Variant {
	public String chr;
	public int pos;
	public String ref;
	public String alt;
}

public class Variants {

	public Hashtable<String, TreeMap<Integer, ArrayList<Variant>>> variants = new Hashtable<String, TreeMap<Integer, ArrayList<Variant>>>();
	
	/**
	 * [start, end] both are inclusive.
	 * 
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public ArrayList<Variant> getVariants (String chr, int start, int end) {
		ArrayList<Variant> selectedVariants = null;
		
		TreeMap<Integer, ArrayList<Variant>> variants = this.variants.get(chr);
		
		return selectedVariants;
	}
}
