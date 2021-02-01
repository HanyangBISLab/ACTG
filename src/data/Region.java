package data;

import java.util.ArrayList;

public class Region implements Comparable<Region> {
	
	public final static int UNDEFINED = -1;
	public final static int CDS = 0;
	public final static int INTRON = 1;
	public final static int UTR5 = 2;
	public final static int UTR3 = 3;
	public final static int NONCODING = 4;
	
	public int regionClass = UNDEFINED;
	public int start = 0;
	public int end = 0;
	
	public StringBuilder nucleotides = null;
	public Region nextNormal = null;
	public Region prevNormal = null;
	
	public Region nextExonSkipping = null;
	public Region prevExonSkipping = null;
	
	public ArrayList<Region> nextVars = new ArrayList<Region>();
	public ArrayList<Region> prevVars = new ArrayList<Region>();
	
	
	
	/**
	 * [ start, end ]: both loci are inclusive.
	 * @param start
	 * @param end
	 */
	public Region (int start, int end) {
		this.start = start;
		this.end = end;
	}

	@Override
	public int compareTo(Region o) {
		if(this.start < o.start ) {
			return -1;
		}else if(this.start > o.start) {
			return 1;
		}
		return 0;
	}
	
}
