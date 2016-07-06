package ExonGraph;

import java.util.ArrayList;

public class KAminoNode extends ArrayList<KAminoNode> {//ArrayList<KAminoNode> implements KNode {
	/**
	 * 
	 */
	private static final long serialVersionUID = -5372755372154803342L;
	int[] NextState;
	int Depth;
	int FailState;
	Pattern MatchList;

	public KAminoNode(int size) {
		NextState = new int[size];
		for(int i=0; i<size; i++) {
			this.NextState[i] = -1;
		}
		Depth = 0;
		MatchList = null;
	}
	
	public int get_next(String a) {
		return get_next(a.charAt(0));
	}
	public int get_next(char a) {
		return this.NextState[a - 'A'];
	}
	public int get_next(int i) {
		return this.NextState[i];
	}
	
	public void set_next(String a, int i) {
		set_next(a.charAt(0), i);
	}
	public void set_next(char a, int i) {
		this.NextState[a - 'A'] = i;
	}
	public void set_next(int i, int j) {
		this.NextState[i] = j;
	}

	public int get_depth() { return this.Depth; }
	public int get_failstate() { return this.FailState; }
	public Pattern get_matchlist() { return this.MatchList; }
	
	public void set_matchlist(Pattern a) {
		if(MatchList == null) MatchList = new Pattern(a);
		else MatchList.add(a);
	}
	public void set_depth(int a) { this.Depth = a; }
	public void set_failstate(int a) { this.FailState = a; }

}
