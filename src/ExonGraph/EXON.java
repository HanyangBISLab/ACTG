package ExonGraph;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

public class EXON implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 4221958241947364562L;
	public boolean isExon = true;
	public boolean state;


	private String seq, head, tail;
	// seq: string contained this exon
	// head: connect a next exon
	// tail: connect a prev exon
	private int start, end, number;
	// start: starting point this exon
	// end: ending point this exon
	// number: this exon's number
	private boolean check; // Visit checking
	private ArrayList<ArrayList<EXON>> next, prev; // About junction variables
	
	boolean[] frame;

	// constructor
	public EXON(String seq, String head, String tail, int start, int end,
			int number, int trans_cnt) {
		this.state = false;
		// Initialize EXON information
		this.start = start;
		this.end = end;
		this.number = number;
		this.seq = seq;
		this.head = head;
		this.tail = tail;
		this.check = false;

		next = new ArrayList<ArrayList<EXON>>();
		prev = new ArrayList<ArrayList<EXON>>();

		// Initialize next and prev
		for (int i = 0; i < trans_cnt; i++) {// for loop 1 start

			ArrayList<EXON> tempNext = new ArrayList<EXON>();
			ArrayList<EXON> tempPrev = new ArrayList<EXON>();
			// tempNext and tempPrev are containing junction information

			// Initialize tempNext and tempPrev
			for (int j = 0; j < ExonGraph.JVALUE; j++) {
				tempNext.add(null);
				tempPrev.add(null);
				
			}

			// Add tempNext and tempPrev into next and prev each other

			next.add(tempNext);
			prev.add(tempPrev);
		}// for loop 1 end
		
		this.frame = new boolean[3];
		for(int i=0;i<3;i++) frame[i] = false;
	}

	
	
	public int getNumber() { return number; }
	public void setNumber(int number) { this.number = number;}

	// Remove methods for a single junction variable(next and prev)
	public void remove_next_junc(int indexNext, int indexJunc) {
		next.get(indexNext).remove(indexJunc);
	}

	public void remove_prev_junc(int indexPrev, int indexJunc) {
		prev.get(indexPrev).remove(indexJunc);
	}

	// Override to remove_next_junc and remove_prev_junc
	// Remove methods for a set of junction variables(next and prev)
	public void remove_next_junc(int indexNext) {
		next.remove(indexNext);
	}

	public void remove_prev_junc(int indexPrev) {
		prev.remove(indexPrev);
	}

	// Get methods (start, end, number, seq, head, tail, check, prev, and next)
	public int get_start() {
		return start;

	}

	public int get_end() {
		return end;
	}

	public int get_number() {
		return number;
	}

	public String get_seq() {
		return seq;
	}

	public String get_head() {
		return head;
	}

	public String get_tail() {
		return tail;
	}

	public boolean get_check() {
		return check;
	}

	// override
	public EXON get_next(int index1, int index2) {
		return next.get(index1).get(index2);
	}

	public ArrayList<EXON> get_next(int index) {
		return next.get(index);
	}
	
	public ArrayList<ArrayList<EXON>> get_next() {
		return next;
	}

	// override
	public EXON get_prev(int index1, int index2) {
		return prev.get(index1).get(index2);
	}
	
	public ArrayList<EXON> get_prev(int index) {
		return prev.get(index);
	}
	
	public ArrayList<ArrayList<EXON>> get_prev() {
		return prev;
	}

	// Set methods (check, next, and prev)
	public void set_check(boolean check) {
		this.check = check;
	}

	// Q. We need to catch the exception when indexPrev or indexJunc exceeds the
	// limit variable?
	// Include a single junction variable into next
	public void set_next_junc(int indexNext, int indexJunc, EXON junc) {
		next.get(indexNext).set(indexJunc, junc);
	}

	// Include a single junction variable into prev
	public void set_prev_junc(int indexPrev, int indexJunc, EXON junc) {
		prev.get(indexPrev).set(indexJunc, junc);
	}

	// Override to set_next_junc
	// Include a set of junction variables into next
	public void set_next_junc(int indexNext, ArrayList<EXON> juncSet) {
		next.set(indexNext, juncSet);
	}

	// Override to set_prev_junc
	// Include a set of junction variables into prev
	public void set_prev_junc(int indexPrev, ArrayList<EXON> juncSet) {
		prev.set(indexPrev, juncSet);
	}
	
	public void set_start(int start){
		this.start = start;
	}
	
	public void set_end(int end){
		this.end = end;
	}
	
	public void set_seq(String seq){
		this.seq = seq;
	}
	
	public void set_head(String head){
		this.head = head;
	}
	
	public void set_tail(String tail){
		this.tail = tail;
	}

	void configure(BufferedReader bufferedReader, boolean strand)
			throws IOException {
		String readString = null;
		readString = bufferedReader.readLine();
		start = Integer.parseInt(readString);
		readString = bufferedReader.readLine();
		end = Integer.parseInt(readString);

		head = bufferedReader.readLine();
		seq = bufferedReader.readLine();
		tail = bufferedReader.readLine();
		// present = Exon_tmp.get(j);

		if (strand == false) {
			seq = new StringBuilder(seq).reverse().toString();
			head = new StringBuilder(head).reverse().toString();
			tail = new StringBuilder(tail).reverse().toString();
		} //else
			//head = null;

	}


	public boolean[] getFrame() {
		return frame;
	}


	public void setFrame(int i, boolean frame) {
		//if(this.frame[i] == frame) System.out.println("asdfasdf");
		this.frame[i] = frame;
	}

	/*
	 * 
	 * 
	 * RefSeqExonGraph Addition
	 */
	// Include a single junction variable into next
	public void set_next_junc_dynamic(int indexNext, int indexJunc, EXON junc) {
		while( indexNext >= next.size()){
			ArrayList<EXON> tmp = new ArrayList<EXON>();
			for( int i=0; i< ExonGraph.JVALUE; i++){
				tmp.add(null);
			}
			next.add(tmp);
		}
		next.get(indexNext).set(indexJunc, junc);
	}

	// Include a single junction variable into prev
	public void set_prev_junc_dynamic(int indexPrev, int indexJunc, EXON junc) {
		while( indexPrev >= prev.size()){
			ArrayList<EXON> tmp = new ArrayList<EXON>();
			for( int i=0; i< ExonGraph.JVALUE; i++){
				tmp.add(null);
			}
			prev.add(tmp);
		}
		prev.get(indexPrev).set(indexJunc, junc);
	}
	//END
}
