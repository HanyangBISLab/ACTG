package ExonGraph;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Stack;

import Environments.Constants;

public class GENE implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -1327534019288345239L;

	final int TAGLENGTH = 10;

	private String gene_id;
	// gene_id: the id of gene
	private int trans_cnt;
	// trans_cnt:
	private boolean strand;
	// strand: the option of implement
	// true = 1 && false = -1
	private ArrayList<EXON> trans;
	public String chrID = null;
	/*
	This size matches with trans_cnt.
	This is because there is possibility to get different CDS region each transcript models.
	*/
	public ArrayList<Transcript> transcriptList = new ArrayList<Transcript>();
	/*
	 * Intron
	 * 
	 */
	public Transcript intronscript = new Transcript();

	// nucleotide sequences along with traversing.
	private StringBuilder nucleotides = new StringBuilder();
	
	class StartPos {
		LinkedList<EXON> startExon;
		int startPos;
		String nucleotide;

		public StartPos(StartPos temp){
			this.startExon = new LinkedList<EXON>();
			this.startExon.addAll(temp.startExon);
			this.startPos = temp.startPos;
			this.nucleotide = new String(temp.nucleotide);
		}

		public StartPos() {
			this.startExon = new LinkedList<EXON>();
			this.startPos = 0;
			this.nucleotide = new String();
		}

		public StartPos(int i) {
			this.startExon = new LinkedList<EXON>();
			this.startPos = i;
			this.nucleotide = new String();
		}
	}

	// trans:

	// constructor
	public GENE(String gene_id, int trans_cnt, boolean strand) {
		// Initialize GENE information
		this.gene_id = gene_id;
		this.trans_cnt = trans_cnt;
		// Q. need to new operation?
		this.strand = strand;
	}

	// Add exon to trans (ArrayList<EXON>)
	public void add_trans(EXON exon) { trans.add(exon);}

	// Override to add_trans
	// Add exon to trans where the index
	public void add_trans(int index, EXON exon) { trans.add(index, exon);}

	// Get methods
	public ArrayList<EXON> get_trans() { return trans; }
	public String get_gene_id() { return gene_id; }
	public int get_trans_cnt() { return trans_cnt; }
	public boolean get_strand() { return strand; }

	// initial_exon_graph
	void InitialExonGraph(BufferedReader bufferedReader) {
		int Exon_cnt;
		EXON present, tmp;

		String readString = null;
		for (int j = 0; j < trans_cnt; j++) {
			try {
				readString = bufferedReader.readLine();
				readString = bufferedReader.readLine();
				if (j == 0) {
					int d = Integer.parseInt(readString);
					strand = (d > 0 ? true : false);
				}

				readString = bufferedReader.readLine();
				Exon_cnt = Integer.parseInt(readString);

				if(strand) present = new EXON(null, null, null, -1, -1, 0, this.trans_cnt);
				else present = new EXON(null, null, null, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1, 0, this.trans_cnt);
				this.trans.add(present);

				for (int k = 0; k < Exon_cnt; k++) {
					tmp = new EXON(null, null, null, 0, 0, 0, this.trans_cnt);
					tmp.configure(bufferedReader, strand);
					if (strand) {
						if(k == 0) tmp.set_head(null);
						present.set_next_junc(j, 0, tmp);
						tmp.set_prev_junc(j, 0, present);
					} else {
						if(k == 0) tmp.set_tail(null);
						tmp.set_next_junc(j, 0, present);
						present.set_prev_junc(j, 0, tmp);
					}
					present = tmp;
				}
				if (strand == false) {
					present.set_head(null);

					tmp = new EXON(null, null, null, -1, -1, 0, this.trans_cnt);
					tmp.set_next_junc(j, 0, present);
					present.set_prev_junc(j, 0, tmp);
					present = tmp;

					this.trans.set(j, present);
				} 
				else {
					present.set_tail(null);

					tmp = new EXON(null, null, null, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1, 0, this.trans_cnt);
					present.set_next_junc(j, 0, tmp);
					tmp.set_prev_junc(j, 0, present);
					present = tmp;
				}
			} catch (IOException e) { e.printStackTrace(); }			
		}
	}

	// By SeungHyuk
	public void add_Alternative() {
		int loop;
		EXON prevExon, nextExon;
		// Temporary EXON object

		for (loop = 0; loop < this.trans_cnt; loop++) {
			prevExon = trans.get(loop).get_next(loop, 0);
			nextExon = prevExon.get_next(loop, 0);

			if(nextExon.get_next(loop, 0) == null) continue;
			
			while (nextExon.get_next(loop, 0).get_start() != Integer.MAX_VALUE-1) {
				prevExon.set_next_junc(loop, ExonGraph.JVALUE - 1, nextExon.get_next(loop, 0));

				if (prevExon.get_next(loop, ExonGraph.JVALUE - 1) != null)
					prevExon.get_next(loop, ExonGraph.JVALUE - 1).set_prev_junc(loop, ExonGraph.JVALUE - 1, prevExon);
				else System.out.println("Error");

				prevExon = nextExon;
				nextExon = nextExon.get_next(loop, 0);
			}
		}
	}

	// By Seunghyuk
	public void add_Junction() {
		int loop1, loop2, loopA1, loopA2, loopA3, tJunc, hJunc;
		String tempString;
		EXON prevExon, nextExon, tempExon, tPrev, hPrev;

		String Consensus5, Consensus3;
		if(this.strand) {
			Consensus5 = "GT";
			Consensus3 = "AG";
		}
		else {
			Consensus3 = "TG";
			Consensus5 = "GA";
		}

		for (loop1 = 0; loop1 < this.trans_cnt; loop1++) {// for loop start
			prevExon = trans.get(loop1);
			tPrev = prevExon;
			nextExon = prevExon.get_next(loop1, 0);
			hPrev = nextExon;

			while (nextExon != null) {// while statement start
				tJunc = hJunc = 0;
				loopA2 = loopA3 = 2;

				// INTRON == 12, It defines upper the source code by private
				// final
				// IF1 statement start
				if (nextExon.get_start() - prevExon.get_end() < 12 || prevExon.get_tail() == null /*|| nextExon.get_head() == null*/) {
					prevExon = nextExon;
					nextExon = nextExon.get_next(loop1, 0);
					tPrev = prevExon;
					hPrev = nextExon;
				} 
				else {// IF1 end, and ELSE1 statement start
					for (loop2 = 1; loop2 < 6; loop2++) {
						tempString = prevExon.get_tail().substring(loop2, loop2 + 2);
						// substring(a,b) -> a<= string < b


						// compareTo(String) -> if equal than return 0
						if (tempString.compareTo(Consensus5) == 0 && loop2 == 3) {// IF2
							// statement
							// start
							tempString = prevExon.get_tail().substring(tJunc , loop2);
							tempExon = new EXON(null, null, null, 0, 0, 0, this.trans_cnt);
							tempExon.set_start(tPrev.get_end() + 1);
							tempExon.set_end(tPrev.get_end() + loop2 - tJunc);
							tempExon.set_seq(tempString);

							// Set prev and next each
							// (prev) tPrev <-> tempExon (next)
							tPrev.set_next_junc(loop1, 1, tempExon);
							tempExon.set_prev_junc(loop1, 1, tPrev);

							// The next Exon of tPrev is connecting to the next
							// of tempExon
							// (prev) tempExon <-> nextExon (next)
							tempExon.set_next_junc(loop1, loopA2, nextExon);
							nextExon.set_prev_junc(loop1, loopA2, tempExon);

							loopA2++;
							tJunc = loop2;
							tPrev = tempExon;
						}// IF2 statement end

						//ERROR1 nullpoint
						tempString = nextExon.get_head().substring(5 - loop2, 7 - loop2);
						if (tempString.compareTo(Consensus3) == 0 && loop2 == 3) {// IF3
							// statement
							// start
							tempString = nextExon.get_head().substring(7 - loop2 - hJunc, 7 - 2 * hJunc);
							tempExon = new EXON(null, null, null, 0, 0, 0, this.trans_cnt);
							tempExon.set_end(hPrev.get_start() - 1);
							tempExon.set_start(tempExon.get_end()- (loop2 - hJunc) + 1);
							tempExon.set_seq(tempString);

							hPrev.set_prev_junc(loop1, 1, tempExon);
							tempExon.set_next_junc(loop1, 1, hPrev);

							tempExon.set_prev_junc(loop1, loopA3, prevExon);
							prevExon.set_next_junc(loop1, loopA3, tempExon);
							loopA3++;
							hJunc = loop2;
							hPrev = tempExon;
						}// IF3 statement end
					}// for loop end : for(loop2=1; loop2<6; loop2++)
					prevExon = nextExon;
					nextExon = nextExon.get_next(loop1, 0);
					tPrev = prevExon;
					hPrev = nextExon;
				}// ELSE1 statement end

			}// while statement end
		}// for loop end : for(loop1=0; loop1<this.trans_cnt; loop1++)

		
		for (loop1 = 0; loop1 < this.trans_cnt; loop1++) {
			
			prevExon = trans.get(loop1);
			while (prevExon != null) {
				//loopA2 = loopA3 = 2;
				// edge瑜� 3踰덉쑝濡� 怨좎젙
				loopA2 = loopA3 = 3;
				
				hJunc = tJunc = 0;
				if (prevExon.get_end() - prevExon.get_start() < 11) {
					if (prevExon.get_head() != null) prevExon.set_head(null);
					if (prevExon.get_tail() != null) prevExon.set_tail(null);
					prevExon = prevExon.get_next(loop1, 0);
					continue;
				}

				if (prevExon.get_head() != null) {
					// edge瑜� 3踰덉쑝濡� 怨좎젙
					//while (prevExon.get_prev(loop1,0).get_next(loop1, loopA2) != null) loopA2++;

					for (loop2 = 0; loop2 < 4; loop2++) {
						tempString = prevExon.get_seq().substring(loop2, loop2 + 2);
						if (tempString.compareTo(Consensus3) == 0 && loop2 == 1) {
							tempString = prevExon.get_seq().substring(hJunc, loop2 +2);

							tempExon = new EXON(null, null, null, 0, 0, 0, this.trans_cnt);
							tempExon.set_start(prevExon.get_start() + hJunc);
							tempExon.set_end(tempExon.get_start() + loop2 - hJunc + 2 - 1);
							tempExon.set_seq(tempString);

							tempExon.set_next_junc(loop1,0, prevExon);

							for (loopA1 = 0; loopA1 < ExonGraph.JVALUE; loopA1++) {
								tempExon.set_prev_junc(loop1,loopA1, prevExon.get_prev(loop1, loopA1));
								if (tempExon.get_prev(loop1,loopA1) != null) tempExon.get_prev(loop1,loopA1).set_next_junc(loop1, loopA1, tempExon);
								prevExon.set_prev_junc(loop1, loopA1, null);
							}
							tempExon.get_next(loop1,0).set_prev_junc(loop1, 0, tempExon);

							if (hJunc == 0) prevExon.set_prev_junc(loop1,loopA2, tempExon.get_prev(loop1, 0));
							else prevExon.set_prev_junc(loop1,loopA2, tempExon.get_prev(loop1, loopA2-1));
							prevExon.get_prev(loop1,loopA2).set_next_junc(loop1, loopA2, prevExon);
							loopA2++;
							hJunc = loop2 + 2;
						}
					}
				}

				if (prevExon.get_tail() != null) {
					// edge瑜� 3踰덉쑝濡� 怨좎젙
					//while (prevExon.get_next(loop1,0).get_prev(loop1,loopA3) != null) loopA3++;

					for (loop2 = 0; loop2 < 4; loop2++) {
						tempString = prevExon.get_seq().substring(prevExon.get_end() - prevExon.get_start() - loop2-1, prevExon.get_end() - prevExon.get_start() - loop2 +1);
						if (tempString.compareTo(Consensus5) == 0 && loop2 == 1) {
							tempString = prevExon.get_seq().substring(prevExon.get_end() - prevExon.get_start() - loop2-1, prevExon.get_end() - prevExon.get_start() + 1 - tJunc);
							tempExon = new EXON(null, null, null, 0, 0, 0, this.trans_cnt);
							tempExon.set_end(prevExon.get_end() - tJunc);
							tempExon.set_start(tempExon.get_end() - (loop2 + 2 - tJunc) + 1);
							tempExon.set_seq(tempString);

							tempExon.set_prev_junc(loop1,0, prevExon);

							for (loopA1 = 0; loopA1 < ExonGraph.JVALUE; loopA1++) {
								tempExon.set_next_junc(loop1, loopA1, prevExon.get_next(loop1, loopA1)); //NULL ERROR
								if (tempExon.get_next(loop1, loopA1) != null) tempExon.get_next(loop1, loopA1).set_prev_junc(loop1, loopA1,tempExon);
								prevExon.set_next_junc(loop1, loopA1, null);
							}
							tempExon.get_prev(loop1,0).set_next_junc(loop1,0,tempExon);

							if (tJunc == 0) prevExon.set_next_junc(loop1, loopA3, tempExon.get_next(loop1,0));
							else prevExon.set_next_junc(loop1,loopA3, tempExon.get_next(loop1,loopA3-2));
							prevExon.get_next(loop1,loopA3).set_prev_junc(loop1,loopA3, prevExon);
							loopA3++;
							tJunc = loop2 + 2;
						}
					}
				}

				if (prevExon.get_head() != null) prevExon.set_head(null);
				if (prevExon.get_tail() != null) prevExon.set_tail(null);

				if (hJunc != 0 || tJunc != 0) {
					String temp;
					int start, end;
					if (hJunc != 0) start = prevExon.get_prev(loop1,0).get_end() + 1;
					else start = prevExon.get_start();
					if (tJunc != 0) end = prevExon.get_next(loop1,0).get_start() - 1;
					else end = prevExon.get_end();
					temp = new String(prevExon.get_seq().substring(start - prevExon.get_start(), start - prevExon.get_start() + end	- start+1));
					prevExon.set_seq(temp);
					prevExon.set_start(start);
					prevExon.set_end(end);
				}
				prevExon = prevExon.get_next(loop1,0);
			}
		}
	}// End add_junction

	// ExonGraphDiv By Hyunwoo
	void ExonGraphDiv()
	{
		int i, j, k;
		EXON tmp, pre_exon, pre_compare, tmp2;
		String str;
		int compare_min, compare_min_pos, compare_now_pos;
		int exon_min, exon_min_pos, exon_now_pos;

		//------------------ 占쎈―�뜐�뜝�룞�삕域밟뫁�굲占쎌쥙�쓡野껁깷�쐻占쏙옙占� 癲ル슪�맋�몭�씛�삕�뜝占�----------------- //
		for(i=0;i<this.trans_cnt;i++) {
			pre_exon = trans.get(i);
			while(pre_exon != null) {
				for(j=0;j<this.trans_cnt;j++) {
					if(i == j) continue;
					pre_compare = trans.get(j);
					while(pre_compare != null) {
						if(pre_exon.get_end() < pre_compare.get_start()) break;
						//pre_exon
						else if(pre_exon.get_start() == pre_compare.get_start()) {
							if(pre_exon.get_end() > pre_compare.get_end()) {
								str = new String(pre_exon.get_seq().substring(pre_compare.get_end()+1 - pre_exon.get_start()));
								tmp = new EXON(str, null, null, pre_compare.get_end()+1, pre_exon.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp.set_next_junc(i, k, pre_exon.get_next(i, k));
									if(tmp.get_next(i, k) != null) tmp.get_next(i, k).set_prev_junc(i, k, tmp);
									pre_exon.set_next_junc(i, k, null);
								}

								pre_exon.set_end(pre_compare.get_end());
								str = new String(pre_exon.get_seq().substring(0, pre_exon.get_end() - pre_exon.get_start() + 1));
								pre_exon.set_seq(str);

								if(tmp.get_next(i, 0) != null || (tmp.get_next(i, 0) == null && tmp.get_next(i, 1) == null)) {
									pre_exon.set_next_junc(i, 0, tmp);
									tmp.set_prev_junc(i, 0, pre_exon);
								}
								if(tmp.get_next(i, 1) != null) {
									pre_exon.set_next_junc(i, 1, tmp);
									tmp.set_prev_junc(i, 1, pre_exon);
								}
								break;
							}
							else if(pre_exon.get_end() == pre_compare.get_end()) break;
							else {
								str = new String(pre_compare.get_seq().substring(pre_exon.get_end() +1 - pre_compare.get_start()));
								tmp = new EXON(str, null, null, pre_exon.get_end() + 1, pre_compare.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp.set_next_junc(j, k, pre_compare.get_next(j, k));
									if(tmp.get_next(j, k) != null) tmp.get_next(j, k).set_prev_junc(j, k, tmp);
									pre_compare.set_next_junc(j, k, null);
								}

								pre_compare.set_end(pre_exon.get_end());
								str = new String(pre_compare.get_seq().substring(0, pre_compare.get_end() - pre_compare.get_start()+1));
								pre_compare.set_seq(str);

								if(tmp.get_next(j, 0) != null || (tmp.get_next(j, 0) == null && tmp.get_next(j, 1) == null)) {
									pre_compare.set_next_junc(j, 0, tmp);
									tmp.set_prev_junc(j, 0, pre_compare);
								}
								if(tmp.get_next(j, 1) != null) {
									pre_compare.set_next_junc(j, 1, tmp);
									tmp.set_prev_junc(j, 1, pre_compare);
								}
								break;
							}
						}
						else if(pre_exon.get_start() < pre_compare.get_start()) {
							if(pre_exon.get_end() > pre_compare.get_end()) {

								str = new String(pre_exon.get_seq().substring(pre_compare.get_end()+1 - pre_exon.get_start()));
								tmp2 = new EXON(str, null, null, pre_compare.get_end() + 1, pre_exon.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp2.set_next_junc(i, k, pre_exon.get_next(i, k));
									if(tmp2.get_next(i, k) != null) tmp2.get_next(i, k).set_prev_junc(i, k, tmp2);
									pre_exon.set_next_junc(i, k, null);
								}

								tmp = new EXON(pre_compare.get_seq(), null, null, pre_compare.get_start(), pre_compare.get_end(), 0, this.trans_cnt);

								pre_exon.set_end(pre_compare.get_start() - 1);
								str = new String(pre_exon.get_seq().substring(0, pre_exon.get_end() - pre_exon.get_start()+1));
								pre_exon.set_seq(str);

								if(tmp2.get_next(i, 0) != null || (tmp2.get_next(i, 0) == null && tmp2.get_next(i, 1) == null)) {
									pre_exon.set_next_junc(i, 0, tmp);
									tmp.set_prev_junc(i, 0, pre_exon);

									tmp.set_next_junc(i, 0, tmp2);
									tmp2.set_prev_junc(i, 0, tmp);
								}
								if(tmp2.get_next(i, 1) != null) {
									pre_exon.set_next_junc(i, 1, tmp);
									tmp.set_prev_junc(i, 1, pre_exon);

									tmp.set_next_junc(i, 1, tmp2);
									tmp2.set_prev_junc(i, 1, tmp);
								}
								break;
							}
							else if(pre_exon.get_end() == pre_compare.get_end()) {
								tmp = new EXON(pre_compare.get_seq(), null, null, pre_compare.get_start(), pre_compare.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp.set_next_junc(i, k, pre_exon.get_next(i, k));
									if(tmp.get_next(i, k) != null) tmp.get_next(i, k).set_prev_junc(i, k, tmp);
									pre_exon.set_next_junc(i, k, null);
								}

								pre_exon.set_end(pre_compare.get_start()-1);
								str = new String(pre_exon.get_seq().substring(0, pre_exon.get_end()- pre_exon.get_start()+1));
								pre_exon.set_seq(str);

								if(tmp.get_next(i, 0) != null || (tmp.get_next(i, 0) == null && tmp.get_next(i, 1) == null)) {
									pre_exon.set_next_junc(i, 0, tmp);
									tmp.set_prev_junc(i, 0, pre_exon);
								}
								if(tmp.get_next(i, 1) != null) {
									pre_exon.set_next_junc(i, 1, tmp);
									tmp.set_prev_junc(i, 1, pre_exon);
								}
								break;
							}
							else {
								str = new String(pre_compare.get_seq().substring(0, pre_exon.get_end() - pre_compare.get_start()+1));
								tmp = new EXON(str, null, null, pre_compare.get_start(), pre_exon.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp.set_next_junc(i, k, pre_exon.get_next(i, k));
									if(tmp.get_next(i, k) != null) tmp.get_next(i, k).set_prev_junc(i, k, tmp);
									pre_exon.set_next_junc(i, k, null);
								}

								pre_exon.set_end(tmp.get_start() - 1);
								str = new String(pre_exon.get_seq().substring(0, pre_exon.get_end() - pre_exon.get_start()+1));
								pre_exon.set_seq(str);

								if(tmp.get_next(i, 0) != null || (tmp.get_next(i, 0) == null && tmp.get_next(i, 1) == null)) {
									pre_exon.set_next_junc(i, 0, tmp);
									tmp.set_prev_junc(i, 0, pre_exon);
								}
								if(tmp.get_next(i, 1) != null) {
									pre_exon.set_next_junc(i, 1, tmp);
									tmp.set_prev_junc(i, 1, pre_exon);
								}

								str = new String(pre_compare.get_seq().substring(tmp.get_end() - pre_compare.get_start() + 1));
								tmp2 = new EXON(str, null, null, tmp.get_end() + 1, pre_compare.get_end(), 0, this.trans_cnt);

								for(k=0;k<ExonGraph.JVALUE;k++) {
									tmp2.set_next_junc(j, k, pre_compare.get_next(j, k));
									if(tmp2.get_next(j, k) != null) tmp2.get_next(j, k).set_prev_junc(j, k, tmp2);
									pre_compare.set_next_junc(j, k, null);
								}

								pre_compare.set_end(tmp2.get_start() - 1);
								str = new String(pre_compare.get_seq().substring(0, pre_compare.get_end() - pre_compare.get_start()+1));
								pre_compare.set_seq(str);

								if(tmp2.get_next(j, 0) != null || (tmp2.get_next(j, 0) == null && tmp2.get_next(j, 1) == null)) {
									pre_compare.set_next_junc(j, 0, tmp2);
									tmp2.set_prev_junc(j, 0, pre_compare);
								}
								if(tmp2.get_next(j, 1) != null) {
									pre_compare.set_next_junc(j, 1, tmp2);
									tmp2.set_prev_junc(j, 1, pre_compare);
								}
								break;
							}
						}
						else {
							compare_now_pos = pre_compare.get_start();
							compare_min = Integer.MAX_VALUE;
							compare_min_pos = -1;
							for(k=0;k<3;k++) {
								if(pre_compare.get_next(j, k) != null) {
									if(pre_compare.get_next(j, k).get_start() < compare_min) {
										compare_min = pre_compare.get_next(j, k).get_start();
										compare_min_pos = k;
									}
								}
							}
							if(compare_min_pos == -1) pre_compare = null;
							else {
								pre_compare = pre_compare.get_next(j, compare_min_pos);
								while(pre_compare.get_prev(j, 1) != null) {
									if(pre_compare.get_prev(j, 1).get_start() > compare_now_pos) pre_compare = pre_compare.get_prev(j, 1);
									else break;
								}
							}
						}
					}
				}
				exon_now_pos = pre_exon.get_start();
				exon_min = Integer.MAX_VALUE;
				exon_min_pos = -1;
				for(k=0;k<3;k++) {
					if(pre_exon.get_next(i, k) != null) {
						if(pre_exon.get_next(i, k).get_start() < exon_min) {
							exon_min = pre_exon.get_next(i, k).get_start();
							exon_min_pos = k;
						}
					}
				}
				if(exon_min_pos == -1) pre_exon = null;
				else {
					pre_exon = pre_exon.get_next(i, exon_min_pos);
					while(pre_exon.get_prev(i, 1) != null) {
						if(pre_exon.get_prev(i, 1).get_start() > exon_now_pos) pre_exon = pre_exon.get_prev(i, 1);
						else break;
					}
				}
			}
		}
	}


	// ExonGraphMerge By Hyunwoo
	public void ExonGraphMerge() {
		int i, j, k;
		EXON pre_exon, pre_compare;
		int compare_min, compare_min_pos, compare_now_pos;
		int exon_min, exon_min_pos, exon_now_pos;
		for (i = 0; i < this.trans_cnt; i++) {
			pre_exon = this.trans.get(i);
			while (pre_exon != null) {
				for (j = 0; j < this.trans_cnt; j++) {
					pre_compare = this.trans.get(j);
					while (pre_compare != null) {
						if (pre_exon.get_start() == pre_compare.get_start()
								&& pre_exon.get_end() == pre_compare.get_end()) {
							/*
							if (pre_exon.get_seq() != null && pre_compare.get_seq() != null) {
								if(pre_exon.get_seq().compareTo(pre_compare.get_seq()) == 0) ;// printf("OK\n");
								else System.out.println("No");
							}*/

							if (pre_exon == pre_compare) break;

							// if(pre_compare->prev[0][j] == NULL &&
							// pre_compare->prev[1][j] == NULL) {
							if (pre_compare == this.trans.get(j)) {
								for (k = 0; k < ExonGraph.JVALUE; k++) {
									pre_exon.set_next_junc(j, k, pre_compare.get_next(j, k));
									pre_exon.set_prev_junc(j, k, pre_compare.get_prev(j, k));
									if (pre_compare.get_next(j, k) != null) pre_compare.get_next(j, k).set_prev_junc(j, k, pre_exon);
									if (pre_compare.get_prev(j, k) != null) pre_compare.get_prev(j, k).set_next_junc(j, k, pre_exon);
								}
								this.trans.set(j, pre_exon);
							}
							else {
								for (k = 0; k < ExonGraph.JVALUE; k++) {
									pre_exon.set_next_junc(j, k, pre_compare.get_next(j, k));
									pre_exon.set_prev_junc(j, k, pre_compare.get_prev(j, k));
									if (pre_compare.get_next(j, k) != null) pre_compare.get_next(j, k).set_prev_junc(j, k, pre_exon);
									if (pre_compare.get_prev(j, k) != null) pre_compare.get_prev(j, k).set_next_junc(j, k, pre_exon);
								}

							}
							break;
						} else {
							compare_now_pos = pre_compare.get_start();
							compare_min = Integer.MAX_VALUE;
							compare_min_pos = -1;
							for (k = 0; k < ExonGraph.JVALUE - 1; k++) {
								if (pre_compare.get_next(j, k) != null) {
									if (pre_compare.get_next(j, k).get_start() < compare_min) {
										compare_min = pre_compare.get_next(j, k).get_start();
										compare_min_pos = k;
									}
								}
							}
							if (compare_min_pos == -1) pre_compare = null;
							else {
								pre_compare = pre_compare.get_next(j,
										compare_min_pos);
								while (pre_compare.get_prev(j, 1) != null) {
									if (pre_compare.get_prev(j, 1).get_start() > compare_now_pos) pre_compare = pre_compare.get_prev(j, 1);
									else break;
								}
							}
						}
					}
				}
				exon_now_pos = pre_exon.get_start();
				exon_min = Integer.MAX_VALUE;
				exon_min_pos = -1;
				for (k = 0; k < ExonGraph.JVALUE - 1; k++) {
					if (pre_exon.get_next(i, k) != null) {
						if (pre_exon.get_next(i, k).get_start() < exon_min) {
							exon_min = pre_exon.get_next(i, k).get_start();
							exon_min_pos = k;
						}
					}
				}
				if (exon_min_pos == -1)
					pre_exon = null;
				else {
					pre_exon = pre_exon.get_next(i, exon_min_pos);
					while (pre_exon.get_prev(i, 1) != null) {
						if (pre_exon.get_prev(i, 1).get_start() > exon_now_pos) pre_exon = pre_exon.get_prev(i, 1);
						else break;
					}
				}
			}
		}
	}

	void RemoveEdges() {
		for(EXON i : trans) this.RemoveEdges(i);
	}

	void RemoveEdges(EXON present) {		
		if (present == null || present.get_check()) return;

		present.set_check(true);
		for (int j = 0; j < trans_cnt; j++) {
			for (int l = 0; l < ExonGraph.JVALUE; l++) {
				for (int k = 0; k < trans_cnt; k++) {
					for (int m = 0; m < ExonGraph.JVALUE; m++) {
						if (k == j && l == m) continue;
						if (present.get_next(j, l) == present.get_next(k, m)) {
							if (m != 0) present.set_next_junc(k, m, null);
							else present.set_next_junc(j, l, null);
						}
						if (present.get_prev(j, l) == present.get_prev(k, m)) {
							if (m != 0) present.set_prev_junc(k, m, null);
							else present.set_prev_junc(j, l, null);
						}
					}
				}
			}
		}

		for(int i=0; i < trans_cnt; i++) {
			for (int p=0; p < ExonGraph.JVALUE; p++) RemoveEdges(present.get_next(i, p));
		}
	}


	

	public void AminoExonSearch(boolean visited, KAminoTree ac_tree) {
			
		//candidate = new ArrayList<Candidate>();
		for(int loop=0 ; loop<this.trans_cnt; loop++) {
			if(strand == true) {
				AminoExonSearch(ac_tree.get_tree(), trans.get(loop), 0, new StartPos(0), visited, 0);
				AminoExonSearch(ac_tree.get_tree(), trans.get(loop), 0, new StartPos(1), visited, 1);
				AminoExonSearch(ac_tree.get_tree(), trans.get(loop), 0, new StartPos(2), visited, 2);
			}
			else {
				AminoExonSearch(ac_tree.get_rtree(), trans.get(loop), 0, new StartPos(0), visited, 0);
				AminoExonSearch(ac_tree.get_rtree(), trans.get(loop), 0, new StartPos(1), visited, 1);
				AminoExonSearch(ac_tree.get_rtree(), trans.get(loop), 0, new StartPos(2), visited, 2);
			}
		}
	}

	int AminoExonSearch(KAminoNode ac_tree, EXON exon ,int prePos, StartPos startPos, boolean visited, int frame) {
		
		this.nucleotides.setLength(0);
		
		Stack<EXON> exonStack = new Stack<EXON>();
		Stack<Integer> treePosStack = new Stack<Integer>();
		Stack<Integer> frameStack = new Stack<Integer>();
		Stack<StartPos> startPosStack = new Stack<StartPos>();
		
		// stack init
		exonStack.add(exon);
		treePosStack.add(prePos);
		startPosStack.add(startPos);
		frameStack.add(frame);
		
		while(!exonStack.isEmpty()) {
			
			EXON targetExon = exonStack.pop();
			int targetPos = treePosStack.pop();
			StartPos targetStartPos = startPosStack.pop();
			int targetFrame = frameStack.pop();
			
			boolean skipExon = false;
			
			int loop1;
			char nucl;
			int textPos, treePos, framePos;
			String exonSeq = targetExon.get_seq();
			LinkedList<EXON> present = null;
			
			if(targetStartPos.startPos == 0 && targetExon.getFrame()[0] == visited && targetStartPos.startExon.size() == 0) skipExon = true;
			if(targetFrame != 0 && targetExon.getFrame()[3-targetFrame] == visited) skipExon = true;
			
			if(skipExon) continue;
			
			framePos = targetStartPos.nucleotide.length();

			textPos=0+targetFrame;
			treePos=targetPos;
			
			// The end of transcript translates one more
			boolean isMoreTranslatable = false;
			if(Constants.TRANSLATION_END_OF_TRANSCRIPT){
				if(this.strand){
					if(targetExon.get_start() == Integer.MAX_VALUE-1){
						if(targetStartPos.nucleotide.length() == 2){
							targetStartPos.nucleotide += "N";
							nucleotides.append("N");
							nucl = Codon.NuclToAmino(targetStartPos.nucleotide);
							treePos = ac_tree.get(treePos).get_next(nucl);
							isMoreTranslatable = true;
						}
					}
				}else{
					if(targetExon.get_start() == -1 && targetFrame == 2){
						targetStartPos.nucleotide = "N";
						nucleotides.append("N");
					}
				}
				
				if(isMoreTranslatable){
					if(ac_tree.get(treePos).get_matchlist() != null){
						Pattern TPM = ac_tree.get(treePos).get_matchlist();
						addMapping(targetStartPos, textPos, TPM);					
					}
				}
			}
			
			
			// The end of transcript translates one more END//
			
			if(exonSeq != null) {
				targetStartPos.startExon.add(targetExon);
				present = targetStartPos.startExon;
				
				while(textPos < exonSeq.length()) {
					nucleotides.append(exonSeq.charAt(textPos));
					
					if(targetStartPos.nucleotide.length() != 2 ) {
						targetStartPos.nucleotide += exonSeq.charAt(textPos);
						textPos++;
					}
					else {
						targetStartPos.nucleotide += exonSeq.charAt(textPos);
						nucl = 'X';
						try{
						if(this.strand) nucl = Codon.NuclToAmino(targetStartPos.nucleotide);
						else nucl = Codon.NuclToAmino_R(targetStartPos.nucleotide);
						}catch(Exception e){
							System.out.println(targetStartPos.nucleotide);
							System.out.println(nucl);
						}
						
						if(ac_tree.get(treePos).get_next(nucl) == ac_tree.get(ac_tree.get(treePos).get_failstate()).get_next(nucl)) {
							targetStartPos.startPos += (ac_tree.get(treePos).get_depth() - ac_tree.get(ac_tree.get(treePos).get_next(nucl)).get_depth() + 1)*3;
							
							while(targetStartPos.startPos >= present.get(0).get_seq().length()) {
								targetStartPos.startPos -= present.get(0).get_seq().length();
								present.remove(0);
								if(present.size() == 0) {
									targetStartPos.nucleotide = "";
									break;
								}
								
								if(present.get(0).getFrame()[(3 - (targetStartPos.startPos%3))%3] == visited && targetStartPos.startPos < present.get(0).get_seq().length()) {
									skipExon = true;
									break;
								}

							}
						}
						
						if(skipExon) break;
						
						treePos = ac_tree.get(treePos).get_next(nucl);

						// GFF�� NextSearch Result�뙆�씪�쓣 留뚮뱶�뒗 猷⑦떞
						if(ac_tree.get(treePos).get_matchlist() != null){
							Pattern TPM = ac_tree.get(treePos).get_matchlist();
							addMapping(targetStartPos, textPos, TPM);
							
						}
						textPos++;
						targetStartPos.nucleotide = "";
					}
				}
			}
			
			if(skipExon) continue;
			
			for(loop1=0; loop1<this.trans_cnt; loop1++){
				for(int jEdge : Environments.Constants.SEARCH_EDGES){
					if(targetExon.get_next(loop1,jEdge) != null) {
						exonStack.add(targetExon.get_next(loop1,jEdge));
						treePosStack.add(treePos);
						
						if(exonSeq != null) {
							if(exonSeq.length() <= targetFrame) {
								startPosStack.add(new StartPos(targetFrame - exonSeq.length()));
								frameStack.add(targetFrame - exonSeq.length());
							}
							else {
								startPosStack.add(new StartPos(targetStartPos));
								frameStack.add(0);
							}
						}
						else {
							startPosStack.add(new StartPos(targetFrame));
							frameStack.add(targetFrame);
						}
					}
				}
			}

			if (targetFrame != 0) targetExon.setFrame(3-targetFrame, visited);
			else targetExon.setFrame(framePos, visited);
		}

		return 0;
	}


	
	//
	//RefSeqExonGraph Addition
	public void set_trans(ArrayList<EXON> trans){
		this.trans = trans;
		
	}
		
	public void set_trans_cnt(int cnt){
		this.trans_cnt = cnt;
	}
	
	public void incre_trans_cnt(){
		this.trans_cnt++;		
	}
	//END
	
	public void set_intron(){
		ArrayList<ExonRangeType> allExon = new ArrayList<ExonRangeType>();
		
		for(int i=0; i<trans_cnt; i++){
			
			ArrayList<ExonRangeType> exonList = transcriptList.get(i).exonList;
			int sizeOfExonList = exonList.size();
			for(int j=0; j<sizeOfExonList; j++){
				allExon.add(exonList.get(j));
			}
		}
		//Sort
		Collections.sort(allExon, new TranscriptComparator());
		
		//Init intron sites (threat it as pseudo)
		intronscript.transcriptID = "INTRON";
		
		int sizeOfallExon = allExon.size();
		int start = allExon.get(0).start;
		int end = allExon.get(0).end;
		for(int i=1; i<sizeOfallExon; i++){
			if(start <= allExon.get(i).start && end >= allExon.get(i).start){
				
				if(end < allExon.get(i).end){
					end = allExon.get(i).end;
				}
				
			}else if(start < allExon.get(i).start){
				//end+1 == allExon.get(i).start �� 媛숈쑝硫�,
				//intron援ш컙�씠 0�씠誘�濡� �젣�쇅�븳�떎.
				if(end+1 != allExon.get(i).start){
					ExonRangeType ERT = new ExonRangeType();
					ERT.start = end+1;
					ERT.end = allExon.get(i).start-1;
					ERT.isCDS = false;
					ERT.isMut = false;
					intronscript.exonList.add(ERT);
				}
				
				start = allExon.get(i).start;
				end = allExon.get(i).end;
			}
		}
		
		//DONE!
	}
	
	//Version 1.08
	public void addMapping(StartPos startPos, int textPos, Pattern TPM){
		int longestPeptideLength = 0;
		for(int i=0; i<=TPM.size(); i++){
			if(i == 0){
				longestPeptideLength = TPM.get_amino().length();
			}else if(longestPeptideLength < TPM.get(i-1).get_amino().length()){
				longestPeptideLength = TPM.get(i-1).get_amino().length();
			}
		}
		
		
		for(int i=0; i<=TPM.size(); i++){
			Pattern curTPM = null;
			StartPos GFFStartPos = new StartPos(startPos);
			LinkedList<EXON> exonList = GFFStartPos.startExon;
			
			
			if(i == 0){
				curTPM = TPM;
			}else{
				// 寃뱀퀜�꽌 留듯븨�릺�뒗 寃쎌슦
				// ex-
				// AAAACC
				//    ACC
				
				curTPM = TPM.get(i-1);
			}
			
			try{
				
				int interval = 0;
				int startLoci = 0;
				int endLoci = 0;
				int ntLength = curTPM.get_amino().length()*3;
				int len = this.nucleotides.length();
				StringBuilder sequence = new StringBuilder();
				if(this.strand) {
					int index = len - ntLength;
					for(; index<len; index++) {
						sequence.append(this.nucleotides.charAt(index));
					}
				} else {
					// reverse-complementary
					int index = len - ntLength;
					for(; index<len; index++) {
						switch(this.nucleotides.charAt(index)) {
						case 'A': sequence.append('T'); break;
						case 'C': sequence.append('G'); break;
						case 'T': sequence.append('A'); break;
						case 'G': sequence.append('C'); break;
						default : sequence.append('N'); break;
						}
						sequence = sequence.reverse();
					}
				}
				
				for(int j=exonList.size()-1; j>=0; j--){
					EXON tempExon = exonList.get(j);
					if(j == exonList.size()-1){
						endLoci = tempExon.get_start() + textPos;
					}else{
						endLoci = tempExon.get_end();
					}
					
					startLoci = tempExon.get_start();
					
					interval += endLoci - startLoci + 1;
					
					if(ntLength <= interval){
						GFFStartPos.startPos = (interval - ntLength);
						for(int k=0; k<j; k++){
							exonList.removeFirst();
						}
						break;
					}
				}
				
				Flat.write(GFFStartPos, textPos, this.strand, this.trans_cnt, curTPM.getOutput(), transcriptList, gene_id, this.chrID, sequence.toString());
			}catch(Exception E){
				
			}
		}
	}
	

}