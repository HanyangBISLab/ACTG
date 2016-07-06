package proteinDB;

import java.util.ArrayList;
import java.util.LinkedList;

import Environments.Constants;
import ExonGraph.Protein;

public class ProteinMappingWithSAVs {
	
	public static final int F = 'F'-'A';
	public static final int L = 'L'-'A';
	public static final int I = 'I'-'A';
	public static final int M = 'M'-'A';
	public static final int V = 'V'-'A';
	public static final int S = 'S'-'A';
	public static final int P = 'P'-'A';
	public static final int T = 'T'-'A';
	public static final int A = 'A'-'A';
	public static final int Y = 'Y'-'A';
	public static final int H = 'H'-'A';
	public static final int Q = 'Q'-'A';
	public static final int N = 'N'-'A';
	public static final int K = 'K'-'A';
	public static final int D = 'D'-'A';
	public static final int E = 'E'-'A';
	public static final int C = 'C'-'A';
	public static final int W = 'W'-'A';
	public static final int R = 'R'-'A';
	public static final int G = 'G'-'A';


	private calculateED cal = null;
	
	public static boolean[][] SNPtable = new boolean[26][26]; 

	public ProteinMappingWithSAVs(){
		createSNPTable();
	}
	
	public void proteinDBMappingWithSAVs(LinkedList<String> queryString, ArrayList<Protein> proteinDB, ProteinDBFlat mappedPeptides){

		for(String queryList : queryString){
			for(int j = 0; j < proteinDB.size(); j++){
				cal = new calculateED(proteinDB.get(j).proteinSeq, queryList);
				int isSNP = cal.calculate();

				if(isSNP == 1){
					mappedPeptides.putMapping(proteinDB.get(j).header, queryList);
				}
			}
		}
	}
	
	private void createSNPTable() {
		for(int i=0; i<SNPtable.length; i++){
			for(int j=0; j<SNPtable[i].length; j++){
				SNPtable[i][j] = false;
			}
		}
		
		SNPtable[F] = new boolean[26];
		SNPtable[F]['S'-'A'] = true;
		SNPtable[F]['Y'-'A'] = true;
		SNPtable[F]['C'-'A'] = true;
		SNPtable[F]['L'-'A'] = true;
		SNPtable[F]['I'-'A'] = true;
		SNPtable[F]['M'-'A'] = true;
		SNPtable[F]['V'-'A'] = true;
		
		SNPtable[L] = new boolean[26]; 
		SNPtable[L]['S'-'A'] = true;
		SNPtable[L]['W'-'A'] = true;
		SNPtable[L]['F'-'A'] = true;
		SNPtable[L]['L'-'A'] = true;
		SNPtable[L]['I'-'A'] = true;
		SNPtable[L]['M'-'A'] = true;
		SNPtable[L]['V'-'A'] = true;
		SNPtable[L]['P'-'A'] = true;
		SNPtable[L]['H'-'A'] = true;
		SNPtable[L]['Q'-'A'] = true;
		SNPtable[L]['R'-'A'] = true;

		SNPtable[I] = new boolean[26];
		SNPtable[I]['F'-'A'] = true;
		SNPtable[I]['L'-'A'] = true;
		SNPtable[I]['M'-'A'] = true;
		SNPtable[I]['V'-'A'] = true;
		SNPtable[I]['T'-'A'] = true;
		SNPtable[I]['N'-'A'] = true;
		SNPtable[I]['K'-'A'] = true;
		SNPtable[I]['S'-'A'] = true;
		SNPtable[I]['R'-'A'] = true;

		SNPtable[M] = new boolean[26];
		SNPtable[M]['F'-'A'] = true;
		SNPtable[M]['L'-'A'] = true;
		SNPtable[M]['I'-'A'] = true;
		SNPtable[M]['V'-'A'] = true;
		SNPtable[M]['T'-'A'] = true;
		SNPtable[M]['K'-'A'] = true;
		SNPtable[M]['R'-'A'] = true;

		SNPtable[V] = new boolean[26];
		SNPtable[V]['F'-'A'] = true;
		SNPtable[V]['L'-'A'] = true;
		SNPtable[V]['I'-'A'] = true;
		SNPtable[V]['M'-'A'] = true;
		SNPtable[V]['A'-'A'] = true;
		SNPtable[V]['D'-'A'] = true;
		SNPtable[V]['E'-'A'] = true;
		SNPtable[V]['G'-'A'] = true;

		SNPtable[S] = new boolean[26];
		SNPtable[S]['F'-'A'] = true;
		SNPtable[S]['L'-'A'] = true;
		SNPtable[S]['Y'-'A'] = true;
		SNPtable[S]['C'-'A'] = true;
		SNPtable[S]['W'-'A'] = true;
		SNPtable[S]['P'-'A'] = true;
		SNPtable[S]['T'-'A'] = true;
		SNPtable[S]['A'-'A'] = true;
		SNPtable[S]['N'-'A'] = true;
		SNPtable[S]['T'-'A'] = true;
		SNPtable[S]['I'-'A'] = true;
		SNPtable[S]['R'-'A'] = true;
		SNPtable[S]['G'-'A'] = true;

		SNPtable[P] = new boolean[26]; 
		SNPtable[P]['L'-'A'] = true;
		SNPtable[P]['H'-'A'] = true;
		SNPtable[P]['Q'-'A'] = true;
		SNPtable[P]['R'-'A'] = true;
		SNPtable[P]['S'-'A'] = true;
		SNPtable[P]['T'-'A'] = true;
		SNPtable[P]['A'-'A'] = true;

		SNPtable[T] = new boolean[26];
		SNPtable[T]['I'-'A'] = true;
		SNPtable[T]['M'-'A'] = true;
		SNPtable[T]['N'-'A'] = true;
		SNPtable[T]['K'-'A'] = true;
		SNPtable[T]['S'-'A'] = true;
		SNPtable[T]['R'-'A'] = true;
		SNPtable[T]['P'-'A'] = true;
		SNPtable[T]['A'-'A'] = true;

		SNPtable[A] = new boolean[26];
		SNPtable[A]['V'-'A'] = true;
		SNPtable[A]['D'-'A'] = true;
		SNPtable[A]['E'-'A'] = true;
		SNPtable[A]['G'-'A'] = true;
		SNPtable[A]['S'-'A'] = true;
		SNPtable[A]['P'-'A'] = true;
		SNPtable[A]['T'-'A'] = true;

		SNPtable[Y] = new boolean[26];
		SNPtable[Y]['F'-'A'] = true;
		SNPtable[Y]['L'-'A'] = true;
		SNPtable[Y]['S'-'A'] = true;
		SNPtable[Y]['C'-'A'] = true;
		SNPtable[Y]['W'-'A'] = true;
		SNPtable[Y]['H'-'A'] = true;
		SNPtable[Y]['Q'-'A'] = true;
		SNPtable[Y]['N'-'A'] = true;
		SNPtable[Y]['K'-'A'] = true;
		SNPtable[Y]['D'-'A'] = true;
		SNPtable[Y]['E'-'A'] = true;

		SNPtable[H] = new boolean[26];
		SNPtable[H]['L'-'A'] = true;
		SNPtable[H]['P'-'A'] = true;
		SNPtable[H]['R'-'A'] = true;
		SNPtable[H]['Y'-'A'] = true;
		SNPtable[H]['Q'-'A'] = true;
		SNPtable[H]['N'-'A'] = true;
		SNPtable[H]['K'-'A'] = true;
		SNPtable[H]['D'-'A'] = true;
		SNPtable[H]['E'-'A'] = true;

		SNPtable[Q] = new boolean[26];
		SNPtable[Q]['L'-'A'] = true;
		SNPtable[Q]['P'-'A'] = true;
		SNPtable[Q]['R'-'A'] = true;
		SNPtable[Q]['Y'-'A'] = true;
		SNPtable[Q]['H'-'A'] = true;
		SNPtable[Q]['N'-'A'] = true;
		SNPtable[Q]['K'-'A'] = true;
		SNPtable[Q]['D'-'A'] = true;
		SNPtable[Q]['E'-'A'] = true;

		SNPtable[N] = new boolean[26];
		SNPtable[N]['I'-'A'] = true;
		SNPtable[N]['T'-'A'] = true;
		SNPtable[N]['S'-'A'] = true;
		SNPtable[N]['Y'-'A'] = true;
		SNPtable[N]['H'-'A'] = true;
		SNPtable[N]['Q'-'A'] = true;
		SNPtable[N]['K'-'A'] = true;
		SNPtable[N]['D'-'A'] = true;
		SNPtable[N]['E'-'A'] = true;

		SNPtable[K] = new boolean[26];
		SNPtable[K]['I'-'A'] = true;
		SNPtable[K]['M'-'A'] = true;
		SNPtable[K]['T'-'A'] = true;
		SNPtable[K]['R'-'A'] = true;
		SNPtable[K]['Y'-'A'] = true;
		SNPtable[K]['H'-'A'] = true;
		SNPtable[K]['Q'-'A'] = true;
		SNPtable[K]['N'-'A'] = true;
		SNPtable[K]['D'-'A'] = true;
		SNPtable[K]['E'-'A'] = true;

		SNPtable[D] = new boolean[26];
		SNPtable[D]['V'-'A'] = true;
		SNPtable[D]['A'-'A'] = true;
		SNPtable[D]['G'-'A'] = true;
		SNPtable[D]['Y'-'A'] = true;
		SNPtable[D]['H'-'A'] = true;
		SNPtable[D]['Q'-'A'] = true;
		SNPtable[D]['N'-'A'] = true;
		SNPtable[D]['K'-'A'] = true;
		SNPtable[D]['E'-'A'] = true;

		SNPtable[E] = new boolean[26]; 
		SNPtable[E]['V'-'A'] = true;
		SNPtable[E]['A'-'A'] = true;
		SNPtable[E]['G'-'A'] = true;
		SNPtable[E]['Y'-'A'] = true;
		SNPtable[E]['H'-'A'] = true;
		SNPtable[E]['Q'-'A'] = true;
		SNPtable[E]['N'-'A'] = true;
		SNPtable[E]['K'-'A'] = true;
		SNPtable[E]['D'-'A'] = true;

		SNPtable[C] = new boolean[26];
		SNPtable[C]['F'-'A'] = true;
		SNPtable[C]['S'-'A'] = true;
		SNPtable[C]['Y'-'A'] = true;
		SNPtable[C]['W'-'A'] = true;
		SNPtable[C]['R'-'A'] = true;
		SNPtable[C]['S'-'A'] = true;
		SNPtable[C]['G'-'A'] = true;
		
		SNPtable[W] = new boolean[26]; 
		SNPtable[W]['L'-'A'] = true;
		SNPtable[W]['S'-'A'] = true;
		SNPtable[W]['C'-'A'] = true;
		SNPtable[W]['R'-'A'] = true;
		SNPtable[W]['S'-'A'] = true;
		SNPtable[W]['G'-'A'] = true;

		SNPtable[R] = new boolean[26]; 
		SNPtable[R]['L'-'A'] = true;
		SNPtable[R]['P'-'A'] = true;
		SNPtable[R]['H'-'A'] = true;
		SNPtable[R]['Q'-'A'] = true;
		SNPtable[R]['C'-'A'] = true;
		SNPtable[R]['W'-'A'] = true;
		SNPtable[R]['S'-'A'] = true;
		SNPtable[R]['G'-'A'] = true;
		SNPtable[R]['I'-'A'] = true;
		SNPtable[R]['M'-'A'] = true;
		SNPtable[R]['T'-'A'] = true;
		SNPtable[R]['K'-'A'] = true;

		SNPtable[G] = new boolean[26]; 
		SNPtable[G]['V'-'A'] = true;
		SNPtable[G]['A'-'A'] = true;
		SNPtable[G]['D'-'A'] = true;
		SNPtable[G]['E'-'A'] = true;
		SNPtable[G]['C'-'A'] = true;
		SNPtable[G]['W'-'A'] = true;
		SNPtable[G]['R'-'A'] = true;
		SNPtable[G]['S'-'A'] = true;
	}
}

class calculateED {

	private String sequence = null;
	private String query = null;

	public calculateED(String sequence, String query) {
		this.sequence = sequence;
		this.query = query;
	}

	public int calculate() {

		if(Constants.IS_I_SAME_WITH_L){
			sequence = sequence.replaceAll("I", "L");
			query = query.replaceAll("I", "L");
		}
		
		int count = 0;
		char[] seq = sequence.toCharArray();
		char[] que = query.toCharArray();

		char sequenceChar = 0;
		char queryChar = 0;


		for (int i = 0; i < seq.length - que.length + 1; i++) {
			count = 0;
			sequenceChar = seq[i];
			for (int j = 0; j < que.length; j++) {
				if (seq[i + j] != que[j]){
					count++;
					sequenceChar = seq[i + j];
					queryChar = que[j];
				}

				if (count > 1)
					break;
			}
			if (count == 1) {
				int check = checkSNP(sequenceChar, queryChar);
				if(check== 1){
					return 1;
				}
			}
		}

		return 0;
	}
	
	private int checkSNP(char sequenceChar, char queryChar) {
		if(ProteinMappingWithSAVs.SNPtable[sequenceChar-'A'][queryChar-'A']){
			return 1;
		}
		return 0;
		
	}

}
