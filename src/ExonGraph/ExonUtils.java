package ExonGraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import ExonGraph.GENE.StartPos;

class MutationStructure_ extends MutationStructure{
	int endPos;
}

public class ExonUtils {

	public static final int DEL_START_INDEX = AddMutation.MEDIAN - AddMutation.MUT_TYPE_NUM;
	public static final int SNP_START_INDEX = AddMutation.MEDIAN;
	public static final int INS_START_INDEX = AddMutation.MEDIAN + AddMutation.MUT_TYPE_NUM;
	
	private static ArrayList<Integer> jncStartList = null;
	private static ArrayList<Integer> jncEndList = null;
	
	/**
	 * A, C, G, T로 이루어진 String, S를 받으면 생물학적 보수연산을 하여 반환한다.
	 * A, C, G, T를 제외한 문자는 보수연산이 존재하지 않으므로 기존 문자를 반환한다.
	 * 
	 * @param fastaSeq 문자열
	 * @return 기존 S에 보수연산을 한 String
	 * @see {@link https://en.wikipedia.org/wiki/FASTA_format}
	 * 
	 */
	public static String complement(String fastaSeq) {
		StringBuilder complementedText = new StringBuilder();
		Hashtable<Character, Character> complementHashTable = new Hashtable<>();
		int length = fastaSeq.length();
		
		complementHashTable.put('A', 'T');
		complementHashTable.put('T', 'A');
		complementHashTable.put('C', 'G');
		complementHashTable.put('G', 'C');
		
		
		for(int i=0; i<length; i++){
			// A, C, T, G 이 외의 코드가 들어갈 가능성이 있음
			// EX> N, K, ... 
			// https://en.wikipedia.org/wiki/FASTA_format 에서 Sequence representation 항목을 참조
			if(complementHashTable.get(fastaSeq.charAt(i)) != null){
				complementedText.append(complementHashTable.get(fastaSeq.charAt(i)));
			}else{
				complementedText.append(fastaSeq.charAt(i));
			}
			
		}
		
		return complementedText.toString();
	}
	

	/**
	 * getGtfAttr
	 * GTF에서 attribute 항목을 Parsing
	 * 
	 */
	public static String getGtfAttr(String[] attr, String tag){
		
		for(String _s : attr){
			if(_s.contains(tag)){
				return _s.replaceAll("[\"\\s]|"+tag, "");
			}
		}
		
		return null;
	}
	
	/**
	 * exon에서 발생한 이벤트를 체크한다.
	 * 
	 * 
	 */

	public static void variationInit(){
		jncStartList = new ArrayList<Integer>();
		jncEndList = new ArrayList<Integer>();
	}
	
	public static ArrayList<MutationStructure_> getVariation(StartPos startPos, int exonIndex, int trans_cnt, boolean strand){
		EXON exon = startPos.startExon.get(exonIndex);
		ArrayList<MutationStructure_> MSList = new ArrayList<MutationStructure_>();
		ArrayList<ExonStructure> prevExonz = getPrevExon(exon, trans_cnt);
		EXON prevExon = null;
		HashMap<Integer, Boolean> checkEvent = new HashMap<Integer, Boolean>();
		
		for(ExonStructure es : prevExonz){
			prevExon = es.exon;
			
			if(checkEvent.get(es.out) == null){
				checkEvent.put(es.out, true);
			}else{
				continue;
			}
			
			MutationStructure_ MS = new MutationStructure_();
			
			String var = null;
			
			switch(es.out){
			// DEL
			// : mark as -1
			case DEL_START_INDEX: case DEL_START_INDEX+1:  case DEL_START_INDEX+2: 	
				if(var == null){ var = "del"; MS.mutation = -1;}
			// SNP
			// : mark as 0
			case SNP_START_INDEX: case SNP_START_INDEX+1: case SNP_START_INDEX+2:
				if(var == null){ var = "snp"; MS.mutation = 0; }
			// INS
			// : mark as 1
			case INS_START_INDEX: case INS_START_INDEX+1: case INS_START_INDEX+2:
				if(var == null){ var = "ins"; MS.mutation = 1; }
			MS.pos = exon.get_start();
			MS.substitution = exon.get_seq();
			
			// 어디서 왔는지, 기존이 무엇인지만 알면 됨
			int iList[] = {0, 1, 2, 3, ExonGraph.JVALUE-1};
			for(int i : iList){
				if(prevExon.get_next(es.trans_cnt,i) != null && prevExon.get_next(es.trans_cnt,i).get_start() == exon.get_start()){
					MS.origin = prevExon.get_next(es.trans_cnt, i).get_seq();
					break;
				}
			}
			
			
			MSList.add(MS);
			
			break;
			
			default: break;
			}
			
		}

		// JNC : mark as -2
		// ALT : mark as +2
		
		if(exonIndex != 0){
			prevExon = startPos.startExon.get(exonIndex-1);
			
			// Alternative Splicing
			for(int i=0; i<trans_cnt; i++){
				if(prevExon.get_next(i, ExonGraph.JVALUE - 1) == exon){
					MutationStructure_ MS = new MutationStructure_();
					MS.pos = prevExon.get_end();
					MS.substitution = exon.get_seq();
					MS.endPos = exon.get_start();
					MS.mutation = 2;
					MSList.add(MS);
					break;
				}
			} 
		
			
			// Junction Variation
			// 1, 2, 3 이 junction variation과 연관된 에지 번호임
			int jncList[] = {1, 2, 3};
			int decisionMask = 0b0000_0;
			boolean isJnc = false;
			
			
			// Decision Tree Search
			for(int ti=0; ti<trans_cnt; ti++){
				for(int jnc : jncList){
					if(prevExon.get_next(ti, jnc) == exon){
						isJnc = true;
						if(jnc == 1){
							decisionMask |= 0b0000_1;
						}
					}
				}
				
				if(prevExon.get_prev(ti, 0) != null){
					decisionMask |= 0b0001_0;
				}
				
				// A` Condition
				if(prevExon.get_next(ti, 0) != null){
					boolean isAPrime = false;

					for(int ti_=0; ti_<trans_cnt; ti_++){
						if(prevExon.get_prev(ti_, 1) != null){
							for(int ti__=0; ti__<trans_cnt; ti__++){
								if(prevExon.get_prev(ti_, 1).get_next(ti__, 0) == exon){
									isAPrime = true;
								}
							}	
						}
						
						if(prevExon.get_prev(ti_, 0) != null){
							for(int ti__=0; ti__<trans_cnt; ti__++){
								if(prevExon.get_prev(ti_, 0).get_next(ti__, 0) == exon){
									isAPrime = true;
								}
							}	
						}
					}
					
					for(int ti_=0; ti_<trans_cnt; ti_++){
						if(prevExon.get_next(ti_, 3) == exon){
							isAPrime = false;
						}
					}
					
					if(!isAPrime){
						decisionMask |= 0b0010_0;
					}
				}
				
				if(exon.get_prev(ti, 0) != null){
					decisionMask |= 0b0100_0;
				}
				
				if(prevExon.get_next(ti, 0) != null
						&&
					(prevExon.get_next(ti, 0).get_start() == prevExon.get_end()+1)){
					
					EXON intraExon = prevExon.get_next(ti, 0);
					for(int ti_=0; ti_<trans_cnt; ti_++){
						if(intraExon.get_next(ti_, 0) == exon){
							decisionMask |= 0b1000_0;
						}
					}
				}
			}
			
			
			if(isJnc){
				//MASK CONFIRMATION START
				/*System.out.println("MASK: "+decisionMask);
				if((decisionMask & 0b0001_1) == 0b0000_1){
					// case: B`
					System.out.println("B`");
				}else if((decisionMask & 0b0001_1) == 0b0001_1){
					// case: A
					System.out.println("A");
				}else if((decisionMask & 0b0010_1) == 0b0000_0){
					// case: A`
					System.out.println("A`");
				}else if((decisionMask & 0b0110_1) == 0b0010_0){
					// case: B
					System.out.println("B");					
				}else if((decisionMask & 0b1110_1) == 0b0110_0){
					// case: C
					System.out.println("C");
				}else if((decisionMask & 0b1110_1) == 0b1110_0){
					// case: D
					System.out.println("D");
				}*/
				// -MASK CONFIRMATION END
				
				MutationStructure_ MS = new MutationStructure_();
				
				if((decisionMask & 0b0001_1) == 0b0000_1){
					// case: B`
					MS.pos = prevExon.get_start();
					MS.endPos = exon.get_start();
					
					MS.substitution = prevExon.get_seq() + String.valueOf(exon.get_seq().charAt(0));
					MS.origin = getPadding(prevExon.get_seq().length())+String.valueOf(exon.get_seq().charAt(0));
					
				}else if((decisionMask & 0b0001_1) == 0b0001_1){
					// case: A
					MS.pos = prevExon.get_end();
					MS.endPos = exon.get_end();
					
					MS.origin = prevExon.get_seq();
					MS.substitution = String.valueOf(MS.origin.charAt(MS.origin.length()-1)) + exon.get_seq();
					MS.origin = String.valueOf(MS.origin.charAt(MS.origin.length()-1)) + getPadding(exon.get_seq().length());
					
				}else if((decisionMask & 0b0010_1) == 0b0000_0){
					// case: A`
					EXON originExon = null;
					for(int ti=0; ti<trans_cnt; ti++){
						if(prevExon.get_prev(ti, 1) != null){
							originExon = prevExon.get_prev(ti, 1);
							break;
						}
						
						if(prevExon.get_prev(ti, 0) != null){
							originExon = prevExon.get_prev(ti, 0);
							break;
						}
					}
					
					MS.pos = originExon.get_end();
					MS.endPos = prevExon.get_end();
					
					MS.origin = originExon.get_seq();
					MS.substitution = String.valueOf(MS.origin.charAt(MS.origin.length()-1)) + prevExon.get_seq();
					MS.origin = String.valueOf(MS.origin.charAt(MS.origin.length()-1)) + getPadding(prevExon.get_seq().length());
					
				}else if((decisionMask & 0b0110_1) == 0b0010_0){
					// case: B
					EXON originExon = null;
					for(int ti=0; ti<trans_cnt; ti++){
						if(exon.get_next(ti, 1) != null){
							originExon = exon.get_next(ti, 1);
							break;
						}
						
						if(exon.get_next(ti, 0) != null){
							originExon = exon.get_next(ti, 0);
							break;
						}
					}
					
					MS.pos = exon.get_start();
					MS.endPos = originExon.get_start();
					
					MS.substitution = exon.get_seq() + String.valueOf(originExon.get_seq().charAt(0));
					MS.origin = getPadding(exon.get_seq().length()) + String.valueOf(originExon.get_seq().charAt(0));
					
				}else if((decisionMask & 0b1110_1) == 0b0110_0){
					// case: C
					EXON originExon = null;
					for(int ti=0; ti<trans_cnt; ti++){
						if(exon.get_prev(ti, 0) != null){
							originExon = exon.get_prev(ti, 0);
							break;
						}
					}
					
					MS.pos = originExon.get_start();
					MS.endPos = exon.get_start();
					
					MS.origin = originExon.get_seq() + String.valueOf(exon.get_seq().charAt(0));
					MS.substitution = getPadding(originExon.get_seq().length()) + String.valueOf(exon.get_seq().charAt(0));
					
					
				}else if((decisionMask & 0b1110_1) == 0b1110_0){
					// case: D
					EXON originExon = null;
					for(int ti=0; ti<trans_cnt; ti++){
						if(prevExon.get_next(ti, 0) != null){
							originExon = prevExon.get_next(ti, 0);
							break;
						}
					}
					
					MS.pos = prevExon.get_end();
					MS.endPos = originExon.get_end();
					
					MS.origin = prevExon.get_seq();
					MS.origin = String.valueOf(MS.origin.charAt(MS.origin.length()-1)) + originExon.get_seq();
					MS.substitution = String.valueOf(prevExon.get_seq().charAt(prevExon.get_seq().length()-1)) + getPadding(originExon.get_seq().length());
					
					
				}				
				MS.mutation = -2;
				
				// jnc의 경우는 중복문제가 있음
				// 따라서 중복을 제거해야함
				boolean isDel = false;
				for(int i=0; i<jncStartList.size(); i++){
					if(MS.pos == jncStartList.get(i) && MS.endPos == jncEndList.get(i)){
						isDel = true;
					}
				}
				
				if(!isDel){
					jncStartList.add(MS.pos);
					jncEndList.add(MS.endPos);
					MSList.add(MS);
				}
			}
			
			
		}
		
		return MSList;
	}
	
	
	public static String getPadding(int len){
		
		StringBuilder SB = new StringBuilder();
		
		for(int i=0; i<len; i++){
			SB.append("-");
		}
		
		return SB.toString();
	}
	
	/**
	 * 
	 */
	
	public static ArrayList<ExonStructure> getPrevExon(EXON exon, int trans_cnt){
		ArrayList<ExonStructure> prevExon = new ArrayList<ExonStructure>();
		
		for(int i=0; i<trans_cnt; i++){
			for(int j=0; j<ExonGraph.JVALUE; j++){
				if(exon.get_prev(i,j) != null){
					ExonStructure es = new ExonStructure();
					es.exon = exon.get_prev(i,j);
					es.in = j;
					es.trans_cnt = i;
					
					prevExon.add(es);
				}
			}
		}
		
		for(ExonStructure es : prevExon){
			for(int i=0; i<trans_cnt; i++){
				for(int j=0; j<ExonGraph.JVALUE; j++){
					if(es.exon.get_next(i,j) == exon){
						es.out = j;
					}
				}
			}
		}
		
		return prevExon;
	}
	
	/*
	 * nextExon
	 * 
	 * 
	 */
	public static ArrayList<ExonStructure> getNextExon(EXON exon, int trans_cnt){
		ArrayList<ExonStructure> nextExon = new ArrayList<ExonStructure>();
		
		for(int i=0; i<trans_cnt; i++){
			for(int j=0; j<ExonGraph.JVALUE; j++){
				if(exon.get_next(i,j) != null){
					ExonStructure es = new ExonStructure();
					es.exon = exon.get_next(i,j);
					es.in = j;
					es.trans_cnt = i;
					
					nextExon.add(es);
				}
			}
		}
		
		for(ExonStructure es : nextExon){
			for(int i=0; i<trans_cnt; i++){
				for(int j=0; j<ExonGraph.JVALUE; j++){
					if(es.exon.get_prev(i,j) == exon){
						es.out = j;
					}
				}
			}
		}
		
		return nextExon;
	}
}

