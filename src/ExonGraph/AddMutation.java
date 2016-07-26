package ExonGraph;

/*
 * AddMutation | ver 1.0.0
 * Writer: progistar
 * 
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import Environments.Constants;

class MutationStructure{
	//mutation numbering
	//mutation == -1	->	Deletion
	//mutation == 1	->	Insertion
	//mutation == 0	->	SNP
	int pos, mutation;
	
	// origin : 기존 문자열
	// substitution : 뮤테이션 문자열
	String origin, substitution;
	
}

class LargeScaleDeletion extends MutationStructure {
	EXON startE = null;
	
	// Large Deletion의 실제 시작과 끝
	int endPos;
	
}

// Exon이 걸쳐있는 범위 정보를 갖고 있는 클래스
class ExonRange implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5861622815037956658L;
	// start : exon 내에서 시작 위치
	// end : exon 내에서 끝 위치
	// mut : true == 뮤테이션, false == 노말
	int start;
	int end;
	boolean isMut;
}

// Exon의 연결 상태를 갖고 있는 클래스
class ExonStructure {
	//
	EXON exon;
	int trans_cnt;
	
	// in : InEdge
	// out : OutEdge
	int in;
	int out;
}

public class AddMutation {
	
	// RefSeq Homo_sapiens.GRCh37.71.gtf & Ensembl GTF 기준 
	// 엑손 평균 길이는 252개
	private int BIN_SIZE = 252;
	public static final int MUT_TYPE_NUM = 3;
	public static final int MEDIAN = 7;
	public static final int EACH_MUT_MAX_SIZE = 3;
	private Hashtable<Integer, ArrayList<MutationStructure>> MSHash = new Hashtable<Integer, ArrayList<MutationStructure>>();
	private Hashtable<Integer, ArrayList<LargeScaleDeletion>> LargeScaleDeletionHash = new Hashtable<Integer, ArrayList<LargeScaleDeletion>>();
	private int curTransCount = 0;
	
	// LargeScaleDeletion 처리관련 변수
	
	// Path way 
	private String VCFs[] = null;
			
	private boolean strand;
	private boolean visit;
	private String CHR_IDENTIFIER = null;
	
	private int trans_cnt = 0;
	
	/**
	 * VCFs folder를 입력
	 * 
	 * @param CHR_IDENTIFIER
	 * @param VCFPATH
	 * @throws IOException
	 */
	public AddMutation( String CHR_IDENTIFIER, String VCFPATH) throws IOException{
		int length = 0;
						
		File[] fileList = new File(VCFPATH).listFiles();
		if(fileList == null){
			fileList = new File[1];
			fileList[0] = new File(VCFPATH);
		}
		
		length = fileList.length;
		this.VCFs = new String[length];
		
		for(int i=0; i<length; i++){
			this.VCFs[i] = fileList[i].getPath();
		}
		
		this.CHR_IDENTIFIER = CHR_IDENTIFIER;
		
		// 해당 염색체의 뮤테이션 구조를 만듦
		
		ConstructMutationStructure();
	}
	
	
	/**
	 * VCF 파일의 뮤테이션 정보를 읽고 자료구조를 만듦
	 * 
	 */
	private void ConstructMutationStructure() throws IOException{
		
		for(String s : VCFs){
			
			
			File f =new File(s);
			FileReader FR = new FileReader(f);
			BufferedReader BR = new BufferedReader(FR);
			String line;
			
			// VCF 파일을 읽고 구조를 만드는 루틴
			while((line = BR.readLine()) != null){
				String[] lineSplit = null;
				
				// VCF 레코드는 항목이 10개 이상임
				lineSplit = line.split("\t");
				if(lineSplit.length < 10){
					Logger.getLogger("AddMutation").log(Level.WARNING, "There is a missing value caused by wrong VCF format.");
					continue; 
				}
				
				
				// 두 염색체 번호가 다르면 넘어감
				if(!lineSplit[0].equalsIgnoreCase(CHR_IDENTIFIER)){
					continue;
				}
				
				//Consists mutation structure MS.
				ArrayList<MutationStructure> MSList = new ArrayList<MutationStructure>();
				
				MutationStructure MS = new MutationStructure();
				MSList.add(MS);
				
				int pos = MS.pos = Integer.parseInt(lineSplit[1]);
				String origin = MS.origin = lineSplit[3];
				String substitution = MS.substitution = lineSplit[4];
				
				if(substitution.contains(",")){
					String[] substitutions = substitution.split(",");
					
					MS.substitution  = substitutions[0];
					
					for(int i=1; i<substitutions.length; i++){
						MutationStructure ms = new MutationStructure();
						ms.substitution = substitutions[i]; 
						MSList.add(ms);
					}
				}
				
				
				// Mutation 종류를 결정
				
				//SNP
				if(MS.origin.length() == MS.substitution.length()){
					MS.mutation = 0;
				}
				//Deletion
				else if(MS.origin.length() > MS.substitution.length()){
					MS.mutation = -1;
				}
				//Insertion
				else{
					MS.mutation = 1;
				}
				
				int mutation = MS.mutation;
				int BIN_NUM = MS.pos/BIN_SIZE;
				if(MSHash.get(BIN_NUM) == null){
					MSHash.put(BIN_NUM, new ArrayList<MutationStructure>());
				}
				
				for(int i=0; i<MSList.size(); i++){
					MutationStructure ms = MSList.get(i);
					ms.origin = origin;
					ms.pos = pos;
					ms.mutation = mutation;
					
					MSHash.get(BIN_NUM).add(MSList.get(i));
				}
				
				
			}
			
			BR.close();
			FR.close();
		}
		
		// 정렬
		if(MSHash.size() == 0) return;
		Iterator<Integer> keys = (Iterator)MSHash.keys();
		
		
		while(keys.hasNext()){
			ArrayList<MutationStructure> MSList = MSHash.get(keys.next());
			
			Collections.sort(MSList, new MSComparator());
			
			// 중복 제거
			int pivot = 0;
			for(int i=1; i<MSList.size(); i++){
				if(MSList.get(i).pos == MSList.get(pivot).pos &&
					MSList.get(i).mutation == MSList.get(pivot).mutation &&
					MSList.get(i).substitution.equals(MSList.get(pivot).substitution)){
					MSList.remove(i);
					i -= 1;
				}else{
					pivot = i;
				}
			}
		}

	}
	
	public void addToGene(GENE gene){
		
		
		if(MSHash.size() == 0){
			return;
		}
		
		if(gene.get_trans().get(0).get_check())
			visit = false;
		else
			visit = true;
		
		int transSize = gene.get_trans_cnt();
		
		for(int i=0; i<transSize; i++){
			strand = gene.get_strand();
			trans_cnt = gene.get_trans_cnt();
			curTransCount = i;
			
			// Large Scale Deletion이 발견되면 isLargeSacleDeletion = true 임
			addToExon(gene.get_trans().get(i));
			
		}
	}
	
	private void addToExon(EXON exon){
		if(visit == exon.get_check()){
			return;
		}
		exon.set_check(visit);
		
		// Transcript의 끝
		if(exon.get_end() == Integer.MAX_VALUE-1){
			return;
		}
		
		//extract Mutation group ranged from exon-start to exon-end
		ArrayList<MutationStructure> msStructure = extractGroup(exon.get_start(), exon.get_end());
		ArrayList<ExonStructure> prevExon = ExonUtils.getPrevExon(exon, trans_cnt);
		ArrayList<ExonStructure> nextExon = ExonUtils.getNextExon(exon, trans_cnt);
		
		if(msStructure.size() != 0){
			//extract normal-exon and mutation-exon range from exon-start to exon-end
			//It uses msStructrue.
			ArrayList<ExonRange> ER = extractRange(exon.get_start(), exon.get_end(), msStructure);
			
			//makeMid Exonz using ER && Connect them each other
			EXON startExon = makeExonz(exon, ER, prevExon, nextExon );
			
			//add mutation to midz
			startExon = addMut(exon, startExon, msStructure);
			
			
		}

		// exon이 바뀌지 않음
		//ArrayList<ExonStructure> nexts = getNextExon(exon);
		
		for(ExonStructure es : nextExon){
			addToExon(es.exon);
		}
	}
	
	
	private ArrayList<MutationStructure> extractGroup(int start, int end){
		ArrayList<MutationStructure> MS = new ArrayList<MutationStructure>();
		int binStart = start/BIN_SIZE;
		int binEnd = end/BIN_SIZE;
		
		for(int binIndex=binStart; binIndex<=binEnd; binIndex++){
			ArrayList<MutationStructure> MSList = MSHash.get(binIndex);
			if(MSList != null){
				for(MutationStructure ms : MSList){
					if(ms.pos >= start && ms.pos <= end){
						
						if(ms.pos + ms.origin.length() - 1 > end){
							// Large Scale Deletion일 경우 제외
							// HandleLargeScaleDeletion 메소드에서 따로 처리함
							// TODO: LargeScaleDeletion
							System.out.println("TODO: LargeScaleDeletion");
							/*System.out.println(ms.pos);*/
						}else{
							MS.add(ms);
						}
						
					}else if(ms.pos > end){
						break;
					}
				}
			}
			
		}
		
		return MS;
	}

	private ArrayList<ExonRange> extractRange(int start, int end, ArrayList<MutationStructure> msStructure){
		ArrayList<ExonRange> exonRange = new ArrayList<ExonRange>();
		int norStart = start;
		int mutStart = msStructure.get(0).pos, mutEnd = msStructure.get(0).pos + msStructure.get(0).origin.length() - 1;
		int mutSize = msStructure.size();
		
		
		for(int i= 1; i<mutSize; i++){
			MutationStructure ms = msStructure.get(i);
			
			int msStart = ms.pos;
			int msEnd = ms.pos + ms.origin.length() - 1;
			
			if(mutStart <= msStart && mutEnd >= msStart){
				mutEnd = mutEnd > msEnd ? mutEnd : msEnd;
			}else{
				NorOrMut:
					while(true){
						
						// 뮤테이션일 경우
						if(norStart == mutStart){
							
							exonRange.add(getExonRange(mutStart, mutEnd, true));
							
							norStart = mutEnd + 1;
							mutStart = msStart;
							mutEnd = msEnd;
							
							break NorOrMut;
						}
						// 노말일 경우
						else{
							
							exonRange.add(getExonRange(norStart, mutStart-1, false));
							// 이 경우 반드시 norStart == mutStart를 만족함
							// 따라서 다음 루프에서는 상위 if문에 해당되어 루프가 종료됨							
							norStart = mutStart;
 
						}
					}
			}
			
		} // End For Loop
		
		if(norStart == mutStart){
			exonRange.add(getExonRange(mutStart, mutEnd, true));
		}else if(norStart < mutStart){
			exonRange.add(getExonRange(norStart, mutStart-1, false));
			exonRange.add(getExonRange(mutStart, mutEnd, true));
		}
		
		if(mutEnd != end){
			exonRange.add(getExonRange(mutEnd+1, end, false));
		}
		
		
		/* Debugging Code
		 * It prints reads*/
		/*for(ExonRange er : exonRange){
			System.out.println(er.start + "-" +er.end +"-" +er.isMut);
		}
		System.out.println("--");*/
		
		return exonRange;
	}
	
	private ExonRange getExonRange(int start, int end, boolean mutation){
		ExonRange returnRange = new ExonRange();
		
		returnRange.end = end;
		returnRange.start = start;
		returnRange.isMut = mutation;
		
		return returnRange;
	}
	
	private EXON makeExonz(EXON exon, ArrayList<ExonRange> ER, ArrayList<ExonStructure> prevES, ArrayList<ExonStructure> nextES){
		EXON startExon = new EXON("START", null, null, exon.get_start(), exon.get_end(), -1, trans_cnt);
		ArrayList<EXON> midzList = new ArrayList<EXON>();
		EXON newExon = null;
		
		int exonStart = exon.get_start();
		int exonNumber = exon.get_number();
		String exonSeq = exon.get_seq();
		int exonSeqLen = exonSeq.length();
		
		//makeExonz
		for(ExonRange er : ER){
			
			int start = er.start;
			int end = er.end;
			int textPos = start - exonStart;
			
			//mutation
			if(er.isMut){
				
				while(exonStart + textPos <= end){
					
					newExon = new EXON(String.valueOf(exonSeq.charAt(textPos)), null, null, exonStart+textPos, exonStart+textPos, exonNumber, trans_cnt);
					newExon.set_check(visit);
					
					midzList.add(newExon);
					
					textPos ++;
				}
				
			}
			//normal
			else{
				newExon = new EXON(exonSeq.substring(start - exonStart, end - exonStart + 1), null, null, start, end, exonNumber, trans_cnt);
				newExon.set_check(visit);
				midzList.add(newExon);
				
			}
			
		}
		
		startExon.set_next_junc(0, 0, midzList.get(0));
		
		for(ExonStructure es : nextES){
			midzList.get(midzList.size()-1).set_next_junc(es.trans_cnt, es.in, es.exon);
			es.exon.set_prev_junc(es.trans_cnt, es.out, midzList.get(midzList.size()-1));
		}
		
		for(ExonStructure es : prevES){
			midzList.get(0).set_prev_junc(es.trans_cnt, es.in, es.exon);
			es.exon.set_next_junc(es.trans_cnt , es.out, midzList.get(0));
		}
		
		//connect !
		for(int i=1; i<midzList.size(); i++){
			
			EXON left = midzList.get(i-1);
			EXON right = midzList.get(i);
			
			// 뮤테이션 연결은 transcript 0만 이용해서 연결함.
			left.set_next_junc(0, 0, right);
			right.set_prev_junc(0, 0, left);
			
		}
		
		
		
		/* Debugging Code
		 * It prints connecting exons
		 * 
		 */ 
		/*newExon = startExon.get_next(0, 0);
		
		while(true){
			System.out.println(newExon.get_start()+"-"+newExon.get_end()+":"+newExon.get_seq());
			
			if(newExon.get_end() == exon.get_end()){
				break;
			}
			
			for(ExonStructure es : prevES){
				newExon = newExon.get_next(es.trans_cnt, 0);
				break;
			}
		}
		
		System.out.println("-------------");*/
		
		return startExon;
	}
	
	private EXON addMut(EXON exon, EXON startExon, ArrayList<MutationStructure> msStructures){
		
				
		int[] mutCount = new int[MUT_TYPE_NUM];
		
			
		int prevMsPos = -1;
		EXON tempExon = startExon.get_next(0, 0);
		for(MutationStructure ms : msStructures){
			if(prevMsPos != ms.pos){
				prevMsPos = ms.pos;
				
				// 뮤테이션 카운트 초기화
				for(int i=0; i<MUT_TYPE_NUM; i++){
					mutCount[i] = 0;
				}
			}
			int mutType = ms.mutation + 1;
			
			// 동일 장소에서 종류별로 각 3개를 초과하여 적용될 수 없음
			// ex>
			// snp, snp, ins, ins 가능
			// snp, snp, snp, snp 불가능 (snp가 3개를 초과)
			// snp, snp, snp, ins, ins, ins, del, del, del 가능
			// del, del, del, ins, del 불가능 (del이 3개를 초과)
			if(mutCount[mutType] == 3){
				continue;
			}
			
			//locate to the ms position

			
			while(!(tempExon.get_start() == ms.pos)){
				tempExon = tempExon.get_next(0, 0);
			}
			
			ArrayList<ExonStructure> tempPrevExon = ExonUtils.getPrevExon(tempExon, trans_cnt);
			ArrayList<ExonStructure> tempNextExon = ExonUtils.getNextExon(tempExon, trans_cnt);
			
			
			//deletion
			if(ms.mutation == -1){
				EXON endE = tempExon;
				EXON startE = tempExon;
				int msEndPos = ms.pos + ms.origin.length();
				int startPos = ms.pos + ms.substitution.length();
				
				//startPos
				while(!(startE.get_start() == startPos)){
					startE = startE.get_next(0, 0);
				}startE = startE.get_prev(0, 0);
				//endPos
				endE = startE;
				
				while(!(endE.get_start() >= msEndPos-1)){
					endE = endE.get_next(0, 0);
				}
				
				tempNextExon = ExonUtils.getNextExon(endE, trans_cnt);
				
				for(ExonStructure es : tempNextExon){
					es.exon.set_prev_junc(es.trans_cnt, MEDIAN + EACH_MUT_MAX_SIZE*ms.mutation + mutCount[mutType], startE);
					startE.set_next_junc(es.trans_cnt, MEDIAN + EACH_MUT_MAX_SIZE*ms.mutation + mutCount[mutType], es.exon);
				}
				
				
				mutCount[mutType] ++;
			}
			
			//insertion OR SNP
			else{
				String substitution;
				
				if(!strand){
					substitution = ExonUtils.complement(ms.substitution);
				}else{
					substitution = ms.substitution;
				}
				
				EXON mutExon = new EXON(substitution, null, null, ms.pos, ms.pos, exon.get_number(), trans_cnt);
				
				for(ExonStructure es : tempPrevExon){
					es.exon.set_next_junc(es.trans_cnt, MEDIAN + EACH_MUT_MAX_SIZE*ms.mutation + mutCount[mutType], mutExon);
					mutExon.set_prev_junc(es.trans_cnt, es.in, es.exon);
				}
				
				for(ExonStructure es: tempNextExon){
					es.exon.set_prev_junc(es.trans_cnt, MEDIAN + EACH_MUT_MAX_SIZE*ms.mutation + mutCount[mutType], mutExon);
					mutExon.set_next_junc(es.trans_cnt, es.out, es.exon);
				}
				
				
				mutCount[mutType] ++;
			}
			
		}
		
		
		return startExon;
	}

	

	public void setBinSize(int BIN_SIZE){
		this.BIN_SIZE = BIN_SIZE;
	}
	
	public int getBinSize(){
		return this.BIN_SIZE;
	}
}


class MSComparator implements Comparator<MutationStructure>{
	 
    @Override
    public int compare(MutationStructure ms1, MutationStructure ms2) {
        return ms1.pos > ms2.pos ? 1 : -1;
    }
} 