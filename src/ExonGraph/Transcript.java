package ExonGraph;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

class ExonRangeType extends ExonRange implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -2521571711161357871L;
	// start, end, mut에 대한 정보는 ExonRange에 있음.
	// 본 클래스에서는 mut 대신  isCDS를 사용함.
	public boolean isCDS = false;
	
}

public class Transcript implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3255582802045066217L;
	public ArrayList<ExonRangeType> exonList;
	public String transcriptID;
	private boolean isPseudo;
	
	// 생성자
	public Transcript() {
		exonList = new ArrayList<ExonRangeType>();
		transcriptID = null;
		isPseudo = true;
	}
	
	public boolean isPseudo(){
		return isPseudo;
	}
	
	public void reverseExon(){
		
		Collections.sort(exonList, new TranscriptComparator());
		
	}
	
	public void addExon(String[] GTFLineSplit, boolean strand){
		String[] attr = GTFLineSplit[8].split(";");

		ExonRangeType curERT = new ExonRangeType();
		curERT.start =  Integer.parseInt(GTFLineSplit[3]);
		curERT.end =  Integer.parseInt(GTFLineSplit[4]);
		String type = GTFLineSplit[2]; 
		
		// CDS-exon, exon-CDS 둘 중 어떤 순서도 나올 수 있음.
		// 따라서 두 순서에 대한 처리를 고려함.
		
		if(type.equalsIgnoreCase("CDS")){
			curERT.isCDS = true;	isPseudo = false;
		}
		
		// exon size가 0이면 비교대상이 없으므로, 바로 추가함.
		if(exonList.size() == 0){
			transcriptID= ExonUtils.getGtfAttr(attr,"transcript_id");
			exonList.add(curERT);
			return;
		}
		
		int ERTIndex = exonList.size() - 1;
		ExonRangeType prevERT = exonList.get(ERTIndex);

		// 둘다 exon이면 바로 추가하고 종료.
		if(!(curERT.isCDS || prevERT.isCDS)){
			exonList.add(curERT);
			return;
		}
		
		// prevERT = exon 이고 curERT = CDS 가 되게끔 유지.
		if(prevERT.isCDS){
			prevERT = curERT;
			curERT = exonList.get(ERTIndex);
		}
		
		// 만약, exon과 CDS의 범위가 겹치지 않으면 한 세트의 exon-CDS가 아니므로 넘어감.
		if(!(prevERT.start <= curERT.start && prevERT.end >= curERT.end)){
			if(exonList.get(ERTIndex).isCDS){
				exonList.add(prevERT);
			}else{
				exonList.add(curERT);
			}
			
			return;
		}
		
		// exon과 CDS가 겹치는 범위를 고려.
		// case1: exon과 CDS가 완전히 겹치는 경우.
		//	exon: -----------
		//	CDS : -----------
		if(prevERT.start == curERT.start && prevERT.end == curERT.end){
			exonList.get(ERTIndex).isCDS = true;
			return;
		}
		
		// 마지막 엑손 범위 정보를 제거함.
		exonList.remove(ERTIndex);
		
		// case2: exon이 CDS보다 먼저 시작할 경우
		//	exon: ----------
		//	CDS :   --------
		ExonRangeType forwardERT = null;
		if(prevERT.start < curERT.start){
			forwardERT = new ExonRangeType();
			forwardERT.start = prevERT.start;
			forwardERT.end = curERT.start - 1;
		}
		
		// case3: exon이 CDS보다 늦게 끝나는 경우
		//	exon: ----------
		//	CDS : --------
		ExonRangeType backwardERT = null;
		if(prevERT.end > curERT.end){
			backwardERT = new ExonRangeType();
			backwardERT.start = curERT.end + 1;
			backwardERT.end = prevERT.end;
		}
		
		if(strand){
			if(forwardERT != null){
				exonList.add(forwardERT);
			}exonList.add(curERT);
			if(backwardERT != null){
				exonList.add(backwardERT);
			}
		}else{
			if(backwardERT != null){
				exonList.add(backwardERT);
			}exonList.add(curERT);
			if(forwardERT != null){
				exonList.add(forwardERT);
			}
		}
		
	}
	
}

class TranscriptComparator implements Comparator<ExonRangeType>{
	 
    @Override
    public int compare(ExonRangeType ERT1, ExonRangeType ERT2) {
    	if(ERT1.start == ERT2.start){
    		return 0;
    	}
    	
        return ERT1.start > ERT2.start ? 1 : -1;
    }
} 