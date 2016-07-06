package ExonGraph;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.logging.Level;
import java.util.logging.Logger;

import Environments.Constants;
import ExonGraph.GENE.StartPos;

public class GFF {

	private static FileWriter GFFWriter = null;
	private static BufferedWriter GFFBW = null;
	
	public static ArrayList<String> temporaryGFF = null;
	/**
	 * NextSearch의 결과를 GFF 포맷으로 파일을 작성한다.
	 * @throws IOException 
	 * 
	 */
	public static void write(StartPos startPos, int textPos, boolean strand, int trans_cnt) throws IOException{
		LinkedList<EXON> connectedExonList = startPos.startExon;
		
		int sizeOfExon = connectedExonList.size();
		// 파일 오픈
		open(Constants.OUTPUT_VSG_GFF);
		temporaryGFF = new ArrayList<String>();
		//
		ExonUtils.variationInit();
		for(int i=0; i<sizeOfExon; i++){
			EXON curExon = connectedExonList.get(i);
			StringBuilder GFFSB = new StringBuilder();
			// seqName 부터 차레대로 기입
			// 0: seqName
			GFFSB.append(ExonGraph.Chrom)
			.append("\t");
			// 1: source
			GFFSB.append("ACTG")
			.append("\t");
			// 2: feature
			if(curExon.isExon){
				GFFSB.append("exon").append("\t");
			}else{
				GFFSB.append("intron").append("\t");
			}
			// 3: start
			if(i == 0)
				GFFSB.append(curExon.get_start() + startPos.startPos);
			else
				GFFSB.append(curExon.get_start());
			GFFSB.append("\t");
			// 4: end
			if(i == sizeOfExon-1)
				GFFSB.append(curExon.get_start() + textPos);
			else
				GFFSB.append(curExon.get_end());
			GFFSB.append("\t");
			// 5: score
			GFFSB.append(".")
			.append("\t");
			// 6: strand
			if(strand){
				GFFSB.append("+");
			}else{
				GFFSB.append("-");
			}
			GFFSB.append("\t");
			
			// 7: frame
			GFFSB.append(".")
			.append("\t");
			
			temporaryGFF.add(GFFSB.toString());
			
			// 8: ID in attributes
			//GFFSB.append("ID=").append(IDIndex);
			
			// variation에 대한 정보는 추가하지 않음
			// variation 정보를 추가할 필요가 있으면, 해당 주석을 제거하면 됨
			/*ArrayList<MutationStructure_> variation = ExonUtils.getVariation(startPos, i, trans_cnt, strand);
			
			if(variation.size() != 0)
				GFFSB.append(MS_toString(variation));*/
			
		}
				
	}
	
	/**
	 * Deprecated
	 * 
	 * @param MS_
	 * @return
	 */
	private static String MS_toString(ArrayList<MutationStructure_> MS_){
		StringBuilder str = new StringBuilder();
		
		str.append(";");
		
		int msSize= MS_.size();
		for(int i=0; i<msSize; i++){
			MutationStructure_ ms = MS_.get(i);
			// 구분자 & 
			if(i != 0){ str.append("&"); }
			
			// 뮤테이션 타입을 기입.
			switch(ms.mutation){
			case -2: str.append("jnc("); break; 
			case -1: str.append("del("); break;
			case 0: str.append("snv("); break;
			case 1: str.append("ins("); break;
			case 2: str.append("alt("); break;
			default: Logger.getLogger("GFF").log(Level.SEVERE, "You should check your code related to GFF Maker.");
			}
			
			// Chr 번호 기입
			str.append(ExonGraph.Chrom).append(":");
			
			// start Position 기입
			str.append(ms.pos).append(":");
			
			// 뮤테이션 타입마다 표시하는 게 다름
			if(ms.mutation >= -1 && ms.mutation <= 1){
				str.append(ms.origin).append("|").append(ms.substitution);
			}else{
				str.append(ms.endPos);
			}
			
			// 마감자
			str.append(")");
		}
		
		return str.toString();
	}
	
	private static void open(String fileName) throws IOException{
		if(GFFWriter == null || GFFBW == null){
			GFFWriter = new FileWriter(fileName);
			GFFBW = new BufferedWriter(GFFWriter);
		}
	}
	
	
	
	public static void close() throws IOException{
		int GFFIDIndex = 0;
		
		ArrayList<Result> resultList = new ArrayList<Result>(ResultMapping.ResultTable.values());
		
		for(Result result : resultList){
			GFFIDIndex++;
			for(String gff : result.GFFResult){
				GFFBW.append(gff+"ID="+GFFIDIndex+"\n");
			}
			
		}
		
		
		if(GFFWriter != null && GFFBW != null){
			GFFBW.close();
			GFFWriter.close();
		}
	}
}