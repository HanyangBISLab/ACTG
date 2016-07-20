package ExonGraph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.logging.Level;
import java.util.logging.Logger;

import Environments.Constants;
import ExonGraph.GENE.StartPos;

public class Flat {

	private static final int UTR_5 = 3;
	private static final int UTR_3 = 2;
	private static final int PSEUDO = -1;
	private static final int CDS = 4;
	private static final int FS = 1;
	private static final int UNKNOWN = -2;
	private static int theBestMappingType = 0;
	private static int theCurrMappingType = 0;
	
	private static File valOutput = null;
	private static FileWriter FWVal = null;
	private static BufferedWriter BWVal = null;
	
	private static File logOutput = null;
	private static FileWriter FWLog = null;
	private static BufferedWriter BWLog = null;
	private static String geneID_ = null;

	private static Result results = null;
	
	
	
	public static void write(StartPos startPos, int textPos, boolean strand, int trans_cnt, String valText, ArrayList<Transcript> transcriptList, String geneID) throws IOException {
		
		theBestMappingType = UNKNOWN;
		geneID_ = geneID;
		
		// 최초 1회 파일을 생성함
		if (valOutput == null) {
			String defaultValName = Constants.OUTPUT_VSG_FLAT;
			String defaultLogName = Constants.OUTPUT_VSG_LOG;
			valOutput = new File(defaultValName);
			valOutput.createNewFile();
			FWVal = new FileWriter(valOutput);
			BWVal = new BufferedWriter(FWVal);

			logOutput = new File(defaultLogName);
			logOutput.createNewFile();
			FWLog = new FileWriter(logOutput);
			BWLog = new BufferedWriter(FWLog);
			
			// Index를 추가함
			BWVal.append("GFFID\tPeptide\tGeneID\tAttribute\n");
			BWLog.append("GFFID\tGeneID\tTranscriptID\tAttribute\n");
		}

		
		int sizeOfExon = startPos.startExon.size();

		StringBuilder eventText = new StringBuilder();
		StringBuilder normalText = new StringBuilder();
		int count = 0;
		int jncGap = 0;
		boolean isJncOrAlt = false;
		// jnc 변이에 대한 중복제거를 위함
		
		
		// Event Extractor
		ExonUtils.variationInit();
		for (int i = 0; i < sizeOfExon; i++) {
			ArrayList<MutationStructure_> variation = ExonUtils.getVariation(startPos, i, trans_cnt, strand);
						
			if (variation.size() != 0) {
				count++;
				eventText.append(MS_toString(variation));
				for(int j=0; j<variation.size(); j++){
					// jnc의 경우 변형된 nucleotide 수 만큼 패널티를 줌
					MutationStructure_ var = variation.get(j);
					if(var.mutation == -2){
						
						for(int k=0; k<var.origin.length(); k++){
							if(var.origin.charAt(k) == '-'){
								jncGap++;
							}
						}
						
						isJncOrAlt = true;
					}
					
					if(var.mutation == 2){
						isJncOrAlt = true;
					}
				}
			}
		}
		
		if(eventText.length() != 0){
			eventText.deleteCharAt(eventText.lastIndexOf(";"));
		}

		results = new Result();
		normalText.append(valText).append("\t");
		// GENEID
		normalText.append(String.valueOf(geneID_)).append("\t");
		// CDSRegion
		if(!isJncOrAlt)
			normalText.append(extractCDSRegion(transcriptList, startPos, textPos, strand, jncGap, count));
		// Event Information
		normalText.append(eventText.toString());
		
		StringBuilder resultKey = new StringBuilder(normalText.toString());
		// chr기입
		results.chr = GFF.temporaryGFF.get(0).split("\t")[0];
		for(String gff : GFF.temporaryGFF){
			resultKey.append(gff);
		}
		if(ResultMapping.ResultTable.get(resultKey.toString()) == null){
			
			results.flatResult = normalText.toString();
			results.GFFResult = GFF.temporaryGFF;			
			ResultMapping.ResultTable.put(resultKey.toString(), results);
			
		}
		
		
	}

	//TODO: FS method
	// True => FS
	private static boolean isFS(Transcript transcript, int matchingStart, int matchingEnd, boolean strand){
		
		ArrayList<ExonRangeType> ERT = transcript.exonList;
		int sizeOfERT = ERT.size();
		int indexOfERT = 0;
		int frame = 0;
		int cumulativePos = 0;
		
		
		// find mapped CDS to the firstMappedExon
		if(strand){
			for( ; indexOfERT<sizeOfERT; indexOfERT++){
				if(!ERT.get(indexOfERT).isCDS) continue;
				
				ExonRangeType curERT = ERT.get(indexOfERT);
				
				if(curERT.start <= matchingStart && curERT.end >= matchingStart){
					frame = Math.abs(cumulativePos + curERT.start - matchingStart);
					break;
				}else{
					cumulativePos += (curERT.end-curERT.start+1);
				}
			}
		}else{
			indexOfERT = sizeOfERT - 1;
			for( ; indexOfERT > -1; indexOfERT--){
				if(!ERT.get(indexOfERT).isCDS) continue;
				
				ExonRangeType curERT = ERT.get(indexOfERT);
				
				if(curERT.start <= matchingEnd && curERT.end >= matchingEnd){
					frame = Math.abs(cumulativePos + curERT.end - matchingEnd);
					break;
				}else{
					cumulativePos += (curERT.end-curERT.start+1);
				}
			}
		}
		
		// Check Frame using 3 product rule
		
		String frameString = String.valueOf(frame);
		int lenOfString = frameString.length();
		frame = 0;
		
		for(int i=0; i<lenOfString; i++){
			frame += (frameString.charAt(i)-48);
		}
		
		
		
		if(frame%3 == 0){
			return false;
		}
		
		return true;
	}
	
	private static String extractCDSRegion(ArrayList<Transcript> transcriptList, StartPos startPos, int textPos, boolean strand, int jncGap, int count) throws IOException {
		
		LinkedList<EXON> connectedExonList = startPos.startExon;
		int matchingLength = 0;
		int matchingStart = 0;
		int matchingEnd = 0;
		
		// CDS정보를 추출함
		EXON lastExon = null;
		int sizeOfExonList = connectedExonList.size();
		
		// Matching Region에 대한 정보를 저장
		// 짝수 인덱스: 시작 위치
		// 홀수 인덱스: 끝 위치
		int[] matchingRegion = new int[sizeOfExonList*2];
		
		
		for (int i=0; i<sizeOfExonList; i++) {
			lastExon = connectedExonList.get(i);
			
			
			if(sizeOfExonList == 1){
				matchingLength += textPos - startPos.startPos + 1;
				matchingStart = lastExon.get_start() + startPos.startPos;
				matchingEnd = lastExon.get_start() + textPos;
				
				matchingRegion[i*2] = lastExon.get_start() + startPos.startPos;
				matchingRegion[i*2+1] = lastExon.get_start() + textPos;
				
			}
			else if( i == 0){
				matchingStart = lastExon.get_start() + startPos.startPos;
				matchingLength += lastExon.get_end() - matchingStart + 1;
				
				matchingRegion[i*2] = lastExon.get_start() + startPos.startPos;
				matchingRegion[i*2+1] = lastExon.get_end();
			}
			else if(sizeOfExonList-1 == i){
				matchingEnd = lastExon.get_start() + textPos;
				matchingLength += textPos + 1;
				
				matchingRegion[i*2] = lastExon.get_start();
				matchingRegion[i*2+1] = lastExon.get_start() + textPos;
			}else{
				matchingLength += lastExon.get_end() - lastExon.get_start() + 1;
				
				matchingRegion[i*2] = lastExon.get_start();
				matchingRegion[i*2+1] = lastExon.get_end();
			}
			
		}

		// 각 transcript의 CDS구간이 Shared Peptide와 많이 겹치면 점수가 높음
		int score[] = new int[transcriptList.size()];
		
		for(int i=0; i<score.length; i++){
			score[i] = 0;
		}

		int transListIndex = 0;

		// transcript와 매칭된 패턴 사이의 similarity를 계산
		for (Transcript T : transcriptList) {
			ArrayList<ExonRangeType> ERT = T.exonList;
			boolean isMatchingRegion = false;
			int containScore = 0;
			
			for (int siteIndex = 0; siteIndex < ERT.size(); siteIndex++) {
				ExonRangeType curRange = ERT.get(siteIndex);				

				for(int i=0; i<sizeOfExonList; i++){
					// CDS-Start ---- matchingStart, matchingEnd ---- CDS-End
					if (curRange.start <= matchingRegion[i*2]  && curRange.end >= matchingRegion[i*2+1]) {
						if (curRange.isCDS) {
							score[transListIndex] += matchingRegion[i*2+1] - matchingRegion[i*2] +1;
						}
						containScore += matchingRegion[i*2+1] - matchingRegion[i*2] +1;
						isMatchingRegion = true;
					}
					// CDS-Start ---- matchingStart ---- CDS-End
					else if (curRange.start <= matchingRegion[i*2] && curRange.end >= matchingRegion[i*2]) {
						if (curRange.isCDS) {
							score[transListIndex] += curRange.end - matchingRegion[i*2] +1;
						}
						containScore += curRange.end - matchingRegion[i*2] +1;
						isMatchingRegion = true;
					}
					// CDS-Start ---- matchingEnd ---- CDS-End
					else if (curRange.start <= matchingRegion[i*2+1] && curRange.end >= matchingRegion[i*2+1]) {
						if (curRange.isCDS) {
							score[transListIndex] += matchingRegion[i*2+1] - curRange.start +1;
						}
						containScore += matchingRegion[i*2+1] - curRange.start +1;
						isMatchingRegion = true;
					}
					// matchingStart ---- CDS ---- matchingEnd
					else if (matchingRegion[i*2] <= curRange.start && matchingRegion[i*2+1] >= curRange.end) {
						if (curRange.isCDS) {
							score[transListIndex] += curRange.end - curRange.start +1;
						}
						containScore += curRange.end - curRange.start +1;
						isMatchingRegion = true;
					}
				}
				
			}
			
			
			// matchingRegion이 포함되지 않으면 -1로 가장 낮은 점수 부여
			if(!isMatchingRegion || (containScore + jncGap != matchingLength) ){
				score[transListIndex] = UNKNOWN;				
			}
			
			

			// transcript: ACGTCG-ACGTGT-ACGGAT가 있다고 가정하면,
			//                TCG-ACG   -ACG 와 같은 방식으로 맵핑되면 안됨
			// 즉, 완벽하게 exon에 맵핑되어야 함
			if(score[transListIndex] != UNKNOWN){
				for(int i=0; i<sizeOfExonList; i++){
					if(i+1 < sizeOfExonList && matchingRegion[i*2 + 1] +1 < matchingRegion[i*2 +2]){
						int prevEnd = matchingRegion[i*2 +1];
						int nextStart = matchingRegion[i*2 +2];
						boolean isPrevEnd = false;
						boolean isNextStart = false;
						for (int siteIndex = 0; siteIndex < ERT.size(); siteIndex++) {
							ExonRangeType curRange = ERT.get(siteIndex);
							
							if(curRange.end == prevEnd) isPrevEnd = true;
							if(curRange.start == nextStart) isNextStart = true;
							
							if(prevEnd < curRange.start && nextStart > curRange.start){
								score[transListIndex] = UNKNOWN;
								break;
							}
						}
						
						if(!(isPrevEnd && isNextStart)){
							score[transListIndex] = UNKNOWN;
							break;
						}
					}
				}
			}
			
			if(T.isPseudo() && score[transListIndex] != UNKNOWN){
				score[transListIndex] = PSEUDO;
			}
			
			transListIndex++;
			
		}

		// 가장 점수가 높은 transcript를 선정
		int bestScoreTranscriptIndex = 0;
		for (int i = 1; i < transcriptList.size(); i++) {
			if (score[i] > score[bestScoreTranscriptIndex]) {
				bestScoreTranscriptIndex = i;
			}
		}

		// CDS start and end (THE BEST CASE)
		String cdsInfo = null;
		for(int i= 0; i< transcriptList.size(); i++){
			if(score[i] == score[bestScoreTranscriptIndex]){
				cdsInfo = getRegionMappingType(transcriptList.get(i), score[i], matchingLength, matchingStart, matchingEnd, strand, true);				
			}
		}
		
		// LOG (NOT THE BEST CASE)
		for(int i=0; i<transcriptList.size(); i++){
			if(theBestMappingType == CDS || theBestMappingType == FS) break;
			
			// 모델이 없는 경우는 제외
			if(score[i] != UNKNOWN && score[i] != PSEUDO){
				
				// 5`utr과 3`utr의 경우만 Log로 출력함
				String mappingInfo = getRegionMappingType(transcriptList.get(i), score[i], matchingLength, matchingStart, matchingEnd, strand, false);
				
				if(mappingInfo.length() == 0){
					continue;
				}
				if(theCurrMappingType == theBestMappingType && 
						score[i] == score[bestScoreTranscriptIndex]){
					continue;
				}
				
				StringBuilder logs = new StringBuilder();
				logs.append(String.valueOf(geneID_)).append("\t");
				logs.append(String.valueOf(transcriptList.get(i).transcriptID)).append("\t");
				logs.append(mappingInfo);
				
				results.logResult.add(logs.toString()+"\n");
				
			}
		}
		
		if(cdsInfo.length() != 0 && count != 0){
			cdsInfo += ";";
		}
		return cdsInfo;
	}

	private static String getRegionMappingType(Transcript transcript, int score, int matchingLength, int matchingStart, int matchingEnd, boolean strand, boolean isBest){
		
		ArrayList<ExonRangeType> exonz = transcript.exonList;
		String cdsInfo = null;
		theCurrMappingType = 0;
		
		if (score == PSEUDO) {
			cdsInfo = "pseudo " + (ExonGraph.Chrom.replace("chr", "")) + " - - -";
			theCurrMappingType = PSEUDO;
		}
		// fully CDS!
		else if(score == matchingLength){
			
			if(isFS(transcript, matchingStart, matchingEnd, strand)){
				theCurrMappingType = FS;
			}else{
				theCurrMappingType = CDS;
			}
		}
		// 어떠한 경우에도 포함되지 않음
		// 즉, 모델에 없는 경우
		// unknown으로 표기
		// Intron의 경우는 unknown대신 intron으로 표기
		else if(score == UNKNOWN){
			if(Constants.IS_INTRON){
				cdsInfo = "intron";
			}else{
				cdsInfo = "unknown";
			}
			
			theCurrMappingType = UNKNOWN;
		}
		// 무조건 UTR!
		else {
			
			cdsInfo = "utr " + (ExonGraph.Chrom.replace("chr", "")) + " ";

			for (ExonRangeType ERT : exonz) {
				if(ERT.isCDS){
					
					if(ERT.start > matchingStart){
						if(strand){
							theCurrMappingType = UTR_5;
						}else{
							theCurrMappingType = UTR_3;
						}
						
					}else if(ERT.end < matchingEnd){
						if(strand){
							theCurrMappingType = UTR_3;
						}else{
							theCurrMappingType = UTR_5;
						}
						
					}
					
					break;
				}
			}
			int cdsStart = 0;
			int cdsEnd = 0;
			
			for(ExonRangeType ERT : exonz){
				if(ERT.isCDS){
					if(cdsStart == 0){
						cdsStart = ERT.start;
						cdsEnd = ERT.end;
					}else{
						cdsEnd = ERT.end;
					}
					
				}
			}
			
			if(theCurrMappingType == UTR_3){
				if(strand){
					cdsInfo += (cdsEnd) + " - x";
				}else{
					cdsInfo += (cdsStart) + " - x";
				}
			}else if(theCurrMappingType == UTR_5){
				if(strand){
					cdsInfo += (cdsStart) + " x -";
				}else{
					cdsInfo += (cdsEnd) + " x -";
				}
			}
			
		}
		if(isBest){
			theBestMappingType = theCurrMappingType > theBestMappingType ? theCurrMappingType : theBestMappingType;
			if(theBestMappingType == CDS){
				cdsInfo = "";
			}else if(theBestMappingType == FS){
				cdsInfo = "FS";
			}
		}else{
			
			if(theCurrMappingType != UTR_3 && theCurrMappingType != UTR_5){
				return "";
			}
			
		}
		
		return cdsInfo;
	}
	
	private static String MS_toString(ArrayList<MutationStructure_> MS_) {
		StringBuilder str = new StringBuilder();

		int msSize = MS_.size();
		for (int i = 0; i < msSize; i++) {
			MutationStructure_ ms = MS_.get(i);

			// 뮤테이션 타입을 기입.
			switch (ms.mutation) {
			case -2:
				str.append("jnc ");
				break;
			case -1:
				str.append("del ");
				break;
			case 0:
				str.append("snv ");
				break;
			case 1:
				str.append("ins ");
				break;
			case 2:
				str.append("alt ");
				break;
			default:
				Logger.getLogger("GFF").log(Level.SEVERE, "You should check your code related to GFF Maker.");
			}

			// Chr 번호 기입
			str.append(ExonGraph.Chrom.replace("chr", "")).append(" ");

			// start Position 기입
			str.append(ms.pos).append(" ");

			// 뮤테이션 타입마다 표시하는 게 다름
			if (ms.mutation >= -2 && ms.mutation <= 1) {
				str.append(ms.origin).append(" ").append(ms.substitution);
			} else {
				str.append("- -");
			}
			
			// 구분자 ;
			str.append(";");
						
		}

		return str.toString();
	}

	public static void close() throws IOException {
		int GFFIDIndex = 0;
		
		ArrayList<Result> resultList = new ArrayList<Result>(ResultMapping.ResultTable.values());
		
		for(Result result : resultList){
			GFFIDIndex++;
			BWVal.append(GFFIDIndex+"\t"+result.flatResult+"\n");
			
			for(String log : result.logResult){
				BWLog.append(GFFIDIndex+"\t"+log+"\n");
			}
		}
		
		
		if (BWVal != null) {
			BWVal.close();
			FWVal.close();
			
			BWLog.close();
			FWLog.close();
		}

	}

}
