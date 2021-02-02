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
	private static final int NONCODING = -1;
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
	private static String CHR = null;
	
	
	public static synchronized void write(StartPos startPos, int textPos, boolean strand, int trans_cnt, 
			String valText, ArrayList<Transcript> transcriptList, String geneID, String chrID, String nucleotides) throws IOException {
		GFF.write(startPos, textPos, strand, trans_cnt, chrID);
		
		CHR = chrID;
		theBestMappingType = UNKNOWN;
		geneID_ = geneID;
		
		// 理쒖큹 1�쉶 �뙆�씪�쓣 �깮�꽦�븿
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
			
			// Index瑜� 異붽��븿
			BWVal.append("GFFID\tPeptide\tNucleotide\tGeneID\tAttribute\n");
			BWLog.append("GFFID\tGeneID\tNucleotide\tTranscriptID\tAttribute\n");
		}

		
		int sizeOfExon = startPos.startExon.size();

		StringBuilder eventText = new StringBuilder();
		StringBuilder normalText = new StringBuilder();
		int count = 0;
		int jncGap = 0;
		boolean isJncOrAlt = false;
		// jnc 蹂��씠�뿉 ���븳 以묐났�젣嫄곕�� �쐞�븿
		
		
		// Event Extractor
		ExonUtils.variationInit();
		for (int i = 0; i < sizeOfExon; i++) {
			ArrayList<MutationStructure_> variation = ExonUtils.getVariation(startPos, i, trans_cnt, strand);
						
			if (variation.size() != 0) {
				count++;
				eventText.append(MS_toString(variation));
				for(int j=0; j<variation.size(); j++){
					// jnc�쓽 寃쎌슦 蹂��삎�맂 nucleotide �닔 留뚰겮 �뙣�꼸�떚瑜� 以�
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
		// Nucleotide Sequences
		normalText.append(nucleotides).append("\t");
		// GENEID
		normalText.append(String.valueOf(geneID_)).append("\t");
		// CDSRegion
		if(!isJncOrAlt)
			normalText.append(extractCDSRegion(transcriptList, startPos, textPos, strand, jncGap, count));
		// Event Information
		normalText.append(eventText.toString());
		
		StringBuilder resultKey = null;
		
		if(Constants.EQUIVALENT_TEST_USING_GENOMIC_LOCATION){
			resultKey = new StringBuilder();
		}else{
			resultKey = new StringBuilder(normalText.toString());
		}
		
		// chr湲곗엯
		results.chr = GFF.temporaryGFF.get(0).split("\t")[0];
		for(String gff : GFF.temporaryGFF){
			
			if(Constants.EQUIVALENT_TEST_USING_GENOMIC_LOCATION){
				String[] spliter = gff.split("\t");
				resultKey.append(spliter[3]).append("\t").append(spliter[4]).append("\t");
			}else{
				resultKey.append(gff);
			}
			
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
					frame = Math.abs(cumulativePos + matchingStart - curERT.start);
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
		
		// CDS�젙蹂대�� 異붿텧�븿
		EXON lastExon = null;
		int sizeOfExonList = connectedExonList.size();
		
		// Matching Region�뿉 ���븳 �젙蹂대�� ���옣
		// 吏앹닔 �씤�뜳�뒪: �떆�옉 �쐞移�
		// ���닔 �씤�뜳�뒪: �걹 �쐞移�
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

		// 媛� transcript�쓽 CDS援ш컙�씠 Shared Peptide�� 留롮씠 寃뱀튂硫� �젏�닔媛� �넂�쓬
		int score[] = new int[transcriptList.size()];
		
		for(int i=0; i<score.length; i++){
			score[i] = 0;
		}

		int transListIndex = 0;

		// transcript�� 留ㅼ묶�맂 �뙣�꽩 �궗�씠�쓽 similarity瑜� 怨꾩궛
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
			
			
			// matchingRegion�씠 �룷�븿�릺吏� �븡�쑝硫� -1濡� 媛��옣 �궙�� �젏�닔 遺��뿬
			if(!isMatchingRegion || (containScore + jncGap != matchingLength) ){
				score[transListIndex] = UNKNOWN;				
			}
			
			

			// transcript: ACGTCG-ACGTGT-ACGGAT媛� �엳�떎怨� 媛��젙�븯硫�,
			//                TCG-ACG   -ACG �� 媛숈� 諛⑹떇�쑝濡� 留듯븨�릺硫� �븞�맖
			// 利�, �셿踰쏀븯寃� exon�뿉 留듯븨�릺�뼱�빞 �븿
			/*if(score[transListIndex] != UNKNOWN){
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
			}*/
			
			if(T.isPseudo() && score[transListIndex] != UNKNOWN){
				score[transListIndex] = NONCODING;
			}
			
			transListIndex++;
			
		}

		// 媛��옣 �젏�닔媛� �넂�� transcript瑜� �꽑�젙
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
			
			// 紐⑤뜽�씠 �뾾�뒗 寃쎌슦�뒗 �젣�쇅
			if(score[i] != UNKNOWN && score[i] != NONCODING){
				
				// 5`utr怨� 3`utr�쓽 寃쎌슦留� Log濡� 異쒕젰�븿
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
		if (score == NONCODING) {
			cdsInfo = "noncoding(chr" + (CHR.replace("chr", "")) + ")";
			theCurrMappingType = NONCODING;
		}
		// fully CDS!
		else if(score == matchingLength){
			
			if(isFS(transcript, matchingStart, matchingEnd, strand)){
				theCurrMappingType = FS;
			}else{
				theCurrMappingType = CDS;
			}
		}
		// �뼱�뼚�븳 寃쎌슦�뿉�룄 �룷�븿�릺吏� �븡�쓬
		// 利�, 紐⑤뜽�뿉 �뾾�뒗 寃쎌슦
		// unknown�쑝濡� �몴湲�
		// Intron�쓽 寃쎌슦�뒗 unknown���떊 intron�쑝濡� �몴湲�
		else if(score == UNKNOWN){
			if(Constants.IS_INTRON){
				cdsInfo = "exon-extension";
			}else{
				cdsInfo = "unknown";
			}
			
			theCurrMappingType = UNKNOWN;
		}
		// 臾댁“嫄� UTR!
		else {
			
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
					cdsInfo = "3`-UTR(chr" + CHR.replace("chr", "")+ ":g."+(cdsEnd) + ")";
				}else{
					cdsInfo = "3`-UTR(chr" + CHR.replace("chr", "")+ ":g."+(cdsStart) + ")";
				}
			}else if(theCurrMappingType == UTR_5){
				if(strand){
					cdsInfo = "5`-UTR(chr" + CHR.replace("chr", "")+ ":g."+(cdsStart) + ")";
				}else{
					cdsInfo = "5`-UTR(chr" + CHR.replace("chr", "")+ ":g."+(cdsEnd) + ")";
				}
			}
			
		}
		if(isBest){
			theBestMappingType = theCurrMappingType > theBestMappingType ? theCurrMappingType : theBestMappingType;
			if(theBestMappingType == CDS){
				cdsInfo = "CDS";
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

			// 裕ㅽ뀒�씠�뀡 ���엯�쓣 湲곗엯.
			switch (ms.mutation) {
			case -2:
				str.append("jnc");
				break;
			case -1:
				str.append("del");
				break;
			case 0:
				str.append("snv");
				break;
			case 1:
				str.append("ins");
				break;
			case 2:
				str.append("exon-skipping");
				break;
			default:
				Logger.getLogger("GFF").log(Level.SEVERE, "You should check your code related to GFF Maker.");
			}

			// Chr 踰덊샇 湲곗엯
			str.append("(chr").append(CHR.replace("chr", "")).append(":g.");

			// start Position 湲곗엯
			str.append(ms.pos);

			// 裕ㅽ뀒�씠�뀡 ���엯留덈떎 �몴�떆�븯�뒗 寃� �떎由�
			if (ms.mutation >= -2 && ms.mutation <= 1) {
				str.append(ms.origin).append(">").append(ms.substitution);
			}str.append(")");
			
			// 援щ텇�옄 ;
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
