package sixFrameTranslation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;

import ExonGraph.Codon;

public class acSearch {
	
	public int GFFID = 1;
	public String chr = null;
	public SFTOutput sftOut = null;
	public long totalFileSize = 0;
	public int totalProcess = 0;
	
	public void search(AhoCorasick ROOT, String fastaFilePath, String fastaName, long thisFileSize) throws IOException{
		if(totalProcess == 0){
			System.out.println("0 %");
		}
		
		if(sftOut == null){
			sftOut = new SFTOutput();
		}
		
		chr = fastaName.substring(0, fastaName.lastIndexOf("."));
		if(!chr.contains("chr")){
			chr = "chr"+chr;
		}
		
		FileReader FR= new FileReader(fastaFilePath);
		BufferedReader BR = new BufferedReader(FR);
		
		
		String line = null;
		int i,j,loc=0;
		String codonChar;
		
		LinkedList<String> match = new LinkedList<String>();

		StringBuilder SB = new StringBuilder();
		
		while((line = BR.readLine()) != null){
			if(line.startsWith(">"))
				continue;
			SB.append(line);
		}
		
		// Translation
		int length = SB.length();
		
		/*StringBuilder firstTranslation = new StringBuilder();
		StringBuilder secondTranslation = new StringBuilder();
		StringBuilder thirdranslation = new StringBuilder();
		
		int turn = 0;
		for(i=0; i<length-3; i++){
			codonChar = SB.substring(i, i+3);
			
			switch(turn){
			case 0: firstTranslation.append(Codon.NuclToAmino(codonChar)); turn++; break;
			case 1: secondTranslation.append(Codon.NuclToAmino(codonChar)); turn++; break;
			case 2: thirdranslation.append(Codon.NuclToAmino(codonChar)); turn = 0; break;
			}
		}
		
		turn = 0;
		firstTranslation.setLength(0);
		secondTranslation.setLength(0);
		thirdranslation.setLength(0);
		
		firstTranslation.trimToSize();
		secondTranslation.trimToSize();
		thirdranslation.trimToSize();
		
		System.out.println("-");
		for(i=length; i>2; i--){
			codonChar = SB.substring(i-3, i);
			
			switch(turn){
			case 0: firstTranslation.append(Codon.NuclToAmino_R_C(codonChar)); turn++; break;
			case 1: secondTranslation.append(Codon.NuclToAmino_R_C(codonChar)); turn++; break;
			case 2: thirdranslation.append(Codon.NuclToAmino_R_C(codonChar)); turn = 0; break;
			}
		}*/
		
		for(i = 0 ; i+5 <= length ; i = i + 3){
			
			//	f_1-frame
			codonChar = SB.substring(i, i+3);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 0);
			for(j=0; j<match.size(); j++){
				loc = i+4 - match.get(j).length() * 3 ;
				addOutput(loc, i+3, match.get(j), true, 0);
			}
			//	f_2-frame
			codonChar = SB.substring(i+1, i+4);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 1);
			for(j=0; j<match.size(); j++){
				loc = i+5 - match.get(j).length() * 3 ;
				addOutput(loc, i+4, match.get(j), true, 1);
			}
			//	f_3-frame
			codonChar = SB.substring(i+2, i+5);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 2);		
			for(j=0; j<match.size(); j++){
				loc = i+6 - match.get(j).length() * 3 ;
				addOutput(loc, i+5, match.get(j), true, 2);
			}
			
		}
		
		
		//	Exceptional Handling********************
		if(length>i+3){
			//	f_1-frame
			codonChar = SB.substring(i, i+3);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 0);	
			for(j=0; j<match.size(); j++){
				loc = i+4 - match.get(j).length() * 3 ;
				addOutput(loc, i+3, match.get(j), true, 0);
			}
			//	f_2-frame
			codonChar = SB.substring(i+1, i+4);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 1);
			for(j=0; j<match.size(); j++){
				loc = i+5 - match.get(j).length() * 3 ;
				addOutput(loc, i+4, match.get(j), true, 1);
			}
		}
		else if(length==i+3){
			//	f_1-frame
			codonChar = SB.substring(i, i+3);
			match = ROOT.searchString(Codon.NuclToAmino(codonChar), 0);
			for(j=0; j<match.size(); j++){
				loc = i+4 - match.get(j).length() * 3 ;
				addOutput(loc, i+3, match.get(j), true, 0);
			}
		}
		//******************************************
		
		
		ROOT.curState[0] = ROOT.ROOT;
		ROOT.curState[1] = ROOT.ROOT;
		ROOT.curState[2] = ROOT.ROOT;
		
		for(i = length ; i >= 5 ; i = i - 3){
			
			
			//	r_1-frame
			codonChar = SB.substring(i-3, i);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 0);
			for(j=0; j<match.size(); j++){
				loc = i-3 + match.get(j).length() * 3 ;
				addOutput(i-2, loc, match.get(j), false, 0);
			}
			//	r_2-frame
			codonChar = SB.substring(i-4, i-1);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 1);
			for(j=0; j<match.size(); j++){
				loc = i-4 + match.get(j).length() * 3 ;
				addOutput(i-3, loc, match.get(j), false, 1);
			}
			//	r_3-frame
			codonChar = SB.substring(i-5, i-2);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 2);
			for(j=0; j<match.size(); j++){
				loc = i-5 + match.get(j).length() * 3 ;
				addOutput(i-4, loc, match.get(j), false, 2);
			}
		}
		
		
		//	Exceptional Handling********************
		if(i>3){
			//	r_1-frame
			codonChar = SB.substring(i-3, i);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 0);
			for(j=0; j<match.size(); j++){
				loc = i-3 + match.get(j).length() * 3 ;
				addOutput(i-2, loc, match.get(j), false, 0);
			}
			//	r_2-frame
			codonChar = SB.substring(i-4, i-1);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 1);
			for(j=0; j<match.size(); j++){
				loc = i-4 + match.get(j).length() * 3 ;
				addOutput(i-3, loc, match.get(j), false, 1);
			}
		}
		else if(i==3){
			//	r_1-frame
			codonChar = SB.substring(i-3, i);
			match = ROOT.searchString(Codon.NuclToAmino_R_C(codonChar), 0);
			for(j=0; j<match.size(); j++){
				loc = i-3 + match.get(j).length() * 3 ;
				addOutput(i-2, loc, match.get(j), false, 0);
			}
		}
		//******************************************
		totalProcess += (int)100*((double)thisFileSize/(double)totalFileSize);
		System.out.println(totalProcess+" %");
		
		
		BR.close();
		
	}

	public void addOutput(int start, int end, String peptide, boolean strand, int frame) throws IOException{
		if(sftOut == null){
			System.out.println("ERROR");
			return;
		}
		
		sftOut.writerFlat(GFFID, peptide, chr);
		sftOut.writeGFF(chr, start, end, strand, frame, GFFID);
		GFFID++;
	}
	
	public void done() throws IOException{
		if(sftOut != null){
			sftOut.close();
		}
	}
}