package sixFrameTranslation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import Environments.Constants;

public class SFTOutput {
	
	public FileWriter FWFlat=null, FWGFF=null;
	public BufferedWriter BWFlat=null, BWGFF=null;
	
	public SFTOutput() throws IOException{
		FWFlat = new FileWriter(Constants.OUTPUT_SFT_FLAT);
		FWGFF = new FileWriter(Constants.OUTPUT_SFT_GFF);
		
		BWFlat = new BufferedWriter(FWFlat);
		BWGFF = new BufferedWriter(FWGFF);
		
		BWFlat.append("GFFID\tPeptide\tChr");
		BWFlat.newLine();
	}
	
	public void writerFlat(int GFFID, String peptide, String chr) throws IOException{		
		BWFlat.append(GFFID+"\t"+peptide+"\t"+chr.replaceAll("chr", ""));
		BWFlat.newLine();
	}
	
	public void writeGFF(String chr, int start, int end, boolean strand, int frame, int GFFID) throws IOException{
		String strandString = "-";
		if(strand){
			strandString = "+";
		}
		
		BWGFF.append(chr+"\tACTG\t.\t"+start+"\t"+end+"\t.\t"+strandString+"\t"+frame+"\t"+"ID="+GFFID);
		BWGFF.newLine();
	}
	
	public void close() throws IOException{
		if(BWFlat != null){
			BWFlat.close(); FWFlat.close();
		}
		if(BWGFF != null){
			BWGFF.close(); FWGFF.close();
		}
	}
}
