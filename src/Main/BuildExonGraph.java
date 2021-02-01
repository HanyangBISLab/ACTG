package Main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import Environments.Codon;
import Environments.Constants;
import Thread.ExonGraphGFT;
import data.GTF;
import data.KeywordTree;
import data.Pattern;

public class BuildExonGraph {
	
	public static int TheNumberOfSuccess = ExonGraphGFT.TheNumberOfSuccess;
	public static int TheNumberOfTasks = 0;
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length != 1){
			System.exit(-1);
		}
		
		String params = args[0];
		
		try {
			Constants.readParams(params, Constants.CONSTRUCTION_PHASE);
		} catch (ParserConfigurationException | SAXException e) {
			// TODO Auto-generated catch block
			Environments.Error.exitError(Environments.Error.PARAM_ERROR, "XML_PARSER_ERROR");
		}
		long startTime= System.currentTimeMillis();

		File[] fastaFileList = new File(Constants.REFERENCE_GENOME_PATH).listFiles();
		File[] gtfFileListCheck = new File(Constants.GTF_PATH).listFiles();
		File[] gtfFileList = new File[gtfFileListCheck.length];
		
		for(int i=0; i<fastaFileList.length; i++){
			String fastaName = fastaFileList[i].getName().toUpperCase();
			if(!fastaName.contains(".FA")){
				Environments.Error.exitError(Environments.Error.FORMAT_ERROR, "FASTA_FILE_HAS_TO_BE_.FA");
			}
			
			String identifier = fastaName.substring(0, fastaName.lastIndexOf("."));
			
			for(int j=0; j<gtfFileListCheck.length; j++){
				String gtfName = gtfFileListCheck[j].getName().toUpperCase();
				if(!gtfName.contains(".GTF")){
					Environments.Error.exitError(Environments.Error.FORMAT_ERROR, "GTF_FILE_HAS_TO_BE_.GTF");
				}
				
				if(gtfName.substring(0, gtfName.lastIndexOf(".")).equalsIgnoreCase(identifier)){
					gtfFileList[i] = gtfFileListCheck[j];
					break;
				}
			}
		}
		
		if(gtfFileListCheck.length == 1){
			gtfFileList[0] = gtfFileListCheck[0];
		}
		
		
		// load gtf
		
		// build keyword tree
		ArrayList<Pattern> patterns = loadPatterns("C:\\Users\\progi\\Desktop\\Projects\\ACTG\\test.txt");
		KeywordTree ktree = new KeywordTree(patterns);
		Codon.Mapping();
		for(File file : gtfFileList) {
			GTF gtf = new GTF(file);
			gtf.setSequence(fastaFileList);
			gtf.find(ktree);
		}
		
		System.out.println("Construction Time for Exon Graph : " + (System.currentTimeMillis()-startTime)/1000 + " Sec" );

	
	}
	
	public static ArrayList<Pattern> loadPatterns (String fileName) throws IOException {
		ArrayList<Pattern> patterns = new ArrayList<Pattern>();
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			Pattern pattern = new Pattern();
			pattern.aaSequence = line;
			patterns.add(pattern);
		}
		BR.close();
		
		return patterns;
	}
	
}
