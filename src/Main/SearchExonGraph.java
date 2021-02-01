package Main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedList;
import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;

import Environments.Codon;
import Environments.Constants;
import ExonGraph.ExonGraph;
import ExonGraph.ExonGraphGF;
import ExonGraph.Flat;
import ExonGraph.GFF;
import ExonGraph.KAminoNode;
import ExonGraph.KAminoTree;
import ExonGraph.Pattern;
import ExonGraph.Protein;
import proteinDB.ProteinDBFlat;
import proteinDB.ProteinMappingWithSAVs;
import sixFrameTranslation.AhoCorasick;
import sixFrameTranslation.acSearch;

public class SearchExonGraph {

	
	public static void main(String[] args) throws IOException, ClassNotFoundException{
		if(args.length == 2){
			System.setOut(new PrintStream(new FileOutputStream(args[1])));
		}
		
		long startTime= System.currentTimeMillis();
		
		String params = args[0];
		
		try {
			Constants.readParams(params, Constants.SEARCH_PHASE);
		} catch (ParserConfigurationException | SAXException e) {
			// TODO Auto-generated catch block
			Environments.Error.exitError(Environments.Error.PARAM_ERROR, "XML_PARSER_ERROR");
		}
		
		String peptideListPath = Constants.PEPTIDE_LIST_PATH;
		ExonGraphGF EG = null;
		
		if(Constants.IS_VSG){
			EG = loadGraph(Constants.GRAPH_PATH);
			
			if(EG == null){
				Environments.Error.exitError(Environments.Error.FACTAL_ERROR, "GRAPH_IS_INVALID");
			}
		}

		// Codon Table
		Codon.Mapping();
		
		
		KAminoTree tmp = new KAminoTree(50, 30);
		FileReader peptideListFR = new FileReader(peptideListPath);
		BufferedReader peptideListBR = new BufferedReader(peptideListFR);
		LinkedList<String> peptideList = new LinkedList<String>();
		String peptide = null;
		while((peptide = peptideListBR.readLine()) != null) {
			
			// FORMAT CHECKING
			char[] peptideValidation = peptide.toCharArray();
			for(int i=0; i<peptideValidation.length; i++){
				
				if(peptideValidation[i] < 'A' || peptideValidation[i] > 'Z'){
					Environments.Error.exitError(Environments.Error.FORMAT_ERROR, "PEPTIDE_FORMAT_IS_INVALID");
				}
			}
			
			peptideList.add(peptide);
			Pattern pat1 = null;
			Pattern pat2 = null;
			if(Environments.Constants.IS_I_SAME_WITH_L){
				pat1 = new Pattern(peptide.replaceAll("I", "L"), peptide, peptide);
				pat2 = new Pattern(new StringBuilder(peptide.replaceAll("I", "L")).reverse().toString(), new StringBuilder(peptide).reverse().toString(), peptide);
			}else{
				pat1 = new Pattern(peptide, peptide, peptide);
				pat2 = new Pattern(new StringBuilder(peptide).reverse().toString(), new StringBuilder(peptide).reverse().toString(), peptide);
			}
			
			
			tmp.BuildKeyword(pat1, tmp.Tree);
			tmp.BuildKeyword(pat2, tmp.RTree);
		}
		
		peptideListBR.close();
		peptideListFR.close();
		
		tmp.BuildFlink(); tmp.BuildOlink(); tmp.NFAtoDFA();
		
		// ProteinDB �꽌移섎�� 怨좊젮�븯硫�, ProteinDB�꽌移섎�� �닔�뻾�븯怨� �궎�썙�뱶 �듃由щ�� �깉濡� 留뚮벀
		if(Constants.IS_PROTEIN_DB){
			System.out.println("Mapping on the protein database");
			ProteinDBFlat mappedPeptides = extractUnmappedPeptideListFromProteinDB(tmp.Tree, peptideList);
			
			tmp = new KAminoTree(50, 30);
			
			// Remove already mapped peptides
			for(int i=0; i<peptideList.size(); i++){
				String pep = peptideList.get(i);
				for(String p : mappedPeptides.seqList){
					if(pep.equalsIgnoreCase(p)){
						peptideList.remove(i--);
						break;
					}
				}
			}
			
			if(Constants.IS_VSG){
				for(String originalPeptide : peptideList){
					String replacedPeptide = null;
					if(Environments.Constants.IS_I_SAME_WITH_L){
						replacedPeptide = originalPeptide.replaceAll("I", "L");
					}else{
						replacedPeptide = originalPeptide;
					}
					
					Pattern pat1 = new Pattern(replacedPeptide, originalPeptide, originalPeptide);
					Pattern pat2 = new Pattern(new StringBuilder(replacedPeptide).reverse().toString(), new StringBuilder(originalPeptide).reverse().toString(), originalPeptide);
					tmp.BuildKeyword(pat1, tmp.Tree);
					tmp.BuildKeyword(pat2, tmp.RTree);
				}
				
				tmp.BuildFlink(); tmp.BuildOlink(); tmp.NFAtoDFA();
			}
			System.out.println("Mapping on the protein database is complete");
		}
		
		if(Constants.IS_SFT){
			AhoCorasick ROOT = new AhoCorasick();
			int pepSize = peptideList.size();
			for(int i=0; i<pepSize; i++){
				String peptideSeq = peptideList.get(i);
				if(Constants.IS_I_SAME_WITH_L){
					peptideSeq = peptideSeq.replaceAll("I", "L");
				}
				
				ROOT.addPattern(peptideSeq, peptideList.get(i));
			}
			ROOT.buildFailState();
			
			//TODO: Init SFT
			System.out.println("Mapping on the six-frame translation database");
			File[] fastaFileList = new File(Constants.REFERENCE_GENOME_PATH).listFiles();
			acSearch ac = new acSearch();
			long totalSize = 0;
			
			for(int i=0; i<fastaFileList.length; i++){
				String fastaName = fastaFileList[i].getName().toUpperCase();
				if(!fastaName.contains(".FA")){
					Environments.Error.exitError(Environments.Error.FORMAT_ERROR, "FASTA_FILE_HAS_TO_BE_.FA");
				}
				totalSize += fastaFileList[i].getTotalSpace();
			}
			
			ac.totalFileSize = totalSize;
			for(int i=0; i<fastaFileList.length; i++){
				ac.search(ROOT, fastaFileList[i].getPath(), fastaFileList[i].getName(), fastaFileList[i].getTotalSpace());
			}
			
			//TODO: Search on SFT
			ac.done();
			System.out.println("Mapping on the six-frame translation database is complete");
		}
	
		if(EG != null){
			System.out.println("Mapping on the variant splice graph");
			ExonGraph.nowStatus = !ExonGraph.nowStatus;
			EG.AnimoExonSearch(ExonGraph.nowStatus, tmp);
			System.out.println("Mapping on the variant splice graph is complete");
		}
		
		System.out.println("Elapsed Time : " + (System.currentTimeMillis()-startTime)/1000 + " Sec" );
		
		GFF.close();
		Flat.close();
	}
	
	public static ExonGraphGF loadGraph(String fileName) throws ClassNotFoundException, IOException{
		return new ExonGraphGF(fileName);
	}
	
	
	/**
	 * protein fasta �뙆�씪�쓣 �씫怨�, ac_tree(peptide list�젙蹂대�� ���옣)�� 留듯븨�릺�뒗 
	 * peptides瑜� 諛섑솚�븿.
	 * 
	 * @param ac_tree
	 * @return
	 * @throws IOException
	 * 	 * 
	 */

	public static ProteinDBFlat extractUnmappedPeptideListFromProteinDB(KAminoNode ac_tree, LinkedList<String> peptideList) throws IOException{
		
		File[] proteinDBList = new File(Constants.PROTEIN_DB_PATH).listFiles();
		if(proteinDBList == null){
			proteinDBList = new File[1];
			proteinDBList[0] = new File(Constants.PROTEIN_DB_PATH);
		}
		
		ProteinDBFlat mappedPeptides = new ProteinDBFlat();
		ProteinDBFlat mappedSAVPeptides = null;
		if(Constants.IS_SAV) mappedSAVPeptides = new ProteinDBFlat();
		ProteinMappingWithSAVs PSWS = null;
		if(Constants.IS_SAV) PSWS = new ProteinMappingWithSAVs();
		
		for(File proteinDBFile : proteinDBList){
			FileReader FR = new FileReader(proteinDBFile);
			BufferedReader BR = new BufferedReader(FR);			
			String line = null;
			StringBuilder proteinSeq = new StringBuilder();
			int index = 0;
			ArrayList<Protein> proteinDB = new ArrayList<Protein>();
			
			// ProteinDB瑜� �씫�쓬
			while((line = BR.readLine())!=null){
				if(line.startsWith(">")){
					if(proteinSeq.length() != 0){
						proteinDB.get(index++).proteinSeq = proteinSeq.toString();
						proteinSeq.delete(0, proteinSeq.length());
					}
					
					proteinDB.add(new Protein());
					proteinDB.get(index).header = line;
				}else{
					proteinSeq.append(line);
				}
			}
			// DB 留덉�留� �씪�씤 異붽�
			if(proteinSeq.length() != 0){
				proteinDB.get(index).proteinSeq = proteinSeq.toString();
				proteinSeq.delete(0, proteinSeq.length());
			}
			
			proteinDBMapping(ac_tree, proteinDB, mappedPeptides);
			if(PSWS != null) PSWS.proteinDBMappingWithSAVs(peptideList, proteinDB, mappedSAVPeptides);
			
			
			BR.close();
			FR.close();
		}
		
		mappedPeptides.writeMapping(Constants.OUTPUT_PROTEIN_NORMAL);
		if(mappedSAVPeptides != null){
			mappedSAVPeptides.writeMapping(Constants.OUTPUT_PROTEIN_SAV);
			// Merge SAV and NOR mapping results into one
			for(String seqSAV : mappedSAVPeptides.seqList){
				for(String headerSAV : mappedSAVPeptides.headerMapping.get(seqSAV)){
					mappedPeptides.putMapping(headerSAV, seqSAV);
				}
			}
		
		}
		return mappedPeptides;
	}
	
	public static void proteinDBMapping(KAminoNode ac_tree, ArrayList<Protein> proteinDB, ProteinDBFlat mappedPeptides){
		for(Protein protein : proteinDB){
			int treePos = 0;
			int index = 0;
			String proteinSeq = null;
			String proteinHeader = protein.header;
			
			if(Environments.Constants.IS_I_SAME_WITH_L){
				proteinSeq = protein.proteinSeq.replaceAll("I", "L"); 
			}else{
				proteinSeq = protein.proteinSeq;
			}
			
			int proteinSeqLen = proteinSeq.length();
			while(index < proteinSeqLen){
				char nucl = proteinSeq.charAt(index);
				if(ac_tree.get(treePos).get_next(nucl) == ac_tree.get(ac_tree.get(treePos).get_failstate()).get_next(nucl)) {

				}
				
				treePos = ac_tree.get(treePos).get_next(nucl);
				
				if(ac_tree.get(treePos).get_matchlist() != null){
					Pattern TPM = ac_tree.get(treePos).get_matchlist();
					
					for(int i=0; i<=TPM.size(); i++){
						String mappedPeptide = null;
						
						if(i == 0){
							mappedPeptide = TPM.get_origin();
						}else{
							mappedPeptide = TPM.get(i-1).get_origin();
						}
						mappedPeptides.putMapping(proteinHeader, mappedPeptide);
						
					}
					
				}
				
				index++;
			}
		}
		
		return;
	}
}
