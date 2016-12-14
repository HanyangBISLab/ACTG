package Environments;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import ExonGraph.ExonGraphGF;
import Thread.ExonGraphGFT;

public class Constants {
	public static final String VERSION = "1.10";
	
	public static String PEPTIDE_LIST_PATH = null;
	public static String PROTEIN_DB_PATH = null;
	public static String VCF_PATH = null;
	public static String REFERENCE_GENOME_PATH = null;
	public static String GTF_PATH = null;
	public static String GRAPH_PATH = null;
	
	public static String GRAPH_OUTPUT_PATH = null;
	public static String OUTPUT_ROOT = null;
	public static String OUTPUT_PROTEIN_NORMAL = null;
	public static String OUTPUT_PROTEIN_SAV = null;
	public static String OUTPUT_VSG_FLAT = null;
	public static String OUTPUT_VSG_GFF = null;
	public static String OUTPUT_VSG_LOG = null;
	public static String OUTPUT_SFT_FLAT = null;
	public static String OUTPUT_SFT_GFF = null;
	
	public static final int CONSTRUCTION_PHASE = 0;
	public static final int SEARCH_PHASE = 1;
	
	public static int[] SEARCH_EDGES = null;
	public static boolean IS_MUTATION = false;
	public static boolean IS_ALTERNATIVE = false;
	public static boolean IS_JUNCTION = false;
	
	public static boolean IS_PROTEIN_DB = false;
	public static boolean IS_VSG = false;
	public static boolean IS_SFT = false;
	
	public static boolean IS_I_SAME_WITH_L = false;
	public static boolean IS_INTRON = false;
	public static boolean IS_SAV = false;
	
	public static String MAPPING_METHOD = null;
	
	public static boolean TRANSLATION_END_OF_TRANSCRIPT = true;
	public static boolean EQUIVALENT_TEST_USING_GENOMIC_LOCATION = false;
	
	//EDGE INFO
	public static final int INTRON_EDGE = 13;
	
	//THREAD INFO
	public static int THE_NUMBER_OF_THREADS = 3;
	
	public static int readParams(String params, int case_) throws IOException, ParserConfigurationException, SAXException{
		System.out.println("ACTG v"+VERSION);
		System.out.println();
		
		File file = new File(params);
		

		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document doc = builder.parse(file);
		doc.getDocumentElement().normalize();
		
		NodeList nl = null;
		Node node = null;
		String type = null;

		int theNumberOfSearchEdges = ExonGraphGF.JVALUE;
		
		
		nl = doc.getElementsByTagName("TheNumberOfThreads");
		THE_NUMBER_OF_THREADS = Integer.parseInt(nl.item(0).getTextContent());
		
		System.out.println("<System setting>");
		System.out.println(" The number of threads: "+THE_NUMBER_OF_THREADS);
			
		if(case_ == SEARCH_PHASE){
			
			//MAPPING METHOD
			nl = doc.getElementsByTagName("MappingMethod");
			node = nl.item(0);
			MAPPING_METHOD = node.getTextContent();
			if(MAPPING_METHOD.contains("P")){
				IS_PROTEIN_DB = true;
			}
			if(MAPPING_METHOD.contains("V")){
				IS_VSG = true;
			}
			if(MAPPING_METHOD.contains("S")){
				IS_SFT = true;
			}
			
			//IL
			nl = doc.getElementsByTagName("ILSame");
			node  = nl.item(0);
			if(node.getTextContent().equalsIgnoreCase("yes")){
				IS_I_SAME_WITH_L = true;
			}
			
			//PROTEIN DB
			if(IS_PROTEIN_DB){
				nl = doc.getElementsByTagName("SAVs");
				node = nl.item(0);
				if(node.getTextContent().equalsIgnoreCase("yes")){
					IS_SAV = true;
				}
			}
			
			
			//VARIANT SPLICE GRAPH
			if(IS_VSG){
				nl = doc.getElementsByTagName("JunctionVariation");
				node = nl.item(0);
				
				if(node.getTextContent().equalsIgnoreCase("yes")){
					IS_JUNCTION = true;
				}else{
					theNumberOfSearchEdges -= 3;
				}
				
				nl = doc.getElementsByTagName("ExonSkipping");
				node = nl.item(0);

				if(node.getTextContent().equalsIgnoreCase("yes")){
					IS_ALTERNATIVE = true;
				}else{
					theNumberOfSearchEdges -= 1;
				}
				
				nl = doc.getElementsByTagName("Mutation");
				node = nl.item(0);

				if(node.getTextContent().equalsIgnoreCase("yes")){
					IS_MUTATION = true;
				}else{
					theNumberOfSearchEdges -= 9;
				}
											
				nl = doc.getElementsByTagName("IntronMapping");
				node = nl.item(0);

				if(node.getTextContent().equalsIgnoreCase("yes")){
					IS_INTRON = true;
				}else{
					theNumberOfSearchEdges -=1;
				}
				
			}
		}
		
		// 이벤트를 고려한 Edges 배열
		SEARCH_EDGES = new int[theNumberOfSearchEdges];
		int edgeIndex = 0;
		for(int i=0; i<ExonGraphGF.JVALUE; i++){
			if(i == ExonGraphGF.JVALUE-1 && IS_ALTERNATIVE){
				SEARCH_EDGES[edgeIndex++] = i;
			}else if( (i > 0 && i < 4) && IS_JUNCTION){
				SEARCH_EDGES[edgeIndex++] = i;
			}else if( (i > 3 && i < ExonGraphGF.JVALUE-2) && IS_MUTATION){
				SEARCH_EDGES[edgeIndex++] = i;
			}else if(i == 0){
				SEARCH_EDGES[edgeIndex++] = i;
			}else if(i == INTRON_EDGE && IS_INTRON){
				SEARCH_EDGES[edgeIndex++] = i;
			}
		}
						
		nl = doc.getElementsByTagName("Input");
		
		if(case_ == CONSTRUCTION_PHASE){
			for(int i=0; i<nl.getLength(); i++){
				node = nl.item(i);
				type = node.getAttributes().getNamedItem("type").getNodeValue();
				
				if(type.equalsIgnoreCase("transcriptome")){
					GTF_PATH = node.getTextContent();
				}else if(type.equalsIgnoreCase("referenceGenome")){
					REFERENCE_GENOME_PATH = node.getTextContent();
				}
				
			}
		}else if(case_ == SEARCH_PHASE){
			for(int i=0; i<nl.getLength(); i++){
				node = nl.item(i);
				type = node.getAttributes().getNamedItem("type").getNodeValue();
				
				if(IS_VSG && type.equalsIgnoreCase("graphFile")){
					GRAPH_PATH = node.getTextContent();
				}else if(IS_PROTEIN_DB && type.equalsIgnoreCase("proteinDB")){
					PROTEIN_DB_PATH = node.getTextContent();
				}else if(type.equalsIgnoreCase("peptideList")){
					PEPTIDE_LIST_PATH = node.getTextContent();
				}else if(IS_MUTATION && type.equalsIgnoreCase("mutation")){
					VCF_PATH = node.getTextContent();
				}else if(IS_SFT && type.equalsIgnoreCase("referenceGenome")){
					REFERENCE_GENOME_PATH = node.getTextContent();
				}
			}
		}
		
		
		nl = doc.getElementsByTagName("Output");
		
		if(case_ == CONSTRUCTION_PHASE){
			for(int i=0; i<nl.getLength(); i++){
				node = nl.item(i);
				type = node.getAttributes().getNamedItem("type").getNodeValue();
				
				if(type.equalsIgnoreCase("GraphFile")){
					GRAPH_OUTPUT_PATH = node.getTextContent();
				}
			}
		}else if(case_ == SEARCH_PHASE){
			node = nl.item(0);
			OUTPUT_ROOT = node.getTextContent();
			File pepFile = new File(PEPTIDE_LIST_PATH);
			if(!pepFile.exists()){
				Error.exitError(Error.NO_SUCH_A_FILE, "PEPTIDE_LIST_PATH");
			}
			
			String peptideListName = pepFile.getName();
			peptideListName = peptideListName.substring(0, peptideListName.lastIndexOf("."));
			
			OUTPUT_PROTEIN_NORMAL = OUTPUT_ROOT+"/NOR_"+peptideListName+".txt";
			OUTPUT_PROTEIN_SAV = OUTPUT_ROOT+"/SAV_"+peptideListName+".txt";
			OUTPUT_VSG_FLAT = OUTPUT_ROOT+"/VSG_"+peptideListName+".flat";
			OUTPUT_VSG_GFF = OUTPUT_ROOT+"/VSG_"+peptideListName+".gff";
			OUTPUT_VSG_LOG = OUTPUT_ROOT+"/VSG_"+peptideListName+".log";
			OUTPUT_SFT_FLAT = OUTPUT_ROOT+"/SFT_"+peptideListName+".flat";
			OUTPUT_SFT_GFF = OUTPUT_ROOT+"/SFT_"+peptideListName+".gff";
			
		}
		
		
		//Validation
		File isFile = null;
		if(case_ == CONSTRUCTION_PHASE){
			
			if(REFERENCE_GENOME_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "GENOME_REFERENCE_PATH");
			}else{
				isFile = new File(REFERENCE_GENOME_PATH);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, REFERENCE_GENOME_PATH);
				}
				
				if(!isFile.isDirectory()){
					Error.exitError(Error.IS_NOT_A_DIRECTORY, REFERENCE_GENOME_PATH);
				}
				
				if(isFile.list().length == 0){
					Error.exitError(Error.IS_EMPTY, REFERENCE_GENOME_PATH);
				}
			}
			
			if(GTF_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "GTF_PATH");
			}else{
				isFile = new File(GTF_PATH);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, GTF_PATH);
				}
				
				if(!isFile.isDirectory()){
					Error.exitError(Error.IS_NOT_A_DIRECTORY, GTF_PATH);
				}
				
				if(isFile.list().length == 0){
					Error.exitError(Error.IS_EMPTY, GTF_PATH);
				}
			}
			
			if(GRAPH_OUTPUT_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "GRAPH_OUTPUT_PATH");
			}
			
			
		}else if(case_ == SEARCH_PHASE){
			if(PEPTIDE_LIST_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "PEPTIDE_LIST_PATH");
			}else{
				isFile = new File(PEPTIDE_LIST_PATH);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, PEPTIDE_LIST_PATH);
				}
			}
			
			if(IS_PROTEIN_DB && PROTEIN_DB_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "PROTEIN_DB_PATH");
			}else if(IS_PROTEIN_DB){
				isFile = new File(PROTEIN_DB_PATH);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, PROTEIN_DB_PATH);
				}				
				if(isFile.isDirectory() && isFile.list().length == 0){
					Error.exitError(Error.IS_EMPTY, PROTEIN_DB_PATH);
				}
			}
			
			if(IS_VSG && GRAPH_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "GRAPH_PATH");
			}else if(IS_VSG){
				isFile = new File(GRAPH_PATH);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, GRAPH_PATH);
				}
			}
			
			if(IS_MUTATION && VCF_PATH == null){
				Error.exitError(Error.PARAM_ERROR, "VCF_PATH");
			}else if(IS_MUTATION){
				isFile = new File(VCF_PATH);
				if(!isFile.exists() ){
					Error.exitError(Error.NO_SUCH_A_FILE, VCF_PATH);
				}
				if(isFile.isDirectory() && isFile.list().length == 0){
					Error.exitError(Error.IS_EMPTY, VCF_PATH);
				}
				
			}
			
			if(OUTPUT_ROOT == null){
				Error.exitError(Error.PARAM_ERROR, "GFF_OUTPUT_PATH");
			}else{
				isFile = new File(OUTPUT_ROOT);
				if(!isFile.exists()){
					Error.exitError(Error.NO_SUCH_A_FILE, "OUTPUT_ROOT");
				}
			}
		}
		
		if(case_ == SEARCH_PHASE) printEnv();
		
		return 0;
	}
	
	public static void printEnv(){
		System.out.println("Environment");
		System.out.println("\tMapping method: "+MAPPING_METHOD);
		System.out.println("\tPeptide list: "+PEPTIDE_LIST_PATH);
		if(IS_I_SAME_WITH_L){
			System.out.println("\tIsoleucine is equivalent to leucine");
		}
		
		if(IS_PROTEIN_DB){
			System.out.println("Protein database mapping");
			System.out.println("\tProtein database: "+PROTEIN_DB_PATH);
			
			if(IS_SAV){
				System.out.println("\tSingle amino-acids variantion");
			}
		}
		
		if(IS_VSG){
			System.out.println("Variant splice graph mapping");
			if(IS_JUNCTION){
				System.out.println("\tJunction variation");
			}
			
			if(IS_ALTERNATIVE){
				System.out.println("\tExon skipping");
			}
			
			if(IS_INTRON){
				System.out.println("\tExon-extension");
			}
			
			if(IS_MUTATION){
				System.out.println("\tSingle nucleotide variation");
				System.out.println("\t\tVCF: "+VCF_PATH);
			}
		}
		
		if(IS_SFT){
			System.out.println("Six-frame translation mapping");
			System.out.println("\tReference genome: "+REFERENCE_GENOME_PATH);
		}
		
	}
}
