package ExonGraph;

import java.io.PrintWriter;
import java.util.HashSet;

public class ExonGraph {
	static public boolean oneMode = true;
	static public int JVALUE = 15; // The maximum number of junction variables
	/**
	 * 0~3: JNC
	 * 4~6: SNV
	 * 7~9: INS
	 * 10~12: DEL
	 * 13: INTRON
	 * 14: ALT
	 * 
	 */
	static public int CHR_NUM = 24;
	static public int CLS; // Parameter or Static ??
	private boolean initialize;
	//RefSeqExonGraph Addition
	static protected Chromosome chromosome[];
	//END
	//private boolean nowStatus = true;
	static public boolean nowStatus = false;
	static public int NodeNumber = 0;
	
	static public PrintWriter err;
	static public PrintWriter err1;
	int XXX = 27;
	
	static public String Chrom;
	static public String GeneID;
	
	static public HashSet<String> BioType = new HashSet<String>();
	
	
	//RefSeqExonGraph Addition
	public ExonGraph(){
		if(this.getClass() == ExonGraphGF.class){
			if(chromosome == null)
				chromosome = new Chromosome[CHR_NUM];
			
		}else
			chromosome = new Chromosome[CHR_NUM];
	}
	//END
	
	
	public boolean get_initialize() { return initialize; }
	
	
	public Chromosome[] getChromosome() {
		return chromosome;
	}


	
	public void AnimoExonSearch(boolean visited, KAminoTree tree) {
		//statistics
		int totalGenes = 0;
		for(int i=0; i<CHR_NUM; i++){
			totalGenes += chromosome[i].get_genecnt();
		}
		
		totalGenes /= 100;
		
		if(totalGenes == 0){
			totalGenes = 1;
		}
		
		int count = 0;
		int ratio = 0;
		int interval = 5;
		for(int i=0; i<CHR_NUM; i++){
			
			for(int j=0;j<chromosome[i].get_genecnt();j++){
				Chrom = chromosome[i].get_name(); GeneID = new String(chromosome[i].get_gene(j).get_gene_id());
				
				//if(!GeneID.equalsIgnoreCase("ENSG00000124783")) continue;
				
				chromosome[i].get_gene(j).AminoExonSearch(visited, tree);
				count ++;
				if(ratio + interval <= count / totalGenes){
					ratio = count/totalGenes;
					System.out.println(ratio + " %");
				}
			}
		}
	}
	
}