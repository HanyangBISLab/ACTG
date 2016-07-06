package ExonGraph;

import java.io.Serializable;
import java.util.ArrayList;


public class Chromosome implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 7212090891662502901L;
	private int Gene_cnt;
	private String name;	
	private GENE[] GeneArray;
	
	public Chromosome(String readString) {
		this.name = readString;
		Gene_cnt = 0;
		GeneArray = null;
	
	}

	//TODO:: Chromosome
	//RefSeqExonGraph Addition
	public Chromosome(ArrayList<GENE> genes, String name){
		this.Gene_cnt = genes.size();
		this.name = name;
		
		GeneArray = new GENE[Gene_cnt];
		
		for(int i=0; i<Gene_cnt; i++){
			GeneArray[i] = genes.get(i);
		}
		
	}
	//END
	
	public Chromosome(String name, String geneCnt){
		this.Gene_cnt = Integer.parseInt(geneCnt);
		this.name = name;
		
		GeneArray = new GENE[this.Gene_cnt];
	}

	public GENE get_gene(int i) { return GeneArray[i]; }
	public String get_name() { return name;}
	public int get_genecnt() { return Gene_cnt; }

	public void setGene(int i, GENE gene) { this.GeneArray[i] = gene; }
}
