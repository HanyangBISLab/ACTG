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
	public GENE[] GeneArray;
	//TODO:: Chromosome
	public Chromosome(ArrayList<GENE> genes, String name){
		this.Gene_cnt = genes.size();
		this.name = name;
		
		GeneArray = new GENE[Gene_cnt];
		
		for(int i=0; i<Gene_cnt; i++){
			GeneArray[i] = genes.get(i);
		}
	}
	
	public String get_name() { return name;}
	public int get_genecnt() { return Gene_cnt; }

}
