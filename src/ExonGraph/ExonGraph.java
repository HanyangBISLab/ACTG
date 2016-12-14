package ExonGraph;

import java.io.PrintWriter;
import java.util.HashSet;

import Environments.Constants;
import Thread.GeneT;

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
	
	static public HashSet<String> BioType = new HashSet<String>();
	
	
	public ExonGraph(){
		if(this.getClass() == ExonGraphGF.class){
			if(chromosome == null)
				chromosome = new Chromosome[CHR_NUM];
			
		}else
			chromosome = new Chromosome[CHR_NUM];
	}
	
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
		int geneDefaultInterval = 500;
		boolean complete = false;
		for(int i=0; i<CHR_NUM; i++){
			int geneCnt = chromosome[i].get_genecnt();
			for(int j=0;j<geneCnt;){
				
				if(Thread.activeCount() >= Constants.THE_NUMBER_OF_THREADS){
					chromosome[i].GeneArray[j].chrID = chromosome[i].get_name();
					chromosome[i].GeneArray[j].AminoExonSearch(visited, tree);
					count ++;
					j++;
				}else{
					int startIndex = j;
					int endIndex = j+geneDefaultInterval;
					if(endIndex > geneCnt){
						endIndex = geneCnt;
					}
					
					GeneT geneT = new GeneT(chromosome[i], startIndex, endIndex, visited, tree);
					geneT.start();
					j = endIndex;
					count+= endIndex - startIndex;
				}
				
				if(ratio + interval <= count / totalGenes){
					ratio = count/totalGenes;
					
					
					if(ratio > 100){
						System.out.println("100 %");
						complete = true;
					}else{
						System.out.println(ratio + " %");
						if(ratio == 100){
							complete = true;
						}
					}
					
				}
			}
		}
		
		while(Thread.activeCount() != 1){
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		if(!complete){
			System.out.println("100 %");
		}
	}
	
}