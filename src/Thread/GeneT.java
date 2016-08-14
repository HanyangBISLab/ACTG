package Thread;

import ExonGraph.Chromosome;
import ExonGraph.GENE;
import ExonGraph.KAminoTree;

public class GeneT extends Thread{
	Chromosome chr = null;
	public boolean visited;
	public KAminoTree tree = null;
	public int startIndex = 0;
	public int endIndex = 0;
	
	public GeneT(Chromosome chr, int startIndex, int endIndex, boolean visited, KAminoTree tree){
		this.chr = chr;
		this.visited = visited;
		this.tree = tree;
		this.startIndex = startIndex;
		this.endIndex = endIndex;
	}
	
	public void run(){
		for(int index = startIndex; index < endIndex; index++){
			chr.GeneArray[index].chrID = chr.get_name();
			chr.GeneArray[index].AminoExonSearch(visited, tree);
		}
		
	}
}
