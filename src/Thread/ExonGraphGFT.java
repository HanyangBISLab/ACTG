package Thread;

import java.io.File;
import java.io.IOException;

import ExonGraph.ExonGraphGF;

public class ExonGraphGFT extends Thread{
	public static int TheNumberOfSuccess = 0;
	
	public File gtfFile = null;
	public File fastaFile = null;
	public ExonGraphGF EG = null;
	public int cls = 0;
	public boolean isDone = false; 
	
	public ExonGraphGFT(ExonGraphGF EG, File gtfFile, File fastaFile, int cls){
		this.gtfFile = gtfFile;
		this.fastaFile = fastaFile;
		this.cls = cls;
		this.EG = EG;
	}
	
	public void run(){
		try {
			EG = new ExonGraphGF(gtfFile, fastaFile, cls);
			TheNumberOfSuccess++;
			isDone = true;
			System.out.println(fastaFile.getName()+" is done");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
