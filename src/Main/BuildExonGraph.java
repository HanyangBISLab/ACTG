package Main;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.OutputStream;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import Environments.Constants;
import ExonGraph.ExonGraphGF;
import Thread.ExonGraphGFT;

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
		// fasta가 들어있는 폴더
		
		// LargeScale Deletion

		File[] fastaFileList = new File(Constants.REFERENCE_GENOME_PATH).listFiles();
		// gtf가 들어있는 폴더
		// LargeScale Deletion

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
		
		//!! fasta 파일 이름(확장자제외)과 gtf에 적혀있는 chrName을 동일하게 할 것
		//   ex> 1.fa 이면  gtf 파일안에 적혀있는 chrName은 1
		//       X.fa 이면 gtf 파일안에 적혀있는 chrName은 X
				
		ExonGraphGF.CHR_NUM = fastaFileList.length;
		TheNumberOfTasks = fastaFileList.length;
		ExonGraphGF EG = null;
		
		ExonGraphGFT[] exonGraphGFTList = new ExonGraphGFT[ExonGraphGFT.TheNumberOfThreads-1];
		if(gtfFileList.length == fastaFileList.length){
			// GTF 파일 수와 Fasta 파일 수가 같을 경우
			// 즉, 1.fa  -- 1.gtf
			//    2.fa  -- 2.gtf
			//    ....
			//    10.fa  -- 10.gtf
			if(fastaFileList.length == 1){
				EG = new ExonGraphGF(gtfFileList[0], fastaFileList[0], -1);
				ExonGraphGFT.TheNumberOfSuccess++;
				System.out.println(fastaFileList[0].getName()+" is done");
			}
			
			for(int fIndex=0; fIndex < fastaFileList.length-1; fIndex++){
				if(fIndex == 0){
					EG = new ExonGraphGF(gtfFileList[fastaFileList.length-1], fastaFileList[fastaFileList.length-1], -1);
					ExonGraphGFT.TheNumberOfSuccess++;
					System.out.println(fastaFileList[fastaFileList.length-1].getName()+" is done");
				}
				
				if(Thread.activeCount() != ExonGraphGFT.TheNumberOfThreads){
					
					for(int i=0; i<ExonGraphGFT.TheNumberOfThreads-1; i++){
						if(exonGraphGFTList[i] == null || exonGraphGFTList[i].isDone){
							exonGraphGFTList[i] = new ExonGraphGFT(EG, gtfFileList[fIndex], fastaFileList[fIndex], -1);
							exonGraphGFTList[i].start();
							break;
						}
					}
					
				}else{
					EG = new ExonGraphGF(gtfFileList[fIndex], fastaFileList[fIndex], -1);
					ExonGraphGFT.TheNumberOfSuccess++;
					System.out.println(fastaFileList[fIndex].getName()+" is done");
				}
				
			}
		}else if(gtfFileList.length == 1){
			// 하나의 GTF 파일
			for(int fIndex=0; fIndex < fastaFileList.length-1; fIndex++){
				if(fIndex == 0){
					EG = new ExonGraphGF(gtfFileList[0], fastaFileList[fastaFileList.length-1], -1);
					ExonGraphGFT.TheNumberOfSuccess++;
					System.out.println(fastaFileList[fastaFileList.length-1].getName()+" is done");
				}
				
				if(Thread.activeCount() != ExonGraphGFT.TheNumberOfThreads){
					
					for(int i=0; i<ExonGraphGFT.TheNumberOfThreads-1; i++){
						if(exonGraphGFTList[i] == null || exonGraphGFTList[i].isDone){
							exonGraphGFTList[i] = new ExonGraphGFT(EG, gtfFileList[0], fastaFileList[fIndex], -1);
							exonGraphGFTList[i].start();
							break;
						}
					}
					
				}else{
					EG = new ExonGraphGF(gtfFileList[0], fastaFileList[fIndex], -1);
					ExonGraphGFT.TheNumberOfSuccess++;
					System.out.println(fastaFileList[fIndex].getName()+" is done");
				}
			}
		}
		
		
		for(int i=0; i<exonGraphGFTList.length; i++){
			if(exonGraphGFTList[i] != null){
				if(exonGraphGFTList[i].getState() != Thread.State.TERMINATED){
					exonGraphGFTList[i].join();
				}
			}
			
		}

		createGraphFile(EG, Constants.GRAPH_OUTPUT_PATH);
		System.out.println("Construction Time for Exon Graph : " + (System.currentTimeMillis()-startTime)/1000 + " Sec" );

	
	}
	
	public static void createGraphFile(ExonGraphGF EG, String fileName) throws IOException{
		OutputStream outFile = new FileOutputStream(fileName);
		OutputStream buffer = new BufferedOutputStream(outFile);
		ObjectOutput output = new ObjectOutputStream(buffer);
		
		output.writeObject(EG.getChromosome());
		output.close();
		outFile.close();
		
	}
	
}
