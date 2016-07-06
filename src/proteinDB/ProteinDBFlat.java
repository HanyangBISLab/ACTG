package proteinDB;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.LinkedList;

public class ProteinDBFlat {
	public Hashtable<String, LinkedList<String>> headerMapping;
	public LinkedList<String> seqList;
	
	public ProteinDBFlat(){
		headerMapping = new Hashtable<String, LinkedList<String>>();
		seqList = new LinkedList<String>();
	}
	
	public void putMapping(String header, String pepSeq){
		LinkedList<String> mappedList = headerMapping.get(pepSeq);
		
		if(mappedList == null){
			mappedList = new LinkedList<String>();
			seqList.add(pepSeq);
		}
		
		mappedList.add(header);
		headerMapping.put(pepSeq, mappedList);
	}
	
	public void writeMapping(String fileName){
		FileWriter FW = null;
		try {
			FW = new FileWriter(fileName);
			BufferedWriter BW = new BufferedWriter(FW);
			
			LinkedList<String> headers = null;
			for(String pepSeq : seqList){
				
				headers = headerMapping.get(pepSeq);
				BW.append(pepSeq).append("\n");
				for(String header : headers){
					BW.append(header).append("\n");
				}
				BW.append("\n");
				
			}
			
			BW.close(); FW.close();
		} catch (IOException e) {}
		
	}
}
