package sixFrameTranslation;

import java.util.LinkedList;
//import java.util.Hashtable;
import java.util.LinkedList;


public class AhoCorasick {
	
	public state ROOT = null;
	public static final int sizeOfChild = 26;
	public state[] curState = null;
	
	public AhoCorasick() {
		ROOT = new state(null);
		ROOT.failState = ROOT;
		curState = new state[3];
		
		curState[0] = ROOT;
		curState[1] = ROOT;
		curState[2] = ROOT;
	}
	
	public void addPattern(String str, String output){
		state curState = ROOT;
		String pat = str;
		int patLength = pat.length();
		
		
		for(int i=0; i<patLength; i++){
			int characterIndex = pat.charAt(i) - 'A';			
			if(curState.nextState[characterIndex] == null){
				curState.nextState[characterIndex] = new state(String.valueOf(pat.charAt(i)));
			}
			curState = curState.nextState[characterIndex];
		}
		curState.addOutput(output);
	}	
	
	public void buildFailState(){	
		state curState = ROOT;		
		LinkedList<state> queue = new LinkedList<state>();
		
		for(int i=0; i<sizeOfChild; i++){
			if(curState.nextState[i] != null){	
				queue.add(curState.nextState[i]);
				curState.nextState[i].failState = ROOT;
			}
		}
		
		
		while(!queue.isEmpty()){	
			curState = queue.get(0);	
			queue.remove(0);
			
			for(int i=0; i<sizeOfChild; i++){
				if(curState.nextState[i] == null){	
					continue;						
				}
				
				queue.add(curState.nextState[i]);	
				
				if(curState.failState.nextState[i] != null){	
					curState.nextState[i].failState = curState.failState.nextState[i];
				}else{											
					state failState = curState.failState;
					
					while(true){								
										
						if(failState.nextState[i] != null){		
							curState.nextState[i].failState = failState.nextState[i];
							break;
						}else if(failState == ROOT){			
							curState.nextState[i].failState = ROOT;
							break;
						}
						failState = failState.failState;
					}
				}
			}
		}
		
		buildOutput();
	}
	
	public void buildOutput(){
		state curState = ROOT;
		LinkedList<state> queue = new LinkedList<state>();
		
		for(int i=0; i<sizeOfChild; i++){
			if(curState.nextState[i] != null){
				queue.add(curState.nextState[i]);
			}
		}
		
		while(!queue.isEmpty()){	
			curState = queue.get(0);
			queue.remove(0);
			
			for(int i=0; i<sizeOfChild; i++){
				if(curState.nextState[i] == null){
					continue;
				}
				
				queue.add(curState.nextState[i]);
				
				curState.nextState[i].addOutput(curState.nextState[i].failState.getOutput());
			}
		}
	}
	
		
	public LinkedList<String> searchString(char transCodon, int pos){
		LinkedList<String> match = new LinkedList<String>();
		
		int characterIndex = transCodon - 'A';
		if(curState[pos].nextState[characterIndex] != null){
			curState[pos] = curState[pos].nextState[characterIndex];
			for(String t : curState[pos].output){
				match.add(t);
			}
		}
		else{
			while(true){
				curState[pos] = curState[pos].failState;

				if(curState[pos].nextState[characterIndex] != null){
					curState[pos] = curState[pos].nextState[characterIndex];
					for(String t : curState[pos].output){
						match.add(t);
					}
					break;
				}
				else if(curState[pos] == ROOT){
					curState[pos] = ROOT;
					break;
				}
			}
		}
		return match;
	}
}

class state {
	public String character = null;
	public state[] nextState = null;
	public state failState = null;
	public LinkedList<String> output = null;
	
	public state(String character){
		this.nextState = new state[AhoCorasick.sizeOfChild];
		for(int i=0; i<AhoCorasick.sizeOfChild; i++){
			nextState[i] = null;
		}

		this.character = character;
		this.output = new LinkedList<String>();
	}
	
	public void addOutput(String output){
		this.output.add(output);
	}
	
	public void addOutput(LinkedList<String> outputList){
		for(int i=0; i<outputList.size(); i++){
			this.output.add(outputList.get(i));
		}
	}
	
	public LinkedList<String> getOutput(){
		return this.output;
	}
}
