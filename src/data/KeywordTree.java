package data;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

public class KeywordTree {

	public Node root = null;
	public ArrayList<Pattern> patterns = null;
	
	public KeywordTree (ArrayList<Pattern> patterns) {
		root = new Node();
		root.nt = '-';
		
		this.patterns = patterns;
		System.out.println("Building prefix-trie...");
		int patternIndex = 0;
		for(Pattern pattern : patterns) {
			Node curNode = root;
			int len = pattern.aaSequence.length();
			for(int i=0; i<len; i++) {
				char nt = pattern.aaSequence.charAt(i);
				if(nt == 'I' || nt == 'L') nt = 'J';
				
				if(curNode.nexts.get(nt) != null) {
					curNode = curNode.nexts.get(nt);
				} else {
					Node newNode = new Node();
					newNode.nt = nt;
					curNode.nexts.put(nt, newNode);
					newNode.parent = curNode;
					curNode = newNode;
				}
				
				if(i == len-1) {
					curNode.outputs.add(patternIndex);
				}
			}
			patternIndex++;
		}
		
		buildFailureLink();
		buildOutputLink();
		
	}
	
	public void buildFailureLink () {
		LinkedList<Node> queue = new LinkedList<Node>();
		queue.add(root);
		
		// build failure links
		root.failNode = root;
		root.parent = root;
		Node curNode = queue.poll();
		System.out.println("Building failure link...");
		while(curNode != null) {
			Node pNode = curNode.parent;
			if(pNode == root) {
				curNode.failNode = root;
			} else {
				boolean isConnected = false;
				while(pNode != root) {
					pNode = pNode.failNode;
					Node fNode = pNode.nexts.get(curNode.nt);
					if(fNode != null) {
						curNode.failNode = fNode;
						isConnected = true;
						break;
					}
				}
				if(!isConnected) curNode.failNode = root;
			}
			
			if(curNode.nexts.size() != 0) {
				Iterator<Character> nts = (Iterator<Character>) curNode.nexts.keys();
				while(nts.hasNext()) {
					queue.add(curNode.nexts.get(nts.next()));
				}
			}
			curNode = queue.poll();
		}
	}
	
	public void buildOutputLink () {
		System.out.println("Building output link...");
		LinkedList<Node> queue = new LinkedList<Node>();
		queue.add(root);
		Node curNode = queue.poll();
		while(curNode != null) {
			Node fNode = curNode.failNode;
			if(fNode.outputs.size() != 0) {
				curNode.outputs.addAll(fNode.outputs);
			}
			
			if(curNode.nexts.size() != 0) {
				Iterator<Character> nts = (Iterator<Character>) curNode.nexts.keys();
				while(nts.hasNext()) {
					queue.add(curNode.nexts.get(nts.next()));
				}
			}
			curNode = queue.poll();
		}
	}
}
