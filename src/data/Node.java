package data;

import java.util.ArrayList;
import java.util.Hashtable;

public class Node {

	public char nt;
	public ArrayList<Integer> outputs = new ArrayList<Integer>();
	public Node parent;
	public Node failNode;
	public Hashtable<Character, Node> nexts = new Hashtable<Character, Node>();
	
}
