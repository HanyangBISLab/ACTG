package ExonGraph;

import java.util.ArrayList;

public class KAminoTree {
	
	public Pattern Set;
	private int PSize;
	public KAminoNode Tree;
	public KAminoNode RTree;
	public int sigma;
	
	public KAminoTree(int size, int sig) {
		Set = null;
		Tree = new KAminoNode(sig);
		Tree.add(new KAminoNode(sig));
		RTree = new KAminoNode(sig);
		RTree.add(new KAminoNode(sig));
		PSize = size;
		sigma = sig;
	}
	
	public void BuildKeyword(Pattern a, KAminoNode Tree)
	{
		int i;
		int now_state;

		//System.out.println(a.get_amino());
		now_state = 0;
		int StringSize = a.get_amino().length();
		for(i=0;i<StringSize;i++) {
			if(Tree.get(now_state).get_next(a.get_amino().charAt(i)) == -1) { // 현재 노드에 해당 엣지가 존재하지 않은 경우
				Tree.add(new KAminoNode(sigma));
				int state_cnt = Tree.size() -1;
				Tree.get(now_state).set_next(a.get_amino().charAt(i), state_cnt);
				Tree.get(state_cnt).set_depth(Tree.get(now_state).get_depth() + 1);
				now_state = state_cnt;
			}
			else {	// 현재 노드에 해당 엣지가 존재하는 경우
				now_state = Tree.get(now_state).get_next(a.get_amino().charAt(i));
			}
			if(i == StringSize-1) Tree.get(now_state).set_matchlist(a);
		}
	}
	
	public KAminoNode get_tree() { return this.Tree; }
	public KAminoNode get_rtree() { return this.RTree; }
	
	public void BuildFlink() {
		BuildFlink(Tree);
		BuildFlink(RTree);
	}
	public void BuildOlink() { BuildOlink(Tree);BuildOlink(RTree); }
	
	public void BuildFlink(KAminoNode Tree) {
		ArrayList<KAminoNode> queue;
		KAminoNode tmp;
		int i, flink;

		queue = new ArrayList<KAminoNode>();
		
		
		queue.add(Tree.get(0));
		tmp = queue.get(0);
		queue.remove(0);
		
		tmp.set_failstate(0);
		
		for(i=0;i<this.sigma;i++) {
			if(tmp.get_next(i) != -1) {
				Tree.get(tmp.get_next(i)).set_failstate(0);
				queue.add(Tree.get(tmp.get_next(i)));
			}
		}
		
		while(queue.size() != 0) {
			tmp = queue.get(0);
			queue.remove(0);

			for(i=0;i<this.sigma;i++) {
				if(tmp.get_next(i) != -1) {
					queue.add(Tree.get(tmp.get_next(i)));
					flink = tmp.get_failstate();
					while(flink != 0) {
						if(Tree.get(flink).get_next(i) != -1) break;
						flink = Tree.get(flink).get_failstate();
					}
					if(Tree.get(flink).get_next(i) != -1) Tree.get(tmp.get_next(i)).set_failstate(Tree.get(flink).get_next(i));
					else Tree.get(tmp.get_next(i)).set_failstate(flink);
				}
			}
		}
	}
	
	public void BuildOlink(KAminoNode Tree) {
		ArrayList<KAminoNode> queue;
		KAminoNode tmp;
		int i, flink;

		queue = new ArrayList<KAminoNode>();
		
		
		queue.add(Tree.get(0));
		tmp = queue.get(0);
		queue.remove(0);
		
		for(i=0;i<this.sigma;i++) {
			if(tmp.get_next(i) != -1) {
				queue.add(Tree.get(tmp.get_next(i)));
			}
		}
		
		while(queue.size() != 0) {
			tmp = queue.get(0);
			queue.remove(0);

			for(i=0;i<this.sigma;i++) {
				if(tmp.get_next(i) != -1) {
					queue.add(Tree.get(tmp.get_next(i)));
				}
			}
			flink = tmp.get_failstate();
			if(Tree.get(flink).MatchList != null) {
				if(tmp.MatchList != null) {
					tmp.MatchList.add(Tree.get(flink).MatchList);
					tmp.MatchList.addAll(Tree.get(flink).MatchList);
				}
				else {
					tmp.MatchList = Tree.get(flink).MatchList;
				}
			}
		}
	}

	public void NFAtoDFA() {
		NFAtoDFA(Tree);
		NFAtoDFA(RTree);
	}
	
	public void NFAtoDFA(KAminoNode Tree) {
		ArrayList<KAminoNode> queue;
		KAminoNode tmp;
		int i;
	
		queue = new ArrayList<KAminoNode>();
		
		queue.add(Tree.get(0));
		tmp = queue.get(0);
		queue.remove(0);

		for(i=0;i<this.sigma;i++) {
			if(tmp.get_next(i) != -1) queue.add(Tree.get(tmp.get_next(i)));
			else tmp.set_next(i, 0);
		}
	
		while(queue.size() != 0) {
			tmp = queue.get(0);
			queue.remove(0);
			
			for(i=0;i<this.sigma;i++) {
				if(tmp.get_next(i) == -1) tmp.set_next(i, Tree.get(tmp.get_failstate()).get_next(i));
				else queue.add(Tree.get(tmp.get_next(i)));
			}
		}
	}
}
