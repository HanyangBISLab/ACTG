package data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

import Environments.Codon;

public class Transcript {

	public ArrayList<Region> regions = null;
	public String transcriptID = null;
	public String geneID = null;
	public String geneName = null;
	public String transcriptType = null;
	public String chr = null;
	public boolean strand = true;
	
	
	public Transcript (String chr, String transcriptID, String geneID, String geneName, String transcriptType, boolean strand) {
		this.regions = new ArrayList<Region>();
		this.transcriptID = transcriptID;
		this.geneID = geneID;
		this.geneName = geneName;
		this.transcriptType = transcriptType;
		this.strand = strand;
	}
	
	/**
	 * Define CDS/UTR/NONCODING classes
	 * 
	 */
	public void refineTranscript () {
		Collections.sort(this.regions);
		int sizeOfRegions = regions.size();
		boolean isCDS = false;
		ArrayList<Region> exons = new ArrayList<Region>();
		ArrayList<Region> cdss = new ArrayList<Region>();
		for(int i=0; i<sizeOfRegions; i++) {
			if(this.regions.get(i).regionClass == Region.CDS) {
				isCDS = true;
				cdss.add(this.regions.get(i));
			} else {
				exons.add(this.regions.get(i));
			}
		}
		
		// if there is no CDS class, then whole exons are defined as non-coding region
		if(!isCDS) {
			for(int i=0; i<sizeOfRegions; i++) {
				this.regions.get(i).regionClass = Region.NONCODING;
			}
		} else {
			// UTR and CDS should be separated properly.
			ArrayList<Region> newRegion = new ArrayList<Region>();
			
			int exonIndex = 0;
			int cdsIndex = 0;
			while(exonIndex != exons.size() && cdsIndex != cdss.size()) {
				Region exon = exons.get(exonIndex);
				Region cds = cdss.get(cdsIndex);
				
				// inclusive exon and cds
				if(exon.start <= cds.start && exon.end >= cds.end) {
					if(exon.start == cds.start && exon.end == cds.end) newRegion.add(cds);
					else {
						if(exon.start == cds.start) {
							newRegion.add(cds);
							exon.start = cds.end+1;
							newRegion.add(exon);
						} else if(exon.end == cds.end) {
							exon.end = cds.start-1;
							newRegion.add(exon);
							newRegion.add(cds);
						} else {
							System.out.println("Region Error was detected!");
						}
					}
					exonIndex++;
					cdsIndex++;
				} else {
					newRegion.add(exon);
					exonIndex++;
				}
			}
			
			while(exonIndex != exons.size()) {
				Region exon = exons.get(exonIndex);
				newRegion.add(exon);
				exonIndex++;
			}
			
			this.regions = newRegion;
			Collections.sort(this.regions);
			
			// 5'-UTR and 3'-UTR?
			boolean turn = false;
			for(int i=0; i<this.regions.size(); i++) {
				if(this.regions.get(i).regionClass == Region.UNDEFINED) {
					if(this.strand) {
						if(turn) {
							this.regions.get(i).regionClass = Region.UTR3;
						} else {
							this.regions.get(i).regionClass = Region.UTR5;
						}
					} else {
						if(turn) {
							this.regions.get(i).regionClass = Region.UTR5;
						} else {
							this.regions.get(i).regionClass = Region.UTR3;
						}
					}
				} else {
					turn = true;
				}
			}
		}
		
		// Introns
		ArrayList<Region> newRegion = new ArrayList<Region>();
		newRegion.add(this.regions.get(0));
		for(int i=0; i<this.regions.size()-1; i++) {
			Region lRegion = this.regions.get(i);
			Region rRegion = this.regions.get(i+1);
			
			// is not adjacent?
			
			if(lRegion.end+1 != rRegion.start) {
				Region intron = new Region(lRegion.end+1, rRegion.start-1);
				intron.regionClass = Region.INTRON;
				newRegion.add(intron);
			}
			newRegion.add(rRegion);
		}
		Collections.sort(newRegion);
		this.regions = newRegion;
		
		// linear DAG
		for(int i=0; i<this.regions.size()-1; i++) {
			Region lRegion = this.regions.get(i);
			Region rRegion = this.regions.get(i+1);
			lRegion.nextNormal = rRegion;
			rRegion.prevNormal = lRegion;
		}
	}
	
	public void setExonSkipping () {
		for(int i=0; i<this.regions.size()-1; i++) {
			Region lRegion = this.regions.get(i);
			for(int j=i+1; j<this.regions.size(); j++) {
				Region rRegion = this.regions.get(j);
				// CDS and Noncoding only
				if(lRegion.regionClass == Region.CDS && rRegion.regionClass == Region.CDS) {
					lRegion.nextExonSkipping = rRegion;
					rRegion.prevExonSkipping = lRegion;
					break;
				} else if(lRegion.regionClass == Region.NONCODING && rRegion.regionClass == Region.NONCODING) {
					lRegion.nextExonSkipping = rRegion;
					rRegion.prevExonSkipping = lRegion;
					break;
				}
			}
		}
	}
	
	public void setVariants (Variants variants) {
		
		//TODO variants..!
		for(int i=0; i<this.regions.size(); i++) {
			Region region = this.regions.get(i);
			// sorted variants
			ArrayList<Variant> selectedVariants = variants.getVariants(chr, region.start, region.end);
			
			ArrayList<Region> varRegions = new ArrayList<Region>();
			for(Variant variant : selectedVariants) {
				if(region.start == variant.pos) {
					
				} else if(region.start < variant.pos && variant.pos < region.end) {
					
				} else if(region.end == variant.pos) {
					
				}
			}
			
		}
		
	}
	
	public int find (KeywordTree keywordTree) {
		int count = 0;
		int size = this.regions.size();
		StringBuilder nts = new StringBuilder();
		
		Node curNode = keywordTree.root;
		LinkedList<Region> path = new LinkedList<Region>();
		path.add(this.regions.get(0)); // stack
		if(this.strand) {
			// TODO!
			while(!path.isEmpty()) {
				Region region = path.pop();
				
			}
			
			for(int i=0; i<size; i++) {
				Region region = this.regions.get(i);
				if(region.regionClass == Region.CDS) {
					int len = region.nucleotides.length();
					for(int j=0; j<len; j++) {
						nts.append(region.nucleotides.charAt(j));
						
						if(nts.length() == 3) {
							char aa = Codon.NuclToAmino(nts.toString());
							if(aa == 'I' || aa == 'L') aa = 'J';
							while(curNode.nexts.get(aa) == null) {
								if(curNode == keywordTree.root) {
									break;
								}
								curNode = curNode.failNode;
							}
							
							if(curNode.nexts.get(aa) != null) curNode = curNode.nexts.get(aa);
							if(curNode.outputs.size() != 0) {
								System.out.println("FIND-");
								for(int patternIndex : curNode.outputs) {
									System.out.println(keywordTree.patterns.get(patternIndex).aaSequence +"\t" +this.transcriptID+"\t"+this.geneID+"\t"+this.transcriptType);
								}
							}
							
							nts.setLength(0);
						}
					}
				}
			}
		} else {
			
		}
		
		return count;
	}
}
