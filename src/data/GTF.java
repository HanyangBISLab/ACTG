package data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

public class GTF {

	public Hashtable<String, ArrayList<Transcript>> transcriptsPerChr = new Hashtable<String, ArrayList<Transcript>>();
	
	public GTF (File file) throws IOException {
		readFile(file);
	}
	
	public void find (KeywordTree keywordTree) {
		if(transcriptsPerChr.size() == 0) return;
		Iterator<String> chrs = (Iterator<String>) transcriptsPerChr.keys();
		while(chrs.hasNext()) {
			String chr = chrs.next();
			ArrayList<Transcript> transcripts = transcriptsPerChr.get(chr);
			
			for(Transcript transcript : transcripts) {
				transcript.find(keywordTree);
			}
		}
	}
	
	public void setSequence (File[] fastaFiles) throws IOException {
		for(File file : fastaFiles) {
			String chr = file.getName().split("\\.")[0];
			ArrayList<Transcript> transcripts = transcriptsPerChr.get(chr);
			StringBuilder sequences = new StringBuilder();
			if(transcripts != null) {
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;
				
				while((line = BR.readLine()) != null) {
					if(line.startsWith(">")) continue;
					sequences.append(line);
				}
				
				BR.close();
				
				for(Transcript transcript : transcripts) {
					for(Region region : transcript.regions) {
						region.nucleotides = new StringBuilder(sequences.substring(region.start-1, region.end).toUpperCase());
						if(!transcript.strand) {
							int len = region.nucleotides.length();
							for(int i=0; i<len; i++) {
								switch(region.nucleotides.charAt(i)) {
								case 'A': region.nucleotides.setCharAt(i, 'T'); break;
								case 'C': region.nucleotides.setCharAt(i, 'G'); break;
								case 'G': region.nucleotides.setCharAt(i, 'C'); break;
								case 'T': region.nucleotides.setCharAt(i, 'A'); break;
								default: break;
								}
							}
						}
					}
				}
			}
			
			
		}
	}
	
	private void readFile (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Hashtable<String, Transcript> regionMap = new Hashtable<String, Transcript>();
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith("#")) continue; // skip header
			
			String[] fields = line.split("\t");
			
			String chr = fields[0];
			String type = fields[2];
			if(type.equalsIgnoreCase("exon") || type.equalsIgnoreCase("cds")) {
				String[] attr = fields[8].split(";");
				String transcriptID = getGtfAttr(attr, "transcript_id");
				String key = chr+"_"+transcriptID;
				
				// set locus information
				boolean strand = fields[6].equalsIgnoreCase("+") ? true : false;
				int start = Integer.parseInt(fields[3]);
				int end = Integer.parseInt(fields[4]);
				String geneID = getGtfAttr(attr, "gene_id");
				String geneName = getGtfAttr(attr, "gene_name");
				String transcriptType = getGtfAttr(attr, "transcript_type");
				
				Region region = new Region(start, end);
				if(type.equalsIgnoreCase("cds")) region.regionClass = Region.CDS;
				
				Transcript transcript = regionMap.get(key);
				if(transcript == null) transcript = new Transcript(chr, transcriptID, geneID, geneName, transcriptType, strand);
				transcript.regions.add(region);
				regionMap.put(key, transcript);
			}
		}
		
		BR.close();
		
		Iterator<String> chr_transcriptIDs = (Iterator<String>) regionMap.keys();
		while(chr_transcriptIDs.hasNext()) {
			String key = chr_transcriptIDs.next();
			String[] chr_transcriptID = key.split("_");
			String chr = chr_transcriptID[0];
			
			Transcript transcript = regionMap.get(key);
			transcript.refineTranscript(); // define regional information such as CDS/ UTR/ NonCoding /Intron ..
			
			ArrayList<Transcript> transcripts = this.transcriptsPerChr.get(chr);
			if(transcripts == null) transcripts = new ArrayList<Transcript>();
			transcripts.add(transcript);
			this.transcriptsPerChr.put(chr, transcripts);
		}
	}
	
	/**
	 * Convert gtf attribute
	 * 
	 * @param attr
	 * @param tag
	 * @return
	 */
	public String getGtfAttr(String[] attr, String tag){
		
		for(String _s : attr){
			if(_s.contains(tag)){
				return _s.replaceAll("[\"\\s]|"+tag, "");
			}
		}
		
		return null;
	}
}
