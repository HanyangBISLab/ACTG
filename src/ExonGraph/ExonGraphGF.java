package ExonGraph;
import java.io.BufferedInputStream;
/*
 * ExonGraphGF | ver 1.0.0
 * Writer: progistar
 * 
 * Description:
 * ExonGraph瑜� GTF�뙆�씪怨� FASTA�뙆�씪�쓣 �씠�슜�븯�뿬 �깮�꽦
 * Fasta�뙆�씪�쓽 �씠由꾩� chrNum.fa �삎�떇
 * ex>
 * 1.fa
 * 2.fa
 * ...
 * 22.fa
 * x.fa
 * y.fa
 * 
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.logging.Level;
import java.util.logging.Logger;

import Environments.Constants;


public class ExonGraphGF extends ExonGraph implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -996012475149260761L;
	//
	private RandomAccessFile RAFFasta = null;
	private int fastaDescLength = 0;
	private int whiteSpaceInterval = 0;
	
	// �뿼�깋泥� 媛쒖닔瑜� ���옣
	// ExonGraphGF �겢�옒�뒪 媛� 怨듭쑀
	static private int CHR_COUNT = 0;
	private int CUR_CHR_COUNT = 0;
	
	/*
	 * gtfFileLineInfo
	 * 0: seqid (Chr ID)
	 * 1: source
	 * 2: type (Gene, Exon)
	 * 3: start
	 * 4: end
	 * 5: score
	 * 6: strand
	 * 7: phase
	 * 8: attributes (ID)
	 */
	
	// �쑀�쟾�옄瑜� ArrayList濡� ���옣
	private ArrayList<GENE> genes = null;
	// �듃�옖�뒪�겕由쏀듃瑜� ArrayList濡� ���옣
	private ArrayList<EXON> trans = null;
	
	// Junction Variation�쓣 怨좊젮踰붿쐞
	// Exon�쓽 �븵, �뮘 7 媛쒖쓽 臾몄옄瑜� 怨좊젮
	private final int HEAD_LENGTH = 7;
	private final int TAIL_LENGTH = 7;
	
	public ExonGraphGF (String GF) throws IOException, ClassNotFoundException{
		InputStream file = new FileInputStream(GF);
		InputStream buffer = new BufferedInputStream(file);
		ObjectInput input = new ObjectInputStream (buffer);
		
		
		System.out.println("Initializing variant splice graph");
		
		// Inits
		chromosome = (Chromosome[])input.readObject();
		CHR_COUNT = chromosome.length;
		ExonGraph.CHR_NUM = CHR_COUNT;
		
		input.close();
		file.close();
		
		// Mutation 泥섎━
		
		System.out.println("The Number of Chromosomes: "+CHR_COUNT);
		
		if(Constants.IS_MUTATION){
			for(int index=0; index<CHR_COUNT; index++){
				Chromosome curChr = chromosome[index];
								
				String CHR_IDENTIFIER = curChr.get_name().replace("chr", "");
				AddMutation addMutation = new AddMutation(CHR_IDENTIFIER, Constants.VCF_PATH);
				
				int geneMaxIndex = curChr.get_genecnt();
				
				for(int geneIndex=0; geneIndex<geneMaxIndex; geneIndex++){
					addMutation.addToGene(curChr.GeneArray[geneIndex]);
				}
			}
		}
		
		
		System.out.println("The variant splice graph is successfully loaded");
	}

	public ExonGraphGF (File gtfFile, File fastaFile, int cls) throws IOException{
		// �뿼�깋泥� 媛쒖닔 利앷�
		CHR_COUNT++;
		CUR_CHR_COUNT = CHR_COUNT;
		// ExonGraph �깮�꽦�떆,
		// Junction Variation, Alternative Splicing, 洹몃━怨�
		// 硫붾え由� 理쒖쟻�솕瑜� �쐞�븳 �샃�뀡
		CLS = cls;
		
		File fasta = fastaFile;
		File gtf = gtfFile;
		
		// �빐�떦 �뙆�씪�씠 �뾾�쓣 寃쎌슦, ExonGraph瑜� �깮�꽦�븯吏� 紐삵븿
		if(!fasta.exists() || !gtf.exists()){
			System.out.println("File is not existed!");
			return;
		}
		
		genes = new ArrayList<GENE>();
		
		// �듃�옖�뒪�겕由쏀듃 媛쒖닔 珥덇린�솕
		int trans_cnt = 1;
		
		// 
		String prevGeneID = null;
		
		// �뿼�깋泥� �떇蹂꾩옄
		// Fasta File �씠由꾩씠 �떇蹂꾩옄濡� �궗�슜
		String CHR_IDENTIFIER = fastaFile.getName().substring(0, fastaFile.getName().lastIndexOf("."));
		
		GENE gene = null;
		EXON exon = null;
		EXON prevExon_ = null;
		
		String transcript_id = null;
		
		FileReader FRGtf = new FileReader(gtf);		
		BufferedReader BRGtf = new BufferedReader(FRGtf);
		RAFFasta = new RandomAccessFile(fastaFile, "r");
		
		// > �뵒�뒪�겕由쎌뀡 湲몄씠瑜� 援ы븿
		DescLength:
			for(;;){
				byte[] readBuf = new byte[256];
				int readLine = RAFFasta.read(readBuf);
				
				fastaDescLength = new String(readBuf).indexOf('\n');
				
				if(fastaDescLength != -1 || readLine == -1){
					RAFFasta.seek(0);
					// CRLF
					fastaDescLength+=1;
					break DescLength;
				}
			}
		
		// whitespace�쓽 媛꾧꺽�쓣 援ы븿
		// whitespace媛� 256byte瑜� �꽆�뼱 媛� 寃쎌슦, bufCount瑜� 怨꾩궛�븯�뿬 whitespaceInterval�쓣 援ы븿
		boolean isSet = false;
		int bufCount = 0;
		WhiteSpaceInterval:
			for(;;){
				if(!isSet)
					RAFFasta.seek(fastaDescLength);
				isSet = true;
				byte[] readBuf = new byte[256];
				int readLine = RAFFasta.read(readBuf);
				whiteSpaceInterval = (new String(readBuf).indexOf('\n'));
				
				
				if(whiteSpaceInterval != -1 || readLine == -1){
					whiteSpaceInterval += bufCount * 256;
					break WhiteSpaceInterval;
				}else{
					bufCount ++;
				}
			}
		
		boolean GeneStrand = true;
		MakeExonGraph:
			while(true){
				// GTF �븳 以꾩쓣 �씫�쓬
				String line=BRGtf.readLine();
				// line�뿉 null�씪 寃쎌슦�뒗 �뙆�씪�쓽 �걹�엫
				// �빐�떦 寃쎌슦, �몢 媛�吏� 媛��뒫�꽦�씠 議댁옱
				// case1: �빐�떦 �뿼�깋泥� �젙蹂닿� �븳 以꾩씠�씪�룄 泥섎━�맂 寃쎌슦 (gene蹂��닔�뒗 null�씠 �븘�떂)
				// �씠 寃쎌슦 留덉�留� �쑀�쟾�옄 �젙蹂대�� 泥섎━�빐�빞�븿
				//
				// case2: �빐�떦 �뿼�깋泥� �젙蹂닿� �븳 以꾨룄 泥섎━媛� �븞�맂 寃쎌슦 (gene蹂��닔�뒗 null)
				// �씠 寃쎌슦 �쑀�쟾�옄 �젙蹂닿� �쟾�� �뾾�뒗 �긽�깭�씪�꽌 泥섎━�븷 寃껋씠 �뾾�쓬
				
				// 蹂� if臾몄뿉 �뱾�뼱媛� 寃쎌슦, 諛섎뱶�떆 醫낅즺�맖 
				if(line == null){
	
					// case2 �뿉 �빐�떦
					if(gene == null){
						break MakeExonGraph;
					}
					
					// case1 �뿉 �빐�떦
					int gene_trans_cnt = gene.get_trans_cnt();
					for(int i=0; i<gene_trans_cnt; i++){
						gene.transcriptList.get(i).reverseExon();
						EXON e =trans.get(i);
						while(e.get_start() != Integer.MAX_VALUE-1){
							if(i+1 != gene_trans_cnt){
								e.set_next_junc_dynamic(gene_trans_cnt-1,0, null);
								e.set_prev_junc_dynamic(gene_trans_cnt-1, 0, null);
							}
							e = e.get_next(i, 0);
						}
							
						e.get_prev(i, 0).set_tail(null);
							trans.get(i).get_next(i,0).set_head(null);
						
						if(i+1 != gene_trans_cnt){
							e.set_next_junc_dynamic(gene_trans_cnt-1,0, null);
							e.set_prev_junc_dynamic(gene_trans_cnt-1, 0, null);
						}
					}
					
					if (cls == -1 || cls == 1) genes.get(genes.size()-1).add_Alternative();
					if (cls == -1 || cls == 0) genes.get(genes.size()-1).add_Junction();
					if (cls != 2) genes.get(genes.size()-1).ExonGraphDiv();
					if (cls != 2) genes.get(genes.size()-1).ExonGraphMerge();
					if(cls != 2)
					genes.get(genes.size()-1).RemoveEdges();
					genes.get(genes.size()-1).set_intron();
					addIntron(genes.get(genes.size()-1));
					
					break MakeExonGraph;
				}
				
				// �씠 寃쎌슦濡� �꽆�뼱�삤硫�, line�� null�씠 �븘�땶寃� 蹂댁옣�맖
	
				//Meta info region
				// Meta Info�쓽 寃쎌슦�뒗 �쟾泥섎━�븷 �젙蹂닿� �뾾�쓬 (ver 1.0.0 �꽕怨� 湲곗�)
				if(line.startsWith("#")){
					continue;
				}
	
				String[] record = line.split("\t");
					
				// GTF�쓽 �빆紐⑹� 珥� 9媛쒕줈 援ъ꽦
				// �룷留룹뿉 留욎� �븡�쑝硫� �빐�떦 �씪�씤�� 臾댁떆
				
				if(record.length != 9){
					Logger.getLogger("ExonGraphGF").log(Level.WARNING, "There is a missing value caused by wrong GTF format.");
					continue;
				}
				
				
				// �쁽�옱 �떎猷⑤뒗 �뿼�깋泥� �젙蹂대쭔 嫄몃윭�깂
				// ex>
				// �뿼�깋泥� 1踰덈쭔 �떎猷щ떎硫�, 1踰덉뿉 ���븳 GTF �씪�씤留� �뿰�궛
				
				if( !record[0].equalsIgnoreCase(CHR_IDENTIFIER)){ continue; }
	
				
				// ���엯�씠 exon�씪 寃쎌슦
				String feature = record[2];
				if(feature.equalsIgnoreCase("exon") || feature.equalsIgnoreCase("CDS")){
					String[] attr = record[8].split(";");
					
					// �꽭 媛�吏� 寃쎌슦媛� �엳�쓬
					// case1: �쁽�옱 exon�씠 媛��옣 泥섏쓬 exon�씤 寃쎌슦
					// gene怨� trans瑜� 珥덇린�솕 �븯怨�, exon�쓣 異붽��븿
					//
					// case2: �쁽�옱 exon�쓽 TI媛� �씠�쟾 exon�쓽 TI�� �떎瑜� 寃쎌슦
					// 	case2-1: �쁽�옱 gene�씠 �씠�쟾 gene怨� �떎瑜� 寃쎌슦
					//	�씠�쟾 gene�쓣 genes�뿉 異붽��븯怨�, case1怨� 媛숈� 猷⑦떞�쓣 �닔�뻾�븿
					//
					//	case2-2: �쁽�옱 gene�씠 �씠�쟾 gene怨� 媛숈쓣 寃쎌슦
					//	�씠�쟾 gene�뿉 �깉濡쒖슫 �듃�옖�뒪�겕由쏀듃瑜� 異붽��븯怨�, exon�쓣 異붽��븿
					//
					// case3: �쁽�옱 exon�쓽 TI媛� �씠�쟾 exon�쓽 TI�� �룞�씪�븳 寃쎌슦
					// �씠�쟾 gene怨� �씠�쟾 trans�뿉 exon�쓣 異붽��븿
					
					// case3�쓽 寃쎌슦
					if(transcript_id != null && (transcript_id.equalsIgnoreCase(ExonUtils.getGtfAttr(attr,"transcript_id")))){
					
						gene.transcriptList.get(gene.get_trans_cnt() - 1).addExon(record, GeneStrand);
					}
					// case1, 2�쓽 寃쎌슦
					// case1怨� 2瑜� 臾띠� �씠�쑀�뒗 case1怨� case2媛� gene怨� trans瑜� 珥덇린�솕�븯�뒗 猷⑦떞�씠 以묐났�릺湲� �븣臾몄엫
					else{
						
						//Different Transcripts
						
						transcript_id = ExonUtils.getGtfAttr(attr, "transcript_id");
						
						// TI�뒗 exon援щ텇�뿉 �엳�뼱�꽌 諛섎뱶�떆 �븘�슂�븳 �슂�냼�엫
						// TI�뿉 ���빐�꽌 missing value媛� 諛쒖깮�븯硫�, 臾댁떆
						if(transcript_id == null){
							continue;
						}
						
						//Create Gene
						exon = null;
						prevExon_ = null;
						
						//Gene ID
						String GeneID = ExonUtils.getGtfAttr(attr, "gene_id");
						//case2-1�뿉 �빐�떦
						if(gene != null
								&& !prevGeneID.equalsIgnoreCase(GeneID)){

							int gene_trans_cnt = gene.get_trans_cnt();
							for(int i=0; i<gene_trans_cnt; i++){
								gene.transcriptList.get(i).reverseExon();
								EXON e =trans.get(i);
																
								while(e.get_start() != Integer.MAX_VALUE-1){
										
									if(i+1 != gene_trans_cnt){
										e.set_next_junc_dynamic(gene_trans_cnt-1,0, null);
										e.set_prev_junc_dynamic(gene_trans_cnt-1, 0, null);
									}
									e = e.get_next(i, 0);
								}
								
									e.get_prev(i, 0).set_tail(null);
									trans.get(i).get_next(i,0).set_head(null);
								
								
								if(i+1 != gene_trans_cnt){
									e.set_next_junc_dynamic(gene_trans_cnt-1,0, null);
									e.set_prev_junc_dynamic(gene_trans_cnt-1, 0, null);
								}
							}
							
							prevGeneID = null;
							
	
							if (cls == -1 || cls == 1) genes.get(genes.size()-1).add_Alternative();
							if (cls == -1 || cls == 0) genes.get(genes.size()-1).add_Junction();
							if (cls != 2) genes.get(genes.size()-1).ExonGraphDiv();
							if (cls != 2) genes.get(genes.size()-1).ExonGraphMerge();
							if(cls != 2)
							genes.get(genes.size()-1).RemoveEdges();
							genes.get(genes.size()-1).set_intron();
							addIntron(genes.get(genes.size()-1));
						}
						//Strand
						if(record[6].equalsIgnoreCase("-")){
							GeneStrand = false;
						}else{
							GeneStrand = true;
						}
						
						// case1怨� case2-1�뿉 �빐�떦
						// gene怨� trans瑜� 珥덇린�솕
						if(prevGeneID == null){
							gene = new GENE(GeneID, trans_cnt, GeneStrand);
							
							trans = new ArrayList<EXON>();
							
							gene.set_trans(trans);
							genes.add(gene);
						}
						
						// case2-2�뿉 �빐�떦
						if(prevGeneID != null && prevGeneID.equalsIgnoreCase(GeneID)){
							gene.incre_trans_cnt();
						}
						prevGeneID = GeneID;
						
						// exon�쓽 �떆�옉怨� �걹�쓣 留뚮벀
						exon =  new EXON(null, null, null, -1, -1, 0, gene.get_trans_cnt());
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, new EXON(null, null, null, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1, 0, gene.get_trans_cnt()));
						exon.get_next(gene.get_trans_cnt()-1, 0).set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_ = exon;
												
						gene.transcriptList.add(new Transcript());
						gene.transcriptList.get(gene.get_trans_cnt() - 1).addExon(record, GeneStrand);
						trans.add(exon);
					}
					
					// CDS�씤 寃쎌슦�뿉�뒗 CDS�젙蹂대�� exon�뿉 ���옣�븯吏� �븡�쓬.
					if(feature.equalsIgnoreCase("CDS")){
						continue;
					}
					
					// exon�쓽 �떆�옉怨� �걹�쓽 genome position
					int start =  Integer.parseInt(record[3]);
					int end =  Integer.parseInt(record[4]);
					
					// start蹂대떎 end媛� �겕嫄곕굹 媛숈븘�빞 �븿
					// �씠 �궗�빆�씠 蹂댁옣�릺�뼱�빞 addmutation怨� 媛숈� �뿰�궛�쓣 �닔�뻾�븷 �븣, �삱諛붾Ⅴ寃� �옉�룞�븿
					// dirty data�쓽 �븳 醫낅쪟�엫
					if(start > end){
						int temp = start;
						start = end;
						end = temp;
					}
					
					// Fasta�뿉�꽌 exon�쓽 臾몄옄�뿴�쓣 �씫�쓬
					String[] head_seq_tail = getFastaString(start, end, gene.get_strand());
					
					// exon �깮�꽦
					exon = new EXON(head_seq_tail[1], head_seq_tail[0], head_seq_tail[2], start, end, 0, gene.get_trans_cnt());
					
					exon.isExon = true;
	
					 
					boolean direction = false;
					
					if(prevExon_.get_start() < exon.get_start()){
						direction = true;
					}
					
					// start瑜� 湲곗��쑝濡� �삤由꾩감�닚 �젙�젹�씠 �릺�뼱�빞 �븿
					// direction�씠 true�씤 寃쎌슦�뒗 prevExon_ �떎�쓬�뿉 �쁽�옱 exon�쓣 �뿰寃고븯硫� �맖
					// direction�씠 false�씤 寃쎌슦�뒗 prevExon_ �씠�쟾�뿉 �쁽�옱 exon�쓣 �뿰寃고븯硫� �맖
					
					// direction�씠 true�씠硫� �빆�긽 �씠 猷⑦떞�쑝濡� �뿰寃�
					if(direction || prevExon_.get_seq() == null){
						exon.set_prev_junc(gene.get_trans_cnt()-1, 0, prevExon_);
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, prevExon_.get_next(gene.get_trans_cnt()-1, 0));
						prevExon_.get_next(gene.get_trans_cnt()-1, 0).set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_.set_next_junc(gene.get_trans_cnt()-1, 0, exon);
	
					}
					// direction�씠 false�씠硫� 珥덇린 �뿰寃곗쓣 �젣�쇅�븯怨� �빆�긽 �씠 猷⑦떞�쑝濡� �뿰寃�
					// 珥덇린 �뿰寃곗� �몢 寃쎌슦 �룞�씪�븯湲� �븣臾몄엫
					else if(prevExon_.get_seq() != null){
	
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, prevExon_);
						exon.set_prev_junc(gene.get_trans_cnt()-1, 0, prevExon_.get_prev(gene.get_trans_cnt()-1, 0));
						prevExon_.get_prev(gene.get_trans_cnt()-1, 0).set_next_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_.set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
					}
	
					
					prevExon_ = exon;
	
				}
			}
		
		// �뙆�씪 �떕湲�
		
		BRGtf.close();
		FRGtf.close();
		RAFFasta.close();
		
		//Chromosome �깮�꽦
		if(!CHR_IDENTIFIER.startsWith("chr")){
			CHR_IDENTIFIER = "chr"+CHR_IDENTIFIER;
		}
		chromosome[CUR_CHR_COUNT-1] = new Chromosome(genes, CHR_IDENTIFIER);
	}
	
	private String[] getFastaString(int start, int end, boolean strand) throws IOException{
		
		int interval = end-start+1;

		int diff = 0;
		if((start-HEAD_LENGTH)%whiteSpaceInterval==0) { diff = 1; }
		
		RAFFasta.seek(start -1 - HEAD_LENGTH   + fastaDescLength   +  (start-HEAD_LENGTH)/whiteSpaceInterval - diff);

		byte[] tempByte = new byte[interval +(interval+HEAD_LENGTH+TAIL_LENGTH)/whiteSpaceInterval + 1     + HEAD_LENGTH + TAIL_LENGTH];
		RAFFasta.read(tempByte);
		
		String fastaSeq = new String(tempByte);
		fastaSeq = fastaSeq.replaceAll("\n", "");
		
		
		String[] returnString = new String[3];
		fastaSeq = fastaSeq.toUpperCase();
				
		if(!strand){
			fastaSeq = ExonUtils.complement(fastaSeq);
		}
		
		returnString[0] = fastaSeq.substring(0, HEAD_LENGTH); //HEAD
		returnString[1] = fastaSeq.substring(HEAD_LENGTH, interval + HEAD_LENGTH); //SEQ
		returnString[2] = fastaSeq.substring(interval+HEAD_LENGTH, interval+HEAD_LENGTH+TAIL_LENGTH); //TAIL
		
		
		return returnString;
		
	}
	
	private void addIntron(GENE gene) throws IOException{
		ArrayList<ExonRangeType> ERTList = gene.intronscript.exonList;
		ArrayList<EXON> trans = gene.get_trans();
		Hashtable<Integer, EXON> frontIntronHash = new Hashtable<Integer, EXON>();
		Hashtable<Integer, EXON> backIntronHash = new Hashtable<Integer, EXON>();
		int trans_cnt = gene.get_trans_cnt();
		
		if(ERTList.size()!=0){
			// Intron Hash瑜� 留뚮벀
			for(ExonRangeType ERT : ERTList){
				// Fasta�뿉�꽌 exon�쓽 臾몄옄�뿴�쓣 �씫�쓬
				String[] head_seq_tail = getFastaString(ERT.start, ERT.end, gene.get_strand());
				
				// exon �깮�꽦
				EXON exon = new EXON(head_seq_tail[1], null, null, ERT.start, ERT.end, 0, gene.get_trans_cnt());
				exon.isExon = false;
				
				frontIntronHash.put(ERT.start-1, exon);
				backIntronHash.put(ERT.end+1, exon);
			}
			
			for(EXON exon : trans){
				exonTraversal(exon, trans_cnt, frontIntronHash, backIntronHash);
			}
		}
		
		
	}
	
	private void exonTraversal(EXON exon, int trans_cnt, Hashtable<Integer, EXON> frontIntronHash, Hashtable<Integer, EXON> backIntronHash){
		
		if(exon.state){
			return;
		}else{
			exon.state = true;
		}
		
		int start = exon.get_start();
		int end = exon.get_end();
		
		
		
		
		if(frontIntronHash.get(end) != null){
			EXON intron = frontIntronHash.get(end);
			exon.set_next_junc(0, Constants.INTRON_EDGE, intron);
			intron.set_prev_junc(0, Constants.INTRON_EDGE, exon);
			
			
			intron.state = true;
		}
		if(backIntronHash.get(start) != null){
			EXON intron = backIntronHash.get(start);
			intron.set_next_junc(0, Constants.INTRON_EDGE, exon);
			exon.set_prev_junc(0, Constants.INTRON_EDGE, intron);
			
			intron.state = true;
		}
		
		
		for(int i=0; i<trans_cnt; i++){
			for(int j=0; j<ExonGraph.JVALUE; j++){
				if(exon.get_next(i, j) != null){
					exonTraversal(exon.get_next(i, j), trans_cnt, frontIntronHash, backIntronHash);
				}
			}
		}
	}
}
