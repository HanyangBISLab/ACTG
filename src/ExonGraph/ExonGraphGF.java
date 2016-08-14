package ExonGraph;
import java.io.BufferedInputStream;
/*
 * ExonGraphGF | ver 1.0.0
 * Writer: progistar
 * 
 * Description:
 * ExonGraph를 GTF파일과 FASTA파일을 이용하여 생성
 * Fasta파일의 이름은 chrNum.fa 형식
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
	
	// 염색체 개수를 저장
	// ExonGraphGF 클래스 간 공유
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
	
	// 유전자를 ArrayList로 저장
	private ArrayList<GENE> genes = null;
	// 트랜스크립트를 ArrayList로 저장
	private ArrayList<EXON> trans = null;
	
	// Junction Variation을 고려범위
	// Exon의 앞, 뒤 7 개의 문자를 고려
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
		
		// Mutation 처리
		
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
		// 염색체 개수 증가
		CHR_COUNT++;
		CUR_CHR_COUNT = CHR_COUNT;
		// ExonGraph 생성시,
		// Junction Variation, Alternative Splicing, 그리고
		// 메모리 최적화를 위한 옵션
		CLS = cls;
		
		File fasta = fastaFile;
		File gtf = gtfFile;
		
		// 해당 파일이 없을 경우, ExonGraph를 생성하지 못함
		if(!fasta.exists() || !gtf.exists()){
			System.out.println("File is not existed!");
			return;
		}
		
		genes = new ArrayList<GENE>();
		
		// 트랜스크립트 개수 초기화
		int trans_cnt = 1;
		
		// 
		String prevGeneID = null;
		
		// 염색체 식별자
		// Fasta File 이름이 식별자로 사용
		String CHR_IDENTIFIER = fastaFile.getName().substring(0, fastaFile.getName().lastIndexOf("."));
		
		GENE gene = null;
		EXON exon = null;
		EXON prevExon_ = null;
		
		String transcript_id = null;
		
		FileReader FRGtf = new FileReader(gtf);		
		BufferedReader BRGtf = new BufferedReader(FRGtf);
		RAFFasta = new RandomAccessFile(fastaFile, "r");
		
		// > 디스크립션 길이를 구함
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
		
		// whitespace의 간격을 구함
		// whitespace가 256byte를 넘어 갈 경우, bufCount를 계산하여 whitespaceInterval을 구함
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
				// GTF 한 줄을 읽음
				String line=BRGtf.readLine();
				// line에 null일 경우는 파일의 끝임
				// 해당 경우, 두 가지 가능성이 존재
				// case1: 해당 염색체 정보가 한 줄이라도 처리된 경우 (gene변수는 null이 아님)
				// 이 경우 마지막 유전자 정보를 처리해야함
				//
				// case2: 해당 염색체 정보가 한 줄도 처리가 안된 경우 (gene변수는 null)
				// 이 경우 유전자 정보가 전혀 없는 상태라서 처리할 것이 없음
				
				// 본 if문에 들어갈 경우, 반드시 종료됨 
				if(line == null){
	
					// case2 에 해당
					if(gene == null){
						break MakeExonGraph;
					}
					
					// case1 에 해당
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
				
				// 이 경우로 넘어오면, line은 null이 아닌게 보장됨
	
				//Meta info region
				// Meta Info의 경우는 전처리할 정보가 없음 (ver 1.0.0 설계 기준)
				if(line.startsWith("#")){
					continue;
				}
	
				String[] lineSplit = line.split("\t");
					
				// GTF의 항목은 총 9개로 구성
				// 포맷에 맞지 않으면 해당 라인은 무시
				
				if(lineSplit.length != 9){
					Logger.getLogger("ExonGraphGF").log(Level.WARNING, "There is a missing value caused by wrong GTF format.");
					continue;
				}
				
				
				// 현재 다루는 염색체 정보만 걸러냄
				// ex>
				// 염색체 1번만 다룬다면, 1번에 대한 GTF 라인만 연산
				
				if( !lineSplit[0].equalsIgnoreCase(CHR_IDENTIFIER)){ continue; }
	
				
				// 타입이 exon일 경우
				String feature = lineSplit[2];
				if(feature.equalsIgnoreCase("exon") || feature.equalsIgnoreCase("CDS")){
					String[] attr = lineSplit[8].split(";");
					
					// 세 가지 경우가 있음
					// case1: 현재 exon이 가장 처음 exon인 경우
					// gene과 trans를 초기화 하고, exon을 추가함
					//
					// case2: 현재 exon의 TI가 이전 exon의 TI와 다른 경우
					// 	case2-1: 현재 gene이 이전 gene과 다를 경우
					//	이전 gene을 genes에 추가하고, case1과 같은 루틴을 수행함
					//
					//	case2-2: 현재 gene이 이전 gene과 같을 경우
					//	이전 gene에 새로운 트랜스크립트를 추가하고, exon을 추가함
					//
					// case3: 현재 exon의 TI가 이전 exon의 TI와 동일한 경우
					// 이전 gene과 이전 trans에 exon을 추가함
					
					// case3의 경우
					if(transcript_id != null && (transcript_id.equalsIgnoreCase(ExonUtils.getGtfAttr(attr,"transcript_id")))){
					
						gene.transcriptList.get(gene.get_trans_cnt() - 1).addExon(lineSplit, GeneStrand);
					}
					// case1, 2의 경우
					// case1과 2를 묶은 이유는 case1과 case2가 gene과 trans를 초기화하는 루틴이 중복되기 때문임
					else{
						
						//Different Transcripts
						
						transcript_id = ExonUtils.getGtfAttr(attr, "transcript_id");
						
						// TI는 exon구분에 있어서 반드시 필요한 요소임
						// TI에 대해서 missing value가 발생하면, 무시
						if(transcript_id == null){
							continue;
						}
						
						//Create Gene
						exon = null;
						prevExon_ = null;
						
						//Gene ID
						String GeneID = ExonUtils.getGtfAttr(attr, "gene_id");
						//case2-1에 해당
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
						if(lineSplit[6].equalsIgnoreCase("-")){
							GeneStrand = false;
						}else{
							GeneStrand = true;
						}
						
						// case1과 case2-1에 해당
						// gene과 trans를 초기화
						if(prevGeneID == null){
							gene = new GENE(GeneID, trans_cnt, GeneStrand);
							
							trans = new ArrayList<EXON>();
							
							gene.set_trans(trans);
							genes.add(gene);
						}
						
						// case2-2에 해당
						if(prevGeneID != null && prevGeneID.equalsIgnoreCase(GeneID)){
							gene.incre_trans_cnt();
						}
						prevGeneID = GeneID;
						
						// exon의 시작과 끝을 만듦
						exon =  new EXON(null, null, null, -1, -1, 0, gene.get_trans_cnt());
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, new EXON(null, null, null, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1, 0, gene.get_trans_cnt()));
						exon.get_next(gene.get_trans_cnt()-1, 0).set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_ = exon;
												
						gene.transcriptList.add(new Transcript());
						gene.transcriptList.get(gene.get_trans_cnt() - 1).addExon(lineSplit, GeneStrand);
						trans.add(exon);
					}
					
					// CDS인 경우에는 CDS정보를 exon에 저장하지 않음.
					if(feature.equalsIgnoreCase("CDS")){
						continue;
					}
					
					// exon의 시작과 끝의 genome position
					int start =  Integer.parseInt(lineSplit[3]);
					int end =  Integer.parseInt(lineSplit[4]);
					
					// start보다 end가 크거나 같아야 함
					// 이 사항이 보장되어야 addmutation과 같은 연산을 수행할 때, 올바르게 작동함
					// dirty data의 한 종류임
					if(start > end){
						int temp = start;
						start = end;
						end = temp;
					}
					
					// Fasta에서 exon의 문자열을 읽음
					String[] head_seq_tail = getFastaString(start, end, gene.get_strand());
					
					// exon 생성
					exon = new EXON(head_seq_tail[1], head_seq_tail[0], head_seq_tail[2], start, end, 0, gene.get_trans_cnt());
					
					exon.isExon = true;
	
					 
					boolean direction = false;
					
					if(prevExon_.get_start() < exon.get_start()){
						direction = true;
					}
					
					// start를 기준으로 오름차순 정렬이 되어야 함
					// direction이 true인 경우는 prevExon_ 다음에 현재 exon을 연결하면 됨
					// direction이 false인 경우는 prevExon_ 이전에 현재 exon을 연결하면 됨
					
					// direction이 true이면 항상 이 루틴으로 연결
					if(direction || prevExon_.get_seq() == null){
						exon.set_prev_junc(gene.get_trans_cnt()-1, 0, prevExon_);
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, prevExon_.get_next(gene.get_trans_cnt()-1, 0));
						prevExon_.get_next(gene.get_trans_cnt()-1, 0).set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_.set_next_junc(gene.get_trans_cnt()-1, 0, exon);
	
					}
					// direction이 false이면 초기 연결을 제외하고 항상 이 루틴으로 연결
					// 초기 연결은 두 경우 동일하기 때문임
					else if(prevExon_.get_seq() != null){
	
						exon.set_next_junc(gene.get_trans_cnt()-1, 0, prevExon_);
						exon.set_prev_junc(gene.get_trans_cnt()-1, 0, prevExon_.get_prev(gene.get_trans_cnt()-1, 0));
						prevExon_.get_prev(gene.get_trans_cnt()-1, 0).set_next_junc(gene.get_trans_cnt()-1, 0, exon);
						prevExon_.set_prev_junc(gene.get_trans_cnt()-1, 0, exon);
					}
	
					
					prevExon_ = exon;
	
				}
			}
		
		// 파일 닫기
		
		BRGtf.close();
		FRGtf.close();
		RAFFasta.close();
		
		//Chromosome 생성
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
			// Intron Hash를 만듦
			for(ExonRangeType ERT : ERTList){
				// Fasta에서 exon의 문자열을 읽음
				String[] head_seq_tail = getFastaString(ERT.start, ERT.end, gene.get_strand());
				
				// exon 생성
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
