package Environments;

import java.util.HashMap;

public class Codon {
	static HashMap<String, Character> NuclToAmino;
	static HashMap<String, Character> NuclToAminoRev_R;
	static HashMap<String, Character> NuclToAminoRev_R_C;

	static public void Mapping() {
		NuclToAmino = new HashMap<String, Character>();

		/// --------------------------------- N to A -----------------------------------------//

		Character tmp = 'F';
		NuclToAmino.put("TTT", tmp);NuclToAmino.put("TTC", tmp);

		tmp = 'L';
		NuclToAmino.put("TTA", tmp);NuclToAmino.put("TTG", tmp);NuclToAmino.put("CTT", tmp);
		NuclToAmino.put("CTC", tmp);NuclToAmino.put("CTA", tmp);NuclToAmino.put("CTG", tmp);

		if(Constants.IS_I_SAME_WITH_L) tmp = 'L';
		else tmp = 'I';
		NuclToAmino.put("ATT", tmp);NuclToAmino.put("ATC", tmp);NuclToAmino.put("ATA", tmp);

		tmp = 'M';
		NuclToAmino.put("ATG", tmp);

		tmp = 'V';
		NuclToAmino.put("GTT", tmp);NuclToAmino.put("GTC", tmp);NuclToAmino.put("GTA", tmp);
		NuclToAmino.put("GTG", tmp);

		tmp = 'S';
		NuclToAmino.put("TCT", tmp);NuclToAmino.put("TCC", tmp);NuclToAmino.put("TCA", tmp);
		NuclToAmino.put("TCG", tmp);NuclToAmino.put("AGT", tmp);NuclToAmino.put("AGC", tmp);

		tmp = 'P';
		NuclToAmino.put("CCT", tmp);NuclToAmino.put("CCC", tmp);NuclToAmino.put("CCA", tmp);
		NuclToAmino.put("CCG", tmp);

		tmp = 'T';
		NuclToAmino.put("ACT", tmp);NuclToAmino.put("ACC", tmp);NuclToAmino.put("ACA", tmp);
		NuclToAmino.put("ACG", tmp);

		tmp = 'A';
		NuclToAmino.put("GCT", tmp);NuclToAmino.put("GCC", tmp);NuclToAmino.put("GCA", tmp);
		NuclToAmino.put("GCG", tmp);

		tmp = 'Y';
		NuclToAmino.put("TAT", tmp);NuclToAmino.put("TAC", tmp);

		tmp = 'X';
		NuclToAmino.put("TAA", tmp);NuclToAmino.put("TGA", tmp);NuclToAmino.put("TAG", tmp);

		tmp = 'H';
		NuclToAmino.put("CAT", tmp);NuclToAmino.put("CAC", tmp);

		tmp = 'Q';
		NuclToAmino.put("CAA", tmp);NuclToAmino.put("CAG", tmp);

		tmp = 'N';
		NuclToAmino.put("AAT", tmp);NuclToAmino.put("AAC", tmp);

		tmp = 'K';
		NuclToAmino.put("AAA", tmp);NuclToAmino.put("AAG", tmp);

		tmp = 'D';
		NuclToAmino.put("GAT", tmp);NuclToAmino.put("GAC", tmp);

		tmp = 'E';
		NuclToAmino.put("GAA", tmp);NuclToAmino.put("GAG", tmp);

		tmp = 'C';
		NuclToAmino.put("TGT", tmp);NuclToAmino.put("TGC", tmp);

		tmp = 'W';
		NuclToAmino.put("TGG", tmp);

		tmp = 'R';
		NuclToAmino.put("CGT", tmp);NuclToAmino.put("CGC", tmp);NuclToAmino.put("CGA", tmp);
		NuclToAmino.put("CGG", tmp);NuclToAmino.put("AGA", tmp);NuclToAmino.put("AGG", tmp);

		tmp = 'G';
		NuclToAmino.put("GGT", tmp);NuclToAmino.put("GGC", tmp);NuclToAmino.put("GGA", tmp);
		NuclToAmino.put("GGG", tmp);


		// The end of down stream exon in transcript
		tmp = 'S';
		NuclToAmino.put("TCN", tmp);
		
		tmp = 'L';
		NuclToAmino.put("CTN", tmp);
		
		tmp = 'P';
		NuclToAmino.put("CCN", tmp);
		
		tmp = 'R';
		NuclToAmino.put("CGN", tmp);
		
		tmp = 'T';
		NuclToAmino.put("ACN", tmp);
		
		tmp = 'V';
		NuclToAmino.put("GTN", tmp);
		
		tmp = 'A';
		NuclToAmino.put("GCN", tmp);
		
		tmp = 'G';
		NuclToAmino.put("GGN", tmp);

		/*
		 * Reverse Table
		 */


		NuclToAminoRev_R = new HashMap<String, Character>();

		/// --------------------------------- N to A for Reverse -----------------------------------------//
		tmp = 'F';
		NuclToAminoRev_R.put("TTT", tmp);NuclToAminoRev_R.put("CTT", tmp);

		tmp = 'L';
		NuclToAminoRev_R.put("ATT", tmp);NuclToAminoRev_R.put("GTT", tmp);NuclToAminoRev_R.put("TTC", tmp);
		NuclToAminoRev_R.put("CTC", tmp);NuclToAminoRev_R.put("ATC", tmp);NuclToAminoRev_R.put("GTC", tmp);
		
		if(Constants.IS_I_SAME_WITH_L) tmp = 'L';
		else tmp = 'I';
		NuclToAminoRev_R.put("TTA", tmp);NuclToAminoRev_R.put("CTA", tmp);NuclToAminoRev_R.put("ATA", tmp);

		tmp = 'M';
		NuclToAminoRev_R.put("GTA", tmp);

		tmp = 'V';
		NuclToAminoRev_R.put("TTG", tmp);NuclToAminoRev_R.put("CTG", tmp);NuclToAminoRev_R.put("ATG", tmp);
		NuclToAminoRev_R.put("GTG", tmp);

		tmp = 'S';
		NuclToAminoRev_R.put("TCT", tmp);NuclToAminoRev_R.put("CCT", tmp);NuclToAminoRev_R.put("ACT", tmp);
		NuclToAminoRev_R.put("GCT", tmp);NuclToAminoRev_R.put("TGA", tmp);NuclToAminoRev_R.put("CGA", tmp);

		tmp = 'P';
		NuclToAminoRev_R.put("TCC", tmp);NuclToAminoRev_R.put("CCC", tmp);NuclToAminoRev_R.put("ACC", tmp);
		NuclToAminoRev_R.put("GCC", tmp);

		tmp = 'T';
		NuclToAminoRev_R.put("TCA", tmp);NuclToAminoRev_R.put("CCA", tmp);NuclToAminoRev_R.put("ACA", tmp);
		NuclToAminoRev_R.put("GCA", tmp);

		tmp = 'A';
		NuclToAminoRev_R.put("TCG", tmp);NuclToAminoRev_R.put("CCG", tmp);NuclToAminoRev_R.put("ACG", tmp);
		NuclToAminoRev_R.put("GCG", tmp);

		tmp = 'Y';
		NuclToAminoRev_R.put("TAT", tmp);NuclToAminoRev_R.put("CAT", tmp);

		tmp = 'X';
		NuclToAminoRev_R.put("AAT", tmp);NuclToAminoRev_R.put("AGT", tmp);NuclToAminoRev_R.put("GAT", tmp);

		tmp = 'H';
		NuclToAminoRev_R.put("TAC", tmp);NuclToAminoRev_R.put("CAC", tmp);

		tmp = 'Q';
		NuclToAminoRev_R.put("AAC", tmp);NuclToAminoRev_R.put("GAC", tmp);

		tmp = 'N';
		NuclToAminoRev_R.put("TAA", tmp);NuclToAminoRev_R.put("CAA", tmp);

		tmp = 'K';
		NuclToAminoRev_R.put("AAA", tmp);NuclToAminoRev_R.put("GAA", tmp);

		tmp = 'D';
		NuclToAminoRev_R.put("TAG", tmp);NuclToAminoRev_R.put("CAG", tmp);

		tmp = 'E';
		NuclToAminoRev_R.put("AAG", tmp);NuclToAminoRev_R.put("GAG", tmp);

		tmp = 'C';
		NuclToAminoRev_R.put("TGT", tmp);NuclToAminoRev_R.put("CGT", tmp);

		tmp = 'W';
		NuclToAminoRev_R.put("GGT", tmp);

		tmp = 'R';
		NuclToAminoRev_R.put("TGC", tmp);NuclToAminoRev_R.put("CGC", tmp);NuclToAminoRev_R.put("AGC", tmp);
		NuclToAminoRev_R.put("GGC", tmp);NuclToAminoRev_R.put("AGA", tmp);NuclToAminoRev_R.put("GGA", tmp);

		tmp = 'G';
		NuclToAminoRev_R.put("TGG", tmp);NuclToAminoRev_R.put("CGG", tmp);NuclToAminoRev_R.put("AGG", tmp);
		NuclToAminoRev_R.put("GGG", tmp);
				
		// The end of down stream exon in transcript
		tmp = 'S';
		NuclToAmino.put("NCT", tmp);
		
		tmp = 'L';
		NuclToAmino.put("NTC", tmp);
		
		tmp = 'P';
		NuclToAmino.put("NCC", tmp);
		
		tmp = 'R';
		NuclToAmino.put("NGC", tmp);
		
		tmp = 'T';
		NuclToAmino.put("NCA", tmp);
		
		tmp = 'V';
		NuclToAmino.put("NTG", tmp);
		
		tmp = 'A';
		NuclToAmino.put("NCG", tmp);
		
		tmp = 'G';
		NuclToAmino.put("NGG", tmp);

		
		NuclToAminoRev_R_C = new HashMap<String, Character>();

		/// --------------------------------- N to A for Reverse (complementary) -----------------------------------------//
		//simply changed "NuclToAminoRev" to be complementary
		
		tmp = 'F';
		NuclToAminoRev_R_C.put("AAA", tmp);NuclToAminoRev_R_C.put("GAA", tmp);

		tmp = 'L';
		NuclToAminoRev_R_C.put("TAA", tmp);NuclToAminoRev_R_C.put("CAA", tmp);NuclToAminoRev_R_C.put("AAG", tmp);
		NuclToAminoRev_R_C.put("GAG", tmp);NuclToAminoRev_R_C.put("TAG", tmp);NuclToAminoRev_R_C.put("CAG", tmp);
		NuclToAminoRev_R_C.put("NAG", tmp);
		
		if(Constants.IS_I_SAME_WITH_L) tmp = 'L';
		else tmp = 'I';
		NuclToAminoRev_R_C.put("AAT", tmp);NuclToAminoRev_R_C.put("GAT", tmp);NuclToAminoRev_R_C.put("TAT", tmp);

		tmp = 'M';
		NuclToAminoRev_R_C.put("CAT", tmp);

		tmp = 'V';
		NuclToAminoRev_R_C.put("AAC", tmp);NuclToAminoRev_R_C.put("GAC", tmp);NuclToAminoRev_R_C.put("TAC", tmp);
		NuclToAminoRev_R_C.put("CAC", tmp);NuclToAminoRev_R_C.put("NAC", tmp);

		tmp = 'S';
		NuclToAminoRev_R_C.put("AGA", tmp);NuclToAminoRev_R_C.put("GGA", tmp);NuclToAminoRev_R_C.put("TGA", tmp);
		NuclToAminoRev_R_C.put("CGA", tmp);NuclToAminoRev_R_C.put("ACT", tmp);NuclToAminoRev_R_C.put("GCT", tmp);
		NuclToAminoRev_R_C.put("NGA", tmp);

		tmp = 'P';
		NuclToAminoRev_R_C.put("AGG", tmp);NuclToAminoRev_R_C.put("GGG", tmp);NuclToAminoRev_R_C.put("TGG", tmp);
		NuclToAminoRev_R_C.put("CGG", tmp);NuclToAminoRev_R_C.put("NGG", tmp);

		tmp = 'T';
		NuclToAminoRev_R_C.put("AGT", tmp);NuclToAminoRev_R_C.put("GGT", tmp);NuclToAminoRev_R_C.put("TGT", tmp);
		NuclToAminoRev_R_C.put("CGT", tmp);NuclToAminoRev_R_C.put("NGT", tmp);

		tmp = 'A';
		NuclToAminoRev_R_C.put("AGC", tmp);NuclToAminoRev_R_C.put("GGC", tmp);NuclToAminoRev_R_C.put("TGC", tmp);
		NuclToAminoRev_R_C.put("CGC", tmp);NuclToAminoRev_R_C.put("NGC", tmp);

		tmp = 'Y';
		NuclToAminoRev_R_C.put("ATA", tmp);NuclToAminoRev_R_C.put("GTA", tmp);

		tmp = 'X';
		NuclToAminoRev_R_C.put("TTA", tmp);NuclToAminoRev_R_C.put("TCA", tmp);NuclToAminoRev_R_C.put("CTA", tmp);

		tmp = 'H';
		NuclToAminoRev_R_C.put("ATG", tmp);NuclToAminoRev_R_C.put("GTG", tmp);

		tmp = 'Q';
		NuclToAminoRev_R_C.put("TTG", tmp);NuclToAminoRev_R_C.put("CTG", tmp);
		
		tmp = 'N';
		NuclToAminoRev_R_C.put("ATT", tmp);NuclToAminoRev_R_C.put("GTT", tmp);
		
		tmp = 'K';
		NuclToAminoRev_R_C.put("TTT", tmp);NuclToAminoRev_R_C.put("CTT", tmp);
		
		tmp = 'D';
		NuclToAminoRev_R_C.put("ATC", tmp);NuclToAminoRev_R_C.put("GTC", tmp);
		
		tmp = 'E';
		NuclToAminoRev_R_C.put("TTC", tmp);NuclToAminoRev_R_C.put("CTC", tmp);

		tmp = 'C';
		NuclToAminoRev_R_C.put("ACA", tmp);NuclToAminoRev_R_C.put("GCA", tmp);

		tmp = 'W';
		NuclToAminoRev_R_C.put("CCA", tmp);

		tmp = 'R';
		NuclToAminoRev_R_C.put("ACG", tmp);NuclToAminoRev_R_C.put("GCG", tmp);NuclToAminoRev_R_C.put("TCG", tmp);
		NuclToAminoRev_R_C.put("CCG", tmp);NuclToAminoRev_R_C.put("TCT", tmp);NuclToAminoRev_R_C.put("CCT", tmp);
		NuclToAminoRev_R_C.put("NCG", tmp);

		tmp = 'G';
		NuclToAminoRev_R_C.put("ACC", tmp);NuclToAminoRev_R_C.put("GCC", tmp);NuclToAminoRev_R_C.put("TCC", tmp);
		NuclToAminoRev_R_C.put("CCC", tmp);NuclToAminoRev_R_C.put("NCC", tmp);
	}

	public static Character NuclToAmino(String a) {
		if(NuclToAmino.get(a.toString()) == null)
			return 'X';
		else
			return NuclToAmino.get(a.toString());
	}

	public static Character NuclToAmino_R(String a) {
		if(NuclToAminoRev_R.get(a.toString()) == null)
			return 'X';
		else
			return NuclToAminoRev_R.get(a.toString());
	}
	
	public static Character NuclToAmino_R_C(String a) {
		if(NuclToAminoRev_R_C.get(a.toString()) == null)
			return 'X';
		else
			return NuclToAminoRev_R_C.get(a.toString());
	}
}
