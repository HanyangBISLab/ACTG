package ExonGraph;

import java.util.ArrayList;

@SuppressWarnings("serial")
public class Pattern extends ArrayList<Pattern> {
	
	private String Nucleotide;
	private String Amino;
	private double p_flank, n_flank;
	private String output;
	private String Origin;
	
	public Pattern(String amino, String nucl, double p, double n) {
		this.Amino = amino;
		this.Nucleotide = nucl;
		this.p_flank = p;
		this.n_flank = n;
	}
	
	public Pattern(String amino, String Origin, String nucl, double p, double n, String output) {
		this.Amino = amino;
		this.Nucleotide = nucl;
		this.p_flank = p;
		this.n_flank = n;
		this.output = output;
		this.Origin = Origin;
	}
	
	public Pattern(Pattern tmp) {
		this.Amino = tmp.get_amino();
		this.Nucleotide = tmp.get_nucleotide();
		this.p_flank = tmp.get_pflank();
		this.n_flank = tmp.get_nflank();
		this.output = tmp.output;
		this.Origin = tmp.Origin;
	}
	
	public double get_pflank() { return this.p_flank; }
	public double get_nflank() { return this.n_flank; }
	public String get_amino() { return this.Amino; }
	public String get_nucleotide() { return this.Nucleotide; }
	public String get_origin() { return this.Origin; }
	
	public void set_nucleotide(String str) { this.Nucleotide = str; }
	public void set_Amino(String str) { this.Amino = str; }
	public void set_pflank(double tmp) { this.p_flank = tmp; }
	public void set_nflank(double tmp) { this.n_flank = tmp; }
	public void set_Origin(String str) { this.Origin = str; }
	
	public String getOutput() {
		return output;
	}
	public void setOutput(String output) {
		this.output = new String(output);
	}
	
}
