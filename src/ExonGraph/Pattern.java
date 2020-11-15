package ExonGraph;

import java.util.ArrayList;

@SuppressWarnings("serial")
public class Pattern extends ArrayList<Pattern> {
	
	private String Nucleotide;
	private String Amino;
	private String output;
	private String Origin;
	
	public Pattern(String amino, String nucl) {
		this.Amino = amino;
		this.Nucleotide = nucl;
	}
	
	public Pattern(String amino, String Origin, String nucl, String output) {
		this.Amino = amino;
		this.Nucleotide = nucl;
		this.output = output;
		this.Origin = Origin;
	}
	
	public Pattern(Pattern tmp) {
		this.Amino = tmp.get_amino();
		this.Nucleotide = tmp.get_nucleotide();
		this.output = tmp.output;
		this.Origin = tmp.Origin;
	}
	
	public String get_amino() { return this.Amino; }
	public String get_nucleotide() { return this.Nucleotide; }
	public String get_origin() { return this.Origin; }
	
	public void set_nucleotide(String str) { this.Nucleotide = str; }
	public void set_Amino(String str) { this.Amino = str; }
	public void set_Origin(String str) { this.Origin = str; }
	
	public String getOutput() {
		return output;
	}
	public void setOutput(String output) {
		this.output = new String(output);
	}
	
}
