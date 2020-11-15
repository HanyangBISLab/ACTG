package ExonGraph;

import java.util.ArrayList;

@SuppressWarnings("serial")
public class Pattern extends ArrayList<Pattern> {
	
	private String Amino;
	private String output;
	private String Origin;
	
	public Pattern(String amino) {
		this.Amino = amino;
	}
	
	public Pattern(String amino, String Origin, String output) {
		this.Amino = amino;
		this.output = output;
		this.Origin = Origin;
	}
	
	public Pattern(Pattern tmp) {
		this.Amino = tmp.get_amino();
		this.output = tmp.output;
		this.Origin = tmp.Origin;
	}
	
	public String get_amino() { return this.Amino; }
	public String get_origin() { return this.Origin; }
	
	public void set_Amino(String str) { this.Amino = str; }
	public void set_Origin(String str) { this.Origin = str; }
	
	public String getOutput() {
		return output;
	}
	public void setOutput(String output) {
		this.output = new String(output);
	}
	
}
