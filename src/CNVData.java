
public class CNVData {
    int startPos;
    int endPos;
    int chrom;
    String name;
	public CNVData(int chr, int start, int end, int index) {
		startPos = start;
		endPos=end;
		chrom = chr;
		name = ""+index;
	}

}

