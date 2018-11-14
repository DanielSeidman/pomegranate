
public class SVData {
    int startPos;
    int endPos;
    int chrom;
    String name;
	public SVData(int chr, int start, int end, int index) {
		startPos = start;
		endPos=end;
		chrom = chr;
		name = ""+index;
	}

}
