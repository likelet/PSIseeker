/*
 * This class encodes the basic output format of PSI seeker with following columns 
 chr position base supporCountInTreat totalCountInTreat supporCountControl   totalCountControl Pvalue adjustpvalue
 */

package pisseeker;

import pisseeker.pub.FisherExactTest;



/**
 *
 * @author zhaoqi
 */
public class PSIout {
    
    private String chr;
    
    private int position; // 1 based  position
    private String strand;// strand information
    private boolean strandB;//boolean format for strand
    private char exbase=' ';//base at previous position
    private char base= ' '; // nuclic acid at candited position
    private String readsString ="";
    private int supporCountInTreat=0;
    
    private int totalCountInTreat=0;
    
    private int supporCountControl=0;
    
    private int totalCountControl = 0;
    
    private double Pvalue=1;//fisher test Pvalue
    
    private double adjustP = 1;//ajust Pvalue
    
    
    
    

    public PSIout(String chr, int position) {
        this.chr = chr;
        this.position = position;
        
    }

    public PSIout(String chr, int position, boolean strand) {
        this.chr = chr;
        this.position = position;
        if(strand){
            this.strand = "-";
        }else{
            this.strand = "+";
        }
        this.strandB=strand;
        
    }
    
    
   

    public boolean isStrandB() {
        return strandB;
    }

    public void setStrandB(boolean strandB) {
        this.strandB = strandB;
    }


    
    /**
     * Get the value of adjustP
     *
     * @return the value of adjustP
     */
    public double getAdjustP() {
        return adjustP;
    }

    /**
     * Set the value of adjustP
     *
     * @param adjustP new value of adjustP
     */
    public void setAdjustP(double adjustP) {
        this.adjustP = adjustP;
    }


    /**
     * Get the value of Pvalue
     *
     * @return the value of Pvalue
     */
    public double getPvalue() {
        return Pvalue;
    }

    /**
     * Set the value of Pvalue
     *
     * @param Pvalue new value of Pvalue
     */
    public void setPvalue(double Pvalue) {
        this.Pvalue = Pvalue;
    }


    /**
     * Get the value of totalCountControl
     *
     * @return the value of totalCountControl
     */
    public int getTotalCountControl() {
        return totalCountControl;
    }

    /**
     * Set the value of totalCountControl
     *
     * @param totalCountControl new value of totalCountControl
     */
    public void setTotalCountControl(int totalCountControl) {
        this.totalCountControl = totalCountControl;
    }

    

    /**
     * Get the value of supporCountControl
     *
     * @return the value of supporCountControl
     */
    public int getSupporCountControl() {
        return supporCountControl;
    }

    /**
     * Set the value of supporCountControl
     *
     * @param supporCountControl new value of supporCountControl
     */
    public void setSupporCountControl(int supporCountControl) {
        this.supporCountControl = supporCountControl;
    }


    /**
     * Get the value of totalCountInTreat
     *
     * @return the value of totalCountInTreat
     */
    public int getTotalCountInTreat() {
        return totalCountInTreat;
    }

    /**
     * Set the value of totalCountInTreat
     *
     * @param totalCountInTreat new value of totalCountInTreat
     */
    public void setTotalCountInTreat(int totalCountInTreat) {
        this.totalCountInTreat = totalCountInTreat;
    }


    /**
     * Get the value of supporCountInTreat
     *
     * @return the value of supporCountInTreat
     */
    public int getSupporCountInTreat() {
        return supporCountInTreat;
    }

    /**
     * Set the value of supporCountInTreat
     *
     * @param supporCountInTreat new value of supporCountInTreat
     */
    public void setSupporCountInTreat(int supporCountInTreat) {
        this.supporCountInTreat = supporCountInTreat;
    }


    /**
     * Get the value of base
     *
     * @return the value of base
     */
    public char getBase() {
        return base;
    }

    /**
     * Set the value of base
     *
     * @param base new value of base
     */
    public void setBase(char base) {
        this.base = base;
    }

    public String getReadsString() {
        return readsString;
    }

    public void setReadsString(String readsString) {
        this.readsString = readsString;
    }


    /**
     * Get the value of position
     *
     * @return the value of position
     */
    public int getPosition() {
        return position;
    }

    /**
     * Set the value of position
     *
     * @param position new value of position
     */
    public void setPosition(int position) {
        this.position = position;
    }

    /**
     * Get the value of chr
     *
     * @return the value of chr
     */
    public String getChr() {
        return chr;
    }

    /**
     * Set the value of chr
     *
     * @param chr new value of chr
     */
    public void setChr(String chr) {
        this.chr = chr;
    }

    public char getExbase() {
        return exbase;
    }

    public void setExbase(char exbase) {
        this.exbase = exbase;
    }

    
    
    
    
    public void add1supporCountInTreat(){
        this.supporCountInTreat++;
    }
    
    public void add1totalCountInTreat(){
        this.totalCountInTreat++;
    }
    public void add1supporCountControl(){
        this.supporCountControl++;
    }
    public void add1totalCountControl(){
        this.totalCountControl++;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }
    
    
    @Override

    public int hashCode() {

//        System.out.println("------hashCode----");

        return 1;

    }

    @Override

    public boolean equals(Object obj) {

//        System.out.println("------equals----");

        if (obj instanceof PSIout) {

            PSIout po = (PSIout) obj;

            if (this.position==(po.getPosition())) {

                return true;

            }

        }

        return false;

    }
    
    //TEST attributive
    private double enrichmentScore=1;

    public double getEnrichmentScore() {
        return enrichmentScore;
    }

    public void setEnrichmentScore(double enrichmentScore) {
        this.enrichmentScore = enrichmentScore;
    }
    
        private double poissonP;

    /**
     * Get the value of poissonP
     *
     * @return the value of poissonP
     */
    public double getPoissonP() {
        return poissonP;
    }

    /**
     * Set the value of poissonP
     *
     * @param poissonP new value of poissonP
     */
    public void setPoissonP(double poissonP) {
        this.poissonP = poissonP;
    }

    private double swapFDR=1;

    /**
     * Get the value of swapFDR
     *
     * @return the value of swapFDR
     */
    public double getSwapFDR() {
        return swapFDR;
    }

    /**
     * Set the value of swapFDR
     *
     * @param swapFDR new value of swapFDR
     */
    public void setSwapFDR(double swapFDR) {
        this.swapFDR = swapFDR;
    }

    
    @Override
    public String toString() {
        return chr + "\t" + position + "\t" +exbase+"\t"+ base + "\t"+readsString+"\t"+strand+"\t" + supporCountInTreat + "\t" + totalCountInTreat +"\t" + supporCountControl +"\t" + totalCountControl + "\t" + Pvalue + "\t" + adjustP ;
    }
    
    public String toString2() {
        return chr + "\t" + position + "\t" +exbase+"\t"+ base + "\t"+readsString+"\t"+strand+"\t" + supporCountInTreat + "\t" + totalCountInTreat +"\t" + supporCountControl +"\t" + totalCountControl + "\t" + Pvalue + "\t" + adjustP+"\t"+this.enrichmentScore ;
    }
    

    
    
}
