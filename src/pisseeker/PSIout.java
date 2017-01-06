/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

/**
 *
 * @author Administrator
 * @since 2017-1-6
 * @coding time 11:24:56
 * @author Qi Zhao
 */
public class PSIout {

    private String chr;
    private int position;
    private String strand;
    private boolean strandB;
    private char exbase = ' ';
    private char base = ' ';
    private String readsString = "";
    private int supporCountInTreat = 0;
    private int totalCountInTreat = 0;
    private int supporCountControl = 0;
    private int totalCountControl = 0;
    private double FisherPvalue = 1.0D;
    private double fihserAdjustP = 1.0D;
    private double enrichmentScore = 1.0D;
    private double swapFDR = 1.0D;
    private String firstCigar;
    private double backRatio=0;
    private double lamda[] =new double[3]; //local lamda list with each conrespponding lamdaBG,lamda 1k, lamda 5k
    private double PoissonP=1;

    /**
     * Get the value of PoissonP
     *
     * @return the value of PoissonP
     */
    public double getPoissonP() {
        return PoissonP;
    }

    /**
     * Set the value of PoissonP
     *
     * @param PoissonP new value of PoissonP
     */
    public void setPoissonP(double PoissonP) {
        this.PoissonP = PoissonP;
    }


    public PSIout(String chr, int position) {
        this.chr = chr;
        this.position = position;
    }

    public PSIout(String chr, int position, boolean strand) {
        this.chr = chr;
        this.position = position;
        if (strand) {
            this.strand = "-";
        } else {
            this.strand = "+";
        }
        this.strandB = strand;
    }

    public boolean isStrandB() {
        return this.strandB;
    }

    public void setStrandB(boolean strandB) {
        this.strandB = strandB;
    }

    public double getFihserAdjustP() {
        return this.fihserAdjustP;
    }

    public void setFihserAdjustP(double fihserAdjustP) {
        this.fihserAdjustP = fihserAdjustP;
    }

    public double getFisherPvalue() {
        return this.FisherPvalue;
    }

    public void setFisherPvalue(double FisherPvalue) {
        this.FisherPvalue = FisherPvalue;
    }

    public int getTotalCountControl() {
        return this.totalCountControl;
    }

    public void setTotalCountControl(int totalCountControl) {
        this.totalCountControl = totalCountControl;
    }

    public int getSupporCountControl() {
        return this.supporCountControl;
    }

    public void setSupporCountControl(int supporCountControl) {
        this.supporCountControl = supporCountControl;
    }

    public int getTotalCountInTreat() {
        return this.totalCountInTreat;
    }

    public void setTotalCountInTreat(int totalCountInTreat) {
        this.totalCountInTreat = totalCountInTreat;
    }

    public int getSupporCountInTreat() {
        return this.supporCountInTreat;
    }

    public void setSupporCountInTreat(int supporCountInTreat) {
        this.supporCountInTreat = supporCountInTreat;
    }

    public char getBase() {
        return this.base;
    }

    public void setBase(char base) {
        this.base = base;
    }

    public String getReadsString() {
        return this.readsString;
    }

    public void setReadsString(String readsString) {
        this.readsString = readsString;
    }

    public int getPosition() {
        return this.position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public String getChr() {
        return this.chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public char getExbase() {
        return this.exbase;
    }

    public void setExbase(char exbase) {
        this.exbase = exbase;
    }

    public void add1supporCountInTreat() {
        this.supporCountInTreat += 1;
    }

    public void add1totalCountInTreat() {
        this.totalCountInTreat += 1;
    }

    public void add1supporCountControl() {
        this.supporCountControl += 1;
    }

    public void add1totalCountControl() {
        this.totalCountControl += 1;
    }

    public String getStrand() {
        return this.strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public int hashCode() {
        return 1;
    }

    public boolean equals(Object obj) {
        if ((obj instanceof PSIout)) {
            PSIout po = (PSIout) obj;

            if (this.position == po.getPosition()) {
                return true;
            }

        }

        return false;
    }

    public double getEnrichmentScore() {
        return this.enrichmentScore;
    }

    public void setEnrichmentScore(double enrichmentScore) {
        this.enrichmentScore = enrichmentScore;
    }

   

    public double getSwapFDR() {
        return this.swapFDR;
    }

    public void setSwapFDR(double swapFDR) {
        this.swapFDR = swapFDR;
    }

    public String getFirstCigar() {
        return this.firstCigar;
    }

    public void setFirstCigar(String firstCigar) {
        this.firstCigar = firstCigar;
    }

    public void calBackRatio(){
        this.backRatio=(double)supporCountControl/(double)totalCountControl;
    }

    public double getBackRatio() {
        return backRatio;
    }
    

    public double[] getLamda() {
        return lamda;
    }

    public void setLamda(double[] lamda) {
        this.lamda = lamda;
    }
    
    
    public String toString() {
        return this.chr + "\t" + this.position + "\t" + this.exbase + "\t" + this.base + "\t" + this.readsString + "\t" + this.strand + "\t" + this.supporCountInTreat + "\t" + this.totalCountInTreat + "\t" + this.supporCountControl + "\t" + this.totalCountControl + "\t" + this.FisherPvalue + "\t" + this.fihserAdjustP;
    }

    public String toString2() {
        return this.chr + "\t" + this.position + "\t" + this.exbase + "\t" + this.base + "\t" + this.firstCigar + "\t" + this.readsString + "\t" + this.strand + "\t" + this.supporCountInTreat + "\t" + this.totalCountInTreat + "\t" + this.supporCountControl + "\t" + this.totalCountControl + "\t" + this.FisherPvalue + "\t" + this.fihserAdjustP + "\t" + this.enrichmentScore;
    }
}
