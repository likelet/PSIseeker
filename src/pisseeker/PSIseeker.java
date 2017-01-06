package pisseeker;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import pisseeker.pub.DNAsequenceProcess;
import pisseeker.pub.FisherExactTest;

public final class PSIseeker {

    private final SamReader srt;
    private final SamReader src;
    public int Thread = 3;
    public String tbam;
    public String cbam;
    public IndexedFastaSequenceFile Indexgenomefile;

    private HashMap<String, HashSet<PSIout>> positiveResultMap = new HashMap();
    private HashMap<String, HashSet<PSIout>> negativeResultMap = new HashMap();
    
    private HashMap<String, IntervalTree>  positiveBackGroundMap = new HashMap<String, IntervalTree>();// result in positve strand 
    private HashMap<String, IntervalTree>  negativeBackGroundMap = new HashMap<String, IntervalTree>();// result in negative strand;
    

    private HashMap<String, IntervalTree> positiveTreeControlMap = new HashMap();
    private HashMap<String, IntervalTree> negativeTreeControlMap = new HashMap();
    private HashMap<String, IntervalTree> positiveTreeTreatMap = new HashMap();
    private HashMap<String, IntervalTree> negativeTreeTreatMap = new HashMap();
    
    private HashMap<String,Double> positiveGloableLamda=new HashMap<String,Double> ();
    private HashMap<String,Double> negativeGloableLamda=new HashMap<String,Double> ();
    private HashMap<String,Double> positiveGloableNormalizedFactor=new HashMap<String,Double> ();
    private HashMap<String,Double> negativeGloableNormalizedFactor=new HashMap<String,Double> ();
    public int filternumber = 2;
    public double filterTreatRatio = 0.1D;
    public double enrichmentThreshold = 0D;

    private HashSet<String> chrlist = new HashSet();

    public PSIseeker(String tbam, String cbam, String genomefile)
            throws IOException {
        this.cbam = cbam;
        this.tbam = tbam;
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        this.Indexgenomefile = new IndexedFastaSequenceFile(new File(genomefile));
        this.srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile1)
                .index(new File(bamfile1
                        .getAbsolutePath() + ".bai")));

        this.src = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile2)
                .index(new File(bamfile2
                        .getAbsolutePath() + ".bai")));
    }

    public void initializeTreatIntervalTree() {
        SAMRecordIterator itertemp = srt.iterator();
        //get chrome information
        SAMFileHeader fileHeader = srt.getFileHeader();
        SAMSequenceDictionary seqDic = fileHeader.getSequenceDictionary();
        List<SAMSequenceRecord> seqRecList = seqDic.getSequences();
        for (SAMSequenceRecord seqRec : seqRecList) {
            //store chr information
            this.chrlist.add(seqRec.getSequenceName());
        }
        System.out.println("Initializing initializeTreatIntervalTree  and site infomation");
        while (itertemp.hasNext()) {
            SAMRecord sitem = (SAMRecord) itertemp.next();
            if ((sitem.getReadUnmappedFlag())
                    || (sitem.getDuplicateReadFlag())
                    || (sitem.getNotPrimaryAlignmentFlag())
                    || (sitem.getReadFailsVendorQualityCheckFlag())
                    ) {
                continue;
            }
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();
            String chrome = sitem.getReferenceName();
            if (strand) {
                if (!this.negativeTreeTreatMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    this.negativeTreeTreatMap.put(chrome, potree);
                } else {
                    ((IntervalTree) this.negativeTreeTreatMap.get(chrome)).put(start, end, sitem);
                }
            } else if (!strand) {
                if (!this.positiveTreeTreatMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    this.positiveTreeTreatMap.put(chrome, potree);
                } else {
                    ((IntervalTree) this.positiveTreeTreatMap.get(chrome)).put(start, end, sitem);
                }
            }
        }
        itertemp.close();
    }

    public void initializeControlIntervalTree() {
        System.out.println("Initializing initializeControlIntervalTree  and site infomation");

        SAMRecordIterator itertemp = this.src.iterator();
        while (itertemp.hasNext()) {
            SAMRecord sitem = (SAMRecord) itertemp.next();
            if ((sitem.getReadUnmappedFlag())
                    || (sitem.getDuplicateReadFlag())
                    || (sitem.getNotPrimaryAlignmentFlag())
                    || (sitem.getReadFailsVendorQualityCheckFlag())) {
                continue;
            }
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();
            String chrome = sitem.getReferenceName();
            if (!this.chrlist.contains(chrome)) {
                continue;
            }
            if (strand) {
                if (!this.negativeTreeControlMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    this.negativeTreeControlMap.put(chrome, potree);
                } else {
                    ((IntervalTree) this.negativeTreeControlMap.get(chrome)).put(start, end, sitem);
                }
            } else if (!strand) {
                if (!this.positiveTreeControlMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    this.positiveTreeControlMap.put(chrome, potree);
                } else {
                    ((IntervalTree) this.positiveTreeControlMap.get(chrome)).put(start, end, sitem);
                }
            }
        }

        itertemp.close();
    }

    public void process() throws IOException {
        initializeTreatIntervalTree();
        initializeControlIntervalTree();
     
        try {
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator iterator = this.chrlist.iterator(); iterator.hasNext();) {
                String chr = (String) iterator.next();
                if ((this.positiveTreeTreatMap.get(chr) == null) && (this.positiveTreeControlMap.get(chr) == null)) {
                    continue;
                }
                runPSIseekerthread = new runPSIseekerThread(chr, false, this.positiveTreeTreatMap.get(chr),  this.positiveTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            for (Iterator iterator = this.chrlist.iterator(); iterator.hasNext();) {
                String chr = (String) iterator.next();
                if ((this.negativeTreeTreatMap.get(chr) == null) && (this.negativeTreeControlMap.get(chr) == null)) {
                    continue;
                }
                runPSIseekerthread = new runPSIseekerThread(chr, true,  this.negativeTreeTreatMap.get(chr), this.negativeTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            pool.shutdown();
            pool.awaitTermination(9223372036854775807L, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            Logger.getLogger(PSIseeker.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        CalculateClass cc=new CalculateClass(positiveResultMap, 
                negativeResultMap, 
                positiveBackGroundMap, 
                negativeBackGroundMap, 
                positiveGloableLamda, 
                negativeGloableLamda,
                positiveGloableNormalizedFactor,
                negativeGloableNormalizedFactor,
                Thread,  chrlist) ;
        
    }

    public synchronized void addnegativeResultMap(String chrom, HashSet<PSIout> psioutSET) {
        System.out.println("put\tN\t" + chrom + "\t" + psioutSET.size());
        this.negativeResultMap.put(chrom, psioutSET);
    }

    public synchronized void addnpositiveResultMap(String chrom, HashSet<PSIout> psioutSET) {
        System.out.println("put\tP\t" + chrom + "\t" + psioutSET.size());
        this.positiveResultMap.put(chrom, psioutSET);
    }

    public synchronized void addNegativeBackGroundMap(String chrom, IntervalTree backtree) {
//        System.out.println("Candidate sites at positiveB strand in \t" + chrom + "\t" + backtree.size());
        this.negativeBackGroundMap.put(chrom, backtree);
    }

    public synchronized void addPositiveBackGroundMap(String chrom, IntervalTree backtree) {
//        System.out.println("Candidate sites at positiveB strand in \t" + chrom + "\t" + backtree.size());
        this.positiveBackGroundMap.put(chrom, backtree);
    }
     public synchronized void addPositiveNormalizedFactor(String chrom,double factor){
         System.out.println("NormalizedFactor at positive strand in \t" + chrom + "\t" + factor);
        this.positiveGloableNormalizedFactor.put(chrom, factor);
    }
    public synchronized void addNegativeNormalizedFactor(String chrom,double factor){
        System.out.println("NormalizedFactor at Negative strand in \t" + chrom + "\t" + factor);
        this.negativeGloableNormalizedFactor.put(chrom, factor);
    }
    public synchronized void addPositiveGlobalLamda(String chrom,double lam){
         System.out.println("Lamda at positive strand in \t" + chrom + "\t" + lam);
        this.positiveGloableLamda.put(chrom, lam);
    }
    public synchronized void addNegativeGlobalLamda(String chrom,double lam){
        System.out.println("Lamda at Negative strand in \t" + chrom + "\t" + lam);
        this.negativeGloableLamda.put(chrom, lam);
    }
    
    
    
    public void print(String fileout) throws IOException {
        System.out.println("Writing out..");
        FileWriter fw = new FileWriter(fileout);
        ArrayList templist = new ArrayList();
        ArrayList plist = new ArrayList();
        fw.append("chr\tposition\texbase\tbase\trepresentSeq\tstrand\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\tenrichmentScore\n");
        for (Iterator chrit = this.chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            Iterator poit;
            if (this.negativeResultMap.get(chr) != null) {
                HashSet negpomap = (HashSet) this.negativeResultMap.get(chr);
                System.out.println("Write\tN\t" + chr + "\t" + negpomap.size());
                for (poit = negpomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
                    fw.append(pso.toString2() + "\n");
                    fw.flush();
                }
            }

            if (this.positiveResultMap.get(chr) != null) {
                HashSet pospomap = (HashSet) this.positiveResultMap.get(chr);
                System.out.println("Write\tP\t" + chr + "\t" + pospomap.size());
                for (poit = pospomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
                    fw.append(pso.toString2() + "\n");
                    fw.flush();
                }
            }
        }
        Iterator poit;
        fw.close();
    }

    public int getFilternumber() {
        return this.filternumber;
    }

    public void setFilternumber(int filternumber) {
        this.filternumber = filternumber;
    }

    public void fishertest(PSIout pis) {
        FisherExactTest test = new FisherExactTest();
        int a = pis.getSupporCountControl();
        int b = pis.getTotalCountControl();
        if (a == 0) {
            a = 1;
        }
        if (b == 0) {
            b = 1;
        }
        pis.setFisherPvalue(test.getTwoTailP(pis.getSupporCountInTreat(), pis.getTotalCountInTreat() - pis.getSupporCountInTreat(), a, b - a));
    }

    class runPSIseekerThread implements Runnable {

        private String subchr;
        private IntervalTree intervalTreat;
        private IntervalTree intervalControl;
        private boolean strand;

        public runPSIseekerThread(String subchr, boolean strand, IntervalTree intervalTreat, IntervalTree intervalControl) {
            this.subchr = subchr;
            this.intervalTreat = intervalTreat;
            this.intervalControl = intervalControl;
            this.strand = strand;
        }

        public void run() {
            HashSet psioutSet = new HashSet();
            HashSet<Integer> posset = new HashSet<Integer>();//storage background information
            //positve strand 
            double tamplam=0;
            if (strand) {
                for (Iterator poit = this.intervalTreat.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();

                    PSIout psi = new PSIout(subchr, sr.getAlignmentEnd(), this.strand);

                    if (!psioutSet.add(psi)) {
                        continue;
                    }
                    Iterator tempit1 = intervalTreat.overlappers(psi.getPosition(), psi.getPosition());
                    Iterator tempit2 = intervalControl.overlappers(psi.getPosition(), psi.getPosition());
                    CountReadsTreat(psi, strand, tempit1);
                    CountReadsControl(psi, strand, tempit2);

                    if ((psi.getSupporCountInTreat() < PSIseeker.this.filternumber) || (psi.getTotalCountControl() < PSIseeker.this.filternumber)) {
                        psioutSet.remove(psi);
                    } else if ((double)psi.getSupporCountInTreat() / (double)psi.getTotalCountInTreat() < PSIseeker.this.filterTreatRatio) {
                        psioutSet.remove(psi);
                    } else {
                        double a = psi.getSupporCountControl();
                        double b = psi.getTotalCountControl();
                        if (psi.getSupporCountControl() == 0) {
                            a = 1.0D;
                        }
                        if (psi.getTotalCountControl() == 0) {
                            b = 1.0D;
                        }
                        double tempa = (double)psi.getSupporCountInTreat() / (double)psi.getTotalCountInTreat() / (a / b);
                        psi.setEnrichmentScore(tempa);
                        if (tempa < PSIseeker.this.enrichmentThreshold) {
                            psioutSet.remove(psi);
                        } else {
                            psi.setFirstCigar(sr.getCigarString());
                            psi.setReadsString(sr.getReadString());

                            psi.setExbase(Character.toUpperCase(DNAsequenceProcess.getChargefBase(PSIseeker.this.Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition()))));
                        }
                    }
                }

                PSIseeker.this.addnegativeResultMap(subchr, psioutSet);
                
                //get background reads coverage and stored in a interval tree
                IntervalTree potree = new IntervalTree();
                for (Iterator poit = intervalControl.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    if (!posset.add(sr.getAlignmentEnd())) {
                        continue;
                    }
                    PSIout psi = new PSIout(subchr, sr.getAlignmentEnd());
                    Iterator<IntervalTree.Node> tempit2 = intervalControl.overlappers(sr.getAlignmentEnd(), sr.getAlignmentEnd());
                    CountReadsControl(psi, strand, tempit2);
                    if (psi.getTotalCountControl() < filternumber) {
                        continue;
                    }
                    if (((double) psi.getSupporCountControl() / (double) psi.getTotalCountControl()) < filterTreatRatio) {
                        continue;
                    }
                    potree.put(psi.getPosition(), psi.getPosition() + 1, psi);
//                    psi.calBackRatio();
                    tamplam+=psi.getSupporCountControl();
                }
//                intervalControl.clear();
                PSIseeker.this.addNegativeBackGroundMap(subchr, potree);
                double factor=(double)intervalTreat.size()/(double)intervalControl.size();
                PSIseeker.this.addNegativeNormalizedFactor(subchr,factor);
                PSIseeker.this.addNegativeGlobalLamda(subchr, tamplam/potree.size()*factor);
                
            } else if (!this.strand) {                             //positve strand 
                for (Iterator poit = this.intervalTreat.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    PSIout psi = new PSIout(this.subchr, sr.getAlignmentStart(), this.strand);
                    if (!psioutSet.add(psi)) {
                        continue;
                    }
                    Iterator tempit1 = this.intervalTreat.overlappers(psi.getPosition(), psi.getPosition());
                    Iterator tempit2 = this.intervalControl.overlappers(psi.getPosition(), psi.getPosition());
                    CountReadsTreat(psi, this.strand, tempit1);
                    CountReadsControl(psi, this.strand, tempit2);

                    if ((psi.getSupporCountInTreat() < PSIseeker.this.filternumber) || (psi.getTotalCountControl() < PSIseeker.this.filternumber)) {
                        psioutSet.remove(psi);
                    } else if ((double)psi.getSupporCountInTreat() / (double)psi.getTotalCountInTreat() < PSIseeker.this.filterTreatRatio) {
                        psioutSet.remove(psi);
                    } else {
                        double a = psi.getSupporCountControl();
                        double b = psi.getTotalCountControl();
                        if (psi.getSupporCountControl() == 0) {
                            a = 1.0D;
                        }
                        if (psi.getTotalCountControl() == 0) {
                            b = 1.0D;
                        }
                        double tempa = (double)psi.getSupporCountInTreat() / (double)psi.getTotalCountInTreat() / (a / b);
                        psi.setEnrichmentScore(tempa);
                        if (tempa < PSIseeker.this.enrichmentThreshold) {
                            psioutSet.remove(psi);
                        } else {
                            psi.setReadsString(sr.getReadString());
                            psi.setFirstCigar(sr.getCigarString());

                            psi.setExbase(Character.toUpperCase(PSIseeker.this.Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition() - 2)));
                        }
                    }
                }

                PSIseeker.this.addnpositiveResultMap(this.subchr, psioutSet);
                

                //get background reads coverage and stored in a interval tree
                IntervalTree potree = new IntervalTree();
                for (Iterator poit = intervalControl.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    if (!posset.add(sr.getAlignmentStart())) {
                        continue;
                    }
                    PSIout bt = new PSIout(subchr, sr.getAlignmentStart());
                    Iterator<IntervalTree.Node> tempit2 = intervalControl.overlappers(sr.getAlignmentEnd(), sr.getAlignmentEnd());
                    CountReadsControl(bt, strand, tempit2);
                    if (bt.getTotalCountControl() < filternumber) {
                        continue;
                    }
                    if (((double) bt.getSupporCountControl() / (double) bt.getTotalCountControl()) < filterTreatRatio) {
                        continue;
                    }
                    potree.put(bt.getPosition(), bt.getPosition() + 1, bt);
                    tamplam+=bt.getSupporCountControl();
                }
//                intervalControl.clear();
                PSIseeker.this.addPositiveBackGroundMap(subchr, potree);
                double factor=(double)intervalTreat.size()/(double)intervalControl.size();
                PSIseeker.this.addPositiveNormalizedFactor(subchr,factor);
                PSIseeker.this.addPositiveGlobalLamda(subchr, tamplam/potree.size()*factor);
                
            }
        }

        public void CountReadsTreat(PSIout po, boolean strand, Iterator<IntervalTree.Node> readsIt) {
            
            for (Iterator iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = (IntervalTree.Node) iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();

                if (strand) {
                    if (iscoverNegative(po.getPosition(), sr)) {
                        po.add1supporCountInTreat();
                    }
                } else if (iscoverPositve(po.getPosition(), sr)) {
                    po.add1supporCountInTreat();
                }
                if (isCoverAtMcigar(po.getPosition(), sr)) {
                    po.add1totalCountInTreat();
                }
            }
//            System.out.println(po.getSupporCountInTreat()+"\t"+po.getTotalCountInTreat());
        }

        public void CountReadsControl(PSIout po, boolean strand, Iterator<IntervalTree.Node> readsIt) {
            for (Iterator iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = (IntervalTree.Node) iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();
                if (strand) {
                    if (iscoverNegative(po.getPosition(), sr)) {
                        po.add1supporCountControl();
                    }
                } else if (iscoverPositve(po.getPosition(), sr)) {
                    po.add1supporCountControl();
                }
//                po.add1totalCountControl();
                if (isCoverAtMcigar(po.getPosition(), sr)) {
                    po.add1totalCountControl();
                }
            }
        }

        public boolean iscoverNegative(int pos, SAMRecord sr) {
            return sr.getAlignmentEnd() == pos;
        }

        public boolean iscoverPositve(int pos, SAMRecord sr) {
            return sr.getAlignmentStart() == pos;
        }

        public boolean isCoverAtMcigar(int pos, SAMRecord sr) {
            Cigar cgar = sr.getCigar();
            int currentlength = pos - sr.getAlignmentStart();
            List elist = cgar.getCigarElements();
            boolean iscover = true;
            for (Iterator it = elist.iterator(); it.hasNext();) {
                CigarElement ce = (CigarElement) it.next();
                int le = ce.getLength();
                String o = ce.getOperator().name();
                if ((currentlength < le) && (o.equals("M"))) {
                    break;
                }
                if ((currentlength < le) && (o.equals("N"))) {
                    iscover = false;
                    break;
                }
                currentlength -= le;
            }

            return iscover;
        }
    }

    
    
    public static void main(String[] args) throws IOException {
//        PSIseeker ps = new PSIseeker("E:\\迅雷下载\\3.bam", "E:\\迅雷下载\\4.bam", "E:\\迅雷下载\\dm6.fa");
        PSIseeker ps = new PSIseeker("/Users/zhaoqi/NetBeansProjects/SMULTQ02-3_chr4.bam", "/Users/zhaoqi/NetBeansProjects/SMULTQ02-4_chr4.bam", "/Users/zhaoqi/NetBeansProjects/dm6.fa");

        ps.process();
        ps.print("out.txt");
        
//           IntervalTree potree = new IntervalTree();
//        potree.put(1, 2, 3);
//        potree.put(4, 5, 5);
//        potree.put(6, 6, 6);
//        potree.put(7, 7, 7);
//        potree.put(8, 8, 8);
//        potree.overlappers(4, 4);
//        for (Iterator it = potree.overlappers(4, 8); it.hasNext();) {
//            IntervalTree.Node test = (IntervalTree.Node) it.next();
//            System.out.println(test.getValue());
//        }
        
        
    }

}
