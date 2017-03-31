/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
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

/**
 *
 * @author Administrator
 * @since 2016-12-27
 * @coding time 16:57:48
 * @author Qi Zhao
 */
public final class PSIseekerInterval {

    private final SamReader srt;
    private final SamReader src;
    private final IndexedFastaSequenceFile Indexgenomefile;

    //storage chrome and position information
    private HashMap<String, HashSet<PSIout>> positiveResultMap = new HashMap<String, HashSet<PSIout>>();// result in positve strand 
    private HashMap<String, HashSet<PSIout>> negativeResultMap = new HashMap<String, HashSet<PSIout>>();// result in negative strand;
    
    private HashMap<String, IntervalTree>  positiveBackGroundMap = new HashMap<String, IntervalTree>();// result in positve strand 
    private HashMap<String, IntervalTree>  negativeBackGroundMap = new HashMap<String, IntervalTree>();// result in negative strand;
    
    //storage reads
    private HashMap<String, IntervalTree> positiveTreeControlMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> negativeTreeControlMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> positiveTreeTreatMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> negativeTreeTreatMap = new HashMap<String, IntervalTree>();
    public int filternumber = 2;//at least 2 read support
    public double filterTreatRatio = 2;//at least 2 read support
    public double enrichmentThreshold = 1;//at least fold change
    public int Thread=3;// parellel compute thread bumber 
    public HashSet<String> chrlist = new HashSet<String>();// storage chrlist infromation

    /**
     * @param args the command line arguments
     */
    //tbam and cbam are both sorted bamfile 
    public PSIseekerInterval(String tbam, String cbam, String genomefile) throws IOException {
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        this.Indexgenomefile = new IndexedFastaSequenceFile(new File(genomefile));;
        srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile1).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
        src = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile2).
                index(new File(bamfile2.getAbsolutePath() + ".bai"))
        );

    }

    public void initializeTreatIntervalTree() {
        //first cycle initialized chromesome and positions
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
            if (sitem.getReadUnmappedFlag()) {
                continue;
            }
            if (sitem.getDuplicateReadFlag()) {
                continue;
            }
            if (sitem.getNotPrimaryAlignmentFlag()) {
                continue;
            }
            if (sitem.getReadFailsVendorQualityCheckFlag()) {
                continue;
            }
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            if (strand) {
                if (!negativeTreeTreatMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    negativeTreeTreatMap.put(chrome, potree);
                } else {
                    negativeTreeTreatMap.get(chrome).put(start, end, sitem);
                }
            } else if (!strand) {
                if (!positiveTreeTreatMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    positiveTreeTreatMap.put(chrome, potree);

                } else {
                    positiveTreeTreatMap.get(chrome).put(start, end, sitem);
                }
            }
        }
        itertemp.close();
    }

    public void initializeControlIntervalTree() {
        System.out.println("Initializing initializeControlIntervalTree  and site infomation");
        //first cycle initialized chromesome and positions
        SAMRecordIterator itertemp = src.iterator();
        while (itertemp.hasNext()) {
            SAMRecord sitem = (SAMRecord) itertemp.next();
            if (sitem.getReadUnmappedFlag()) {
                continue;
            }
            if (sitem.getDuplicateReadFlag()) {
                continue;
            }
            if (sitem.getNotPrimaryAlignmentFlag()) {
                continue;
            }
            if (sitem.getReadFailsVendorQualityCheckFlag()) {
                continue;
            }
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            if (!this.chrlist.contains(chrome)) {
                continue;
            }
            if (strand) {
                if (!negativeTreeControlMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    negativeTreeControlMap.put(chrome, potree);
                } else {
//                    po.setExbase(Indexgenomefile.getSequence(chrome).getBaseString().charAt(end+1));
                    negativeTreeControlMap.get(chrome).put(start, end, sitem);
                }
            } else if (!strand) {
                if (!positiveTreeControlMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    positiveTreeControlMap.put(chrome, potree);
                } else {
                    positiveTreeControlMap.get(chrome).put(start, end, sitem);

                }
            }
        }
        itertemp.close();
    }

    //run analysis parallel
    public void process() throws IOException {
        this.initializeTreatIntervalTree();
        this.initializeControlIntervalTree();

//        System.out.println("Start Counting reads in parallel mode ");
        try {
            // run paralla
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);//Creat a new thread pool
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if (positiveTreeTreatMap.get(chr) == null && positiveTreeControlMap.get(chr) == null) {
                    continue;
                }
                runPSIseekerthread = new runPSIseekerThread(chr, false, positiveTreeTreatMap.get(chr), positiveTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if (negativeTreeTreatMap.get(chr) == null && negativeTreeControlMap.get(chr) == null) {
                    continue;
                }
                runPSIseekerthread = new runPSIseekerThread(chr, true, negativeTreeTreatMap.get(chr), negativeTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            // write out 
        } catch (InterruptedException ex) {
            Logger.getLogger(PSIseekerInterval.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    class runPSIseekerThread implements Runnable {

        private String subchr;
        private IntervalTree intervalTreat;
        private IntervalTree intervalControl;
        private boolean strand; // true for negative; true for positive

        public runPSIseekerThread(String subchr, boolean strand, IntervalTree intervalTreat, IntervalTree intervalControl) {
            this.subchr = subchr;
            this.intervalTreat = intervalTreat;
            this.intervalControl = intervalControl;
            this.strand = strand;
        }

        @Override
        public void run() {
            HashSet<PSIout> psioutSet = new HashSet<PSIout>();
            HashSet<Integer> possetS = new HashSet<Integer>();
            HashSet<Integer> posset = new HashSet<Integer>();//storage background information
//            System.out.println("Start counting reads from " + ToolsforCMD.print_ansi_PURPLE(subchr));

            if (strand) {
                for (Iterator poit = this.intervalTreat.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();

                    PSIout psi = new PSIout(this.subchr, sr.getAlignmentEnd(), this.strand);

                    if (!psioutSet.add(psi)) {
                        continue;
                    }
                    Iterator tempit1 = this.intervalTreat.overlappers(psi.getPosition(), psi.getPosition());
                    Iterator tempit2 = this.intervalControl.overlappers(psi.getPosition(), psi.getPosition());
                    CountReadsTreat(psi, this.strand, tempit1);
                    CountReadsControl(psi, this.strand, tempit2);

                    if ((psi.getSupporCountInTreat() < PSIseekerInterval.this.filternumber) || (psi.getTotalCountControl() < PSIseekerInterval.this.filternumber)) {
                        psioutSet.remove(psi);
                    } else if (psi.getSupporCountInTreat() / psi.getTotalCountInTreat() < PSIseekerInterval.this.filterTreatRatio) {
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
                        double tempa = psi.getSupporCountInTreat() / psi.getTotalCountInTreat() / (a / b);
                        psi.setEnrichmentScore(tempa);
                        if (tempa < PSIseekerInterval.this.enrichmentThreshold) {
                            psioutSet.remove(psi);
                        } else {
                            psi.setFirstCigar(sr.getCigarString());
                            psi.setReadsString(sr.getReadString());

                            psi.setExbase(Character.toUpperCase(DNAsequenceProcess.getChargefBase(PSIseekerInterval.this.Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition()))));
                        }
                    }
                }

                PSIseekerInterval.this.addnegativeResultMap(this.subchr, psioutSet);
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
                }
//                intervalControl.clear();
                PSIseekerInterval.this.addNegativeBackGroundMap(subchr, potree);
            } else if (!strand) {
//                System.out.println(intervalTreat.size());
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

                    if ((psi.getSupporCountInTreat() < PSIseekerInterval.this.filternumber) || (psi.getTotalCountControl() < PSIseekerInterval.this.filternumber)) {
                        psioutSet.remove(psi);
                    } else if (psi.getSupporCountInTreat() / psi.getTotalCountInTreat() < PSIseekerInterval.this.filterTreatRatio) {
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
                        double tempa = psi.getSupporCountInTreat() / psi.getTotalCountInTreat() / (a / b);
                        psi.setEnrichmentScore(tempa);
                        if (tempa < PSIseekerInterval.this.enrichmentThreshold) {
                            psioutSet.remove(psi);
                        } else {
                            psi.setReadsString(sr.getReadString());
                            psi.setFirstCigar(sr.getCigarString());

                            psi.setExbase(Character.toUpperCase(PSIseekerInterval.this.Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition() - 2)));
                        }
                    }
                }

                PSIseekerInterval.this.addnpositiveResultMap(this.subchr, psioutSet);

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
                }
//                intervalControl.clear();
                PSIseekerInterval.this.addPositiveBackGroundMap(subchr, potree);
                
            }

        }

        //Count reads from treated bamfile
        public void CountReadsTreat(PSIout po, boolean strand, Iterator<IntervalTree.Node> readsIt) {

            for (Iterator<IntervalTree.Node> iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();
//                if(po.getPosition()==414041){
//                        System.out.println(sr.getAlignmentStart()+"\t"+sr.getAlignmentEnd()+"\t"+sr.getCigarString()+"\t"+isCoverAtMcigar(po.getPosition(), sr));
//                    }
                if (strand) {
                    if (this.iscoverNegative(po.getPosition(), sr)) {
                        po.add1supporCountInTreat();
                    }
                } else if (this.iscoverPositve(po.getPosition(), sr)) {
                    po.add1supporCountInTreat();
                }
                //exclude reads gap span position like "20M2000N20M "
                if (isCoverAtMcigar(po.getPosition(), sr)) {
                    po.add1totalCountInTreat();
                } else {
//                        System.out.println(po.getPosition()+"\t"+sr.getCigarString());
                }
            }

        }

        //Count reads from Control bamfile
        public void CountReadsControl(PSIout po, boolean strand, Iterator<IntervalTree.Node> readsIt) {
            for (Iterator<IntervalTree.Node> iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();
                if (strand) {
                    if (this.iscoverNegative(po.getPosition(), sr)) {
                        po.add1supporCountControl();
                    }
                } else if (this.iscoverPositve(po.getPosition(), sr)) {
                    po.add1supporCountControl();
                }
                //exclude reads gap span position like "20M2000N20M "
//                po.add1totalCountControl();
                if (isCoverAtMcigar(po.getPosition(), sr)) {
                    po.add1totalCountControl();
                }
            }

        }

        public boolean iscoverNegative(int pos, SAMRecord sr) {
            if (sr.getAlignmentEnd() == pos) {
                return true;
            } else {
                return false;
            }
        }

        public boolean iscoverPositve(int pos, SAMRecord sr) {
            if (sr.getAlignmentStart() == pos) {
                return true;
            } else {
                return false;
            }
        }

        public boolean isCoverAtMcigar(int pos, SAMRecord sr) {
            Cigar cgar = sr.getCigar();
            int currentlength = pos - sr.getAlignmentStart();
            List<CigarElement> elist = cgar.getCigarElements();
            boolean iscover = true;
            for (Iterator it = elist.iterator(); it.hasNext();) {
                CigarElement ce = (CigarElement) it.next();
                int le = ce.getLength();
                String o = ce.getOperator().name();
                if (currentlength <le && o.equals("M")) {
                    break;
                } else if (currentlength <le && o.equals("N")) {
                    iscover = false;
                    break;
                } else {
                    currentlength = currentlength - le;
                }
            }
//            System.out.println("finish");
            return iscover;
        }

    }

    //for parallel synchronized method 
    public synchronized void addnegativeResultMap(String chrom, HashSet<PSIout> psioutSET) {
        System.out.println("Candidate sites at positive strand in \t" + chrom + "\t" + psioutSET.size());
        this.negativeResultMap.put(chrom, psioutSET);
    }

    public synchronized void addnpositiveResultMap(String chrom, HashSet<PSIout> psioutSET) {
        System.out.println("Candidate sites at positive strand in \t" + chrom + "\t" + psioutSET.size());
        this.positiveResultMap.put(chrom, psioutSET);
    }

    public synchronized void addNegativeBackGroundMap(String chrom, IntervalTree backtree) {
        System.out.println("Candidate sites at positiveB strand in \t" + chrom + "\t" + backtree.size());
        this.negativeBackGroundMap.put(chrom, backtree);
    }

    public synchronized void addPositiveBackGroundMap(String chrom, IntervalTree backtree) {
        System.out.println("Candidate sites at positiveB strand in \t" + chrom + "\t" + backtree.size());
        this.positiveBackGroundMap.put(chrom, backtree);
    }
    
    
    public void print(String fileout) throws IOException {
        System.out.println("Writing out..");
        FileWriter fw = new FileWriter(fileout);
        ArrayList<PSIout> templist = new ArrayList<PSIout>();
        ArrayList<Double> plist = new ArrayList<Double>();
        fw.append("chr\tposition\texbase\tbase\trepresentSeq\tstrand\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\tenrichmentScore\n");
        for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            if (negativeResultMap.get(chr) != null) {

                HashSet<PSIout> negpomap = negativeResultMap.get(chr);
                for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
//                    fishertest(pso);
                    fw.append(pso.toString2() + "\n");
                    fw.flush();
//                    templist.add(pso);
                }
            }
            if (positiveResultMap.get(chr) != null) {
                HashSet<PSIout> pospomap = positiveResultMap.get(chr);
                for (Iterator poit = pospomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
//                    fishertest(pso);
                    fw.append(pso.toString2() + "\n");
                    fw.flush();
                }
            }
        }
        //FDR calculation && write out
//        ArrayList<Double> fdrlist = new FDR(plist).getFdrDoublelist();

        fw.close();
    }
    // why not work？
    public void fishertest(PSIout pis) {
        FisherExactTest test = new FisherExactTest();
        int a = pis.getSupporCountControl();
        if (a == 0) {
            a = 1;
        }
        if(pis.getTotalCountInTreat() - pis.getSupporCountInTreat()<0){
            System.out.println(pis.getPosition());
        }
        pis.setFisherPvalue(test.getTwoTailP(pis.getSupporCountInTreat(), pis.getTotalCountInTreat() - pis.getSupporCountInTreat(), a, pis.getTotalCountControl() - a));
    }

    public HashMap<String, HashSet<PSIout>> getPositiveResultMap() {
        return positiveResultMap;
    }

    public HashMap<String, HashSet<PSIout>> getNegativeResultMap() {
        return negativeResultMap;
    }

    
    
    
    
    public static void main(String[] args) throws IOException {
//       PSIseekerInterval ps = new PSIseekerInterval("E:\\迅雷下载\\3.bam", "E:\\迅雷下载\\4.bam", "E:\\迅雷下载\\dm6.fa");
        PSIseekerInterval ps = new PSIseekerInterval("E:\\迅雷下载\\SMULTQ02-3_chr4.bam", "E:\\迅雷下载\\SMULTQ02-4_chr4.bam", "E:\\迅雷下载\\dm6.fa");
        ps.process();
        ps.print("E:\\迅雷下载\\out.txt");
    }
}
