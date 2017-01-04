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
    public int Thread = 3;
    public String tbam;
    public String cbam;
    public IndexedFastaSequenceFile Indexgenomefile;
    private double controlLibsize = 0;
    private double treatLibsize = 0;

    //storage chrome and position information
    private HashMap<String, HashSet<PSIout>> positiveResultMap = new HashMap<String, HashSet<PSIout>>();// result in positve strand 
    private HashMap<String, HashSet<PSIout>> negativeResultMap = new HashMap<String, HashSet<PSIout>>();// result in negative strand;

    private HashMap<String, IntervalTree> positiveTreeControlMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> negativeTreeControlMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> positiveTreeTreatMap = new HashMap<String, IntervalTree>();
    private HashMap<String, IntervalTree> negativeTreeTreatMap = new HashMap<String, IntervalTree>();
//    private HashMap<Integer,PSIout> pomap=new HashMap();
    public int filternumber = 2;//at least 2 read support
    public HashSet<String> chrlist = new HashSet<String>();

    /**
     * @param args the command line arguments
     */
    //tbam and cbam are both sorted bamfile 
    public PSIseekerInterval(String tbam, String cbam, String genomefile) throws IOException {
        this.cbam = cbam;
        this.tbam = tbam;
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

        this.initializeTreatIntervalTree();
        this.initializeControlIntervalTree();
        this.process();

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
             if (sitem.getReadUnmappedFlag())  continue;
            if (sitem.getDuplicateReadFlag())  continue;
            if (sitem.getNotPrimaryAlignmentFlag())  continue;
            if (sitem.getReadFailsVendorQualityCheckFlag()) continue;
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            treatLibsize++;
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
            if (sitem.getReadUnmappedFlag())   continue;
            if (sitem.getDuplicateReadFlag())  continue;
            if (sitem.getNotPrimaryAlignmentFlag())   continue;
            if (sitem.getReadFailsVendorQualityCheckFlag())  continue;
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            if (!this.chrlist.contains(chrome)) {
                continue;
            }
            controlLibsize++;
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

//        System.out.println("Start Counting reads in parallel mode ");
        try {
            // run paralla
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);//Creat a new thread pool
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if (positiveTreeTreatMap.get(chr) == null&&positiveTreeControlMap.get(chr) == null)   continue;
                runPSIseekerthread = new runPSIseekerThread(chr, false,positiveTreeTreatMap.get(chr), positiveTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if (negativeTreeTreatMap.get(chr) == null&&negativeTreeControlMap.get(chr) == null)  continue;
                runPSIseekerthread = new runPSIseekerThread(chr, true,negativeTreeTreatMap.get(chr), negativeTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            // write out 
        } catch (InterruptedException ex) {
            Logger.getLogger(PSIseekerParallel.class.getName()).log(Level.SEVERE, null, ex);
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
            this.strand=strand;
        }

        @Override
        public void run() {
            HashSet<PSIout> psioutSet=new HashSet<PSIout> ();
//            System.out.println("ok");
//            System.out.println("Start counting reads from " + ToolsforCMD.print_ansi_PURPLE(subchr));
            if (strand) {
                for (Iterator poit = intervalTreat.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    PSIout psi = new PSIout(subchr, sr.getAlignmentEnd(), strand);
                    if (!psioutSet.add(psi)) continue;
                    Iterator<IntervalTree.Node> tempit1 = intervalTreat.overlappers(psi.getPosition(), psi.getPosition());
                    Iterator<IntervalTree.Node> tempit2 = intervalControl.overlappers(psi.getPosition(), psi.getPosition());
                    CountReadsTreat(psi, tempit1);
                    CountReadsControl(psi, tempit2);
                    if (psi.getSupporCountInTreat() < filternumber || psi.getTotalCountControl() < filternumber) {
                        psioutSet.remove(psi);
                    }else{
                        psi.setBase(new DNAsequenceProcess().getReverseComplimentary(sr.getReadString()).charAt(0));
                        psi.setReadsString(sr.getReadString());
                        psi.fishertest();
                    }
                }
//                System.out.println("ADD N");
                addnegativeResultMap(subchr,psioutSet);
            } else {
                for (Iterator poit = intervalTreat.iterator(); poit.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) poit.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    PSIout psi = new PSIout(subchr, sr.getAlignmentStart(), strand);
                    if (!psioutSet.add(psi)) continue;
                    Iterator<IntervalTree.Node> tempit1 = intervalTreat.overlappers(psi.getPosition(), psi.getPosition());
                    Iterator<IntervalTree.Node> tempit2 = intervalControl.overlappers(psi.getPosition(), psi.getPosition());
                    CountReadsTreat(psi, tempit1);
                    CountReadsControl(psi, tempit2);
                    if (psi.getSupporCountInTreat() < filternumber || psi.getTotalCountControl() < filternumber) {
                        psioutSet.remove(psi);
                    }else{
                        psi.setBase(sr.getReadString().charAt(0));
                        psi.setReadsString(sr.getReadString());
                        psi.fishertest();
                    }
                }
                addnpositiveResultMap(subchr,psioutSet);
            }

        }

        //Count reads from treated bamfile
        public void CountReadsTreat(PSIout po, Iterator<IntervalTree.Node> readsIt) {
//            System.out.println("yes");
            for (Iterator<IntervalTree.Node> iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();
                boolean strand = sr.getReadNegativeStrandFlag();
                if (strand) {
                    if (this.iscoverNegative(po.getPosition(), sr)) 
                        po.add1supporCountInTreat();
                } else {
                    if (this.iscoverPositve(po.getPosition(), sr)) 
                        po.add1supporCountInTreat();
                }
                 //exclude reads gap span position like "20M2000N20M "
                if (isCoverAtMcigar(po.getPosition(), sr)) {
                    po.add1totalCountInTreat();
                }
            }
            
            
            
        }

        //Count reads from Control bamfile
        public void CountReadsControl(PSIout po, Iterator<IntervalTree.Node> readsIt) {
//            System.out.println("NO");
            for (Iterator<IntervalTree.Node> iterator = readsIt; iterator.hasNext();) {
                IntervalTree.Node nd = iterator.next();
                SAMRecord sr = (SAMRecord) nd.getValue();
                boolean strand = sr.getReadNegativeStrandFlag();
                    if (strand) {
                        if (this.iscoverNegative(po.getPosition(), sr)) {
                            po.add1supporCountControl();
                        }
                    } else {
                        if (this.iscoverPositve(po.getPosition(), sr)) {
                            po.add1supporCountControl();
                        }
                    }
                  //exclude reads gap span position like "20M2000N20M "
                po.add1totalCountControl();
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
                int l = ce.getLength();
                String o = ce.getOperator().name();
                if (currentlength <= l && o.equals("M")) {
                    break;
                } else if (currentlength <= l && !o.equals("M")) {
                    iscover = false;
                    break;
                }

            }
            return iscover;
        }
    }

    //for parallel
    public synchronized void  addnegativeResultMap(String chrom,HashSet<PSIout> psioutSET ){
        System.out.println("put\tN\t"+chrom+"\t"+psioutSET.size() );
        this.negativeResultMap.put(chrom, psioutSET);
    }
    public synchronized void  addnpositiveResultMap(String chrom,HashSet<PSIout> psioutSET ){
        System.out.println("put\tP\t"+chrom+"\t"+psioutSET.size() );
        this.positiveResultMap.put(chrom, psioutSET);
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
                System.out.println("Write\tN\t" + chr + "\t" + negpomap.size());
                for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
                    fw.append(pso.toString() + "\n");
                    fw.flush();
//                    templist.add(pso);
                }
            }
            if (positiveResultMap.get(chr) != null) {
                HashSet<PSIout> pospomap = positiveResultMap.get(chr);
                System.out.println("Write\tP\t"+chr+"\t"+pospomap.size() );
                for (Iterator poit = pospomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
                    fw.append(pso.toString() + "\n");
                    fw.flush();
                }
            }
        }
        //FDR calculation && write out
//        ArrayList<Double> fdrlist = new FDR(plist).getFdrDoublelist();
      

        fw.close();
    }

    public int getFilternumber() {
        return filternumber;
    }

    public void setFilternumber(int filternumber) {
        this.filternumber = filternumber;
    }

    public static void main(String[] args) throws IOException {
        PSIseekerInterval ps = new PSIseekerInterval("E:\\迅雷下载\\3.bam", "E:\\迅雷下载\\4.bam", "E:\\迅雷下载\\dm6.fa");
        ps.print("out.txt");

//        IntervalTree potree = new IntervalTree();
//        potree.put(1, 2, 3);
//        potree.put(4, 5, 4);
//        potree.overlappers(4, 4);
//        for (Iterator it =potree.overlappers(4, 4); it.hasNext();) {
//            IntervalTree.Node test=(IntervalTree.Node) it.next();
//            System.out.println(test.getValue());
//        }
    }
}
