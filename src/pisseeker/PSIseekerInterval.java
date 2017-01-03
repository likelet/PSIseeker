/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

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
import pisseeker.MultipleCorrection.FDR;
import pisseeker.pub.DNAsequenceProcess;
import pisseeker.pub.ToolsforCMD;

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
    public int Thread = 1;
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
            if (sitem.getReadUnmappedFlag()) continue;
            if (sitem.getDuplicateReadFlag()) continue;
            if (sitem.getNotPrimaryAlignmentFlag())continue;
            if (sitem.getReadFailsVendorQualityCheckFlag())continue;
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            treatLibsize++;
            if (strand) {
                if (!negativeTreeTreatMap.containsKey(chrome)) {
                    //samrecord
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    negativeTreeTreatMap.put(chrome, potree);
                    // candidate site
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    negativeResultMap.put(chrome, pomap);
                } else {
                    negativeTreeTreatMap.get(chrome).put(start, end, sitem);
                    // candidate site
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setReadsString(sitem.getReadString());
                    negativeResultMap.get(chrome).add(po);
                }
            } else if (!strand) {
                if (!positiveTreeTreatMap.containsKey(chrome)) {
                    IntervalTree potree = new IntervalTree();
                    potree.put(start, end, sitem);
                    positiveTreeTreatMap.put(chrome, potree);
                    //positive
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, start, strand);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    positiveResultMap.put(chrome, pomap);

                } else {
                    positiveTreeTreatMap.get(chrome).put(start, end, sitem);
                    //positive
                    PSIout po = new PSIout(chrome, start, strand);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setReadsString(sitem.getReadString());
                    positiveResultMap.get(chrome).add(po);
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

        System.out.println("Start Counting reads in parallel mode ");
        try {
            // run paralla
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);//Creat a new thread pool
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if(positiveTreeTreatMap.get(chr)==null)continue;
                if(positiveTreeControlMap.get(chr)==null)continue;
                runPSIseekerthread = new runPSIseekerThread(chr, negativeResultMap.get(chr), positiveTreeTreatMap.get(chr), positiveTreeControlMap.get(chr));
                pool.submit(runPSIseekerthread);
            }
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                if(negativeTreeTreatMap.get(chr)==null){
                    continue;
                }
                if(negativeTreeControlMap.get(chr)==null){
                    continue;
                }
                runPSIseekerthread = new runPSIseekerThread(chr, positiveResultMap.get(chr), negativeTreeTreatMap.get(chr), negativeTreeControlMap.get(chr));
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
        private HashSet<PSIout> pomap;
        private IntervalTree intervalTreat;
        private IntervalTree intervalControl;

        public runPSIseekerThread(String subchr, HashSet<PSIout> pomap, IntervalTree intervalTreat, IntervalTree intervalControl) {
            this.subchr = subchr;
            this.pomap = pomap;
            this.intervalTreat = intervalTreat;
            this.intervalControl = intervalControl;
        }

        @Override
        public void run() {
//            System.out.println("Start counting reads from " + ToolsforCMD.print_ansi_PURPLE(subchr));
            for (Iterator poit = pomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                Iterator<SAMRecord> tempit1 = intervalTreat.overlappers(pso.getPosition(), pso.getPosition());
                
                Iterator<SAMRecord> tempit2 = intervalControl.overlappers(pso.getPosition(), pso.getPosition());
                CountReadsTreat(pso, tempit1);
                CountReadsControl(pso, tempit2);
            }

        }

        //Count reads from treated bamfile
        public void CountReadsTreat(PSIout po, Iterator<SAMRecord> readsIt) {
            for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
                System.out.println("yes");
                SAMRecord sr = iterator.next();
                boolean strand = sr.getReadNegativeStrandFlag();
                if (strand == po.isStrandB()) {
                    if (strand) {
                        if (this.iscoverNegative(po.getPosition(), sr)) {
                            po.add1supporCountInTreat();
                        }
                        po.add1totalCountInTreat();

                    } else {
                        if (this.iscoverPositve(po.getPosition(), sr)) {
                            po.add1supporCountInTreat();
                        }
                        po.add1totalCountInTreat();
                    }

                }

            }
        }

        //Count reads from Control bamfile
        public void CountReadsControl(PSIout po, Iterator<SAMRecord> readsIt) {
            for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
                SAMRecord sr = iterator.next();
                boolean strand = sr.getReadNegativeStrandFlag();
                if (strand == po.isStrandB()) {
                    if (strand) {
                        if (this.iscoverNegative(po.getPosition(), sr)) {
                            po.add1supporCountControl();
                        }
                        po.add1totalCountControl();
                    } else {
                        if (this.iscoverPositve(po.getPosition(), sr)) {
                            po.add1supporCountControl();
                        }
                        po.add1totalCountControl();
                    }
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
    }

    public void print(String fileout) throws IOException {
        System.out.println("Writing out..");
        FileWriter fw = new FileWriter(fileout);
        ArrayList<PSIout> templist = new ArrayList<PSIout>();
        ArrayList<Double> plist = new ArrayList<Double>();
        fw.append("chr\tposition\tbase\tstrand\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\n");

        for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            if(negativeResultMap.get(chr)!=null){
                HashSet<PSIout> negpomap = negativeResultMap.get(chr);

                for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
                    //remove reads less than 2 
//                    if (pso.getSupporCountInTreat() < this.filternumber) {
//                        continue;
//                    }
//                    //remove reads 
//                    if (pso.getTotalCountControl() < this.filternumber) {
//                        continue;
//                    }
//                    if (pso.getSupporCountInTreat() / pso.getTotalCountInTreat() < pso.getSupporCountControl() / pso.getTotalCountControl()) {
//                        continue;
//                    }
                    pso.fishertest();
                    plist.add(pso.getPvalue());
                    templist.add(pso);
                }
            }
            if (positiveResultMap.get(chr) != null) {
                HashSet<PSIout> pospomap = positiveResultMap.get(chr);
                for (Iterator poit = pospomap.iterator(); poit.hasNext();) {
                    PSIout pso = (PSIout) poit.next();
//                    if (pso.getSupporCountInTreat() < this.filternumber) {
//                        continue;
//                    }
//                    //remove reads 
//                    if (pso.getTotalCountControl() < this.filternumber) {
//                        continue;
//                    }
//                    // remove site show lower propotion in treatment
//                    if (pso.getSupporCountInTreat() / pso.getTotalCountInTreat() < pso.getSupporCountControl() / pso.getTotalCountControl()) {
//                        continue;
//                    }
                    pso.fishertest();
                    plist.add(pso.getPvalue());
                    templist.add(pso);
                }
            }
        }
        //FDR calculation && write out
        ArrayList<Double> fdrlist = new FDR(plist).getFdrDoublelist();
        for (int i = 0; i < templist.size(); i++) {
            PSIout psi = templist.get(i);

            if (psi.getStrand().endsWith("+")) {
                psi.setBase(Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition() - 2));
            } else {
                psi.setBase(Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition()));
            }
            psi.setAdjustP(fdrlist.get(i));
            fw.append(psi.toString() + "\t" + "\n");
        }
        fw.flush();
        fw.close();
    }

    public int getFilternumber() {
        return filternumber;
    }

    public void setFilternumber(int filternumber) {
        this.filternumber = filternumber;
    }

    public static void main(String[] args) throws IOException {
//        PSIseekerInterval ps = new PSIseekerInterval("E:\\迅雷下载\\SMULTQ02-3_chr4.bam", "E:\\迅雷下载\\SMULTQ02-4_chr4.bam", "E:\\迅雷下载\\dm6.fa");
//        ps.process();
//        ps.print("out.txt");
    
        IntervalTree potree = new IntervalTree();
        potree.put(1, 2, 3);
        potree.put(4, 5, 4);
        potree.overlappers(4, 4);
        for (Iterator it =potree.overlappers(4, 4); it.hasNext();) {
            
            System.out.println(it.next());
        }
        

    }
}
