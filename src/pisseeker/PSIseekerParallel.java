/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

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
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import pisseeker.MultipleCorrection.FDR;
import pisseeker.pub.ToolsforCMD;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import pisseeker.pub.DNAsequenceProcess;

/**
 *
 * @author zhaoqi
 */
public class PSIseekerParallel {

    private final SamReader srt;
    public  int Thread=1;
    public String tbam;
    public String cbam;
    public String genomefile;
//    private libsize
    
    //storage chrome and position information
    private HashMap<String, HashSet<PSIout>> positiveResultMap = new HashMap<String, HashSet<PSIout>>();// result in positve strand 
    private HashMap<String, HashSet<PSIout>> negativeResultMap = new HashMap<String, HashSet<PSIout>>();// result in negative strand;
//    private HashMap<Integer,PSIout> pomap=new HashMap();
    public int filternumber = 5;//at least 2 read support
    public HashSet<String> chrlist = new HashSet<String>();

    /**
     * @param args the command line arguments
     */
    //tbam and cbam are both sorted bamfile 
    
    
    public PSIseekerParallel(String tbam, String cbam,String genomefile) throws IOException {
        this.cbam=cbam;
        this.tbam=tbam;
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        this.genomefile=genomefile;
        srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile1).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
       
        this.initialize();
        
        
        
    }
    public PSIseekerParallel(String tbam, String cbam,String genomefile, String out) throws IOException {
        this.cbam=cbam;
        this.tbam=tbam;
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        this.genomefile=genomefile;
        srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile1).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
    
        this.initialize();
    }
    

    public void initialize(){
        // in one specific chrome
        ArrayList<SAMRecord> samlistTreat = new ArrayList<SAMRecord>();
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
       
        System.out.println("Initializing candidate site infomation");
       
        int count=0;
        while (itertemp.hasNext()) {
            SAMRecord sitem = (SAMRecord) itertemp.next();
            if (sitem.getReadUnmappedFlag())continue;
            if (sitem.getDuplicateReadFlag()) continue;
            if (sitem.getNotPrimaryAlignmentFlag()) continue;
            if (sitem.getReadFailsVendorQualityCheckFlag()) continue;
            int end = sitem.getAlignmentEnd();
            int start = sitem.getAlignmentStart();
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand)
            String chrome = sitem.getReferenceName();
            count++;
            if (strand) {
                if (!negativeResultMap.containsKey(chrome)) {
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    negativeResultMap.put(chrome, pomap);
                } else {
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setReadsString(sitem.getReadString());
                    negativeResultMap.get(chrome).add(po);
                }
            } else if (!strand) {
                if (!positiveResultMap.containsKey(chrome)) {
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, start, strand);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    positiveResultMap.put(chrome, pomap);
                } else {
                    PSIout po = new PSIout(chrome, start, strand);
//                    System.out.println(chrome);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setReadsString(sitem.getReadString());
                    positiveResultMap.get(chrome).add(po);
                }
            }
        }
        itertemp.close();
        
    }
    
    //run analysis parallel
    public void process() throws IOException {

        System.out.println("Start Counting reads in parrallel mode ");
        try {
            // run paralla
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);//Creat a new thread pool
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                runPSIseekerthread = new runPSIseekerThread(chr, negativeResultMap.get(chr), tbam, cbam);
                pool.submit(runPSIseekerthread);
            }
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                runPSIseekerthread = new runPSIseekerThread(chr, positiveResultMap.get(chr), tbam, cbam);
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
        private SamReader subsrt;
        private SamReader subsrc;

        public runPSIseekerThread(String subchr, HashSet<PSIout> pomap, String tbam, String cbam) {
            this.subchr = subchr;
            this.pomap = pomap;
            File bamfile1 = new File(tbam);
            File bamfile2 = new File(cbam);
            this.subsrt = SamReaderFactory.makeDefault().open(
                    SamInputResource.of(bamfile1).
                    index(new File(bamfile1.getAbsolutePath() + ".bai"))
            );
            this.subsrc = SamReaderFactory.makeDefault().open(
                    SamInputResource.of(bamfile2).
                    index(new File(bamfile2.getAbsolutePath() + ".bai"))
            );
        }

        @Override
        public void run() {
            System.out.println("Start counting reads from " + ToolsforCMD.print_ansi_PURPLE(subchr));
            for (Iterator poit = pomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                SAMRecordIterator tempit1 = subsrt.queryOverlapping(subchr, pso.getPosition(), pso.getPosition());
                SAMRecordIterator tempit2 = subsrc.queryOverlapping(subchr, pso.getPosition(), pso.getPosition());
                CountReadsTreat(pso, tempit1);
                CountReadsControl(pso, tempit2);
            }
            
            try {
                subsrt.close();
                subsrc.close();
            } catch (IOException ex) {
                Logger.getLogger(PSIseekerParallel.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        //Count reads from treated bamfile
        public void CountReadsTreat(PSIout po, SAMRecordIterator readsIt) {
            for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
                SAMRecord sr = iterator.next();
               //exclude bad reads
                if (sr.getReadUnmappedFlag())continue;
                if (sr.getDuplicateReadFlag()) continue;
                if (sr.getNotPrimaryAlignmentFlag()) continue;
                if (sr.getReadFailsVendorQualityCheckFlag()) continue;
                
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
            readsIt.close();
        }

        //Count reads from Control bamfile
        public void CountReadsControl(PSIout po, SAMRecordIterator readsIt) {
            for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
                SAMRecord sr = iterator.next();
                if (sr.getReadUnmappedFlag())continue;
                if (sr.getDuplicateReadFlag()) continue;
                if (sr.getNotPrimaryAlignmentFlag()) continue;
                if (sr.getReadFailsVendorQualityCheckFlag()) continue;
                
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
            readsIt.close();
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
        fw.append("chr\tposition\texbase\tbase\trepresentSeq\tstrand\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\tenrichmentScore\n");

        for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            HashSet<PSIout> negpomap = negativeResultMap.get(chr);
            if(negpomap==null) continue;
            for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                //remove reads less than 2 
                if (pso.getSupporCountInTreat() < this.filternumber) {
                    continue;
                }
                //remove reads 
                if (pso.getTotalCountControl() < this.filternumber) {
                    continue;
                }
                 pso.fishertest();
                double a=pso.getSupporCountControl();
                double b=pso.getTotalCountControl();
                if (pso.getSupporCountControl() == 0) {
                    a = 1;
                }
                if (pso.getTotalCountControl() == 0) {
                    b = 1;
                }
//                System.out.println(pso.toString2());
                double tempa=((double)pso.getSupporCountInTreat() / (double)pso.getTotalCountInTreat()) / (a / b);
//                System.out.println(tempa);
                pso.setEnrichmentScore(tempa);
                
                plist.add(pso.getPvalue());
                templist.add(pso);
            }
            HashSet<PSIout> pospomap = positiveResultMap.get(chr);
            if(negpomap==null) continue;
            for (Iterator poit = pospomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                if (pso.getSupporCountInTreat() < this.filternumber) {
                    continue;
                }
                //remove reads 
                if (pso.getTotalCountControl() < this.filternumber) {
                    continue;
                }
                // remove site show lower propotion in treatment
//                if (pso.getSupporCountInTreat() / pso.getTotalCountInTreat() < pso.getSupporCountControl() / pso.getTotalCountControl()) {
//                    continue;
//                }
                // avoiding zero devider
                pso.fishertest();
                double a=pso.getSupporCountControl();
                double b=pso.getTotalCountControl();
                if (pso.getSupporCountControl() == 0) {
                    a = 1;
                }
                if (pso.getTotalCountControl() == 0) {
                    b = 1;
                }
//                System.out.println(pso.toString2());
                double tempa=((double)pso.getSupporCountInTreat() / (double)pso.getTotalCountInTreat()) / (a / b);
//                System.out.println(tempa);
                pso.setEnrichmentScore(tempa);
                

                plist.add(pso.getPvalue());
                templist.add(pso);
            }
        }
        //FDR calculation && write out
        IndexedFastaSequenceFile Indexgenomefile=new IndexedFastaSequenceFile (new File(genomefile));;
        ArrayList<Double> fdrlist = new FDR(plist).getFdrDoublelist();
        for (int i = 0; i < templist.size(); i++) {
            PSIout psi = templist.get(i);
            
            if(psi.getStrand().endsWith("+")){
                psi.setExbase(Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition()-2));
            }else{
                psi.setExbase(DNAsequenceProcess.getChargefBase(Indexgenomefile.getSequence(psi.getChr()).getBaseString().charAt(psi.getPosition())));
            }
            psi.setAdjustP(fdrlist.get(i));
            fw.append(psi.toString2() + "\t" + "\n");
        }
        fw.flush();
        fw.close();
        try {
            //release  memory resources
            Indexgenomefile.close();
        } catch (IOException ex) {
            Logger.getLogger(PSIseekerParallel.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public int getFilternumber() {
        return filternumber;
    }

    public void setFilternumber(int filternumber) {
        this.filternumber = filternumber;
    }

    public static void main(String[] args) throws IOException {
//        PSIseekerParallel ps=new PSIseekerParallel("E:\\迅雷下载\\SMULTQ02-3_chr4.bam", "E:\\迅雷下载\\SMULTQ02-4_chr4.bam","E:\\迅雷下载\\dm6.fa" );
//        ps.process();
//        ps.print("out.txt");
//        System.out.println("ABCD".charAt(1));
//        File bamfile1 = new File("E:\\迅雷下载\\SMULTQ02-3.clean.fq.gz.Aligned.sortedByCoord.out.bam");
//        String genomefile = "E:\\迅雷下载\\dm6.fa";
        SamReader srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of("E:\\迅雷下载\\SMULTQ02-3_chr4.bam").
                index(new File(new File("E:\\迅雷下载\\SMULTQ02-3_chr4.bam").getAbsolutePath() + ".bai")
        ));
        SAMRecordIterator tempit1=srt.queryOverlapping("chr4", 0, 0);
        for (Iterator iterator = tempit1; iterator.hasNext();) {
            SAMRecord next = (SAMRecord) iterator.next();
//            System.out.print(next.getCigar()+"\t");
                if (next.getReadUnmappedFlag())continue;
                if (next.getDuplicateReadFlag()) continue;
                if (next.getNotPrimaryAlignmentFlag()) continue;
                if (next.getReadFailsVendorQualityCheckFlag()) continue;
                if(next.getMappingQuality()<=2) continue;
            List<CigarElement> ce=next.getCigar().getCigarElements(); 
            
            int le1=0;
            for (int i = 0; i < ce.size(); i++) {
                int le=ce.get(i).getLength();
            String op=ce.get(i).getOperator().toString();
                if (op.equals("S")) {
                    le1 += le;
                }
            }
            if(le1!=0){
//                System.out.println((double)le1/next.getReadLength()+"\t"+(next.getReadLength()-le1) );
                System.out.println((next.getReadLength()-le1) );
            }
            

//
        }
////
//
//        IndexedFastaSequenceFile genomeFile = new IndexedFastaSequenceFile(new File(genomefile));
//        System.out.println(genomeFile.getSequence("chr2L").getBaseString().charAt(1));
////        System.out.println(genomeFile.getSequenceDictionary().getSequence("chr2L"));
//SAMRecordIterator tempit1=srt.iterator();
//        
//        for (Iterator iterator = tempit1; iterator.hasNext();) {
//            SAMRecord next = (SAMRecord) iterator.next();
//            System.out.println(next.getReferencePositionAtReadPosition(1));
////
//        }
    }
}
