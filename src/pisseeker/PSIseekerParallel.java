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
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
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
    public  int Thread=4;
    public String tbam;
    public String cbam;
    public IndexedFastaSequenceFile Indexgenomefile;
    //storage chrome and position information
    private HashMap<String, HashSet<PSIout>> positiveResultMap = new HashMap<String, HashSet<PSIout>>();// result in positve strand 
    private HashMap<String, HashSet<PSIout>> negativeResultMap = new HashMap<String, HashSet<PSIout>>();// result in negative strand;
//    private HashMap<Integer,PSIout> pomap=new HashMap();
    public int filternumber = 2;//at least 2 read support
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
        this.Indexgenomefile=new IndexedFastaSequenceFile (new File(genomefile));;
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
        this.Indexgenomefile=new IndexedFastaSequenceFile (new File(genomefile));;
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
        SAMRecordIterator iter = srt.iterator();
        //get chrome information
        SAMFileHeader fileHeader = srt.getFileHeader();
        SAMSequenceDictionary seqDic = fileHeader.getSequenceDictionary();
        List<SAMSequenceRecord> seqRecList = seqDic.getSequences();
        for (SAMSequenceRecord seqRec : seqRecList) {
            //store chr information
            this.chrlist.add(seqRec.getSequenceName());
        }
        System.out.println("Initializing candidate site infomation");
        while (iter.hasNext()) {
            SAMRecord sitem = (SAMRecord) iter.next();
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
            boolean strand = sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand).

            String chrome = sitem.getReferenceName();
            if (strand) {
                if (!negativeResultMap.containsKey(chrome)) {
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setExbase(Indexgenomefile.getSequence(chrome).getBaseString().charAt(end+1));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    negativeResultMap.put(chrome, pomap);
                } else {
                    PSIout po = new PSIout(chrome, end, strand);
                    po.setBase(new DNAsequenceProcess().getReverseComplimentary(sitem.getReadString()).charAt(0));
                    po.setReadsString(sitem.getReadString());
                    po.setExbase(Indexgenomefile.getSequence(chrome).getBaseString().charAt(end+1));
                    negativeResultMap.get(chrome).add(po);
                }
            } else if (!strand) {
                if (!positiveResultMap.containsKey(chrome)) {
                    HashSet<PSIout> pomap = new HashSet<PSIout>();
                    PSIout po = new PSIout(chrome, start, strand);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setExbase(Indexgenomefile.getSequence(chrome).getBaseString().charAt(start-1));
                    po.setReadsString(sitem.getReadString());
                    pomap.add(po);
                    positiveResultMap.put(chrome, pomap);
                } else {
                    PSIout po = new PSIout(chrome, start, strand);
//                    System.out.println(chrome);
                    po.setBase(sitem.getReadString().charAt(0));
                    po.setExbase(Indexgenomefile.getSequence(chrome).getBaseString().charAt(start-1));
                    po.setReadsString(sitem.getReadString());
                    positiveResultMap.get(chrome).add(po);
                }
            }

        }
        iter.close();
    }
    
    public void process() throws IOException {
        

        System.out.println("Start Counting reads in parrallel mode ");
        try {
            // run paralla
            ExecutorService pool = Executors.newFixedThreadPool(this.Thread);//creat a new thread pool
            
            runPSIseekerThread runPSIseekerthread = null;
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                runPSIseekerthread=new runPSIseekerThread(chr,negativeResultMap.get(chr),tbam,cbam);
                pool.submit(runPSIseekerthread);
            }
            for (Iterator<String> iterator = chrlist.iterator(); iterator.hasNext();) {
                String chr = iterator.next();
                runPSIseekerthread=new runPSIseekerThread(chr,positiveResultMap.get(chr),tbam,cbam);
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

        public void CountReadsTreat(PSIout po, SAMRecordIterator readsIt) {
            for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
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
            readsIt.close();

        }

        public void CountReadsControl(PSIout po, SAMRecordIterator readsIt) {
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
        fw.append("chr\tposition\tbase\tstrand\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\n");

        for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();

            HashSet<PSIout> negpomap = negativeResultMap.get(chr);
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
                if (pso.getSupporCountInTreat() / pso.getTotalCountInTreat() < pso.getSupporCountControl() / pso.getTotalCountControl()) {
                    continue;
                }
                pso.fishertest();
                plist.add(pso.getPvalue());
                templist.add(pso);

            }
            HashSet<PSIout> pospomap = positiveResultMap.get(chr);
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
                if (pso.getSupporCountInTreat() / pso.getTotalCountInTreat() < pso.getSupporCountControl() / pso.getTotalCountControl()) {
                    continue;
                }
                pso.fishertest();
                plist.add(pso.getPvalue());
                templist.add(pso);

            }
        }
        //FDR calculation && write out
        ArrayList<Double> fdrlist = new FDR(plist).getFdrDoublelist();
        for (int i = 0; i < templist.size(); i++) {
            PSIout psi = templist.get(i);
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
//        new PSIseekerParallel("test/SMULTQ02-1.clean.fq.gz.Aligned.sortedByCoord.out.bam", "test/SMULTQ02-2.clean.fq.gz.Aligned.sortedByCoord.out.bam", "test/result.txt","2");

        File bamfile1 = new File("E:\\迅雷下载\\SMULTQ02-4.clean.fq.gz.Aligned.sortedByCoord.out.bam");
        String genomefile="E:\\迅雷下载\\dm6.fa";
          SamReader srt = SamReaderFactory.makeDefault().referenceSequence(new File(genomefile)).open(
                SamInputResource.of(bamfile1).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
//
  
    IndexedFastaSequenceFile  genomeFile=new IndexedFastaSequenceFile (new File(genomefile));
        System.out.println(genomeFile.getSequence("chr2L").getBaseString().charAt(0));
//        System.out.println(genomeFile.getSequenceDictionary().getSequence("chr2L"));
//SAMRecordIterator tempit1=srt.iterator();
//        
//        for (Iterator iterator = tempit1; iterator.hasNext();) {
//            SAMRecord next = (SAMRecord) iterator.next();
//            System.out.println(next.getReferencePositionAtReadPosition(1));
////
//        }
    }
}
