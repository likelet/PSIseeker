/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import pisseeker.MultipleCorrection.FDR;

/**
 *
 * @author Administrator
 */
public class PISseeker {
//    private String gaifile="";

    private final SamReader srt;
    private final SamReader src;
    private int readsRegion = 500;

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
    public PISseeker(String tbam, String cbam, String out) throws IOException {

        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        srt = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile1).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
        src = SamReaderFactory.makeDefault().open(
                SamInputResource.of(bamfile2).
                index(new File(bamfile1.getAbsolutePath() + ".bai"))
        );
        this.process();
        this.print(out);

    }

    public void process() throws IOException {

//    
        // in one specific chrome
        ArrayList<SAMRecord> samlistTreat = new ArrayList<SAMRecord>();
//        String chrome="chr2L";
//        Iterator iter = srt.query(chrome,0,0,true);
//        System.out.println("Get record from record "+chrome );
//        
        HashSet<Integer> hashposition = new HashSet<Integer>();

        //first cycle initialized chromesome and positions
        SAMRecordIterator iter = srt.iterator();

        System.out.println("Start first reading cycle for position infomation");

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
            if (chrlist.add(chrome)) {
                HashSet<PSIout> pomap = new HashSet<PSIout>();
                if (strand) {
                    PSIout po = new PSIout(chrome, end, strand);
                    pomap.add(po);
                    negativeResultMap.put(chrome, pomap);
                } else {
                    PSIout po = new PSIout(chrome, start, strand);
                    pomap.add(po);
                    positiveResultMap.put(chrome, pomap);
                }

            } else if (strand) {
                PSIout po = new PSIout(chrome, end, strand);
                negativeResultMap.get(chrome).add(po);
            } else {
                PSIout po = new PSIout(chrome, start, strand);
                positiveResultMap.get(chrome).add(po);
            }
        }
        iter.close();

        System.out.println("Start Counting reads");
        for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            HashSet<PSIout> negpomap = negativeResultMap.get(chr);
            for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                SAMRecordIterator tempit1 = srt.queryContained(chr, pso.getPosition() - readsRegion, pso.getPosition() + readsRegion);
                SAMRecordIterator tempit2 = src.queryContained(chr, pso.getPosition() - readsRegion, pso.getPosition() + readsRegion);
                CountReadsTreat(pso,tempit1);
                CountReadsControl(pso, tempit2);
            }
            HashSet<PSIout> pospomap = positiveResultMap.get(chr);
            for (Iterator poit = pospomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                SAMRecordIterator tempit1 = srt.queryContained(chr, pso.getPosition() - readsRegion, pso.getPosition() + readsRegion);
                SAMRecordIterator tempit2 = src.queryContained(chr, pso.getPosition() - readsRegion, pso.getPosition() + readsRegion);
                CountReadsTreat(pso,tempit1);
                CountReadsControl(pso, tempit2);
            }
        }

    }

    public void CountReadsTreat(PSIout po, SAMRecordIterator readsIt) {
        for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
            SAMRecord sr = iterator.next();
            boolean strand = sr.getReadNegativeStrandFlag();
            if (!strand && po.isStrandB()) {
                continue;
            }
            if (iscover(po.getPosition(), sr).equals("onsite")) {
                po.add1supporCountInTreat();
            }
            if (iscover(po.getPosition(), sr).equals("cover")) {
                po.add1totalCountInTreat();
            }
            readsIt.close();

        }
        
    }
    
    public void CountReadsControl(PSIout po, Iterator readsIt) {
        for (Iterator<SAMRecord> iterator = readsIt; iterator.hasNext();) {
            SAMRecord sr = iterator.next();
            boolean strand = sr.getReadNegativeStrandFlag();
            if (!strand && po.isStrandB()) {
                continue;
            }
            if (iscover(po.getPosition(), sr).equals("onsite")) {
                po.add1supporCountControl();
            }
            if (iscover(po.getPosition(), sr).equals("cover")) {
                po.add1totalCountControl();
            }

        }
    }
    
    
    
    public void print(String fileout) throws IOException {
        System.out.println("Writing out..");
        FileWriter fw = new FileWriter(fileout);
        ArrayList<PSIout> templist = new ArrayList<PSIout>();
        ArrayList<Double> plist = new ArrayList<Double>();
        fw.append("chr\tposition\tbase\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\n");

          for (Iterator chrit = chrlist.iterator(); chrit.hasNext();) {
            String chr = (String) chrit.next();
            HashSet<PSIout> negpomap = negativeResultMap.get(chr);
            for (Iterator poit = negpomap.iterator(); poit.hasNext();) {
                PSIout pso = (PSIout) poit.next();
                if (pso.getSupporCountInTreat() < this.filternumber) {
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
            fw.append(psi.toString() + "\t" + "\r\n");
        }

        fw.flush();
        fw.close();
    }

    public String iscover(int pos, SAMRecord sr) {

        if (sr.getAlignmentEnd() == pos) {
            return "onsite";
        } else if (sr.getAlignmentEnd() > pos && sr.getAlignmentStart() <= pos) {
            return "cover";
        } else if (sr.getAlignmentEnd() < pos) {
            return "down";
        } else {
            return "up";
        }
    }
}
