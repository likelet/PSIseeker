/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.SortedMap;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 *
 * @author Administrator
 */
public class PISseeker {
//    private String gaifile="";
    private  SAMFileReader srt;
    private  SAMFileReader src;
    
    //storage chrome and position information
    private HashMap<String,LinkedHashMap <Integer,PSIout>> PositiveResultMap=new HashMap<String,LinkedHashMap <Integer,PSIout>>();// result in positve strand 
    private HashMap<String,LinkedHashMap <Integer,PSIout>> negativeResultMap=new HashMap<String,LinkedHashMap <Integer,PSIout>>();// result in negative strand;
//    private HashMap<Integer,PSIout> pomap=new HashMap();
    public int filternumber=2;//at least 2 read support
    public HashSet<String> chrlist=new HashSet<String>();
    /**
     * @param args the command line arguments
     */
    //tbam and cbam are both sorted bamfile 
    public PISseeker(String  tbam, String cbam,String out) throws IOException{
        
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
         srt = new SAMFileReader(bamfile1, new File(bamfile1.getAbsolutePath() + ".bai"));
         src = new SAMFileReader(bamfile2, new File(bamfile2.getAbsolutePath() + ".bai"));
        this.process();
        this.print(out);
        
    }
   
    public void process(){
        
//    
      
      // in one specific chrome
        ArrayList<SAMRecord> samlistTreat=new ArrayList<SAMRecord> ();
//        String chrome="chr2L";
//        Iterator iter = srt.query(chrome,0,0,true);
//        System.out.println("Get record from record "+chrome );
//        
        HashSet<Integer> hashposition = new HashSet<Integer>();
        
        //first cycle initialized chromesome and positions
        Iterator iter = srt.iterator();
        System.out.println("Start first reading cycle for position infomation");
        ArrayList<Integer> positionlist=new ArrayList<Integer>();
        while (iter.hasNext()) {
            SAMRecord sitem = (SAMRecord) iter.next();
            if(sitem.getReadUnmappedFlag()) continue;
            if(sitem.getDuplicateReadFlag()) continue;
            if(sitem.getNotPrimaryAlignmentFlag()) continue;
            if(sitem.getReadFailsVendorQualityCheckFlag()) continue;
            int end = sitem.getAlignmentEnd();
            int start= sitem.getAlignmentStart();
            boolean strand=sitem.getReadNegativeStrandFlag();//strand of the query (false for forward; true for reverse strand).
           
            String  chrome=sitem.getReferenceName();
            if(chrlist.add(chrome)){
                LinkedHashMap <Integer,PSIout> pomap =new LinkedHashMap <Integer,PSIout>() ;
                if(strand){
                    PSIout po = new PSIout(chrome,end,strand);
                    pomap.put(end, po);
                    negativeResultMap.put(chrome, pomap);
                }else{
                    PSIout po = new PSIout(chrome,start,strand);
                    pomap.put(start, po);
                    negativeResultMap.put(chrome, pomap);
                }
                
            }else{
                if(resultMap.get(chrome).keySet().contains(key)) continue;
                PSIout po = new PSIout(chrome,end,strand);
                resultMap.get(chrome).put(key, po);
            }
            samlistTreat.add(sitem);
        }
//        srt.close();
        System.out.println("sorting");
        Collections.sort(positionlist);
        //second cycles
        System.out.println("Start second traversing cycle counting tags from treatment");
        Iterator iter1 = samlistTreat.iterator();
        int startm = 0;
        int endm = positionlist.size();
        while (iter1.hasNext()) {
            SAMRecord sitem = (SAMRecord) iter1.next();
            endm = positionlist.indexOf(sitem.getAlignmentEnd());
            for (int i = startm; i <=endm; i++) {
                int end = positionlist.get(i);
                if (this.iscover(end, sitem).equals("cover")) {
                    pomap.get(end).add1totalCountInTreat();
                } else if (this.iscover(end, sitem).equals("onsite")) {
                    pomap.get(end).add1totalCountInTreat();
                    pomap.get(end).add1supporCountInTreat();
                } else if (this.iscover(end, sitem).equals("up")) {
                    startm++;
                }
            }
        }
        Iterator iter2 = src.query(chrome,0,0,true);
        System.out.println("Start second traversing cycle counting tags from input library");
        startm = 0;
        endm = positionlist.size();
        while (iter2.hasNext()) {
            SAMRecord sitem = (SAMRecord) iter2.next();
            if(sitem.getReadUnmappedFlag()) continue;
            if(sitem.getDuplicateReadFlag()) continue;
            if(sitem.getNotPrimaryAlignmentFlag()) continue;
            if(sitem.getReadFailsVendorQualityCheckFlag()) continue;
            endm=positionlist.indexOf(sitem.getAlignmentEnd());
            for (int i = startm; i <=endm; i++) {
                int end=positionlist.get(i);
                if(this.iscover(end, sitem).equals("cover")){
                        pomap.get(end).add1totalCountControl();
                    }else if(this.iscover(end, sitem).equals("onsite")){
                        pomap.get(end).add1totalCountControl();
                        pomap.get(end).add1supporCountControl();
                    }else if(this.iscover(end, sitem).equals("up")){
                        startm++;
                    }
            } 
        }

    }
    
    public void print(String fileout) throws IOException {
        System.out.println("Writing out..");
        FileWriter fw = new FileWriter(fileout);
        ArrayList<PSIout> templist=new  ArrayList<PSIout>();
        ArrayList<Double> plist=new ArrayList<Double>();
        fw.append("chr\tposition\tbase\tsupporCountInTreat\ttotalCountInTreat\tsupporCountControl\ttotalCountControl\tPvalue\tadjustP\n");

        for (Iterator it = pomap.keySet().iterator(); it.hasNext();) {
            int pos = (int) it.next();
            PSIout psi = pomap.get(pos);

//            psi.fishertest();
            if (psi.getSupporCountInTreat() < this.filternumber) {
                continue;
            }
            psi.fishertest();
            plist.add(psi.getPvalue());
            templist.add(psi);

        }
        //FDR calculation && write out
        ArrayList<Double> fdrlist=new FDR(plist).getFdrDoublelist();
        for (int i = 0; i < templist.size(); i++) {
            PSIout psi=templist.get(i);
            psi.setAdjustP(fdrlist.get(i));
            fw.append(psi.toString() +"\t"+ "\r\n");
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
