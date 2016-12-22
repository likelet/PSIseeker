/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
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
    private HashMap<String,ArrayList<PSIout>> pomap=new HashMap();
    /**
     * @param args the command line arguments
     */
    //tbam and cbam are both sorted bamfile 
    public PISseeker(String  tbam, String cbam){
        
        File bamfile1 = new File(tbam);
        File bamfile2 = new File(cbam);
        SAMFileReader srt = new SAMFileReader(bamfile1, new File(bamfile1.getAbsolutePath() + ".bai"));
        SAMFileReader src = new SAMFileReader(bamfile2, new File(bamfile2.getAbsolutePath() + ".bai"));
    }
   
    public void process(){
        LinkedList<SAMRecord> ls1=new LinkedList<SAMRecord>();
         LinkedList<SAMRecord> ls2=new LinkedList<SAMRecord>();
        HashSet<String> hashposition=new HashSet<String>();
         HashSet<String> chrsome=new HashSet<String>();
         
        Iterator iter1 = srt.iterator();
        Iterator iter2 = src.iterator();
        
        while (iter1.hasNext()) {
            SAMRecord sitem = (SAMRecord) iter1.next();
            String chrom = sitem.getReferenceName();
            int start = sitem.getAlignmentStart();
            int end = sitem.getAlignmentEnd();
            if (hashposition.add(chrom + end)) {
                PSIout po = new PSIout(chrom, end, "");
                if (chrsome.add(chrom)) {
                    ArrayList<PSIout> psilist = new ArrayList<PSIout>();
                    psilist.add(po);
                    pomap.put(chrom, psilist);
                } else {
                    pomap.get(chrom).add(po);
                }
            }else{
                int index=pomap.get(chrom).size()-1;//array last one
                pomap.get(chrom).get(index).add1supporCountInTreat();
            }       
        }
    }
}
