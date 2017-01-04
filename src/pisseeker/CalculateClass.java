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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Qi Zhao
 * @since 2017-1-3
 * @coding time 15:23:26
 * @author Qi Zhao
 */
public class CalculateClass {
    private double lamda;// lamda parameter for poisson distribution
    private ArrayList<PSIout> psilist=new ArrayList<PSIout>();

    public CalculateClass(HashSet<PSIout> psiset) {
        double libtreat1=0;
        double libtreat2=0;
        for (Iterator it = psiset.iterator(); it.hasNext();) {
            PSIout psi = (PSIout) it.next();
//            libtreat1+=psi.get
            psilist.add(psi);
        }
        
    }
    
    
    //get cigar reads with N covered posision N
   
   
    
    public static void main(String[] args) {
        SAMRecord sr =new SAMRecord(new SAMFileHeader());
        sr.setAlignmentStart(100);
        sr.setCigarString("23M100N20M");
//        System.out.println(CalculateClass.getMstart(223, sr));
    }
}
