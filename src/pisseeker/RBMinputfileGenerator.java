/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

import htsjdk.samtools.util.IntervalTree;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;

/**
 *
 * @author zhaoqi
 */
public class RBMinputfileGenerator {
    
    
    
    public  IntervalTree hashset2resultTree(HashSet<PSIout> psiHashset){
        IntervalTree potree=new IntervalTree(); 
        for (Iterator it = psiHashset.iterator(); it.hasNext();) {
            PSIout psi = (PSIout) it.next();
             potree.put(psi.getPosition(), psi.getPosition(), psi);
        }
        return potree;
    }
    
    
    // flanking position number should be 1(3),2(5),3(7),4(9)...10
    public static ArrayList<String> RBMinputfileGeneratorSingleChrome(IntervalTree resultTree, int range) {
        ArrayList<String> strlist = new ArrayList<String>();
//        Iterator it=resultTree.iterator();
        LinkedList<Double> plist = new LinkedList<Double>();
        for (Iterator it = resultTree.iterator(); it.hasNext();) {
            PSIout next = (PSIout) it.next();
            plist.add(next.getPoissonP());
        }

        //get poisson P aroung the each position
        //avoid edge numbers
        for (int i = 0; i < range; i++) {
            plist.addFirst(0.0);
        }
        for (int i = 0; i < range; i++) {
            plist.addLast(0.0);
        }
        String str = "";
        for (int i = 0; i < plist.size(); i++) {
//            if(plist.get(i)==0) continue;
            if (plist.get(i + range) == 0) {
                continue;
            }
            str = "";
            for (int j = i; j < i + 2 * range; j++) {
                str = str + plist.get(j) + "\t";
            }
            if (str != "") {
                strlist.add(str.substring(0, str.length() - 1));// remove last tab 
            }
        }

        return strlist;
    }
    
  
}
