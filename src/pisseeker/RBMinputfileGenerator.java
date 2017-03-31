/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

import htsjdk.samtools.util.IntervalTree;
import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * @author zhaoqi
 */
public class RBMinputfileGenerator {
    
    public static String RBMinputfileGeneratorSingleChrome(IntervalTree resultTree){
        String str="";
//        Iterator it=resultTree.iterator();
        ArrayList<Double> plist=new ArrayList<Double>();
        for (Iterator it = resultTree.iterator(); it.hasNext();) {
            PSIout next = (PSIout) it.next();
            plist.add(next.getPoissonP()); 
        }
        
        return str;
    }
}
