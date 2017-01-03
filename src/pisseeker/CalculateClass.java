/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

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
    
    
    
    
}
