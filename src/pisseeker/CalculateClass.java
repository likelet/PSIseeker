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
import htsjdk.samtools.util.IntervalTree;
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
import org.apache.commons.math3.distribution.PoissonDistribution;

/**
 *
 * @author Qi Zhao
 * @since 2017-1-3
 * @coding time 15:23:26
 * @author Qi Zhao
 */
public class CalculateClass {

    private HashMap<String, HashSet<PSIout>> positiveResultMap;// result in positve strand 
    private HashMap<String, HashSet<PSIout>> negativeResultMap;// result in negative strand;

    private HashMap<String, IntervalTree> positiveBackGroundMap;// result in positve strand 
    private HashMap<String, IntervalTree> negativeBackGroundMap;// result in negative strand;

    private HashMap<String, Double> positiveGloableLamda;
    private HashMap<String, Double> negativeGloableLamda;
    private int thread;
    private HashSet<String> chrlist;

    public CalculateClass(HashMap<String, HashSet<PSIout>> positiveResultMap, HashMap<String, HashSet<PSIout>> negativeResultMap, HashMap<String, IntervalTree> positiveBackGroundMap, HashMap<String, IntervalTree> negativeBackGroundMap, HashMap<String, Double> positiveGloableLamda, HashMap<String, Double> negativeGloableLamda, int thread, HashSet<String> chrlist) {
        this.positiveResultMap = positiveResultMap;
        this.negativeResultMap = negativeResultMap;
        this.positiveBackGroundMap = positiveBackGroundMap;
        this.negativeBackGroundMap = negativeBackGroundMap;
        this.positiveGloableLamda = positiveGloableLamda;
        this.negativeGloableLamda = negativeGloableLamda;
        this.thread = thread;
        this.chrlist = chrlist;
        this.process();
    }

    
    
    public void process(){
         try {
            ExecutorService pool = Executors.newFixedThreadPool(this.thread);
//            PSIseeker.getPoissonPthread runPSIseekerthread = null;
            for (Iterator iterator =chrlist .iterator(); iterator.hasNext();) {
                String chr = (String) iterator.next();
                if ((this.positiveResultMap.get(chr) == null) && (this.positiveBackGroundMap.get(chr) == null)) {
                    continue;
                }
                getPoissonPthread getP= new getPoissonPthread(chr, this.positiveResultMap.get(chr),  this.positiveBackGroundMap.get(chr),this.positiveGloableLamda.get(chr));
                pool.submit(getP);
            }
            for (Iterator iterator = this.chrlist.iterator(); iterator.hasNext();) {
                String chr = (String) iterator.next();
                if ((this.negativeResultMap.get(chr) == null) && (this.negativeBackGroundMap.get(chr) == null)) {
                    continue;
                }
               getPoissonPthread getP= new getPoissonPthread(chr, this.negativeResultMap.get(chr),  this.negativeBackGroundMap.get(chr),this.negativeGloableLamda.get(chr));
               pool.submit(getP);
            }
            pool.shutdown();
            pool.awaitTermination(9223372036854775807L, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            Logger.getLogger(PSIseeker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }


    
    
    class getPoissonPthread implements Runnable {
        private String chr;
        private HashSet<PSIout> subResultMap;
        private IntervalTree subIntervalMap;
        private double globalLam;

        public getPoissonPthread(String chr, HashSet<PSIout> subResultMap, IntervalTree subIntervalMap, double globalLam) {
            this.chr = chr;
            this.subResultMap = subResultMap;
            this.subIntervalMap = subIntervalMap;
            this.globalLam = globalLam;
        }
        
        

        public void run() {
            for (Iterator it = subResultMap.iterator(); it.hasNext();) {
                PSIout psi = (PSIout) it.next();
                double lamd5k=getAverageLamda(subIntervalMap.overlappers(psi.getPosition()-2500, psi.getPosition()+2500));
                double lamd1k=getAverageLamda(subIntervalMap.overlappers(psi.getPosition()-500, psi.getPosition()+500));
                double lamda=Math.max(lamd1k, lamd5k);
                lamda=Math.max(lamda, globalLam);
                PoissonDistribution pdis=new PoissonDistribution(lamda*50);
                double currentRatio=(double)psi.getSupporCountInTreat()/(double)psi.getTotalCountInTreat();
                psi.setPoissonP(1-pdis.cumulativeProbability((int)(currentRatio*50)));
                
            }
        }
        
        //Get average enrichment ratio of a background interval 
        public double getAverageLamda(Iterator iterm){
            double tamplam=0;
            double count=0;
            for (Iterator it = iterm; it.hasNext();) {
                IntervalTree.Node nd = (IntervalTree.Node) it.next();
                PSIout ps = (PSIout) nd.getValue();
                count++;
                tamplam+=ps.getBackRatio();
            }
            return tamplam/count;
            
        }

    }
    
    //get cigar reads with N covered posision N
   
    
 
}
