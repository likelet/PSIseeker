/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

/**
 *
 * @author Administrator
 * @since 2017-1-5
 * @coding time 16:11:08
 * @author Qi Zhao
 */
public class basicTerm {
    private String chr;
    private int position;
    private int supporCountControl=0;
    private int totalCountControl=0;
    private double backRatio;

    public basicTerm(String chr, int pos) {
        this.chr = chr;
        this.position = pos;
    }
    
    public void add1supporCountControl(){
        this.supporCountControl++;
    }
    public void add1totalCountControl(){
        this.totalCountControl++;
    }

    
    
    
    public void getE(){
        if(supporCountControl==0)
         backRatio=1/(double)totalCountControl;   
    }

    public String getChr() {
        return chr;
    }

    public int getPosition() {
        return position;
    }

    public int getSupporCountControl() {
        return supporCountControl;
    }

    public int getTotalCountControl() {
        return totalCountControl;
    }

    public double getbackRatio() {
        return backRatio;
    }
    
    
   
    
}
