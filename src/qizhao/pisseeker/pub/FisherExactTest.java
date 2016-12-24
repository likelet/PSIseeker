package pisseeker.pub;


import java.math.BigDecimal;
//import java.math.BigInteger;
import java.math.MathContext;

/*
 * 适用于样本数小于5的情况下2X2表的精确检验
 * 
 * 
 */

/**
 *
 * @author ZHAO Qi
 * @version jdk 1.7.0
 */
public class FisherExactTest {
    
     private MathContext mc=new MathContext(7);
     private double[] log;
    
     
      public double getTwoTailP(int a,int b,int c,int d){
        double p1=getLeftTailP(a,b,c,d);
        double p2=getRightTailP(a,b,c,d);
        
        //System.out.println(p1+"\t"+p2);
        
        return 2*Math.min(Math.min(p1,p2),0.5);
    }
     
//    public double getTwoTailP(int a,int b,int c,int d){
//        //左尾概率
//        BigDecimal p1 = BigDecimal.valueOf(0);
//        //右尾概率
//        BigDecimal p2 = BigDecimal.valueOf(0);
//        
//        int af,bf,cf,df;
//        //重新排列2X2表
//        int temp2=Math.min(Math.min(a+b, c+d), Math.min(a+c, b+d));
//        
//        int temp = Math.min(Math.min(a, b), Math.min(c, d));
//        if(temp==a){
//            af=a;
//            bf=b;
//            cf=c;
//            df=d;
//        }else if(temp==b){
//            af=b;
//            bf=a;
//            cf=d;
//            df=c;
//        }else if(temp==c){
//            af=c;
//            bf=d;
//            cf=a;
//            df=b;
//        }else{
//            af=d;
//            bf=c;
//            cf=b;
//            df=a;
//        }
//        //System.out.println(af+" "+bf+" "+cf+" "+df);
//        //BigDecimal Pv=getPSingle(af,af,bf,cf,df);
//        for (int i = 0; i <= af; i++) {
//            //System.out.println(getPSingle(i,af,bf,cf,df)+"\t ok");
//                p1=p1.add(getPSingle(i,af,bf,cf,df),mc);
//            
//        }
//        
//        for (int i = af; i <=temp2; i++) {
//            p2=p2.add(getPSingle(i,af,bf,cf,df),mc);
//        }
//        return 2*Math.min(Math.min(p1.doubleValue(),p2.doubleValue()),0.5);
//    }
    //计算左尾概率
    public double getLeftTailP(int a,int b,int c,int d){
        // p1 = BigDecimal.valueOf(0);
        double p1=0;
        int temp2=Math.min(Math.min(a+b, c+d), Math.min(a+c, b+d));
        //System.out.println("temp2"+temp2);
        int arr[]=this.reArrange(a, b, c, d); 
        
        for (int i = 0; i <=arr[0]; i++) {
            p1=p1+getPSingle(i,arr[0],arr[1],arr[2],arr[3]);
        }
        //System.out.println(p1.);
        return p1;
    }
    
    //计算右尾
    public double getRightTailP(int a,int b,int c,int d){
        //BigDecimal p1 = BigDecimal.valueOf(0);
        double p1=0;
        int temp2=Math.min(Math.min(a+b, c+d), Math.min(a+c, b+d));
        //System.out.println("temp2== "+temp2);
        int arr[]=this.reArrange(a, b, c, d); 
        //System.out.println(arr[0]+arr[1]+arr[2]+arr[3]);
        //System.out.println(arr[0]);
        for (int i = arr[0]; i <=temp2; i++) {
            //System.out.println(getPSingle(i,arr[0],arr[1],arr[2],arr[3]));
            p1+=getPSingle(i,arr[0],arr[1],arr[2],arr[3]);
            //System.out.println(p1.doubleValue());
            //System.out.println(p1);
        }
        
        return p1;
    }
    
    //返回构成每个表的精确概率（符合超几何分布）
    public double getPSingle(int x,int af,int bf,int cf, int df){
       
        BigDecimal P;
        int n1,n2,m1,m2=0;
        n1=af+bf;
        n2=cf+df;
        m1=af+cf;
        m2=bf+df;
        int n=af+bf+cf+df;
        
        
        this.log = new double[n + 1];
         log[0] = 0.0;
        //log[i] = logi + log(i - 1) + log(i - 2) + ... + log1
        for (int i = 1; i <= n; i++) {
            log[i] = log[i - 1] + Math.log(i);
        }
        
        double pr = (log[n1] + log[n2] + log[m1] + log[m2]) - (log[n1-x]+log[m1-x]+log[m2-n1+x]+ log[af + bf + cf + df]+log[x]);
        return Math.exp(pr);
        //System.out.println(x+" "+(n1-x)+" "+(m1-x)+" "+(m2-n1+x)+" "+n);
        //System.out.println(mutiT(n)+"\t"+mutiT(x)+"\t"+mutiT(n1-x)+"\t"+mutiT(m1-x)+"\t"+mutiT(m2-n1+x));
//        BigDecimal x1=mutiT(n1).divide(mutiT(n),mc);
//        BigDecimal x2=mutiT(n2).divide(mutiT(n1-x),mc);
//        BigDecimal x3=mutiT(m1).divide(mutiT(m1-x),mc);
//        BigDecimal x4=mutiT(m2).divide(mutiT(m2-n1+x).multiply(mutiT(x)),mc);
//        P=x1.multiply(x2.multiply(x3.multiply(x4)));
       //System.out.println(P.doubleValue());
        //int af,bf,cf,df=0;
        
        //return P;
    }
    
    //计算阶乘
//    public BigDecimal mutiT(int i){
//         BigDecimal sum=BigDecimal.valueOf(1);
//        if(i==0||i==1){
//            return BigDecimal.valueOf(1);
//        }else{
//            for (int j = 1; j < i+1; j++) {
//                sum=sum.multiply(BigDecimal.valueOf(j));
//            }
//            return sum;
//        }
//        
//    }
//    
    //重排2联表
    public int[]reArrange(int a,int b,int c,int d){
        int[] arr=new int[4];
        
        //int af,bf,cf,df;
     
        //int temp2=Math.min(Math.min(a+b, c+d), Math.min(a+c, b+d));
        int temp = Math.min(Math.min(a, b), Math.min(c, d));
        if(temp==a){
            arr[0]=a;
            arr[1]=b;
            arr[2]=c;
            arr[3]=d;
        }else if(temp==b){
            arr[0]=b;
            arr[1]=a;
            arr[2]=d;
            arr[3]=c;
        }else if(temp==c){
            arr[0]=c;
            arr[1]=d;
            arr[2]=a;
            arr[3]=b;
        }else{
            arr[0]=d;
            arr[1]=c;
            arr[2]=b;
            arr[3]=a;
        }
        //System.out.println(arr[);
        return arr;
    }
    
    public static void main(String arg[]) {
        FisherExactTest test=new FisherExactTest();
        //System.out.println(test.mutiT(5));
        //System.out.println(test.getPSingle(2, 2,23,5,30));
        System.out.println(test.getPSingle(0,2,23,5,30));
       System.out.println(test.getTwoTailP(2,23,5,30));
    }
}
