/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pisseeker;

import java.io.FileNotFoundException;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import pisseeker.pub.FunctionClass;
import pisseeker.pub.ToolsforCMD;

/**
 *
 * @author Qi Zhao
 */
public class Console {

    public static void main(String[] args) throws IOException, FileNotFoundException, ParserConfigurationException {
        long start = System.currentTimeMillis();
        String version = "0.0.1";
//        System.out.println();

        if (args.length == 0) {

            System.out.println(ToolsforCMD.print_ansi_PURPLE(ToolsforCMD.getDAtoolstr()));
            System.out.println(ToolsforCMD.print_ansi_PURPLE("PSIseeker is a java-based command tool for identifying pseudouridylic acid site in RNA from NGS data, "
                    + "PSIseeker applied a simple statistic method, Fisher' Exact Test, to inference the candicated site that occuring pseudouridylic event. "
                    + "It requires sorted bamfiles from both treatment libarary and control library as input to perform comprehensive analysis. "
                    + "The current version is " + version) + "\r\n");

            System.out.println("Please input args\n Type "
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar ")
                    + ToolsforCMD.print_ansi_CYAN("-h")
                    + " for help\r\n");
        } else if (args[0].endsWith("-h")) {

            System.out.println(ToolsforCMD.print_ansi_YELLOW("Calling sites with input library :   \r\n\t\t")
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar -run")
                    + ToolsforCMD.print_ansi_CYAN(" <TreatBam(sorted)> <ControlBan(sorted)> <outputfile> [optional parameters] "));
            System.out.println(ToolsforCMD.print_ansi_YELLOW("Calling sites with input library(Parallel Mode) :   \r\n\t\t")
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar -runMthread")
                    + ToolsforCMD.print_ansi_CYAN(" <TreatBam(sorted)> <ControlBan(sorted)> <outputfile> <genome> <threadnum> [optional parameters] "));

            System.out.println(ToolsforCMD.print_ansi_WHITE("Extra paramters for options"));
            System.out.println(ToolsforCMD.print_ansi_RED("\r\n\t\t-covT\t ")
                    + ToolsforCMD.print_ansi_YELLOW("User defined minimum reads threshold for supporting truncated site in treatment and qualified reads in control,DEFAULT 2\r\n"));

        } else if (args[0].endsWith("-run")) {
            PSIseeker psiseeker = new PSIseeker(args[1], args[2]);
            if (FunctionClass.getArgsParameter(args, "-covT") != null) {
                psiseeker.filternumber = Integer.parseInt(FunctionClass.getArgsParameter(args, "-covT"));
                System.out.println(psiseeker.getFilternumber() + "minimun reads threshold is applied for analysis");
            }
            psiseeker.process();
            psiseeker.print(args[3]);
        } else if (args[0].endsWith("-runMthread")) {
            PSIseekerParallel psiseeker = new PSIseekerParallel(args[1], args[2],args[3]);
            if (FunctionClass.getArgsParameter(args, "-covT") != null) {
                psiseeker.filternumber = Integer.parseInt(FunctionClass.getArgsParameter(args, "-covT"));
                System.out.println(psiseeker.getFilternumber() + "minimun reads threshold is applied for analysis");
            }
            psiseeker.process();
            psiseeker.print(args[4]);
        }
        
        
        //get run time
        long end = System.currentTimeMillis();

        System.out.println("Total running time is " + ToolsforCMD.print_ansi_YELLOW((end - start) + "s"));

    }
}
