/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

import java.io.FileNotFoundException;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
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
            System.out.println(ToolsforCMD.print_ansi_PURPLE("A java-based software for identifying pseudouridylic acid site in RNA from NGS data, version " + version) + "\r\n");

            System.out.println("Please input args\n Type "
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar ")
                    + ToolsforCMD.print_ansi_CYAN("-h")
                    + " for help\r\n");
        } else if (args[0].endsWith("-h")) {
            
            System.out.println(ToolsforCMD.print_ansi_YELLOW("Calling sites with input library :   \r\n\t\t")
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar -run")
                    + ToolsforCMD.print_ansi_CYAN(" <TreatBam(sorted)> <ControlBan(sorted)> <outputfile> "));
            System.out.println(ToolsforCMD.print_ansi_YELLOW("Calling sites with input library(Parallel Mode) :   \r\n\t\t")
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar -runMthread")
                    + ToolsforCMD.print_ansi_CYAN(" <TreatBam(sorted)> <ControlBan(sorted)> <outputfile> <threadnum> "));
        }else if (args[0].endsWith("-run")) {
            new PSIseeker(args[1],args[2],args[3]);
        }else if (args[0].endsWith("-runMthread")) {
            new PSIseekerParallel(args[1],args[2],args[3]);
        }
          long end = System.currentTimeMillis();
        System.out.println("Total running time is " + ToolsforCMD.print_ansi_YELLOW((end - start)  + "s"));
    }
}
