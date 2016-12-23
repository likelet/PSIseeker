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
        long start = System.nanoTime();
        String version = "0.0.1";
//        System.out.println();

        if (args.length == 0) {

            System.out.println(ToolsforCMD.print_ansi_PURPLE(ToolsforCMD.getDAtoolstr()));
            System.out.println(ToolsforCMD.print_ansi_PURPLE("Java-based Data Analysis tool for biological data process, version " + version) + "\r\n");

            System.out.println("Please input args\n Type "
                    + ToolsforCMD.print_ansi_GREEN("java -jar PSIseeker.jar ")
                    + ToolsforCMD.print_ansi_CYAN("-h")
                    + " for help\r\n");
        } else if (args[0].endsWith("-h")) {
        }
    }
}
