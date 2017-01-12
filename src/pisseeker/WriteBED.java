/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pisseeker;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.RuntimeIOException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 *
 * @author Administrator
 * @since 2017-1-12
 * @coding time 14:33:32
 * @author Qi Zhao
 */
public class WriteBED {

    public static void WriteBED(HashMap<String, IntervalTree> intervalMap, File bedout) {

        IOUtil.assertFileIsWritable(bedout);
        try {
            final BufferedWriter out = IOUtil.openFileForBufferedWriting(bedout);
            int score = 500;
            for (Iterator it = intervalMap.keySet().iterator(); it.hasNext();) {
                String chr = (String) it.next();
                if(!chr.equalsIgnoreCase("XI")) continue;
                IntervalTree temptree = intervalMap.get(chr);
                Iterator tempit1 = temptree.iterator();
                for (Iterator it2 = tempit1; it2.hasNext();) {
                    IntervalTree.Node nd = (IntervalTree.Node) it2.next();
                    SAMRecord sr = (SAMRecord) nd.getValue();
                    final String strand = sr.getReadNegativeStrandFlag() ? "-" : "+";
                    final List<?> fields = CollectionUtil.makeList(sr.getReferenceName(), nd.getStart(), nd.getEnd(), sr.getReadName(), score, strand);
                    out.append(fields.stream().map(String::valueOf).collect(Collectors.joining("\t")));
                    out.newLine();

                }
            }
            out.close();
        } catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }

    }
}
