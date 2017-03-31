package pisseeker.MultipleCorrection;

import org.datavec.api.records.reader.RecordReader;
import org.datavec.api.records.reader.impl.csv.CSVRecordReader;
import org.datavec.api.split.FileSplit;
import org.deeplearning4j.datasets.datavec.RecordReaderDataSetIterator;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.layers.RBM;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Restricted Boltzmann Machine Pipeline
 *
 * Based on Hinton et al.'s work
 *
 * Great reference:
 * http://www.iro.umontreal.ca/~lisa/publications2/index.php/publications/show/239
 *
 * Created by redfish on 2017/3/30.
 */
public class RBMPipeline {
    /**
     * Initial RBM model and Pre-train the model
     *
     * @param nIn  the visible input numbers
     * @param nOut the hidden output numbers
     * @param trainFilePath the filepath of training file
     * @return  MultiLayerNetwork the trained model
    **/
    public MultiLayerNetwork modelInitialandTrain(Integer nIn, Integer nOut,String trainFilePath) throws IOException, InterruptedException {
        //initial the recordReader
        RecordReader recordReader = new CSVRecordReader(0,",");
        recordReader.initialize(new FileSplit(new File(trainFilePath)));
        DataSetIterator iter2 = new RecordReaderDataSetIterator(recordReader,5000);
        //initial the model
        MultiLayerConfiguration conf =  new NeuralNetConfiguration.Builder()
                .seed(123)
                .iterations(1)
                .list()
                .layer(0, new RBM.Builder()
                        .preTrainIterations(1)
                        .nIn(nIn)
                        .nOut(nOut)
                        .build()
                )
                .pretrain(true).backprop(false)
                .build();
        MultiLayerNetwork model = new MultiLayerNetwork(conf);
        model.init();
        //Pre-train the model
        model.fit(iter2);
        return model;
    }

    /**
     * Initial RBM model and Pre-train the model
     *
     * @param model  the Pre-trained model
     * @param evaluationFilepath the filepath for evaluation
     * @param resultPath the filepath of result
     **/
    public void getEvaluationFile( MultiLayerNetwork model,String evaluationFilepath,String resultPath) throws IOException, InterruptedException {
        FileOutputStream resultOut = new FileOutputStream(resultPath,true);
        StringBuffer temp_out = new StringBuffer();
        RecordReader recordReader2 = new CSVRecordReader(0,",");
        recordReader2.initialize(new FileSplit(new File(evaluationFilepath)));
        DataSetIterator test = new RecordReaderDataSetIterator(recordReader2,1);
        while(test.hasNext()){
            DataSet currentDataSet = test.next();
            INDArray testdata = currentDataSet.getFeatureMatrix();
            //model.output mean get the activations of input
            INDArray output = model.output(testdata, Layer.TrainingMode.TEST);
            for(int j = 0; j < output.rows(); j++)
            {
                temp_out.append(output.getDouble(j,0)+"\n");
            }
        }
        recordReader2.close();
        resultOut.write(temp_out.toString().getBytes("utf-8"));
        resultOut.close();
    }
    public static void main(String[] args){

    }
}
