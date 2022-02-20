import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;


public class OMP {

    public static void main(String args[]) {


        int numVars =16609144;
        int numSamples = 332500;
//        int targetcol =0;
//        double[] resids = new double[numSamples];
//        ArrayList<Feature> vars = new ArrayList<Feature>(1000001);
////        ArrayList<Feature> varstrans = new ArrayList<Feature>();
//        for(int j = 0;j<numVars;j++){
////            System.out.print(j+",");
//            vars.add(new Feature(0,j));
//        }

//        int[] cc = new int[numVars];
//        Arrays.fill(cc,0);
        int counter = 0;
        double number;
        String filepath = "./webspam_wc_normalized_unigram.svm";

        ArrayList<Feature> vars = new ArrayList<Feature>(numVars);
        for(int j = 0;j<numVars;j++){
//            System.out.print(j+",");
            vars.add(new Feature(0,j));
        }
        counter = 0;
        try {
//            BufferedReader br = new BufferedReader(new FileReader(new File(filepath)));
            LineIterator it = FileUtils.lineIterator(new File(filepath),"UTF-8");
            String CurrentLine = "";
//            String[] splitline = new String[numVars]; need this
            String[] s;
            String splitline[];
            int i = 0;
            int feature;

            while (it.hasNext() && counter <numSamples) {
                CurrentLine = it.nextLine();
//                i++;
//                if(i<6000000)continue;
                if (i == -1) {
                    i = 0;
                    continue;
                }

//                System.arraycopy(CurrentLine.split(","),0,splitline,0,numVars); need this



//                System.out.println(CurrentLine);
                splitline = CurrentLine.split(" ");
//                if(splitline.length<numVars)continue;
//                CurrentLine = "";
//                CurrentLine = null;


//                for (int k = 0; k < numVars; k++) { need this
//                System.out.println("len:"+splitline.length);
                for (int k = 0; k < splitline.length; k++){
//                    number = Double.parseDouble(splitline[k]); need this
//                    if (k != targetcol) { need this
                    if(k!=0){
                        s = splitline[k].split(":",2);
                        feature = Integer.parseInt(s[0]);
                        feature = feature-1;
                        number = Double.parseDouble(s[1]);
//                        System.out.println(cc[feature]);
                        if(number!=0) {

                            vars.get(feature).data.add(number); // k instead of feature
                            vars.get(feature).rows.add(counter);


                            vars.get(feature).num_of_features++;
                            vars.get(feature).sum_for_corr += number;
                            vars.get(feature).square_sum_for_corr += number * number;



                        }
                    } else {
//
                            number = Double.parseDouble(splitline[k]);
                            if(number == -1)number = 0;
//                            target[counter] = number;
                            vars.get(numVars-1).data.add(number);
                            vars.get(numVars-1).num_of_features++;
//                        }

                    }

                }

                counter++;

                System.out.println("iter:"+counter);
            }
            it.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


        double[] target = new double[numSamples];
        try {
            FileWriter w = new FileWriter("targetcol");

            for (int i = 0; i < numSamples; i++) {
                target[i] = vars.get(numVars - 1).data.get(i);
                w.write(target[i]+"\n");
            }
            w.close();
        }catch (IOException E){
            E.printStackTrace();
        }

//        for(int c = 0;c<50;c++) {
//            double y[] = new double[numSamples];
//            ArrayList<Double> cf = new ArrayList<>();
//            ArrayList<Integer> svar = new ArrayList<>();
//
//            Arrays.fill(y,0);
//            filepath = "coef"+(c+1);
//            try {
//
//                LineIterator it = FileUtils.lineIterator(new File(filepath), "UTF-8");
//                String CurrentLine = "";
//
//                String[] s;
//                String splitline[];
//                int i = 0;
//                int feature;
//
//                while (it.hasNext()) {
//                    CurrentLine = it.nextLine();
//
//                    if (i == -1) {
//                        i = 0;
//                        continue;
//                    }
//
//
//                    splitline = CurrentLine.split(" ");
//                    if(i==(c+1)){
//                        cf.add(-Double.parseDouble(splitline[0]));
//                    }else{
//                        cf.add(-Double.parseDouble(splitline[0]));
//                        svar.add(Integer.parseInt(splitline[1]));
//                    }
//
//
//
//                    i++;
//
//                    System.out.println("iter:" + i);
//                }
//                System.out.println(cf);
//                it.close();
//                FileWriter w = new FileWriter("y"+(c+1));
//                for(int l = 0;l<cf.size()-1;l++) {
//                    for(int m = 0;m<vars.get(svar.get(l)).rows.size();m++){
//                        y[vars.get(svar.get(l)).rows.get(m)]+= vars.get(svar.get(l)).data.get(m)* cf.get(l);
//                    }
//                }
//                for(int l = 0;l<numSamples;l++){
//                    y[l]+=cf.get(cf.size()-1);
//                    y[l] = Math.exp(y[l])/(1+Math.exp(y[l]));
//                    w.write(y[l]+"\n");
//                }
//
//                w.close();
//
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }

        for(int j = 0;j<numVars;j++){
            vars.get(j).rows.trimToSize();
            vars.get(j).data.trimToSize();
        }


        double starttime = System.nanoTime();
        OMPRet omp = OMPRet.sparseomp(vars,numVars,numSamples,target,"logistic",3,50);

        System.out.println(omp.selectedvars);

        double endtime = System.nanoTime();
        double totaltime = (endtime-starttime)/1000000000;
        System.out.println("total:"+totaltime);
    }
}
