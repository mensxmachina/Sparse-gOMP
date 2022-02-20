import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

public class Feature {
    public int id,num_of_features;
    public double correlation = 0,sum_for_corr = 0,square_sum_for_corr = 0,sum_with_target=0;
    public ArrayList<Double> data ;
    public ArrayList<Integer> rows ;
//    public ArrayList<Double> standardized_data = new ArrayList<Double>();
    double mean = 0,sd = 0;
    boolean standardized = false;
    boolean calculatesd = false;
    boolean calculatemean = false;
//    public ArrayList<Double> scores = new ArrayList<Double>();
//    public ArrayList<Double> betas = new ArrayList<Double>();
//    public HashMap<Integer,Integer> hashrows = new HashMap<Integer, Integer>();
    public Feature(int num_of_features,int id){
        this.num_of_features = num_of_features;
        this.id = id;
        data = new ArrayList<Double>(num_of_features);
        rows = new ArrayList<Integer>(num_of_features);
    }

    public ArrayList<Double> getData() {
        return data;
    }

    public ArrayList<Integer> getRows() {
        return rows;
    }

    public double getCorrelation() {
        return correlation;
    }

    public double getSquare_sum_for_corr() {
        return square_sum_for_corr;
    }

    public double getSum_for_corr() {
        return sum_for_corr;
    }

    public double getSum_with_target() {
        return sum_with_target;
    }

    public int getId() {
        return id;
    }

    public int getNum_of_features() {
        return num_of_features;
    }

    public double calculatemean(int num_of_samples){
        if(calculatemean == true)return this.mean;
        this.mean = 0;
        for(int i = 0;i<data.size();i++){
            this.mean+= data.get(i);
        }
        this.mean = this.mean/(double)num_of_samples;
        calculatemean = true;
        return this.mean;
    }

    public double calculatesd(int num_of_samples){
        if(calculatesd == true)return this.sd;
        this.sd = 0;
        double sum = 0;
        for(int i = 0;i<data.size();i++){
            sum+= data.get(i);
        }
//        if(mean==0){
//             mean = calculatemean(num_of_samples);
//        }
        double mean = sum/num_of_samples;
        for(int i = 0;i<data.size();i++){
            this.sd += Math.pow(data.get(i) - mean, 2);
        }
        this.sd += Math.pow(0-mean,2)*(num_of_samples-data.size());
//        this.sd+= Math.pow(0 - mean, 2)*(num_of_samples-data.size());
        this.sd = (double) Math.sqrt(this.sd/(num_of_samples));
        calculatesd = true;
        return this.sd;
    }

//    public void standardize_data(int num_of_samples){
//        if(standardized == false) {
//            if (mean == 0) {
//                mean = calculatemean(num_of_samples);
//            }
//            if (sd == 0) {
//                sd = calculatesd(num_of_samples);
//            }
//            for (int i = 0; i < data.size(); i++) {
//
//                standardized_data.add((data.get(i) - mean) / sd);
////            sum_for_corr+=standardized_data.get(i);
////            square_sum_for_corr+=standardized_data.get(i) * standardized_data.get(i);
//            }
//            standardized = true;
//        }
//    }
}
