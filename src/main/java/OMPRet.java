import org.omg.CORBA.OMGVMCID;

import java.lang.reflect.Array;
import java.util.*;

public class OMPRet {
    public ArrayList<Integer> selectedvars;
    public ArrayList<Boolean> selectedvarsbool;
    ArrayList<ArrayList<Integer>> equiv;
    public OMPRet(){
        selectedvars = new ArrayList<Integer>();
        selectedvarsbool = new ArrayList<Boolean>();
        equiv = new ArrayList<ArrayList<Integer>>();
    }
    static double correlationCoefficient(Feature X,
                                         double Y[], int n, Feature target)
    {

        double sum_X = 0, sum_Y = 0, sum_XY = 0;
        double squareSum_X = 0, squareSum_Y = 0;

        sum_Y = target.sum_for_corr;
        squareSum_Y = target.square_sum_for_corr;
        sum_X = X.sum_for_corr;
        squareSum_X = X.square_sum_for_corr;
//        for (int i = 0; i < n; i++)
//        {
//            sum_Y = sum_Y + Y[i];


        for(int i = 0;i<X.rows.size();i++) {
//            sum_X = sum_X + X.data.get(i);
            // sum of X[i] * Y[i].
            sum_XY = sum_XY + X.data.get(i) * Y[X.rows.get(i)];

            // sum of square of array elements.
//            squareSum_X = squareSum_X + X.data.get(i) * X.data.get(i);

        }
//            }else {
//                sum_X+= (0-X.mean)/X.sd;
//                sum_XY+= ((0-X.mean)/X.sd)* Y[i];
//                squareSum_X = squareSum_X + (0-X.mean)/X.sd*(0-X.mean)/X.sd;
//            }
//            squareSum_Y = squareSum_Y + Y[i] * Y[i];
//        }
        // use formula for calculating correlation
        // coefficient.
        double corr = (double)(n * sum_XY - sum_X * sum_Y)/
                (double)(Math.sqrt((n * squareSum_X -
                        sum_X * sum_X) * (n * squareSum_Y -
                        sum_Y * sum_Y)));

        return corr;
    }

    static double calculateBIC(double res,int dof,int n,String test){
        if(test.equals("linear") ) {
            return (double) (n * Math.log(2 * Math.PI / n * res) + n + dof * Math.log(n));
        }else if( test.equals("logistic")){
            return (double) (res + dof * Math.log(n));
        }else{
            return 0;
        }
    }


    public static OMPRet sparseomp(ArrayList<Feature> vars, int numVars, int numSamples, double[] target, String test, double DBIC, int maxF){
        double mean = vars.get(numVars-1).calculatemean(numSamples);
        OMPRet ret = new OMPRet();
//        double[] resids = new double[numSamples];
        int selectedvar = -1;
        double prevDev = 0;
        double currBIC = 0;
        double prevBIC = 0;
        OMPTest testres = new OMPTest(numSamples,1);
        for(int i = 0;i<numSamples;i++){
            testres.resids[i] = target[i] - mean;
            vars.get(numVars-1).sum_for_corr+=testres.resids[i];
            vars.get(numVars-1).square_sum_for_corr+=testres.resids[i]*testres.resids[i];
        }
        double maxcorr = 0;
        for(int i = 0;i<numVars-1;i++){
//                vars.get(i).standardize_data(numSamples);

            vars.get(i).correlation = Math.abs(correlationCoefficient(vars.get(i), testres.resids, numSamples,vars.get(numVars-1)));

            if(vars.get(i).correlation>maxcorr){
                maxcorr = vars.get(i).correlation;
                selectedvar = i;
            }
        }
        if(test.equals("linear")) {
            for (int i = 0; i < numSamples; i++) {
                prevDev += Math.pow(testres.resids[i], 2);
            }

        }else if (test.equals("logistic")){
            prevDev = (double) (-2*(numSamples*vars.get(numVars-1).mean*Math.log(vars.get(numVars-1).mean)+(numSamples-numSamples*vars.get(numVars-1).mean)*Math.log(1-vars.get(numVars-1).mean)));
        }
        prevBIC = calculateBIC(prevDev, 1, numSamples,test);
        ret.selectedvars.add(selectedvar);
        int iters = 0;





        while(ret.selectedvars.size()<numVars-1){
            System.out.println("iter: "+iters);
            if(test.equals("linear")){
//                    testres = OMPTest.LinearTestImproved(target,vars,varstrans,numVars,numSamples,ret.selectedvars,prevDev);
            }else if(test.equals("logistic")){
                testres = OMPTest.LogisticTestImproved(target,vars,numVars,numSamples,ret.selectedvars,prevDev);
            }

            currBIC = calculateBIC(testres.currDev,ret.selectedvars.size()+1,numSamples,test);
            System.out.println(currBIC);
            System.out.println(prevBIC);

            if(((prevBIC-currBIC)<DBIC|| (ret.selectedvars.size()>maxF) || (ret.selectedvars.size() == numVars))){
                ret.selectedvars.remove((Object)selectedvar);

                System.out.println("vgike");

//                backstepres = backwardstep(target,targetarr,dataset,numVars,numSamples,test,DBIC,results.selectedvars,results.equiv);
//                results.selectedvars = backstepres.selvars;
//                results.equiv = backstepres.equiv;
                return ret;
            }

            prevDev = testres.currDev;
            prevBIC = currBIC;
            maxcorr = 0;
            vars.get(numVars-1).sum_for_corr = 0;
            vars.get(numVars-1).square_sum_for_corr = 0;
            for(int i = 0;i<numSamples;i++){

                vars.get(numVars-1).sum_for_corr+= testres.resids[i];
                vars.get(numVars-1).square_sum_for_corr+= testres.resids[i]*testres.resids[i];
            }

            for(int i = 0;i<numVars-1;i++){

                if(!ret.selectedvars.contains(i)) {
                    vars.get(i).correlation = Math.abs(correlationCoefficient(vars.get(i), testres.resids, numSamples, vars.get(numVars - 1)));

                    if (vars.get(i).correlation > maxcorr) {
                        maxcorr = vars.get(i).correlation;
                        selectedvar = i;
                    }
                }
            }

            ret.selectedvars.add(selectedvar);
            iters++;
        }
        return ret;
    }

}
