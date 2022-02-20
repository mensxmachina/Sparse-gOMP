import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;


public class OMPTest {

    static ArrayList<Feature> X;
    static ArrayList<ArrayList<HashMap<Integer,Double>>> XmultX = new ArrayList<ArrayList<HashMap<Integer, Double>>>();
//    static ArrayList<double[]> lastcols = new ArrayList<double[]>();
//    static ArrayList<double[]> lastrows = new ArrayList<double[]>();
    static double[] prevcoef;
    double []resids;
    double currDev;
    static boolean hasconverged;
    int i;
    double[] be;
    static double negLL;
//    static double[] pleasedontcrush;
    public OMPTest(int numSamples,int numofselvars){
        this.resids = new double[numSamples];
        this.currDev = 0;
        i = 0;
        this.be = new double[numofselvars];
        prevcoef = new double[1];
        prevcoef[0] = 0;
//        pleasedontcrush = new double[numSamples];
    }

//    static double sparse_dot_prod(Feature X, Feature Y,int numSamples){
//
//        double sum_XY = 0;
//        int counter = 0;
//        double standardX = (0-X.mean)/X.sd;
//        double standardY = (0-Y.mean)/Y.sd;
//        ArrayList<Integer> common = new ArrayList<Integer>(X.hashrows.keySet());
//        ArrayList<Integer> rem = new ArrayList<Integer>(Y.hashrows.keySet());
//        common.retainAll(rem);
//
//        counter+= common.size();
//        int Xindex;
//        int Yindex;
//        for(int i = 0;i<common.size();i++){
//            Xindex = X.hashrows.get(common.get(i));
//            Yindex = Y.hashrows.get(common.get(i));
//            sum_XY+= X.standardized_data.get(Xindex)*Y.standardized_data.get(Yindex);
//        }
//        ArrayList<Integer> uncommon = new ArrayList<Integer>(X.hashrows.keySet());
//
//        uncommon.removeAll(rem);
//
//        counter+=uncommon.size();
//        for(int i = 0;i<uncommon.size();i++){
//            Xindex = X.hashrows.get(uncommon.get(i));
//
//            sum_XY+= standardY*X.standardized_data.get(Xindex);
//        }
//        uncommon = new ArrayList<Integer>(Y.hashrows.keySet());
//        rem = new ArrayList<Integer>(X.hashrows.keySet());
//        uncommon.removeAll(rem);
//
//        counter+=uncommon.size();
//        for(int i = 0;i<uncommon.size();i++){
//            Yindex = Y.hashrows.get(uncommon.get(i));
//            sum_XY+= standardX*Y.standardized_data.get(Yindex);
//        }
//        sum_XY += standardY*standardX*(numSamples-counter);
//
//        return sum_XY;
//
//    }
//
//    static double dot_prod_ones(Feature X,int numSamples){
//        double sumXY = 0;
//        int dataiter = 0;
//        double standardX = (0-X.mean)/X.sd;
//        for(int i = 0;i<numSamples;i++){
//            if(X.hashrows.containsKey(i)){
//                sumXY+=X.standardized_data.get(dataiter);
//
//                dataiter++;
//            }else{
//                sumXY+= standardX;
//
//            }
//        }
//
//        return sumXY;
//    }
//
//    static double sparse_dense_dot(ArrayList<Feature> vars, ArrayList<Integer> selvars, double[] M, int currcol){
//        double dotprod = 0;
//        for(int i = 0;i<selvars.size()+1;i++){
//            if(i == 0){
//                dotprod+=M[i];
//            }else{
//                if(vars.get(selvars.get(i-1)).hashrows.containsKey(currcol)){
//                    int index = vars.get(selvars.get(i-1)).hashrows.get(currcol);
//                    dotprod+= vars.get(selvars.get(i-1)).standardized_data.get(index)*M[i];
//                }else{
//                    dotprod+= ((0-vars.get(selvars.get(i-1)).mean)/vars.get(selvars.get(i-1)).sd) * M[i];
//                }
//            }
//        }
//
//        return dotprod;
//    }


//    static double[][] mat_mult_sparse(double[][] M, ArrayList<Feature> X, ArrayList<Integer> selvars, ArrayList<Feature> vars, int numSamples){
//        double[][] A = new double[M.length][numSamples];
//
//
//        for(int i = 0;i<M.length;i++){
//            for(int j = 0;j<numSamples;j++){
//
//
//                    A[i][j] = sparse_dense_dot(vars,selvars, M[i],j);
//
//            }
//        }
//        return A;
//    }

    public static double[][] invert(double a[][])
    {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i)
            b[i][i] = 1;

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                            -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }

// Method to carry out the partial-pivoting Gaussian
// elimination.  Here index[] stores pivoting order.

    public static void gaussian(double a[][], int index[])
    {
        int n = index.length;
        double c[] = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            double c1 = 0;
            for (int j=0; j<n; ++j)
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i)
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1)
                {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }


//    public static OMPTest LinearTestImproved(double[] target, ArrayList<Feature> vars,ArrayList<Feature>varstrans,int numVars,int numSamples,ArrayList<Integer> selvars,double prevDev){
//        OMPTest retval = new OMPTest(numSamples,selvars.size());
//
//
//
////        System.out.println(selvars);
////        int lastselvar = selvars.get(selvars.size()-1);
////        ArrayList<Integer> sortedselvars = new ArrayList<Integer>(selvars);
////        Collections.sort(sortedselvars);
////        System.out.println(sortedselvars);
////        int ind = sortedselvars.indexOf(lastselvar);
////        System.out.println(ind);
//        vars.get(selvars.get(selvars.size()-1)).standardize_data(numSamples);
//
//        X.add(vars.get(selvars.get(selvars.size()-1)));
////        X.add(ind,vars.get(selvars.get(selvars.size()-1)));
//
//
//
//
//
//        double lastrow = 0;
//        double[] llrow = new double[selvars.size()+1];
//        double[] lastcol = new double[selvars.size()];
//        double starttime = System.nanoTime();
//        for(int i = 0;i<selvars.size()+1;i++){
//            if(i==0){
//                lastcol[i] = dot_prod_ones(vars.get(selvars.get(selvars.size()-1)),numSamples);
//
//                llrow[i] = lastcol[i];
//
//            }else if(i<selvars.size()) {
//
//                lastcol[i] = sparse_dot_prod(vars.get(selvars.get(i-1)), vars.get(selvars.get(selvars.size() - 1)), numSamples);
//
//                llrow[i] = lastcol[i];
//
//            }else{
//                lastrow = sparse_dot_prod(vars.get(selvars.get(selvars.size() - 1)),vars.get(selvars.get(selvars.size() - 1)), numSamples);
//
//                llrow[i] = lastrow;
//            }
//        }
//        double endtime = System.nanoTime();
//        double totaltime = (endtime-starttime)/1000000000;
//        System.out.println("sparse dot prods:"+totaltime);
//
//        lastcols.add(lastcol);
//        lastrows.add(llrow);
//        // AN VALUME STO SOSTO INDEX TO LASTCOL KAI ROW TOTE MPOREI TA PROTA CALCULATIONS NA EXUN PIO POLA
//        // STOIXEIA APTA TELEYTEA LOGO TOY SORTING KAI TOTE O XTX GINETAI OPOS NANAI, DEN ISXIEI PIA TO
//        // ADD STO TELOS TOY PINAKA
//        int lastlength = 1;
//        int counter = 0;
//        double[] temp;
//        double[] temp2;
//        double[][] XTX_array = new double[selvars.size()+1][selvars.size()+1];
//
////        for(int i = 0;i<selvars.size()+1;i++){
////            for(int j = 0;j<selvars.size()+1;j++){
////                if(i == 0 && j == 0){
////                    XTX_array[0][0] = numSamples;
////                }else if(i == 0){
////                    XTX_array[i][j] = dot_prod_ones(X.get(j-1),numSamples);
////                }else if(j == 0){
////                    XTX_array[i][j] = dot_prod_ones(X.get(i-1),numSamples);
////                }
////                else{
////                    XTX_array[i][j] = sparse_dot_prod(X.get(i-1),X.get(j-1),numSamples);
////                }
////            }
////        }
//
//
//        XTX_array[0][0] = numSamples;
//        starttime = System.nanoTime();
//        while(counter < lastcols.size()){
//
//            temp = lastrows.get(counter).clone();
//            temp2 = lastcols.get(counter).clone();
//            for(int i = 0;i<temp2.length;i++){
//                XTX_array[i][lastlength] = temp2[i];
//            }
//            for(int i = 0;i<temp.length;i++){
//                XTX_array[lastlength][i] = temp[i];
//            }
//
//            counter++;
//            lastlength = temp.length;
//        }
//        endtime = System.nanoTime();
//        totaltime = (endtime-starttime)/1000000000;
//        System.out.println("xtx creation:"+totaltime);
//
////        for(int i = 0;i<XTX_array.length;i++){
////            System.out.println(Arrays.toString(XTX_array[i]));
////        }
//
//
//
//
//        RealMatrix XTX_ = new Array2DRowRealMatrix(XTX_array);
//
//
////        System.out.println(XTX_);
//        starttime = System.nanoTime();
//        RealMatrix XTXinv = new CholeskyDecomposition(XTX_).getSolver().getInverse();
//        endtime = System.nanoTime();
//        totaltime = (endtime-starttime)/1000000000;
//        System.out.println("xtxinv cre:"+totaltime);
//
//
//        double[][] XTX_inv = XTXinv.getData();
////        double[][] XTX_inv = invert(XTX_array);
//
////        System.out.println("--");
////        for(int i = 0;i<XTX_array.length;i++){
////            System.out.println(Arrays.toString(XTX_inv[i]));
////        }
//
//        starttime = System.nanoTime();
//        double[][] XTXX = mat_mult_sparse(XTX_inv,X,selvars,vars,numSamples);
//        endtime = System.nanoTime();
//        totaltime = (endtime-starttime)/1000000000;
//        System.out.println("mat mult sparse:"+totaltime);
//
//
//
//
//
//
//        double[] bhat;
//        RealMatrix XTXX_ = new Array2DRowRealMatrix(XTXX);
//
//        starttime = System.nanoTime();
//        double[] target2 = new double[numSamples];
//        for(int i = 0;i<numSamples;i++){
//            target2[i] = target[i];
//        }
//
//        bhat = XTXX_.operate(target2);
//
//
////        System.out.println("bhat");
////        System.out.println(Arrays.toString(bhat));
//
//
//        double[] XopB = new double[numSamples];
//        double sum;
//        for(int i = 0;i<numSamples;i++){
//            sum = 0;
//            for(int j = 0;j<selvars.size()+1;j++){
//                if(j == 0){
//                    sum += bhat[j];
//                }else{
//                    if(X.get(j-1).hashrows.containsKey(i)){
//                        int index = X.get(j-1).hashrows.get(i);
//                        sum+= X.get(j-1).standardized_data.get(index)*bhat[j];
//                    }else{
//                        double mean = (double) (((0-X.get(j-1).mean)/X.get(j-1).sd)*bhat[j]);
//                        sum+= mean;
//                    }
//                }
//            }
//            XopB[i] = sum;
//
//        }
//        endtime = System.nanoTime();
//        totaltime = (endtime-starttime)/1000000000;
//        System.out.println("bhat operations:"+totaltime);
//        double[] resids2 = new double[numSamples];
//        for(int i = 0;i<numSamples;i++){
//            retval.resids[i] = target[i] - XopB[i];
//            resids2[i] = retval.resids[i];
//
//        }
//        RealVector res = new ArrayRealVector(resids2);
//
//        retval.currDev = (double) res.dotProduct(res);
////        System.out.println(retval.currDev);
//
//
//        return retval;
//    }

//    public static  OMPTest LinearTest(double[] target, ArrayList<Feature> vars,int numVars,int numSamples,ArrayList<Integer> selvars,double prevDev){
//        OMPTest retval = new OMPTest(numSamples,selvars.size());
//        int maxrows = 0;
//        int currrows = 0;
//
//        for(int i = 0;i<selvars.size();i++){
//            currrows = vars.get(selvars.get(i)).rows.get(vars.get(selvars.get(i)).rows.size()-1);
//            if(maxrows<currrows){
//                maxrows = currrows;
//            }
//        }
//        double[][] X = new double[maxrows][selvars.size()+1];
//        double[] newtarget = new double[maxrows];
//        double[] multresult = new double[maxrows];
//        for(int i = 0;i<maxrows;i++){
//            X[i][0] = 1;
//            newtarget[i] = target[i];
//        }
//        int dataiter;
//
//        for(int i = 0;i<selvars.size();i++){
//            dataiter = 0;
//            vars.get(selvars.get(i)).standardize_data(numSamples);
//            for(int j = 0;j<maxrows;j++){
//                if(vars.get(selvars.get(i)).rows.contains(j)){
//                    X[j][i+1] = vars.get(selvars.get(i)).standardized_data.get(dataiter);
//                    dataiter++;
//                }else{
//                    X[j][i+1] = (0-vars.get(selvars.get(i)).mean)/vars.get(selvars.get(i)).sd;
//                }
//            }
//        }
//
////        for(int i = 0;i<10;i++){
////            for(int j = 0;j<selvars.size()+1;j++){
////                System.out.println(X[i][j]);
////            }
////            System.out.println();
////        }
//        double sum = 0;
//        try{
//            RealMatrix X2 = new Array2DRowRealMatrix(X);
//            DecompositionSolver solver = new SingularValueDecomposition(X2).getSolver();
//            RealVector target2 = new ArrayRealVector(newtarget);
//            RealVector betas = solver.solve(target2);
//            System.out.println(betas);
//            double[] betas2 = betas.toArray();
//            for (int i = 0; i < maxrows; i++) {
//                for (int j = 0; j < selvars.size()+1; j++) {
//                    sum += X[i][j] * betas2[j];
//                }
//                multresult[i] = sum;
//
//                sum = 0;
//            }
//            for (int i = 0; i < maxrows; i++) {
//                retval.resids[i] = target[i] - multresult[i];
//
//                retval.currDev += Math.pow(retval.resids[i], 2);
//            }
//            for(int i = maxrows;i<numSamples;i++){
//                retval.resids[i] = target[i]-vars.get(numVars-1).mean;
//                retval.currDev += Math.pow(retval.resids[i], 2);
//            }
//
//        }catch (Exception e){
//            System.out.println(e);
//            double [] betas2 = new double[selvars.size()];
//            for(int i = 0;i<selvars.size();i++){
//                betas2[i] = 0;
//            }
//            retval.currDev = prevDev;
//            for(int i = 0;i<numSamples;i++){
//                retval.resids[i] = Math.pow(10,6);
//            }
//        }
//
//        return retval;
//    }

    public static double [][] cross_x_y(double [][] x,double [][] y, int xcols,int xrows,int ycols){
        double[][] f = new double[xcols][ycols];
        double sum  = 0;
        for(int i =0;i<ycols;i++){
            for(int j = 0;j<xcols;j++){
                for(int k = 0;k<xrows;k++){
                    sum+= y[k][i] * x[k][j];
                }
                f[j][i] = sum;
                sum = 0;
            }
        }
        return f;
    }

    public static double calcDevRes(double [] p, double[] y, double[] expyhat) {
        int psize = p.length;
        double s = 0;
        for (int i = 0; i < psize; i++) {
//            System.out.println(p[i]);
            if (y[i] == 1.0) {
                if (p[i] == 0.0) {
                    s += expyhat[i];
                } else {
                    s += Math.log(p[i]);
                }
            } else {
                if (p[i] == 1.0) {
                    s += expyhat[i];
                } else {
                    s += Math.log(1 - p[i]);
                }
            }

        }
        return s;
    }

    public static  OMPTest LogisticTestImproved(double[] target, ArrayList<Feature> vars,int numVars,int numSamples,ArrayList<Integer> selvars,double prevDev) {

        X = vars;

        double[] coef = new double[prevcoef.length+1];

        System.arraycopy(prevcoef, 0, coef, 0, prevcoef.length - 1);

        coef[coef.length - 1] = prevcoef[prevcoef.length - 1];

        OMPTest ret =  fit(target,1.0E-4,100,1.0E-4,0.5,vars,selvars,coef);
        prevcoef = coef;
        return ret;
    }

//    public static double computeNegativeLogLikelihood(double[] coef, double[] y, ArrayList<Feature> X) {
//        double LL = 0;
//        final int nsamples = y.length;
//        final int constant = X.size();
//
//        for (int i = 0; i < nsamples; ++i) {
//            double wx = 0;
//            for (int n = 0; n < (X.size()); ++n) {
//                if(X.get(n).hashrows.containsKey(i)){
//                    int index = X.get(n).hashrows.get(i);
//                    wx += X.get(n).data.get(index) * coef[n];
//                }else{
//                    wx += ((0-X.get(n).mean)/X.get(n).sd) * coef[n];
//                }
//            }
//            // constant
//            wx += coef[constant];
////            LL += (1 - y[i]) * wx - FastMath.log(1 + FastMath.exp(wx));
//            LL += (1 - y[i]) * wx - ((wx <= 30) ? Math.log(1 + FastMath.exp(wx)) : wx);
//        }
//
//        return -LL;
//    }

    public static double computeNegativeLogLikelihoodimproved(double[] coef, double[] y,ArrayList<Integer> selvars) {
        double LL = 0;
        final int nsamples = y.length;
        final int constant = selvars.size();
        double[] wx =  new double[nsamples];
        Arrays.fill(wx,0);
//        double sumb0 = 0;
//        for(int i = 0;i<coef.length;i++){
//            if(i == coef.length-1) {
//                for(int j = 0;j<coef.length-1;j++){
//                    sumb0 += coef[j]*X.get(selvars.get(j)).calculatemean(nsamples)/X.get(selvars.get(j)).calculatesd(nsamples);
//                }
//                coef[i] = -coef[i]+sumb0;
//            }else{
//
//                coef[i] = -coef[i]*X.get(selvars.get(i)).calculatesd(nsamples);
//            }
//        }
        for(int n = 0;n<selvars.size();n++){
            for(int i = 0;i<X.get(selvars.get(n)).rows.size();i++){
                wx[X.get(selvars.get(n)).rows.get(i)]+=X.get(selvars.get(n)).data.get(i)*coef[n];
            }
        }
        for(int i = 0;i<nsamples;i++){
            wx[i]+= coef[constant];
            LL += (1 - y[i]) * wx[i] - ((wx[i] <= 30) ? Math.log(1 + FastMath.exp(wx[i])) : wx[i]);
        }
        return -LL;
    }


    public static double computeNegativeLogLikelihoodimproved(double[] coef, double t, double[] delta, double[] y, ArrayList<Integer> selvars) {
        double LL = 0;
        final int nsamples = y.length;
        final int constant = selvars.size();
        double[] wx =  new double[nsamples];
//        double sumb0 = 0;
//        for(int i = 0;i<coef.length;i++){
//            if(i == coef.length-1) {
//                for(int j = 0;j<coef.length-1;j++){
//                    sumb0 += coef[j]*X.get(selvars.get(j)).calculatemean(nsamples)/X.get(selvars.get(j)).calculatesd(nsamples);
//                }
//                coef[i] = -coef[i]+sumb0;
//            }else{
//
//                coef[i] = -coef[i]*X.get(selvars.get(i)).calculatesd(nsamples);
//            }
//        }
        Arrays.fill(wx,0);
        for(int n = 0;n<selvars.size();n++){
            for(int i = 0;i<X.get(selvars.get(n)).rows.size();i++){
                wx[X.get(selvars.get(n)).rows.get(i)]+=X.get(selvars.get(n)).data.get(i)*(coef[n]+t*delta[n]);
            }
        }
        for(int i = 0;i<nsamples;i++){
            wx[i]+= (coef[constant] + t * delta[constant]);
            LL += (1 - y[i]) * wx[i] - ((wx[i] <= 30) ? Math.log(1 + FastMath.exp(wx[i])) : wx[i]);
        }
        return -LL;
    }

//    public static double computeNegativeLogLikelihood(double[] coef, double t, double[] delta, double[] y, ArrayList<Feature> X) {
//        double LL = 0;
//        final int nsamples = y.length;
//        final int constant = X.size();//maybe xoris to +1
//
//
//        for (int i = 0; i < nsamples; ++i) {
//            double wx = 0;
//            for (int n = 0; n < (X.size()); ++n) {
//                if(X.get(n).hashrows.containsKey(i)){
//                    int index = X.get(n).hashrows.get(i);
//                    wx += X.get(n).data.get(index) * (coef[n] + t * delta[n]);
//                }else{
//                    wx += ((0-X.get(n).mean)/X.get(n).sd) * (coef[n] + t * delta[n]);
//                }
//
//            }
//            // constant
//            wx += (coef[constant] + t * delta[constant]);
////            LL += (1 - y[i]) * wx - FastMath.log(1 + FastMath.exp(wx));
//            LL += (1 - y[i]) * wx - ((wx <= 30) ? Math.log(1 + FastMath.exp(wx)) : wx);
//        }
//
//        return -LL;
//    }

//    public static double[] vectormult (Feature X,Feature Y,int numsamples){
//        ArrayList<Integer> allsamples = new ArrayList<Integer>();
//        for(int i = 0;i<numsamples;i++){
//            allsamples.add(i);
//        }
//        allsamples.removeAll(X.rows);
//        allsamples.removeAll(Y.rows);
//        double standardX = (0-X.mean)/X.sd;
//        double standardY = (0-Y.mean)/Y.sd;
//        double multXY = standardY*standardX;
//        double[] ret = new double[numsamples];
//        for(int i = 0;i<numsamples;i++){
//            if(allsamples.contains(i)){
//                ret[i] = multXY;
//            }else if(X.hashrows.containsKey(i) && Y.hashrows.containsKey(i)){
//                int indexX = X.hashrows.get(i);
//                int indexY= Y.hashrows.get(i);
//                ret[i] = X.data.get(indexX)*Y.data.get(indexY);
//
//            }else if(X.hashrows.containsKey(i)){
//                int indexX = X.hashrows.get(i);
//                ret[i] = X.data.get(indexX)*standardY;
//
//            }else if (Y.hashrows.containsKey(i)){
//                int indexY= Y.hashrows.get(i);
//                ret[i] = Y.data.get(indexY)*standardX;
//            }
//        }
//        return ret;
//    }
//
//    public static double[] badvectormult (Feature X,Feature Y,int numsamples){
//
//        double standardX = (0-X.mean)/X.sd;
//        double standardY = (0-Y.mean)/Y.sd;
//        double multXY = standardY*standardX;
//        double[] ret = new double[numsamples];
//        for(int i = 0;i<numsamples;i++){
//            if(X.hashrows.containsKey(i) && Y.hashrows.containsKey(i)){
//                int indexX = X.hashrows.get(i);
//                int indexY= Y.hashrows.get(i);
//                ret[i] = X.data.get(indexX)*Y.data.get(indexY);
//
//            }else if(X.hashrows.containsKey(i)){
//                int indexX = X.hashrows.get(i);
//                ret[i] = X.data.get(indexX)*standardY;
//
//            }else if (Y.hashrows.containsKey(i)){
//                int indexY= Y.hashrows.get(i);
//                ret[i] = Y.data.get(indexY)*standardX;
//            }else{
//                ret[i] = multXY;
//            }
//        }
//        return ret;
//    }


    public static HashMap<Integer,Double> badvectormultimproved (int X1,int Y,int numsamples){

//        double standardX = (0-X.mean)/X.sd;
//        double standardY = (0-Y.mean)/Y.sd;
//        double multXY = standardY*standardX;
        HashMap<Integer,Double> ret = new HashMap<>();
//        Arrays.fill(ret,0);
//        int size = 0;
//        Feature searchhere,searchfrom;
        HashMap<Integer, Double> Xhash = new HashMap<Integer, Double>();
//        HashMap<Integer, Double> Yhash = new HashMap<Integer, Double>();

        //TODO NA TA VALO SE ENA HASHMAP KAI AN IPARXEI IDI TO KEY TOTE TO VALUE NA GINETAI O POLLPLASIASMOS METAKSI TOUS
        for(int i = 0;i<X.get(X1).rows.size();i++){
            Xhash.put(X.get(X1).rows.get(i),X.get(X1).data.get(i));
        }
//        for(int i = 0;i<X.get(Y).rows.size();i++) {
//            Yhash.put(X.get(Y).rows.get(i), X.get(Y).data.get(i));
//        }
//        Set<Integer> keys =  Xhash.keySet();
        int row;
        double num;
        for(int i = 0;i<X.get(Y).rows.size();i++){
//            System.out.println("i:"+i);
            row = X.get(Y).rows.get(i);
            num = X.get(Y).data.get(i);
            if(Xhash.containsKey(row)){
//                System.out.println("in");
//                int indexX = X.hashrows.get(i);
//                int indexY= Y.hashrows.get(i);
                ret.put(row, Xhash.get(row)*num);
            }

        }
        return ret;
    }




    public static  OMPTest fit(double[] y, double tol, int maxit, double alpha, double beta, ArrayList<Feature> vars, ArrayList<Integer> selvars, double[] coef) {
        final int nsamples = y.length;
        final int constant = selvars.size();
        OMPTest ret = new OMPTest(nsamples,selvars.size()+1);
//        for(int i = 0;i<X.size();i++){
//            System.out.println(X.get(i).id);
////            System.out.println(X.get(i).data);
//        }

//        this.coef = new double[X.length + 1];
        final double[] gradient = new double[constant + 1];
        final ArrayRealVector gradientVector = new ArrayRealVector(gradient, false);

        final double[][] hessian = new double[constant + 1][constant + 1];
        final Array2DRowRealMatrix hessianMatrix = new Array2DRowRealMatrix(hessian, false);

        final double[] p = new double[nsamples];
        final double[] pp = new double[nsamples];
//        double totaltimepre = 0;
//        double totaltimegrad = 0;
//        double hessiantimeg = 0;
//        double ludtime = 0;
//        double compneglltime = 0;
        double lambda;
        int it = 0;
        negLL = computeNegativeLogLikelihoodimproved(coef, y,selvars);
        double nextNegLL;
//        ArrayList<ArrayList<double[]>> XmultX = new ArrayList<ArrayList<double[]>>();
//        for (int i1 = 0; i1 < X.size(); ++i1) {
//            Feature Xi1 = X.get(i1);
//            ArrayList<double[]> temp = new ArrayList<double[]>();
//            for (int i2 = i1; i2 < X.size(); ++i2) {
////                System.out.println(i1+" "+i2);
//                Feature Xi2 = X.get(i2);
//
//                temp.add(badvectormult(Xi1, Xi2, nsamples));
//
//            }
//            XmultX.add(temp);
//
//        }
        //TODO O XMULTX MPOREI NA KRATAEI MESA ARRAYLISTS KAI OXI DOUBLE PINAKES AFOU POLLA APTOUS PINAKES EINAI ZEROS
        int Xi2 = selvars.get(selvars.size()-1);
        for (int i1 = 0; i1 < selvars.size(); ++i1) {
            int Xi1 = selvars.get(i1);
            if(i1 == selvars.size()-1) {
                ArrayList<HashMap<Integer,Double>> temp = new ArrayList<>();

                temp.add(badvectormultimproved(Xi1, Xi2, nsamples));
                XmultX.add(temp);
            }else {

                XmultX.get(i1).add(badvectormultimproved(Xi1, Xi2, nsamples));
            }









        }


//        System.out.println(XmultX.get(0).get(0)[0]);
//        System.out.println("-----");
        double wx[] = new double[nsamples];

        do {
//            double starttime = System.nanoTime();
            // Do pre-computations
            Arrays.fill(wx,0);
            for(int n = 0;n<selvars.size();n++){
                for(int i = 0;i<X.get(selvars.get(n)).rows.size();i++){
                    wx[X.get(selvars.get(n)).rows.get(i)]+=X.get(selvars.get(n)).data.get(i) * coef[n];
                }
            }
            for(int i = 0;i<nsamples;i++){
                wx[i]+= coef[constant];
                p[i] = (double) (1 - (1 / (1 + FastMath.exp(wx[i]))));
                pp[i] = p[i] * (1 - p[i]);
            }
//            for (int i = 0; i < nsamples; ++i) {
//                double wx = 0;
//                for (int n = 0; n < X.size(); ++n) {
//                    if(X.get(n).hashrows.containsKey(i)){
//                        int index = X.get(n).hashrows.get(i);
//                        wx += X.get(n).data.get(index) * coef[n];
//                    }else{
//                        wx += ((0-X.get(n).mean)/X.get(n).sd) * coef[n];
//                    }
//                }
//                // constant
//                wx += coef[constant];
//
//                p[i] = 1 - (1 / (1 + FastMath.exp(wx)));
//                pp[i] = p[i] * (1 - p[i]);
//            }
//            System.out.println("Pi");
//            System.out.println(Arrays.toString(p));

            
//            double endtime = System.nanoTime();
//            double totaltime = (endtime-starttime)/1000000000;
//            totaltimepre+=totaltime;



            // Compute Gradient
//            starttime = System.nanoTime();
            for (int n = 0; n < selvars.size(); ++n) {
//                double standardX = (0-X.get(n).mean)/X.get(n).sd;
                gradient[n] = 0;
                for (int i = 0; i < X.get(selvars.get(n)).rows.size(); ++i) {

                    gradient[n] += (1 - y[X.get(selvars.get(n)).rows.get(i)] - p[X.get(selvars.get(n)).rows.get(i)]) * X.get(selvars.get(n)).data.get(i);

                }

            }
            // Constant
            gradient[constant] = 0;
            for (int i = 0; i < nsamples; ++i) {
                gradient[constant] += (1 - y[i] - p[i]);
            }
//            endtime = System.nanoTime();
//            totaltime = (endtime-starttime)/1000000000;
//            totaltimegrad +=totaltime;

            // Compute Hessian
//            starttime = System.nanoTime();
            for (int i1 = 0; i1 < XmultX.size(); ++i1) {
//                Feature Xi1 = X.get(i1);
                int i3 = i1;
                for (int i2 = 0; i2 < XmultX.get(i1).size(); ++i2) {
//                    Feature Xi2 = X.get(i2);

                    double s = 0;
//                    double[] XX = badvectormult(Xi1,Xi2,nsamples);
//                    System.out.println(XX[0]);
//                    System.out.println(i1+" "+i2);
                    HashMap<Integer,Double> XX = XmultX.get(i1).get(i2);
                    Set<Integer> keys =  XX.keySet();
                    for (int i:keys) {

                        s -= pp[i] * XX.get(i);
                    }
                    hessian[i3][i1] = hessian[i1][i3] = s;
                    i3++;
                }
            }
            // Constant with others
            for (int j = 0; j < selvars.size(); ++j) {
                double s = 0;
                for (int i = 0; i < X.get(selvars.get(j)).rows.size(); ++i) {

                    s -= pp[X.get(selvars.get(j)).rows.get(i)] * X.get(selvars.get(j)).data.get(i);


                }
                hessian[j][constant] = hessian[constant][j] = s;
            }
            // Constant with constant
            double s = 0;
            for (int i = 0; i < nsamples; ++i) {
                s -= pp[i];
            }
            hessian[constant][constant] = s;
//            endtime = System.nanoTime();
//            totaltime = (endtime-starttime)/1000000000;
//            hessiantimeg+=totaltime;
            // Compute step
            double[] delta;
            try {
//                starttime = System.nanoTime();
                delta = new LUDecomposition(hessianMatrix).getSolver().solve(gradientVector.mapMultiply(-1)).toArray();
//                endtime = System.nanoTime();
//                totaltime = (endtime-starttime)/1000000000;
//                ludtime = totaltime;
            } catch (Exception e) {
                // If Hessian is not invertible, we switch to simple gradient descent
                delta = gradientVector.toArray();


            }

            lambda = 0;
            for (int i = 0; i < gradient.length; ++i) {
                lambda += gradient[i] * delta[i];
            }

            // Backtracking line search
            double t = 1;
//            starttime = System.nanoTime();
            // Because we minimize, f(x + tDx) > f(x) + a*t*Dx*f(x)'
            while ((nextNegLL = computeNegativeLogLikelihoodimproved(coef, t, delta, y,selvars)) > negLL + alpha * t * lambda) {
                t *= beta;
            }
//            endtime = System.nanoTime();
//            totaltime = (endtime-starttime)/1000000000;
//            compneglltime+=totaltime;

            // Update delta and lambda in case t != 1
            if (t != 1) {
                for (int i = 0; i < delta.length; ++i) {
                    delta[i] *= t;
                }
                lambda = 0;
                for (int i = 0; i < gradient.length; ++i) {
                    lambda += gradient[i] * delta[i];
                }
            }

            // Update coefficients
//            System.out.println(Arrays.toString(coef));
            for (int i = 0; i < coef.length; ++i) {

                coef[i] += delta[i];

            }

            // Compute negative log-likelihood
            negLL = nextNegLL;//computeNegativeLogLikelihood(coef, y, X);

            ++it;
//            System.out.println("it:"+it);
            // Check conditions
            if (Double.isInfinite(negLL) || Double.isNaN(negLL)) {
                return null;
            }
        } while (lambda / 2 >= tol && it < maxit);
        if (lambda / 2 < tol) {
            hasconverged = true;
        }
//        double sumb0 = 0;
//        for(int i = 0;i<coef.length;i++){
//            if(i == coef.length-1) {
//                for(int j = 0;j<coef.length-1;j++){
//                    sumb0 += coef[j]*X.get(selvars.get(j)).calculatemean(nsamples)/X.get(selvars.get(j)).calculatesd(nsamples);
//                }
//                coef[i] = -coef[i]+sumb0;
//            }else{
//
//                coef[i] = -coef[i]*X.get(selvars.get(i)).calculatesd(nsamples);
//            }
//        }

        ret.currDev = 2*negLL;
        try {
            FileWriter w = new FileWriter("coef"+selvars.size());
            for(int i = 0;i<coef.length;i++){
                if(i<coef.length-1) {
                    w.write(coef[i] + " " + selvars.get(i) + "\n");
                }else{
                    w.write(coef[i] + "\n");
                }
            }
            w.close();
        }catch (IOException e){
            e.printStackTrace();
        }

        for (int i = 0; i < nsamples; i++) {
            ret.resids[i] = p[i] - 1 + y[i];

        }

        return ret;
//        System.out.println("Finished after iteration " + it);
    }

//    public static  OMPTest LogisticTest(double[] target, ArrayList<Feature> vars,int numVars,int numSamples,ArrayList<Integer> selvars,double prevDev) {
//        OMPTest retval = new OMPTest(numSamples,selvars.size()+1);
//        double my = vars.get(numVars-1).mean;
//        double d1 = (double) (numSamples * my * Math.log(my) + (numSamples - numSamples * my) * Math.log(1 - my));
//        double d2;
//        double lambda[];
//        int maxrows = 0;
//        int currrows = 0;
//        vars.get(numVars-1).betas.clear();
//        vars.get(numVars-1).scores.clear();
//        for(int i = 0;i<selvars.size();i++){
//            currrows = vars.get(selvars.get(i)).rows.get(vars.get(selvars.get(i)).rows.size()-1);
//            if(maxrows<currrows){
//                maxrows = currrows;
//            }
//        }
//
//
//        double[][] X= new double[maxrows][selvars.size()+1];
//        for(int i=0;i<maxrows;i++){
//            for(int j = 0;j<selvars.size()+1;j++){
//                X[i][j] = 0;
//            }
//        }
//        for(int i = 0;i<maxrows;i++){
//            X[i][0] = 1;
//        }
//
//        int dataiter;
//        for(int i = 0;i<selvars.size();i++){
//            dataiter = 0;
//            vars.get(selvars.get(i)).standardize_data(numSamples);
//            for(int j = 0;j<maxrows;j++){
//                if(vars.get(selvars.get(i)).rows.contains(j)){
//                    X[j][i+1] = vars.get(selvars.get(i)).standardized_data.get(dataiter);
//                    dataiter++;
//                }else{
//                    X[j][i+1] = (0-vars.get(selvars.get(i)).mean)/vars.get(selvars.get(i)).sd;
//                }
//            }
//        }
//        double be[] = new double[selvars.size()+1];
//        be[0] = Math.log(my) - Math.log(1 - my);
//
//        double[] y2 = new double[numSamples];
//        for (int i = 0; i < numSamples; i++) {
//            y2[i] = target[i] - my;
//        }
//        double sum = 0.0;
//        double[] der = new double[selvars.size()+1];
//        for (int i = 0; i < selvars.size()+1; i++) {
//            for (int j = 0; j < maxrows; j++) {
//                sum += X[j][i] * y2[j];
//            }
//            //TODO ADD FOR NUMSAMPLES AND NOT MAXROWS
//            der[i] = sum;
//            sum = 0.0;
//        }
//        double[][] der2;
//        double[][] x2 = new double[maxrows][selvars.size()+1];
//        for (int i = 0; i < maxrows; i++) {
//            for (int j = 0; j < selvars.size()+1; j++) {
//                x2[i][j] = X[i][j] * my * (1 - my);
//            }
//        }
//        der2 = cross_x_y(X, x2,  selvars.size()+1, maxrows,  selvars.size()+1);
//        RealMatrix deriv2 = new Array2DRowRealMatrix(der2);
//        DecompositionSolver solver = new SingularValueDecomposition(deriv2).getSolver();
//        RealVector deriv = new ArrayRealVector(der);
//        RealVector u2 = solver.solve(deriv);
//        double[] u = u2.toArray();
//
//        lambda = new double[ selvars.size()+1];
//        for (int i = 0; i <  selvars.size()+1; i++) {
//            lambda[i] = u[i] * der[i];
//
//            be[i] += u[i];
//            vars.get(numVars-1).betas.add(be[i]);
//        }
//        double[] p = new double[maxrows];
//        RealMatrix newdata = new Array2DRowRealMatrix(X);
//        double[] ykapelo = newdata.operate(be);
//
//        double[] expyhat = new double[maxrows];
//        for (int i = 0; i < maxrows; i++) {
//            expyhat[i] = Math.exp(-ykapelo[i]);
//            p[i] = 1 / (1 + expyhat[i]);
//            vars.get(numVars-1).scores.add(p[i]);
//        }
//        d2 = calcDevRes(p, target, expyhat);
//        double beta = 0.5;
//        double ta;
//        double[] B = new double[ selvars.size()+1];
//        double[] nextB;
//        Arrays.fill(B, 0.0);
//        double tol = 1e-09;
//        int i = 2;
//        double[] yminusp;
//        double[] resids;
//        double pavrg;
//        double[][] xmultW = new double[maxrows][selvars.size()+1];
//        double[] W = new double[maxrows];
//        Arrays.fill(W, 0.0);
//        for (; (d2 - d1 > tol) && (i < 100); ++i) {
//            vars.get(numVars-1).betas.clear();
//            vars.get(numVars-1).scores.clear();
//            d1 = d2;
//            yminusp = new double[maxrows];
//            for (int k = 0; k < maxrows; k++) {
//                yminusp[k] = target[k] - p[k];
//            }
//            sum = 0.0;
//            der = new double[ selvars.size()+1];
//            for (int k = 0; k <  selvars.size()+1; k++) {
//                for (int j = 0; j < maxrows; j++) {
//                    sum += X[j][k] * yminusp[j];
//                }
//                der[k] = sum;
//                sum = 0.0;
//            }
//            for(int k = 0;k<maxrows;k++){
//                W[k] = p[k]*(1-p[k]);
//            }
//            for(int k = 0;k<selvars.size()+1;k++){
//                for(int j = 0;j<maxrows;j++){
//                    xmultW[j][k] = X[j][k]*W[j];
//                }
//            }
//            der2 = cross_x_y(X,xmultW,selvars.size()+1,maxrows,selvars.size()+1);
//            deriv2 = new Array2DRowRealMatrix(der2);
//            solver = new SingularValueDecomposition(deriv2).getSolver();
//            deriv = new ArrayRealVector(der);
//            u2 = solver.solve(deriv);
//            u = u2.toArray();
//            lambda = new double[ selvars.size()+1];
//
//            for (int k = 0; k <  selvars.size()+1; k++) {
//                lambda[k] = u[k] * der[k];
//                be[k] += u[k];
//                vars.get(numVars-1).betas.add(be[k]);
//            }
//            ta = 1/beta;
//            nextB = new double[selvars.size()+1];
//            do{
//                vars.get(numVars-1).betas.clear();
//                vars.get(numVars-1).scores.clear();
//                ta = ta* beta;
//                for(int k = 0;k<selvars.size()+1;k++){
//                    nextB[k] = be[k]+ta*u[k];
//                }
//                Arrays.fill(ykapelo, 0.0);
//                for (int k = 0; k < maxrows; k++) {
//                    for (int j = 0; j <  selvars.size()+1; j++) {
//                        ykapelo[k] += X[k][j] * nextB[j];
//                    }
//                }
//                for (int k = 0; k < maxrows; k++) {
//                    expyhat[k] = Math.exp(-ykapelo[k]);
//                    p[k] = 1 / (1 + expyhat[k]);
//                    vars.get(numVars-1).scores.add(p[k]);
//                }
//                d2 = calcDevRes(p,target,expyhat);
//            }while (d2-d1>1e-4*ta*lambda[0] && ta>1e-09);
//            for(int k = 0;k<selvars.size()+1;k++){
//                be[k] = nextB[k];
//            }
//        }
//
//        retval.i = i;
//        retval.currDev =-2.0 * d2;
//        vars.get(numVars-1).betas.clear();
//        for( i = 0;i<selvars.size()+1;i++){
//            retval.be[i] = be[i];
//            vars.get(numVars-1).betas.add(be[i]);
//        }
//        resids = new double[numSamples];
//        for( i = 0;i<maxrows;i++){
//            retval.resids[i] = target[i] - p[i];
//        }
//        for( i = maxrows;i<numSamples;i++){
//            retval.resids[i] = target[i] - vars.get(numVars-1).mean;
//        }
//
//
//
//
//        return retval;
//    }
}
