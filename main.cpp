#include <iostream>
#include <fstream>
#include <cmath>
#include <bitset>
#include <cstdint>


// on testing
#define TESTCASES 100000
#define BIG_ERROR_PERCENT 10
#define SHORTTEST 10
#define EPSILON 0.0000000000000000001 //this one is kinda arbitrary

// on the lognum itself
#define BASE 2
#define INTBITS 4
#define FRACBITS 4
#define W_BITS (INTBITS + FRACBITS)

using namespace std;
static double logPrecision, minLogVal, maxLogVal, nearZeroRealVal, minRealVal, maxRealVal;

void setup();
void verifyLogs();
void fullRangeTestSetAddition();
void fullRangeTestSetMultiply();
void smallTestSet();
void deltaComparisonTestSet();
void percentErrorAnalysis(string goldpath, string testpath);
void MSE_Analysis(string goldpath, string testpath);

struct logString {
    string bits;
    double aslog;
};

string closestLog(double);
string int2bin(int a);
string sint2bin(int b);
void verify2bins();

int main() {
    setup();
    srand(545); // Do not change this value. I set it to 545. Please leave that.

    verifyLogs(); // works last time I checked
    verify2bins();// works last time I checked

    return 0;
}

void setup() {
    cout << "\nGiven " << W_BITS << " bits total (" << INTBITS << " are int, "
                << FRACBITS << " are frac)" <<" for log_base = " << BASE << "\n\n";
    // calculate log ranges
    logPrecision = pow(BASE,-FRACBITS);
    minLogVal = -1.0 * pow(BASE,INTBITS-1);
    maxLogVal = -1.0 * minLogVal - logPrecision;

    // corresponding real range
    maxRealVal = pow(BASE,maxLogVal);
    minRealVal = -1.0*pow(BASE,abs(minLogVal));
    nearZeroRealVal = pow(BASE,minLogVal);

    // the printed values (and stored) are a function of printing precision

    cout << "Log resolution is: " << logPrecision << endl;
    cout << "Highest positive log value before overflow is: " << maxLogVal << endl;
    cout << "Most negative log value before overflow is: " << minLogVal << endl;
    cout << "(shown with double floating point precision)\n" << endl;

    cout << "Range of real values representable is: (" << minRealVal << ", " << maxRealVal << ").\n";
    cout << "Smallest (abs) real value before 0 overflow is: " << nearZeroRealVal << endl;
    cout << "(shown with double floating point precision)\n" << endl;
}
void fullRangeTestSetAddition(){
   cout << "Generating REAL and LOG inputs across (1/2) input range for addition" << endl;
   cout << "\tThis guarantees that no overflow will occur" << endl;
   string file1name("real_addition_inputs.txt"),file2name("log_addition_inputs.txt"),
                    file3name("golden_real_addition_results.txt");
   ofstream realInputs(file1name);
   ofstream logInputs(file2name);
   ofstream goldOutputs(file3name);

   // we need to pick a + b to add that will not cause overflow
    double maxAddend = 1.0*maxRealVal/BASE;

    for (int i = 0; i < TESTCASES; ++i) {
        // get integer part
        double a = rand() % (int)maxAddend; double b = rand() % (int)maxAddend;
        // get fractional part
        a -= (double)rand()/RAND_MAX; b -= (double)rand()/RAND_MAX;
        // make sure we get some negatives in there
        if (i%2) {a *= -1.0;} if (i%3) {b *= -1.0;}

        realInputs << a << endl;
        realInputs << b << endl;
        goldOutputs << (a+b) << endl;

        string aLog = closestLog(a); string bLog = closestLog(b);
        logInputs << aLog << endl;
        logInputs << bLog << endl;
    }

    cout << "Three files created: " << file1name << " , " << file2name << " , and " << file3name << endl;
    realInputs.close();
    logInputs.close();
    goldOutputs.close();
}
void fullRangeTestSetMultiply() {
    cout << "Generating REAL and LOG inputs across sqrt(input range) for multiplication" << endl;
    cout << "\tThis guarantees that no overflow will occur" << endl;
    string file1name("real_multiplication_inputs.txt"),file2name("log_multiplication_inputs.txt"),
            file3name("golden_real_multiplication_results.txt");
    ofstream realInputs(file1name);
    ofstream logInputs(file2name);
    ofstream goldOutputs(file3name);

    // we need to pick a + b to add that will not cause overflow
    double maxFactor = sqrt(maxRealVal);

    for (int i = 0; i < TESTCASES; ++i) {
        // get integer part
        double a = rand() % (int)maxFactor; double b = rand() % (int)maxFactor;
        // get fractional part
        a -= (double)rand()/RAND_MAX; b -= (double)rand()/RAND_MAX;
        // make sure we get some negatives in there
        if (i%2) {a *= -1.0;} if (i%3) {b *= -1.0;}

        realInputs << a << endl;
        realInputs << b << endl;
        goldOutputs << (a*b) << endl;

        string aLog = closestLog(a); string bLog = closestLog(b);
        logInputs << aLog << endl;
        logInputs << bLog << endl;
    }

    cout << "Three files created: " << file1name << " , " << file2name << " , and " << file3name << endl;
    realInputs.close();
    logInputs.close();
    goldOutputs.close();
}
void smallTestSet() {
    cout << "Generating REAL and LOG inputs across small range for add/multiply" << endl;
    cout << "\tThese tests squeeze x into range -1 < x_real < +1" << endl;
    cout << "\tThis test set should be used to ensure sufficient precision with a small set of operands.\n";
    cout << "\tIf insufficient precision: add frac. bits OR decrease base." << endl;
    string file1name("real_small_inputs.txt"),file2name("log_small_inputs.txt"),
            file3name("golden_real_multiplication_results.txt"),file4name("golden_real_addition_results.txt");
    ofstream realInputs(file1name);
    ofstream logInputs(file2name);
    ofstream goldMultOutputs(file3name);
    ofstream goldAddOutputs(file4name);

    for (int i = 0; i < TESTCASES; ++i) {
        // DON'T get any integer part
        double a = 0.0; double b = 0.0;
        // get fractional part
        a += (double)rand()/RAND_MAX; b += (double)rand()/RAND_MAX;

        // make sure we get some negatives in there
        if (i%2) {a *= -1.0;} if (i%3) {b *= -1.0;}

        realInputs << a << endl;
        realInputs << b << endl;
        goldMultOutputs<< (a*b) << endl;
        goldAddOutputs << (a+b) << endl;

        string aLog = closestLog(a); string bLog = closestLog(b);
        logInputs << aLog << endl;
        logInputs << bLog << endl;
    }

    cout << "Four files created: " << file1name << " , " << file2name << " , "
                                    << file3name << " and " << file4name << endl;
    realInputs.close();
    logInputs.close();
    goldAddOutputs.close();
    goldMultOutputs.close();
}
void deltaComparisonTestSet() {
    cout << "Generating REAL and LOG inputs across (1/2) input range for DELTA testing" << endl;
    cout << "\tTesting done via log addition - guarantee that no overflow will occur" << endl;
    string file1name("real_delta_PLUS_inputs.txt"),file2name("log_delta_PLUS_inputs.txt"),
            file3name("golden_real_delta_PLUS_results.txt");
    string file4name("real_delta_MINUS_inputs.txt"),file5name("log_delta_MINUS_inputs.txt"),
            file6name("golden_real_delta_MINUS_results.txt");
    ofstream realPLUSInputs(file1name); ofstream realMINUSInputs(file4name);
    ofstream logPLUSInputs(file2name);  ofstream logMINUSInputs(file5name);
    ofstream goldPLUSOutputs(file3name);ofstream goldMINUSOutputs(file6name);

    // we need to pick a + b to add that will not cause overflow
    double maxAddend = 1.0*maxRealVal/BASE;

    for (int i = 0; i < TESTCASES; ++i) {
        double a,b;
        // 1. Pick a number |a| in the addition range - okay to change a by adding an integer, multiplying by -1, etc.
        a = (double)rand()/RAND_MAX; // add "+ rand() % (int)maxAddend" to move a across the range
        // 2. Calculate A_fixed (float here);
        double logA = log(a)/log(BASE);
        // 3. Calculate B_fixed (float here) by moving (up to) 1 away in log domain
        double eps = (double)rand()/RAND_MAX;
        double logB = logA + eps;
        // 4. Go back to the number |b|
        b = pow(BASE,logB);

        // Exercise a DeltaPlus test and a DeltaMinus using a+b, a-b
        double bplus(b),bminus(b);
        if (i % 2) {a *= -1.0;}

        if (a<0) {
            if (bplus>0) {bplus *= -1.0;}  // delta plus means b should be negative too
            if (bminus<0) {bminus *= -1.0;}// delta minus means b would be positive then
        } else { // a > 0 so make b > 0, bminus < 0
            if (bplus<0) {bplus *= -1.0;} // delta plus means b should be positive too
            if (bminus>0) {bminus *= -1.0;}//delta minus means b should be negative then
        }

        realPLUSInputs << a << endl;
        realMINUSInputs << a << endl;
        realPLUSInputs << bplus << endl;
        realMINUSInputs << bminus << endl;

        goldPLUSOutputs << (a+bplus) << endl;
        goldMINUSOutputs << (a+bminus) << endl;

        string aLog = closestLog(a); string bPlusLog = closestLog(bplus); string bMinusLog = closestLog(bminus);
        logPLUSInputs << aLog << endl;
        logMINUSInputs << aLog << endl;
        logPLUSInputs << bPlusLog << endl;
        logMINUSInputs << bMinusLog << endl;
    }

    cout << "Six files created: " << file1name << " , " << file2name << " , " << file3name <<
    " , " << file4name << " , " << file5name << " , and " << file6name << endl;
    realPLUSInputs.close();
    realMINUSInputs.close();
    logPLUSInputs.close();
    logMINUSInputs.close();
    goldPLUSOutputs.close();
    goldMINUSOutputs.close();
}

/*
 * YOUR responsibility to verify input files (same length, only doubles, etc)
 *      and of course that you're not comparing the test N's output to test M's golden
 * */
void percentErrorAnalysis(string goldpath, string testpath) {
    ifstream gold(goldpath), test(testpath);
    int testCount(0); // by the end, this should = TESTCASES,or you changed something
    int largeErrorCt(0);
    double numGold,numTest; // corresponding values
    double totalPercentError(0.0);
    while (!gold.eof() && !test.eof()) {
        ++testCount;
        gold >> numGold;
        test >> numTest;
        double samplePercentError = 100.0 * abs(numGold-numTest)/numGold;
        if (samplePercentError > BIG_ERROR_PERCENT) {++largeErrorCt;}
        totalPercentError += samplePercentError;
    }

    double avgPercentError = totalPercentError/testCount;
    cout << "In " << testCount << " cases, saw " << largeErrorCt <<
            " cases with significant error ( >" << BIG_ERROR_PERCENT << "%)\n";
    printf("Average Percent Error: %2.6f",avgPercentError);

}

/*
 * YOUR responsibility to verify input files (same length, only doubles, etc)
 *      and of course that you're not comparing the test N's output to test M's golden
 * */
void MSE_Analysis(string goldpath, string testpath) {
    ifstream gold(goldpath), test(testpath);
    int testCount(0); // by the end, this should = TESTCASES,or you changed something
    int largeErrorCt(0);
    double numGold,numTest; // corresponding values
    double diff(0.0), totalPercentError(0.0);
    double MSE(0.0);
    while (!gold.eof() && !test.eof()) {
        ++testCount;
        gold >> numGold;
        test >> numTest;
        diff = abs(numGold-numTest);
        double samplePercentError = 100.0 * diff/numGold;
        if (samplePercentError > BIG_ERROR_PERCENT) {++largeErrorCt;}
        MSE += diff*diff;
    }

    double RMSD = sqrt(MSE);

    cout << "In " << testCount << " cases, saw " << largeErrorCt <<
         " cases with significant error ( >" << BIG_ERROR_PERCENT << "%)\n";
    printf("MSE of %d ops: %.7f", MSE);
    printf("\nRMSD of %d ops: %.7f",RMSD);
}


/*
 * Code you wrote while you were tired
 * Are you *sure* it works...?
 * Also you didn't call Lily or Meg this weekend
 * */
string closestLog(double x) {
    // the 'true' logval... can be positive (|x| > 1) or negative ( -1 < x < 1)
    double logfloat = log(abs(x))/log(BASE);
    double optionHigh(0.0),optionLow(0.0);

    // Right away, need to get rid of the zero case. By now, logfloat is inf and this error will propogate
    if (x == 0) {
        return to_string(minLogVal);
    }

    if (logfloat >= 0) { // |x| > 1
        // cutoff for overshoot is okay: logfloat is within 1 resolution of next whole #
        if (abs(ceil(logfloat)-logfloat) < logPrecision) {
            // easy, now your two options are
            // the whole # L, or L - precision
            optionHigh = ceil(logfloat);
            optionLow = optionHigh - logPrecision;
        } else {
            // normal procedure
            // 1. Get a close approximation for integer part
            int logfixedint = floor(logfloat);

            // 2. Get the fraction bits
            double fracpart = logfloat - logfixedint; // positive bc logfixedint <= logfloat

            // 2a. Shift them over
            int lowfrac = floor(pow(BASE,FRACBITS) * fracpart);
            int highfrac = lowfrac + 1;

            // 3. A good guess
            optionLow = 1.0*logfixedint + pow(BASE,-FRACBITS)*lowfrac;
            optionHigh = 1.0*logfixedint + pow(BASE,-FRACBITS)*highfrac;
            //optionHigh = optionLow + logPrecision;
        }

    } else { // -1 < x < +1
        // cutoff for overshoot is better: logfloat w/n 1 resolution of next whole #
        if (abs(logfloat - ceil(logfloat)) < logPrecision) {
            // easy, now your two options are
            // the whole # L, or L - precision
            optionHigh = ceil(logfloat);
            optionLow = optionHigh - logPrecision;
        } else {
            // normal procedure
            // 1. Get a close approximation for integer part
            int logfixedint = floor(logfloat);

            // 2. Get the fraction bits
            double fracpart = logfloat - logfixedint; // positive bc logfixedint <= logfloat

            // 2a. Shift them over
            double logfrac = floor(pow(BASE,FRACBITS) * fracpart);

            // 3. Possible fraction parts
            int lowfrac = floor(logfrac);
            int highfrac = ceil(logfrac);

            // 4. Leading to possible guesses
            optionLow = 1.0*logfixedint + pow(BASE,-FRACBITS)*lowfrac;
            optionHigh = 1.0*logfixedint + pow(BASE,-FRACBITS)*highfrac; // optionLow + logPrecision;
        }
    }

    if (optionLow == -0) {optionLow = 0;}
    if (optionHigh == -0) {optionHigh = 0;}


    double LowtoReal = pow(BASE,optionLow); double HightoReal = pow(BASE,optionHigh);
    double absXreal = abs(x);

    if (absXreal < LowtoReal || absXreal > HightoReal) {
        cout << "Error type 1 occured on x = " << x << " with lowerLog = "
            << optionLow << " and higherLog = " << optionHigh << endl;
    }
    //if (abs(optionLow - optionHigh) > logPrecision + EPSILON) { // also works
    if (abs(optionLow + logPrecision - optionHigh) > EPSILON) { // error on x_real = 0 if not yet caught
        cout << "Error type 2 occured on x = " << x << " with lowerLog = "
             << optionLow << " and higherLog = " << optionHigh << endl;
    }

    // DON'T just pick the closer log value to float log(|x|)
    //      if (abs(optionLow-logfloat) < abs(optionHigh-logfloat)) {

    // DO pick the log value that gives the 2^(logval) closest to |x|
    if (abs(LowtoReal-x) < abs(HightoReal-x)) {
        // Used to show how we decided in console view
//        cout << "For x = " << x << " choose LOWER: " << optionLow << " which gives: " << pow(BASE,optionLow) << endl;
        return ("Option low but first convert to lognum bit format");
    }


// Used to show how we decided in console view
//    cout << "For x = " << x << " choose UPPER: " << optionHigh  << " which gives: " << pow(BASE,optionHigh) << endl;
    return ("Option high but first convert to lognum bit format");

}

/*
 * Prints console message if any of the log procedures failed
 * */
void verifyLogs() {
    for (int i = 0; i < TESTCASES; ++i) {
        int num = rand() % (int)(maxRealVal - 1*minRealVal);
        num -= (int)(maxRealVal);
        closestLog(num);
    }
}


/*
 * Turn the fractional part of the lognum into a bitvector, padded out to FRACBITS length
 * takes in an unsigned integer - there is NO way this should be interpretted as a negative
 * */
string int2bin(int a) {
    return bitset<FRACBITS>(a).to_string();
}

/*
 * Turn the integer part of the lognum into a bitvector, padded accordingly
 * takes in a signed integer. Need to check sign and deal with negative crap
 * */
string sint2bin(int b) {
    if (b >= 0) {
        return bitset<INTBITS>(b).to_string();
    }
    // ahh, now the fun case
    // bitset takes in an UNSIGNED number as the ctor, so
    // we need to do this manually. The python "bin()" function would be nice
    // alternatively we could use AP_INT<> type, but that's another big #include for not much


    // 1. Get the positive version
    int bPos = -1*b;

    // 2. Get binary representation
    auto x = bitset<INTBITS>(bPos);

    // 3. Get 1's complement
    x.flip();

    // 4. Get 2's complement by ripple adding 1...
    bool carry = x[0];  // , carry(1,x) = x
    x[0] = ~x[0];       // sum(1,x) = ~x
    for (int i = 1; i < INTBITS; ++i) {
        bool sum = x[i] ^ carry; // sum(x,y) = x ^ y, carry(x,y) = x&y
        carry = x[i] & carry; // propogate the carry if this bit is a 1
        x[i] = sum;
    }
    return x.to_string();
}

void verify2bins() {
    // Case 1, integer part.  Likely place for a screw up
        // Given INTBITS = x, smallest number is -2^(x-1)
        //                    largest number is  -1 + 2^(x-1)
        // Width of range is 2^x
    cout << "Begin signed integer testing. Hope this works\n";
    int biggestSigned = -1 + pow(2,INTBITS-1);
    int smallestSigned = -1*pow(2,INTBITS-1);
    int intRange = pow(2,INTBITS);
    int halfIntRange = intRange/2;

//    cout << "\nDecimal " << smallestSigned << " as bitvector : " << sint2bin(smallestSigned);
//    cout << "\nDecimal " << biggestSigned << " as bitvector : " << sint2bin(biggestSigned);
//    cout << "\nDecimal " << 0 << " as bitvector : " << sint2bin(0);
    for (int i = smallestSigned; i <= biggestSigned; ++i) {
//        int num = rand() % halfIntRange;
//        if(i%2) {num *= -1;}
        cout << "\nDecimal " << i << " as bitvector : " << sint2bin(i);
    }

    cout << "\n\nBegin unsigned integer testing. This should work\n";
    // Case 2, fractional part.  Unlikely to incur errors
        // Given FRACBITS - y, smallest number is 0
        //                     largest number is -1+2^y
        // width of range is 2^t
    int biggestU = -1 + pow(2,FRACBITS);
    int smallestU = 0;
    int fracRange = pow(2,FRACBITS);
//    cout << "\nDecimal " << biggestU << " as bitvector : " << sint2bin(biggestU);
//    cout << "\nDecimal " << smallestU << " as bitvector : " << sint2bin(smallestU);
    for (int i = smallestU; i <= biggestU; ++i) {
//        int num = rand() % fracRange;
        cout << "\nDecimal " << i << " as bitvector : " << int2bin(i);
    }
}