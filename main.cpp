#include <iostream>
#include <fstream>
#include <cmath>
#include <bitset>
#include "ConversionEngine.cpp"


// on testing
#define TESTCASES 100000
#define BIG_ERROR_PERCENT 1000
#define SHORTTEST 10
#define EPSILON 0.0000000000000000001 //this one is kinda arbitrary
#define SEEDVAL 545 // Do not change this. This val is 545. Do not change it

// on the lognum itself
#define BASE 2
#define INTBITS 4
#define FRACBITS 4
#define WBITS (INTBITS + FRACBITS)

using namespace std;
static double logPrecision, minLogVal, maxLogVal, nearZeroRealVal, minRealVal, maxRealVal;

// should be first function called in main
void setup();

// builds your test set
void fullRangeTestSetAddition();
void fullRangeTestSetMultiply();
void smallTestSet();
void deltaComparisonTestSet();

// analyzes your results (need real NOT LOG results)
void percentErrorAnalysis(string goldpath, string testpath);
void MSE_Analysis(string goldpath, string testpath);


// conversion from log<-->real
// required to construct input cases that have 0 conversion error
double convertToDouble(lognum L) {
    double d = pow(BASE,L.getLogval());
    if(L.getSignBit()) {return d;}
    return (-1.0 * d);
}

int convertToInt(lognum L) {
    double d = pow(BASE,L.getLogval());
    if (L.getSignBit()) {
        return d;
    }
    return (int)(-1*d);
}

// better quantizing - good to go
string logQuantize(double);
// unsigned bitvector of length FRACBITS
string int2bin(int a);
// 2's c bitvector of length INTBITS
string sint2bin(int b);
// 2's c bitvector of length WBITS
string sint3bin(int b);

// Run us once each time you change the configuration
// console output will let you know if any issues
void verifyLogs();
void verify2bins();
void verifyElog();
void verifyQuant();

// convention, "1" for abs(x) > 0, "0" otherwise
string sgn(double a) {
    string s = a > 0 ? "1" : "0";
    return s;
}

// Functional but deprecated
string closestLog(double); // Griffin will want closestLog, since he needs bit vectors for VHDL
double bestLog(double); // Matt will want bestLog, since he needs numeric types for C++

// unused for now
void InputFileWrite(ofstream &RealInputs, ofstream &LogInputs, double x, double y);

int main() {
    setup();          // Always run setup first, this resets for any changes you made.

/*
  // Run us for each setup
    verifyLogs();  // goes through best log, which is the "thinking" logic...
    verify2bins(); // you need this to get correct bit vector's (2's complement fixed point)
    verifyElog();
    verifyQuant();

*/

/*
    double maxAddend = 1.0*maxRealVal/BASE;
    for (int i = 0; i < 4; ++i) {
        double x = rand() % (int) maxAddend; double y = rand() % (int)maxAddend;
        x -= (double)rand()/RAND_MAX; y -= (double)rand()/RAND_MAX;
        if (i % 2) {x *= -1.0;} if (i % 3) {y *= -1.0;}
        cout << "x: " << x << "\t y: " << y << "\n";
        cout << "Logs: " << logQuantize(x) << "\t" << logQuantize(y) << "\n\n";
    }
*/

/*
 *    // we need to pick a + b to add that will not cause overflow
    double maxAddend = 1.0*maxRealVal/BASE;

    for (int i = 0; i < TESTCASES; ++i) {
        // get integer part
        double a = rand() % (int)maxAddend; double b = rand() % (int)maxAddend;
        // get fractional part
        a -= (double)rand()/RAND_MAX; b -= (double)rand()/RAND_MAX;
        // make sure we get some negatives in there
        if (i%2) {a *= -1.0;} if (i%3) {b *= -1.0;}
 */

    //fullRangeTestSetAddition();

    // This is an int, not a double, so that that the to_string gives us "BASE_2_" not "BASE_2.0000_"
        // you can switch it to double if using a non integer base...
    int ourbase = (2);
    cout << "Base " << ourbase << " results:\n";
    string ending = "_BASE_" + to_string(ourbase) + "_Intb_" + to_string(INTBITS) + "_Fracb_" + to_string(FRACBITS) + ".txt";

    string goldR("../golden_real_addition_results");
    string testR("../real_addition_test_results");
    string gold(goldR + ending), test(testR + ending);
    cout << "File 1: " << gold << endl;
    cout << "File 1: " << test << endl;

    percentErrorAnalysis(gold,test);
    cout << endl;
    MSE_Analysis(gold,test);

    return 0;
}

void setup() {
    cout << "\nGiven " << WBITS << " bits total (" << INTBITS << " are int, "
                << FRACBITS << " are frac)" <<" for log_base = " << BASE << "\n\n";
    // calculate log ranges
    logPrecision = pow(BASE,-FRACBITS);
    minLogVal = -1.0 * pow(BASE,INTBITS-1);
    maxLogVal = -1.0 * minLogVal - logPrecision;

    // corresponding real range
    maxRealVal = pow(BASE,maxLogVal);
    minRealVal = -1.0*pow(BASE,maxLogVal);
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
   srand(SEEDVAL); // Do not change this value. I set it to 545. Please leave that.
   cout << "Generating REAL and LOG inputs across (1/2) input range for addition" << endl;
   cout << "\tThis guarantees that no overflow will occur" << endl;
   string file1name("../real_addition_inputs"),file2name("../log_addition_inputs"),
                    file3name("../golden_real_addition_results");
   string ending = "_BASE_" + to_string(BASE) + "_Intb_" + to_string(INTBITS) + "_Fracb_" + to_string(FRACBITS) + ".txt";
   file1name += ending; file2name += ending; file3name += ending;
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

        // this finds the closest representable value and replaces each input with that
        // this way we only use inputs that have a perfect conversion to lognum format
        lognum aL = toLogNum(a);
        lognum bL = toLogNum(b);
        a = convertToDouble(aL);
        b = convertToDouble(bL);

        realInputs << a << endl;
        realInputs << b << endl;
        goldOutputs << (a+b) << endl;

        string aLog = logQuantize(a); string bLog = logQuantize(b);
        logInputs << aLog << endl;
        logInputs << sgn(a) << endl;
        logInputs << bLog << endl;
        logInputs << sgn(b) << endl;
    }

    cout << "Three files created: " << file1name << " , " << file2name << " , and " << file3name << endl;
    realInputs.close();
    logInputs.close();
    goldOutputs.close();
}
void fullRangeTestSetMultiply() {
    srand(SEEDVAL); // Do not change this value. I set it to 545. Please leave that.
    cout << "Generating REAL and LOG inputs across sqrt(input range) for multiplication" << endl;
    cout << "\tThis guarantees that no overflow will occur" << endl;
    string file1name("../real_multiplication_inputs"),file2name("../log_multiplication_inputs"),
            file3name("../golden_real_multiplication_results");
    string ending = "_BASE_" + to_string(BASE) + "_Intb_" + to_string(INTBITS) + "_Fracb_" + to_string(FRACBITS) + ".txt";
    file1name += ending; file2name += ending; file3name += ending;
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


        // this finds the closest representable value and replaces each input with that
        // this way we only use inputs that have a perfect conversion to lognum format
        lognum aL = toLogNum(a);
        lognum bL = toLogNum(b);
        a = convertToDouble(aL);
        b = convertToDouble(bL);

        realInputs << a << endl;
        realInputs << b << endl;
        goldOutputs << (a*b) << endl;

        string aLog = logQuantize(a); string bLog = logQuantize(b);
        logInputs << aLog << endl;
        logInputs << sgn(a) << endl;
        logInputs << bLog << endl;
        logInputs << sgn(b) << endl;
    }

    cout << "Three files created: " << file1name << " , " << file2name << " , and " << file3name << endl;
    realInputs.close();
    logInputs.close();
    goldOutputs.close();
}
void smallTestSet() {
    srand(SEEDVAL); // Do not change this value. I set it to 545. Please leave that.
    cout << "Generating REAL and LOG inputs across small range for add/multiply" << endl;
    cout << "\tThese tests squeeze x into range -1 < x_real < +1" << endl;
    cout << "\tThis test set should be used to ensure sufficient precision with a small set of operands.\n";
    cout << "\tIf insufficient precision: add frac. bits OR decrease base." << endl;
    string file1name("../real_small_inputs"),file2name("../log_small_inputs"),
            file3name("../golden_real_multiplication_results"),file4name("../golden_real_addition_results");
    string ending = "_BASE_" + to_string(BASE) + "_Intb_" + to_string(INTBITS) + "_Fracb_" + to_string(FRACBITS) + ".txt";
    file1name += ending; file2name += ending; file3name += ending; file4name += ending;
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


        // this finds the closest representable value and replaces each input with that
        // this way we only use inputs that have a perfect conversion to lognum format
        lognum aL = toLogNum(a);
        lognum bL = toLogNum(b);
        a = convertToDouble(aL);
        b = convertToDouble(bL);

        realInputs << a << endl;
        realInputs << b << endl;
        goldMultOutputs<< (a*b) << endl;
        goldAddOutputs << (a+b) << endl;

        string aLog = logQuantize(a); string bLog = logQuantize(b);
        logInputs << aLog << endl;
        logInputs << sgn(a) << endl;
        logInputs << bLog << endl;
        logInputs << sgn(b) << endl;
    }

    cout << "Four files created: " << file1name << " , " << file2name << " , "
                                    << file3name << " and " << file4name << endl;
    realInputs.close();
    logInputs.close();
    goldAddOutputs.close();
    goldMultOutputs.close();
}

// This function is not necessarily compatible with the 3/1/21 update
void deltaComparisonTestSet() {
    srand(SEEDVAL); // Do not change this value. I set it to 545. Please leave that.
    cout << "Generating REAL and LOG inputs across (1/2) input range for DELTA testing" << endl;
    cout << "\tTesting done via log addition - guarantee that no overflow will occur" << endl;
    string file1name("../real_delta_PLUS_inputs"),file2name("../log_delta_PLUS_inputs"),
            file3name("../golden_real_delta_PLUS_results");
    string file4name("../real_delta_MINUS_inputs"),file5name("../log_delta_MINUS_inputs"),
            file6name("../golden_real_delta_MINUS_results");
    string ending = "_BASE_" + to_string(BASE) + "_Intb_" + to_string(INTBITS) + "_Fracb_" + to_string(FRACBITS) + ".txt";
    file1name += ending; file2name += ending; file3name += ending;
    file4name += ending; file5name += ending; file6name += ending;

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

        string aLog = logQuantize(a);
        string bPlusLog = logQuantize(bplus);
        string bMinusLog = logQuantize(bminus);

        logPLUSInputs << aLog << endl; logPLUSInputs << sgn(a) << endl;
        logMINUSInputs << aLog << endl; logMINUSInputs << sgn(a) << endl;
        logPLUSInputs << bPlusLog << endl; logPLUSInputs << sgn(bplus) << endl;
        logMINUSInputs << bMinusLog << endl; logMINUSInputs << sgn(bminus) << endl;
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
    double numGold(0.0),numTest(0.0); // corresponding values
    double totalPercentError(0.0);
    while (!gold.eof() && !test.eof()) {
        if (testCount > (TESTCASES - 5)){
            // stop here
            cout.flush();
        }
        ++testCount;
        gold >> numGold;
        test >> numTest;
        double samplePercentError = 100.0 * abs((numGold-numTest)/numGold);
        if (samplePercentError > BIG_ERROR_PERCENT) {
            ++largeErrorCt;
        }
        totalPercentError += samplePercentError;
    }

    double avgPercentError = totalPercentError/testCount;
    cout << "In " << testCount << " cases, saw " << largeErrorCt <<
            " cases with significant error ( >" << BIG_ERROR_PERCENT << "%)\n";
    printf("Average Percent Error: %2.6f",avgPercentError);
    cout << "\n\n";

}

/*
 * Do not use this with delta function outputs. Those tests are biased and will
 * not be assessed accurately with MSE.  Use %error one.
 *
 * YOUR responsibility to verify input files (same length, only doubles, etc)
 *      and of course that you're not comparing the test N's output to test M's golden
 * */
void MSE_Analysis(string goldpath, string testpath) {
    ifstream gold(goldpath), test(testpath);
    int testCount(0); // by the end, this should = TESTCASES,or you changed something
    int largeErrorCt(0);
    double numGold(0.0),numTest(0.0); // corresponding values
    double diff(0.0), totalPercentError(0.0);
    double MSE(0.0);
    while (!gold.eof() && !test.eof()) {
        ++testCount;
        gold >> numGold;
        test >> numTest;
        diff = abs(numGold-numTest);
        double samplePercentError = 100.0 * abs(diff/numGold);
        if (samplePercentError > BIG_ERROR_PERCENT) {++largeErrorCt;}
        MSE += diff*diff;
    }

    double RMSD = sqrt(MSE);

    cout << "In " << testCount << " cases, saw " << largeErrorCt <<
         " cases with significant error ( >" << BIG_ERROR_PERCENT << "%)\n";
    printf("MSE of %d ops: %.7f", testCount,MSE);
    printf("\nRMSD of %d ops: %.7f",testCount,RMSD);
    printf("\n\tside note: if you got infinite percent error that just "
           "means a result =  0 test case came up. Generate new tests\n\n");
}

/*
 * I fully tested this function when the output for x_fixed was a double
 * I really don't think there is any issue now that I changed the output
 * to be a string.  The function returns at different points, but
 * inputs should still take the same paths they did before.
 *
 * That method lives on as "double bestLog(double)" just in case
 * if you need to find the procedure, go through bestLog, it will be easier
 * and it had better push the 5 cases through the same way:
 *      x = 1 case - dealt with up front
 *     |x| > 1
 *          - tricky
 *          - normal
 *     |x| < 1
 *          - tricky
 *          - normal
 *
 * */
string closestLog(double x) {
    string sgn = x > 0 ? "1": "0";
    double absXreal = abs(x);
    // the 'true' logval... can be positive (|x| > 1) or negative ( -1 < x < 1)
    double logfloat = log(abs(x))/log(BASE);
    double optionHigh(0.0),optionLow(0.0);

    // These will become bit vectors that we concatenate
    string intstr,fracstr,totalstr;

    // Flag 1: potential for overshoot w/ |a| > 1
    // Flag 2: potential for overshoot with |a| < 1
    bool flag1(false),flag2(false);

    // Right away, need to get rid of the zero case. By now, logfloat is inf and this error will propogate
    if (x == 0) { // should look like 10...0_0..0
        intstr = sint2bin((int)round(minLogVal));  // minLogVal was really an int stored as a double
        fracstr = int2bin(0);
        totalstr = intstr + fracstr;
        return totalstr;
    }

    // If we can't represent the number accurately, saturate to get as close as we can
    if (logfloat > maxLogVal) { // abs(x) is outside of our range (too big)
        intstr = sint2bin(maxLogVal);
        for (int i = 0; i < FRACBITS; ++i) {
            fracstr += "1";
        }
        totalstr = intstr + fracstr;
        return totalstr;
    } else if (logfloat < minLogVal) { //abs(x) requires greater precision than we have
        intstr = sint2bin(minLogVal);
        fracstr = int2bin(0);
        totalstr = intstr + fracstr;
        return  totalstr;
    }

    // need these variable names the same for all remaining paths
    int logfixedint, lowfrac,highfrac,thefrac;
    double fracpart;

    if (logfloat >= 0) { // |x| > 1
        // cutoff for overshoot is okay: logfloat is within 1 resolution of next whole #
        if (abs(ceil(logfloat)-logfloat) < logPrecision) {
            // easy, now your two options are
            // the whole # L, or L - precision   ex. logfloat  = 2.99999
            optionHigh = ceil(logfloat);       //ex. optionHI  = 3
            optionLow = optionHigh - logPrecision;// optionLOW = 3 - logprecision
            flag1 = true;
        } else {
            // normal procedure
            // 1. Get a close approximation for integer part
            logfixedint = floor(logfloat);

            // 2. Get the fraction bits
            fracpart = logfloat - logfixedint; // positive bc logfixedint <= logfloat

            // 2a. Shift them over
            lowfrac = floor(pow(BASE,FRACBITS) * fracpart);
            highfrac = lowfrac + 1;

            // 3. A good guess
            optionLow = 1.0*logfixedint + pow(BASE,-FRACBITS)*lowfrac;
            optionHigh = 1.0*logfixedint + pow(BASE,-FRACBITS)*highfrac;
            //optionHigh = optionLow + logPrecision; // equivalent

            if (optionLow == optionHigh) {
                // Indicates that the double logfrac was exactly an integer...
                // Unlikely, but possible.
                // If you see this error more than once in a blue moon, something is wrong
                cout << "Error type 0 occured on x = " << x << " with lowerLog = "
                     << optionLow << " and higherLog = " << optionHigh << endl;
            }
        }

    } else { // -1 < x < +1
        // cutoff for overshoot is better: logfloat w/n 1 resolution of next whole #
        if (abs(logfloat - ceil(logfloat)) < logPrecision) {
            // easy, now your two options are
            // the whole # L, or L - precision             ex. logfloat  = -2.01
            optionHigh = ceil(logfloat);                // ex. optionHI  = -2
            optionLow = optionHigh - logPrecision;      // ex. optionLO = - 2 - log precision
            flag2 = true;
        } else {
            // normal procedure
            // 1. Get a close approximation for integer part
            logfixedint = floor(logfloat);

            // 2. Get the fraction bits
            fracpart = logfloat - logfixedint; // positive bc logfixedint <= logfloat

            // 2a. Shift them over
            double logfrac = (pow(BASE,FRACBITS) * fracpart);

            // 3. Possible fraction parts
            lowfrac = floor(logfrac);
            highfrac = ceil(logfrac);

            // 4. Leading to possible guesses
            optionLow = 1.0*logfixedint + pow(BASE,-FRACBITS)*lowfrac;
            optionHigh = 1.0*logfixedint + pow(BASE,-FRACBITS)*highfrac; // optionLow + logPrecision;

            if (optionLow == optionHigh) {
                // Indicates that the double logfrac was exactly an integer...
                // Unlikely, but possible.
                // If you see this error more than once in a blue moon, something is wrong
                cout << "Error type 0 occured on x = " << x << " with lowerLog = "
                     << optionLow << " and higherLog = " << optionHigh << endl;
            }
        }
    }

    // Doubles support +0 and -0.  We don't really want -0
    if (optionLow == -0) {optionLow = 0;}
    if (optionHigh == -0) {optionHigh = 0;}


    double LowtoReal = pow(BASE,optionLow); double HightoReal = pow(BASE,optionHigh);

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

    if (flag1 || flag2) {
        if (flag1) {
            // Take upper branch if 2^5.00 gives a better approximation than 2^4.99
            if (abs(HightoReal - absXreal) <= abs(LowtoReal - absXreal)) {
                logfixedint = (int)round(optionHigh);// this should already be an int, round to take care of float behavior
                thefrac = 0; highfrac = thefrac;
            } else {
                // this branch indicates frac int should be maxed out
                logfixedint = (int)round(optionHigh) - 1;
                thefrac = (int)round(pow(BASE, FRACBITS)) - 1; lowfrac = thefrac;
            }
        } else {
            // Flag 2
            // Take upper branch if 2^-5.00 gives a better approximation than 2^-5.01
            if(abs(HightoReal - absXreal) <= abs(LowtoReal-absXreal)) {
                logfixedint = (int)round(optionHigh); // this should already be an int, rount to take care of float behavior
                thefrac = 0; highfrac = thefrac;
            } else {
                logfixedint = (int)round(optionHigh) - 1;
                thefrac = (int)round(pow(BASE, FRACBITS)) - 1; lowfrac = thefrac;
            }

        }
        intstr = sint2bin(logfixedint);
        fracstr = int2bin(thefrac);
        totalstr = intstr + fracstr;
        return totalstr;
    }


    if (abs(LowtoReal-absXreal) < abs(HightoReal-absXreal)) {
        // Used to show how we decided in console view
//        cout << "For x = " << x << " choose LOWER: " << optionLow << " which gives: " << pow(BASE,optionLow) << endl;

        intstr = sint2bin(logfixedint);
        fracstr = int2bin(lowfrac);
        totalstr = intstr + fracstr;
        return totalstr;
    }


// Used to show how we decided in console view
//    cout << "For x = " << x << " choose UPPER: " << optionHigh  << " which gives: " << pow(BASE,optionHigh) << endl;
    intstr = sint2bin(logfixedint);
    fracstr = int2bin(highfrac);
    totalstr = intstr + fracstr;
    return totalstr;
}

/*
 * Prints console message if any of the log procedures failed
 * */
void verifyLogs() {
    cout << "Testing bestLog finder (closestLog should do the same)" << endl;
    cout << "If no errors occurred, ";
    for (int i = 0; i < TESTCASES; ++i) {
        int num = rand() % (int)(maxRealVal - 1*minRealVal);
        num -= (int)(maxRealVal);
        bestLog(num);
    }
    cout << "this sentence is uninterrupted.\n" << endl;
}


/*
 * Turn the fractional part of the lognum into a bitvector, padded out to FRACBITS length
 * takes in an unsigned integer - there is NO way this should be interpreted as a negative
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
    cout << "Begin signed integer testing. Hope this works";
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

    cout << "\n\nBegin unsigned integer testing. This should work";
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

    cout << "\n\nCompleted binary 2's complement fixed point conversion testing\n\n";
}

/*
 *  *
 *    1    x = 1 case - dealt with up front
 *    2_   |x| > 1
 *     2a    - tricky
 *     2b    - normal
 *    3_   |x| < 1
 *     3a    - tricky
 *     3b    - normal
 * */
double bestLog(double x) {
    // the 'true' logval... can be positive (|x| > 1) or negative ( -1 < x < 1)
    double logfloat = log(abs(x))/log(BASE);
    double optionHigh(0.0),optionLow(0.0);
    // Right away, need to get rid of the zero case. By now, logfloat is inf and this error will propogate
    if (x == 0) {
        // CASE 1 ----------------------------------------------
        return minLogVal;
    }

    if (logfloat > maxLogVal) { // abs(x) very big, this is the best we can do
        return maxLogVal;
    } else if (logfloat < minLogVal) { // abs(x) very small, this is the best we can do
        return minLogVal;
    }

     if (logfloat >= 0) { // |x| > 1
        // cutoff for overshoot is okay: logfloat is within 1 resolution of next whole #
        if (abs(ceil(logfloat)-logfloat) < logPrecision) {
        // CASE 2a ----------------------------------------------
            // easy, now your two options are
            // the whole # L, or L - precision
            optionHigh = ceil(logfloat);
            optionLow = optionHigh - logPrecision;
        } else {
        // CASE 2b ----------------------------------------------
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

            if (optionLow == optionHigh) {
                // Indicates that the double logfrac was exactly an integer...
                // Unlikely, but possible.
                // If you see this error more than once in a blue moon, something is wrong
                cout << "Error type 0 occured on x = " << x << " with lowerLog = "
                     << optionLow << " and higherLog = " << optionHigh << endl;
            }
        }
    } else { // -1 < x < +1
        // cutoff for overshoot is better: logfloat w/n 1 resolution of next whole #
        if (abs(logfloat - ceil(logfloat)) < logPrecision) {
        // CASE 3a ----------------------------------------------
            // easy, now your two options are
            // the whole # L, or L - precision
            optionHigh = ceil(logfloat);
            optionLow = optionHigh - logPrecision;
        } else {
        // CASE 3b ----------------------------------------------
            // normal procedure
            // 1. Get a close approximation for integer part
            int logfixedint = floor(logfloat);
            // 2. Get the fraction bits
            double fracpart = logfloat - logfixedint; // positive bc logfixedint <= logfloat
            // 2a. Shift them over
            double logfrac = (pow(BASE,FRACBITS) * fracpart);
            // 3. Possible fraction parts
            int lowfrac = floor(logfrac);
            int highfrac = ceil(logfrac);
            // 4. Leading to possible guesses
            optionLow = 1.0*logfixedint + pow(BASE,-FRACBITS)*lowfrac;
            optionHigh = 1.0*logfixedint + pow(BASE,-FRACBITS)*highfrac; // optionLow + logPrecision;

            if (optionLow == optionHigh) {
                // Indicates that the double logfrac was exactly an integer...
                // Unlikely, but possible.
                // If you see this error more than once in a blue moon, something is wrong
                cout << "Error type 0 occured on x = " << x << " with lowerLog = "
                     << optionLow << " and higherLog = " << optionHigh << endl;
            }
        }
    }

    // Get rid of negative 0
    if (optionLow == -0) {optionLow = 0;}
    if (optionHigh == -0) {optionHigh = 0;}

    // At this point, all paths have been executed
    // We are only detecting errors and then selecting the correct higher/lower
    // We are not changing higher/lower anymore. Those are fixed by this point


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
    if (abs(LowtoReal-absXreal) < abs(HightoReal-absXreal)) {
        // Used to show how we decided in console view
//        cout << "For x = " << x << " choose LOWER: " << optionLow << " which gives: " << pow(BASE,optionLow) << endl;
        return (optionLow);
    }
// Used to show how we decided in console view
//    cout << "For x = " << x << " choose UPPER: " << optionHigh  << " which gives: " << pow(BASE,optionHigh) << endl;
    return (optionHigh);

}

/*
 * Cleaner quantization method, as suggested by Dr. Beerel 10/5
 * Much easier than the others to understand...
 *
 *
 * */
string logQuantize(double x) {
    double absXreal = abs(x);
    double logfloat = log(absXreal)/log(BASE);
    double logshift = logfloat * pow(BASE,FRACBITS);
    double shiftedminLogVal = minLogVal * pow(BASE,FRACBITS);
    string lognum;

    // saturation:
    if (absXreal < 1 && logshift < shiftedminLogVal) {
       // Since abs|x| ~= 0 and log|x| << minimum logval, saturate and return
        logshift = shiftedminLogVal;
        lognum = sint3bin(shiftedminLogVal);
        return lognum;
    }

    if (x == 0) {
        // Chose ZERO bound lognum: = minLogVal
        lognum = "1";
        for (int i = 0; i < (INTBITS + FRACBITS - 1); ++i) {
            lognum = lognum + "0";
        }
        return lognum;
    }

    int upper = (int)ceil(logshift);
    double upperboundlog = 1.0 * upper * pow(BASE,-FRACBITS);
    double upperboundreal = pow(BASE,upperboundlog);

    int lower = (int)floor(logshift);
    double lowerboundlog = 1.0 * lower * pow(BASE,-FRACBITS);
    double lowerboundreal = pow(BASE,lowerboundlog);

    string lowerlognum = sint3bin(lower);
    string upperlognum = sint3bin(upper);

    if (abs(upperboundreal-absXreal) < abs(lowerboundreal-absXreal)) {
        // upper bound better
        lognum = upperlognum;
        return lognum;
    }

    lognum = lowerlognum;
    return lognum;
}

string sint3bin(int b) {
    if (b >= 0) {
        return bitset<WBITS>(b).to_string();
    }
    // ahh, now the fun case
    // bitset takes in an UNSIGNED number as the ctor, so
    // we need to do this manually. The python "bin()" function would be nice
    // alternatively we could use AP_INT<> type, but that's a big #include


    // 1. Get the positive version
    int bPos = -1*b;

    // 2. Get binary representation
    auto x = bitset<WBITS>(bPos);

    // 3. Get 1's complement
    x.flip();

    // 4. Get 2's complement by ripple adding 1...
    bool carry = x[0];  // , carry(1,x) = x
    x[0] = ~x[0];       // sum(1,x) = ~x

    // starting at LSB + 1, move through
    for (int i = 1; i < WBITS ; ++i) {
        bool sum = x[i] ^ carry; // sum(x,y) = x ^ y, carry(x,y) = x&y
        carry = x[i] & carry; // propogate the carry if this bit is a 1
        x[i] = sum;
    }
    return x.to_string();
}

// Run once with each config to test correctness
void verifyElog() {
    logQuantize(0); cout << "\n\n";
    string s1,s2;
    for (int i = 0; i < SHORTTEST; ++i) {
        double j = double(rand())/RAND_MAX;
        if (i%2) {j*=-1;}
        logQuantize(j);

        j+=(rand() % (int)maxRealVal);
        if (i%4) {j*=-1;}
        logQuantize(j);
    }
}

// Run once with each config to test correctness
void verifyQuant() {
    for (int i = 0; i < TESTCASES; ++i) {
        double j = double(rand())/RAND_MAX;
        if (i%2) {j*=-1;}
        if (logQuantize(j) != closestLog(j)) {
            cout << "Error - closetLog produced: " << closestLog(j) << "\n\n";
        }

        j+=(rand() % (int)maxRealVal);
        if (i%4) {j*=-1;}
        if (logQuantize(j) != closestLog(j)) {
            cout << "Error - closetLog produced: " << closestLog(j) << "\n\n";
        }
    }
}

// This could be used, but isn't
void InputFileWrite(ofstream &realInputs, ofstream &logInputs, double a, double b) {
    realInputs << a << endl;
    realInputs << b << endl;

    string aLog = logQuantize(a); string bLog = logQuantize(b);
    logInputs << aLog << endl;
    logInputs << sgn(a) << endl;
    logInputs << bLog << endl;
    logInputs << sgn(b) << endl;
}