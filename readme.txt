Overview of Shared Log System's "main.cpp"
	This file exists to give generate tests for log domain operations.
	The main method() controls what functions are executed, namely:
		1 - Sets up the log-precision for chosen base and bit width
		2 - Verifies quantization and fixed point conversion are successful
		3 - Generate test inputs in real domain and log domains
		    and golden outputs in real domain. Stored as .txt files
		4 - Analysis of generated outpputs with golden outputs

How to use:
	1. Set parameters "BASE, "INTBITS" and "FRACBITS"

	2. In the main function, call these functions in this order:
	   main() {
		   setup();
		   verifyLogs();
		   verify2bins();
		   verifyElog();
	   } 
	   
	   Compile & run, then check console output to confirm setup is successful.
	
	3. In the main function, keep only the call to setup and add calls
	   for the desired tests. ex: addition and delta tests:
	   main() {
	   	setup();
	   	fullRangeTestSetAddition();
	   	deltaComparisonTest();
	   } 

	   Compile & run.

	4. Check your working directory for the generated test input/gold output files.
	   Use these with your "working engine" to create test outputs.
	   The test output files MUST contain 1 (decimal) result per line,
	   and should contain the same number of lines as the gold output file.  If not,
	   congrats, you've found a bug!

	5. Noting the names and paths to your test outputs, run selected analysis in
	   main function by using appropriate function. ex:
	   main() {
	   	percentErrorAnalysis("golden1.txt","mytest1.txt");
		MSE_analysis("golden2.txt","mytest2.txt");
	   }

	   Compile and run. Outputs will be displayed to the console, but not saved.
	   
Things to change & how to change them:
	#define MYVAR X
		These can be used to change the base, precision, number of test cases
		and tolerable error cutoffs.  Do not change the seed value.

	How things get printed (to say, your input.txt file)
		Go to the function which generates that test.  This is one of:
			fullRangeTestSetAddition(), fullRangeTestSetMultiply(), smallTestSet()
			or deltaComparisonTest().

I		In that function, switch the lines that look like:
		"realInputs << a << endl;" to change real input.txt
		"logInputs << aLog << endl;" to change loginput.txt

		If you change one of these lines, you should probably change all of the ones
		that write to the same file. (If you change realinputs << a << endl, then
		you should probably change realinputs << b << end)

		NOTE: if you change this for one function, you'll have to change it in the
		other three test generating functions for consistency.  If we move beyond
		four test functions, I'll find a better way to deal with this.

	What range of numbers are being tested
		Probably don't do this.  The test functions select the widest range of 
		numbers that will not cause overflow errors.

		Each test function has a variable (something like "maxAddend", for addition)
		that can be set to determine the input range.

	Sign of numbers being tested
		Erase any line that looks like "if (i%2) {a*= -1.0;}"

	Test only integer numbers
		Erase any line that looks like "a -= (double) rand()/RAND_MAX"

	Where the .txt files are stored
		In each test function, file1name
		Just like with how things get printed, changes to one function
		will not carry over to the other functions, so you may need to
		repeat the change for each function you need.
