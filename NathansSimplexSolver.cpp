//from https://www.youtube.com/watch?v=_5T70cAXDJ0 Stopped at 2:08

#include<pybind11/pybind11.h>
#include<pybind11/stl.h>

#include<iostream> //cerr
#include<string> //string
#include<vector> //vector
#include<sstream> //stringstream
#include<limits> //numeric_limits
#include<chrono> //time_point
#include<utility> //pair
#include<map> //map
#include<cmath> //signbit

#include <unistd.h> //sleep

#define LOG_ERROR 0
#define LOG_WARN 1
#define LOG_INFO 2

//from https://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            log(LOG_ERROR,message); \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)


namespace py = pybind11;

using Scalar = double;

float DummyFunc(float arg1, float arg2) {
    return arg1 + arg2;
}

typedef std::chrono::high_resolution_clock clock_;
typedef std::chrono::duration<double, std::ratio<1> > second_;

class SimplexSolver {
protected:
    Scalar* tableauBody = NULL;
    size_t numVar = 0;
    size_t numBaseVar = 0;
    size_t numConstr = 0;
    std::string* varNames = NULL;
    std::map<std::string, size_t> varIDs;
    //int* basisVars = NULL; DEPRECATED SINCE WE'RE GOING TO SCAN THROUGH AT THE END TO CHECK FOR DEGENERACY.
    size_t tableauHeight = 0;
    size_t tableauWidth = 0;
    bool modelEngaged = false;
    std::chrono::time_point<clock_> tic_time;
    Scalar SCALAR_INF = 0;
    size_t SIZE_T_INF = 0;
    size_t bottomRowStart = 0;
    size_t bottomRowEnd = 0;
    size_t lastColStart = 0;
    size_t lastColEnd = 0;
    bool pointerOwnership = true;
    bool enableVariableNames = true;
    bool auxProblemSolved = false;
    bool maximizationProblem = true;

    size_t maxIter = 0;
    Scalar maxTime = 0;

    bool modelOptimized;

    bool enableLiveUpdates = false;
    size_t liveUpdateIter = 0;
    Scalar liveUpdateTime = 0;

    std::vector<std::pair<size_t, std::string>> logs;

    static const int UNBOUNDED_FLAG = 001;
    static const int INFEASIBLE_FLAG = 002;

    static const size_t CONSTR_EQ = 0;
    static const size_t CONSTR_LEQ = 1;
    static const size_t COSNTR_GEQ = 2;
    
public:
    SimplexSolver() {
        SCALAR_INF = std::numeric_limits<Scalar>::max();
        SIZE_T_INF = std::numeric_limits<size_t>::max();

        maxIter = SIZE_T_INF;
        maxTime = SCALAR_INF;
    }

    ~SimplexSolver() {
        DisengageModel();
    }

    void log(size_t level, const std::string& message) {
        logs.push_back(std::make_pair(level, message));
    }

    std::vector<std::pair<size_t,std::string>> GetLogs() {
        return logs;
    }

    void EngageModel(
        size_t& initNumVar,
        size_t& initNumConstr,
        std::vector<std::string>& initVarNames,
        std::vector<Scalar>& initA,
        std::vector<Scalar>& initB,
        std::vector<Scalar>& initC,
        std::vector<size_t>& constrTypes,
        bool initMaximization
        ) {
        if (modelEngaged) {
            DisengageModel();
        }

        ASSERT(initNumVar == initVarNames.size(),"ERROR! The number of variable names provided does not match the specified number of variables.");
        ASSERT(initA.size() == initNumVar * initNumConstr,"ERROR! The size of the specified A matrix does not match the number of variables and constraints specified.");
        ASSERT(initB.size() == initNumConstr,"ERROR! The size of the specified b vector does not match the specified number of constraints.");
        ASSERT(initC.size() == initNumVar + 1,"ERROR! The size of the specified c vector does not match the specified number of variables.");
        ASSERT(constrTypes.size() == initNumConstr,"ERROR! The size of the specified constraint type vector does not match the specified number of constraints.");

        maximizationProblem = initMaximization;

        //First, we need to deal with the different constraint types
        size_t numLEQ = 0;
        size_t numGEQ = 0;
        size_t numEQ = 0;
        for (size_t i = 0; i < initNumConstr; i++) {
            if (constrTypes[i] == CONSTR_EQ) {
                numLEQ++;
                numGEQ++;
                numEQ++;
            }
            else if (constrTypes[i] == CONSTR_LEQ) {
                numLEQ++;
            }
            else {
                numGEQ++;
            }
        }

        numVar = initNumVar + numLEQ + 2 * numGEQ;
        numBaseVar = initNumVar;
        numConstr = initNumConstr + numEQs;

        if (numGEQ > 0) {
            //This indicates that we'll need to sovle the auxilary problem first and then the normal problem
            //Solving the auxilary problem involves including the initial objective with it's function as a constraint.
            auxProblemSolved = false;
            numVar += 1;
            numConstr += 1;
        }
        else {
            //This indicates that each of the LEQ's slack vars form the initial basis. We don't need to solve the auxilary problem
            auxProblemSolved = true;
        }

        tableauHeight = numConstr + 1;
        tableauWidth = numVar + 1;
        tableauBody = new Scalar[tableauWidth * tableauHeight];


        enableVariableNames = true;
        varNames = new std::string[numVar];
        //First, load the original variables
        for (size_t i = 0; i < initNumVar; i++) {
            varNames[i] = initVarNames[i];
            varIDs[varNames[i]] = i;
        }
        //Now create the slack variables. Every constraint has a slack variable so they'll just be indexed in sequence starting after numBaseVar.
        std::string baseSlackVarName = "SLACK_VAR_";
        for (size_t i = 0; i < numConstr; i++) {
            varNames[numBaseVar + i] = baseSlackVarName + std::to_string(i);
        }
        //Now create the artificial variables. These will be indexed in sequence following the slack variables.
        std::string baseArtVarNAme = "ARTIF_VAR_";
        for (size_t i = 0; i < numGEQ; i++) {
            varNames[numBaseVar + numConstr + i] = baseArtVarNAme + std::to_string(i);
        }
        //If artificial variables are present, the last variable will the objective of the original problem.
        if (numGEQ > 0) {
            varNames[numBaseVar + numConstr + numGEQ] = "ORIG_OBJ_FUNC";
        }

        //Now populate the talbleau body.
        //First, the original A matrix.
        size_t numEQEncountered = 0;
        for (size_t i = 0; i < initA.size(); i++) {
            //First, find the 2-D indices of the original A matrix.
            size_t ii = i / numVar; //Constraint Index
            size_t jj = i - ii * numVar; //Variable Index

            //Now map those 2-D indices to the 2-D indices of the tableau.
            //jj remains the same since original variables come first in variable order
            //ii remains the same except it is repeated if a constraint is an equality constraint. 
            bool thisIsAnEQConstr = constrTypes[ii] == CONSTR_EQ;

            //Thus, we need to shift ii forward the number of equality constraints we already encountered.
            ii += numEQEncountered;

            size_t j = ii * tableauWidth + jj;
            tableauBody[j] = initA[i];

            if (thisIsAnEQConstr) {
                //Repeat this value on the next row if this is an EQConstr.
                ii++;
                j = ii * tableauWidth + jj;
                tableauBody[j] = initA[i];
            }

            //Only iterate numEQEncountered after we've totally finished a given row.
            if ((jj == initNumVar - 1) && (thisIsAnEQConstr)) {
                numEQEncountered++;
            }
        }

        //Now Fill in the Constant Column
        numEQEncountered = 0;
        for (size_t i = 0; i < numConstr; i++) {
            //First, find the index of the original. This is i.
            //Then, map that index onto the indices of the tableau.
            size_t j = (i + numEQEncountered) * tableauWidth
            tableauBody[j] = initB[i];

            //If this constr is an EQ constr, repeat it's value on the next tableau row.
            bool thisIsAnEQConstr = constrTypes[ii] == CONSTR_EQ;
            if (thisIsAnEQConstr) {
                j += tableauWidth;
                tableauBody[j] = initB[i];
            }
        }

        //Now fill in the slack variable values.
        //LEQ constraints have a 1 in the diagonal and 0s everywhere else
        TODO

        //GEQ constraints have a -1 in the diagonal and 0s everywhere else
        TODO

        //Now fill in the art. var. values.
        //LEQ constraints have zeros everywhere
        TODO

        //GEQ constraints have a 1 in the "diagonal" and zeros everywhere else
        TODO

        //Now fill in the objective function row.
        //This behaves differently if we need to solve the auxilary problem first or not.
        bottomRowStart = tableauWidth * numConstr;
        bottomRowEnd = bottomRowStart + numVar;

        lastColStart = tableauWidth - 1;
        lastColEnd = tableauHeight * tableauWidth;

        if (auxProblemSolved) {
            //The bottom row will be the original objective function.

            //The first entries will correspond to the original values. (But Negative)
            for (size_t i = 0; i < numBaseVar; i++) {
                tableauBody[i + bottomRowStart] = -initC[i];
            }
            //The obj value will be put in the appropriate place.
            tableauBody[bottomRowEnd-1] = initC[numBaseVar-1];

            //The rest of the obj function row will be zeros.
            for (size_t i = 0; i < numLEQ; i++) {
                tableauBody[i + bottomRowStart + numBaseVar] = 0;
            }
        }
        else {
            //The bottom-most row will be the auxilary objective function.

            //We want to minimize the art. vars. but since that objective function and associated minimiation is problematic, we'll solve for the equivalent objective coefs.
            //  For an example, see 10:10 of https://www.youtube.com/watch?v=-RtEzwfMqxk&list=PLg2tfDG3Ww4vyVtIvTUY2JaOZDbQStcsb&index=8
            //The first coefs will be the sum of the A matrix columns
            TODO

            //The Slack variable columns will be -1 if the associated constraint is a leq constraint, 0 otherwise.
            TODO

            //The art. var. columns will be 0
            TODO

            //The original objective variable will have a value of zero
            TODO

            //The constant column will be the sum of the original constant column
            TODO


            //The next row up will be the original objective function. 
            //The first columns correspond to the original variables. (Made negative since P - Coefs = Const)
            TODO

            //The constant column is the same.
            TODO

            //The original objective variable has a coef of 1.
            TODO

            //All other variables have a coef of 0.
            TODO
        }

        pointerOwnership = true;

        modelOptimized = false;
        modelEngaged = true;
    }

    void EngageModel(size_t initNumVar, size_t initNumConstr, Scalar* initTableauBody, std::string* initVarNames = NULL, bool transferOwnership = false) {
        if (modelEngaged) {
            DisengageModel();
        }

        tableauBody = initTableauBody;
        numVar = initNumVar;
        numConstr = initNumConstr;
        varNames = initVarNames;
        if (initVarNames == NULL) {
            enableVariableNames = false;
        }
        else {
            for (size_t i = 0; i < numVar; i++) {
                varIDs[varNames[i]] = i;
            }
        }

        tableauHeight = numConstr + 1;
        tableauWidth = numVar + 1;

        bottomRowStart = tableauWidth * numConstr;
        bottomRowEnd = bottomRowStart + numVar;

        lastColStart = tableauWidth - 1;
        lastColEnd = tableauHeight * tableauWidth;

        pointerOwnership = transferOwnership;

        modelEngaged = true;
        modelOptimized = false;
    }

    void DisengageModel() {
        if (modelEngaged) {
            if (pointerOwnership) {
                if (enableVariableNames) {
                    delete[] varNames;
                }
                delete[] tableauBody;
                if (basisVars != NULL) {
                    delete[] basisVars;
                }
                delete[] slackRelationships;
            }

            tableauBody = NULL;
            numVar = 0;
            numBaseVar = 0;
            numConstr = 0;
            varNames = NULL;
            basisVars = NULL;
            tableauHeight = 0;
            tableauWidth = 0;
            varIDs.clear();
            pointerOwnership = false;
            enableVariableNames = false;
            numSlackRelationships = 0;

            modelOptimized = false;
            modelEngaged = false;
        }
    }

    void setMaxIter(size_t newMaxIter) {
        maxIter = newMaxIter;
    }

    void setMaxTime(Scalar newMaxTime) {
        maxTime = newMaxTime;
    }

    std::string TableauToString() {
        return TableauToTerminalString();
    }

    std::string TableauToMarkdownString() {
        if (!enableVariableNames) {
            varNames = new std::string[numVar];
            for (size_t i = 0; i < numVar; i++) {
                varNames[i] = std::to_string(i);
            }
        }

        std::stringstream ss;
        ss << "$\\begin{array}{c|";
        for (size_t i = 0; i < numVar; i++) {
            ss << "c";
        }
        ss << "|c}\nBASIS & ";
        for (size_t i = 0; i < numVar; i++) {
            ss << varNames[i] << " & ";
        }
        ss << "CONSTANT \\\\\n\\hline ";
        for (size_t i = 0; i < numConstr; i++) {
            ss << varNames[basisVars[i]];
            for (size_t j = 0; j < tableauWidth; j++) {
                size_t ii = i * tableauWidth + j;

                ss << " & " << tableauBody[ii];
            }
            ss << "\\\\\n";
        }
        ss << "\\hline ";
        for (size_t j = 0; j < tableauWidth; j++) {
            size_t ii = numConstr * tableauWidth + j;

            ss << " & " << tableauBody[ii];
        }
        ss << "\\\\\n\\end{array}$";

        if (!enableVariableNames) {
            delete[] varNames;
        }

        return ss.str();
    }

    std::string TableauToCSVString() {
        if (!enableVariableNames) {
            varNames = new std::string[numVar];
            for (size_t i = 0; i < numVar; i++) {
                varNames[i] = std::to_string(i);
            }
        }

        std::stringstream ss;
        ss << "BASIS, |, ";
        for (size_t i = 0; i < numVar; i++) {
            ss << varNames[i] << ", ";
        }
        ss << "|, CONSTANT \n";

        for (size_t i = 0; i < numVar+2; i++) {
            ss << "-, ";
        }
        ss << "\n";

        for (size_t i = 0; i < numConstr; i++) {
            ss << varNames[basisVars[i]] << ", |";
            for (size_t j = 0; j < tableauWidth; j++) {
                size_t ii = i * tableauWidth + j;

                if (j == tableauWidth-2) {
                    ss << ", |";
                }
                ss << ", " << tableauBody[ii];
            }
            ss << "\n";
        }
        
        for (size_t i = 0; i < numVar+2; i++) {
            ss << "-, ";
        }
        ss << "\n , , |";

        for (size_t j = 0; j < tableauWidth; j++) {
            size_t ii = numConstr * tableauWidth + j;

            if (j == tableauWidth-2) {
                    ss << ", |";
                }

            ss << ", " << tableauBody[ii];
        }
        ss << "\n";

        if (!enableVariableNames) {
            delete[] varNames;
        }

        return ss.str();
    }

    std::string FormatString(const std::string& str, size_t width) {
        size_t numSpaces;
        if (width < str.length()) {
            numSpaces = 0; //Otherwise the will silently overflow since size_t can't be negative.
        }
        else {
            numSpaces = width - str.length();
        }
        size_t numFirstHalf = numSpaces / 2;
        size_t numSecondHalf = numSpaces - numFirstHalf;

        std::stringstream ss;
        for (size_t i = 0; i < numFirstHalf; i++) {
            ss << " ";
        }
        ss << str;
        for (size_t i = 0; i < numSecondHalf; i++) {
            ss << " ";
        }
        return ss.str();
    }

    std::string FormatScalar(const Scalar& num, size_t width, size_t precision = 6) {
        std::stringstream ss;
        ss.precision(precision);
        bool addNeg = false;
        if (std::signbit(num)) {
            if (num == 0) {
                ss << -num;
                addNeg = false;
            }
            else {
                ss << -num;
                addNeg = true;
            }
        }
        else{
            ss << num;
        }

        std::string str = FormatString(ss.str(),width);
        if (addNeg) {
            size_t startIndex = -1;
            for (size_t i = 0; i < str.length(); i++) {
                if (str.at(i) != ' ') {
                    startIndex = i;
                    break;
                }
            }

            str.replace(startIndex-1,1,"-");
        }
        return str;
    }

    std::string FormatSizeT(size_t num, size_t width) {
        std::stringstream ss;
        ss << num;
        std::string str = FormatString(ss.str(),width);

        return str;
    }

    std::string TableauToTerminalString() {
        if (!enableVariableNames) {
            varNames = new std::string[numVar];
            for (size_t i = 0; i < numVar; i++) {
                varNames[i] = std::to_string(i);
            }
        }
        
        size_t cellWidth = 10;
        for (size_t i = 0; i < numVar; i++) {
            if (varNames[i].length() > cellWidth) {
                cellWidth = varNames[i].length() + 2;
            }
        }
        for (size_t i = 0; i < tableauHeight * tableauWidth; i++) {
            std::stringstream ss;
            ss.precision(6);
            ss << tableauBody[i];
            std::string str = ss.str();
            if (str.length() + 2 > cellWidth) {
                cellWidth = str.length() + 2;
            }
        }

        const std::string vertBar = " | ";
        const std::string horiBar(cellWidth * (tableauWidth + 1) + 2 * vertBar.length(), '-');

        std::stringstream ss;

        ss << FormatString("BASIS",cellWidth) << vertBar;
        for (size_t i = 0; i < numVar; i++) {
            ss << FormatString(varNames[i],cellWidth);
        }
        ss << vertBar << "CONSTANT\n" << horiBar << "\n";

        for (size_t i = 0; i < numConstr; i++) {
            ss << FormatString(varNames[basisVars[i]],cellWidth) << vertBar;
            for (size_t j = tableauWidth * i; j < tableauWidth * (i+1) - 1; j++) {
                ss << FormatScalar(tableauBody[j],cellWidth);
            }
            ss << vertBar << FormatScalar(tableauBody[tableauWidth * (i+1) - 1],cellWidth) << "\n";
        }

        ss << horiBar << "\n";
        for (size_t i = 0; i < cellWidth; i++) {
            ss << " ";
        }
        ss << vertBar;

        for (size_t j = tableauWidth * (tableauHeight-1); j < tableauWidth * tableauHeight - 1; j++) {
            ss << FormatScalar(tableauBody[j],cellWidth);
        }
        ss << vertBar << FormatScalar(tableauBody[tableauWidth * tableauHeight - 1],cellWidth) << "\n";

        
        if (!enableVariableNames) {
            delete[] varNames;
        }

        return ss.str();
    }

    bool LessThanZero(const Scalar& val, const Scalar& tollerance = 1e-7) {
        return (val + tollerance) < 0.0;
    }

    bool GreaterThanZero(const Scalar& val, const Scalar& tollerance = 1e-7) {
        return (val - tollerance) > 0.0;
    }

    bool EqualsZero(const Scalar& val, const Scalar& tollerance = 1e-7) {
        if (val < 0) {
            return -val < tollerance;
        }
        else{
            return val < tollerance;
        }
    }

    void tic() {
        tic_time = clock_::now();
    }

    Scalar toc() {
        return std::chrono::duration_cast<second_>(clock_::now() - tic_time).count();
    }

    bool GetOptimalityStatus(bool maximize) {
        bool optimal = true;
        if (maximize) {
            for (size_t i = bottomRowStart; i < bottomRowEnd; i++) {
                if (LessThanZero(tableauBody[i])) {
                    optimal = false;
                    break;
                }
            }
        }
        else {
            for (size_t i = bottomRowStart; i < bottomRowEnd; i++) {
                if (GreaterThanZero(tableauBody[i])) {
                    optimal = false;
                    break;
                }
            }
        }
        return optimal;
    }

    size_t ComputePivotColumn(bool maximize) {
        size_t columnIndex = 0;
        if (maximize) {
            Scalar minBottomVal = tableauBody[bottomRowStart];
            size_t i = bottomRowStart + 1;
            for (size_t j = 1; j < numVar; i++, j++) {
                if (tableauBody[i] < minBottomVal) {
                    minBottomVal = tableauBody[i];
                    columnIndex = j;
                }
            }
        }
        else {
            Scalar maxBottomVal = tableauBody[bottomRowStart];
            size_t i = bottomRowStart + 1;
            for (size_t j = 1; j < numVar; i++, j++) {
                if (tableauBody[i] > maxBottomVal) {
                    maxBottomVal = tableauBody[i];
                    columnIndex = j;
                }
            }
        }

        return columnIndex;
    }

    size_t ComputePivotRow(size_t pivotCol) {
        size_t rowIndex = 0;

        size_t lastColIndex = lastColStart;
        size_t pivotColIndex = pivotCol;
        Scalar minComparisonVal;
        if (LessThanZero(tableauBody[pivotColIndex])) { // Since we're only interested in the smallest positive value. Note that we're only looking at the pivot value since the constant column will allways be positive.
            minComparisonVal = SCALAR_INF; 
        }
        else {
            minComparisonVal = tableauBody[lastColIndex] / tableauBody[pivotColIndex];
        }

        lastColIndex += tableauWidth;
        pivotColIndex += tableauWidth;
        for (size_t j = 1; j < numConstr; j++, lastColIndex += tableauWidth, pivotColIndex += tableauWidth) {
            if (EqualsZero(tableauBody[pivotColIndex])) {
                continue; //This will return an infinite result which, under no circumstances, is the best.
            }
            if (LessThanZero(tableauBody[pivotColIndex])) { // Since we're only interested in the smallest positive value. Note that we're only looking at the pivot value since the constant column will allways be positive.
                continue;
            }
            
            Scalar comparisonVal = tableauBody[lastColIndex] / tableauBody[pivotColIndex];

            if (comparisonVal < minComparisonVal) {
                rowIndex = j;
                minComparisonVal = comparisonVal;
            }
        }

        if (minComparisonVal == SCALAR_INF) {
            throw UNBOUNDED_FLAG;
        }

        return rowIndex;
    }

    void PerformPivot(const size_t& pivotCol, const size_t& pivotRow) {
        //Step 1: Divide the pivot row by the pivot value
        Scalar pivotValue = tableauBody[tableauWidth * pivotRow + pivotCol];

        size_t startVal = tableauWidth * pivotRow;
        size_t endVal = startVal + tableauWidth;

        for (size_t i = startVal; i < endVal; i++) {
            tableauBody[i] /= pivotValue;
        }

        //Step 2: Row reduce all the other rows so that the pivot value (which is now 1) is the only non-zero value in the pivot column
        for (size_t row_i = 0; row_i < tableauHeight; row_i++) {
            if (row_i == pivotRow) {
                continue;
            }

            Scalar rowMultiplier = tableauBody[row_i * tableauWidth + pivotCol];
            for (size_t col_i = 0; col_i < tableauWidth; col_i++) {
                tableauBody[row_i * tableauWidth + col_i] -= rowMultiplier * tableauBody[pivotRow * tableauWidth + col_i];
            }
        }

        //Step 3: Swap in the appropriate basis valriable
        basisVars[pivotRow] = pivotCol;
    }

    Scalar GetVariableValue(std::string varName) {
        ASSERT(enableVariableNames,"Error! You cannot get the variable values by name when variable names are not enabled. Use GetVariableValue(size_t index) instead.");

        if (!modelOptimized) {
            log(LOG_WARN,"Attempting to access the value of \"" + varName + "\" within a model that has not yet been optimized.");
        }
        size_t varIndex;
        try {
            varIndex = varIDs.at(varName);
        }
        catch (const std::out_of_range& err) {
            ASSERT(false,"Error! Attempting to access the value of a variable that is not registered as part of this model.");
        }

        size_t basisIndex;
        bool indexFound = false;
        for (size_t i = 0; i < numConstr; i++) {
            if (basisVars[i] == varIndex) {
                basisIndex = i;
                indexFound = true;
                break;
            }
        }

        if (!indexFound) {
            return 0.0;
        }
        else {
            return tableauBody[(basisIndex+1) * tableauWidth - 1];
        }
    }

    Scalar GetVariableValue(size_t varIndex) {
        size_t basisIndex;
        bool indexFound = false;
        for (size_t i = 0; i < numConstr; i++) {
            if (basisVars[i] == varIndex) {
                basisIndex = i;
                indexFound = true;
                break;
            }
        }

        if (!indexFound) {
            return 0.0;
        }
        else {
            return tableauBody[(basisIndex+1) * tableauWidth - 1];
        }
    }

    Scalar GetObjectiveValue() {
        return tableauBody[tableauWidth * tableauHeight - 1];
    }

    void SetLiveUpdateSettings(size_t newLiveUpdateIter, Scalar newLiveUpdateTime) {
        enableLiveUpdates = true;
        liveUpdateIter = newLiveUpdateIter;
        liveUpdateTime = newLiveUpdateTime;
    }

    void PrintUpdate(size_t itr, Scalar time) {
        std::cout << "Itr: " << FormatSizeT(itr,10) << " Time: " << FormatScalar(time,15) << " Objective Value: " << GetObjectiveValue() << '\n';
    }

    size_t SolveInitProblemWithSolver() {
        //NOTES ON STUFF TO FIX:
        // For each Constraint, add a artif. var. If the coef. if neg. put -1 as the coef for the art. var.
        // If the coef. is positive put 1 as the art. var. coef.
        // Then, for the obj. func. values sum up the non-art. var. coefs and multiply by negatvie one as shown here: https://www.youtube.com/watch?v=-RtEzwfMqxk&list=PLg2tfDG3Ww4vyVtIvTUY2JaOZDbQStcsb&index=8

        if (basisVars == NULL) {
            ASSERT(pointerOwnership,"Error! It's unclear how to deal with a NULL basis that we don't have ownership of.");
            basisVars = new int[numConstr];
        }
        otherSolver.setMaxIter(maxIter);
        otherSolver.setMaxTime(maxTime);

        //Add an artificial variable to each constraint
        size_t numArtificialVars = numConstr;

        size_t numAuxVars = numArtificialVars;
        size_t otherNumVars = numVar + numAuxVars;
        Scalar* otherTableau = new Scalar[(otherNumVars + 1) * (numConstr + 1)];

        size_t otherTableauHeight = tableauHeight;
        size_t otherTableauWidth = otherNumVars + 1;

        //Load the other Tableau Coefs
        Scalar sumConsts = 0;
        for (size_t i = 0; i < numConstr; i++) {
            //Constant Value:
            size_t myIndex = (i+1) * tableauWidth - 1;
            size_t otherIndex = (i+1) * otherTableauWidth - 1;
            Scalar constVal = tableauBody[myIndex];
            otherTableau[otherIndex] = constVal;
            sumConsts += constVal;
            
            //Regular Var Coefs
            for (size_t j = 0; j < numVar; j++) {
                size_t myIndex = i * tableauWidth + j;
                size_t otherIndex = i * otherTableauWidth + j;
                otherTableau[otherIndex] = tableauBody[myIndex];
            }
        }
        //Artificial Variable Coefs
        for (size_t i = 0; i < numArtificialVars; i++) {
            size_t varI = numVar + i;
            for (size_t j = 0; j < numConstr; j++) {
                if (i == j) {
                    otherTableau[j * otherTableauWidth + varI] = -1;
                }
                else {
                    otherTableau[j * otherTableauWidth + varI] = 0;
                }
            }
        }

        //Now go through and specify new objective row (to minimize all auxilary variables)
        for (size_t i = 0; i < numVar; i++) {
            otherTableau[otherTableauWidth * numConstr + i] = 0;
        }
        for (size_t i = numVar; i < otherNumVars; i++) {
            otherTableau[otherTableauWidth * numConstr + i] = 1;
        }
        otherTableau[otherTableauHeight * otherTableauWidth - 1] = 0;

        //Now specify a initial basis.
        //Since the Aux valriables can be whatever, We'll specify that the original variables must be zero.
        //  i.e. The initial basis must consist of just the artificial variables.
        for (size_t i = 0; i < numArtificialVars; i++) {
            basisVars[i] = numVar + i;
        }

        //Engage Model
        std::string* otherVarNames = NULL;
        if (enableVariableNames_InitProblem && enableVariableNames) {
            otherVarNames = new std::string[otherNumVars];
            for (size_t i = 0; i < numVar; i++) {
                otherVarNames[i] = varNames[i];
            }
            for (size_t i = 0; i < numArtificialVars; i++) {
                otherVarNames[i + numVar] = "ARTIF_VAR[" + std::to_string(i) + "]";
            }
        }
        otherSolver.EngageModel(otherNumVars,numConstr,otherTableau,numSlackRelationships,slackRelationships,basisVars,otherVarNames,false);
        if (enableLiveUpdates) {
            otherSolver.SetLiveUpdateSettings(liveUpdateIter, liveUpdateTime);
        }

        //Perform Forced Pivot
        std::cout << otherSolver.TableauToTerminalString() << "\n";
        size_t forcePivotCol = 0;
        size_t forcePivotRow = otherSolver.ComputePivotRow(forcePivotCol);
        otherSolver.PerformPivot(forcePivotCol,forcePivotRow);

        //Sove regularly
        size_t otherIterations = otherSolver.Solve();

        if (enableVariableNames_InitProblem && enableVariableNames) {
            delete[] otherVarNames;
        }

        //Copy over logs
        std::vector<std::pair<size_t,std::string>> otherLogs = otherSolver.GetLogs();
        for (size_t i = 0; i < otherLogs.size(); i++) {
            logs.push_back(std::make_pair(otherLogs[i].first,"AUXILARY PROBLEM: " + otherLogs[i].second));
        }

        //Check Feasibility
        //In order to be feasible, the feasibility variable and all artificial variables must be absent from the basis.
        bool feasible = true;
        for (size_t i = 0; i < numConstr; i++) {
            if (basisVars[i] >= numVar) {
                feasible = false;
                break;
            }
        }

        delete[] otherTableau;

        if (!feasible) {
            throw INFEASIBLE_FLAG;
        }

        //The otherSolver acted directly on this solver's basisVars object.
        //  That, and the fact that all the non artificial/feasibility variables lie after the original variables
        //  makes it so that the computed basis is valid for the original model as well. There is no need to
        //  copy it over.

        

        return otherIterations;
    }

    size_t SolveCurrentSetup(bool maximize, size_t& iterationNum, size_t& liveUpdateIterCounter, size_t& liveUpdateTimeCounter) {
        //Make sure model is engaged before proceeding
        ASSERT(modelEngaged,"ERROR! Solver cannot execute if there is no model engaged.");


        bool solutionIsOptimal = GetOptimalityStatus(maximize);

        
        size_t exitCode = 0;

        std::cout << TableauToTerminalString() << "\n";
        while (!solutionIsOptimal) {
            if (iterationNum > maxIter) {
                exitCode = 1;
                log(LOG_INFO,"Maximum number of Iterations Exceeded.");
                break;
            }
            iterationNum++;
            Scalar currentTime = toc();
            if (currentTime > maxTime) {
                exitCode = 2;
                log(LOG_INFO,"Maximum solve time exceeded.");
                break;
            }

            if (enableLiveUpdates) {
                liveUpdateIterCounter++;
                if (liveUpdateIterCounter == liveUpdateIter) {
                    PrintUpdate(iterationNum,currentTime);
                    liveUpdateIterCounter = 0;
                }
                if (currentTime > liveUpdateTime * liveUpdateTimeCounter) {
                    PrintUpdate(iterationNum,currentTime);
                    liveUpdateTimeCounter++;
                }
            }

            size_t pivotCol = ComputePivotColumn(maximize);
            size_t pivotRow;
            try {
                pivotRow = ComputePivotRow(pivotCol);
            }
            catch (const int& e) {
                switch (e) {
                    case UNBOUNDED_FLAG:
                        log(LOG_ERROR,"The Model was proven to be unbounded.");
                    default:
                        throw e;
                }
            }
            PerformPivot(pivotCol,pivotRow);

            solutionIsOptimal = GetOptimalityStatus(maximize);

            std::cout << TableauToTerminalString() << "\n";
        }

        if (exitCode == 0) {
            log(LOG_INFO,"Optimal solution found.");
            modelOptimized = true;
        }
        log(LOG_INFO,std::to_string(iterationNum) + " iterations");
        log(LOG_INFO,std::to_string(toc()) + " seconds");
        log(LOG_INFO,"Best Found Objective Function Value: " + std::to_string(GetObjectiveValue()));
    }

    size_t Solve() {
        size_t iterationNum = 0;
        size_t liveUpdateIterCounter = 0;
        size_t liveUpdateTimeCounter = 0;

        if (!auxProblemSolved) {
            SolveCurrentSetup(false,iterationNum,liveUpdateIterCounter,liveUpdateTimeCounter);

            //Assert that the auxilary problem was solved to optimality. Otherwise abort.
            TODO

            //Chop off the Auxilary Objective Function Row and all artificial variable / objective function columns
            //  In the Tableau
            TODO

            //  And in the variable name list
            TODO

            //Update all the tableau size constants (numVar, tableaWidth, etc.)
            TODO
        }

        //Now that the auxilary problem is sovled, the tableau is ready to be solved normally.
        SolveCurrentSetup(maximizationProblem,iterationNum,liveUpdateIterCounter,liveUpdateTimeCounter);
    }
};



PYBIND11_MODULE(nathans_Simplex_solver_py, handle) {
    handle.doc() = "The C++ implementation of Nathan's Simplex Solver.\n To operate, instantiate a SimplexSolver object by passing the following as arguments: numVar, numConstr, varNames, the A matrix (flattened using A.flatten() or A.flatten('C')), the b vector, the vector.\nThen specify any neccecary solver options.\nThen execute the solving operation using the Solve() method.";
    handle.def("DummyFunc", &DummyFunc);

    py::class_<SimplexSolver>(handle, "SimplexSolver")
        .def(py::init<>())
        .def("EngageModel", static_cast<void (SimplexSolver::*)(size_t&, size_t&,std::vector<std::string>&,std::vector<std::pair<size_t,size_t>>&,std::vector<Scalar>&,std::vector<Scalar>&,std::vector<Scalar>&)>(&SimplexSolver::EngageModel))
        .def("Solve", &SimplexSolver::Solve)
        .def("TableauToString", &SimplexSolver::TableauToString)
        .def("setMaxIter",&SimplexSolver::setMaxIter)
        .def("setMaxTime",&SimplexSolver::setMaxTime)
        .def("GetLogs",&SimplexSolver::GetLogs)
        .def("GetVariableValue", static_cast<Scalar (SimplexSolver::*)(std::string)>(&SimplexSolver::GetVariableValue))
        .def("GetVariableValue", static_cast<Scalar (SimplexSolver::*)(size_t)>(&SimplexSolver::GetVariableValue))
        .def("GetObjectiveValue",&SimplexSolver::GetObjectiveValue)
        .def("SetLiveUpdateSettings",&SimplexSolver::SetLiveUpdateSettings);
}