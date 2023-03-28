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

template <typename T>
struct BasicTableau {
public:
    T* tableauBody = NULL;
    size_t numVar = 0;
    size_t numConstr = 0;
    
    size_t tableauHeight = 0;
    size_t tableauWidth = 0;
    bool pointerOwnership = true;

    struct Iterator {
    private:
        size_t index;
        T* tableauBody;
        size_t iterDistance;
        size_t end;

    public:
        Iterator(T* bodyRef, size_t initIndex, size_t initDist, size_t initEnd) {
            index = initIndex;
            tableauBody = bodyRef;
            iterDistance = initDist;
            end = initEnd;
        }

        T& operator*() {
            return tableauBody[index];
        }

        Iterator& operator++() {
            index += iterDistance;
            return *this;
        }

        friend bool operator== (const Iterator& a, const Iterator& b) {
            return a.index == b.index;
        }
        friend bool operator!= (const Iterator& a, const Iterator& b) {
            return a.index != b.index;
        }

        bool isEnd() {
            return index < end;
        }
    };

    BasicTableau(size_t initNumVar, size_t initNumConstr) {
        numVar = initNumVar;
        numConstr = initNumConstr;
        
        tableauHeight = numConstr + 1;
        tableauWidth = numVar + 1;

        tableauBody = new T[tableauHeight * tableauWidth];
        pointerOwnership = true;
    }

    BasicTableau(size_t initNumVar, size_t initNumConstr, T* initTableauBody) {
        tableauBody = initTableauBody;
        numVar = initNumVar;
        numConstr = initNumConstr;
        tableauHeight = numConstr + 1;
        tableauWidth = numVar + 1;
        pointerOwnership = false;
    }

    BasicTableau() {
        tableauBody = NULL;
        numVar = 0;
        numConstr = 0;
        tableauHeight = 0;
        tableauWidth = 0;
        pointerOwnership = false;
    }

    ~BasicTableau() {
        if (pointerOwnership) {
            delete[] tableauBody;
        }
    }

    T& operator[](size_t row, size_t col) {
        //FIXME: For production code, comment out these lines.
        ASSERT(row < tableauHeight,"Out of bounds error!");
        ASSERT(col < tableauWidth,"Out of bounds error!");

        return tableauBody[row * tableauWidth + col];
    }

    Iterator ColumnIterator(size_t col) {
        return Iterator(tableauBody, col, tableauWidth, tableauHeight * tableauWidth);
    }
    Iterator RowIterator(size_t row) {
        return Iterator(tableauBody, row * tableauWidth, 1, (row + 1) * tableauWidth);
    }

    T* getRawPtr() {
        return tableauBody;
    }
};

struct Tableau: public BasicTableau<Scalar> {
public:
    std::string* varNames = NULL;
    std::map<std::string, size_t> varIDs;
    bool enableVariableNames = true;

    Tableau(size_t initNumVar, size_t initNumConstr, std::vector<std::string> initVarNames = {}): BasicTableau<Scalar>(initNumVar, initNumConstr) {
        if (initVarNames.size() > 0) {
            enableVariableNames = true;
            varNames = new std::string[initNumVar];
            for (size_t i = 0; i < initNumVar; i++) {
                varNames[i] = initVarNames[i];
                varIDs[varNames[i]] = i;
            }
        }
        else {
            enableVariableNames = false;
            varNames = NULL;
        }

    }

    Tableau(size_t initNumVar, size_t initNumConstr, Scalar* initTableauBody): BasicTableau<Scalar>(initNumVar, initNumConstr, initTableauBody) {
        varNames = NULL;
        enableVariableNames = false;
    }

    Tableau() {
        varNames = NULL;
        enableVariableNames = false;
    }

    ~Tableau() {
        if (enableVariableNames) {
            varIDs.clear();
            delete[] varNames;
        }
    }

    void AddVarNames(std::vector<std::string> initVarNames) {
        if (enableVariableNames) {
            delete[] varNames;
            varIDs.clear();
        }
        enableVariableNames = true;
        varNames = new std::string[numVar];
        for (size_t i = 0; i < numVar; i++) {
            varNames[i] = initVarNames[i];
            varIDs[varNames[i]] = i;
        }
    }

    void AddVarNames(std::string* initVarNames) {
        if (enableVariableNames) {
            delete[] varNames;
            varIDs.clear();
        }
        enableVariableNames = true;
        varNames = initVarNames;
        for (size_t i = 0; i < numVar; i++) {
            varIDs[varNames[i]] = i;
        }
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
}

class SimplexSolver {
protected:
    Tableau tableau;
    BasicTableau<FormatingRelationship> initialTableauInfo;
    BasicTableau<FormatingRelationship> sfTableauInfo;
    
    size_t numBaseVar = 0;

    bool modelEngaged = false;
    std::chrono::time_point<clock_> tic_time;
    Scalar SCALAR_INF = 0;
    size_t SIZE_T_INF = 0;
    
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

    Scalar* shiftValues;
    bool* invertStatus;
    int* augmented_OldToNew_Index;
    int* augmented_NewToOld_Index;
    int* duplicatedConstr_OldToNew_Index;
    int* duplicatedConstr_NewToOld_Index;
    

    
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
        std::vector<Scalar> initA, //This is intentionally passed by value since we change some of it's value when converting to standard form.
        std::vector<Scalar> initB, //This is intentionally passed by value since we change some of it's value when converting to standard form.
        std::vector<Scalar> initC, //This is intentionally passed by value since we change some of it's value when converting to standard form.
        std::vector<size_t> constrTypes, //This is intentionally passed by value since we change some of it's value when converting to standard form.
        std::vector<std::pair<Scalar,Scalar>> initBounds, //This is intentionally passed by value since we change some of it's value when converting to standard form.
        std::vector<std::pair<bool,bool>>& initBoundActivations,
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

        //Re-write the problem into a stardard form talbeau

        //1: Compute Tableau Size
        //This is involve adding a number of additional variables and constraints.
        //So to start, I'll calculate the number of additional variables and constraints that I need.
        numBaseVar = initNumVar;

        //1.1: Now deal with the variable bounds
        shiftValues = new Scalar[initNumVar];
        invertStatus = new bool[initNumVar];
        augmented_OldToNew_Index = new int[initNumVar];
        size_t numAugmentedVariables = 0;
        size_t numActiveUpperBounds = 0;

        for (size_t i = 0; i < initNumVar; i++) {
            shiftValues[i] = 0;
            invertStatus[i] = false;
            augmented_OldToNew_Index[i] = -1;

            Scalar lowerBound = initBounds[i].first;
            Scalar upperBound = initBounds[i].second;

            bool hasLowerBound = initBoundActivations[i].first;
            bool hasUpperBound = initBoundActivations[i].second;

            if (hasLowerBound) {
                //We'll shift the definition of this value so that it's lower bound is zero.
                //So we don't have to add any constraints or variables.
                /*
                Example:
                    -5*x1 + 2*x2 <= 2
                    2 <= x1 <= 3

                    Becomes
                    -5*x1' + 2*x2 <= -8
                    0 <= x1' <= 1

                    Here, the shift value would be -2 since x1' = x1 - 2

                    Note that we'll also need to multiply this whole equation by -1 in order to maintain standard form.
                */
                if (lowerBound != 0) {
                    for (size_t j = 0; j < initNumConstr; j++) {
                        initB[j] += lowerBound * initA[j * initNumVar + i];
                    }

                    initC[initNumVar] += lowerBound * initC[i];
                    
                    shiftValues[i] = -lowerBound;
                    initBounds[i].first = 0;
                }

                if (hasUpperBound) {
                    //We'll have to add a constraint to specify this upper bound.
                    numActiveUpperBounds++;
                    initBounds[i].second -= lowerBound; //This must still be postiive so it's still properly formated to be put as a constraint later.
                }
                else {
                    //Standard form assumes no upper bound, so we don't need to do anything here.
                }
            }
            else {
                if (hasUpperBound) {
                    //We'll simulate a lower bound by multiplying by negative 1 and then shifting the new lower bound to zero.
                    //So we don't have to add any constraints or variables.

                    /*
                    Example:
                    3*x1 + 2*x2 <= 1
                    x <= 3

                    Becomes
                    -3*x1' + 2*x2 <= -8
                    x1' >= 0

                    Here, the invertStatus would be true and the shift value would be 3 since x1' = -x1 + 3

                    Note that we don't need to adda constraint here since x1' is already in standard form.
                    Also note that we'll need to multiply this whole equation by -1 in order to maintain standard form.
                    */
                    initC[initNumVar] -= upperBound * initC[i];
                    initC[i] *= -1;
                    for (size_t j = 0; j < initNumConstr; j++) {
                        initB[j] -= upperBound * initA[j * initNumVar + i];
                        initA[j * initNumVar + i] *= -1;
                    }
                    shiftValues[i] = upperBound;
                    invertStatus[i] = true;
                    initBounds[i].first = 0;
                    initBoundActivations[i].first = true;
                    initBoundActivations[i].second = false;
                }
                else {
                    //There are no bounds to shift. So we'll have to create a new variable with the same but opposite coefficients to simulate the negative behavior of this variable.
                    //But since both variables will have bounds at zero, we don't need to add any constraints.
                    augmented_OldToNew_Index[i] = initNumVar + numAugmentedVariables;
                    numAugmentedVariables++;
                }
            }
        }

        //1.1.1: Some of the variable shifts may have caused some b vector values to become negative. If this is the case, we must multiply the entire constraint by negative one and change the inequality constraitn type.
        for (size_t j = 0; j < initNumConstr; j++) {
            if (initB[j] < 0) {
                for (size_t i = 0; i < initNumVar; i++) {
                    initA[j * initNumVar + i] *= -1;
                }
                initB[j] *= -1

                if (constrTypes[j] == CONSTR_LEQ) {
                    constrTypes[j] = COSNTR_GEQ;
                }
                else if (constrTypes[j] == COSNTR_GEQ) {
                    constrTypes[j] = CONSTR_LEQ;
                }
            }
        }


        //1.2: Now deal with the different constraint types in the original A matrix
        duplicatedConstr_OldToNew_Index = new int*[initNumConstr];
        size_t numLEQ = 0;
        size_t numGEQ = 0;
        size_t numEQ = 0;
        for (size_t i = 0; i < initNumConstr; i++) {
            duplicatedConstr_OldToNew_Index[i] = -1;
            if (constrTypes[i] == CONSTR_EQ) {
                duplicatedConstr_OldToNew_Index[i] = initNumConstr + numEQ;
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

        //Each LEQ has a slack var
        //Each GEQ has a surplus var and an artificial var
        //Each augmented variable must be added
        //Each active upper bound will have a slack var.
        size_t numSlackSurplus = initNumConstr + numEQ + numActiveUpperBounds;
        size_t numArtificialVars = numGEQ;

        size_t numAdditionalVars = numAugmentedVariables + numSlackSurplus + numArtificialVars;
        size_t numAdditionalConstrs = numEQ + numActiveUpperBounds;
        
        if (numGEQ > 0) {
            //This indicates that we'll need to sovle the auxilary problem first and then the normal problem
            //Solving the auxilary problem involves including the initial objective with it's function as a constraint.
            auxProblemSolved = false;
            numAdditionalVars++;
            numAdditionalConstrs++;
        }
        else {
            //This indicates that each of the LEQ's slack vars form the initial basis. We don't need to solve the auxilary problem
            auxProblemSolved = true;
        }

        size_t numVar = initNumVar + numAdditionalVars;
        
        size_t numConstr = initNumConstr + numAdditionalConstrs;

        tableau = Tableau(numVar,numConstr);
        

        //2: Now generate the variable names
        enableVariableNames = true;
        varNames = new std::string[numVar]; //The deleting of this array is handled by the tableau.
        //2.1: First, load the original variables
        for (size_t i = 0; i < initNumVar; i++) {
            varNames[i] = initVarNames[i];
            varIDs[varNames[i]] = i;
        }
        //2.2: Now create the augmented negative variables from when a variable does not have a bounds.
        for (size_t i = 0; i < initNumVar; i++) {
            if (augmented_OldToNew_Index[i] != -1) {
                varNames[augmented_OldToNew_Index[i]] = "-(" + varNames[i] + ")";
            }
        }

        //2.2: Now create the slack variables. Every constraint has a slack variable so they'll just be indexed in sequence starting after numBaseVar.
        std::string baseSlackVarName = "SLACK_VAR_";
        for (size_t i = 0; i < numSlackSurplus; i++) {
            varNames[numBaseVar + numAugmentedVariables + i] = baseSlackVarName + std::to_string(i);
        }
        //2.3: Now create the artificial variables. These will be indexed in sequence following the slack variables.
        std::string baseArtVarNAme = "ARTIF_VAR_";
        for (size_t i = 0; i < numGEQ; i++) {
            varNames[numBaseVar + numAugmentedVariables + numConstr + i] = baseArtVarNAme + std::to_string(i);
        }
        //2.4: If artificial variables are present, the last variable will the objective of the original problem.
        if (numGEQ > 0) {
            varNames[numBaseVar + numAugmentedVariables + numConstr + numGEQ] = "ORIG_OBJ_FUNC";
        }

        tableau.AddVarNames(varNames);


        augmented_NewToOld_Index = new int[numVar];
        for (size_t i = 0; i < numVar; i++) {
            augmented_NewToOld_Index[i] = -1;
        }
        for (size_t oldI = 0; oldI < initNumVar; oldI++) {
            if (augmented_OldToNew_Index[oldI] != -1) {
                augmented_NewToOld_Index[augmented_OldToNew_Index[oldI]] = oldI;
            }
        }

        duplicatedConstr_NewToOld_Index = new int[numConstr];
        for (size_t i = 0; i < numConstr; i++) {
            duplicatedConstr_NewToOld_Index[i] = -1;
        }
        for (size_t oldI = 0; oldI < initNumConstr; oldI++) {
            if (duplicatedConstr_OldToNew_Index[oldI] != -1) {
                duplicatedConstr_NewToOld_Index[duplicatedConstr_OldToNew_Index[oldI]] = oldI;
            }
        }


        //3: Now populate the talbleau body.

        /* The tableau will look like this:
            |~ Original Vars ~|~ Augmented Vars ~|~ Slack/Surplus Vars ~|~ Artificial Vars ~|~ Original Objective Var ~|~ Constant Column ~|      Explaination:
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        A        |        A'        |         +- I         |         I         |            0             |         b         |  Original Model Equations with appropraite augmented, slack, surplus, and artificial variables added.
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        A''      |        A'''      |         +- I'        |         I'        |            0             |         b'        |  Duplicated Model Equations for replacing each EQ constraint with a LEQ and a GEQ constraint.
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        I        |        0         |           I          |         0         |            0             |       Bounds      |  Variable Bound Enforcement
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        c        |        c'        |           0          |         0         |            1             |        c[-1]      |  Original Objective Function
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |      sum(A)     |     sum(A')      |       -1 or 0        |         0         |            0             |       sum(c)      |  Auxilary Objective Function
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|

            Where A, b, and c are the original A, b, and c Matrix / vectors.
            A' and c' are the negative A and c coefficients for each augmented variable.
            I indicates the identity matrix.
            A'', A''', I', and b' are the appropriate variations of A, I, and b for each equality constraint that was duplicated into a GEQ and LEQ constraint.

            I'll handle each region below. Regions are labeled as follows:
            |~ Original Vars ~|~ Augmented Vars ~|~ Slack/Surplus Vars ~|~ Artificial Vars ~|~ Original Objective Var ~|~ Constant Column ~|
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        11       |        12        |           13         |         14        |            15            |         16        |
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        21       |        22        |           23         |         24        |            25            |         26        |
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        31       |        32        |           33         |         34        |            35            |         36        |
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        41       |        42        |           43         |         44        |            45            |         46        |
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
            |        51       |        52        |           53         |         54        |            55            |         56        |
            |-----------------|------------------|----------------------|-------------------|--------------------------|-------------------|
        */

        //3.1 Original Model Eqations
        //3.1.1 REGION 11, Original A matrix (accounting for shifted and inverted variables)
        for (size_t i = 0; i < initA.size(); i++) {
            //Cast the one-dimensional index, i, into two dimensional indices ii and jj
            size_t ii = i / initNumVar;
            size_t jj = i - ii;

            tableau[ii,jj] = initA[i]; //Shift and invert operations were written to initA earlier.
        }

        //3.1.2 REGION 12, Augmented variables model equation coefs.
        size_t jStart = initNumVar;
        size_t jStop = initNumVar + numAugmentedVariables; //numEQ added since this region (along with region 23 should be square)
        for (size_t newJ = jStart; newJ < jStop; j++) {
            size_t oldJ = (size_t) augmented_NewToOld_Index[newJ];
            for (size_t i = 0; i < initNumConstr; i++) {
                tableau[i,newJ] = -tableau[i,oldJ]; //Negative since the augmented variable (being a positive Scalar) represents the negative values of the original variable.
            }
        }

        //3.1.3 REGION 13, Slack/Surplus variables in the original model equations
        jStart = initNumVar + numAugmentedVariables;
        jStop = jStart + numSlackSurplus; //numEQ added since this region (along with region 23 should be square)
        for (size_t i = 0; i < initNumConstr; i++) {
            for (size_t j = jStart; j < jStop; j++) { 
                tableau[i,j] = 0; //One repeated "assignment" operation is faster than n "if" operations. So I'll just assign all zeros and then re-asign the one non-zero value.
            }

            if (constrTypes[i] == COSNTR_GEQ) {
                tableau[i,jStart + i] = -1; 
            }
            else {
                //The LEQ produced from each EQ will be listed in the original position. The GEQ produced will be placed in Regions 21-26.
                tableau[i,jStart + i] = 1;
            }
        }

        //3.1.4 REGION 14, Artificial variables in the original model equations.
        
        jStart = initNumVar + numAugmentedVariables + numSlackSurplus;
        jStop = jStart + numArtificialVars;
        size_t numGEQEncountered = 0;
        for (size_t i = 0; i < initNumConstr; i++) {
            for (size_t j = jStart; j < jStop; j++) {
                tableau[i,j] = 0; //One repeated "assignment" operation is faster than n "if" operations. So I'll just assign all zeros and then re-asign the one non-zero value.
            }

            if (constrTypes[i] == COSNTR_GEQ) {
                tableau[i,jStart + numGEQEncountered] = 1; 
                numGEQEncountered++;
            }
        }

        //3.1.5 REGION 15, Original Objective Function coefs in the original model equations.
        if (!auxProblemSolved) {
            size_t j = initNumVar + numAugmentedVariables + numSlackSurplus+ numArtificialVars;
            for (size_t i = 0; i < initNumConstr; i++) {
                tableau[i,j] = 0;
            }
        }

        //3.1.6 REGION 16, Constants in the original model equations (adjusted by variable shifts and inverts, handled earlier)
        size_t j = tableauWidth - 1;
        for (size_t i = 0; i < initNumConstr; i++) {
            tableau[i,j] = initC[i];
        }

        //3.2 Duplated Equality Constraints
        //3.2.1 REGION 21 & 22, Coefs of the original and augmented variables in the model equality constraints (adjusted for variable shifts and inverts) that are duplicated.
        size_t iStart = initNumConstr;
        size_t iStop = initNumConstr + numEQ;
        for (size_t newI = iStart; newI < iStop; newI++) {
            size_t oldI = (size_t) duplicatedConstr_NewToOld_Index[newI];
            for (size_t j = 0; j < initNumVar + numAugmentedVariables; j++) {
                tableau[newI,j] = tableau[oldI,j];
            }
        }

        //3.2.2 REGION 23, Surplus Variables for each of the duplicated equality constraints.
        jStart = initNumVar + numAugmentedVariables;
        jStop = initNumVar + numAugmentedVariables + initNumConstr + numEQ;
        for (size_t i = iStart; i < iStop; i++) {
            for (size_t j = jStart; j < jStop; j++) {
                tableau[i,j] = 0;
            }
            tableau[i,jStart + i] = -1; //Since all of these duplicated constraints are GEQ constraints.
        }

        //3.2.3 REGION 24, Artificial variables for each of the duplicated equality constraints.
        //Same iStart and iStop as REGION 23
        jStart = initNumVar + numAugmentedVariables + initNumConstr + numEQ;
        jStop = jStart + numArtificialVars;
        for (size_t i = iStart; i < iStop; i++) {
            for (size_t j = jStart; j < jStop; j++) {
                tableau[i,j] = 0;
            }
            tableau[i,jStart + numGEQEncountered] = 1; //Since all of these constraints are GEQ constraints.
            numGEQEncountered++;
        }

        //3.2.4 REGION 25 Original objective function variable in the duplicated equality constraints.
        if (!auxProblemSolved) {
            size_t j = initNumVar + initNumConstr + numEQ + numArtificialVars;
            for (size_t i = iStart; i < iStop; i++) {
                tableau[i,j] = 0;
            }
        }

        //3.2.5 REGION 26 Constant column for the duplicated equality constraints.
        j = tableauWidth - 1;
        for (size_t newI = iStart; newI < iStop; newI++) {
            size_t oldI = (size_t) duplicatedConstr_NewToOld_Index[newI];
            talbeau[newI,j] = tableau[oldI,j];
        }

        //3.3 Variable Bound Enforcement (REGIONs 31-36)
        size_t numUpperBoundsEncountered = 0;
        iStart = initNumConstr + numEQ;
        jStart = initNumVar + numAugmentedVariables + initNumConstr + numEQ;
        for (size_t i = 0; i < initNumVar; i++) { 
            //Remember, an augmented variable will only exist if the variable (after shifting and inverting) has no upper bound. So we don't need to inlude augmented variables here.
            //Additionally, after shifting, inteverting, and/or augmenting each variable, all variables' lower bound is 0 which is assumed by standard form.
            //So we don't need to explicitly state any lower bounds here.
            if (initBoundActivations[i].second) { //If this variable has an upper bound.
                size_t rowIndex = iStart + numUpperBoundsEncountered;
                auto rowItr = tableau.RowIterator(rowIndex);
                while (!rowItr.isEnd()) {
                    *rowItr = 0;
                    rowItr++;
                }

                Scalar upperBound = initBounds[i].second;
                tableau[rowIndex,i] = 1;

                size_t slackVarIndex = jStart + numUpperBoundsEncountered;
                tableau[rowIndex,slackVarIndex] = 1;

                talbeau[rowIndex,tableau.tableauWidth - 1] = upperBound;
            }
        }

        //3.4 Orginal Objective Function (REGIONs 41-46)
        

        

        


        //EVERYTHING FROM HERE DOWN IS OLD CODE.


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
        //  LEQ constraints have a 1 in the diagonal and 0s everywhere else
        //  GEQ constraints have a -1 in the diagonal and 0s everywhere else
        //Also we'll fill in the artificial variable values.
        //  LEQ constraints have zeros everywhere
        //  GEQ constraints have a 1 in the "diagonal" and zeros everywhere else
        numEQEncountered = 0;
        numGEQEncountered = 0;
        for (size_t i = 0; i < numConstr; i++) {
            size_t ii = i - numEQEncountered;
            size_t constrType = constrTypes[ii];
            if (constrType == CONSTR_EQ) {
                for (size_t iii = 0; iii < numConstr; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + i] = 1;

                for (size_t iii = 0; iii < numArtificialVars; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + numConstr + iii] = 0;
                }

                i++;
                numEQEncountered++;
                
                for (size_t iii = 0; iii < numConstr; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + i] = -1;

                for (size_t iii = 0; iii < numArtificialVars; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + numConstr + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + numConstr + numGEQEncountered] = 1;
                numGEQEncountered++;
            }
            else if (constrType == CONSTR_LEQ) {
                for (size_t iii = 0; iii < numConstr; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + i] = 1;

                for (size_t iii = 0; iii < numArtificialVars; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + numConstr + iii] = 0;
                }
            }
            else {
                for (size_t iii = 0; iii < numConstr; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + iii] = -1;

                for (size_t iii = 0; iii < numArtificialVars; iii++) {
                    tableauBody[tableauWidth * i + numBaseVar + numConstr + iii] = 0;
                }
                tableauBody[tableauWidth * i + numBaseVar + numConstr + numGEQEncountered] = 1;
                numGEQEncountered++;
            }
        }

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
            for (size_t i = 0; i < numBaseVar; i++) {
                Scalar columnSum = 0;
                for (size_t j = 0; j < numConstr; j++) {
                    columnSum += tableauBody[j * tableauWidth + i];
                }

                tableauBody[(tableauHeight - 1) * tableauWidth + i] = columnSum;
            }

            //The Slack variable columns will be -1 if the associated constraint is a geq constraint, 0 otherwise.
            numEQEncountered = 0;
            for (size_t i = 0; i < numConstr; i++) {
                if (i - numEQEncountered == CONSTR_EQ) {
                    tableauBody[(tableauHeight - 1) * tableauWidth + numBaseVar + i] = 0;
                    i++;
                    numEQEncountered++;
                    tableauBody[(tableauHeight - 1) * tableauWidth + numBaseVar + i] = -1;
                }
                else if (i - numEQEncountered == CONSTR_LEQ) {
                    tableauBody[(tableauHeight - 1) * tableauWidth + numBaseVar + i] = 0;
                }
                else {
                    tableauBody[(tableauHeight - 1) * tableauWidth + numBaseVar + i] = -1;
                }
            }

            //The art. var. columns will be 0
            for (size_t i = 0; i < numArtificialVars; i++) {
                tableauBody[(tableauHeight - 1) * tableauWidth + numBaseVar + numConstr + i] = 0;
            }

            //The original objective variable will have a value of zero
            tableauBody[tableauHeight * tableauWidth - 2] = 0;

            //The constant column will be the sum of elements of the original constant column that correspond to GEQ constriants.
            Scalar columnSum = 0;
            numEQEncountered = 0;
            for (size_t i = 0; i < initNumConstr; i++) {
                if (constrTypes[i] == CONSTR_LEQ) {}
                else {
                    columnSum += initC[i];
                }
            }
            tableauBody[tableauHeight * tableaWidth - 1] = columnSum;


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

            delete[] shiftValues;
            delete[] invertStatus;
            delete[] augmented_OldToNew_Index;
            delete[] augmented_NewToOld_Index;
            delete[] duplicatedConstr_OldToNew_Index;
            delete[] duplicatedConstr_NewToOld_Index;
        }
    }

    void setMaxIter(size_t newMaxIter) {
        maxIter = newMaxIter;
    }

    void setMaxTime(Scalar newMaxTime) {
        maxTime = newMaxTime;
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