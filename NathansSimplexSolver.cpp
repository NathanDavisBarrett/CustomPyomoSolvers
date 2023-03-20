//from https://www.youtube.com/watch?v=_5T70cAXDJ0 Stopped at 2:08

#include<pybind11/pybind11.h>
#include<pybind11/stl.h>

#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<map>
#include<limits>
#include<time.h>

//from https://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
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

class SimplexSolver {
private:
    Scalar* tableauBody = NULL;
    size_t numVar = 0;
    size_t numConstr = 0;
    std::string* varNames = NULL;
    int* basisVars = NULL;
    size_t tableauHeight = 0;
    size_t tableauWidth = 0;
    std::map<std::string, std::string> options;
    std::map<std::string, std::vector<std::string>> allowedOptions;
    bool modelEngaged = false;
    time_t tic_time;
    Scalar SCALAR_INF = 0;
    size_t SIZE_T_INF = 0;
    size_t bottomRowStart = 0;
    size_t bottomRowEnd = 0;
    size_t lastColStart = 0;
    size_t lastColEnd = 0;
    
public:
    SimplexSolver() {
        SCALAR_INF = std::numeric_limits<Scalar>::max();
        SIZE_T_INF = std::numeric_limits<size_t>::max();
        SetupOptions();
    }

    ~SimplexSolver() {
        DisengageModel();
    }

    void EngageModel(
        size_t& initNumVar,
        size_t& initNumConstr,
        std::vector<std::string>& initVarNames,
        std::vector<Scalar>& initA,
        std::vector<Scalar>& initB,
        std::vector<Scalar>& initC,
        std::vector<int>& initBasis
        ) {
        if (modelEngaged) {
            DisengageModel();
        }

        ASSERT(initNumVar == initVarNames.size(),"ERROR! The number of variable names provided does not match the specified number of variables.");
        ASSERT(initA.size() == initNumVar * initNumConstr,"ERROR! The size of the specified A matrix does not match the number of variables and constraints specified.");
        ASSERT(initB.size() == initNumConstr,"ERROR! The size of the specified b vector does not match the specified number of constraints.");
        ASSERT(initC.size() == initNumVar + 1,"ERROR! The size of the specified c vector does not match the specified number of variables.");
        ASSERT(initBasis.size() <= initNumConstr,"ERROR! The size of the spcified basis is larger than the number of constraints. This is not allowed.");

        numVar = initNumVar;
        numConstr = initNumConstr;

        varNames = new std::string[numVar];
        for (size_t i = 0; i < numVar; i++) {
            varNames[i] = initVarNames[i];
        }

        tableauHeight = numConstr + 1;
        tableauWidth = numVar + 1;
        tableauBody = new Scalar[tableauWidth * tableauHeight];

        for (size_t i = 0; i < initA.size(); i++) {
            size_t ii = i / numVar;
            size_t jj = i - ii * numVar;

            size_t j = ii * tableauWidth + jj;
            tableauBody[j] = initA[i];
        }

        size_t rightColIndex = tableauWidth - 1;
        for (size_t i = 0; i < numConstr; i++) {
            tableauBody[rightColIndex] = initB[i];
            rightColIndex += tableauWidth;
        }

        size_t bottomRowIndex = tableauWidth * numConstr;
        for (size_t i = 0; i < tableauWidth; i++) {
            if (i == tableauWidth - 1) {
                tableauBody[i + bottomRowIndex] = initC[i];
            }
            else{
                tableauBody[i + bottomRowIndex] = -initC[i];
            }
        }

        basisVars = new int[numConstr];
        for (size_t i = 0; i < numConstr; i++) {
            if (i < initBasis.size()) {
                basisVars[i] = initBasis[i];
            }
            else {
                basisVars[i] = -1;
            }
        }

        bottomRowStart = tableauWidth * numConstr;
        bottomRowEnd = bottomRowStart + numVar;

        lastColStart = tableauWidth - 1;
        lastColEnd = tableauHeight * tableauWidth;

        modelEngaged = true;
    }

    void DisengageModel() {
        if (modelEngaged) {
            delete[] varNames;
            delete[] tableauBody;
            delete[] basisVars;

            tableauBody = NULL;
            numVar = 0;
            numConstr = 0;
            varNames = NULL;
            basisVars = NULL;
            tableauHeight = 0;
            tableauWidth = 0;

            modelEngaged = false;
        }
    }

    void SetupOptions() {
        allowedOptions["Basis_IO_Approach"] = {"MaximizeRC"};
        options["Basis_IO_Approach"] = "MaximizeRC";

        allowedOptions["maxIter"] = {"inf","anyPositiveInt"};
        options["maxIter"] = "inf";

        allowedOptions["maxTime"] = {"inf","anyPositiveFloat"};
        options["maxTime"] = "inf";
    }

    std::string TableauToString() {
        return TableauToMarkdownString();
    }

    std::string TableauToMarkdownString() {
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

        return ss.str();
    }

    std::string TableauToCSVString() {
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

        return ss.str();
    }

    std::map<std::string, std::string> GetCurrentOptions() {
        return options;
    }

    template <typename T>
    T StringToNumeric(const std::string& str) {
        std::stringstream ss(str);
        T result;
        ss >> result;
        return result;
    }

    void CheckOption(const std::string& id, const std::string& op) {
        if (allowedOptions.find(id) == allowedOptions.end()) {
            ASSERT(false,"Error! \"" + id + "\" is not a recognized setting.");
        }
        std::vector<std::string> validOptions = allowedOptions[id];

        bool valid = std::find(validOptions.begin(), validOptions.end(), op) != validOptions.end();
        if (!valid) {
            bool handled = false;
            if (std::find(validOptions.begin(), validOptions.end(), "anyPositiveInt") != validOptions.end()) {
                try {
                    int test = StringToNumeric<int>(op); //FIXME: This might not allways throw an error. For example, it it's a float, then this will be fine. Or, if it's a string that happens to have an int in it, then it's fine.
                    if (test >= 0) {
                        handled = true;
                    }
                }
                catch (const std::invalid_argument& e) {}
            }
            if (std::find(validOptions.begin(), validOptions.end(), "anyPositiveFloat") != validOptions.end()) {
                try {
                    Scalar test = StringToNumeric<Scalar>(op); //FIXME: Similar thing here.
                    if (test >= 0) {
                        handled = true;
                    }
                }
                catch (const std::invalid_argument& e) {}
            }

            if (handled) {
                return;
            }

            std::stringstream message;
            message << "Error! \"" << op << "\" is not a valid option for \"" << id << "\". Valid options are [";
            for (size_t i = 0; i < validOptions.size(); i++) {
                message << "\"" << validOptions[i] << "\"";
                if (i != validOptions.size() - 1) {
                    message << ", ";
                }
            }
            message << "]";
            ASSERT(false, message.str());
        }
    }

    void SetOptions(std::map<std::string, std::string>& newOptions) {
        for (auto i : newOptions) {
            CheckOption(i.first,i.second);
            options[i.first] = i.second;
        }
    }

    bool LessThanZero(const Scalar& val, const Scalar& tollerance = 1e-7) {
        return (val - tollerance) < 0.0;
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
        tic_time = time(NULL);
    }

    Scalar toc() {
        return difftime(time(NULL), tic_time);
    }

    bool GetOptimalityStatus() {
        bool optimal = true;
        for (size_t i = bottomRowStart; i < bottomRowEnd; i++) {
            if (LessThanZero(tableauBody[i])) {
                optimal = false;
                break;
            }
        }
        return optimal;
    }

    size_t ComputePivotColumn_MaximizeRC() {
        size_t columnIndex = 0;
        Scalar minBottomVal = tableauBody[bottomRowStart];

        size_t i = bottomRowStart + 1;
        for (size_t j = 1; j < numVar; i++, j++) {
            if (tableauBody[i] < minBottomVal) {
                minBottomVal = tableauBody[i];
                columnIndex = j;
            }
        }

        return columnIndex;
    }

    size_t ComputePivotRow_MaximizeRC(const size_t& pivotCol) {
        size_t rowIndex = 0;

        size_t lastRowIndex = lastColStart;
        size_t pivotColIndex = pivotCol;
        Scalar minComparisonVal = tableauBody[lastRowIndex] / tableauBody[pivotColIndex];

        if (LessThanZero(minComparisonVal)) { 
            minComparisonVal = SCALAR_INF; // Since we're only interested in the smallest positive value.
        }

        lastRowIndex += tableauWidth;
        pivotColIndex += tableauWidth;
        for (size_t j = 1; j < numConstr; j++, lastRowIndex += tableauWidth, pivotColIndex += tableauWidth) {
            if (EqualsZero(tableauBody[pivotColIndex])) {
                continue; //This will return an infinite result which, under no circumstances, is the best.
            }
            
            Scalar comparisonVal = tableauBody[lastRowIndex] / tableauBody[pivotColIndex];
            

            if (LessThanZero(comparisonVal)) { // Since we're only interested in the smallest positive value.
                continue;
            }

            if (comparisonVal < minComparisonVal) {
                rowIndex = j;
                minComparisonVal = comparisonVal;
            }
        }

        ASSERT(minComparisonVal != SCALAR_INF, "ERROR! No positive comparison values were found in the computation of the pivot row. This indicates that the model is infeasible??");
        ASSERT(SCALAR_INF == SCALAR_INF,"ERROR! This type of scalar comparison is not allowed! I'll need to make an \"IS_INF\" function and use it in the line above"); //FIXME

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
    
    void Solve() {
        //Make sure model is engaged before proceeding
        ASSERT(modelEngaged,"ERROR! Solver cannot execute if there is no model engaged.");

        //Make sure basis is complete. If not, run pre-solve.
        bool basisIsComplete = true;
        for (size_t i = 0; i < numConstr; i++) {
            if (basisVars[i] == -1) {
                basisIsComplete = false;
                break;
            }
        }

        ASSERT(basisIsComplete,"ERROR! I haven't coded how to handle an incomplete basis yet.");

        bool solutionIsOptimal = GetOptimalityStatus();

        size_t maxIter;
        if (options["maxIter"] == "inf") {
            maxIter = SIZE_T_INF;
        }
        else {
            maxIter = StringToNumeric<size_t>(options["maxIter"]);
        }

        Scalar maxTime;
        if (options["maxTime"] == "inf") {
            maxTime = SCALAR_INF;
        }
        else {
            maxTime = StringToNumeric<Scalar>(options["maxTime"]);
        }

        size_t (SimplexSolver::*ComputePivotColumn)() const = NULL;
        size_t (SimplexSolver::*ComputePivotRow)(size_t) const = NULL;


        std::string basisIOApproach = options["Basis_IO_Approach"];
        if (basisIOApproach == "MaximizeRC") {
            ComputePivotColumn = &(this->ComputePivotColumn_MaximizeRC);
            ComputePivotRow = &(this->ComputePivotRow_MaximizeRC);
        }
        else {
            ASSERT(false,"Error! Basis_IO_Approach not recognized. This is a bug.");
        }


        size_t iterationNum = 0;
        tic();
        size_t exitCode = 0;
        while (!solutionIsOptimal) {
            if (iterationNum > maxIter) {
                exitCode = 1;
                break;
            }
            iterationNum++;
            if (toc() > maxTime) {
                exitCode = 2;
                break;
            }

            size_t pivotCol = ComputePivotColumn();
            size_t pivotRow = ComputePivotRow(pivotCol);
            PerformPivot(pivotCol,pivotRow);

        }
    };

};

PYBIND11_MODULE(nathans_Simplex_solver_py, handle) {
    handle.doc() = "The C++ implementation of Nathan's Simplex Solver.\n To operate, instantiate a SimplexSolver object by passing the following as arguments: numVar, numConstr, varNames, the A matrix (flattened using A.flatten() or A.flatten('C')), the b vector, the vector.\nThen specify any neccecary solver options.\nThen execute the solving operation using the Solve() method.";
    handle.def("DummyFunc", &DummyFunc);

    py::class_<SimplexSolver>(handle, "SimplexSolver")
        .def(py::init<>())
        .def("EngageModel", &SimplexSolver::EngageModel)
        .def("Solve", &SimplexSolver::Solve)
        .def("GetCurrentOptions", &SimplexSolver::GetCurrentOptions)
        .def("SetOptions", &SimplexSolver::SetOptions)
        .def("TableauToString", &SimplexSolver::TableauToString);
}