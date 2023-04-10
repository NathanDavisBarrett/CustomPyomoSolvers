#ifndef TASK_FUNCTOR_H
#define TASK_FUNCTOR_H

#include <atomic>
#include <unistd.h>
#include <sstream>

class TaskFunctor {
public:
    bool complete = false;
    bool error = false;
    virtual void Execute() = 0;
};

template <typename T>
class GetMaxMinOfContiguousArray: public TaskFunctor {
public:
    T* baseArray;
    size_t indexStart;
    size_t indexStop;
    size_t stepSize;
    bool max;

    std::atomic<size_t>* indexToBeat;

    GetMaxMinOfContiguousArray() {
        baseArray = NULL;
        indexStart = 0;
        indexStop = 0;
        stepSize = 0;
        max = false;
        indexToBeat = NULL;

        complete = false;
        error = false;
    }

    GetMaxMinOfContiguousArray(T* initBaseArray, std::atomic<size_t>* initIndexToBeat, size_t initIndexStart, size_t initIndexStop, bool initMax = true, size_t initStepSize = 1) {
        baseArray = initBaseArray;
        indexToBeat = initIndexToBeat;

        indexStart = initIndexStart;
        indexStop = initIndexStop;

        stepSize = initStepSize;
        max = initMax;

        complete = false;
        error = false;
    }

    void Execute() {
        if (max) {
            for (size_t i = indexStart; i < indexStop; i += stepSize) {
                if (baseArray[i] > baseArray[*indexToBeat]) {
                    *indexToBeat = i;
                }
            }
        }
        else {
            for (size_t i = indexStart; i < indexStop; i += stepSize) {
                if (baseArray[i] < baseArray[*indexToBeat]) {
                    *indexToBeat = i;
                }
            }
        }
        complete = true;
    }
};

#endif