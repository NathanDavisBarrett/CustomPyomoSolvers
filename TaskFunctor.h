#ifndef TASK_FUNCTOR_H
#define TASK_FUNCTOR_H

#include <atomic>

class TaskFunctor {
public:
    bool complete = false;
    bool error = false;
    virtual void Execute() = 0;
};

void WaitForTasksToBeDone(TaskFunctor* tasks, size_t numTasks) {
    for (size_t i = 0; i < numTasks; i++) {
        while (!tasks[i].complete) {
            //Infinite loop till it's complete.
        }
    }
}

template <typename T>
class GetMaxMinOfContiguousArray: public TaskFunctor {
public:
    T* baseArray;
    size_t indexStart;
    size_t indexStop;
    size_t stepSize;
    bool max;

    std::atomic<size_t>* indexToBeat;

    GetMaxMinOfContiguousArray(T* initBaseArray, std::atomic<size_t>* initIndexToBeat, bool initMax = true, size_t initIndexStart, size_t initIndexStop, size_t initStepSize = 1) {
        baseArray = initBaseArray;
        indexToBeat = initIndexToBeat;

        indexStart = initIndexStart;
        indexStop = initIndexStop;

        stepSize = initStepSize;
        max = initMax;
    }

    void Execute() {
        if (max) {
            for (size_t i = indexStart + 1; i < indexStop; i += stepSize) {
                if (baseArray[i] > baseArray[*indexToBeat]) {
                    *indexToBeat = i;
                }
            }
        }
        else {
            for (size_t i = indexStart + 1; i < indexStop; i += stepSize) {
                if (baseArray[i] < baseArray[*indexToBeat]) {
                    *indexToBeat = i;
                }
            }
        }
    }
}

#endif