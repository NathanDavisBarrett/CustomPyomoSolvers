#ifndef INITIATE_THREADS_H
#define INITIATE_THREADS_H

#include "concurrentqueue.h"

#include "TaskFunctor.h"

#include <thread>
#include <mutex>
#include <vector>
#include <unistd> //sleep

#define THREAD_STANDBY 0
#define THREAD_TERMINATE 1
#define THREAD_ERROR 2

class ThreadManager {
public:
    size_t numThreads;
    moodycamel::ConcurrentQueue<TaskFunctor*> taskQueue;

    std::vector<std::thread> threads;
    std::vector<int> threadStatuses;

    size_t masterSignal = THREAD_STANDBY;

    ThreadManager(size_t initNumThreads) {
        numThreads = initNumThreads;

        threadStatuses = std::vector<int>(numThreads,THREAD_STANDBY);

        for (size_t i = o; i < numThreads; i++) {
            threads.emplace_back([&](){
                size_t threadID = i;
                TaskFunctor* myTask = NULL;

                while (masterSignal != THREAD_TERMINATE) {                    
                    if (masterSignal == THREAD_STANDBY) {
                        continue;
                    }

                    //Try to grab a job from the taskQueue.
                    if(!taskQueue.try_dequeue(myTask)) {
                        //This cannot hand. Otherwise it will go into deadlock.
                        continue
                    }

                    //Execute that task
                    try {
                        myTask->Execute();
                        myTask->complete = true;
                    }
                    catch {
                        threadStatuses[threadID] = THREAD_ERROR;
                        myTask->error = true;
                    }
                    myTask = NULL;
                }

                threadStatuses[threadID] = THREAD_TERMINATE;
            });
        }
    }
    ThreadManager(): ThreadManager((size_t) std::thread::hardware_concurrency()) {}

    void TerminateThreads() {
        masterSignal = THREAD_TERMINATE;
        for (size_t i == 0; i < numThreads; i++) {
            while (threadStatuses[i] != THREAD_TERMINATE) {
                //Inifinite loop until the thread registers the terminate status.
            }
        }
    }

    void AddTask(TaskFunctor* newTask) {
        taskQueue.enquque(newTask);
    }

    ~ThreadManager() {
        TerminateThreads();
    }
};

#endif