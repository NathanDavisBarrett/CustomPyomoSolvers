#ifndef INITIATE_THREADS_H
#define INITIATE_THREADS_H

#include "concurrentqueue.h"

#include "TaskFunctor.h"

#include <thread>
#include <mutex>

#define THREAD_STANDBY 0
#define THREAD_ACTIVE 1
#define THREAD_WAITING 2
#define THREAD_TERMINATE 3
#define THREAD_ERROR 4

void ThreadFunction(moodycamel::ConcurrentQueue<TaskFunctor*>* taskQueue, int* masterSignal , int* status, size_t threadID) {
    TaskFunctor* myTask = NULL;

    while (*masterSignal != THREAD_TERMINATE) {                  
        if (*masterSignal == THREAD_STANDBY) {
            *status = THREAD_STANDBY;
            continue;
        }

        //Try to grab a job from the taskQueue.
        if(!taskQueue->try_dequeue(myTask)) {
            //This cannot hang. Otherwise it will go into deadlock.
            *status = THREAD_WAITING;
            continue;
        }
        
        *status = THREAD_ACTIVE;

        //Execute that task
        try {
            myTask->Execute();
            myTask->complete = true;
        }
        catch (...) {
            *status = THREAD_ERROR;
            myTask->error = true;
        }
        myTask = NULL;
    }

    *status = THREAD_TERMINATE;
}

class ThreadManager {
public:
    size_t numThreads;
    moodycamel::ConcurrentQueue<TaskFunctor*> taskQueue;

    std::thread* threads;
    int* threadStatuses;

    int masterSignal = THREAD_STANDBY;

    ThreadManager(size_t initNumThreads) {
        masterSignal = THREAD_ACTIVE;

        numThreads = initNumThreads;

        threads = new std::thread[numThreads];
        threadStatuses = new int[numThreads];

        for (size_t i = 0; i < numThreads; i++) {
            threadStatuses[i] = THREAD_WAITING;
            threads[i] = std::thread(ThreadFunction,&taskQueue,&masterSignal,&(threadStatuses[i]),i);
            while (threadStatuses[i] != THREAD_WAITING) {
                //Wait for the thread to fully activate.
            }
        }
    }
    ThreadManager(): ThreadManager((size_t) std::thread::hardware_concurrency()) {}

    void PauseAll() {
        masterSignal = THREAD_STANDBY;
    }

    void ResumeAll() {
        masterSignal = THREAD_ACTIVE;
    }

    void TerminateThreads() {
        masterSignal = THREAD_TERMINATE;
        for (size_t i = 0; i < numThreads; i++) {
            while (threadStatuses[i] != THREAD_TERMINATE) {
                //Inifinite loop until the thread registers the terminate status.
            }
            threads[i].join();
        }
    }

    void AddTask(TaskFunctor* newTask) {
        taskQueue.enqueue(newTask);
    }

    ~ThreadManager() {
        TerminateThreads();
        delete[] threads; 
        delete[] threadStatuses;
    }
};

#endif