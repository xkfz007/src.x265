/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com
 *****************************************************************************/

#ifndef X265_THREADPOOL_H
#define X265_THREADPOOL_H

#include "common.h"

namespace x265 {
// x265 private namespace

class ThreadPool;

int getCpuCount();

// Any class that wants to distribute work to the thread pool must
// derive from JobProvider and implement FindJob().
class JobProvider
{
protected:

    ThreadPool   *m_pool;

    JobProvider  *m_nextProvider;
    JobProvider  *m_prevProvider;

public:

    JobProvider(ThreadPool *p) : m_pool(p), m_nextProvider(0), m_prevProvider(0) {}

    virtual ~JobProvider() {}

    void setThreadPool(ThreadPool *p) { m_pool = p; }

    // Register this job provider with the thread pool, jobs are available
    void enqueue();

    // Remove this job provider from the thread pool, all jobs complete
    void dequeue();

    // Worker threads will call this method to find a job.  Must return true if
    // work was completed.  False if no work was available.
    virtual bool findJob(int threadId) = 0;

    // All derived objects that call Enqueue *MUST* call flush before allowing
    // their object to be destroyed, otherwise you will see random crashes involving
    // partially freed vtables and you will be unhappy
    void flush();

    friend class ThreadPoolImpl;
    friend class PoolThread;
};

// Abstract interface to ThreadPool.  Each encoder instance should call
// AllocThreadPool() to get a handle to the singleton object and then make
// it available to their job provider structures (wave-front frame encoders,
// etc).
class ThreadPool
{
protected:

    // Destructor is inaccessable, force the use of reference counted Release()
    ~ThreadPool() {}

    virtual void enqueueJobProvider(JobProvider &) = 0;

    virtual void dequeueJobProvider(JobProvider &) = 0;

public:

    // When numthreads == 0, a default thread count is used. A request may grow
    // an existing pool but it will never shrink.
    static ThreadPool *allocThreadPool(int numthreads = 0);

    static ThreadPool *getThreadPool();

    virtual void pokeIdleThread() = 0;

    // The pool is reference counted so all calls to AllocThreadPool() should be
    // followed by a call to Release()
    virtual void release() = 0;

    virtual int  getThreadCount() const = 0;

    friend class JobProvider;
};
} // end namespace x265

#endif // ifndef X265_THREADPOOL_H
