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

#ifndef X265_WINXP_H
#define X265_WINXP_H

#if defined(_WIN32) && (_WIN32_WINNT < 0x0600) // _WIN32_WINNT_VISTA

#ifdef _MSC_VER
#include <intrin.h> // _InterlockedCompareExchange64
#endif

namespace x265 {
/* non-native condition variable */
typedef struct
{
    CRITICAL_SECTION broadcastMutex;
    CRITICAL_SECTION waiterCountMutex;
    HANDLE semaphore;
    HANDLE waitersDone;
    volatile int waiterCount;
    volatile int bIsBroadcast;
} ConditionVariable;

int WINAPI cond_init(ConditionVariable *cond);
void WINAPI cond_broadcast(ConditionVariable *cond);
void WINAPI cond_signal(ConditionVariable *cond);
BOOL WINAPI cond_wait(ConditionVariable *cond, CRITICAL_SECTION *mutex, DWORD wait);
void cond_destroy(ConditionVariable *cond);

/* map missing API symbols to our structure and functions */
#define CONDITION_VARIABLE          x265::ConditionVariable
#define InitializeConditionVariable x265::cond_init
#define SleepConditionVariableCS    x265::cond_wait
#define WakeConditionVariable       x265::cond_signal
#define WakeAllConditionVariable    x265::cond_broadcast
#define XP_CONDITION_VAR_FREE       x265::cond_destroy

} // namespace x265

#else // if defined(_WIN32) && (_WIN32_WINNT < 0x0600)

#define XP_CONDITION_VAR_FREE(x)

#endif // _WIN32_WINNT <= _WIN32_WINNT_WINXP

#endif // ifndef X265_WINXP_H
