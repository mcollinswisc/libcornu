#ifndef __PROGRESS_H__
#define __PROGRESS_H__

#include <time.h>

/**
 * Update a command-line progress meter.
 *
 * @TODO Calculate decimals within this
 */
void print_progress(time_t started, int progress, int total, int decimals);

#endif // #ifndef __PROGRESS_H__
