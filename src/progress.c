#include <math.h>
#include <stdio.h>

#include "progress.h"

static void print_time(time_t t)
{
	if(t < 2 * 60) {
		printf("%ds", (int)t);
	}
	else if(t < 5 * 60) {
		printf("%dm%ds", (int)(t/60), (int)(t%60));
	}
	else if(t < 2 * 60 * 60) {
		printf("%dm", (int)(t/60));
	}
	else if(t < 24 * 60 * 60) {
		printf("%dh%dm", (int)(t/(60*60)), (int)((t/60)%60));
	}
	else {
		printf("%dd%dh", (int)(t/(24*60*60)), (int)((t/(60*60))%(24)));
	}
}

void print_progress(time_t started, int progress, int total, int decimals)
{
	char percent_format[32];
	time_t now, elapsed, remaining;
	float t;
	int i;

	now = time(NULL);
	elapsed = now - started;

	// Find appropriate percentage format
	snprintf(percent_format, 32, "%%.%df%%%%", decimals);

	// Calculations...
	if(progress > 1) {
    t = ceil(elapsed / (float)progress * (float)(total - progress));
		remaining = (time_t)t;
	}
	else {
		remaining = 0;
	}

	printf(percent_format, 100.0 * progress / total, elapsed);

	// Print elapsed time
	fputs(" (", stdout);
	print_time(elapsed);
	fputs(" elapsed, est. ", stdout);
	print_time(remaining);
	fputs(" remaining)", stdout);

	for(i = 0; i < 20; i++) {
		fputc(' ', stdout);
	}

	fputc('\r', stdout);
	fflush(stdout);
}
