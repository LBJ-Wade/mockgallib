#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "comm.h"
#include "msg.h"


static LogLevel log_level= msg_debug;
static char prefix[8]= "";

void msg_set_loglevel(const LogLevel lv)
{
  log_level= lv;

  msg_printf(msg_verbose, "MPI %d nodes\n", comm_n_nodes());
}

void msg_set_prefix(const char prefix_[])
{
  // Start messages with given prefix such as '#'
  strncpy(prefix, prefix_, 7);
}


void msg_printf(const LogLevel msg_level, char const * const fmt, ...)
{
  if(comm_this_rank() == 0 && msg_level >= log_level) {
    va_list argp;

    va_start(argp, fmt);

    fprintf(stdout, "%s", prefix);
    vfprintf(stdout, fmt, argp);
    fflush(stdout);

    va_end(argp);
  }
}

void msg_abort(char const * const fmt, ...)
{
  va_list argp;

  if(log_level <= msg_fatal) {
    va_start(argp, fmt);
    vfprintf(stdout, fmt, argp); fflush(stdout);
    va_end(argp);
  }

  abort();
}

